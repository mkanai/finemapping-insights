library(ggplot2)
library(dplyr)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

df.in_cs <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("shard_trait.pip001.pop", "tsv.bgz"))
df.af_pop <- rgsutil::read_gsfile(fm_insights$get_analysis_path("af_pop.fin_jpn", "tsv.bgz"))

df.traits <- read_trait_summary()

n_samples <-
  dplyr::mutate(df.traits,
    n_samples = ifelse(
      is.na(n_cases),
      n_samples,
      n_samples * (n_cases / n_controls) * (1 - (n_cases / n_controls))
    )
  ) %>%
  dplyr::select(cohort, trait, n_samples) %>%
  tidyr::pivot_wider(
    id_cols = "trait",
    names_from = "cohort",
    names_prefix = "n_samples.",
    values_from = "n_samples"
  )

all.traits <- df.traits %>%
  dplyr::group_by(trait) %>%
  dplyr::filter(length(unique(cohort)) == 3) %>%
  dplyr::pull(trait) %>%
  unique()

df.in_cs <- dplyr::filter(df.in_cs, trait %in% all.traits) %>%
  dplyr::left_join(df.af_pop)

alpha <- 0.01
threshold <- qchisq(alpha, 1, lower.tail = F)

#######################################################################
# Figure 2. overlapping summary
#######################################################################
assign_status <- function(df,
                          df.traits,
                          cohort1,
                          cohort2,
                          pip_threshold = 0.1) {
  # if pip_theshold < 0, use `in CS` alternatively

  pvalue_threshold <- 5e-8
  shared.traits <-
    dplyr::filter(df.traits, cohort %in% c(cohort1, cohort2)) %>%
    dplyr::group_by(trait) %>%
    dplyr::filter(length(unique(cohort)) == 2) %>%
    .$trait %>%
    unique()

  dplyr::filter(df, trait %in% shared.traits) %>%
    dplyr::left_join(n_samples, by = "trait") %>%
    # susie or finemap failure
    dplyr::filter(!is.na(!!as.symbol(paste0(
      "finemap.pip.", cohort1
    ))) &
      !is.na(!!as.symbol(paste0(
        "susie.pip.", cohort1
      )))) %>%
    dplyr::mutate(
      base = cohort1,
      cohort = cohort2,
      n_samples_base = !!as.symbol(paste0("n_samples.", cohort1)),
      n_samples_cohort = !!as.symbol(paste0("n_samples.", cohort2)),
      pip_base = !!as.symbol(paste0("pip.", cohort1)),
      pip_cohort = !!as.symbol(paste0("pip.", cohort2)),
      pip_cohort_failed = is.na(pip_cohort) & (!is.na(!!as.symbol(
        paste0("finemap.pip.", cohort2)
      )) |
        !is.na(!!as.symbol(
          paste0("susie.pip.", cohort2)
        ))),
      in_cs_base = !!as.symbol(paste0("susie.cs_id.", cohort1)),
      in_cs_cohort = !!as.symbol(paste0("susie.cs_id.", cohort2)),
      in_cs_base = !is.na(in_cs_base) & (in_cs_base > 0),
      in_cs_cohort = !is.na(in_cs_cohort) & (in_cs_cohort > 0),
      pip_cohort = dplyr::case_when(
        pip_threshold > 0 ~ pip_cohort,
        in_cs_cohort ~ 1,
        TRUE ~ -1
      ),
      maf_base = 0.5 - abs(0.5 - !!as.symbol(paste0("af.", cohort1))),
      maf_base = ifelse(is.na(maf_base), 0, maf_base),
      maf_cohort = 0.5 - abs(0.5 - !!as.symbol(paste0("af.", cohort2))),
      maf_cohort = ifelse(is.na(maf_cohort), 0, maf_cohort),
      beta_posterior_base = !!as.symbol(paste0("susie.beta_posterior.", cohort1)),
      var_base = 2 * maf_base * (1 - maf_base),
      var_cohort = 2 * maf_cohort * (1 - maf_cohort),
      power_base = pchisq(
        threshold,
        1,
        ncp = n_samples_base * var_base * (beta_posterior_base**2),
        lower.tail = F
      ),
      power_cohort = pchisq(
        threshold,
        1,
        ncp = n_samples_cohort * var_cohort * (beta_posterior_base**2),
        lower.tail = F
      ),
      pvalue_base = !!as.symbol(paste0("pvalue.", cohort1)),
      pvalue_cohort = !!as.symbol(paste0("pvalue.", cohort2)),
      pip_bin = cut(pip_base, pip_bin_breaks),
      status = factor(
        dplyr::case_when(
          pvalue_cohort < 5e-8 ~ as.character(cut(pip_cohort, pip_bin_breaks2)),
          pvalue_cohort < 5e-8 ~ "P < 5e-8",
          pvalue_cohort < 1e-4 ~ "P < 1e-4",
          pvalue_cohort >= 1e-4 ~ "P >= 1e-4",
          is.na(pvalue_cohort) ~ "NA",
          TRUE ~ "ERROR"
        ),
        levels = c(
          "(0.9,1]",
          "(0.5,0.9]",
          "(0.1,0.5]",
          "(0.01,0.1]",
          "(-Inf,0.01]",
          "P < 5e-8",
          "P < 1e-4",
          "P >= 1e-4",
          "NA"
        )
      ),
      status = forcats::fct_recode(status, "[0,0.01]" = "(-Inf,0.01]"),
      gw_sig = pvalue_cohort < pvalue_threshold &
        !is.na(pvalue_cohort),
      category = factor(
        dplyr::case_when(
          pip_base <= 0.9 | is.na(pip_base) ~ "EXCLUDE",
          (pvalue_cohort < pvalue_threshold) &
            (pip_cohort <= pip_threshold |
              is.na(pip_cohort)) ~ "Fine-mapping uncertainity",
          pip_cohort > pip_threshold ~ "Replication",
          # pvalue_cohort < 5e-8 ~ "Fine-mapping uncertainity",
          pvalue_cohort < 0.01 ~ "Replicated association",
          power_cohort > 0.9 ~ "Heterogeneity", !is.na(pvalue_cohort) ~ "Low power",
          is.na(pvalue_cohort) ~ "Population-specific",
          TRUE ~ "ERROR"
        ),
        levels = c(
          "Replication",
          "Fine-mapping uncertainity",
          "Replicated association",
          "Heterogeneity",
          "Low power",
          "Population-specific",
          "EXCLUDE",
          "ERROR"
        )
      )
    )
}

assign_status_any <- function(df.in_cs, df.traits,
                              pip_threshold = 0.1) {
  purrr::map_dfr(cohorts, function(cohort) {
    other_cohorts <- setdiff(cohorts, cohort)
    data <- purrr::map_dfr(other_cohorts, function(other_cohort) {
      assign_status(df.in_cs,
        df.traits,
        cohort,
        other_cohort,
        pip_threshold = pip_threshold
      ) %>%
        dplyr::filter(category != "EXCLUDE") %>%
        dplyr::mutate(status = category)
    })
  }) %>%
    dplyr::group_by(variant, trait) %>%
    dplyr::summarize(
      base = "Any",
      cohort = "Any",
      count = n(),
      discovery = unique(c(
        c("BBJ")[any(pip.BBJ > 0.9, na.rm = TRUE)],
        c("FG")[any(pip.FG > 0.9, na.rm = TRUE)],
        c("UKBB")[any(pip.UKBB > 0.9, na.rm = TRUE)]
      )),
      discovery = factor(
        ifelse(length(discovery) > 1, "Multi", discovery),
        levels = c(cohorts, "Multi")
      ),
      all_status = stringr::str_c(unique(status), collapse = ","),
      gw_sig = any(gw_sig, na.rm = TRUE),
      status = factor(
        dplyr::case_when(
          any(status == "Replication") ~ "Replication",
          any(status == "Fine-mapping uncertainity") ~ "Fine-mapping uncertainity",
          any(status == "Replicated association") ~ "Replicated association",
          any(status == "Heterogeneity") ~ "Heterogeneity",
          any(status == "Low power") ~ "Low power",
          any(status == "Population-specific") ~ "Population-specific",
          TRUE ~ "ERROR"
        ),
        levels = c(
          "Replication",
          "Fine-mapping uncertainity",
          "Replicated association",
          "Heterogeneity",
          "Low power",
          "Population-specific",
          "ERROR"
        )
      )
    )
}

plot_status_bar <- function(df,
                            plot.any = FALSE,
                            axis.x = TRUE) {
  require(GGally)
  if (plot.any) {
    data <- assign_status_any(df, df.traits)
  } else {
    data <-
      gtools::permutations(length(cohorts), 2, cohorts) %>%
      magrittr::set_colnames(c("cohort", "other_cohort")) %>%
      tibble::as_tibble() %>%
      purrr::pmap_dfr(function(cohort, other_cohort) {
        assign_status(df, df.traits, cohort, other_cohort) %>%
          dplyr::filter(category != "EXCLUDE") %>%
          dplyr::mutate(status = category)
      })
    print(dplyr::group_by(data, base, cohort, category) %>% dplyr::count())
  }

  data <- dplyr::mutate(data,
    cohort = stringr::str_c(base, cohort, sep = " ")
  )

  cohort.levels <-
    dplyr::group_by(data, cohort) %>%
    dplyr::summarize(frac = sum(status == "Replication" &
      gw_sig) / sum(gw_sig)) %>%
    dplyr::arrange(frac) %>%
    dplyr::pull(cohort)

  data <- dplyr::mutate(data, cohort = factor(cohort, levels = cohort.levels))

  data2 <-
    dplyr::group_by(data, cohort, gw_sig) %>%
    dplyr::summarize(count = n()) %>%
    dplyr::group_by(cohort) %>%
    dplyr::mutate(frac = count / sum(count))

  plot_bar <- function(data, data2, axis.y = TRUE) {
    ggplot(data, aes(cohort)) +
      geom_bar(aes(fill = status), position = position_fill(reverse = TRUE)) +
      geom_text(
        aes(
          label = scales::percent(..prop.., accuracy = 1),
          # ..count..
          by = cohort,
          color = status
        ),
        stat = "prop",
        size = 2,
        position = position_fill(0.5, reverse = TRUE)
      ) +
      geom_text(
        aes(cohort, 1, label = scales::percent(frac, accuracy = 1)),
        size = 2,
        hjust = -0.5,
        data = data2
      ) +
      my_theme +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.1), add = 0),
        labels = scales::percent
      ) +
      theme(
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      ) +
      locusviz::or_missing(
        !axis.x,
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank()
        )
      ) +
      locusviz::or_missing(
        !axis.y,
        theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"))
      ) +
      labs(y = "% variant", fill = "Replication status") +
      scale_fill_manual(
        values = c(
          `Replication` = PNWColors::pnw_palette("Bay")[1],
          `Fine-mapping uncertainity` = BuenColors::jdb_palette("brewer_blue")[2],
          # `Replicated association` = "#A3BDCC",
          `Replicated association` = PNWColors::pnw_palette("Bay", 9)[6],
          `Heterogeneity` = PNWColors::pnw_palette("Bay", 9)[5],
          # `Low power` = shades::opacity(PNWColors::pnw_palette("Bay")[3], 0.5),
          `Low power` = "grey90",
          `Population-specific` = "grey70"
        )
      ) +
      scale_color_manual(
        values = c(
          `Replication` = ifelse(any(data2$gw_sig), "white", NA),
          `Fine-mapping uncertainity` = "black",
          `Replicated association` = "black",
          `Heterogeneity` = "black",
          `Low power` = "black",
          `Population-specific` = "black"
        ),
        guide = FALSE,
        na.translate = FALSE
      ) +
      coord_flip()
  }
  return(list(
    plot_bar(dplyr::filter(data, gw_sig), dplyr::filter(data2, gw_sig)),
    plot_bar(
      dplyr::filter(data, !gw_sig),
      dplyr::filter(data2, !gw_sig),
      axis.y = FALSE
    )
  ))
}

plt <- list(
  plot_status_bar(df.in_cs, plot.any = TRUE, axis.x = F),
  plot_status_bar(df.in_cs, plot.any = FALSE)
) %>%
  purrr::flatten() %>%
  purrr::reduce(`+`) + patchwork::plot_layout(
    nrow = 2,
    widths = c(9, 11),
    heights = c(1, 6)
  )
plt

cowplot::save_plot(
  "figures/Fig2/Fig2_overlap_summary_v6.pdf",
  plt,
  base_height = 2.4,
  base_width = 6.7,
  device = cairo_pdf
)

#######################################################################

plot_status_bar_ext <- function(df,
                                pip_threshold,
                                plot.any = FALSE,
                                axis.x = TRUE) {
  require(GGally)
  if (plot.any) {
    data <- assign_status_any(df, df.traits, pip_threshold = pip_threshold)
    cohort.levels <- "Any Any"
  } else {
    data <-
      gtools::permutations(length(cohorts), 2, cohorts) %>%
      magrittr::set_colnames(c("cohort", "other_cohort")) %>%
      tibble::as_tibble() %>%
      purrr::pmap_dfr(function(cohort, other_cohort) {
        assign_status(df, df.traits, cohort, other_cohort, pip_threshold = pip_threshold) %>%
          dplyr::filter(category != "EXCLUDE") %>%
          dplyr::mutate(status = category)
      })
    print(dplyr::group_by(data, base, cohort, category) %>% dplyr::count())

    cohort.levels <-
      gtools::permutations(length(cohorts), 2, cohorts) %>%
      magrittr::set_colnames(c("cohort", "other_cohort")) %>%
      tibble::as_tibble() %>%
      purrr::pmap_dfr(function(cohort, other_cohort) {
        assign_status(df, df.traits, cohort, other_cohort, pip_threshold = 0.1) %>%
          dplyr::filter(category != "EXCLUDE") %>%
          dplyr::mutate(status = category)
      }) %>%
      dplyr::mutate(cohort = stringr::str_c(base, cohort, sep = " ")) %>%
      dplyr::group_by(cohort) %>%
      dplyr::summarize(frac = sum(status == "Replication" &
        gw_sig) / sum(gw_sig)) %>%
      dplyr::arrange(frac) %>%
      dplyr::pull(cohort)
    print(cohort.levels)
  }

  data <- dplyr::mutate(data,
    cohort = stringr::str_c(base, cohort, sep = " ")
  )

  data <- dplyr::mutate(data, cohort = factor(cohort, levels = cohort.levels))

  data2 <-
    dplyr::group_by(data, cohort, gw_sig) %>%
    dplyr::summarize(count = n()) %>%
    dplyr::group_by(cohort) %>%
    dplyr::mutate(frac = count / sum(count))

  plot_bar <- function(data, data2, axis.y = TRUE) {
    ggplot(data, aes(cohort)) +
      geom_bar(aes(fill = status), position = position_fill(reverse = TRUE)) +
      geom_text(
        aes(
          label = scales::percent(..prop.., accuracy = 1),
          # ..count..
          by = cohort,
          color = status
        ),
        stat = "prop",
        size = 2,
        position = position_fill(0.5, reverse = TRUE)
      ) +
      geom_text(
        aes(cohort, 1, label = scales::percent(frac, accuracy = 1)),
        size = 2,
        hjust = -0.5,
        data = data2
      ) +
      my_theme +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.1), add = 0),
        labels = scales::percent
      ) +
      theme(
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      ) +
      locusviz::or_missing(
        !axis.x,
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank()
        )
      ) +
      locusviz::or_missing(
        !axis.y,
        theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"))
      ) +
      labs(y = "% variant", fill = "Replication status") +
      scale_fill_manual(
        values = c(
          `Replication` = PNWColors::pnw_palette("Bay")[1],
          `Fine-mapping uncertainity` = BuenColors::jdb_palette("brewer_blue")[2],
          # `Replicated association` = "#A3BDCC",
          `Replicated association` = PNWColors::pnw_palette("Bay", 9)[6],
          `Heterogeneity` = PNWColors::pnw_palette("Bay", 9)[5],
          # `Low power` = shades::opacity(PNWColors::pnw_palette("Bay")[3], 0.5),
          `Low power` = "grey90",
          `Population-specific` = "grey70"
        )
      ) +
      scale_color_manual(
        values = c(
          `Replication` = ifelse(any(data2$gw_sig), "white", NA),
          `Fine-mapping uncertainity` = "black",
          `Replicated association` = "black",
          `Heterogeneity` = "black",
          `Low power` = "black",
          `Population-specific` = "black"
        ),
        guide = FALSE,
        na.translate = FALSE
      ) +
      coord_flip()
  }
  return(plot_bar(dplyr::filter(data, gw_sig), dplyr::filter(data2, gw_sig)))
}

plt2 <-
  list(
    plot_status_bar_ext(
      df.in_cs,
      plot.any = TRUE,
      axis.x = FALSE,
      pip_threshold = 0.05
    ),
    plot_status_bar_ext(df.in_cs, plot.any = FALSE, pip_threshold = 0.05),
    # in 95% CS
    plot_status_bar_ext(
      df.in_cs,
      plot.any = TRUE,
      axis.x = FALSE,
      pip_threshold = -1
    ),
    plot_status_bar_ext(df.in_cs, plot.any = FALSE, pip_threshold = -1)
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(
    ncol = 2,
    byrow = FALSE,
    heights = c(1, 6)
  )

cowplot::save_plot(
  "figures/ExDataFig2/ExDataFig2_overlap_pip005_in_cs.pdf",
  plt2,
  base_height = 2.4,
  base_width = 6.7,
  device = cairo_pdf
)

#######################################################################

plot_summary <- function(df,
                         cohort,
                         axis.y = TRUE,
                         maf.threshold = 0) {
  require(GGally)
  other_cohorts <- setdiff(cohorts, cohort)
  data <- purrr::map_dfr(other_cohorts, function(other_cohort) {
    assign_status(df, df.traits, cohort, other_cohort) %>%
      dplyr::mutate(status = category)
  })

  data2 <-
    dplyr::filter(data, !is.na(pip_bin)) %>%
    dplyr::group_by(cohort, pip_bin) %>%
    dplyr::summarize(
      count_replicated_association = sum(pvalue_cohort < 0.01, na.rm = T),
      count_sig = sum(pvalue_cohort < 5e-8, na.rm = T),
      count_missing = sum(is.na(pvalue_cohort)),
      frac_sig = count_sig / (count_replicated_association + count_sig),
      mean_pip = sum(pip_cohort[which(pvalue_cohort < 5e-8 &
        !is.na(pip_cohort))]) / sum(pvalue_cohort < 5e-8 &
        !is.na(pip_cohort), na.rm = T),
      frac_missing = count_missing / n()
    )

  data3 <- dplyr::filter(data, pvalue_cohort < 5e-8 &
    !is.na(pip_bin)) %>%
    tidyr::drop_na(pip_cohort)

  print(dplyr::filter(data2, pip_bin == "(0.9,1]"))

  plot_line <- function(y, y_count, ylab, xtext = TRUE) {
    ggplot(data2, aes_string("pip_bin", y, color = "cohort")) +
      geom_line(aes(group = cohort)) +
      geom_point(size = 1) +
      # ggrepel::geom_text_repel(aes_string("pip_bin", y, label = y_count), size = 2, direction = "y", point.size = 1) +
      scale_color_manual(values = cohort_colors[other_cohorts]) +
      scale_y_continuous(
        lim = c(0, 0.7),
        labels = ifelse(
          y != "mean_pip",
          scales::percent_format(accuracy = 1),
          scales::number_format(accuracy = 0.1)
        )
      ) +
      labs(x = sprintf("PIP bin (%s)", cohort), y = ylab) +
      my_theme +
      locusviz::or_missing(
        !axis.y,
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank()
        )
      ) +
      theme(
        legend.position = c(0.05, 1),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        plot.tag.position = c(-0.05, 1),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      locusviz::or_missing(
        !xtext,
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm")
        )
      ) +
      locusviz::or_missing(y != "frac_missing", theme(legend.position = "none"))
  }

  plot_violinhalf <- function() {
    pd <- position_dodge(width = 0.05)
    ggplot() +
      gghalves::geom_half_violin(
        aes(pip_bin, pip_cohort, fill = cohort),
        size = 0.2,
        alpha = 0.5,
        scale = "width",
        side = "l",
        data = dplyr::filter(data3, cohort == other_cohorts[1])
      ) +
      gghalves::geom_half_violin(
        aes(pip_bin, pip_cohort, fill = cohort),
        size = 0.2,
        alpha = 0.5,
        scale = "width",
        side = "r",
        data = dplyr::filter(data3, cohort == other_cohorts[2])
      ) +
      geom_line(
        aes(pip_bin, mean_pip, color = cohort, group = cohort),
        position = pd,
        data = data2
      ) +
      geom_point(
        aes(pip_bin, mean_pip, color = cohort),
        size = 1,
        position = pd,
        data = data2
      ) +
      scale_color_manual(values = cohort_colors) +
      scale_fill_manual(values = cohort_colors, guide = FALSE) +
      scale_y_continuous(
        limits = c(0, 1),
        labels = scales::number_format(accuracy = 0.1)
      ) +
      labs(x = sprintf("PIP bin (%s)", cohort), y = "PIP in a secondary population") +
      my_theme +
      locusviz::or_missing(
        !axis.y,
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank()
        )
      ) +
      theme(
        legend.position = c(0.07, 1),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        plot.tag.position = c(-0.05, 1),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }

  return(list(
    # plot_line("frac_sig", "count_sig", "% GW-sig variants", xtext = F),
    plot_violinhalf()
    # plot_line("frac_missing", "count_missing", "% missing variants")
  ))
}


# plt = purrr::map2(c(
#   plot_summary(df.in_cs, cohort = "BBJ"),
#   plot_summary(df.in_cs, cohort = "FG", axis.y = F),
#   plot_summary(df.in_cs, cohort = "UKBB", axis.y = F)
# ), as.list(matrix(
#   letters[3:8], byrow = T, ncol = 3
# )), function(panel, tag) {
#   panel + labs(tag = tag)
# }) %>% purrr::reduce(`+`) + patchwork::plot_layout(ncol = 3, byrow = F)

df.sim <- rgsutil::read_gsfile(fm_insights$get_analysis_path("ukbb_sim.gamma_gwas.sig.pip", "tsv.bgz")) %>%
  dplyr::mutate(high_pip = prob > 0.9) %>%
  tidyr::drop_na(prob)

df.sim.label <-
  dplyr::group_by(df.sim, high_pip) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    x = ifelse(high_pip, 0.9, 0.5),
    y = 100,
    label = scales::percent(count / sum(count), accuracy = 1)
  )

p_sim_gamma <-
  ggplot(df.sim, aes(prob)) +
  geom_histogram(aes(fill = high_pip), binwidth = 0.02, center = 0.01) +
  geom_text(aes(x, y, label = label), data = df.sim.label, size = 2) +
  my_theme +
  scale_fill_manual(values = c(`FALSE` = "grey50", `TRUE` = "lightblue"), guide = FALSE) +
  labs(x = "PIP", y = "# true causal variants")

plt <-
  list(
    plot_summary(df.in_cs, cohort = "BBJ"),
    plot_summary(df.in_cs, cohort = "FG", axis.y = F),
    plot_summary(df.in_cs, cohort = "UKBB", axis.y = F),
    list(p_sim_gamma)
  ) %>%
  purrr::flatten() %>%
  purrr::reduce(`+`) + patchwork::plot_layout(ncol = 3)
plt

cowplot::save_plot(
  "figures/ExDataFig2/ExDataFig2_pip_distribution.pdf",
  plt,
  base_height = 5.0,
  base_width = 7.02,
  device = cairo_pdf
)