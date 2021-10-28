library(ggplot2)
library(dplyr)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

df.in_cs <-
  rgsutil::read_gsfile(fm_insights$get_merged_results_path("shard_trait.pip001.pop", "tsv.bgz"))
df.af_pop <-
  rgsutil::read_gsfile(fm_insights$get_analysis_path("af_pop.fin_jpn", "tsv.bgz"))

df.traits <- read_trait_summary()

n_samples <-
  dplyr::mutate(df.traits,
                n_samples = ifelse(
                  is.na(n_cases),
                  n_samples,
                  n_samples * (n_cases / n_controls) * (1 - (n_cases / n_controls))
                )) %>%
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
      pip_cohort = dplyr::case_when(pip_threshold > 0 ~ pip_cohort,
                                    in_cs_cohort ~ 1,
                                    TRUE ~ -1),
      maf_base = 0.5 - abs(0.5-!!as.symbol(paste0("af.", cohort1))),
      maf_base = ifelse(is.na(maf_base), 0, maf_base),
      maf_cohort = 0.5 - abs(0.5-!!as.symbol(paste0("af.", cohort2))),
      maf_cohort = ifelse(is.na(maf_cohort), 0, maf_cohort),
      beta_posterior_base = !!as.symbol(paste0("susie.beta_posterior.", cohort1)),
      var_base = 2 * maf_base * (1 - maf_base),
      var_cohort = 2 * maf_cohort * (1 - maf_cohort),
      power_base = pchisq(
        threshold,
        1,
        ncp = n_samples_base * var_base * (beta_posterior_base ** 2),
        lower.tail = F
      ),
      power_cohort = pchisq(
        threshold,
        1,
        ncp = n_samples_cohort * var_cohort * (beta_posterior_base ** 2),
        lower.tail = F
      ),
      pvalue_base = !!as.symbol(paste0("pvalue.", cohort1)),
      pvalue_cohort = !!as.symbol(paste0("pvalue.", cohort2)),
      pip_bin = cut(pip_base, pip_bin_breaks),
      gw_sig = pvalue_cohort < pvalue_threshold &
        !is.na(pvalue_cohort),
      status = factor(
        dplyr::case_when(
          pip_base <= 0.9 | is.na(pip_base) ~ "EXCLUDE",
          (pvalue_cohort < pvalue_threshold) &
            (pip_cohort <= pip_threshold |
               is.na(pip_cohort)) ~ "Fine-mapping non-replication",
          pip_cohort > pip_threshold ~ "Fine-mapping replication",
          pvalue_cohort < 0.01 ~ "Replicated association",
          power_cohort > 0.9 ~ "Non-replicated association",!is.na(pvalue_cohort) ~ "Low power",
          is.na(pvalue_cohort) ~ "Missing variants",
          TRUE ~ "ERROR"
        ),
        levels = c(
          "Fine-mapping replication",
          "Fine-mapping non-replication",
          "Replicated association",
          "Non-replicated association",
          "Low power",
          "Missing variants",
          "EXCLUDE",
          "ERROR"
        )
      )
    )
}

assign_status_any <- function(df.in_cs,
                              df.traits,
                              pip_threshold = 0.1,
                              sim.comparison = FALSE) {
  data =
    purrr::map_dfr(cohorts, function(cohort) {
      other_cohorts <- setdiff(cohorts, cohort)
      purrr::map_dfr(other_cohorts, function(other_cohort) {
        assign_status(df.in_cs,
                      df.traits,
                      cohort,
                      other_cohort,
                      pip_threshold = pip_threshold) %>%
          dplyr::filter(status != "EXCLUDE")
      })
    })
  if (sim.comparison) {
    data = dplyr::mutate(data, trait = stringr::str_c(trait, cohort, sep = ":"))
  }
  dplyr::group_by(data, variant, trait) %>%
    dplyr::summarize(
      base = "Any",
      cohort = "Any",
      count = n(),
      discovery = unique(c(c("BBJ")[any(pip.BBJ > 0.9, na.rm = TRUE)],
                           c("FG")[any(pip.FG > 0.9, na.rm = TRUE)],
                           c("UKBB")[any(pip.UKBB > 0.9, na.rm = TRUE)])),
      discovery = factor(
        ifelse(length(discovery) > 1, "Multi", discovery),
        levels = c(cohorts, "Multi")
      ),
      all_status = stringr::str_c(unique(status), collapse = ","),
      gw_sig = any(gw_sig, na.rm = TRUE),
      status = factor(
        dplyr::case_when(
          any(status == "Fine-mapping replication") ~ "Fine-mapping replication",
          any(status == "Fine-mapping non-replication") ~ "Fine-mapping non-replication",
          any(status == "Replicated association") ~ "Replicated association",
          any(status == "Non-replicated association") ~ "Non-replicated association",
          any(status == "Low power") ~ "Low power",
          any(status == "Missing variants") ~ "Missing variants",
          TRUE ~ "ERROR"
        ),
        levels = c(
          "Fine-mapping replication",
          "Fine-mapping non-replication",
          "Replicated association",
          "Non-replicated association",
          "Low power",
          "Missing variants",
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
          dplyr::filter(status != "EXCLUDE")
      })
    print(dplyr::group_by(data, base, cohort, status) %>% dplyr::count())
  }
  
  data <- dplyr::mutate(data,
                        cohort = stringr::str_c(base, cohort, sep = " "))
  
  cohort.levels <-
    dplyr::group_by(data, cohort) %>%
    dplyr::summarize(frac = sum(status == "Fine-mapping replication" &
                                  gw_sig) / sum(gw_sig)) %>%
    dplyr::arrange(frac) %>%
    dplyr::pull(cohort)
  
  data <-
    dplyr::mutate(data, cohort = factor(cohort, levels = cohort.levels))
  
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
      scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0),
                         labels = scales::percent) +
      theme(
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      ) +
      locusviz::or_missing(!axis.x,
                           theme(axis.title.x = element_blank(),
                                 axis.text.x = element_blank())) +
      locusviz::or_missing(!axis.y,
                           theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"))) +
      labs(y = "% variant", fill = "Replication status") +
      scale_fill_manual(
        values = c(
          `Fine-mapping replication` = PNWColors::pnw_palette("Bay")[1],
          `Fine-mapping non-replication` = BuenColors::jdb_palette("brewer_blue")[2],
          `Replicated association` = PNWColors::pnw_palette("Bay", 9)[6],
          `Non-replicated association` = PNWColors::pnw_palette("Bay", 9)[5],
          `Low power` = "grey90",
          `Missing variants` = "grey70"
        )
      ) +
      scale_color_manual(
        values = c(
          `Fine-mapping replication` = ifelse(any(data2$gw_sig), "white", NA),
          `Fine-mapping non-replication` = "black",
          `Replicated association` = "black",
          `Non-replicated association` = "black",
          `Low power` = "black",
          `Missing variants` = "black"
        ),
        guide = FALSE,
        na.translate = FALSE
      ) +
      coord_flip()
  }
  return(list(
    plot_bar(dplyr::filter(data, gw_sig), dplyr::filter(data2, gw_sig)),
    plot_bar(
      dplyr::filter(data,!gw_sig),
      dplyr::filter(data2,!gw_sig),
      axis.y = FALSE
    )
  ))
}

plt <- list(
  plot_status_bar(df.in_cs, plot.any = TRUE, axis.x = F),
  plot_status_bar(df.in_cs, plot.any = FALSE)
) %>%
  purrr::flatten() %>%
  purrr::reduce(`+`) + patchwork::plot_layout(nrow = 2,
                                              widths = c(9, 11),
                                              heights = c(1, 6))
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
                                axis.x = TRUE,
                                sim.comparison = FALSE) {
  require(GGally)
  if (plot.any) {
    data <-
      assign_status_any(df, df.traits, pip_threshold = pip_threshold)
    cohort.levels <- "Any Any"
  } else {
    data <-
      gtools::permutations(length(cohorts), 2, cohorts) %>%
      magrittr::set_colnames(c("cohort", "other_cohort")) %>%
      tibble::as_tibble() %>%
      purrr::pmap_dfr(function(cohort, other_cohort) {
        assign_status(df, df.traits, cohort, other_cohort, pip_threshold = pip_threshold) %>%
          dplyr::filter(status != "EXCLUDE")
      })
    print(dplyr::group_by(data, base, cohort, status) %>% dplyr::count())
    
    cohort.levels <-
      gtools::permutations(length(cohorts), 2, cohorts) %>%
      magrittr::set_colnames(c("cohort", "other_cohort")) %>%
      tibble::as_tibble() %>%
      purrr::pmap_dfr(function(cohort, other_cohort) {
        assign_status(df, df.traits, cohort, other_cohort, pip_threshold = 0.1) %>%
          dplyr::filter(status != "EXCLUDE")
      }) %>%
      dplyr::mutate(cohort = stringr::str_c(base, cohort, sep = " ")) %>%
      dplyr::group_by(cohort) %>%
      dplyr::summarize(frac = sum(status == "Fine-mapping replication" &
                                    gw_sig) / sum(gw_sig)) %>%
      dplyr::arrange(frac) %>%
      dplyr::pull(cohort)
    print(cohort.levels)
  }
  
  data <- dplyr::mutate(data,
                        cohort = stringr::str_c(base, cohort, sep = " "))
  
  data <-
    dplyr::mutate(data, cohort = factor(cohort, levels = cohort.levels))
  
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
      scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0),
                         labels = scales::percent) +
      theme(
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      ) +
      locusviz::or_missing(!axis.x,
                           theme(axis.title.x = element_blank(),
                                 axis.text.x = element_blank())) +
      locusviz::or_missing(!axis.y,
                           theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"))) +
      labs(y = "% variant", fill = "Replication status") +
      scale_fill_manual(
        values = c(
          `Fine-mapping replication` = PNWColors::pnw_palette("Bay")[1],
          `Fine-mapping non-replication` = BuenColors::jdb_palette("brewer_blue")[2],
          `Replicated association` = PNWColors::pnw_palette("Bay", 9)[6],
          `Non-replicated association` = PNWColors::pnw_palette("Bay", 9)[5],
          `Low power` = "grey90",
          `Missing variants` = "grey70"
        )
      ) +
      scale_color_manual(
        values = c(
          `Fine-mapping replication` = ifelse(any(data2$gw_sig), "white", NA),
          `Fine-mapping non-replication` = "black",
          `Replicated association` = "black",
          `Non-replicated association` = "black",
          `Low power` = "black",
          `Missing variants` = "black"
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
  patchwork::plot_layout(ncol = 2,
                         byrow = FALSE,
                         heights = c(1, 6))

cowplot::save_plot(
  "figures/ExDataFig2/ExDataFig2_overlap_pip005_in_cs.pdf",
  plt2,
  base_height = 2.4,
  base_width = 6.7,
  device = cairo_pdf
)

#######################################################################

plot_violin <- function(df,
                        cohort,
                        axis.y = TRUE,
                        maf.threshold = 0) {
  require(GGally)
  other_cohorts <- setdiff(cohorts, cohort)
  data <- purrr::map_dfr(other_cohorts, function(other_cohort) {
    assign_status(df, df.traits, cohort, other_cohort)
  })
  
  data2 <-
    dplyr::filter(data,!is.na(pip_bin)) %>%
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
  
  pd <- position_dodge(width = 0.05)
  plt =
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
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::number_format(accuracy = 0.1)) +
    labs(x = sprintf("PIP bin (%s)", cohort), y = "PIP in a secondary population") +
    my_theme +
    locusviz::or_missing(!axis.y,
                         theme(axis.title.y = element_blank(),
                               axis.text.y = element_blank())) +
    theme(
      legend.position = c(0.07, 1),
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      plot.tag.position = c(-0.05, 1),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  return(plt)
}


plot_sim_violin <- function(df.in_cs, df.traits, df.sim) {
  data.real = purrr::map_dfr(cohorts, function(cohort) {
    other_cohorts <- setdiff(cohorts, cohort)
    purrr::map_dfr(other_cohorts, function(other_cohort) {
      assign_status(df.in_cs, df.traits, cohort, other_cohort)
    }) %>%
      dplyr::mutate(sec.cohort = cohort,
                    cohort = .env$cohort)
  }) %>%
    dplyr::filter(pvalue_cohort < 5e-8 & pip_bin == "(0.9,1]")
  
  data.real.point <-
    dplyr::group_by(data.real, cohort, sec.cohort) %>%
    dplyr::summarize(mean_pip = mean(pip_cohort, na.rm = T))
  
  data.sim = dplyr::filter(
    df.sim,
    cohort != sim.cohort &
      pip > 0.9 &
      pchisq(sim.chisq_marginal, 1, lower.tail = F) < 5e-8 &
      !is.na(sim.pip)
  )
  data.sim.point = dplyr::group_by(data.sim, cohort, sim.cohort) %>%
    dplyr::summarize(mean_pip = mean(sim.pip)) %>%
    dplyr::rename(sec.cohort = sim.cohort) %>%
    dplyr::mutate(sec.cohort = paste0(sec.cohort, ".sim"))
  
  pd <- position_dodge(width = 1)
  ggplot() +
    gghalves::geom_half_violin(
      aes(
        cohort,
        pip_cohort,
        fill = sec.cohort,
        group = interaction(cohort, sec.cohort)
      ),
      size = 0.2,
      alpha = 0.5,
      width = 1,
      scale = "width",
      side = "l",
      data = data.real
    ) +
    gghalves::geom_half_violin(
      aes(cohort, sim.pip, group = interaction(cohort, sim.cohort)),
      fill = "grey50",
      size = 0.2,
      alpha = 0.5,
      width = 1,
      scale = "width",
      side = "r",
      data = data.sim
    ) +
    geom_point(
      aes(cohort, mean_pip, color = sec.cohort),
      size = 1,
      position = pd,
      data = data.real.point
    ) +
    geom_point(
      aes(cohort, mean_pip, color = sec.cohort),
      size = 1,
      position = pd,
      data = data.sim.point
    ) +
    scale_color_manual(values = c(
      cohort_colors,
      c(
        "BBJ.sim" = "grey20",
        "UKBB.sim" = "grey20",
        "FG.sim" = "grey20"
      )
    )) +
    scale_fill_manual(values = cohort_colors, guide = FALSE) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::number_format(accuracy = 0.1)) +
    labs(x = "Discovery cohort (PIP > 0.9)", y = "PIP in a secondary population\n(left: real, right: simulated)") +
    my_theme +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      plot.tag.position = c(-0.05, 1)
    )
}

df.sim = rgsutil::read_gsfile(fm_insights$get_merged_results_path("sim.fm_only.pip09", "tsv.bgz")) %>%
  dplyr::group_split(locus, alleles, trait, cohort) %>%
  purrr::map_dfr(function(data) {
    n_sim_cohorts = sum(!is.na(data$sim.cohort))
    if (n_sim_cohorts == 2) {
      return(data)
    }
    sim_cohorts = setdiff(cohorts, c(data$cohort[1], unique(data$sim.cohort)))
    dplyr::bind_rows(
      data,
      purrr::map_dfr(seq_along(sim_cohorts), ~ data) %>%
        dplyr::mutate(
          sim.cohort = sim_cohorts,
          sim.chisq_marginal = NA,
          sim.pip = NA
        )
    ) %>% tidyr::drop_na(sim.cohort)
  })

plt <-
  list(
    plot_violin(df.in_cs, cohort = "BBJ"),
    plot_violin(df.in_cs, cohort = "FG", axis.y = F),
    plot_violin(df.in_cs, cohort = "UKBB", axis.y = F),
    plot_sim_violin(df.in_cs, df.traits, df.sim)
  ) %>%
  purrr::reduce(`+`) + patchwork::plot_layout(nrow = 1, widths = c(5, 5, 5, 5))

cowplot::save_plot(
  "figures/ExDataFig2/ExDataFig2_pip_distribution.pdf",
  plt,
  base_height = 2.4,
  base_width = 6.98,
  device = cairo_pdf
)

###################################

plot_sim_status_bar = function(df,
                               df.in_cs,
                               df.traits,
                               plot.any = FALSE,
                               axis.x = TRUE,
                               axis.y = TRUE) {
  data = dplyr::filter(df,
                       cohort != sim.cohort &
                         pip > 0.9) %>%
    dplyr::mutate(
      gw_sig = pchisq(sim.chisq_marginal, 1, lower.tail = F) < 5e-8,
      status = factor(
        dplyr::case_when(
          !gw_sig ~ "Not GW significant",
          sim.pip > 0.1 ~ "Fine-mapping replication",
          TRUE ~ "Fine-mapping non-replication"
        ),
        levels = c("Fine-mapping replication", "Fine-mapping non-replication")
      )
    ) %>%
    dplyr::filter(!(gw_sig & is.na(sim.pip)))
  print(dplyr::filter(data, is.na(sim.pip)))
  if (plot.any) {
    data = dplyr::mutate(data,
                         cohort = factor("Any Any"))
  } else {
    cohort.levels <-
      gtools::permutations(length(cohorts), 2, cohorts) %>%
      magrittr::set_colnames(c("cohort", "other_cohort")) %>%
      tibble::as_tibble() %>%
      purrr::pmap_dfr(function(cohort, other_cohort) {
        assign_status(df.in_cs,
                      df.traits,
                      cohort,
                      other_cohort,
                      pip_threshold = 0.1) %>%
          dplyr::filter(status != "EXCLUDE")
      }) %>%
      dplyr::mutate(cohort = stringr::str_c(base, cohort, sep = " ")) %>%
      dplyr::group_by(cohort) %>%
      dplyr::summarize(frac = sum(status == "Fine-mapping replication" &
                                    gw_sig) / sum(gw_sig)) %>%
      dplyr::arrange(frac) %>%
      dplyr::pull(cohort)
    print(cohort.levels)
    
    data = dplyr::mutate(data,
                         cohort = factor(stringr::str_c(cohort, sim.cohort, sep = " "), levels = cohort.levels))
  }
  
  data2 =
    dplyr::group_by(data, cohort) %>%
    dplyr::summarize(frac = mean(gw_sig))
  
  dplyr::filter(data, gw_sig) %>%
    ggplot(aes(cohort)) +
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
    scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0),
                       labels = scales::percent) +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    locusviz::or_missing(!axis.x,
                         theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank())) +
    locusviz::or_missing(!axis.y,
                         theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"))) +
    labs(y = "% variant", fill = "Replication status") +
    scale_fill_manual(
      values = c(
        `Fine-mapping replication` = BuenColors::jdb_palette("brewer_green")[7],
        `Fine-mapping non-replication` = BuenColors::jdb_palette("brewer_green")[2]
      )
    ) +
    scale_color_manual(
      values = c(
        `Fine-mapping replication` = "white",
        `Fine-mapping non-replication` = "black"
      ),
      guide = FALSE,
      na.translate = FALSE
    ) +
    coord_flip()
}


plt = list(
  plot_status_bar_ext(
    df.in_cs,
    plot.any = TRUE,
    axis.x = FALSE,
    pip_threshold = 0.1,
    sim.comparison = TRUE
  ),
  plot_sim_status_bar(
    df.sim,
    df.in_cs,
    df.traits,
    plot.any = TRUE,
    axis.x = FALSE
  ),
  plot_sim_status_bar(df.sim, df.in_cs, df.traits)
) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1, heights = c(1, 1, 6))

cowplot::save_plot(
  "figures/ExDataFig2/ExDataFig2_sim_overlap.pdf",
  plt,
  base_height = 2.4,
  base_width = 3.44,
  device = cairo_pdf
)

###########################################################################

plot_replication_status_csq = function(df.csq, df.in_cs, df.traits) {
  data = assign_status_any(df.in_cs, df.traits) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(df.csq) %>%
    dplyr::count(status, consequence) %>%
    dplyr::group_by(status) %>%
    dplyr::mutate(total = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(frac = n / total,
                  consequence = ordered(consequence, levels = rev(
                    c(
                      "Non-genic",
                      "CRE",
                      "Promoter",
                      "UTR3",
                      "UTR5",
                      "Synonymous",
                      "Missense",
                      "pLoF"
                    )
                  ))) %>%
    dplyr::arrange(desc(consequence)) %>%
    dplyr::group_by(status) %>%
    dplyr::mutate(pos = cumsum(frac))
  
  label.data <- dplyr::group_by(data, status) %>%
    dplyr::summarize(label = ifelse(
      total[1] >= 1e4,
      scales::label_number_si(accuracy = 1)(total[1]),
      as.character(total[1])
    ))
  
  plt =
    ggplot() +
    geom_bar(aes(
      x = factor(status),
      y = frac,
      fill = consequence
    ),
    stat = "identity",
    data = data) +
    geom_text(
      aes(
        x = factor(status),
        y = 1,
        label = label
      ),
      vjust = -0.5,
      size = 2,
      data = label.data
    ) +
    geom_vline(xintercept = 2.5,
               linetype = 'dashed',
               color = "grey50") +
    scale_fill_manual(values = annot_colors) +
    my_theme +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 8, margin = margin(b = 5.5)),
      plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank()
    ) +
    labs(x = "Replication status",
         y = "Proportion of variants",
         fill = "Annotation") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  return(plt)
}


df.csq = rgsutil::read_gsfile(fm_insights$get_merged_results_path("pip09.vtc", "tsv")) %>%
  dplyr::select(variant, consequence, most_severe, gene_most_severe) %>%
  dplyr::distinct()

plt = plot_replication_status_csq(df.csq, df.in_cs, df.traits)


cowplot::save_plot(
  "figures/ExDataFig2/ExDataFig2_replication_status_csq.pdf",
  plt + theme(plot.margin = margin(0.1, 0.1, 0.2, 0.5, unit = "cm")),
  base_height = 3.05,
  base_width = 3.35,
  device = cairo_pdf
)
