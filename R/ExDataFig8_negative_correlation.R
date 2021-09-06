library(ggplot2)
library(magrittr)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

df.csm <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("fm_only.csm_id", "tsv.bgz")) %>%
  dplyr::select(cohort, trait, region, variant, susie.cs_id, csm_id)

plot_cor_panel <- function(trait,
                           lead_variant,
                           window,
                           highlight_variants,
                           highlight_rect,
                           highlight_hlines,
                           df.csm,
                           cohorts = c("BBJ", "FG", "UKBB")) {
  lead_pos <- locusviz::parse_variant(lead_variant)$position
  chromosome <- locusviz::parse_variant(lead_variant)$chromosome
  if (length(window) == 1) {
    window <- rep(window, 2)
  }
  start <- lead_pos - window[1]
  end <- lead_pos + window[2]

  cs.colors <- dplyr::filter(df.csm, trait == .env$trait) %>%
    dplyr::mutate(locusviz::parse_variant(variant)) %>%
    dplyr::filter(start < position & position < end) %>%
    dplyr::arrange(abs(position - position[which(.$variant == lead_variant)[1]])) %$%
    locusviz::get_cs_color_mapping(.$csm_id, highlight_cs_ids = .$csm_id[which(.$variant == lead_variant)])

  df.v_in_cs <- purrr::map_dfr(cohorts, function(cohort) {
    rgsutil::read_gsfile(fm_insights$get_locus_txt_path(cohort, trait, lead_variant)) %>%
      dplyr::mutate(locusviz::parse_variant(variant)) %>%
      dplyr::filter(susie.cs_id > 0 &
        start < position & position < end) %>%
      dplyr::select(variant, rsid)
  }) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      rsid = dplyr::case_when(
        rsid == "5:176514781_ACAG_A" ~ "rs45608536",
        rsid == "16:56997349:C:CA" ~ "rs5817082",
        rsid == "16:57001254:T:TCACA" ~ "rs12720908",
        rsid == "16:57013743:TA:T" ~ "rs71383213",
        TRUE ~ rsid
      )
    )


  df.ld <- purrr::map_dfr(cohorts, function(cohort) {
    R <- rgsutil::read_gsfile(fm_insights$get_ld_matrix_path(cohort, lead_variant, "50000"),
      header = F
    ) %>%
      as.matrix()

    df <- rgsutil::read_gsfile(fm_insights$get_ld_txt_path(cohort, lead_variant)) %>%
      dplyr::mutate(locusviz::parse_variant(variant)) %>%
      dplyr::filter((lead_pos - 50000) < position &
        position < (lead_pos + 50000)) %>%
      dplyr::mutate(idx = dplyr::row_number()) %>%
      dplyr::select(-rsid) %>%
      dplyr::inner_join(df.v_in_cs)

    df2 <- dplyr::filter(
      df,
      variant %in% highlight_variants
    )

    R <- R[df2$idx, df$idx, drop = F]
    colnames(R) <- df$rsid
    as.data.frame(R) %>%
      dplyr::mutate(row = sprintf("%s: %s", cohort, df2$rsid)) %>%
      tidyr::pivot_longer(cols = !row, names_to = "col") %>%
      dplyr::left_join(dplyr::select(df, rsid, position), by = c(col = "rsid"))
  }) %>%
    dplyr::mutate(col = factor(col, levels = unique(col[order(position)])))


  p1 <-
    ggplot() +
    geom_rect(aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ), fill = "grey90") +
    geom_tile(aes(col, row, fill = value), data = df.ld) +
    geom_text(aes(col, row, label = sprintf("%.1g", value)),
      size = 1,
      data = df.ld
    ) +
    geom_hline(yintercept = highlight_hlines, color = "grey50") +
    geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      color = "grey20",
      linetype = "dotted",
      fill = NA,
      data = highlight_rect
    ) +
    scale_fill_distiller(
      palette = "RdBu",
      limits = c(-1, 1),
      oob = scales::squish,
      direction = 1
    ) +
    coord_equal() +
    my_theme +
    theme(
      legend.position = "right",
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    scale_x_discrete(expand = expansion(0)) +
    scale_y_discrete(
      limits = rev,
      expand = expansion(0)
    ) +
    labs(fill = "r")


  p2 <-
    dplyr::filter(df.csm, trait == .env$trait) %>%
    dplyr::distinct(cohort, variant, csm_id) %>%
    dplyr::right_join(df.v_in_cs) %>%
    dplyr::mutate(locusviz::parse_variant(variant)) %>%
    dplyr::mutate(
      cohort = sprintf("%s: 95%% CS", cohort),
      rsid = factor(rsid, levels = unique(rsid[order(position)])),
      csm_id = factor(csm_id)
    ) %>%
    ggplot(aes(rsid, cohort, fill = csm_id)) +
    geom_rect(aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ), fill = "grey90") +
    geom_tile() +
    coord_equal() +
    my_theme +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      )
    ) +
    scale_x_discrete(expand = expansion(0)) +
    scale_y_discrete(limits = rev, expand = expansion(0)) +
    scale_fill_manual(values = cs.colors)

  plt <- p1 / p2
  return(plt)
}

plot_forest_plot <- function(trait,
                             lead_variant,
                             highlight_variants,
                             cohorts = c("BBJ", "FG", "UKBB")) {
  data <-
    purrr::map_dfr(cohorts, function(cohort) {
      data <-
        rgsutil::read_gsfile(fm_insights$get_locus_txt_path(cohort, trait, lead_variant)) %>%
        dplyr::select(
          cohort,
          variant,
          rsid,
          beta_marginal,
          se_marginal,
          pip,
          susie.beta_posterior,
          susie.sd_posterior
        )
    }) %>%
    dplyr::filter(variant %in% highlight_variants) %>%
    dplyr::mutate(
      shape = !is.na(pip) & pip > 0.1,
      beta_posterior = susie.beta_posterior,
      se_posterior = susie.sd_posterior,
      xmin_marginal = beta_marginal - se_marginal,
      xmax_marginal = beta_marginal + se_marginal,
      xmin_posterior = beta_posterior - se_posterior,
      xmax_posterior = beta_posterior + se_posterior
    )

  pd <- position_dodge(width = 0.5)

  data.rect <-
    dplyr::distinct(data, rsid) %>%
    dplyr::mutate(idx = dplyr::row_number()) %>%
    dplyr::arrange(dplyr::desc(rsid)) %>%
    dplyr::mutate(
      ymin = idx - 0.5,
      ymax = idx + 0.5
    ) %>%
    dplyr::filter(idx %% 2 == 0)

  xmin <- min(c(data$beta_marginal, data$beta_posterior))
  xmax <- max(c(data$beta_marginal, data$beta_posterior))

  plt <- purrr::map(c("marginal", "posterior"), function(x) {
    ggplot() +
      geom_rect(
        aes(
          xmin = -Inf,
          xmax = Inf,
          ymin = ymin,
          ymax = ymax
        ),
        fill = "grey90",
        data = data.rect
      ) +
      geom_point(
        aes_string(
          sprintf("beta_%s", x),
          "rsid",
          color = "cohort",
          shape = "shape"
        ),
        size = 1,
        position = pd,
        data = data
      ) +
      geom_errorbarh(
        aes_string(
          y = "rsid",
          xmin = sprintf("xmin_%s", x),
          xmax = sprintf("xmax_%s", x),
          color = "cohort"
        ),
        height = 0,
        position = pd,
        data = data
      ) +
      scale_y_discrete(limits = rev) +
      scale_color_manual(values = cohort_colors) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4)) +
      my_theme +
      theme(
        legend.position = if (x == "marginal") {
          c(1, 0.05)
        } else {
          "none"
        },
        legend.justification = c(1, 0),
        axis.title.y = element_blank(),
        plot.margin = margin(0.1, 0.1, 0, 0.1, unit = "cm")
      ) +
      labs(x = sprintf("%s beta", stringr::str_to_sentence(x))) +
      guides(
        color = guide_legend(title = "Cohort"),
        shape = guide_legend(title = "PIP > 0.1")
      ) +
      coord_cartesian(xlim = c(xmin, xmax))
  }) %>%
    purrr::reduce(`+`) + patchwork::plot_layout(ncol = 1)
  return(plt)
}


p1 <- plot_multi_locuszoom(
  trait = "Height",
  lead_variant = "5:176509193:C:T",
  rsid = "rs244711",
  window = 20000,
  df.csm = df.csm,
  additional_highlight_pos = 176516631,
  additional_highlight_labels = "rs1966265"
)

p2 <- plot_cor_panel(
  trait = "Height",
  lead_variant = "5:176509193:C:T",
  window = 20000,
  df.csm = df.csm,
  highlight_variants = c("5:176509193:C:T", "5:176516631:G:A"),
  highlight_rect = tibble::tibble(
    xmin = c(0.5, 10.5),
    xmax = c(1.5, 11.5),
    ymin = c(4.5, 4.5),
    ymax = c(6.5, 6.5)
  ),
  highlight_hlines = c(2.5, 4.5)
)
p2[[1]] <- p2[[1]] + theme(legend.position = "none")


p3 <- plot_forest_plot(
  trait = "Height",
  lead_variant = "5:176509193:C:T",
  highlight_variants = c("5:176509193:C:T", "5:176516631:G:A")
)
p3[[1]] <- p3[[1]] + theme(legend.position = "none")


p4 <- plot_multi_locuszoom(
  trait = "HDLC",
  lead_variant = "16:57017662:G:A",
  rsid = "rs1801706",
  window = c(22000, 2500),
  cohorts = c("BBJ", "UKBB"),
  df.csm = df.csm,
  additional_highlight_pos = c(57017292, 57016150),
  additional_highlight_labels = c("rs2303790", "rs5742907")
)
p4[[1]] <- p4[[1]] + theme(legend.position = "none")

p5 <- plot_cor_panel(
  trait = "HDLC",
  lead_variant = "16:57017662:G:A",
  window = c(22000, 2500),
  highlight_variants = c("16:57017662:G:A", "16:57017292:A:G", "16:57016150:G:A"),
  highlight_rect = tibble::tibble(
    xmin = 16.5,
    xmax = 19.5,
    ymin = 2.5,
    ymax = 5.5
  ),
  highlight_hlines = 2.5,
  cohorts = c("BBJ", "UKBB"),
  df.csm = df.csm
)
# p5[[1]] = p5[[1]] + theme(legend.position = "none")

p6 <- plot_forest_plot(
  trait = "HDLC",
  lead_variant = "16:57017662:G:A",
  highlight_variants = c("16:57017662:G:A", "16:57017292:A:G", "16:57016150:G:A"),
  cohorts = c("BBJ", "UKBB")
)

cowplot::save_plot(
  "figures/ExDataFig8/SFig_locuszoom_FGFR4_CETP.pdf",
  p1 | p4,
  base_height = 3,
  base_width = 7.2,
  device = cairo_pdf
)

cowplot::save_plot(
  "figures/ExDataFig8/SFig_ld_FGFR4_CETP.pdf",
  p2 | p5,
  base_height = 2,
  base_width = 7.2,
  device = cairo_pdf
)

cowplot::save_plot(
  "figures/ExDataFig8/SFig_forest_plot_FGFR4_CETP.pdf",
  p3 | p6 + patchwork::plot_layout(guides = "collect"),
  base_height = 2,
  base_width = 4.2,
  device = cairo_pdf
)