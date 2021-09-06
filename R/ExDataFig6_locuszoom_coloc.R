library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

params <- list(
  "15:58723939:G:A" = list(
    trait = "HDLC",
    lead_variant = "15:58723939:G:A",
    rsid = "rs2070895",
    cohorts = c("BBJ", "UKBB"),
    tissue = "Liver",
    symbol = "LIPC",
    gtex_label = "GTEx: Liver (LIPC)",
    tag = "a",
    window = 250000
  ),
  "17:7571752:T:G" = list(
    trait = "SkC",
    trait_label = "Skin cancer",
    lead_variant = "17:7571752:T:G",
    rsid = "rs78378222",
    cohorts = c("FG", "UKBB"),
    tissue = "Skin_Not_Sun_Exposed_Suprapubic",
    symbol = "TP53",
    gtex_label = "GTEx: Skin -\nNot Sun Exposed Suprapubic (TP53)",
    tag = "b",
    window = 50000
  ),
  "1:16505320:A:G" = list(
    trait = "GGT",
    lead_variant = "1:16505320:A:G",
    rsid = "rs1497406",
    cohorts = c("BBJ", "UKBB"),
    tissue = "Liver",
    symbol = "EPHA2",
    gtex_label = "GTEx: Liver (EPHA2)",
    tag = "c",
    window = 100000
  ),
  "3:71771215:T:TG" = list(
    trait = "LOY",
    lead_variant = "3:71771215:T:TG",
    rsid = "rs34778241",
    cohorts = c("BBJ", "UKBB"),
    tissue = "Whole_Blood",
    symbol = "EIF4E3",
    gtex_label = "GTEx: Whole blood (EIF4E3)",
    tag = "d",
    window = 100000
  )
)

plts <-
  purrr::map(params, function(p) {
    trait <- p$trait
    trait_label <- ifelse(!is.null(p$trait_label), p$trait_label, trait)
    lead_variant <- p$lead_variant
    rsid <- p$rsid
    cohorts <- p$cohorts
    tissue <- p$tissue
    symbol <- p$symbol
    gtex_label <- p$gtex_label
    tag <- p$tag
    window <- p$window

    lead_pos <- locusviz::parse_variant(lead_variant)$position
    chromosome <- locusviz::parse_variant(lead_variant)$chromosome
    start <- lead_pos - window
    end <- lead_pos + window

    panels <- purrr::map(c(cohorts, "GTEx"), function(cohort) {
      if (cohort == "GTEx") {
        data <- rgsutil::read_gsfile(fm_insights$get_gtex_locus_txt_path(tissue, symbol, lead_variant))
      } else {
        data <-
          rgsutil::read_gsfile(fm_insights$get_locus_txt_path(cohort, trait, lead_variant))
      }
      data <- dplyr::mutate(data, locusviz::parse_variant(variant)) %>%
        locusviz::preprocess(
          lead_variant = lead_variant,
          beta_col = "beta_marginal",
          se_col = "se_marginal",
          cs_id_col = "susie.cs_id",
          r2_col = "lead_r2"
        )
      p_manhattan <- locusviz::plot_manhattan_panel(
        data,
        highlight_pos = lead_pos,
        xlim = c(start, end),
        plot.loglog_p = FALSE,
        nlog10p_threshold = 0,
        point.size = 0.5,
        point.size2 = 1.5,
        line.size = 0.2,
        title = ifelse(cohort == "GTEx",
          gtex_label,
          paste(cohort, trait, sep = ": ")
        ),
        ggtheme = locusviz::get_default_theme(
          fontsize = 6,
          hide.xtext = TRUE,
          hide.xtitle = TRUE,
          legend.position = if (cohort == cohorts[1]) {
            c(1, 1.1)
          } else {
            "none"
          }
        ),
        rasterize = TRUE
      ) + theme(
        plot.title = element_text(
          hjust = 0.01,
          margin = margin(b = -12),
          size = 6
        ),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.15, units = "cm")
      )

      if (cohort == cohorts[1]) {
        p_manhattan <- p_manhattan +
          geom_text(
            aes(x, y, label = label),
            data = tibble::tibble(
              x = lead_pos,
              y = max(data$nlog10p),
              label = rsid
            ),
            size = 2,
            hjust = -0.1
          ) + labs(tag = tag)
      }

      p_fm <- locusviz::plot_fm_panel(
        data,
        highlight_pos = lead_pos,
        xlim = c(start, end),
        ylim = c(0, 1.1),
        ybreaks = seq(0, 1, length = 2),
        point.size = 0.5,
        point.size2 = 1.5,
        ggtheme = locusviz::get_default_theme(
          fontsize = 6,
          hide.xtext = cohort != "GTEx",
          hide.xtitle = TRUE,
          legend.position = "none"
        ) + theme(plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm")),
        rasterize = TRUE
      )

      return(list(p_manhattan, p_fm))
    }) %>% purrr::flatten()

    plt <-
      append(panels, list(
        locusviz::plot_gene_panel(
          chromosome,
          start,
          end,
          highlight_pos = lead_pos,
          fontsize = 6,
          point.size = 1.5,
          label.size = 1.5,
          arrow.rate = 0.075,
          length = unit(0.05, "cm")
        ) + theme(plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm"))
      )) %>%
      purrr::reduce(`+`) + patchwork::plot_layout(ncol = 1, heights = c(rep(c(1, 0.3), 3), 0.05))
    return(plt)
  })

plt <- (plts[[1]] | plts[[2]]) / (plts[[3]] | plts[[4]])
plt

cowplot::save_plot(
  "figures/ExDataFig6_locuszoom_coloc.pdf",
  plt,
  base_height = 9.6,
  base_width = 7.2,
  device = cairo_pdf
)