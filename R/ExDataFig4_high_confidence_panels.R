library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

####################################################

df.traits <- read_trait_summary()
df.pip09.every <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("shard_trait.pip09.every", "tsv"))
df.pip01.every <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("shard_trait.pip09.any.pip01.every", "tsv"))

df.pip01.every.only <- dplyr::filter(df.pip01.every, !(
  stringr::str_c(trait, variant, sep = "-") %in% stringr::str_c(df.pip09.every$trait, df.pip09.every$variant, sep = "-")
))

length(setdiff(
  unique(df.pip01.every$variant),
  unique(df.pip09.every$variant)
))

length(unique(df.pip01.every$variant))
dplyr::filter(
  df.pip01.every,
  consequence %in% c("pLoF", "Missense", "Synonymous")
) %>%
  dplyr::summarize(n_pairs = n(), n_variants = length(unique(variant)))

dplyr::filter(df.pip01.every, !(consequence %in% c("pLoF", "Missense", "Synonymous"))) %>%
  dplyr::summarize(
    n_pairs = n(),
    n_variants = length(unique(variant))
  )
dplyr::filter(df.pip01.every, !(consequence %in% c("pLoF", "Missense", "Synonymous"))) %>%
  dplyr::distinct(variant, coloc_variant = !is.na(max_pip_coloc)) %>%
  dplyr::filter(coloc_variant) %>%
  nrow()


dplyr::group_by(df.pip01.every, variant, trait) %>%
  dplyr::filter(sum(c(pip.BBJ, pip.FG, pip.UKBB) > 0.1, na.rm = T) == n_pop) %>%
  dplyr::filter(sum(c(pvalue.BBJ, pvalue.FG, pvalue.UKBB) > 5e-8, na.rm = T) > 0) %>%
  nrow()

dplyr::group_by(df.pip01.every, variant, trait) %>%
  dplyr::filter(sum(c(pip.BBJ, pip.FG, pip.UKBB) > 0.1, na.rm = T) == n_pop) %>%
  dplyr::filter(sum(c(pvalue.BBJ, pvalue.FG, pvalue.UKBB) > 1e-5, na.rm = T) > 0) %>%
  View()

plot_maf_across_cohorts <- function(df) {
  plt <-
    dplyr::select(df, variant, starts_with("af.")) %>%
    tidyr::pivot_longer(starts_with("af.")) %>%
    dplyr::mutate(
      cohort = stringr::str_split_fixed(name, "\\.", 2)[, 2],
      value = 0.5 - abs(0.5 - value)
    ) %>%
    tidyr::drop_na() %>%
    ggplot(aes(cohort, value)) +
    gghalves::geom_half_violin(
      aes(fill = cohort),
      position = position_nudge(x = -0.1),
      draw_quantiles = c(0.25, 0.5, 0.75),
      size = 0.3
    ) +
    gghalves::geom_half_violin(
      aes(fill = cohort),
      position = position_nudge(x = 0.1),
      draw_quantiles = c(0.25, 0.5, 0.75),
      side = "r",
      size = 0.3
    ) +
    geom_line(aes(group = variant), alpha = 0.2, size = 0.2) +
    geom_point(size = 0.1, color = "grey50") +
    my_theme +
    theme(
      axis.title.y = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.1),
      plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm")
    ) +
    scale_fill_manual(values = cohort_colors, guide = FALSE) +
    scale_x_discrete(limits = rev) +
    scale_y_log10(
      labels = scales::number_format(accuracy = 1e-6, drop0trailing = TRUE),
      breaks = c(0.001, 0.01, 0.05, 0.1, 0.5),
      limits = c(0.0004, 0.5)
    ) +
    coord_flip() +
    labs(y = "MAF")
  return(plt)
}

plot_consequence_barplot <- function(df) {
  data <- dplyr::distinct(df, variant, consequence, most_severe) %>%
    dplyr::mutate(consequence = factor(consequence, levels = annot_levels), ) %>%
    dplyr::group_by(consequence) %>%
    dplyr::summarize(count = n())

  ylim <- intersect(rev(annot_levels), unique(data$consequence))

  plt <-
    ggplot(data, aes(count, consequence)) +
    geom_col(aes(fill = consequence), position = position_stack()) +
    geom_text(aes(label = count),
      hjust = -0.5,
      size = 2
    ) +
    my_theme +
    theme(
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm")
    ) +
    scale_y_discrete(limits = ylim) +
    scale_fill_manual(values = annot_colors) +
    scale_x_continuous(expand = expansion(c(0, 0.15), 0)) +
    labs(x = "No. variants")
  return(plt)
}

# annotate TSS
txdb <- locusviz::load_txdb("hg19")
df.dist <-
  dplyr::filter(df.pip01.every, !(consequence %in% c("pLoF", "Missense", "Synonymous"))) %>%
  dplyr::distinct(variant, consequence) %>%
  purrr::pmap_dfr(function(variant, consequence) {
    lead_variant <- variant[1]
    print(lead_variant)
    window <- 5000000
    lead_pos <- locusviz::parse_variant(lead_variant)$position
    chromosome <- locusviz::parse_variant(lead_variant)$chromosome
    start <- lead_pos - window
    end <- lead_pos + window

    locusviz::compute_distance_to_gene(
      txdb,
      paste0("chr", chromosome),
      start,
      end,
      ref_position = lead_pos,
      type = "GB"
    ) %>%
      dplyr::mutate(
        variant = lead_variant,
        consequence = consequence[1]
      )
  })

plot_distance_to_nearest_gene <- function(df.dist) {
  dplyr::filter(df.dist) %>%
    dplyr::group_by(variant) %>%
    dplyr::summarize(score = min(score), consequence = consequence[1]) %>%
    dplyr::mutate(consequence = factor(consequence, levels = annot_levels)) %>%
    ggplot(aes(score / 1000, fill = consequence)) +
    geom_histogram(bins = 20) +
    my_theme +
    theme(
      legend.title = element_blank(),
      plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm")
    ) +
    scale_fill_manual(values = annot_colors[setdiff(annot_levels, c("pLoF", "Missense", "Synonymous"))], drop = TRUE) +
    scale_x_continuous(breaks = c(0, 500, 1000, 1500)) +
    labs(x = "Distance to the closest gene (kb)", y = "No. variants")
}

df.ems <- rgsutil::read_gsfile(fm_insights$get_analysis_path("max_ems_normalized", "tsv.bgz")) %>%
  dplyr::filter(!(consequence %in% c("pLoF", "Missense", "Synonymous"))) %>%
  dplyr::mutate(
    high_confidence = variant %in% df.pip01.every$variant,
    max_pip_bin = factor(
      dplyr::case_when(
        high_confidence ~ "High-confidence",
        TRUE ~ as.character(cut(
          max_pip, pip_bin_breaks2
        ))
      ),
      levels = c(
        "(-Inf,0.01]",
        "(0.01,0.1]",
        "(0.1,0.5]",
        "(0.5,0.9]",
        "(0.9,1]",
        "High-confidence"
      ),
      labels = c(
        "[0,0.01]",
        "(0.01,0.1]",
        "(0.1,0.5]",
        "(0.5,0.9]",
        "(0.9,1]",
        "High-confidence"
      )
    )
  )

plot_ems_distribution <- function(df.ems) {
  ggplot(df.ems, aes(max_ems_normalized, max_pip_bin)) +
    ggridges::stat_density_ridges(
      aes(fill = 0.5 - abs(0.5 - stat(ecdf))),
      size = 0.2,
      geom = "density_ridges_gradient",
      calc_ecdf = TRUE
    ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "grey50"
    ) +
    scale_fill_viridis_c(name = "Tail\nprobability", direction = -1) +
    scale_x_log10(
      breaks = c(0.1, 1, 10, 100, 1000),
      label = scales::comma_format(accuracy = 0.0001, drop0trailing = TRUE)
    ) +
    my_theme +
    theme(legend.position = "right") +
    labs(x = "Max. normalized EMS for whole blood", y = "Best PIP bin")
}

plot_ems_high_confidence_enrichment <- function(df.ems) {
  dplyr::filter(df.ems, max_pip_bin == "High-confidence") %>%
    ggplot(aes(consequence, max_ems_normalized, color = consequence)) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "grey50"
    ) +
    stat_summary(fun.data = mean_cl_boot, fun.args = list(conf.int = 0.95)) +
    my_theme +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.tag.position = c(-0.025, 1)
    ) +
    scale_color_manual(values = annot_colors) +
    scale_x_discrete(limits = c("Non-genic", "CRE", "UTR5", "UTR3", "Promoter")) +
    scale_y_log10() +
    labs(y = "Max. normalized EMS for whole blood") +
    coord_cartesian(ylim = c(0.5, 100))
}

plt1 <- plot_maf_across_cohorts(df.pip01.every) + labs(tag = "a")
plt2 <- plot_consequence_barplot(df.pip01.every) + labs(tag = "b")
plt3 <- plot_distance_to_nearest_gene(df.dist) + labs(tag = "c")
plt4 <- plot_ems_distribution(df.ems) + labs(tag = "d")
# plt5 = plot_ems_high_confidence_enrichment(df.ems) + labs(tag = "e")


# plt_combined = (plt1 + plt2 + plt3 + patchwork::plot_layout(nrow = 1)) / (plt4 + plt5 + patchwork::plot_layout(nrow = 1, widths = c(2, 1)))
plt_combined <- (plt1 + plt2 + plt3 + patchwork::plot_layout(nrow = 1)) / (plt4 + patchwork::plot_spacer() + patchwork::plot_layout(nrow = 1, widths = c(2.4, 0.6)))
plt_combined

cowplot::save_plot(
  "figures/ExDataFig4_high_confidence_panels.pdf",
  plt_combined,
  base_height = 4.8,
  base_width = 7.2,
  device = cairo_pdf
)