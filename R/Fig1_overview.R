library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

#########################################################
# Fig. 1a: overview
#########################################################
trait <- "CAD"
lead_variant <- "6:12903957:A:G"
tissue <- "Artery_Tibial"
symbol <- "PHACTR1"
window <- 250000
lead_pos <- locusviz::parse_variant(lead_variant)$position
chromosome <- locusviz::parse_variant(lead_variant)$chromosome
start <- lead_pos - window
end <- lead_pos + window

panels <- purrr::map(c(cohorts, "GTEx"), function(cohort) {
  if (cohort == "GTEx") {
    data <- rgsutil::read_gsfile(fm_insights$get_gtex_locus_txt_path(tissue, symbol, lead_variant))
  } else {
    data <- rgsutil::read_gsfile(fm_insights$get_locus_txt_path(cohort, trait, lead_variant))
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
    title = ifelse(
      cohort == "GTEx",
      "GTEx: Artery - Tibial (PHACTR1)",
      paste(cohort, trait, sep = ": ")
    ),
    ggtheme = locusviz::get_default_theme(
      fontsize = 6,
      hide.xtext = TRUE,
      hide.xtitle = TRUE,
      legend.position = if (cohort == "BBJ") {
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

  if (cohort == "BBJ") {
    p_manhattan <- p_manhattan +
      geom_text(
        aes(x, y, label = label),
        data = tibble::tibble(
          x = lead_pos,
          y = max(data$nlog10p),
          label = "rs9349379"
        ),
        size = 2,
        hjust = -0.1
      )
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
    )
  )) %>%
  purrr::reduce(`+`) + patchwork::plot_layout(ncol = 1, heights = c(rep(c(1, 0.3), 4), 0.05))

cowplot::save_plot(
  "figures/Fig1/Fig1_locuszoom.pdf",
  plt,
  base_height = 3.6,
  base_width = 3.6,
  device = cairo_pdf
)

#########################################################
# Fig. 1b: trait upset
#########################################################

df.traits <- read_trait_summary()

matrix <-
  dplyr::group_by(df.traits, trait) %>%
  dplyr::summarize(cohorts = list(cohort)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cohorts_collapsed = factor(purrr::map_chr(cohorts, function(x) {
    paste0(sort(x), collapse = "-")
  }))) %>%
  dplyr::group_by(cohorts, cohorts_collapsed) %>%
  dplyr::summarize(count = n(), degree = length(cohorts[[1]])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cohorts_collapsed = factor(cohorts_collapsed, levels = cohorts_collapsed[order(degree, count, decreasing = TRUE)])) %>%
  tidyr::unnest(cohorts) %>%
  dplyr::mutate(cohorts = forcats::fct_rev(cohorts))

segment <-
  dplyr::mutate(matrix, cohorts = factor(cohorts, ordered = TRUE)) %>%
  dplyr::group_by(cohorts_collapsed) %>%
  dplyr::summarize(y = min(cohorts), yend = max(cohorts))

p_mat <-
  ggplot() +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill
    ),
    data = tibble::tibble(
      xmin = 0.5,
      xmax = Inf,
      ymin = -0.5 + seq(3),
      ymax = 0.5 + seq(3)
    ),
    fill = c("#F2F2F2", "white", "#F2F2F2")
  ) +
  scale_y_discrete() +
  scale_x_discrete() +
  geom_segment(aes(
    x = cohorts_collapsed,
    xend = cohorts_collapsed,
    y = y,
    yend = yend
  ),
  data = segment
  ) +
  geom_point(aes(cohorts_collapsed, cohorts, color = cohorts), data = matrix) +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 6.4),
    legend.position = "none",
    plot.margin = margin(0, 0.1, 0.1, 0.1, unit = "cm")
  ) +
  scale_color_manual(values = cohort_colors)

p_bar <-
  dplyr::distinct(matrix, cohorts_collapsed, count, degree) %>%
  ggplot(aes(cohorts_collapsed, count)) +
  geom_text(aes(label = count), vjust = -0.5, size = 2) +
  geom_col(aes(fill = factor(degree))) +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = BuenColors::jdb_palette("brewer_blue")[c(3, 5, 7)]) +
  scale_y_continuous(expand = expansion(c(0, 0.05), 0)) +
  labs(y = "No. traits") +
  coord_cartesian(clip = "off")

# p_domains =
#   dplyr::group_by(df.traits, domain) %>%
#   dplyr::summarize(count = length(unique(trait))) %>%
#   ggplot(aes(x = count, y = domain)) +
#   geom_col(aes(fill = domain)) +
#   geom_text(aes(
#     label = scales::comma(count, accuracy = 1),
#     hjust = -0.3
#   ), size = 2) +
#   my_theme +
#   theme(
#     axis.title.y = element_blank(),
#     legend.position = "none",
#     plot.margin = margin(0, 0.1, 0.1, 0.1, unit = "cm")
#   ) +
#   scale_x_continuous(expand = expansion(c(0, 0.2), 0)) +
#   scale_y_discrete(limits = rev(domains)) +
#   scale_fill_manual(values = domain_colors) +
#   labs(x = "No. traits") +
#   coord_cartesian(clip = "off")

plot_maf_dist <- function(df, xlab) {
  if (!("max_pip_color" %in% colnames(df))) {
    df <- dplyr::mutate(df, max_pip_coloc = max_pip)
  }
  df.label <-
    dplyr::group_by(df, cohort) %>%
    dplyr::summarize(
      count = sum(max_pip > 0.1, na.rm = T),
      coloc_count = sum(max_pip_coloc > 0.1, na.rm = T),
      frac_coloc = coloc_count / count,
      coloc_count_maf_gt_005 = sum(max_pip_coloc > 0.1 &
        maf > 0.05, na.rm = T),
      coloc_count_maf_le_005 = sum(max_pip_coloc > 0.1 &
        maf <= 0.05, na.rm = T),
      frac_coloc_maf_gt_005 = coloc_count_maf_gt_005 / coloc_count,
      frac_coloc_maf_le_005 = coloc_count_maf_le_005 / coloc_count
    )

  dplyr::filter(df, max_pip_coloc > 0.1) %>%
    ggplot(aes(maf, cohort)) +
    geom_vline(
      xintercept = 0.05,
      linetype = "dashed",
      color = "grey50"
    ) +
    ggridges::geom_density_ridges(
      aes(fill = cohort),
      stat = "binline",
      bins = 50,
      scale = 0.95,
      draw_baseline = FALSE,
      size = 0.2
    ) +
    geom_text(
      aes(
        x = 0.1,
        y = cohort,
        label = scales::percent(frac_coloc_maf_gt_005, accuracy = 1)
      ),
      vjust = -3,
      size = 2,
      data = df.label
    ) +
    geom_text(
      aes(
        x = 0.03,
        y = cohort,
        label = scales::percent(frac_coloc_maf_le_005, accuracy = 1)
      ),
      vjust = -3,
      size = 2,
      data = df.label
    ) +
    my_theme +
    theme(
      legend.position = "none",
      axis.title.y = element_blank()
    ) +
    scale_y_discrete(limits = rev(cohorts), expand = expansion(0, 0)) +
    scale_fill_manual(values = cohort_colors) +
    scale_x_log10(breaks = c(0.01, 0.05, 0.5)) +
    labs(x = xlab) +
    coord_cartesian(xlim = c(0.001, 0.5))
}



df.coloc_maf <- rgsutil::read_gsfile(
  fm_insights$get_merged_results_path("fm_only.coloc.max_pip_per_vc.max_pip01.af", "tsv.bgz")
) %>%
  dplyr::mutate(maf = 0.5 - abs(0.5 - af), max_pip_coloc = max_pip) %>%
  tidyr::drop_na(maf)

df.coding_maf <- rgsutil::read_gsfile(
  fm_insights$get_merged_results_path("fm_only.max_pip_per_cohort.consequence.coding.max_pip01.af", "tsv.bgz")
) %>%
  dplyr::mutate(maf = 0.5 - abs(0.5 - af), max_pip_coloc = max_pip) %>%
  tidyr::drop_na(maf)

p_coloc_maf <- plot_maf_dist(df.coloc_maf, xlab = "MAF of colocalized variants")
p_coding_maf <- plot_maf_dist(df.coding_maf, xlab = "MAF of coding variants")
p_coloc_maf / p_coding_maf

# df.eqtl_samples = read.table("metadata/eqtl_catalogue_r4_population_assignments.tsv",
#                              T,
#                              sep = "\t") %>%
#   tidyr::pivot_longer(cols = c("EUR", "AFR", "SAS", "EAS", "Unassigned")) %>%
#   dplyr::group_by(name) %>%
#   dplyr::summarize(count = sum(value)) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(
#     prop = count / sum(count),
#     name = forcats::fct_other(name, drop = "Unassigned"),
#     name = forcats::fct_reorder2(
#       name,
#       name,
#       count,
#       .fun = function(x, y) {
#         dplyr::case_when(x == "Other" ~ -Inf,
#                          TRUE ~ as.double(y))
#       },
#       .desc = F
#     )
#   )
#
# p_eqtl_samples =
#   ggplot(df.eqtl_samples, aes(count, name)) +
#   geom_col(aes(fill = name))  +
#   geom_text(aes(label = scales::comma(count), hjust = -0.3), size = 2) +
#   my_theme +
#   theme(
#     axis.title.y = element_blank(),
#     legend.position = "none",
#     plot.margin = margin(0, 0.1, 0.1, 0.1, unit = "cm")
#   ) +
#   scale_x_continuous(
#     expand = expansion(c(0, 0.2), 0),
#     breaks = c(0, 2500, 5000),
#     labels = scales::comma
#   ) +
#   scale_fill_manual(values = c(locusviz::get_gnomad_colors(), "Other" = "grey50")) +
#   labs(x = "No. eQTL donors") +
#   coord_cartesian(clip = "off")

# df.eqtl_tissues = read.table("metadata/eqtl_tissues.txt", T, sep = "\t") %>%
#   dplyr::mutate(study = stringr::str_replace(study, "Schmiedel_2018", "Schmiedel")) %>%
#   dplyr::distinct(study, tissue_label) %>%
#   dplyr::bind_rows(tibble::tibble(study = "GTEx", tissue_label = rep(NA, 49))) %>%
#   dplyr::mutate(study = forcats::fct_lump_n(study, 3),
#                 study = forcats::fct_inorder(study)) %>%
#   dplyr::group_by(study) %>%
#   dplyr::summarize(count = n()) %>%
#   dplyr::ungroup()
#
# p_eqtl_tissue =
#   ggplot(df.eqtl_tissues, aes(count, study)) +
#   geom_col(aes(fill = study))  +
#   geom_text(aes(
#     label = scales::comma(count, accuracy = 1),
#     hjust = -0.3
#   ), size = 2) +
#   my_theme +
#   theme(axis.title.y = element_blank(), legend.position = "none") +
#   scale_x_continuous(expand = expansion(c(0, 0.1), 0)) +
#   scale_fill_manual(values = c("grey50", BuenColors::jdb_palette("Darjeeling")[c(5, 4, 2)])) +
#   labs(x = "No. cell-types/tissues") +
#   coord_cartesian(clip = "off")


plt <-
  purrr::reduce(
    list(
      p_bar + labs(tag = "b") + theme(axis.title.y = element_text(margin = margin(0, -12, 0, 0, unit = "pt"))),
      p_mat,
      p_coloc_maf + labs(tag = "g"),
      p_coding_maf + labs(tag = "h") + theme(plot.margin = margin(0.1, 0.1, 0, 0.1, unit = "cm"))
    ),
    `+`
  ) + patchwork::plot_layout(ncol = 1, heights = c(1, 0.4, 1, 0.8))

cowplot::save_plot(
  "figures/Fig1/Fig1_barplots.pdf",
  plt,
  base_height = 3.8,
  base_width = 1.8,
  device = cairo_pdf
)


#########################################################
# Fig. 1: stats bar plots
#########################################################

df.in_cs <-
  rgsutil::read_gsfile(fm_insights$get_merged_results_path("fm_only.in_cs.consequence", "tsv.bgz"))

fct_bin_pip <- function(pip) {
  pip_bin <- cut(pip, pip_bin_breaks2)
  pip_bin <- forcats::fct_relabel(pip_bin, function(x) {
    x <- stringr::str_remove_all(x, "\\(|\\]")
    interval <- stringr::str_split_fixed(x, ",", 2)
    dplyr::case_when(
      x == "-Inf,0.01" ~ "≤1%",
      x == "0.9,1" ~ ">90%",
      TRUE ~ as.character(
        stringr::str_glue(
          "{as.numeric(interval[,1]) * 100}–{as.numeric(interval[,2]) * 100}%"
        )
      )
    )
  })
  return(pip_bin)
}

theme_stats_bar <- list(
  scale_x_continuous(expand = c(0, 0)),
  theme_void(),
  theme(
    legend.position = "none",
    plot.margin = margin(0.1, 0.2, 0.1, 0.2, unit = "cm"),
    text = element_text(size = 8),
    plot.tag = element_text(
      size = 9.6,
      face = "bold",
      margin = margin(0, 0, 0, -12, unit = "pt")
    ),
    plot.tag.position = c(0, 1),
    plot.title = element_text(size = 8, margin = margin(0, 0, -4, 0, unit = "pt"))
  ),
  coord_cartesian(clip = "off")
)


plot_stats_bars <- function(df,
                            cohort,
                            tag = NULL,
                            text.size = 2,
                            width = 0.5,
                            scale = 0.1,
                            tip_length = 0) {
  df <- dplyr::filter(df, cohort == .env$cohort)
  df.filtered <- dplyr::filter(df, susie.cs_id > 0)

  df.n_cs_per_region <-
    dplyr::group_by(df.filtered, cohort, trait, region) %>%
    dplyr::summarize(n = length(unique(susie.cs_id))) %>%
    dplyr::mutate(
      n = cut(n, c(seq(0, 4), Inf)),
      n = forcats::fct_relabel(n, function(x) {
        stringr::str_remove_all(x, "\\([0-9],|\\]") %>%
          stringr::str_replace("Inf", "≥5")
      })
    ) %>%
    dplyr::group_by(cohort, n) %>%
    dplyr::summarize(count = n())

  df.n_var_per_cs <-
    dplyr::group_by(df.filtered, cohort, trait, region, susie.cs_id) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      n = cut(n, c(0, 1, 5, 10, 20, 50, Inf)),
      n = forcats::fct_relabel(n, function(x) {
        x <- stringr::str_remove_all(x, "\\(|\\]")
        interval <- stringr::str_split_fixed(x, ",", 2)
        dplyr::case_when(
          x == "0,1" ~ "1",
          x == "50,Inf" ~ ">50",
          TRUE ~ as.character(
            stringr::str_glue("{as.numeric(interval[,1]) + 1}–{interval[,2]}")
          )
        )
      })
    ) %>%
    dplyr::group_by(cohort, n) %>%
    dplyr::summarize(count = n())

  df.max_pip_bin_per_cs <-
    dplyr::group_by(df.filtered, cohort, trait, region, susie.cs_id) %>%
    dplyr::summarize(max_pip = max(pip, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      max_pip = ifelse(is.infinite(max_pip), 0, max_pip),
      max_pip_bin = fct_bin_pip(max_pip)
    ) %>%
    dplyr::group_by(cohort, max_pip_bin) %>%
    dplyr::summarize(count = n())

  df.max_pip_bin_per_var <-
    dplyr::filter(df, pip > 0.1) %>%
    dplyr::group_by(cohort, trait) %>%
    dplyr::mutate(pip_bin = fct_bin_pip(pip)) %>%
    tidyr::drop_na(pip_bin) %>%
    dplyr::group_by(cohort, pip_bin) %>%
    dplyr::summarize(count = n())

  # numbers
  n_cs <- dplyr::distinct(df.filtered, cohort, trait, region, susie.cs_id) %>%
    nrow()
  n_region <- dplyr::distinct(df.filtered, cohort, trait, region) %>%
    nrow()
  n_trait <- dplyr::distinct(df.filtered, cohort, trait) %>%
    nrow()

  df.cs_le_50 <- dplyr::group_by(df.filtered, cohort, trait, region, susie.cs_id) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n <= 50)
  n_cs_le_50 <- dplyr::distinct(df.cs_le_50, cohort, trait, region, susie.cs_id) %>%
    nrow()
  n_trait_cs_le_50 <- dplyr::distinct(df.cs_le_50, cohort, trait) %>%
    nrow()

  df.cs_pip_gt_10 <-
    dplyr::group_by(df.filtered, cohort, trait, region, susie.cs_id) %>%
    dplyr::summarize(max_pip = max(pip, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(max_pip > 0.1)

  n_cs_pip_gt_10 <- dplyr::distinct(df.cs_pip_gt_10, cohort, trait, region, susie.cs_id) %>%
    nrow()
  n_trait_cs_pip_gt_10 <- dplyr::distinct(df.cs_pip_gt_10, cohort, trait) %>%
    nrow()

  df.var_pip_gt_10 <-
    dplyr::filter(df, pip > 0.1)
  n_var_pip_gt_10 <- nrow(df.var_pip_gt_10)
  n_trait_var_pip_gt_10 <- dplyr::distinct(df.var_pip_gt_10, cohort, trait) %>%
    nrow()



  p1 <-
    ggplot(df.n_cs_per_region, aes(count, cohort)) +
    geom_col(aes(fill = n),
      width = width,
      position = position_stack(reverse = TRUE)
    ) +
    geom_text(
      aes(
        label = scales::comma(count, accuracy = 1),
        color = n
      ),
      size = text.size,
      position = position_stack(0.5, reverse = TRUE)
    ) +
    geom_text(aes(y = 1.5, label = n),
      size = text.size,
      position = position_stack(0.5, reverse = TRUE)
    ) +
    locusviz::annotate_hrange(
      0,
      n_region,
      0.5,
      scale = 0.15,
      label = stringr::str_glue(
        "{scales::comma(n_cs)} independent 95% CS in {scales::comma(n_region)} regions of {n_trait} traits"
      ),
      tip_length = tip_length,
      text.size = text.size
    ) +
    scale_fill_manual(values = shades::opacity(cohort_colors[cohort], seq(0.1, 1, length.out = 5))) +
    scale_color_manual(values = c(rep("black", 4), "white")) +
    theme_stats_bar +
    labs(title = cohort, tag = tag)


  p2 <- ggplot(df.n_var_per_cs, aes(count, cohort)) +
    geom_col(aes(fill = n),
      width = width,
      position = position_stack(reverse = TRUE)
    ) +
    geom_text(
      aes(
        label = scales::comma(count, accuracy = 1),
        color = n
      ),
      size = text.size,
      position = position_stack(0.5, reverse = TRUE)
    ) +
    geom_text(aes(y = 1.5, label = n),
      size = text.size,
      position = position_stack(0.5, reverse = TRUE)
    ) +
    locusviz::annotate_hrange(
      0,
      n_cs_le_50,
      0.5,
      scale = 0.1,
      label = stringr::str_glue(
        "{scales::comma(n_cs_le_50)} CS of {n_trait_cs_le_50} traits fine-mapped to ≤50 variants"
      ),
      tip_length = tip_length,
      text.size = text.size
    ) +
    scale_fill_manual(values = shades::opacity(cohort_colors[cohort], seq(1, 0.1, length.out = 6))) +
    scale_color_manual(values = c("white", rep("black", 5))) +
    theme_stats_bar

  p3 <-
    ggplot(df.max_pip_bin_per_cs, aes(-count, cohort)) +
    geom_col(aes(fill = max_pip_bin),
      width = width,
      position = position_stack(0.5, reverse = T)
    ) +
    geom_text(
      aes(
        label = scales::comma(count, accuracy = 1),
        color = max_pip_bin
      ),
      size = text.size,
      position = position_stack(0.5, reverse = T)
    ) +
    geom_text(
      aes(y = 1.5, label = max_pip_bin),
      size = text.size,
      position = position_stack(0.5, reverse = TRUE)
    ) +
    locusviz::annotate_hrange(
      -n_cs, -(n_cs - n_cs_pip_gt_10),
      0.5,
      scale = 0.025,
      label = stringr::str_glue(
        "{scales::comma(n_cs_pip_gt_10)} CS of {n_trait_cs_pip_gt_10} traits have a variant with PIP >10%"
      ),
      tip_length = tip_length,
      text.size = text.size
    ) +
    scale_fill_manual(values = shades::opacity(cohort_colors[cohort], seq(0.1, 1, length.out = 5))) +
    scale_color_manual(values = c(rep("black", 4), "white")) +
    theme_stats_bar

  p4 <-
    ggplot(df.max_pip_bin_per_var, aes(-count, cohort)) +
    geom_col(aes(fill = pip_bin),
      width = width,
      position = position_stack(0.5, reverse = T)
    ) +
    geom_text(
      aes(
        label = scales::comma(count, accuracy = 1),
        color = pip_bin
      ),
      size = text.size,
      position = position_stack(0.5, reverse = T)
    ) +
    geom_text(
      aes(y = 1.5, label = pip_bin),
      size = text.size,
      position = position_stack(0.5, reverse = TRUE)
    ) +
    locusviz::annotate_hrange(
      -n_var_pip_gt_10,
      0,
      0.5,
      scale = 0.25,
      label = stringr::str_glue(
        "{scales::comma(n_var_pip_gt_10)} variant-trait pairs have PIP >10%"
      ),
      tip_length = tip_length,
      text.size = text.size
    ) +
    scale_fill_manual(values = shades::opacity(cohort_colors[cohort], seq(0.1, 1, length.out = 5)[3:5])) +
    scale_color_manual(values = c(rep("black", 2), "white")) +
    theme_stats_bar

  purrr::reduce(list(p1, p2, p3, p4), `+`) + patchwork::plot_layout(ncol = 1)
}

plot_stats_bars_coloc <- function(tag = NULL,
                                  text.size = 2,
                                  width = 0.5,
                                  scale = 0.1,
                                  tip_length = 0,
                                  color = "#706060",
                                  color2 = "#606070") {
  df.max_pip_per_gene <- rgsutil::read_gsfile(fm_insights$get_merged_eqtl_path("merged.max_pip_per_gene", "tsv.bgz")) %>%
    dplyr::filter(max_pip > 0.01) %>%
    dplyr::mutate(max_pip_bin = fct_bin_pip(max_pip)) %>%
    dplyr::group_by(max_pip_bin) %>%
    dplyr::summarize(count = n(), cohort = "eQTL")

  df.max_pip_per_var <- rgsutil::read_gsfile(fm_insights$get_merged_eqtl_path("merged.max_pip_per_var", "tsv.bgz")) %>%
    dplyr::filter(max_pip > 0.1) %>%
    dplyr::mutate(max_pip_bin = fct_bin_pip(max_pip)) %>%
    dplyr::group_by(max_pip_bin) %>%
    dplyr::summarize(count = n(), cohort = "eQTL")

  df.max_pip_coloc_per_gt <- rgsutil::read_gsfile(
    fm_insights$get_merged_results_path("fm_only.coloc.max_pip_per_gt", "tsv.bgz")
  ) %>%
    dplyr::filter(max_pip_coloc > 0.1) %>%
    dplyr::mutate(max_pip_bin = fct_bin_pip(max_pip_coloc)) %>%
    dplyr::group_by(max_pip_bin) %>%
    dplyr::summarize(count = n(), cohort = "eQTL")

  df.max_pip_coloc_per_vtg <- rgsutil::read_gsfile(
    fm_insights$get_merged_results_path("fm_only.coloc.max_pip_per_vtg", "tsv.bgz")
  ) %>%
    dplyr::filter(max_pip_coloc > 0.1) %>%
    dplyr::mutate(max_pip_bin = fct_bin_pip(max_pip_coloc)) %>%
    dplyr::group_by(max_pip_bin) %>%
    dplyr::summarize(count = n(), cohort = "eQTL")


  n_gene <- dplyr::pull(df.max_pip_per_gene, count) %>% sum()
  n_gene_gt_10 <- dplyr::filter(
    df.max_pip_per_gene,
    max_pip_bin %in% c("10–50%", "50–90%", ">90%")
  ) %>%
    dplyr::pull(count) %>%
    sum()
  n_var_pip_gt_10 <- dplyr::pull(df.max_pip_per_var, count) %>% sum()

  n_gene_pip_coloc <- dplyr::pull(df.max_pip_coloc_per_gt, count) %>% sum()
  n_gene_pip_coloc_gt_10 <- dplyr::filter(
    df.max_pip_coloc_per_gt,
    max_pip_bin %in% c("10–50%", "50–90%", ">90%")
  ) %>%
    dplyr::pull(count) %>%
    sum()
  n_var_pip_coloc_gt_10 <- dplyr::pull(df.max_pip_coloc_per_vtg, count) %>% sum()

  p1 <-
    ggplot(df.max_pip_per_gene, aes(-count, cohort)) +
    geom_col(aes(fill = max_pip_bin),
      width = width,
      position = position_stack(reverse = TRUE)
    ) +
    geom_text(
      aes(
        label = scales::comma(count, accuracy = 1),
        color = max_pip_bin
      ),
      size = text.size,
      position = position_stack(0.5, reverse = TRUE)
    ) +
    geom_text(
      aes(y = 1.5, label = max_pip_bin),
      size = text.size,
      position = position_stack(0.5, reverse = TRUE)
    ) +
    locusviz::annotate_hrange(
      -n_gene, -(n_gene - n_gene_gt_10),
      0.5,
      scale = 0.1,
      label = stringr::str_glue(
        "{scales::comma(n_gene_gt_10)} genes have an eQTL variant with PIP >10% in ≥1 tissue(s)"
      ),
      tip_length = tip_length,
      text.size = text.size
    ) +
    scale_fill_manual(values = shades::opacity(color, seq(0.1, 1, length.out = 5)[2:5])) +
    scale_color_manual(values = c(rep("black", 3), "white")) +
    theme_stats_bar +
    labs(title = "Colocalization", tag = tag)


  p2 <- ggplot(df.max_pip_per_var, aes(-count, cohort)) +
    geom_col(aes(fill = max_pip_bin),
      width = width,
      position = position_stack(0.5, reverse = T)
    ) +
    geom_text(
      aes(
        label = scales::comma(count, accuracy = 1),
        color = max_pip_bin
      ),
      size = text.size,
      position = position_stack(0.5, reverse = T)
    ) +
    geom_text(
      aes(y = 1.5, label = max_pip_bin),
      size = text.size,
      position = position_stack(0.5, reverse = TRUE)
    ) +
    locusviz::annotate_hrange(
      -n_var_pip_gt_10,
      0,
      0.5,
      scale = 0.2,
      label = stringr::str_glue(
        "{scales::comma(n_var_pip_gt_10)} variants have PIP >10% in ≥1 tissue(s)"
      ),
      tip_length = tip_length,
      text.size = text.size
    ) +
    scale_fill_manual(values = shades::opacity(color, seq(0.1, 1, length.out = 5)[3:5])) +
    scale_color_manual(values = c(rep("black", 2), "white")) +
    theme_stats_bar

  p3 <-
    ggplot(df.max_pip_coloc_per_gt, aes(-count, cohort)) +
    geom_col(aes(fill = max_pip_bin),
      width = width,
      position = position_stack(0.5, reverse = T)
    ) +
    geom_text(
      aes(
        label = scales::comma(count, accuracy = 1),
        color = max_pip_bin
      ),
      size = text.size,
      position = position_stack(0.5, reverse = T)
    ) +
    geom_text(
      aes(y = 1.5, label = max_pip_bin),
      size = text.size,
      position = position_stack(0.5, reverse = TRUE)
    ) +
    locusviz::annotate_hrange(
      -n_gene_pip_coloc, -(n_gene_pip_coloc - n_gene_pip_coloc_gt_10),
      0.5,
      scale = 0.05,
      label = stringr::str_glue(
        "{scales::comma(n_gene_pip_coloc_gt_10)} gene-trait pairs have an eQTL variant with coloc PIP >10% in ≥1 tissue(s)"
      ),
      tip_length = tip_length,
      text.size = text.size
    ) +
    scale_fill_manual(values = shades::opacity(color2, seq(0.1, 1, length.out = 5)[3:5])) +
    scale_color_manual(values = c(rep("black", 2), "white")) +
    theme_stats_bar

  p4 <- ggplot(df.max_pip_coloc_per_vtg, aes(-count, cohort)) +
    geom_col(aes(fill = max_pip_bin),
      width = width,
      position = position_stack(0.5, reverse = T)
    ) +
    geom_text(
      aes(
        label = scales::comma(count, accuracy = 1),
        color = max_pip_bin
      ),
      size = text.size,
      position = position_stack(0.5, reverse = T)
    ) +
    geom_text(
      aes(y = 1.5, label = max_pip_bin),
      size = text.size,
      position = position_stack(0.5, reverse = TRUE)
    ) +
    locusviz::annotate_hrange(
      -n_var_pip_coloc_gt_10,
      0,
      0.5,
      scale = 0.15,
      label = stringr::str_glue(
        "{scales::comma(n_var_pip_coloc_gt_10)} variant-trait-gene triples have coloc PIP >10% in ≥1 tissue(s)"
      ),
      tip_length = tip_length,
      text.size = text.size
    ) +
    scale_fill_manual(values = shades::opacity(color2, seq(0.1, 1, length.out = 5)[3:5])) +
    scale_color_manual(values = c(rep("black", 2), "white")) +
    theme_stats_bar


  purrr::reduce(list(p1, p2, p3, p4), `+`) + patchwork::plot_layout(ncol = 1)
}


plt <- (plot_stats_bars(df.in_cs, "BBJ", tag = "c") |
  plot_stats_bars(df.in_cs, "FG", tag = "d")) /
  (plot_stats_bars(df.in_cs, "UKBB", tag = "e") |
    plot_stats_bars_coloc(tag = "f"))
plt

cowplot::save_plot(
  "figures/Fig1/Fig1_stats_bar.pdf",
  plt,
  base_height = 3.6,
  base_width = 7.2,
  device = cairo_pdf
)