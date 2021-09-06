library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

###############################
# ALP ASGR1 locuszoom
###############################

trait <- "ALP"
lead_variant <- "17:7080316:C:T"
rsid <- "rs55714927"
window <- 150000
lead_pos <- locusviz::parse_variant(lead_variant)$position
chromosome <- locusviz::parse_variant(lead_variant)$chromosome
start <- lead_pos - window
end <- lead_pos + window

p1 <-
  purrr::map(c("BBJ", "UKBB"), function(cohort) {
    data <-
      rgsutil::read_gsfile(fm_insights$get_locus_txt_path(cohort, trait, lead_variant)) %>%
      dplyr::mutate(locusviz::parse_variant(variant)) %>%
      locusviz::preprocess(
        lead_variant = lead_variant,
        beta_col = "beta_marginal",
        se_col = "se_marginal",
        cs_id_col = "susie.cs_id",
        r2_col = "lead_r2"
      )

    p_manhattan <- locusviz::plot_manhattan_panel(
      data,
      highlight_pos = c(lead_pos, 7069412),
      xlim = c(start, end),
      plot.loglog_p = FALSE,
      nlog10p_threshold = 0,
      point.size = 0.5,
      point.size2 = 1.5,
      line.size = 0.2,
      title = paste(cohort, trait, sep = ": "),
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
            label = rsid
          ),
          size = 2,
          hjust = -0.1
        ) + labs(tag = "a") + theme(plot.tag = element_text(size = 9.6))
    } else if (cohort == "UKBB") {
      p_manhattan <- p_manhattan +
        geom_text(
          aes(x, y, label = label),
          data = tibble::tibble(
            x = lead_pos,
            y = data$nlog10p[which(data$variant == lead_variant)],
            label = rsid
          ),
          size = 2,
          hjust = -0.1
        ) +
        geom_text(
          aes(x, y, label = label),
          data = tibble::tibble(
            x = 7069412,
            y = max(data$nlog10p),
            label = "rs186021206\n(r2 = 0.86 to del12)"
          ),
          size = 2,
          hjust = 1.1,
          vjust = 0.89
        ) +
        geom_point(
          aes(x, y),
          shape = 15,
          size = 1.5,
          color = "grey50",
          data = tibble::tibble(
            x = 7080256,
            y = max(data$nlog10p)
          )
        ) +
        geom_text(
          aes(x, y, label = label),
          data = tibble::tibble(
            x = 7080256,
            y = max(data$nlog10p),
            label = "del12 (non-imputed)"
          ),
          size = 2,
          hjust = -0.08
        )
    }

    p_fm <- locusviz::plot_fm_panel(
      data,
      highlight_pos = c(lead_pos, 7069412),
      xlim = c(start, end),
      ylim = c(0, 1.1),
      ybreaks = seq(0, 1, length = 2),
      point.size = 0.5,
      point.size2 = 1.5,
      ggtheme = locusviz::get_default_theme(
        fontsize = 6,
        hide.xtext = cohort != "UKBB",
        hide.xtitle = TRUE,
        legend.position = "none"
      ) + theme(plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm")),
      rasterize = TRUE
    )

    return(list(p_manhattan, p_fm))
  }) %>%
  purrr::flatten() %>%
  append(list(
    locusviz::plot_gene_panel(
      chromosome,
      start,
      end,
      highlight_pos = lead_pos,
      fontsize = 6,
      point.size = 1.5,
      label.size = 1.5,
      length = unit(0.05, "cm")
    )
  )) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1, heights = c(rep(c(1, 0.3), 2), 0.2))

cowplot::save_plot(
  "figures/ExDataFig5/ExDataFig5_locuszoom_ALP_ASGR1.pdf",
  p1,
  base_height = 3,
  base_width = 4.8,
  device = cairo_pdf
)

###############################
# rs55714927 forest plot
###############################
df <- rgsutil::read_gsfile(fm_insights$get_analysis_path("rs55714927", "tsv.bgz")) %>%
  dplyr::mutate(nlog10p = pchisq((beta_marginal / se_marginal)**2,
    1,
    log.p = TRUE,
    lower.tail =
      F
  ) / -log(10))


pd <- position_dodge(width = 0.5)

data <- dplyr::group_by(df, trait) %>%
  dplyr::filter(any(pvalue < 5e-8)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(shape = !is.na(pip) & pip > 0.1)

data.rect <-
  dplyr::distinct(data, trait) %>%
  dplyr::mutate(idx = dplyr::row_number()) %>%
  dplyr::arrange(dplyr::desc(trait)) %>%
  dplyr::mutate(
    ymin = idx - 0.5,
    ymax = idx + 0.5
  ) %>%
  dplyr::filter(idx %% 2 == 0)

p2 <-
  ggplot() +
  geom_rect(aes(
    xmin = -Inf,
    xmax = Inf,
    ymin = ymin,
    ymax = ymax
  ),
  fill = "grey90",
  data = data.rect
  ) +
  geom_point(
    aes(beta_marginal, trait, color = cohort, shape = shape),
    size = 1,
    position = pd,
    data = data
  ) +
  geom_errorbarh(
    aes(
      y = trait,
      xmin = beta_marginal - se_marginal,
      xmax = beta_marginal + se_marginal,
      color = cohort
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
    legend.position = c(1, 0.05),
    legend.justification = c(1, 0),
    axis.title.y = element_blank(),
    plot.margin = margin(0.1, 0.1, 0, 0.1, unit = "cm")
  ) +
  labs(x = "Marginal beta") +
  guides(
    color = guide_legend(title = "Cohort"),
    shape = guide_legend(title = "PIP > 0.1")
  )

cowplot::save_plot(
  "figures/Fig2/SFig_rs55714927_forest.pdf",
  p2,
  base_height = 1.5,
  base_width = 2.4,
  device = cairo_pdf
)