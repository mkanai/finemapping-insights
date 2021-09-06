# reticulate::use_python("~/.anyenv/envs/pyenv/versions/anaconda3-2020.02/bin/python", required = TRUE)
fm_insights <- reticulate::import("fm_insights")

# constants
cohorts <- c("BBJ", "FG", "UKBB")
pip_bin_breaks <- c(0, 0.01, 0.1, 0.5, 0.9, 1.0)
pip_bin_breaks2 <- c(-Inf, 0.01, 0.1, 0.5, 0.9, 1.0)

cohort_colors <- c(
  BBJ = "#df5f9b",
  FG = "#514c9f",
  FinnGen = "#514c9f",
  UKBB = "#006683",
  Any = "black"
)
cohort_shapes <- c(BBJ = 16, FG = 15, UKBB = 18)
# https://github.com/macarthur-lab/gnomad_lof/blob/master/R/constants.R
gnomad_pop_colors <- c(
  afr = "#941494",
  oea = "#108C44",
  nfsee = "#6AA5CD",
  fin = "#002F6C",
  jpn = "#BC002D"
)

domains <- c(
  "Metabolic",
  "Lipids",
  "Cardiovascular",
  "Immunological",
  "Hematopoietic",
  "Hepatic",
  "Renal",
  "Skeletal",
  "Neurological",
  "Behavioral",
  "Other"
)
domain_colors <- c(
  "#4C72AE",
  "#E4812F",
  "#F3DE67",
  "#6BBCCC",
  "#5AA245",
  "#BA2E2C",
  "#8C61B6",
  "#80584E",
  "#CE72BB",
  "#BACD3C",
  "#7F7F7F"
)
names(domain_colors) <- domains

# functional annotations
annot_levels <- c(
  "pLoF",
  "Missense",
  "Synonymous",
  "UTR5",
  "UTR3",
  "Promoter",
  "CRE",
  "Non-genic"
)
annot_colors <- c(
  BuenColors::jdb_palette("brewer_celsius")[c(8, 6)],
  "#AAAAAA",
  BuenColors::jdb_palette("brewer_purple")[c(6, 4)],
  BuenColors::jdb_palette("brewer_blue")[c(8, 6, 3)]
)
names(annot_colors) <- annot_levels


coding_colors <- c(
  pLoF = BuenColors::jdb_palette("brewer_celsius")[8],
  Missense = BuenColors::jdb_palette("brewer_celsius")[6],
  `Missense (damaging)` = BuenColors::jdb_palette("brewer_celsius")[6],
  `Missense (benign)` = "#BFB099",
  `Missense (unknown)` = "grey80",
  Synonymous = "#AAAAAA"
)

clinvar_colors <- c(
  BuenColors::jdb_palette("brewer_violet")[c(8, 4)],
  "grey80",
  "#AAAAAA"
)
names(clinvar_colors) <- c("Pathogenic", "Non-pathogenic", "NA", "Synonymous")



my_theme <-
  BuenColors::pretty_plot(fontsize = 8) +
  BuenColors::L_border() +
  theme(
    plot.background = element_blank(),
    plot.margin = margin(0, 0.1, 0, 0.1, unit = "cm"),
    plot.tag = element_text(face = "bold"),
    plot.tag.position = c(0, 1),
    plot.title = element_text(hjust = 4e-3, margin = margin(b = -12)),
    # legend.position = "none",
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.title = element_text(margin = margin(0, 0, 0, 0)),
    legend.background = element_blank(),
    legend.key.size = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0)
  )

read_trait_summary <- function() {
  path <-
    stringr::str_glue("gs://{fm_insights$bucket}/metadata/trait_summary.txt")
  return(rgsutil::read_gsfile(path))
}


compute_p_het <- function(beta1, se1, beta2, se2) {
  purrr::pmap_dbl(data.frame(beta1, se1, beta2, se2), function(beta1, se1, beta2, se2) {
    inv_se2 <- 1 / (c(se1, se2)**2)
    unnorm_beta <- c(beta1, beta2) * inv_se2
    beta_meta <- sum(unnorm_beta) / sum(inv_se2)
    q_meta <- sum((c(beta1, beta2) - beta_meta)**2 * inv_se2)
    return(pchisq(q_meta, 1, lower.tail = F))
  })
}

plot_multi_locuszoom <- function(trait,
                                 lead_variant,
                                 rsid,
                                 window,
                                 cohorts = c("BBJ", "FG", "UKBB"),
                                 df.csm = NULL,
                                 additional_highlight_pos = NULL,
                                 additional_highlight_labels = NULL,
                                 fm.height = 0.3) {
  lead_pos <- locusviz::parse_variant(lead_variant)$position
  chromosome <- locusviz::parse_variant(lead_variant)$chromosome
  if (length(window) == 1) {
    window <- rep(window, 2)
  }
  start <- lead_pos - window[1]
  end <- lead_pos + window[2]

  highlight_pos <- lead_pos
  if (!is.null(additional_highlight_pos)) {
    highlight_pos <- c(highlight_pos, additional_highlight_pos)
  }

  cs.colors <- dplyr::filter(df.csm, trait == .env$trait & cohort %in% cohorts) %>%
    dplyr::mutate(locusviz::parse_variant(variant)) %>%
    dplyr::filter(chromosome == .env$chromosome & .env$start * 0.99 <= position & position <= .env$end * 1.01) %>%
    dplyr::arrange(abs(position - position[which(.$variant == lead_variant)[1]])) %$%
    locusviz::get_cs_color_mapping(.$csm_id, highlight_cs_ids = .$csm_id[which(.$variant == lead_variant)])

  panels <-
    purrr::flatten(purrr::map(cohorts, function(cohort) {
      data <-
        rgsutil::read_gsfile(fm_insights$get_locus_txt_path(cohort, trait, lead_variant)) %>%
        dplyr::left_join(df.csm) %>%
        dplyr::mutate(locusviz::parse_variant(variant)) %>%
        locusviz::preprocess(
          lead_variant = lead_variant,
          beta_col = "beta_marginal",
          se_col = "se_marginal",
          cs_id_col = "csm_id",
          r2_col = "lead_r2"
        )

      p_manhattan <- locusviz::plot_manhattan_panel(
        data,
        highlight_pos = highlight_pos,
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
          legend.position = if (cohort == cohorts[1]) {
            c(1, 1.1)
          } else {
            "none"
          }
        ) + theme(
          plot.title = element_text(
            hjust = 0.01,
            margin = margin(b = -12),
            size = 6
          ),
          legend.margin = margin(0, 0, 0, 0),
          legend.title = element_text(size = 6),
          legend.key.size = unit(0.15, units = "cm")
        ),
        rasterize = TRUE
      )

      if (cohort == cohorts[1]) {
        label_pos <- data$nlog10p[which(data$position == lead_pos)[1]]
        if (label_pos < 10) {
          label_pos <- label_pos + max(data$nlog10p, na.rm = T) / 10
        }
        p_manhattan <- p_manhattan +
          geom_text(
            aes(x, y, label = label),
            data = tibble::tibble(
              x = lead_pos,
              y = label_pos,
              label = rsid
            ),
            size = 2,
            hjust = -0.1
          )
        if (!is.null(additional_highlight_labels)) {
          p_manhattan <- p_manhattan +
            purrr::pmap(tibble::tibble(position = additional_highlight_pos, label = additional_highlight_labels), function(position, label) {
              label_pos <- data$nlog10p[which(data$position == position)[1]]
              if (label_pos < 10) {
                label_pos <- label_pos + max(data$nlog10p, na.rm = T) / 10
              }
              geom_text(
                aes(x, y, label = label),
                data = tibble::tibble(
                  x = .env$position,
                  y = label_pos,
                  label = .env$label
                ),
                size = 2,
                hjust = -0.1
              )
            })
        }
      }

      p_fm <- locusviz::plot_fm_panel(
        data,
        highlight_pos = highlight_pos,
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
        rasterize = TRUE,
        cs.colors = cs.colors
      )

      return(list(p_manhattan, p_fm))
    }))

  panels <- c(panels, list(
    locusviz::plot_gene_panel(
      chromosome,
      start,
      end,
      highlight_pos = highlight_pos,
      fontsize = 6,
      point.size = 1.5,
      label.size = 1.5,
      length = unit(0.05, "cm")
    )
  ))

  plt <- purrr::reduce(c(panels), `+`) + patchwork::plot_layout(ncol = 1, heights = c(rep(c(1, fm.height), length(cohorts)), 0.1))
  return(plt)
}