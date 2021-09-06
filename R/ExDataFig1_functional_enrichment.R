library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

df.csq <- rgsutil::read_gsfile(
  fm_insights$get_merged_results_path("fm_only.max_pip_per_var.consequence", "tsv.bgz")
) %>%
  dplyr::mutate(consequence = factor(consequence))

df.csq.cohort <- rgsutil::read_gsfile(
  fm_insights$get_merged_results_path("fm_only.max_pip_per_cohort.consequence", "tsv.bgz")
) %>%
  dplyr::mutate(consequence = factor(consequence))

df.csq.annot <- rgsutil::read_gsfile(
  fm_insights$get_merged_results_path("fm_only.max_pip_per_var.consequence.annot", "tsv.bgz")
)

plot_csq_fraction <- function(data, title) {
  data <-
    dplyr::mutate(data, max_pip_bin = cut(max_pip, pip_bin_breaks2)) %>%
    dplyr::filter(!is.na(max_pip_bin)) %>%
    dplyr::count(max_pip_bin, consequence) %>%
    dplyr::group_by(max_pip_bin) %>%
    dplyr::mutate(total = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      frac = n / total,
      max_pip_bin = dplyr::recode_factor(max_pip_bin, "(-Inf,0.01]" = "[0,0.01]"),
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
      ))
    ) %>%
    dplyr::arrange(desc(consequence)) %>%
    group_by(max_pip_bin) %>%
    dplyr::mutate(pos = cumsum(frac))

  label.data <- dplyr::group_by(data, max_pip_bin) %>%
    dplyr::summarize(label = ifelse(
      total[1] >= 1e4,
      scales::label_number_si(accuracy = 1)(total[1]),
      as.character(total[1])
    ))

  ggplot() +
    geom_bar(aes(
      x = factor(max_pip_bin),
      y = frac,
      fill = consequence
    ),
    stat = "identity",
    data = data
    ) +
    geom_text(
      aes(
        x = factor(max_pip_bin),
        y = 1,
        label = label
      ),
      vjust = -0.5,
      size = 2,
      data = label.data
    ) +
    scale_fill_manual(values = annot_colors) +
    my_theme +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8, margin = margin(b = 5.5)),
      plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    locusviz::or_missing(
      title != "BBJ",
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.tag.position = c(-0.05, 1)
      )
    ) +
    labs(
      x = "Best PIP bin per variant",
      y = "Proportion of variants",
      title = title,
      fill = "Annotation"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
}

plot_enrichment <- function(data, title) {
  data <- dplyr::mutate(data, max_pip_bin = cut(max_pip, pip_bin_breaks2)) %>%
    dplyr::filter(max_pip_bin %in% c("(-Inf,0.01]", "(0.9,1]")) %>%
    dplyr::count(max_pip_bin, consequence) %>%
    dplyr::group_by(max_pip_bin) %>%
    dplyr::mutate(
      total = sum(n),
      frac = n / total
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      max_pip_bin = dplyr::recode_factor(max_pip_bin, "(-Inf,0.01]" = "[0,0.01]"),
      consequence = ordered(consequence, levels = annot_levels)
    ) %>%
    tidyr::pivot_wider(
      id_cols = "consequence",
      names_from = "max_pip_bin",
      values_from = c("n", "total")
    ) %>%
    dplyr::group_split(consequence) %>%
    purrr::map_dfr(function(data) {
      m <- with(data, matrix(
        c(
          `total_[0,0.01]` - `n_[0,0.01]`,
          `n_[0,0.01]`,
          `total_(0.9,1]` - `n_(0.9,1]`,
          `n_(0.9,1]`
        ),
        nrow = 2,
        byrow = T
      ))
      measure <- epitools::riskratio(m, method = "boot")$measure
      tibble::tibble(
        consequence = data$consequence,
        enrichment = measure[2, "estimate"],
        lower = measure[2, "lower"],
        upper = measure[2, "upper"]
      )
    })

  cohort <- ifelse(title == "All cohorts", "All", title)
  dplyr::mutate(data, cohort = .env$cohort) %>%
    dplyr::select(cohort, consequence, enrichment, lower, upper) %>%
    write.table(
      sprintf(
        "./tables/STable_enrichment_%s.tsv",
        cohort
      ),
      quote = F,
      row.names = F,
      sep = "\t"
    )

  ggplot(data, aes(color = consequence)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_point(aes(consequence, enrichment)) +
    geom_errorbar(aes(
      x = consequence,
      ymin = lower,
      ymax = upper
    ), width = 0) +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm"),
      legend.title = element_blank()
    ) +
    locusviz::or_missing(
      title != "BBJ",
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.tag.position = c(-0.05, 1)
      )
    ) +
    locusviz::or_missing(title != "All cohorts", theme(legend.position = "none")) +
    scale_color_manual(values = annot_colors) +
    scale_y_log10() +
    coord_cartesian(ylim = c(0.3, 150)) +
    labs(y = "Enrichment", color = "Annotation")
}

panels <-
  c(
    purrr::map(cohorts, function(cohort) {
      dplyr::filter(df.csq.cohort, cohort == .env$cohort) %>%
        plot_csq_fraction(title = ifelse(cohort == "FG", "FinnGen", cohort))
    }),
    list(plot_csq_fraction(df.csq, title = "All cohorts")),
    purrr::map(cohorts, function(cohort) {
      dplyr::filter(df.csq.cohort, cohort == .env$cohort) %>%
        plot_enrichment(title = cohort)
    }),
    list(plot_enrichment(df.csq, "All cohorts"))
  )


################################################

continuous_annots <- c(
  "Backgrd_Selection_Stat",
  "BLUEPRINT_FE_META_TISSUE_DNAMETH_MaxCPP",
  "BLUEPRINT_FE_META_TISSUE_H3K27ac_MaxCPP",
  "BLUEPRINT_FE_META_TISSUE_H3K4me1_MaxCPP",
  "GTEx_FE_META_TISSUE_GE_MaxCPP",
  "ASMC",
  "alleleage",
  "Human_Enhancer_Villar_Species_Enhancer_Count"
)
binary_annots <- setdiff(
  colnames(df.csq.annot),
  c(
    "locus",
    "alleles",
    "max_pip",
    "consequence",
    "clinvar",
    "Coding_UCSC",
    "Promoter_UCSC",
    "UTR_3_UCSC",
    "UTR_5_UCSC",
    "CA_H3K27ac_Ulirsch",
    "DHSmerged_Ulirsch",
    "Roadmap_H3K27ac_Ulirsch",
    continuous_annots
  )
)

df.csq.annot.filt <- dplyr::mutate(df.csq.annot, max_pip_bin = cut(max_pip, pip_bin_breaks2)) %>%
  dplyr::filter(max_pip_bin %in% c("(-Inf,0.01]", "(0.9,1]")) %>%
  dplyr::mutate(max_pip_bin = dplyr::recode_factor(max_pip_bin, "(-Inf,0.01]" = "[0,0.01]"))

plot_binary_enrichment <- function(data, cohort) {
  data <- dplyr::group_split(data, max_pip_bin) %>%
    purrr::map_dfr(function(data) {
      tibble::tibble(
        max_pip_bin = rep(data$max_pip_bin[1], 2 * length(binary_annots)),
        consequence = rep(c("All", "Non-genic"), each = length(binary_annots)),
        annot = rep(binary_annots, 2),
        n = c(colSums(data[, binary_annots]), colSums(data[which(data$consequence == "Non-genic"), binary_annots])),
        total = rep(c(
          nrow(data), sum(data$consequence == "Non-genic")
        ), each = length(binary_annots))
      )
    }) %>%
    tidyr::pivot_wider(
      id_cols = c("consequence", "annot"),
      names_from = "max_pip_bin",
      values_from = c("n", "total")
    ) %>%
    dplyr::group_split(consequence, annot) %>%
    purrr::map_dfr(function(data) {
      m <- with(data, matrix(
        c(
          `total_[0,0.01]` - `n_[0,0.01]`,
          `n_[0,0.01]`,
          `total_(0.9,1]` - `n_(0.9,1]`,
          `n_(0.9,1]`
        ),
        nrow = 2,
        byrow = T
      ))
      measure <- epitools::riskratio(m, method = "boot")$measure
      tibble::tibble(
        consequence = data$consequence[1],
        annot = data$annot[1],
        enrichment = measure[2, "estimate"],
        lower = measure[2, "lower"],
        upper = measure[2, "upper"]
      )
    })

  if (cohort == "All") {
    annot_levels <- annot[order(enrichment[which(consequence == "All")], decreasing = T)]
  } else {
    annot_levels <- read.table("./tables/STable_ldsc_binary_enrichment_All.tsv",
      T,
      sep = "\t"
    ) %>%
      dplyr::filter(consequence == "All") %>%
      dplyr::pull(annot)
  }
  data <- dplyr::mutate(data, annot = ordered(annot, levels = annot_levels))

  dplyr::mutate(data, cohort = .env$cohort) %>%
    dplyr::select(cohort, consequence, annot, enrichment, lower, upper) %>%
    dplyr::arrange(cohort, consequence, annot) %>%
    write.table(
      sprintf(
        "./tables/STable_ldsc_binary_enrichment_%s.tsv",
        cohort
      ),
      quote = F,
      row.names = F,
      sep = "\t"
    )

  pd <- position_dodge(width = 1)

  plt <-
    ggplot(data, aes(color = consequence)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_point(aes(annot, enrichment), position = pd) +
    geom_errorbar(aes(
      x = annot,
      ymin = lower,
      ymax = upper
    ),
    width = 0,
    position = pd
    ) +
    my_theme +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      axis.title.x = element_blank()
    ) +
    scale_color_manual(values = BuenColors::jdb_palette("brewer_celsius")[c(8, 2)]) +
    labs(y = "Enrichment", color = "Variants")
  if (cohort == "All") {
    plt <- plt +
      scale_y_log10(breaks = c(1, 5, 10)) +
      coord_cartesian(ylim = c(0.3, 15))
  } else {
    plt <- plt +
      scale_y_log10(breaks = c(1, 10, 100)) +
      coord_cartesian(ylim = c(0.3, 100))
  }
  if (cohort %in% c("BBJ", "FG")) {
    plt <- plt + theme(axis.text.x = element_blank())
  }
  return(plt)
}

plot_continous_enrichment <- function(data) {
  purrr::map(setdiff(
    continuous_annots,
    "Human_Enhancer_Villar_Species_Enhancer_Count"
  ), function(annot) {
    ggplot(data, aes_string("max_pip_bin", annot)) +
      geom_violin() +
      my_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()
      ) +
      locusviz::stat_summary_irq(color = "black") +
      locusviz::or_missing(annot %in% c("ASMC", "alleleage"), scale_y_log10())
  }) %>% purrr::reduce(`+`) + patchwork::plot_layout(nrow = 2)
}

layout <- "
ABCD
EFGH
IIII
"

plt <- append(panels, list(plot_binary_enrichment(df.csq.annot.filt, cohort = "All"))) %>%
  purrr::reduce(`+`) +
  patchwork::plot_annotation(tag_levels = "a") +
  patchwork::plot_layout(design = layout)

cowplot::save_plot(
  "figures/ExDataFig1/ExDataFig1_functional_enrichment.pdf",
  plt,
  base_height = 8.4,
  base_width = 7.2,
  device = cairo_pdf
)
