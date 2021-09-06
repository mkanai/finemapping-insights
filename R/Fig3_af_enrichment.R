library(dplyr)
library(ggplot2)
library(patchwork)
library(magrittr)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")


plot_hist_and_pip <- function(df,
                              xstr,
                              data_type,
                              xlim,
                              xlab,
                              ylab,
                              fill_color,
                              title_y,
                              breaks = c(0.01, 0.1, 1, 10, 100, 1000)) {
  if (is.null(xlim)) {
    xlim <- range(df[[xstr]], na.rm = T)
    print(xlim)
  }

  data <-
    dplyr::mutate(
      df,
      max_pip_bin = cut(max_pip, pip_bin_breaks2),
      max_pip_bin = dplyr::recode_factor(max_pip_bin, "(-Inf,0.01]" = "[0,0.01]"),
      enrichment_bin = cut(!!as.symbol(xstr), c(-Inf, 10, Inf))
    )

  df.p1.label <- dplyr::group_by(data, enrichment_bin) %>%
    dplyr::summarize(count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      frac = count / sum(count),
      x = 10**c(`(-Inf,10]` = 0.75, `(10,Inf]` = 1.5)[enrichment_bin],
      y = 0
    )

  binwidth <- 0.075
  p1 <-
    ggplot(data, aes_string(xstr)) +
    geom_histogram(aes(fill = cut(..x.., c(-Inf, 1, Inf))), binwidth = binwidth) +
    geom_density(aes(y = ..count.. * binwidth), size = 0.2) +
    geom_vline(
      xintercept = 10,
      linetype = "dashed",
      color = "grey50",
      size = 0.2
    ) +
    geom_text(aes(
      x = x,
      y = y,
      label = scales::percent(frac, accuracy = 1)
    ),
    data = df.p1.label,
    vjust = -4,
    size = 2
    ) +
    scale_x_log10(
      label = scales::comma_format(accuracy = 0.0001, drop0trailing = TRUE),
      lim = xlim,
      breaks = breaks
    ) +
    scale_y_continuous(label = scales::comma_format()) +
    my_theme +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    labs(y = sprintf(
      "# %s variants",
      ifelse(data_type == "exomes", "coding", "non-coding")
    )) +
    scale_fill_manual(
      values = c(
        "(1, Inf]" = fill_color,
        "(-Inf,1]" = "grey50"
      ),
      guide = FALSE
    )

  df.p2.label <-
    tidyr::drop_na(data, max_pip) %>%
    dplyr::group_by(max_pip_bin, enrichment_bin) %>%
    dplyr::summarize(count = n()) %>%
    dplyr::group_by(max_pip_bin) %>%
    dplyr::mutate(
      frac = count / sum(count),
      x = 10**c(`(-Inf,10]` = 0.75, `(10,Inf]` = 1.5)[enrichment_bin]
    ) %>%
    dplyr::ungroup()


  p2 <- tidyr::drop_na(data, max_pip) %>%
    ggplot(aes_string(xstr, "max_pip_bin")) +
    # shading: https://stackoverflow.com/questions/49961582/how-shade-area-under-ggridges-curve
    ggridges::stat_density_ridges(aes(fill = cut(..x.., c(-Inf, 1, Inf))),
      geom = "density_ridges_gradient",
      size = 0.2
    ) +
    geom_vline(
      xintercept = 10,
      linetype = "dashed",
      color = "grey50",
      size = 0.2
    ) +
    geom_text(
      aes(
        x = x,
        y = max_pip_bin,
        label = scales::percent(frac, accuracy = 1)
      ),
      data = df.p2.label,
      vjust = -2,
      size = 2
    ) +
    scale_x_log10(
      label = scales::comma_format(accuracy = 0.0001, drop0trailing = TRUE),
      lim = xlim,
      breaks = breaks
    ) +
    scale_fill_manual(
      values = c(
        "(1, Inf]" = fill_color,
        "(-Inf,1]" = "grey50"
      ),
      guide = FALSE
    ) +
    my_theme +
    labs(
      x = xlab,
      y = ifelse(title_y, ylab, stringr::str_remove(ylab, "Best PIP bin\n"))
    )

  if (!title_y) {
    p1 <- p1 + theme(axis.title.y = element_blank())
    p2 <- p2
  }

  return(list(p1, p2))
}

plot_enrichment_fig <- function(df,
                                pop1,
                                pop2,
                                data_type,
                                xlim_enrichment = NULL,
                                xlim_age = NULL,
                                title_y = TRUE,
                                age_ranges = NULL) {
  enrichment_str <- sprintf("enrichment_%s_%s_pseudo", pop1, pop2)

  # color = ifelse(pop1 == "jpn", gnomad_pop_colors["oea"], gnomad_pop_colors["nfsee"])
  color <- ifelse(pop1 == "jpn", cohort_colors["BBJ"], cohort_colors["FG"])

  plt <- c(
    plot_hist_and_pip(
      df,
      enrichment_str,
      data_type,
      xlim = xlim_enrichment,
      xlab = sprintf(
        "AF Enrichment (%s / %s)",
        stringr::str_to_upper(pop1),
        stringr::str_to_upper(pop2)
      ),
      ylab = sprintf("Best PIP bin (%s)", ifelse(pop1 == "jpn", "BBJ", "FG")),
      fill_color = color,
      title_y = title_y
    ),
    list(
      plot_age_fig(
        df,
        pop1,
        pop2,
        data_type,
        xlim = xlim_age,
        title_y = title_y,
        color = color,
        ranges = age_ranges
      )
    )
  )
  return(plt)
}

plot_maf_fig <- function(df,
                         pop,
                         data_type,
                         xlim = NULL,
                         title_y = TRUE) {
  print(xlim)
  af_str <- sprintf("AF_%s", pop)
  plot_hist_and_pip(
    df %>% dplyr::mutate(maf = 0.5 - abs(0.5 - !!as.symbol(af_str))) %>% dplyr::filter(maf > 0),
    "maf",
    data_type,
    xlim = xlim,
    xlab = sprintf(
      "MAF (%s)",
      stringr::str_to_upper(pop)
    ),
    ylab = sprintf("Max PIP bin\n(%s)", ifelse(pop == "jpn", "BBJ", "FG")),
    fill_color = gnomad_pop_colors[pop],
    title_y = title_y,
    breaks = c(10**seq(-4, -1), 0.5)
  )
}

plot_age_fig <- function(df,
                         pop1,
                         pop2,
                         data_type,
                         xlim = NULL,
                         title_y = TRUE,
                         color = "black",
                         ranges = NULL) {
  breaks <- 10**seq(0, 6)
  minor_breaks <- rep(1:9, length(breaks)) * rep(breaks, each = 9)

  enrichment_breaks <- c(-Inf, 10, Inf)

  enrichment_str <- sprintf("enrichment_%s_%s_pseudo", pop1, pop2)

  df <-
    dplyr::mutate(
      df,
      enrichment_bin = cut(!!as.symbol(enrichment_str), enrichment_breaks),
      enrichment_bin = forcats::fct_recode(enrichment_bin, "> 10" = "(10, Inf]", "≤ 10" = "(-Inf,10]"),
    ) %>%
    tidyr::drop_na(age_mode, enrichment_bin)

  plt <- ggplot() +
    locusviz::or_missing(!is.null(ranges), purrr::map(ranges, function(x) {
      geom_rect(
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax,
          fill = fill
        ),
        data = tibble::tibble(
          xmin = x[1],
          xmax = x[2],
          ymin = -Inf,
          ymax = Inf
        ),
        fill = "#F2F2F2"
      )
    })) +
    geom_hline(
      yintercept = 0.5,
      linetype = "dashed",
      color = "grey50",
      size = 0.1
    ) +
    annotation_logticks(sides = "b", size = 0.2) +
    ggrastr::rasterize(stat_ecdf(aes(age_mode, color = enrichment_bin), data = df),
      dpi = 300
    ) +
    scale_x_log10(
      label = scales::comma_format(accuracy = 0.0001, drop0trailing = TRUE),
      limits = xlim,
      breaks = breaks,
      minor_breaks = minor_breaks
    ) +
    my_theme +
    theme(
      legend.position = c(1, 0.15),
      legend.justification = c(1, 0)
    ) +
    scale_color_manual(
      values = c(
        "> 10" = color,
        "≤ 10" = "grey50"
      ),
      breaks = c("> 10")
    ) +
    labs(
      x = "Allele age (generations)",
      y = "Cumulative distribution",
      color = "AF Enrichment"
    ) +
    locusviz::or_missing(!title_y, theme(axis.title.y = element_blank()))

  return(plt)
}


purrr::map(c("exomes", "genomes"), function(data_type) {
  df_jpn <- rgsutil::read_gsfile(fm_insights$get_analysis_path(
    sprintf("af_enrichment.%s.BBJ.max_pip", data_type),
    "tsv.bgz"
  )) %>%
    dplyr::filter(!is.na(AF_jpn))
  df_fin <- rgsutil::read_gsfile(fm_insights$get_analysis_path(
    sprintf("af_enrichment.%s.FG.max_pip", data_type),
    "tsv.bgz"
  )) %>%
    dplyr::filter(!is.na(AF_fin))

  if (data_type == "genomes") {
    print(table(df_jpn$coding_tagging_jpn))
    print(table(df_fin$coding_tagging_fin))
    df_jpn <- dplyr::filter(df_jpn, !coding_tagging_jpn)
    df_fin <- dplyr::filter(df_fin, !coding_tagging_fin)
  }

  pop2 <- c("nfsee", "njkea")
  xlim <- dplyr::bind_rows(
    df_jpn[, sprintf("enrichment_jpn_%s_pseudo", pop2)],
    df_fin[, sprintf("enrichment_fin_%s_pseudo", pop2)]
  ) %>% range(na.rm = T)

  xlim_enrichment <- c(0.01, 1000)
  xlim_age <- c(50, 50000)

  # enrichment fig

  plt <- purrr::map2(c(
    plot_enrichment_fig(
      df_fin,
      "fin",
      "nfsee",
      data_type,
      xlim_enrichment,
      xlim_age,
      age_ranges = NULL # list(c(100, 120), c(200, 400))
    ),
    plot_enrichment_fig(
      df_jpn,
      "jpn",
      "njkea",
      data_type,
      xlim_enrichment,
      xlim_age,
      title_y = FALSE,
      age_ranges = NULL # list(c(3000, 12000) / 25, c(8000, 9000) / 25)
    )
  ), as.list(matrix(
    letters[1:6],
    byrow = T, nrow = 3
  )), function(panel, tag) {
    panel + labs(tag = tag)
  }) %>% purrr::reduce(`+`) + patchwork::plot_layout(ncol = 2, byrow = F)

  if (data_type == "exomes") {
    fname <- sprintf("figures/Fig3.pdf", data_type)
  } else {
    fname <- sprintf("figures/ExDataFig9_af_enrichment_%s.pdf", data_type)
  }

  cowplot::save_plot(
    fname,
    plt,
    base_height = 4.8,
    base_width = 4.8,
    device = cairo_pdf
  )

  # maf fig
  # xlim = c(1e-4, 0.5)
  # p_maf = plot_maf_fig(df_fin, "fin", data_type, xlim) |
  #   plot_maf_fig(df_jpn, "jpn", data_type, xlim)
  # p_maf
  # cowplot::save_plot(
  #   sprintf('figures/Fig3/SFig_maf_%s.pdf', data_type),
  #   p_maf,
  #   base_height = 2.4,
  #   base_width = 4.8,
  #   device = cairo_pdf
  # )
})

###########################
# fisher's exact test

dplyr::mutate(df_fin,
  max_pip_bin = cut(max_pip, pip_bin_breaks2),
  enrichment_bin = cut(enrichment_fin_nfsee_pseudo, c(-Inf, 10, Inf))
) %$%
  table(max_pip_bin, enrichment_bin) %>%
  as.matrix() %>%
  magrittr::extract(c(1, 5), 1:2) %>%
  fisher.test()


dplyr::mutate(df_jpn,
  max_pip_bin = cut(max_pip, pip_bin_breaks2),
  enrichment_bin = cut(enrichment_jpn_njkea_pseudo, c(-Inf, 10, Inf))
) %$%
  table(max_pip_bin, enrichment_bin) %>%
  as.matrix() %>%
  magrittr::extract(c(1, 5), 1:2) %>%
  fisher.test() %>%
  magrittr::extract2("p.value")


#########################
# AFE > 10 and < 0.1 ratio

df_fin_jpn <-
  purrr::map_dfr(c("exomes", "genomes"), function(data_type) {
    df_jpn <- rgsutil::read_gsfile(fm_insights$get_analysis_path(
      sprintf("af_enrichment.%s.BBJ.max_pip", data_type),
      "tsv.bgz"
    )) %>%
      dplyr::filter(!is.na(AF_jpn)) %>%
      dplyr::mutate(
        data_type = .env$data_type,
        pop = "jpn",
        enrichment_pseudo = enrichment_jpn_njkea_pseudo
      )
    df_fin <- rgsutil::read_gsfile(fm_insights$get_analysis_path(
      sprintf("af_enrichment.%s.FG.max_pip", data_type),
      "tsv.bgz"
    )) %>%
      dplyr::filter(!is.na(AF_fin)) %>%
      dplyr::mutate(
        data_type = .env$data_type,
        pop = "fin",
        enrichment_pseudo = enrichment_fin_nfsee_pseudo
      )

    return(
      dplyr::bind_rows(df_jpn, df_fin) %>%
        dplyr::select(data_type, pop, max_pip, enrichment_pseudo)
    )
  })

dplyr::group_by(df_fin_jpn, pop) %>%
  dplyr::summarize(
    count_enrichment_10 = sum(enrichment_pseudo > 10),
    count_enrichment_01 = sum(enrichment_pseudo < 0.1),
    frac_enrichment_10 = count_enrichment_10 / n(),
    frac_enrichment_01 = count_enrichment_01 / n()
  )


# fisher's test
dplyr::group_split(df_fin_jpn, pop, data_type) %>%
  purrr::map_dfr(function(data) {
    pvalue <-
      dplyr::mutate(data,
        max_pip_bin = cut(max_pip, pip_bin_breaks2),
        enrichment_bin = cut(enrichment_pseudo, c(-Inf, 10, Inf))
      ) %$%
      table(max_pip_bin, enrichment_bin) %>%
      as.matrix() %>%
      magrittr::extract(c(1, 5), 1:2) %>%
      fisher.test() %>%
      magrittr::extract2("p.value")
    tibble::tibble(
      pop = data$pop[1],
      data_type = data$data_type[1],
      pvalue = pvalue
    )
  })

# enrichment for coding with AFE > 10
dplyr::mutate(df_fin_jpn, enrichment_bin = cut(enrichment_pseudo, c(-Inf, 10, Inf))) %>%
  dplyr::count(pop, data_type, enrichment_bin) %>%
  tidyr::pivot_wider(
    id_cols = "pop",
    names_from = c("data_type", "enrichment_bin"),
    values_from = "n"
  ) %>%
  dplyr::group_split(pop) %>%
  purrr::map_dfr(function(data) {
    m <- with(data, matrix(
      c(
        `genomes_(-Inf,10]`,
        `exomes_(-Inf,10]`,
        `genomes_(10, Inf]`,
        `exomes_(10, Inf]`
      ),
      nrow = 2,
      byrow = T
    ))
    measure <- epitools::riskratio(m, method = "boot")$measure
    tibble::tibble(
      pop = data$pop[1],
      enrichment = measure[2, "estimate"],
      lower = measure[2, "lower"],
      upper = measure[2, "upper"]
    )
  })