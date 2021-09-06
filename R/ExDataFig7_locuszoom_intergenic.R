library(ggplot2)
library(magrittr)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

df.pip01.every <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("shard_trait.pip09.any.pip01.every", "tsv"))
df.csm <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("fm_only.csm_id", "tsv.bgz")) %>%
  dplyr::select(cohort, trait, region, variant, susie.cs_id, csm_id)

txdb <- locusviz::load_txdb("hg19")
df.dist <-
  dplyr::distinct(df.pip01.every, variant) %>%
  purrr::pmap_dfr(function(variant) {
    lead_variant <- variant[1]
    window <- 5000000
    lead_pos <- locusviz::parse_variant(lead_variant)$position
    chromosome <- locusviz::parse_variant(lead_variant)$chromosome
    start <- lead_pos - window
    end <- lead_pos + window

    locusviz::compute_distance_to_gene(txdb, paste0("chr", chromosome), start, end, ref_position = lead_pos) %>%
      dplyr::mutate(lead_variant = variant)
  })

df.intergenic <-
  dplyr::group_by(df.dist, lead_variant) %>%
  dplyr::filter(score == min(score)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(score > 250000) %>%
  dplyr::rename(variant = lead_variant)

write.table(df.intergenic, "./data/intergenic_variants.dist.txt", quote = F, row.names = F, sep = "\t")

v_intergenic <- dplyr::pull(df.intergenic, lead_variant)

v_intergenic_to_plot <- c(
  "8:128077146:G:A",
  "1:199010721:C:T",
  "2:227326633:A:T",
  "10:9329035:AT:A"
)

window_intergenic <- c(
  "1:199010721:C:T" = 300000,
  "2:227326633:A:T" = 300000,
  "8:128077146:G:A" = 350000,
  "10:9329035:AT:A" = 450000
)
tags <- c(
  "8:128077146:G:A" = "a",
  "1:199010721:C:T" = "b",
  "2:227326633:A:T" = "c",
  "10:9329035:AT:A" = "d"
)

plts <-
  dplyr::filter(
    df.pip01.every,
    variant %in% v_intergenic_to_plot &
      trait %in% c("MCV", "Height", "PrC", "Ca")
  ) %>%
  dplyr::select(trait, variant, rsid) %>%
  purrr::pmap(function(trait, variant, rsid) {
    if (variant %in% c("2:227326633:A:T", "8:128077146:G:A")) {
      cohorts <- c("UKBB", "FG")
    } else {
      cohorts <- c("BBJ", "UKBB")
    }

    plt <-
      plot_multi_locuszoom(
        trait = trait,
        lead_variant = variant,
        rsid = rsid,
        window = window_intergenic[variant],
        cohorts = cohorts,
        df.csm = df.csm,
        fm.height = 0.3
      )
    plt[[1]] <- plt[[1]] + labs(tag = tags[variant])
    return(plt)
  })

plt <- (plts[[3]] | plts[[1]]) / (plts[[2]] | plts[[4]])
plt

cowplot::save_plot(
  "figures/ExDataFig7_locuszoom_intergenic.pdf",
  plt,
  base_height = 7.2,
  base_width = 7.2,
  device = cairo_pdf
)