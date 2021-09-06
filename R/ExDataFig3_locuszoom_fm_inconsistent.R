library(ggplot2)
library(magrittr)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

df.in_cs <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("shard_trait.pip001.pop", "tsv.bgz"))
# df.cs_max_regions = rgsutil::read_gsfile(fm_insights$get_merged_results_path("cs_max_regions", "tsv.bgz"))
df.inconsistent.in <- rgsutil::read_gsfile(
  fm_insights$get_merged_results_path("shard_trait.fm_inconsistent.variant_list", "tsv")
)
df.csm <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("fm_only.csm_id", "tsv.bgz")) %>%
  dplyr::select(cohort, trait, region, variant, susie.cs_id, csm_id)

df.inconsistent <-
  dplyr::select(df.inconsistent.in, trait, variant) %>%
  dplyr::distinct() %>%
  dplyr::left_join(df.in_cs)

variants_to_plot <- c("2:219357858:A:G", "4:45164637:T:G", "11:2165576:G:A", "7:19248278:G:T")

window <- c(
  "11:2165576:G:A" = 50000,
  "7:19248278:G:T" = 200000,
  "4:45164637:T:G" = 50000,
  "2:219357858:A:G" = 200000
)
tags <- c(
  "11:2165576:G:A" = "a",
  "7:19248278:G:T" = "b",
  "4:45164637:T:G" = "c",
  "2:219357858:A:G" = "d"
)
rsids <- c(
  "11:2165576:G:A" = "rs35506085",
  "7:19248278:G:T" = "rs17140875",
  "4:45164637:T:G" = "rs1996023",
  "2:219357858:A:G" = "rs495855"
)

plts <-
  dplyr::filter(df.inconsistent, variant %in% variants_to_plot & trait != "BW") %>%
  dplyr::group_by(trait, variant) %>%
  dplyr::summarize(cohort = list(
    c("BBJ")[!is.na(pvalue.BBJ)],
    c("FG")[!is.na(pvalue.FG)],
    c("UKBB")[!is.na(pvalue.UKBB)]
  )) %>%
  tidyr::unnest(cohort) %>%
  dplyr::summarize(cohorts = list(cohort)) %>%
  purrr::pmap(function(trait, variant, cohorts) {
    plt <-
      plot_multi_locuszoom(
        trait = trait,
        lead_variant = variant,
        rsid = rsids[variant],
        window = window[variant],
        cohorts = cohorts,
        df.csm = df.csm,
        fm.height = 0.3
      )
    plt[[1]] <- plt[[1]] + labs(tag = tags[variant])
    return(plt)
  })

plt <- (plts[[2]] | plts[[4]]) / (plts[[1]] | plts[[3]])

cowplot::save_plot(
  "figures/ExDataFig3_locuszoom_fm_inconsistent.pdf",
  plt,
  base_height = 9.6,
  base_width = 7.2,
  device = cairo_pdf
)
