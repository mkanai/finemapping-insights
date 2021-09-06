library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

prettify_most_severe <- function(most_severe) {
  stringr::str_remove(most_severe, "_variant$") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_sentence() %>%
    stringr::str_replace("utr$", "UTR")
}

prettify_array <- function(array) {
  stringr::str_replace_all(array, "[\\[\\]\"]", "") %>%
    stringr::str_replace_all(",", ", ")
}

combine_values <- function(bbj,
                           fg,
                           ukbb,
                           format = NULL) {
  if (!is.null(format)) {
    formatter <- function(x) {
      sprintf(format, x)
    }
  } else {
    formatter <- function(x) {
      x
    }
  }
  dplyr::bind_cols(x = bbj, y = fg, z = ukbb) %>%
    purrr::pmap_chr(function(x, y, z) {
      purrr::keep(list(x, y, z), function(x) {
        !is.na(x)
      }) %>%
        purrr::map_chr(formatter) %>%
        stringr::str_c(collapse = ", ")
    })
}

get_cohorts <- function(bbj, fg, ukbb) {
  dplyr::bind_cols(
    x = !is.na(bbj),
    y = !is.na(fg),
    z = !is.na(ukbb)
  ) %>%
    purrr::pmap_chr(function(x, y, z) {
      stringr::str_c(c("BBJ", "FG", "UKBB")[c(x, y, z)], collapse = ", ")
    })
}

################################################################################
# Table 1: Extremely population-enriched coding variants
################################################################################
df <- rgsutil::read_gsfile(fm_insights$get_analysis_path("af_enrichment.over5", "tsv"))

tbl1 <-
  dplyr::filter(
    df,
    data_type == "exomes" &
      consequence %in% c("pLoF", "Missense") &
      max_pip > 0.9 & enrichment_pseudo > 10
  ) %>%
  dplyr::mutate(
    cohort = ifelse(pop == "jpn", "BBJ", "FG"),
    most_severe = prettify_most_severe(most_severe),
    pip09_traits = prettify_array(pip09_traits),
    AF_ref = ifelse(is.infinite(enrichment), NA, AF_ref),
    AF_pop = sprintf("%.2g", AF_pop),
    AF_ref = sprintf("%.2g", AF_ref),
    enrichment = sprintf("%.1f", enrichment),
    max_pip = sprintf("%.2f", max_pip)
  )

tbl1 <- dplyr::bind_rows(
  dplyr::filter(tbl1, cohort == "BBJ"),
  dplyr::filter(tbl1, cohort == "FG")
)

dplyr::select(
  tbl1,
  cohort,
  variant,
  rsid,
  gene_most_severe,
  most_severe,
  AF_pop,
  AF_ref,
  enrichment,
  max_pip,
  pip09_traits
) %>%
  write.table(
    "./tables/Table1_enriched_coding.tsv",
    quote = F,
    row.names = F,
    sep = "\t"
  )


################################################################################
# STable: High-PIP pairs
################################################################################
df <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("pip09.vtc", "tsv"))

stbl <- dplyr::mutate(df,
  tissue_max_pip_coloc = prettify_array(tissue_max_pip_coloc)
) %>%
  dplyr::select(-susie.alpha)
write.table(
  stbl,
  "./tables/STable_high_pip_pairs.tsv",
  quote = F,
  row.names = F,
  sep = "\t"
)

################################################################################
# STable: Merged CS summary
################################################################################
df <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("csm.summary", "tsv"))

stbl <- dplyr::mutate(df,
  cohorts = prettify_array(cohorts),
  max_pip_variant = prettify_array(max_pip_variant)
)

write.table(
  stbl,
  "./tables/STable_CS_summary.tsv",
  quote = F,
  row.names = F,
  sep = "\t"
)

################################################################################
# STable: Coloc
################################################################################
df <- rgsutil::read_gsfile(
  fm_insights$get_merged_results_path("fm_only.coloc.max_pip_per_vtg.pip_coloc01", "tsv")
)

stbl <- dplyr::mutate(
  df,
  tissue_max_pip_coloc = prettify_array(tissue_max_pip_coloc),
  eqtl_study = prettify_array(eqtl_study),
  eqtl_symbol = prettify_array(eqtl_symbol),
  eqtl_pip = prettify_array(eqtl_pip)
)

write.table(
  stbl,
  "./tables/STable_coloc.tsv",
  quote = F,
  row.names = F,
  sep = "\t"
)


################################################################################
# STable: High-confidence putative causal coding/non-coding variants
################################################################################
df <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("shard_trait.pip09.any.pip01.every", "tsv"))
df.pip09 <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("shard_trait.pip09.every", "tsv"))

stbl <-
  dplyr::mutate(
    df,
    most_severe = prettify_most_severe(most_severe),
    tissue_max_pip_coloc = prettify_array(tissue_max_pip_coloc),
    cohorts = get_cohorts(pvalue.BBJ, pvalue.FG, pvalue.UKBB)
  )


stbl.eur_specific <- dplyr::filter(stbl, n_pop == 2 &
  !is.na(pvalue.FG) &
  !is.na(pvalue.UKBB))
stbl.coding <- dplyr::filter(stbl, consequence %in% c("pLoF", "Missense", "Synonymous"))
stbl.non_coding <- dplyr::filter(stbl, !(consequence %in% c("pLoF", "Missense", "Synonymous")))

nrow(df)
length(unique(df$variant))
nrow(stbl.eur_specific)
length(unique(stbl.eur_specific$variant))
nrow(stbl.coding)
length(unique(stbl.coding$variant))
nrow(stbl.non_coding)
length(unique(stbl.non_coding$variant))
nrow(df.pip09)
length(unique(df.pip09$variant))

dplyr::filter(stbl.non_coding, !is.na(max_pip_coloc)) %>%
  dplyr::distinct(variant) %>%
  nrow()


dplyr::select(
  stbl.coding,
  variant,
  rsid,
  trait,
  cohorts,
  gene_most_severe,
  most_severe,
  consequence,
  max_pip,
  clinvar,
  max_pip_coloc,
  tissue_max_pip_coloc,
  dplyr::starts_with("af."),
  dplyr::starts_with("beta_marginal."),
  dplyr::starts_with("se_marginal."),
  dplyr::starts_with("pvalue."),
  dplyr::starts_with("pip."),
  dplyr::starts_with("susie.beta_posterior."),
  dplyr::starts_with("susie.sd_posterior.")
) %>%
  write.table(
    "./tables/STable_high_confidence_coding.tsv",
    quote = F,
    row.names = F,
    sep = "\t"
  )

dplyr::select(
  stbl.non_coding,
  variant,
  rsid,
  trait,
  cohorts,
  gene_most_severe,
  most_severe,
  consequence,
  max_pip,
  clinvar,
  max_pip_coloc,
  tissue_max_pip_coloc,
  Ancient_Sequence_Age_Human_Enhancer:CA_H3K27ac_Ulirsch
) %>%
  write.table(
    "./tables/STable_high_confidence_non_coding.tsv",
    quote = F,
    row.names = F,
    sep = "\t"
  )

################################################################################
# STable: High-confidence intergeneic variants
################################################################################
df.intergenic <- read.table("./data/intergenic_variants.dist.txt", T) %>%
  dplyr::transmute(
    variant = variant,
    closest_gene = gene,
    distance = score
  )

df <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("shard_trait.pip09.any.pip01.every", "tsv")) %>%
  dplyr::filter(variant %in% df.intergenic$variant) %>%
  dplyr::group_by(variant) %>%
  dplyr::summarize(
    rsid = rsid[1],
    pip09_traits = stringr::str_c(trait, collapse = ", ")
  )

stbl <- dplyr::inner_join(df.intergenic, df)
dplyr::select(
  stbl,
  variant,
  rsid,
  closest_gene,
  distance,
  pip09_traits
) %>%
  write.table(
    "./tables/STable_high_confidence_intergenic.tsv",
    quote = F,
    row.names = F,
    sep = "\t"
  )

################################################################################
# STable: Extremely population-enriched non-coding variants
################################################################################
df <- rgsutil::read_gsfile(fm_insights$get_analysis_path("af_enrichment.over5", "tsv"))

stbl <-
  dplyr::filter(
    df,
    data_type == "genomes" &
      max_pip > 0.9 & enrichment_pseudo > 10
  ) %>%
  dplyr::mutate(
    cohort = ifelse(pop == "jpn", "BBJ", "FG"),
    most_severe = prettify_most_severe(most_severe),
    pip09_traits = prettify_array(pip09_traits),
    AF_ref = ifelse(is.infinite(enrichment), NA, AF_ref)
  )

dplyr::select(
  stbl,
  cohort,
  variant,
  rsid,
  gene_most_severe,
  most_severe,
  AF_pop,
  AF_ref,
  enrichment,
  max_pip,
  pip09_traits
) %>%
  write.table(
    "./tables/STable_enriched_non_coding.tsv",
    quote = F,
    row.names = F,
    sep = "\t"
  )