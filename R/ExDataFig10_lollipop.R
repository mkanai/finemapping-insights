library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

# filter to in-CS
df.in_cs <-
  rgsutil::read_gsfile(fm_insights$get_merged_results_path("fm_only.in_cs.consequence", "tsv.bgz")) %>%
  dplyr::mutate(variant = locusviz::variant_str2(locus, alleles)) %>%
  dplyr::mutate(locusviz::parse_variant(variant))

df.nonsyn =
  dplyr::filter(df.in_cs,
                consequence %in% c("pLoF", "Missense") &
                  pip > 0.1)

#########################################################
# ExData Fig 10a. constraint
#########################################################
contraint = rgsutil::read_gsfile(
  'gs://gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz'
)

df.contraint =
  dplyr::rename(df.nonsyn, gene = gene_most_severe) %>%
  dplyr::group_by(gene, consequence) %>%
  dplyr::summarize(
    max_pip = max(pip),
    moset_severe_consequence = factor(
      dplyr::case_when(
        consequence == "pLoF" ~ "pLoF",
        consequence == "Missense" ~ "Missense",
        TRUE ~ NA_character_
      ),
      levels = annot_levels
    )
  ) %>%
  dplyr::inner_join(contraint %>% dplyr::select(
    gene,
    pLI,
    oe_lof_upper,
    constraint_flag,
    starts_with("oe_lof_upper")
  )) %>%
  dplyr::filter(!is.na(oe_lof_upper_bin)) %>%
  dplyr::mutate(finemapped = !is.na(max_pip))

label_function = function(x) {
  y = 10 * x + 5
  ifelse(y %% 20 == 0, paste0(y, '%'), "")
}

p1 =
  ggplot(df.contraint, aes(oe_lof_upper_bin)) +
  geom_bar(aes(fill = moset_severe_consequence)) +
  geom_text(
    aes(label = ..count.., group = moset_severe_consequence),
    stat = "count",
    position = position_stack(0.5, reverse = FALSE),
    size = 2
  ) +
  scale_fill_manual(values = annot_colors[c("pLoF", "Missense")]) +
  my_theme +
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.title = element_blank()) +
  scale_x_continuous(labels = label_function, breaks=seq(-0.5, 9.5, 1), limits=c(-0.5, 9.5)) +
  labs(x = "LOEUF decile",
       y = "# genes with fine-mapped\nnonsynonymous variants (PIP > 0.1)",
       tag = "a")
p1

#########################################################
# ExData Fig 10b. lollipop ABCG2
#########################################################
hgvsp = rgsutil::read_gsfile(fm_insights$get_analysis_path("hgvsp", "tsv.bgz"))
clinvar = rgsutil::read_gsfile('gs://xfinemap/annotation/clinvar/variant_summary_20200330.formatted.txt.bgz')

df.in_cs <-
  rgsutil::read_gsfile(fm_insights$get_merged_results_path("fm_only.in_cs.consequence", "tsv.bgz")) %>%
  dplyr::mutate(locusviz::parse_variant(locus),
                variant = locusviz::variant_str2(locus, alleles)) %>%
  dplyr::left_join(hgvsp)

p_APOB =
  dplyr::mutate(df.in_cs,
                label = ifelse(!is.na(hgvsp), hgvsp, rsid)) %>%
  dplyr::filter(pip > 0.1 &
                  gene_most_severe == "APOB") %>%
  dplyr::group_by(cohort, trait, susie.cs_id) %>%
  dplyr::filter(consequence %in% c("pLoF", "Missense", "Synonymous")) %>%
  dplyr::ungroup() %>%
  locusviz:::plot_lollipop(
    gene_symbol = "APOB",
    clinvar = head(clinvar, 0),
    point_colors = cohort_colors,
    point_shapes = c(BBJ = 18, FG = 18, UKBB = 18),
    color_by_cohort = TRUE,
    omit_spacer = TRUE,
    gene_col = 'grey90'
  )

################################################################

p_ABCG2 =
  dplyr::mutate(df.in_cs,
                label = ifelse(!is.na(hgvsp), hgvsp, rsid)) %>%
  dplyr::filter(pip > 0.01 &
                  gene_most_severe == "ABCG2" &
                  trait %in% c("Gout", "UA")) %>%
  dplyr::group_by(cohort, trait, susie.cs_id) %>%
  dplyr::filter(consequence %in% c("pLoF", "Missense", "Synonymous")) %>%
  dplyr::ungroup() %>%
  locusviz:::plot_lollipop(
    gene_symbol = "ABCG2",
    clinvar = clinvar,
    point_colors = cohort_colors,
    point_shapes = c(BBJ = 18, FG = 18, UKBB = 18),
    color_by_cohort = TRUE,
    omit_spacer = TRUE,
    gene_col = 'grey90'
  )

################################################################

p_EPX =
  dplyr::mutate(df.in_cs,
                label = ifelse(!is.na(hgvsp), hgvsp, rsid)) %>%
  dplyr::filter(pip > 0.1 & trait == "Eosino") %>%
  dplyr::group_by(cohort, trait, region, susie.cs_id) %>%
  dplyr::filter(consequence %in% c("pLoF", "Missense", "Synonymous") |
                  (rsid == "rs536070968")) %>%
  dplyr::ungroup() %>%
  locusviz:::plot_lollipop(
    gene_symbol = c("EPX", "LPO"),
    clinvar = clinvar,
    point_colors = cohort_colors,
    point_shapes = c(BBJ = 18, FG = 18, UKBB = 18),
    color_by_cohort = TRUE,
    omit_spacer = TRUE,
    gene_col = 'grey90'
  )


layout = "
AACCCCC
BBBBBBB
DDDDDDD
"

p_APOB[[1]] = p_APOB[[1]] + labs(tag = "b")
p_ABCG2[[1]] = p_ABCG2[[1]] + labs(tag = "c")
p_EPX[[1]] = p_EPX[[1]] + labs(tag = "d")

plt = p1 + p_APOB + p_ABCG2 + p_EPX +
  patchwork::plot_layout(design = layout, heights = c(1.5, 3, 2))
plt

cowplot::save_plot(
  "figures/ExDataFig10_lollipop.pdf",
  plt,
  base_height = 9.6,
  base_width = 7.2,
  device = cairo_pdf
)
