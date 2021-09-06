library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

#########################################################
# Fig. 4: allelic series summary
#########################################################
df.traits <- read_trait_summary()
shared.traits <- df.traits %>%
  dplyr::group_by(trait) %>%
  dplyr::filter(length(unique(cohort)) > 1) %>%
  .$trait %>%
  unique()

# filter to in-CS
df.in_cs <-
  rgsutil::read_gsfile(fm_insights$get_merged_results_path("fm_only.in_cs.consequence", "tsv.bgz")) %>%
  dplyr::mutate(variant = locusviz::variant_str2(locus, alleles)) %>%
  dplyr::mutate(locusviz::parse_variant(variant))

#########################################################
# a) # cohorts x # nonsyn variants per gene
#########################################################
df.nonsyn <-
  dplyr::filter(df.in_cs, consequence %in% c("pLoF", "Missense") &
    pip > 0.1) %>%
  dplyr::mutate(variant = locusviz::variant_str2(locus, alleles))

dplyr::group_by(df.nonsyn, variant) %>%
  dplyr::summarize(
    rsid = rsid[1],
    consequence = consequence[1],
    most_severe = most_severe[1],
    gene_most_severe = gene_most_severe[1],
    max_pip = max(pip),
    pip01_traits = stringr::str_c(trait, collapse = ", ")
  ) %>%
  write.table(
    "./tables/STable_nonsynonymous_coding.txt",
    quote = F,
    row.names = F,
    sep = "\t"
  )

count.per_gene <-
  dplyr::group_by(df.nonsyn, gene_most_severe) %>%
  dplyr::summarize(
    n_variants = length(unique(variant)),
    n_cohorts = length(unique(cohort))
  ) %$%
  table(n_variants, n_cohorts) %>%
  tibble::as_tibble() %>%
  dplyr::rename(count = n) %>%
  dplyr::mutate(
    n_variants = factor(as.numeric(n_variants)),
    frac = count / sum(count)
  )


# No. unique pLoF/missense variants (PIP > 0.1)
print(length(unique(df.nonsyn$variant)))
# mapped on to # genes
print(length(unique(df.nonsyn$gene_most_severe)))

# > 1 variants on the same gene
dplyr::mutate(count.per_gene, n_variants = as.numeric(n_variants)) %>%
  dplyr::filter(n_variants > 1) %>%
  dplyr::pull(count) %>%
  sum() %>%
  print()

# > 1 variants and > 1 cohorts
dplyr::mutate(count.per_gene, n_variants = as.numeric(n_variants)) %>%
  dplyr::filter(n_variants > 1 & n_cohorts > 1) %>%
  dplyr::pull(count) %>%
  sum() %>%
  print()


plt.a <-
  dplyr::mutate(count.per_gene, fill = cut(count, c(-Inf, 0, 3, 10, 50, 100, 1000))) %>%
  ggplot(aes(n_cohorts, n_variants)) +
  geom_tile(aes(fill = fill), color = "white") +
  geom_text(aes(
    label = sprintf(
      "%d\n(%s)",
      count,
      scales::percent_format(accuracy = 0.1, drop0trailing = TRUE)(frac)
    ),
    color = fill
  ), size = 2) +
  geom_rect(
    aes(
      xmin = 1.5,
      xmax = 3.5,
      ymin = 3.5,
      ymax = 8.5
    ),
    fill = NA,
    color = "grey50",
    linetype = "dotted",
    size = 0.5
  ) +
  my_theme +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm")
  ) +
  labs(x = "# cohorts", y = "# pLoF/missense variants per gene\n(Best PIP > 0.1)") +
  scale_fill_manual(
    values = c(
      "grey90",
      BuenColors::jdb_palette("brewer_blue")[8],
      BuenColors::jdb_palette("brewer_celsius")[c(3, 5, 7)],
      BuenColors::jdb_palette("brewer_red")[7]
    )
  ) +
  scale_color_manual(values = c("black", "white", rep("black", 3), "white"))
plt.a

#########################################################
# b) # nonsyn variants per gene
#########################################################
gene_count_threshold <- 3
df.gene_count <-
  dplyr::select(df.nonsyn, variant, gene_most_severe, cohort) %>%
  dplyr::distinct() %>%
  dplyr::group_split(gene_most_severe, variant) %>%
  purrr::map_dfr(~ {
    cohort <- unique(.$cohort)
    x <- c()
    if (length(cohort) == 3) {
      x <- "All"
    } else {
      for (pop in cohorts) {
        if (pop %in% cohort) {
          x <- c(x, pop)
        }
      }
      x <- paste(x, collapse = ",")
    }

    tibble::tibble(
      gene_most_severe = .$gene_most_severe[1],
      cohort = x
    )
  }) %>%
  dplyr::group_by(gene_most_severe) %>%
  dplyr::filter(n() > gene_count_threshold &
    length(unique(cohort)) > 1) %>%
  dplyr::group_by(gene_most_severe, cohort) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::mutate(cohort = factor(
    cohort,
    levels = c(cohorts, "BBJ,UKBB", "BBJ,FG", "FG,UKBB", "All")
  ))

gene_order <-
  dplyr::group_by(df.gene_count, gene_most_severe) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::arrange(count) %>%
  .$gene_most_severe

df.gene_count <-
  dplyr::mutate(df.gene_count,
    gene_most_severe = factor(gene_most_severe, levels = gene_order)
  )

plt.b <-
  ggplot(df.gene_count, aes(count, gene_most_severe)) +
  geom_col(aes(fill = cohort),
    position = position_stack(reverse = T)
  ) +
  my_theme +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 14, 2)) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "italic"),
    legend.position = c(1, 0.05),
    legend.justification = c(1, 0),
    plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm")
  ) +
  labs(x = "# pLoF/missense variants\n(Best PIP > 0.1)", fill = "Cohort") +
  scale_fill_manual(
    values = c(
      cohort_colors,
      `BBJ,UKBB` = PNWColors::pnw_palette("Bay")[3],
      `FG,UKBB` = BuenColors::jdb_palette("samba_night")[2],
      All = "#777788"
    )
  )
plt.b

#########################################################
# c) allelic series summary
#########################################################
df.in_cs.coding <-
  df.in_cs %>%
  dplyr::filter(consequence %in% c("pLoF", "Missense") &
    pip > 0.1) %>%
  dplyr::mutate(to_idx = dplyr::row_number())

gr.in_cs.coding <-
  GenomicRanges::makeGRangesFromDataFrame(
    df.in_cs.coding,
    seqnames.field = "chromosome",
    start.field =  "position",
    end.field = "position"
  )


window <- 50000
coding <- c("pLoF", "Missense", "Synonymous")

# Not
df.in_cs.non_coding <-
  dplyr::filter(df.in_cs, susie.cs_id > 0) %>%
  # dplyr::filter(!(most_severe %in% c("intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant"))) %>%
  dplyr::group_by(cohort, trait, region, susie.cs_id, chromosome) %>%
  dplyr::filter(!any(consequence %in% coding) &
    max(pip, na.rm = T) > 0.1) %>%
  dplyr::summarize(
    position = position[which.max(pip)],
    start = position - window,
    end = position + window,
    variant = variant[which.max(pip)],
    most_severe = stringr::str_c(unique(most_severe), collapse = ","),
    consequence = stringr::str_c(unique(consequence), collapse = ",")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(from_idx = dplyr::row_number())

gr.in_cs.non_coding <-
  GenomicRanges::makeGRangesFromDataFrame(
    df.in_cs.non_coding,
    seqnames.field = "chromosome",
    start.field =  "start",
    end.field = "end"
  )

hits <-
  GenomicRanges::findOverlaps(gr.in_cs.non_coding, gr.in_cs.coding) %>%
  tibble::as_tibble() %>%
  dplyr::rename(from_idx = queryHits, to_idx = subjectHits) %>%
  dplyr::left_join(
    dplyr::select(
      df.in_cs.non_coding,
      cohort,
      trait,
      region,
      susie.cs_id,
      consequence,
      most_severe,
      from_idx,
      variant
    ) %>%
      dplyr::rename(
        cohort1 = cohort,
        trait1 = trait,
        region1 = region,
        susie.cs_id1 = susie.cs_id,
        consequence1 = consequence,
        most_severe1 = most_severe,
        variant1 = variant
      ),
    by = "from_idx"
  ) %>%
  dplyr::left_join(
    dplyr::select(
      df.in_cs.coding,
      cohort,
      trait,
      region,
      pip,
      susie.cs_id,
      consequence,
      most_severe,
      gene_most_severe,
      to_idx,
      variant
    ) %>%
      dplyr::rename(
        cohort2 = cohort,
        trait2 = trait,
        region2 = region,
        susie.cs_id2 = susie.cs_id,
        consequence2 = consequence,
        most_severe2 = most_severe,
        variant2 = variant
      ),
    by = "to_idx"
  ) %>%
  dplyr::filter(trait1 == trait2) %>%
  dplyr::mutate(distance = abs(
    locusviz::parse_variant(variant1)$position - locusviz::parse_variant(variant2)$position
  ))

count.total.coding <-
  dplyr::distinct(df.in_cs.coding, cohort, trait, gene_most_severe) %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarize(count = n())

count.any.total.coding <-
  dplyr::distinct(df.in_cs.coding, trait, gene_most_severe) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::mutate(cohort = "Any")

total_vec <- tibble::deframe(count.total.coding)

xpop.coding <-
  df.in_cs %>%
  dplyr::filter(consequence %in% coding) %>%
  dplyr::mutate(variant = locusviz::variant_str2(locus, alleles)) %>%
  dplyr::group_split(trait, gene_most_severe) %>%
  purrr::map_dfr(function(data) {
    if (all(!(data$pip > 0.1 &
      data$consequence %in% c("pLoF", "Missense")), na.rm = T)) {
      return(NULL)
    }
    purrr::pmap_dfr(tidyr::crossing(cohort1 = cohorts, cohort2 = cohorts), function(cohort1, cohort2) {
      # nonsyn coding variants with PIP > 0.1 in the first pop
      x <-
        dplyr::filter(
          data,
          cohort == cohort1 &
            pip > 0.1 &
            consequence %in% c("pLoF", "Missense")
        )
      # coding CS with the same variant in the second pop
      cs_ids <-
        dplyr::filter(data, cohort == cohort2 &
          variant %in% x$variant) %>%
        .$susie.cs_id %>%
        unique()
      y <-
        dplyr::filter(data, cohort == cohort2 &
          !(susie.cs_id %in% cs_ids) & pip > 0.1)
      y_nonsyn <-
        dplyr::filter(y, consequence %in% c("pLoF", "Missense"))
      y_syn <-
        dplyr::filter(
          y,
          consequence %in% c("Synonymous") &
            !(susie.cs_id %in% unique(y_nonsyn$susie.cs_id))
        )

      count_nonsyn <-
        ifelse(cohort1 == cohort2,
          (nrow(x) + nrow(y_nonsyn)) > 1,
          nrow(x) > 0 & nrow(y_nonsyn) > 0
        )
      count_syn <- nrow(x) > 0 & nrow(y_syn) > 0
      tibble::tibble(
        trait = data$trait[1],
        gene_most_severe = data$gene_most_severe[1],
        cohort1 = cohort1,
        cohort2 = cohort2,
        consequence = c("Nonsynonymous", "Synonymous"),
        count = c(count_nonsyn, count_syn),
        x = stringr::str_c(x$variant, collapse = ","),
        y = purrr::map_chr(list(y_nonsyn, y_syn), function(z) {
          stringr::str_c(z$variant, collapse = ",")
        })
      )
    })
  })

count.coding <-
  dplyr::group_by(xpop.coding, cohort1, cohort2, consequence) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::mutate(frac = count / total_vec[cohort1])

count.non_coding <-
  dplyr::group_by(hits, cohort1, trait1, gene_most_severe, cohort2) %>%
  dplyr::summarize(count = n() > 0) %>%
  dplyr::group_by(cohort1, cohort2) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::mutate(frac = count / total_vec[cohort1])

count.any.coding <-
  dplyr::group_by(xpop.coding, trait, gene_most_severe, consequence) %>%
  dplyr::summarize(count = any(count)) %>%
  dplyr::group_by(consequence) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::mutate(
    cohort1 = "Any",
    cohort2 = "Any",
    frac = count / count.any.total.coding$count[1]
  )

count.any.non_coding <-
  dplyr::group_by(hits, trait1, gene_most_severe) %>%
  dplyr::summarize(count = n() > 0) %>%
  dplyr::ungroup() %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::mutate(
    cohort1 = "Any",
    cohort2 = "Any",
    frac = count / count.any.total.coding$count[1]
  )

allelic_series.coding <-
  dplyr::filter(xpop.coding, consequence == "Nonsynonymous") %>%
  dplyr::group_by(trait, gene_most_severe) %>%
  dplyr::summarize(count = any(count)) %>%
  dplyr::filter(count) %>%
  dplyr::select(-count) %>%
  dplyr::mutate(id = stringr::str_c(trait, gene_most_severe, sep = ":"))

# STables
dplyr::filter(xpop.coding, consequence == "Nonsynonymous" &
  count) %>%
  dplyr::group_by(trait, gene_most_severe) %>%
  dplyr::mutate(variants = list(
    stringr::str_split(stringr::str_c(c(x, y), collapse = ","), ",") %>%
      purrr::pluck(1) %>%
      purrr::keep(function(x) {
        x != ""
      })
  )) %>%
  dplyr::summarize(
    n_uniq_variants = length(unique(variants[[1]])),
    cohorts = stringr::str_c(sort(unique(c(
      cohort1, cohort2
    ))), collapse = ", "),
    is_xpop = any(cohort1 != cohort2),
    variants = stringr::str_c(sort(unique(variants[[1]])), collapse = ", ")
  ) %>%
  dplyr::rename(gene = gene_most_severe) %>%
  dplyr::select(
    gene,
    trait,
    n_uniq_variants,
    cohorts,
    is_xpop,
    variants
  ) %>%
  write.table(
    "./tables/STable_allelic_series_coding.tsv",
    quote = F,
    row.names = F,
    sep = "\t"
  )

dplyr::group_by(hits, trait1, gene_most_severe) %>%
  dplyr::summarize(
    cohorts = stringr::str_c(sort(unique(c(
      cohort1, cohort2
    ))), collapse = ", "),
    is_xpop = any(cohort1 != cohort2),
    coding_non_coding_only = stringr::str_c(trait1, gene_most_severe, sep = ":")[1] %in% allelic_series.coding$id,
    n_uniq_nonsyn_variants = length(unique(variant2)),
    n_non_coding_cs = length(unique(
      stringr::str_c(cohort1, region1, susie.cs_id1)
    )),
    closest_distance = min(distance),
    closest_pair_coding = variant2[which.min(distance)],
    closest_pair_non_coding = variant1[which.min(distance)]
  ) %>%
  dplyr::rename(trait = trait1, gene = gene_most_severe) %>%
  dplyr::select(
    gene,
    trait,
    n_uniq_nonsyn_variants,
    n_non_coding_cs,
    cohorts,
    is_xpop,
    coding_non_coding_only,
    closest_distance,
    closest_pair_coding,
    closest_pair_non_coding
  ) %>%
  write.table(
    "./tables/STable_allelic_series_coding_non_coding.tsv",
    quote = F,
    row.names = F,
    sep = "\t"
  )



# xpop allelic series coding
dplyr::mutate(xpop.coding, is_xpop = cohort1 != cohort2) %>%
  dplyr::group_by(trait, gene_most_severe, consequence) %>%
  dplyr::summarize(
    total = any(count),
    any_xpop = any(count[is_xpop]),
    xpop_only = any_xpop & !any(count[!is_xpop])
  ) %>%
  dplyr::group_by(consequence) %>%
  dplyr::summarize(
    uniq_genes = length(unique(gene_most_severe[total])),
    uniq_genes_xpop = length(unique(gene_most_severe[any_xpop])),
    uniq_genes_xpop_only = length(unique(gene_most_severe[xpop_only])),
    total = sum(total),
    any_xpop = sum(any_xpop),
    xpop_only = sum(xpop_only)
  )
# xpop allelic series non-coding
dplyr::mutate(hits, is_xpop = cohort1 != cohort2) %>%
  dplyr::group_by(trait1, gene_most_severe) %>%
  dplyr::summarize(
    total = TRUE,
    any_xpop = sum(is_xpop) > 0,
    xpop_only = any_xpop & sum(!is_xpop) == 0
  ) %>%
  dplyr::ungroup() %>%
  dplyr::summarize(
    uniq_genes = length(unique(gene_most_severe[total])),
    uniq_genes_xpop = length(unique(gene_most_severe[any_xpop])),
    uniq_genes_xpop_only = length(unique(gene_most_severe[xpop_only])),
    total = sum(total),
    any_xpop = sum(any_xpop),
    xpop_only = sum(xpop_only)
  )

# gene-trait pairs only found in non-coding
dplyr::anti_join(
  dplyr::mutate(hits, is_xpop = cohort1 != cohort2) %>%
    dplyr::group_by(trait1, gene_most_severe) %>%
    dplyr::summarize(
      total = TRUE,
      any_xpop = sum(is_xpop) > 0,
      xpop_only = any_xpop & sum(!is_xpop) == 0
    ),
  dplyr::select(allelic_series.coding, -id)
) %>%
  dplyr::ungroup() %>%
  dplyr::summarize(
    uniq_genes = length(unique(gene_most_severe[total])),
    uniq_genes_xpop = length(unique(gene_most_severe[any_xpop])),
    uniq_genes_xpop_only = length(unique(gene_most_severe[xpop_only])),
    total = sum(total),
    any_xpop = sum(any_xpop),
    xpop_only = sum(xpop_only)
  )

a <-
  dplyr::mutate(hits, is_xpop = cohort1 != cohort2) %>%
  dplyr::group_by(trait1, gene_most_severe) %>%
  dplyr::filter(sum(is_xpop) & sum(!is_xpop) == 0) %>%
  dplyr::ungroup()

sort(table(a$gene_most_severe))

p1 <-
  ggplot(count.total.coding, aes(factor(1), factor(cohort, levels = rev(cohorts)))) +
  geom_tile(aes(fill = cohort), color = "white") +
  geom_text(aes(label = sprintf("%d\n(100%%)", count)), size = 2, color = "white") +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(vjust = 0),
    legend.position = "none"
  ) +
  labs(x = "Total\n# gene-trait", y = "pLoF/missense variant (PIP > 0.1)") +
  scale_fill_manual(values = cohort_colors)

p2.nonsyn <-
  dplyr::filter(count.coding, consequence == "Nonsynonymous") %>%
  ggplot(aes(cohort2, factor(cohort1, levels = rev(cohorts)))) +
  geom_tile(aes(fill = ifelse(frac == 0, NA, frac)), color = "white") +
  geom_text(aes(label = sprintf(
    "%d\n(%s)",
    count,
    scales::percent_format(accuracy = 0.1, drop0trailing = TRUE)(frac)
  )), size = 2) +
  my_theme +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "95% CS with pLoF/\nmissense (PIP > 0.1)", y = "pLoF/missense variant (PIP > 0.1)") +
  scale_fill_stepsn(
    colors = BuenColors::jdb_palette("brewer_red", 7),
    na.value = "grey90"
  )

p2.nc <-
  ggplot(count.non_coding, aes(cohort2, factor(cohort1, levels = rev(cohorts)))) +
  geom_tile(aes(fill = ifelse(frac == 0, NA, frac)), color = "white") +
  geom_text(aes(label = sprintf(
    "%d\n(%s)",
    count,
    scales::percent_format(accuracy = 0.1, drop0trailing = TRUE)(frac)
  )), size = 2) +
  my_theme +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Non-coding CS\nwithin 100 kb (PIP > 0.1)") +
  scale_fill_stepsn(
    colors = BuenColors::jdb_palette("brewer_blue", 7),
    na.value = "grey90"
  )


q1 <- ggplot(count.any.total.coding, aes(factor(""), cohort)) +
  geom_tile(aes(fill = cohort), color = "white", fill = "#777788") +
  geom_text(aes(label = sprintf(
    "%s\n(100%%)", scales::comma_format()(count)
  )), size = 2, color = "white") +
  my_theme +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Total\n# gene-trait", y = "pLoF/missense variant (PIP > 0.1)")

q2.nonsyn <-
  dplyr::filter(count.any.coding, consequence == "Nonsynonymous") %>%
  ggplot(aes(cohort2, factor(cohort1, levels = rev(cohorts)))) +
  geom_tile(aes(fill = factor(frac)), color = "white") +
  geom_text(aes(label = sprintf(
    "%d\n(%s)",
    count,
    scales::percent_format(accuracy = 0.1, drop0trailing = TRUE)(frac)
  )), size = 2) +
  my_theme +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "95% CS with pLoF/missense\n(PIP > 0.1)", y = "pLoF/missense variant (PIP > 0.1)") +
  scale_fill_manual(values = BuenColors::jdb_palette("brewer_red", 7)[4]) +
  scale_x_discrete(expand = expansion(add = 0.53))

q2.nc <-
  ggplot(count.any.non_coding, aes(cohort2, factor(cohort1, levels = rev(cohorts)))) +
  geom_tile(aes(fill = factor(frac)), color = "white") +
  geom_text(aes(label = sprintf(
    "%d\n(%s)",
    count,
    scales::percent_format(accuracy = 0.1, drop0trailing = TRUE)(frac)
  )), size = 2) +
  my_theme +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Non-coding CS\nwithin 100 kb (PIP > 0.1)") +
  scale_fill_manual(values = BuenColors::jdb_palette("brewer_blue", 7)[4]) +
  scale_x_discrete(expand = expansion(add = 0.53))

plt.c <- purrr::reduce(
  list(
    q1 + labs(tag = "c"),
    q2.nonsyn,
    q2.nc,
    p1,
    p2.nonsyn,
    p2.nc
  ),
  `+`
) + patchwork::plot_layout(
  nrow = 2,
  heights = c(1, 3),
  widths = c(1, rep(3, 2))
)


plt <-
  (plt.a + labs(tag = "a") + plt.b + labs(tag = "b") + patchwork::plot_layout(widths = c(4, 5))) / plt.c
plt

cowplot::save_plot(
  "figures/Fig4/Fig4_allelic_series_summary.pdf",
  plt,
  base_height = 5.0,
  base_width = 3.5,
  device = cairo_pdf
)
