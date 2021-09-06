import argparse
import hail as hl

from fm_insights.resources import (
    get_analysis_path,
    get_baseline_annotations,
    get_gtex_sqtl_signifpairs_path,
    get_merged_results_path,
    get_n_pop_dict,
    get_results_path,
    get_ref_freq_info_path,
    get_spliceai_path,
    POPS,
)
from fm_insights.utils import register_log, checkpoint_tmp
from itertools import product


def export_flatfile(ht, outname, extra_fields=[]):
    ht = ht.key_by()
    ht = ht.select(
        "variant",
        "rsid",
        "trait",
        "n_pop",
        "max_pip",
        "consequence",
        "most_severe",
        "gene_most_severe",
        "clinvar",
        "max_pip_coloc",
        "tissue_max_pip_coloc",
        "beta_meta",
        "se_meta",
        "pvalue_meta",
        "pvalue_het",
        "af",
        "beta_marginal",
        "se_marginal",
        "pvalue",
        "pip",
        "susie",
        "finemap",
        *get_baseline_annotations(basename=True),
        *extra_fields,
    )

    ht.flatten().export(outname)


def export_pip09_variants(overwrite: bool = False):
    ht = hl.read_table(get_merged_results_path("fm_only.in_cs.consequence"))
    ht = ht.filter(ht.pip > 0.9)
    ht = ht.key_by("locus", "alleles", "trait", "cohort")

    ht2 = hl.read_table(get_merged_results_path("fm_only.coloc.max_pip_per_vtc"))
    ht2 = ht2.filter(ht2.max_pip_coloc > 0.1)
    ht2 = ht2.key_by("locus", "alleles", "trait", "cohort")
    ht2 = ht2.drop("pip", "eqtl")
    ht = ht.join(ht2, "left")

    ht = ht.checkpoint(get_merged_results_path("pip09.vtc"), overwrite=overwrite)
    ht = ht.key_by()
    ht = ht.select(
        "variant",
        "rsid",
        "trait",
        "cohort",
        "consequence",
        "most_severe",
        "gene_most_severe",
        "clinvar",
        "max_pip_coloc",
        "tissue_max_pip_coloc",
        "af",
        "beta_marginal",
        "se_marginal",
        "pvalue",
        "pip",
        "susie",
        "finemap",
        *get_baseline_annotations(basename=True),
    )

    ht.flatten().export(get_merged_results_path("pip09.vtc", "tsv"))


def export_csm_summary(overwrite: bool = False):
    ht = hl.read_table(get_merged_results_path("fm_only"))
    ht = ht.key_by("locus", "alleles", "trait", "cohort")

    ht_csm = hl.import_table(get_merged_results_path("fm_only.csm_id", "tsv.bgz"), impute=True, min_partitions=200)
    ht_csm = ht_csm.annotate(**hl.parse_variant(ht_csm.variant))
    ht_csm = ht_csm.key_by("locus", "alleles", "trait", "cohort")
    ht_csm = ht_csm.select("csm_id")

    ht = ht.join(ht_csm, "inner")
    ht = ht.group_by("trait", "csm_id").aggregate(
        cohorts=hl.agg.collect_as_set(ht.cohort),
        max_pip=hl.agg.max(ht.susie.pip),
        variants=hl.agg.collect_as_set(hl.struct(cohort=ht.cohort, variant=ht.variant, pip=ht.susie.pip)),
    )
    ht = ht.annotate(
        max_pip_variant=ht.variants.filter(lambda x: ht.max_pip == x.pip).variant, n_variants=hl.len(ht.variants)
    )
    ht = ht.checkpoint(get_merged_results_path("csm.summary"), overwrite=overwrite)
    ht.drop("variants").export(get_merged_results_path("csm.summary", "tsv"))


def export_coloc_vtg(overwrite: bool = False):
    ht = hl.read_table(get_merged_results_path("fm_only.coloc.max_pip_per_vtg"))
    ht = ht.filter(ht.max_pip_coloc > 0.1)
    ht = ht.annotate(max_pip_coloc_eqtl=ht.eqtl.filter(lambda x: ht.max_pip_coloc == ht.max_pip * x.pip))
    ht = ht.annotate(
        tissue_max_pip_coloc=ht.max_pip_coloc_eqtl.tissue,
        eqtl_study=ht.max_pip_coloc_eqtl.study,
        eqtl_symbol=ht.max_pip_coloc_eqtl.symbol,
        eqtl_pip=ht.max_pip_coloc_eqtl.pip,
    )

    ht_csq = hl.read_table(get_merged_results_path(f"fm_only.max_pip_per_var.consequence"))
    ht_csq = ht_csq.select("rsid", "consequence", "most_severe", "gene_most_severe")
    ht = ht.annotate(variant=hl.variant_str(ht.locus, ht.alleles))
    ht = ht.key_by("locus", "alleles").join(ht_csq, "left")
    ht = ht.key_by().select(
        "variant",
        "rsid",
        "trait",
        "consequence",
        "most_severe",
        "gene_most_severe",
        "max_pip",
        "max_pip_coloc",
        "tissue_max_pip_coloc",
        "eqtl_study",
        "eqtl_symbol",
        "eqtl_pip",
    )

    ht.export(get_merged_results_path("fm_only.coloc.max_pip_per_vtg.pip_coloc01", "tsv"))


def export_pip09_in_every_pop(overwrite: bool = False):
    n_pops = get_n_pop_dict()
    ht = hl.read_table(get_merged_results_path("shard_trait.pip001"))
    # PIP > 0.1 in multiple populations
    ht = ht.filter(
        (hl.sum(hl.map(lambda x: x > 0.1, ht.values.pip)) > 1)
        & (hl.sum(hl.map(lambda x: hl.int32(hl.is_defined(x)), ht.values.pvalue)) > 1)
    )

    ht2 = hl.read_table(get_merged_results_path("shard_trait.pip001.pop"))
    ht = ht.select("values").join(ht2, "left")
    ht = ht.annotate(max_pip=hl.max(ht.values.pip))

    ht2 = hl.read_table(get_merged_results_path("fm_only.coloc.max_pip_per_vt"))
    ht2 = ht2.filter(ht2.max_pip_coloc > 0.1)
    ht2 = ht2.drop("max_pip", "eqtl")
    ht = ht.join(ht2, "left")

    ht2 = hl.read_table(get_merged_results_path("fm_only.max_pip_per_var.consequence"))
    ht2 = ht2.drop("max_pip")
    ht = ht.key_by("locus", "alleles").join(ht2, "left")

    # PIP > 0.1 in every population analyzed
    ht2 = ht.drop("values").checkpoint(get_merged_results_path("shard_trait.pip01.every"), overwrite=overwrite)
    export_flatfile(ht2, get_merged_results_path("shard_trait.pip01.every", "tsv"))

    # PIP > 0.9 in one pop, PIP > 0.1 in other populations
    ht2 = (
        ht.filter(hl.any(lambda x: x > 0.9, ht.values.pip))
        .drop("values")
        .checkpoint(get_merged_results_path("shard_trait.pip09.any.pip01.every"), overwrite=overwrite)
    )
    export_flatfile(ht2, get_merged_results_path("shard_trait.pip09.any.pip01.every", "tsv"))

    # PIP > 0.9 in every pop
    ht2 = (
        ht.filter(hl.all(lambda x: x > 0.9, ht.values.pip))
        .drop("values")
        .checkpoint(get_merged_results_path("shard_trait.pip09.every"), overwrite=overwrite)
    )
    export_flatfile(ht2, get_merged_results_path("shard_trait.pip09.every", "tsv"))


def export_splicing():
    def annotate(ht):
        ht_spliceai = hl.read_table(get_spliceai_path())
        ht = ht.annotate(spliceai=ht_spliceai[ht.locus, ht.alleles].spliceai)
        ht_sqtl = hl.read_table(get_gtex_sqtl_signifpairs_path())
        ht = ht.transmute(gtex_sqtl_signif_tissues=hl.set(ht_sqtl[ht.locus, ht.alleles].values.tissue))
        ht = ht.annotate(spliceai=ht.spliceai.drop("DP_AG", "DP_AL", "DP_DG", "DP_DL"))

        # splice_consequences = hl.set(["splice_acceptor_variant", "splice_donor_variant", "splice_region_variant"])
        # ht = ht.filter(
        #     splice_consequences.contains(ht.most_severe)
        #     | (ht.spliceai.DS_MAX > 0.5)
        #     | (hl.is_defined(ht.gtex_sqtl_signif_tissues))
        # )

        return ht

    # high-confidence variants
    ht = hl.read_table(get_merged_results_path("shard_trait.pip09.any.pip01.every"))
    export_flatfile(
        annotate(ht),
        get_merged_results_path("shard_trait.pip09.any.pip01.every.splicing", "tsv"),
        extra_fields=["gtex_sqtl_signif_tissues", "spliceai"],
    )

    # population-enriched variants
    ht = hl.read_table(get_analysis_path("af_enrichment.over5"))
    ht = annotate(ht)
    ht.flatten().export(get_merged_results_path("af_enrichment.over5.splicing", "tsv"))

    # variant-trait pairs with PIP > 0.1
    ht = hl.read_table(get_merged_results_path("fm_only.in_cs.consequence"))
    ht = ht.filter(ht.pip > 0.1)
    ht = annotate(ht)
    ht.flatten().export(get_merged_results_path("fm_only.in_cs.consequence.splicing", "tsv"))


def export_enriched_variants(overwrite: bool = False):
    threshold = 5
    ref_pop = {"jpn": "njkea", "fin": "nfsee"}
    cohorts = {"jpn": "BBJ", "fin": "FG"}
    hts = []
    for data_type, pop in product(["exomes", "genomes"], ["jpn", "fin"]):
        ht = hl.read_table(get_analysis_path(f"af_enrichment.{data_type}.{cohorts[pop]}.max_pip"))
        ht = ht.filter((ht[f"enrichment_{pop}_{ref_pop[pop]}_pseudo"] > threshold) & (ht.max_pip > 0.1))

        ht_csq = hl.read_table(get_merged_results_path(f"fm_only.max_pip_per_var.consequence"))
        ht_csq = ht_csq.drop("max_pip")
        ht = ht.join(ht_csq, "left")

        ht2 = hl.read_table(get_merged_results_path("fm_only.coloc.max_pip_per_var"))
        ht2 = ht2.filter(ht2.max_pip_coloc > 0.1)
        ht2 = ht2.drop("max_pip", "eqtl")
        ht = ht.join(ht2, "left")

        ht = ht.select_globals()
        ht = ht.select(
            data_type=data_type,
            pop=pop,
            variant=hl.variant_str(ht.locus, ht.alleles),
            rsid=ht.rsid,
            most_severe=ht.most_severe,
            gene_most_severe=ht.gene_most_severe,
            consequence=ht.consequence,
            AF_cohort=ht.AF_cohort,
            AF_pop=ht[f"AF_{pop}"],
            AF_ref=hl.float64(ht[f"AF_{ref_pop[pop]}"]),
            enrichment=ht[f"enrichment_{pop}_{ref_pop[pop]}"],
            enrichment_pseudo=ht[f"enrichment_{pop}_{ref_pop[pop]}_pseudo"],
            max_pip=ht.max_pip,
            pip09_traits=ht.pip09_traits,
            pip01_traits=ht.pip01_traits,
            clinvar=ht.clinvar,
            max_pip_coloc=ht.max_pip_coloc,
            tissue_max_pip_coloc=ht.tissue_max_pip_coloc,
            trait_max_pip_coloc=ht.trait_max_pip_coloc,
            # need this to preseve the order
            **{key: ht[key] for key in get_baseline_annotations(basename=True)},
        )
        hts.append(ht)

    ht = hts[0].union(*hts[1:], unify=True)
    ht = ht.checkpoint(get_analysis_path(f"af_enrichment.over{threshold}"), overwrite=args.overwrite)
    ht.key_by().drop("locus", "alleles").export(get_analysis_path(f"af_enrichment.over{threshold}", "tsv"))


def export_variant_sumstats(variant, rsid):
    fields = [
        "cohort",
        "variant",
        "rsid",
        "region",
        "beta_marginal",
        "se_marginal",
        "chisq_marginal",
        "pvalue",
        "pip",
        "susie",
        "finemap",
    ]

    v = hl.parse_variant(variant)

    hts = []
    for pop in POPS:
        ht = hl.read_table(get_results_path(pop))
        ht = ht.select(*fields)
        hts.append(ht)

    ht = hts[0].union(*hts[1:])
    ht = ht.filter((ht.locus == v.locus) & (ht.alleles == v.alleles))
    ht.flatten().export(get_analysis_path(rsid, "tsv.bgz"))


def export_coloc_af():
    # compute max_pip_coloc for variant-cohort
    ht_vtc = hl.read_table(get_merged_results_path("fm_only.coloc.max_pip_per_vtc"))
    ht_vc = ht_vtc.group_by("locus", "alleles", "cohort").aggregate(max_pip_coloc=hl.agg.max(ht_vtc.max_pip_coloc))

    # get AF for each biobank
    ht_af = hl.read_table(get_ref_freq_info_path("BBJ_FG_UKBB"))
    ht_af = ht_af.annotate(af=hl.array([hl.struct(cohort=k, af=ht_af.af_imp[k]) for k in ht_af.af_imp]))
    ht_af = ht_af.select("af")
    ht_af = ht_af.explode("af")
    ht_af = ht_af.transmute(cohort=ht_af.af.cohort, af=ht_af.af.af)
    ht_af = ht_af.key_by("locus", "alleles", "cohort")

    ht = hl.read_table(get_merged_results_path("fm_only.max_pip_per_cohort"))
    ht = ht.filter(ht.max_pip > 0.1)
    ht = ht.join(ht_af, "left")
    ht = ht.join(ht_vc, "left")
    ht = ht.select("max_pip", "af", "max_pip_coloc")
    ht.export(get_merged_results_path("fm_only.coloc.max_pip_per_vc.max_pip01.af", "tsv.bgz"))


def export_coding_af():
    coding = hl.set(["pLoF", "Missense", "Synonymous"])
    ht = hl.read_table(get_merged_results_path("fm_only.max_pip_per_cohort.consequence"))
    ht = ht.filter((ht.max_pip > 0.1) & coding.contains(ht.consequence))

    ht_af = hl.read_table(get_ref_freq_info_path("BBJ_FG_UKBB"))
    ht_af = ht_af.annotate(af=hl.array([hl.struct(cohort=k, af=ht_af.af_imp[k]) for k in ht_af.af_imp]))
    ht_af = ht_af.select("af")
    ht_af = ht_af.explode("af")
    ht_af = ht_af.transmute(cohort=ht_af.af.cohort, af=ht_af.af.af)
    ht_af = ht_af.key_by("locus", "alleles", "cohort")

    ht = ht.join(ht_af, "left")
    ht = ht.select("max_pip", "consequence", "af")
    ht.export(get_merged_results_path("fm_only.max_pip_per_cohort.consequence.coding.max_pip01.af", "tsv.bgz"))


def export_hgvsp():
    ht = hl.read_table(get_merged_results_path("fm_only.max_pip_per_var.consequence"))
    ht = ht.filter(hl.is_defined(ht.hgvsp))
    ht = ht.key_by()
    ht = ht.select(variant=hl.variant_str(ht.locus, ht.alleles), hgvsp=ht.hgvsp.split(":")[1])
    ht.export(get_analysis_path("hgvsp", "tsv.bgz"))


def export_ems():
    ht = hl.read_table(get_merged_results_path("fm_only.max_pip_per_var.consequence"))
    ht = ht.filter(hl.is_defined(ht.max_ems_normalized) & hl.is_defined(ht.max_pip))
    ht = ht.key_by()
    ht = ht.annotate(variant=hl.variant_str(ht.locus, ht.alleles))
    ht = ht.select("variant", "max_pip", "max_ems_normalized", "gene_max_ems_normalized", "consequence")
    ht.export(get_analysis_path("max_ems_normalized", "tsv.bgz"))


def export_sim_gamma_pip():
    ht = hl.read_table("gs://xfinemap/ukbb_simulation/results_nomerge/combined/sim1-50.ALL.ht")
    ht = ht.filter(
        (ht.genotype == "Dosage")
        & (ht.ld == "Dosage")
        & (ht.method == "AVERAGE")
        & (ht.maf_threshold == 1e-4)
        & (ht.info_threshold == 0.8)
        & (ht.window_size == 1500000)
    )
    ht.filter(ht.gamma_gwas & (ht.pvalue < 5e-8)).select().export(
        get_analysis_path("ukbb_sim.gamma_gwas.sig.pip", "tsv.bgz")
    )


def main(args):
    export_pip09_variants(overwrite=args.overwrite)
    export_csm_summary(overwrite=args.overwrite)
    export_coloc_vtg(overwrite=args.overwrite)
    export_coloc_af(overwrite=args.overwrite)
    export_pip09_in_every_pop(overwrite=args.overwrite)
    export_enriched_variants(overwrite=args.overwrite)
    export_splicing()
    export_variant_sumstats("17:7080316:C:T", "rs55714927")
    export_variant_sumstats("5:176509193:C:T", "rs244711")
    export_variant_sumstats("5:176516631:G:A", "rs1966265")
    export_variant_sumstats("16:57017662:G:A", "rs1801706")
    export_coloc_af()
    export_coding_af()
    export_hgvsp()
    export_ems()
    export_sim_gamma_pip()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
