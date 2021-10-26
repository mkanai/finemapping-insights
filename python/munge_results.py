import argparse
from fm_insights.resources.annotations import get_ems_path
import hail as hl
import scipy as sp
from fm_insights.resources import (
    POPS,
    get_analysis_path,
    get_clinvar_path,
    get_dbsnp_path,
    get_merged_results_path,
    get_raw_finemap_path,
    get_raw_sumstats_path,
    get_raw_susie_path,
    get_ref_freq_info_path,
    get_results_path,
    get_sample_size_dict,
    get_trait_mapping_dict,
    get_trait_summary,
    get_vep_annot_path,
)
from fm_insights.utils import checkpoint_tmp, annotate_bed, register_log


def bbj_filter_f(ht):
    # remove imputed SVs for now (e.g., <INS:ME:ALU>)
    return ~ht.variant.contains("<")


def fg_select_f(ht):
    expr = {
        "cohort": ht.cohort,
        "trait": ht.trait,
        "region": ht.region,
        "variant": ht.v,
        "variant_b37": ht.variant_b37,
        "prob": ht.prob,
        "cs": ht.cs,
        "mean": ht.mean,
        "sd": ht.sd,
        "need_to_flip_beta": ht.need_to_flip_beta,
    }
    if "low_purity" in ht.row:
        expr["low_purity"] = ht.low_purity
    return expr


def read_sumstats(path, trait_dict, filter_f=None, liftover=False):
    ht = hl.import_table(
        path,
        impute=True,
        min_partitions=3000,
        types={"beta_marginal": hl.tfloat64, "se_marginal": hl.tfloat64, "chisq_marginal": hl.tfloat64},
    )
    if filter_f is not None:
        ht = ht.filter(filter_f(ht))
    if liftover:
        ht = ht.rename({"variant_b37": "variant", "variant": "variant_b38"})
        ht = ht.filter(hl.is_defined(ht.variant))
        ht = ht.annotate(beta_marginal=hl.if_else(ht.need_to_flip_beta, -ht.beta_marginal, ht.beta_marginal))

    ht = ht.annotate(**hl.parse_variant(ht.variant))
    # map trait_cohort to trait
    ht = ht.annotate(trait=trait_dict[ht.trait], trait_cohort=ht.trait)
    ht = ht.key_by("locus", "alleles", "trait")
    return ht


def read_snp(path, method, trait_dict, filter_f=None, select_f=None, liftover=False):
    ht = hl.import_table(path, min_partitions=200, impute=True)
    if filter_f is not None:
        ht = ht.filter(filter_f(ht))
    if select_f is not None:
        ht = ht.select(**select_f(ht))
    if liftover:
        ht = ht.rename({"variant_b37": "variant", "variant": "variant_b38"})
        ht = ht.filter(hl.is_defined(ht.variant))
        ht = ht.annotate(mean=hl.if_else(ht.need_to_flip_beta, -ht.mean, ht.mean))

    ht = ht.annotate(**hl.parse_variant(ht.variant))

    # map trait_cohort to trait
    ht = ht.annotate(trait=trait_dict[ht.trait], trait_cohort=ht.trait)
    ht = ht.key_by("locus", "alleles", "trait")

    # remove variants in MHC region
    MHC_region = hl.interval(hl.locus("6", 25000000), hl.locus("6", 36000000), includes_end=True)
    ht = ht.annotate(in_MHC=MHC_region.contains(ht.locus))
    ht = ht.filter(~ht.in_MHC)

    if "low_purity" in ht.row:
        ht = ht.annotate(cs=hl.if_else((ht.cs == -1) | (ht.low_purity == 1), -1, ht.cs))

    if "alpha1" in ht.row:
        susie_alpha_expr = {"alpha": hl.array([ht[f"alpha{i}"] for i in range(1, 11)])}
    else:
        susie_alpha_expr = {"alpha": hl.null(hl.tarray(hl.tfloat64))}

    ht = ht.select(
        method=method,
        region=ht.region,
        pip=ht.prob,
        cs_id=ht.cs,
        beta_posterior=ht.mean,
        sd_posterior=ht.sd,
        **susie_alpha_expr,
    )
    return ht


def munge_results(pop, trait="*", filter_f=None, select_f=None, liftover=False, sim=False, overwrite=False):
    if not sim:
        trait_dict = get_trait_mapping_dict(pop)
    else:
        ht_trait = get_trait_summary(pop="ALL", min_n_pop=3)
        trait_dict = (
            ht_trait.select("trait_cohort", "trait").to_pandas().set_index("trait_cohort").T.to_dict("records")[0]
        )
        trait_dict = hl.literal({f"{v}_from_{pop}": v for _, v in trait_dict.items() for pop in POPS})

    ht_sumstats = read_sumstats(
        get_raw_sumstats_path(pop, trait, sim=sim), trait_dict, filter_f=filter_f, liftover=liftover,
    )
    if sim:
        ht_sumstats = ht_sumstats.annotate(cohort=ht_sumstats.cohort.split("_")[0])
    ht_sumstats = checkpoint_tmp(ht_sumstats)

    ht_susie = read_snp(
        get_raw_susie_path(pop, trait, sim=sim),
        "SUSIE",
        trait_dict,
        filter_f=filter_f,
        select_f=select_f,
        liftover=liftover,
    )
    ht_susie = checkpoint_tmp(ht_susie)

    ht_finemap = read_snp(
        get_raw_finemap_path(pop, trait, sim=sim),
        "FINEMAP",
        trait_dict,
        filter_f=filter_f,
        select_f=select_f,
        liftover=liftover,
    )
    ht_finemap = checkpoint_tmp(ht_finemap)

    ht = ht_susie.union(ht_finemap)
    ht = checkpoint_tmp(ht)

    ht = ht.collect_by_key()
    ht = ht.transmute(
        method="Average",
        region=ht.values.region[0],
        pip=(
            hl.case()
            .when(hl.len(ht.values) == 1, ht.values.pip[0])
            .when(hl.is_missing(ht.values.pip[0]), ht.values.pip[1])
            .when(hl.is_missing(ht.values.pip[1]), ht.values.pip[0])
            .when(hl.abs(ht.values.pip[0] - ht.values.pip[1]) < 0.05, hl.mean(ht.values.pip))
            .or_missing()
        ),
        finemap=hl.struct(
            pip=ht.values.pip[ht.values.method.index("FINEMAP")],
            cs_id=ht.values.cs_id[ht.values.method.index("FINEMAP")],
            beta_posterior=ht.values.beta_posterior[ht.values.method.index("FINEMAP")],
            sd_posterior=ht.values.sd_posterior[ht.values.method.index("FINEMAP")],
        ),
        susie=hl.struct(
            pip=ht.values.pip[ht.values.method.index("SUSIE")],
            cs_id=ht.values.cs_id[ht.values.method.index("SUSIE")],
            beta_posterior=ht.values.beta_posterior[ht.values.method.index("SUSIE")],
            sd_posterior=ht.values.sd_posterior[ht.values.method.index("SUSIE")],
            alpha=ht.values.alpha[ht.values.method.index("SUSIE")],
        ),
    )
    ht = checkpoint_tmp(ht)

    ht_sumstats = ht_sumstats.join(ht, "left")
    ht_sumstats = ht_sumstats.checkpoint(get_results_path(pop, sim=sim), overwrite=overwrite)


def get_cs_max_regions(ht, L=10):
    ht = ht.group_by(trait_region=hl.delimit([ht.trait, ht.region])).aggregate(
        uniq_cs_id=hl.agg.collect_as_set(ht.susie.cs_id)
    )
    # '-1' + L CSs
    ht = ht.filter(hl.len(ht.uniq_cs_id) == (L + 1))
    trait_regions = ht.aggregate(hl.agg.collect_as_set(ht.trait_region), _localize=False)
    return trait_regions


def merge_results(overwrite: bool = False):
    fields = [
        "cohort",
        "variant",
        "rsid",
        "region",
        "af",
        "beta_marginal",
        "se_marginal",
        "chisq_marginal",
        "pvalue",
        "pip",
        "susie",
        "finemap",
    ]

    ht_trait = get_trait_summary(pop="ALL", min_n_pop=2)
    traits = ht_trait.aggregate(hl.agg.collect_as_set(ht_trait.trait), _localize=False)

    hts = []
    hts2 = []
    for pop in POPS:
        ht = hl.read_table(get_results_path(pop))

        ht2 = ht.filter(get_cs_max_regions(ht).contains(hl.delimit([ht.trait, ht.region])))
        ht2 = ht2.key_by("cohort", "trait", "region").select().distinct()
        ht2 = ht2.checkpoint(get_analysis_path(f"{pop}.cs_max_regions"), overwrite=overwrite)
        hts2.append(ht2)

        ht = ht.annotate(af=hl.if_else(ht.alleles[1] == ht.minorallele, ht.maf, 1 - ht.maf))
        ht = ht.select(*fields)
        hts.append(ht)

    ht = hl.Table.union(*hts)
    ht2 = hl.Table.union(*hts2)
    ht2.export(get_merged_results_path("cs_max_regions", "tsv.bgz"))

    # output only fm merged ht
    ht_fm = ht.filter(hl.is_defined(ht.finemap) | hl.is_defined(ht.susie))
    ht_fm = ht_fm.checkpoint(get_merged_results_path("fm_only"), overwrite=overwrite)

    # filter to shared traits, output, collect
    ht = ht.filter(traits.contains(ht.trait))
    ht = ht.collect_by_key()
    ht = ht.checkpoint(get_merged_results_path("shard_trait"), overwrite=overwrite)

    return ht


def merge_sim_results(overwrite: bool = False):
    ht_trait = get_trait_summary(pop="ALL", min_n_pop=3)
    traits = ht_trait.aggregate(hl.agg.collect_as_set(ht_trait.trait), _localize=False)

    ht_fm = hl.read_table(get_merged_results_path("fm_only"))
    ht_fm = ht_fm.filter(traits.contains(ht_fm.trait))
    ht_fm = ht_fm.select("cohort", "chisq_marginal", "pip")
    ht_fm = ht_fm.key_by("locus", "alleles", "trait", "cohort")
    ht_fm = checkpoint_tmp(ht_fm)

    hts = []
    for pop in POPS:
        ht = hl.read_table(get_results_path(pop, sim=True))
        ht = ht.annotate(
            # to merge with the original results -- for each variant-trait in discovert cohort
            cohort=ht.trait_cohort.split("_")[-1],
            sim=hl.struct(cohort=ht.cohort, chisq_marginal=ht.chisq_marginal, pip=ht.pip),
        )
        ht = ht.select("cohort", "sim")
        ht = ht.key_by("locus", "alleles", "trait", "cohort")
        hts.append(ht)
    ht_sim = hl.Table.union(*hts)
    ht_sim = checkpoint_tmp(ht_sim)

    ht = ht_fm.join(ht_sim, "left")
    ht = ht.filter((ht.pip > 0.01) | (ht.sim.pip > 0.01))
    ht = ht.checkpoint(get_merged_results_path("sim.fm_only.pip001"), overwrite=overwrite)

    ht = ht.filter((ht.pip > 0.9) | (ht.sim.pip > 0.9))
    ht.flatten().export(get_merged_results_path("sim.fm_only.pip09", "tsv.bgz"))


def compute_max_pip(overwrite: bool = False):
    ht = hl.read_table(get_merged_results_path("fm_only"))

    # max_pip per variant-trait
    ht_max_pip = ht.group_by("locus", "alleles", "trait").aggregate(
        max_pip=hl.agg.max(ht.pip), min_pip=hl.agg.min(ht.pip), n_pop=hl.agg.sum(hl.is_defined(ht.pip))
    )
    ht_max_pip = ht_max_pip.checkpoint(get_merged_results_path("fm_only.max_pip"), overwrite=overwrite)

    # max_pip per variant
    ht_max_pip_per_var = ht_max_pip.group_by("locus", "alleles").aggregate(max_pip=hl.agg.max(ht_max_pip.max_pip))
    ht_max_pip_per_var = ht_max_pip_per_var.checkpoint(
        get_merged_results_path("fm_only.max_pip_per_var"), overwrite=overwrite
    )
    # max_pip per cohort
    ht_max_pip_per_cohort = ht.group_by("locus", "alleles", "cohort").aggregate(max_pip=hl.agg.max(ht.pip))
    ht_max_pip_per_cohort = ht_max_pip_per_cohort.checkpoint(
        get_merged_results_path("fm_only.max_pip_per_cohort"), overwrite=overwrite
    )

    # per-cs
    ht_cs = ht.filter(ht.susie.cs_id > 0)
    ht_cs = ht_cs.group_by("cohort", "trait", "region", ht_cs.susie.cs_id).aggregate(
        cs_size=hl.agg.count(), max_pip=hl.agg.max(ht_cs.pip)
    )
    ht_cs = ht_cs.checkpoint(get_merged_results_path("fm_only.cs"), overwrite=overwrite)

    # cs_ids for CS-merging
    ht_cs_id = ht.filter(ht.susie.cs_id > 0)
    ht_cs_id = ht_cs_id.flatten()
    fields = ["cohort", "trait", "variant", "region", "susie.pip", "susie.cs_id"]
    ht_cs_id.select(*fields).export(get_merged_results_path("fm_only.cs_id", "tsv.bgz"))


def annotate(overwrite: bool = False):
    ht = hl.read_table(get_merged_results_path("fm_only.max_pip_per_var"))

    # VEP
    ht_vep = hl.read_table(get_vep_annot_path())
    ht = ht.join(ht_vep, "left")

    # Clinvar
    ht_clinvar = hl.read_table(get_clinvar_path())
    # cf. gnomad.utils.filtering.filter_to_clinvar_pathogenic
    # filter no_assertion
    no_star_assertions = hl.literal(
        {"no_assertion_provided", "no_assertion_criteria_provided", "no_interpretation_for_the_individual_variant"}
    )
    ht_clinvar = ht_clinvar.filter((hl.set(ht_clinvar.info.CLNREVSTAT).intersection(no_star_assertions).length() == 0))
    ht_clinvar = ht_clinvar.select(
        clinvar=hl.case()
        .when(
            ht_clinvar.info.CLNSIG.map(lambda x: x.lower()).map(lambda x: x.contains("pathogenic")).any(lambda x: x),
            "Pathogenic",
        )
        .when(hl.is_defined(ht_clinvar.info.CLNSIGCONF), "Conflicting")
        .default("Non-pathogenic")
    )
    ht = ht.join(ht_clinvar, "left")

    # dbsnp
    ht_dbsnp = hl.read_table(get_dbsnp_path())
    ht = ht.join(ht_dbsnp, "left")

    # ems
    ht_ems = hl.read_table(get_ems_path("Whole_Blood"))
    ht_ems = ht_ems.select("max_ems_normalized", "gene_max_ems_normalized")
    ht = ht.join(ht_ems, "left")
    annots = ht.row.keys()

    # annotate with bed
    ht = annotate_bed(ht)
    annots = sorted(list(ht.row.keys() - annots))

    ht = ht.annotate(
        consequence=(
            hl.case(missing_false=True)
            .when(hl.is_defined(ht.lof) & (ht.lof != "LC"), "pLoF")
            .when(
                (ht.lof == "LC")
                | (ht.consequence_category == "coding_high")
                | (ht.consequence_category == "coding_medium"),
                "Missense",
            )
            .when(ht.consequence_category == "coding_low", "Synonymous")
            .when(ht.most_severe == "3_prime_UTR_variant", "UTR3")
            .when(ht.most_severe == "5_prime_UTR_variant", "UTR5")
            .when(ht.Promoter_UCSC == 1, "Promoter")
            .when(
                (ht.DHSmerged_Ulirsch == 1) & ((ht.Roadmap_H3K27ac_Ulirsch == 1) | (ht.CA_H3K27ac_Ulirsch == 1)), "CRE"
            )
            .default("Non-genic")
        )
    )

    ht = ht.checkpoint(get_merged_results_path("fm_only.max_pip_per_var.consequence"), overwrite=overwrite)
    ht.select("max_pip", "consequence").export(
        get_merged_results_path("fm_only.max_pip_per_var.consequence", "tsv.bgz")
    )

    ht.select("max_pip", "consequence", "clinvar", *annots).export(
        get_merged_results_path("fm_only.max_pip_per_var.consequence.annot", "tsv.bgz")
    )

    # annotate per-cohort
    ht_cohort = hl.read_table(get_merged_results_path("fm_only.max_pip_per_cohort"))
    ht_cohort = ht_cohort.annotate(**ht[ht_cohort.locus, ht_cohort.alleles].drop("max_pip"))
    ht_cohort = ht_cohort.checkpoint(
        get_merged_results_path("fm_only.max_pip_per_cohort.consequence"), overwrite=overwrite
    )
    ht_cohort.select("max_pip", "consequence").export(
        get_merged_results_path("fm_only.max_pip_per_cohort.consequence", "tsv.bgz")
    )

    # annotate in-CS variants too
    ht_cs = hl.read_table(get_merged_results_path("fm_only"))
    ht_cs = ht_cs.filter((ht_cs.susie.cs_id > 0) | (ht_cs.pip > 0.1))
    ht_cs = ht_cs.annotate(**ht[ht_cs.locus, ht_cs.alleles])
    ht_cs = ht_cs.checkpoint(get_merged_results_path("fm_only.in_cs.consequence"), overwrite=overwrite)
    fields = [
        "locus",
        "alleles",
        "trait",
        "cohort",
        "region",
        "rsid",
        "pip",
        "susie.cs_id",
        "susie.beta_posterior",
        "consequence",
        "most_severe",
        "gene_most_severe",
    ]
    ht_cs.flatten().select(*fields).export(get_merged_results_path("fm_only.in_cs.consequence", "tsv.bgz"))


def meta_analyze_beta(ht):
    ht = ht.annotate(n_pop=hl.sum(ht.values.beta_marginal.map(hl.is_defined)))
    ht = ht.annotate(
        unnorm_beta=ht.values.beta_marginal / (ht.values.se_marginal ** 2), inv_se2=1 / (ht.values.se_marginal ** 2),
    )
    ht = ht.annotate(sum_unnorm_beta=hl.sum(ht.unnorm_beta), sum_inv_se2=hl.sum(ht.inv_se2))
    ht = ht.annotate(beta_meta=ht.sum_unnorm_beta / ht.sum_inv_se2, se_meta=hl.sqrt(1 / ht.sum_inv_se2),)
    ht = ht.annotate(
        pvalue_meta=2 * hl.pnorm(-hl.abs(ht.beta_meta / ht.se_meta)),
        q_meta=hl.sum((ht.values.beta_marginal - ht.beta_meta) ** 2 * ht.inv_se2),
    )
    ht = ht.annotate(pvalue_het=hl.if_else(ht.n_pop > 1, hl.pchisqtail(ht.q_meta, ht.n_pop - 1), hl.null(hl.tfloat64)))
    ht = ht.drop("unnorm_beta", "inv_se2", "sum_unnorm_beta", "sum_inv_se2", "q_meta")
    return ht


def compare_pops(overwrite: bool = False):
    ht = hl.read_table(get_merged_results_path("shard_trait"))

    def get_per_pop_values(field, pops=POPS):
        expr = hl.struct(
            **{
                pop: hl.rbind(
                    ht.values.cohort.index(pop), lambda idx: hl.or_missing(hl.is_defined(idx), ht.values[field][idx]),
                )
                for pop in pops
            }
        )
        return expr

    ht = ht.filter(hl.any(lambda x: x > 0.01, ht.values.pip))
    ht = meta_analyze_beta(ht)
    ht = ht.repartition(1000)
    ht = ht.checkpoint(get_merged_results_path("shard_trait.pip001"), overwrite=overwrite)

    fields = ["pip", "cs_id", "beta_posterior", "sd_posterior"]
    ht = ht.annotate(
        values=ht.values.map(
            lambda x: x.annotate(
                **x.susie.flatten().select(*fields).rename({y: f"susie.{y}" for y in fields}),
                **x.finemap.flatten().select(*fields).rename({y: f"finemap.{y}" for y in fields}),
            )
        )
    )
    ht = ht.annotate(
        variant=hl.variant_str(ht.locus, ht.alleles),
        beta_marginal=get_per_pop_values("beta_marginal"),
        se_marginal=get_per_pop_values("se_marginal"),
        chisq_marginal=get_per_pop_values("chisq_marginal"),
        pvalue=get_per_pop_values("pvalue"),
        af=get_per_pop_values("af"),
        pip=get_per_pop_values("pip"),
        # power=get_per_pop_values("power"),
        susie=hl.struct(
            pip=get_per_pop_values("susie.pip"),
            cs_id=get_per_pop_values("susie.cs_id"),
            beta_posterior=get_per_pop_values("susie.beta_posterior"),
            sd_posterior=get_per_pop_values("susie.sd_posterior"),
        ),
        finemap=hl.struct(
            pip=get_per_pop_values("finemap.pip"),
            cs_id=get_per_pop_values("finemap.cs_id"),
            beta_posterior=get_per_pop_values("finemap.beta_posterior"),
            sd_posterior=get_per_pop_values("finemap.sd_posterior"),
        ),
    )
    ht = ht.drop("values")
    ht.describe()

    # annotate imputation af/info
    ht_ref = hl.read_table(get_ref_freq_info_path("BBJ_FG_UKBB"))
    ht = ht.annotate(**ht_ref[ht.locus, ht.alleles])
    ht_qc = hl.read_table("gs://ukb31063/ukb31063.variant_qc.both_sexes_gwas_samples.autosomes.ht")
    ht = ht.annotate(**{"p_value_hwe.UKBB": ht_qc[ht.locus, ht.alleles].variant_qc.p_value_hwe})

    ht = ht.checkpoint(get_merged_results_path("shard_trait.pip001.pop"), overwrite=overwrite)
    ht.flatten().export(get_merged_results_path("shard_trait.pip001.pop", "tsv.bgz"))

    ht = ht.filter(hl.any(lambda x: x > 0.9, [ht.pip.BBJ, ht.pip.FG, ht.pip.UKBB]))
    ht = ht.checkpoint(get_merged_results_path("shard_trait.pip09.pop"), overwrite=overwrite)
    ht.flatten().export(get_merged_results_path("shard_trait.pip09.pop", "tsv.bgz"))


def annotate_power(overwrite: bool = False):
    ht_trait = get_trait_summary(pop="ALL")
    ht_agg = ht_trait.group_by(ht_trait.trait).aggregate(n_pop=hl.agg.count())
    ht_agg = ht_agg.filter(ht_agg.n_pop == 3)
    traits = ht_agg.aggregate(hl.agg.collect_as_set(ht_agg.trait), _localize=False)

    alpha = 5e-8
    threshold = sp.stats.norm.ppf(alpha / 2) ** 2

    ht = hl.read_table(get_merged_results_path("shard_trait"))
    ht = ht.filter(traits.contains(ht.trait))
    n_dict = hl.dict({pop: get_sample_size_dict(pop, neff=True) for pop in POPS})
    ht = ht.annotate(values=ht.values.map(lambda x: x.annotate(N=n_dict[x.cohort][ht.trait])))
    ht = ht.select(
        values=hl.flatten(
            ht.values.map(
                lambda x: hl.rbind(
                    hl.pchisqtail(threshold, 1, ncp=x.N * 2 * x.af * (1 - x.af) * (x.beta_marginal ** 2)),
                    lambda power_base: ht.values.map(
                        lambda y: hl.or_missing(
                            x.cohort != y.cohort,
                            hl.struct(
                                base=x.cohort,
                                cohort=y.cohort,
                                pvalue_base=x.pvalue,
                                pvalue_cohort=y.pvalue,
                                pip_base=x.pip,
                                pip_cohort=y.pip,
                                cs_id_base=x.susie.cs_id,
                                cs_id_cohort=y.susie.cs_id,
                                power_base=power_base,
                                power_cohort=hl.pchisqtail(
                                    threshold, 1, ncp=y.N * 2 * y.af * (1 - y.af) * (x.beta_marginal ** 2)
                                ),
                            ),
                        )
                    ),
                )
            )
        )
    )
    ht = ht.explode("values")
    ht = ht.filter(hl.is_defined(ht.values))
    ht = ht.transmute(**ht.values)
    ht = ht.checkpoint(get_merged_results_path("shard_trait.power"), overwrite=overwrite)

    ht.filter(ht.pvalue_base < 5e-8).export(get_merged_results_path("shard_trait.power.sig_base", "tsv.bgz"))


def main(args):
    if args.bbj:
        munge_results("BBJ", trait=args.trait, filter_f=bbj_filter_f, sim=args.sim, overwrite=args.overwrite)
    if args.ukbb:
        munge_results("UKBB", trait=args.trait, sim=args.sim, overwrite=args.overwrite)
    if args.fg:
        munge_results(
            "FG", trait=args.trait, select_f=fg_select_f, liftover=True, sim=args.sim, overwrite=args.overwrite,
        )
    if args.merge:
        if not args.sim:
            merge_results(overwrite=args.overwrite)
        else:
            merge_sim_results(overwrite=args.overwrite)
    if args.compute_max_pip:
        compute_max_pip(overwrite=args.overwrite)
    if args.annotate:
        annotate(overwrite=args.overwrite)
    if args.compare_pops:
        compare_pops(overwrite=args.overwrite)
    if args.annotate_power:
        annotate_power(overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bbj", action="store_true")
    parser.add_argument("--ukbb", action="store_true")
    parser.add_argument("--fg", action="store_true")
    parser.add_argument("--sim", action="store_true")
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--compute-max-pip", action="store_true")
    parser.add_argument("--annotate", action="store_true")
    parser.add_argument("--compare-pops", action="store_true")
    parser.add_argument("--annotate-power", action="store_true")
    parser.add_argument("--trait", type=str, default="*")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
