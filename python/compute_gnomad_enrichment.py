import argparse

from hail.typecheck import check
from fm_insights.utils.generic import checkpoint_tmp
import hail as hl
from hail.linalg import BlockMatrix
from itertools import product

from gnomad.utils.vep import (
    filter_vep_to_canonical_transcripts,
    get_most_severe_consequence_for_summary,
    CSQ_CODING_HIGH_IMPACT,
    CSQ_CODING_MEDIUM_IMPACT,
    CSQ_CODING_LOW_IMPACT,
)
from gnomad.resources.grch37.gnomad import public_release, coverage
from gnomad.resources.grch37.gnomad_ld import ld_matrix, ld_index
from gnomad.resources.grch37.reference_data import cpg_sites
from fm_insights.utils import register_log
from fm_insights.resources import (
    get_analysis_path,
    get_results_path,
    get_gwas_variants_path,
    get_gem_j_wga_path,
    get_human_genome_dating_path,
)

coding = hl.set(CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT)


def read_gem_j_wga(data_type: str, coverage_threshold=10):
    ht = hl.read_table(get_gem_j_wga_path(vep=True))
    ht = ht.filter(hl.len(ht.filters) == 0)

    ht_coverage = coverage(data_type).versions["2.1"].ht()
    ht_coverage = ht_coverage.filter(ht_coverage.median > coverage_threshold)
    ht = ht.filter(hl.is_defined(ht_coverage[ht.locus]))

    if data_type == "exomes":
        ht = ht.filter(coding.contains(ht.most_severe_csq))
    else:
        ht = ht.filter(hl.is_missing(ht.most_severe_csq) | ~coding.contains(ht.most_severe_csq))

    ht = ht.select(AC_jpn=ht.info.AC, AN_jpn=ht.info.AN)
    return ht


def calc_pseudo_freq(ht, ac, an):
    max_an = ht.aggregate(hl.agg.max(an))
    return hl.case().when(hl.is_missing(ac) | hl.is_missing(an) | (an == 0), 1 / max_an).default((ac + 1) / an)


def calc_enrichment(numer, denom):
    return (
        hl.case()
        .when(hl.is_missing(numer) | (numer == 0), hl.null(hl.tfloat64))
        .when(hl.is_missing(denom) | (denom == 0), float("inf"))
        .default(numer / denom)
    )


def annotate_coding_tagging(ht_genomes, r2_threshold=0.1):
    hts = {"exomes": hl.read_table(get_analysis_path("af_enrichment.exomes.fin_jpn")), "genomes": ht_genomes}

    for pop in ["fin", "jpn"]:
        pop_gnomad = "eas" if pop == "jpn" else pop

        # extract coding/non-coding variant indicies
        ht_idx = ld_index(pop_gnomad).ht()
        idx = {}
        for data_type in ["exomes", "genomes"]:
            ht = hts[data_type].filter(hts[data_type][f"AF_{pop}"] > 0)
            ht = ht.join(ht_idx, "inner")
            idx[data_type] = ht.aggregate(hl.agg.collect_as_set(ht.idx), _localize=False)

        # compute r2 bm
        bm = ld_matrix(pop_gnomad).bm()
        bm = bm ** 2

        # prune to coding/non-coding pairs with r2 > r2_threshold
        entries = bm.entries(keyed=False)
        entries = entries.filter((entries.entry >= r2_threshold) & (entries.i < entries.j))
        entries = entries.filter(
            (idx["exomes"].contains(entries.i) & idx["genomes"].contains(entries.j))
            | (idx["exomes"].contains(entries.j) & idx["genomes"].contains(entries.i))
        )

        # get tagging variant indices
        tagging_idx = entries.aggregate(
            hl.agg.collect_as_set(entries.i).union(hl.agg.collect_as_set(entries.j)), _localize=False
        )
        ht_idx = ht_idx.filter(tagging_idx.contains(ht_idx.idx)).select().select_globals()
        ht_idx = checkpoint_tmp(ht_idx)

        hts["genomes"] = hts["genomes"].annotate(
            **{f"coding_tagging_{pop}": hl.is_defined(ht_idx[hts["genomes"].locus, hts["genomes"].alleles])}
        )

    return hts["genomes"]


def compute_af_enrichment(data_type: str, overwrite: bool = False):
    # load gnomAD
    ht = public_release(data_type).versions["2.1.1"].ht()
    ht = filter_vep_to_canonical_transcripts(ht)
    ht = get_most_severe_consequence_for_summary(ht)

    ht_coverage = coverage(data_type).versions["2.1"].ht()
    ht_coverage = ht_coverage.filter(ht_coverage.median > 10)
    ht = ht.filter(hl.is_defined(ht_coverage[ht.locus]))

    if data_type == "exomes":
        ht = ht.filter(coding.contains(ht.most_severe_csq))
    else:
        ht = ht.filter(hl.is_missing(ht.most_severe_csq) | ~coding.contains(ht.most_severe_csq))

    # merge Japanese AF
    ht_jpn = read_gem_j_wga(data_type)
    ht = ht.join(ht_jpn, "outer")
    ht = ht.annotate(
        AF_jpn=hl.or_missing(hl.is_finite(ht.AN_jpn), ht.AC_jpn / ht.AN_jpn),
        AF_jpn_pseudo=calc_pseudo_freq(ht, ht.AC_jpn, ht.AN_jpn),
    )

    # get gnomad indexes
    idx = {}
    if data_type == "exomes":
        gnomad_pops = ["afr", "fin", "nfe", "nfe_swe", "nfe_est"]
        # gnomad eas_oea = Non-Japanese-Korean East Asian (NJKEA)
        idx["njkea"] = ht.freq_index_dict[f"gnomad_eas_oea"].collect()[0]
    else:
        gnomad_pops = ["afr", "fin", "nfe_nwe", "nfe_onf", "nfe_seu"]
        # no Japanese/Korean genome available in r2.1.1
        idx["njkea"] = ht.freq_index_dict[f"gnomad_eas"].collect()[0]
    idx.update({pop: ht.freq_index_dict[f"gnomad_{pop}"].collect()[0] for pop in gnomad_pops})

    # compute Non-Finnish-Swedish-Estonian-European (NFSEE) AF
    if data_type == "exomes":
        ac_nfsee = hl.sum([ht.freq[idx["nfe"]].AC, -ht.freq[idx["nfe_swe"]].AC, -ht.freq[idx["nfe_est"]].AC])
        an_nfsee = hl.sum([ht.freq[idx["nfe"]].AN, -ht.freq[idx["nfe_swe"]].AN, -ht.freq[idx["nfe_est"]].AN])
    else:
        ac_nfsee = hl.sum([ht.freq[idx["nfe_nwe"]].AC, ht.freq[idx["nfe_onf"]].AC, ht.freq[idx["nfe_seu"]].AC])
        an_nfsee = hl.sum([ht.freq[idx["nfe_nwe"]].AN, ht.freq[idx["nfe_onf"]].AN, ht.freq[idx["nfe_seu"]].AN])
    ht = ht.annotate(
        AF_nfsee=hl.or_missing(hl.is_finite(an_nfsee), ac_nfsee / an_nfsee),
        AF_nfsee_pseudo=calc_pseudo_freq(ht, ac_nfsee, an_nfsee),
    )

    # compute AF and pseudo-AF
    pops = ["afr", "fin", "njkea"]
    af_expr = {f"AF_{pop}": ht.freq[idx[pop]].AF for pop in pops}
    af_pseudo_expr = {
        f"AF_{pop}_pseudo": calc_pseudo_freq(ht, ht.freq[idx[pop]].AC, ht.freq[idx[pop]].AN) for pop in pops
    }
    jpn_nfsee = [x for pop in ["jpn", "nfsee"] for x in [f"AF_{pop}", f"AF_{pop}_pseudo"]]
    ht = ht.select(*jpn_nfsee, **af_expr, **af_pseudo_expr)
    ht = ht.filter((ht.AF_fin > 0) | (ht.AF_jpn > 0))

    # exclude non-coding variants in LD with coding variants
    if data_type == "genomes":
        ht = checkpoint_tmp(ht)
        ht = annotate_coding_tagging(ht)

    # compute enrichment
    pop_comb = list(product(["jpn", "fin"], ["afr", "njkea", "nfsee"]))
    enrichment_expr = {
        f"enrichment_{pop1}_{pop2}": calc_enrichment(ht[f"AF_{pop1}"], ht[f"AF_{pop2}"]) for pop1, pop2 in pop_comb
    }
    enrichment_pseudo_expr = {
        f"enrichment_{pop1}_{pop2}_pseudo": calc_enrichment(ht[f"AF_{pop1}_pseudo"], ht[f"AF_{pop2}_pseudo"])
        for pop1, pop2 in pop_comb
    }
    ht = ht.annotate(**enrichment_expr, **enrichment_pseudo_expr)

    # annotate CpG context
    ht_cpg = cpg_sites.ht()
    ht = ht.annotate(
        variant_type=hl.case()
        .when(hl.is_defined(ht_cpg[ht.locus, ht.alleles]), "CpG")
        .when(hl.is_transition(ht.alleles[0], ht.alleles[1]), "Non-CpG transition")
        .default("Transversion")
    )

    # merge variant age
    ht_age = hl.read_table(get_human_genome_dating_path())
    ht_age = ht_age.select(
        age_mode=ht_age.AgeMode_Jnt,
        age_mean=ht_age.AgeMean_Jnt,
        age_median=ht_age.AgeMedian_Jnt,
        age_lower=ht_age.AgeCI95Lower_Jnt,
        age_upper=ht_age.AgeCI95Upper_Jnt,
        age_qual_score=ht_age.QualScore_Jnt,
    )
    ht = ht.join(ht_age, "left")

    ht = ht.select_globals()
    ht = ht.checkpoint(get_analysis_path(f"af_enrichment.{data_type}.fin_jpn"), overwrite=overwrite)
    ht.export(get_analysis_path(f"af_enrichment.{data_type}.fin_jpn", "tsv.bgz"))


def export_enrichment(data_type: str, cohort: str, overwrite: bool = False):
    ht = hl.read_table(get_analysis_path(f"af_enrichment.{data_type}.fin_jpn"))

    # filter to GWAS variants
    ht_variants = hl.read_table(get_gwas_variants_path(cohort))
    ht = ht.join(ht_variants.select(), "inner")

    # compute max_pip per cohort
    ht_results = hl.read_table(get_results_path(cohort))
    ht_results = ht_results.filter(hl.is_defined(ht_results.pip))
    ht_results = ht_results.annotate(
        AF_cohort=hl.if_else(ht_results.alleles[1] == ht_results.minorallele, ht_results.maf, 1 - ht_results.maf)
    )
    ht_results = ht_results.group_by("locus", "alleles").aggregate(
        AF_cohort=hl.agg.mean(ht_results.AF_cohort),
        max_pip=hl.agg.max(ht_results.pip),
        pip09_traits=hl.agg.filter(ht_results.pip > 0.9, hl.agg.collect_as_set(ht_results.trait)),
        pip01_traits=hl.agg.filter(
            (ht_results.pip <= 0.9) & (ht_results.pip > 0.1), hl.agg.collect_as_set(ht_results.trait)
        ),
    )

    ht = ht.join(ht_results, "left")

    ht = ht.checkpoint(get_analysis_path(f"af_enrichment.{data_type}.{cohort}.max_pip"), overwrite=overwrite)
    ht.export(get_analysis_path(f"af_enrichment.{data_type}.{cohort}.max_pip", "tsv.bgz"))


def main(args):
    for data_type in ["exomes", "genomes"]:
        compute_af_enrichment(data_type, overwrite=args.overwrite)

        export_enrichment(data_type, "BBJ", overwrite=args.overwrite)
        export_enrichment(data_type, "FG", overwrite=args.overwrite)

    # export af_pop
    hts = []
    for data_type in ["exomes", "genomes"]:
        ht = hl.read_table(get_analysis_path(f"af_enrichment.{data_type}.fin_jpn"))
        ht = ht.select(**{"af_pop.BBJ": ht.AF_jpn, "af_pop.FG": ht.AF_fin, "af_pop.UKBB": ht.AF_nfsee})
        hts.append(ht)
    ht = hts[0].union(*hts[1:])
    ht_variants = hl.read_table(get_gwas_variants_path("BBJ_FG_UKBB_with_LOY"))
    ht = ht.join(ht_variants.select(), "inner")
    ht.export(get_analysis_path("af_pop.fin_jpn", "tsv.bgz"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
