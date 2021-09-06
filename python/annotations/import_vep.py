import argparse
import hail as hl

from gnomad.utils.vep import (
    process_consequences,
    vep_or_lookup_vep,
    filter_vep_to_canonical_transcripts,
    get_most_severe_consequence_for_summary,
    CSQ_CODING_HIGH_IMPACT,
    CSQ_CODING_MEDIUM_IMPACT,
    CSQ_CODING_LOW_IMPACT,
    CSQ_NON_CODING,
)
from gnomad.resources.grch37.gnomad import public_release
from fm_insights.utils import register_log
from fm_insights.resources import get_gwas_variants_path, get_vep_annot_path

coding_high = hl.set(CSQ_CODING_HIGH_IMPACT)
coding_medium = hl.set(CSQ_CODING_MEDIUM_IMPACT)
coding_low = hl.set(CSQ_CODING_LOW_IMPACT)
non_coding = hl.set(CSQ_NON_CODING)


def annotate_consequence_category(csq_expr, annot_location="consequence_category"):
    annot_expr = {
        annot_location: hl.case()
        .when(coding_high.contains(csq_expr), "coding_high")
        .when(coding_medium.contains(csq_expr), "coding_medium")
        .when(coding_low.contains(csq_expr), "coding_low")
        .when(non_coding.contains(csq_expr), "non_coding")
        .or_missing()
    }
    return annot_expr


def main(args):
    if args.vep:
        ht = hl.read_table(get_gwas_variants_path("BBJ_FG_UKBB_with_LOY"))

        # load gnomad site info
        gnomad_ht = public_release("genomes").versions["2.1.1"].ht()
        ht = vep_or_lookup_vep(ht, gnomad_ht)

        ht = filter_vep_to_canonical_transcripts(ht)
        ht = process_consequences(ht)
        ht = get_most_severe_consequence_for_summary(ht)
        ht.describe()
        ht = ht.checkpoint(get_vep_annot_path(most_severe=False), overwrite=args.overwrite)

    ht = hl.read_table(get_vep_annot_path(most_severe=False))

    # extract most severe
    ht = ht.select(
        most_severe=hl.if_else(hl.is_defined(ht.most_severe_csq), ht.most_severe_csq, "intergenic_variant"),
        gene_most_severe=ht.vep.worst_csq_for_variant_canonical.gene_symbol,
        lof=ht.vep.worst_csq_for_variant_canonical.lof,
        hgnc_id=ht.vep.worst_csq_for_variant_canonical.hgnc_id,
        hgvsp=ht.vep.worst_csq_for_variant_canonical.hgvsp,
        transcript_id=ht.vep.worst_csq_for_variant_canonical.transcript_id,
        polyphen_prediction=ht.vep.worst_csq_for_variant_canonical.polyphen_prediction,
        polyphen_score=ht.vep.worst_csq_for_variant_canonical.polyphen_score,
        sift_prediction=ht.vep.worst_csq_for_variant_canonical.sift_prediction,
        sift_score=ht.vep.worst_csq_for_variant_canonical.sift_score,
        protein_coding=ht.protein_coding,
    )

    ht = ht.select_globals()
    ht = ht.annotate(**annotate_consequence_category(ht.most_severe))
    ht.describe()

    ht = ht.checkpoint(get_vep_annot_path(most_severe=True), overwrite=args.overwrite)

    # ht = ht.filter(hl.is_defined(ht.transcript_id) & (ht.is_canonical_vep)).key_by()
    # ht = ht.select(ht.transcript_id).key_by("transcript_id").distinct()
    # ht.export(f"{prefix}/bbj_fg_r4_ukb31063_gwas_variants.vep.transcript_id.tsv.bgz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vep", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
