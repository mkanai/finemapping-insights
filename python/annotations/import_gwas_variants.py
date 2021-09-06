import argparse
import hail as hl

from fm_insights.resources import POPS, get_gwas_variants_path, get_raw_sumstats_path
from fm_insights.utils import register_log


def main(args):
    if args.bbj:
        ht = hl.import_table(
            "gs://xfinemap/annotation/gwas_variants/bbj_gwas_variants.txt.bgz", impute=True, min_partitions=1000
        )
        ht = ht.filter(~ht.variant.contains("<"))
        ht = ht.key_by(**hl.parse_variant(ht.variant))
        ht = ht.checkpoint(get_gwas_variants_path("BBJ"), overwrite=args.overwrite)

    if args.fg:
        ht = hl.import_table(get_raw_sumstats_path("FG", trait="BMI_IRN"), impute=True, min_partitions=1000)
        ht = ht.filter(hl.is_defined(ht.variant_b37))
        ht = ht.key_by(**hl.parse_variant(ht.variant_b37)).select()
        ht = ht.checkpoint(get_gwas_variants_path("FG"), overwrite=args.overwrite)

    if args.ukbb:
        ht = hl.import_table(
            "gs://xfinemap/annotation/gwas_variants/ukb31063_gwas_variants.tsv.bgz", impute=True, min_partitions=1000
        )
        ht = ht.key_by(**hl.parse_variant(ht.variant))
        ht = ht.checkpoint(get_gwas_variants_path("UKBB"), overwrite=args.overwrite)

        # add LOY variants
        ht_loy = hl.import_table(get_raw_sumstats_path("UKBB", trait="LOY"), impute=True, min_partitions=100)
        ht_loy = ht_loy.key_by(**hl.parse_variant(ht_loy.variant)).select()
        ht = ht.select().union(ht_loy)
        ht = ht.distinct()
        ht = ht.checkpoint(get_gwas_variants_path("UKBB_with_LOY"), overwrite=args.overwrite)

    if args.merge:
        hts = [hl.read_table(get_gwas_variants_path(pop)).select() for pop in POPS]
        ht = hts[0].union(*hts[1:])
        ht = ht.distinct()
        ht = ht.checkpoint(get_gwas_variants_path("BBJ_FG_UKBB"), overwrite=args.overwrite)

        ht_loy = hl.read_table(get_gwas_variants_path("UKBB_with_LOY"))
        ht = ht.union(ht_loy)
        ht = ht.distinct()
        ht = ht.checkpoint(get_gwas_variants_path("BBJ_FG_UKBB_with_LOY"), overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bbj", action="store_true")
    parser.add_argument("--fg", action="store_true")
    parser.add_argument("--ukbb", action="store_true")
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
