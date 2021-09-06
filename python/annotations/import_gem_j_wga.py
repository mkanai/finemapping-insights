import argparse
import hail as hl
from gnomad.utils.vep import (
    vep_or_lookup_vep,
    filter_vep_to_canonical_transcripts,
    process_consequences,
    get_most_severe_consequence_for_summary,
)
from gnomad.resources.grch37.gnomad import public_release
from fm_insights.resources import get_gem_j_wga_path
from fm_insights.utils import register_log


def main(args):
    if args.import_vcf:
        ht = hl.import_vcf(
            "gs://ukbb-hail-tmp/BBJ_RIKEN_TMM_20200309_*.vqsr_99.5_99.0.region_flag.genotype_filtering.hwe_flag.annotated.vcf.gz",
            force_bgz=True,
            drop_samples=True,
            n_partitions=1000,
        )
        split_ds = hl.split_multi(ht.rows())
        split_ds = split_ds.annotate(
            info=split_ds.info.annotate(**{key: split_ds.info[key][split_ds.a_index - 1] for key in ["AC", "AF"]})
        )
        ht = split_ds.checkpoint(get_gem_j_wga_path(), overwrite=args.ovewrtite)

    if args.vep:
        ht = hl.read_table(get_gem_j_wga_path())
        gnomad_ht = public_release("genomes").versions["2.1.1"].ht()

        ht = vep_or_lookup_vep(ht, gnomad_ht)
        ht = filter_vep_to_canonical_transcripts(ht)
        ht = process_consequences(ht)
        ht = get_most_severe_consequence_for_summary(ht)
        ht = ht.checkpoint(get_gem_j_wga_path(vep=True), overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--import-vcf", action="store_true")
    parser.add_argument("--vep", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
