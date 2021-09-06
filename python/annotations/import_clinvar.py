import argparse
import hail as hl
from gnomad.utils.vep import (
    vep_or_lookup_vep,
    filter_vep_to_canonical_transcripts,
    process_consequences,
    get_most_severe_consequence_for_summary,
)
from gnomad.resources.grch37.gnomad import public_release
from fm_insights.resources import get_clinvar_path
from fm_insights.utils import register_log


def main(args):
    ht = hl.import_vcf(
        "gs://finemapping-insights/annotations/clinvar/clinvar_20210131.vcf.bgz",
        skip_invalid_loci=True,
        min_partitions=100,
        reference_genome="GRCh37",
    ).rows()
    ht = ht.filter(hl.len(ht.alleles) > 1)  # Get around problematic single entry in alleles array in the clinvar vcf

    # vep
    gnomad_ht = public_release("genomes").versions["2.1.1"].ht()
    ht = vep_or_lookup_vep(ht, gnomad_ht)
    ht = filter_vep_to_canonical_transcripts(ht)
    ht = process_consequences(ht)
    ht = get_most_severe_consequence_for_summary(ht)
    ht.describe()
    ht = ht.checkpoint(get_clinvar_path(), overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
