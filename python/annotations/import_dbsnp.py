import argparse
from fm_insights.utils.generic import checkpoint_tmp
import hail as hl
from fm_insights.resources import get_dbsnp_path
from fm_insights.utils import register_log


def main(args):
    recoding = {
        "NC_000001.10": "1",
        "NC_000002.11": "2",
        "NC_000003.11": "3",
        "NC_000004.11": "4",
        "NC_000005.9": "5",
        "NC_000006.11": "6",
        "NC_000007.13": "7",
        "NC_000008.10": "8",
        "NC_000009.11": "9",
        "NC_000010.10": "10",
        "NC_000011.9": "11",
        "NC_000012.11": "12",
        "NC_000013.10": "13",
        "NC_000014.8": "14",
        "NC_000015.9": "15",
        "NC_000016.9": "16",
        "NC_000017.10": "17",
        "NC_000018.9": "18",
        "NC_000019.9": "19",
        "NC_000020.10": "20",
        "NC_000021.8": "21",
        "NC_000022.10": "22",
        "NC_000023.10": "X",
        "NC_000024.9": "Y",
    }
    ht = hl.import_vcf(
        "gs://finemapping-insights/annotations/dbsnp/dbsnp_b154_grch37_all_GCF_000001405.25_20200514.vcf.bgz",
        header_file="gs://finemapping-insights/annotations/dbsnp/dbsnp_b154_grch37_all_GCF_000001405.25_20200514.vcf.header",
        skip_invalid_loci=True,
        min_partitions=400,
        reference_genome="GRCh37",
        contig_recoding=recoding,
    ).rows()
    ht = ht.select("rsid")
    ht = checkpoint_tmp(ht)
    print(ht.count())

    ht = hl.split_multi(ht, permit_shuffle=True)
    ht = checkpoint_tmp(ht)
    print(ht.count())

    ht = ht.distinct()
    ht = checkpoint_tmp(ht)
    print(ht.count())

    ht = ht.checkpoint(get_dbsnp_path(), overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
