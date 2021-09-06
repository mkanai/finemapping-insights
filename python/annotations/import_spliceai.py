import argparse
import hail as hl
from fm_insights.resources import get_spliceai_path
from fm_insights.utils import checkpoint_tmp, register_log


def main(args):
    mt_snv = hl.import_vcf(
        "gs://finemapping-insights/annotations/spliceai/spliceai_scores.masked.snv.hg19.vcf.bgz", min_partitions=1000
    )
    mt_indel = hl.import_vcf(
        "gs://finemapping-insights/annotations/spliceai/spliceai_scores.masked.indel.hg19.vcf.bgz", min_partitions=2000
    )

    ht = mt_snv.rows().union(mt_indel.rows())
    ht = checkpoint_tmp(ht)

    ht = ht.annotate(spliceai=ht.info.SpliceAI[0].split("\\|"))
    ht = ht.transmute(
        spliceai=hl.struct(
            symbol=ht.spliceai[1],
            DS_MAX=hl.max(hl.map(lambda x: hl.float64(x), ht.spliceai[2:6])),
            DS_AG=hl.float64(ht.spliceai[2]),
            DS_AL=hl.float64(ht.spliceai[3]),
            DS_DG=hl.float64(ht.spliceai[4]),
            DS_DL=hl.float64(ht.spliceai[5]),
            DP_AG=hl.int32(ht.spliceai[6]),
            DP_AL=hl.int32(ht.spliceai[7]),
            DP_DG=hl.int32(ht.spliceai[8]),
            DP_DL=hl.int32(ht.spliceai[9]),
        )
    )
    ht = ht.drop("rsid", "qual", "filters", "info")
    ht.describe()

    ht.write(get_spliceai_path(), overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
