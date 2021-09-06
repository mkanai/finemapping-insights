import argparse
import hail as hl
from fm_insights.resources import get_gtex_sqtl_signifpairs_path
from fm_insights.utils import liftover, register_log


def main(args):
    ht = hl.import_table(
        "gs://finemapping-insights/annotations/GTEx_Analysis_v8_sQTL/*.v8.sqtl_signifpairs.txt.bgz",
        min_partitions=100,
        source_file_field="source",
    )
    ht = ht.transmute(tissue=ht.source.split(r"\/")[-1].split(r"\.")[0])
    ht = ht.key_by(**hl.parse_variant(ht.variant_id.replace("_b38", "").replace("_", ":"), reference_genome="GRCh38"))
    ht = ht.collect_by_key()
    ht = liftover(ht)
    ht.describe()

    ht.write(get_gtex_sqtl_signifpairs_path(), overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
