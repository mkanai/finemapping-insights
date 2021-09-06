import argparse
import hail as hl
from fm_insights.resources import get_ems_path
from fm_insights.utils import checkpoint_tmp, liftover, register_log


def main(args):
    ht_symbol = hl.import_table(
        "gs://finemapping-insights/eqtl_catalogue/gene_counts_Ensembl_96_phenotype_metadata.tsv.bgz", impute=True
    )
    ht_symbol = ht_symbol.key_by(ensgene=ht_symbol.gene_id)
    ht_symbol = ht_symbol.select(symbol=ht_symbol.gene_name)
    ht_symbol = ht_symbol.cache()

    ht = hl.read_table(f"gs://expression-modifier-score/public/ht/ems_{args.tissue}.ht")
    ht = ht.annotate(**hl.parse_variant(ht.v.replace("_b38", "").replace("_", ":"), reference_genome="GRCh38"))
    ht = ht.annotate(ensgene=ht.g.split("\\.")[0])
    ht = ht.annotate(symbol=ht_symbol[ht.ensgene].symbol)
    ht = ht.select("locus", "alleles", "ems_normalized", "symbol")
    ht = ht.group_by("locus", "alleles").aggregate(
        max_ems_normalized=hl.agg.max(ht.ems_normalized),
        x=hl.agg.collect(hl.struct(symbol=ht.symbol, ems_normalized=ht.ems_normalized)),
    )
    ht = ht.annotate(gene_max_ems_normalized=ht.x.filter(lambda x: ht.max_ems_normalized == x.ems_normalized).symbol[0])
    ht = ht.drop("x")
    ht.describe()
    ht = checkpoint_tmp(ht)

    ht = liftover(ht)
    ht.write(get_ems_path(args.tissue), overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--tissue", type=str, default="Whole_Blood")
    args = parser.parse_args()

    register_log()

    main(args)
