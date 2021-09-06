import argparse
import hail as hl
from fm_insights.resources import get_human_genome_dating_path
from fm_insights.utils import register_log


def main(args):
    ht = hl.import_table(
        "gs://finemapping-insights/annotations/human_genome_dating/atlas.chr*.csv.bgz",
        delimiter=",",
        missing=".",
        comment="#",
        impute=True,
        min_partitions=200,
    )
    ht = ht.transmute(
        locus=hl.locus(hl.str(ht.Chromosome), ht.Position), alleles=[ht.AlleleRef, ht.AlleleAlt], rsid=ht.VariantID
    )
    ht = ht.key_by(ht.locus, ht.alleles)

    # qual score > 0.5 and ref == ancestral allele
    ht = ht.filter((ht.QualScore_Jnt > 0.5) & (hl.is_missing(ht.AlleleAnc) | (ht.alleles[0] == ht.AlleleAnc)))

    ht = ht.collect_by_key()
    ht = ht.annotate(
        values=hl.case()
        .when(hl.len(ht.values) == 1, ht.values[0])
        .default(
            hl.rbind(
                ht.values.DataSource.index("Combined"),
                lambda idx: hl.if_else(
                    hl.is_defined(idx),
                    ht.values[idx],
                    ht.values.filter(lambda x: x.QualScore_Jnt == hl.max(ht.values.QualScore_Jnt))[0],
                ),
            )
        )
    )
    ht = ht.transmute(**ht.values)
    ht.describe()
    ht = ht.checkpoint(get_human_genome_dating_path(), overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
