import argparse
import hail as hl

from fm_insights.resources import POPS, get_ref_freq_info_path
from fm_insights.utils import register_log, liftover


def main(args):
    if args.bbj:
        # ht = hl.import_table(
        #     "gs://finemapping-insights/annotations/ref_freq_info/ALT_freq_1000Gp3v5_W01-W05.tsv.bgz",
        #     impute=True,
        #     min_partitions=100,
        # )
        # ht = ht.filter(~ht.variant.contains("<"))
        # ht = ht.key_by(**hl.parse_variant(ht.variant))
        # ht = ht.rename({"AF": "AF_reference"})

        ht = hl.import_table(
            "gs://finemapping-insights/annotations/ref_freq_info/BBJ.info.txt.bgz", impute=True, min_partitions=100,
        )
        ht = ht.filter(~ht.v.contains("<"))
        ht = ht.key_by(**hl.parse_variant(ht.v))
        ht = ht.rename({"ALT_Frq": "AF", "Rsq": "info"})
        ht = ht.select("AF", "info")

        ht = ht.checkpoint(get_ref_freq_info_path("BBJ"), overwrite=args.overwrite)

    if args.fg:
        ht = hl.import_table(
            "gs://finemapping-insights/annotations/ref_freq_info/FG_r6_info_all.txt.bgz",
            impute=True,
            min_partitions=100,
        )
        ht_frq = hl.import_table(
            "gs://finemapping-insights/annotations/ref_freq_info/FG_r6_freq.txt.bgz", impute=True, min_partitions=100
        )
        ht = ht.key_by("ID").join(ht_frq.key_by("ID"), "left")
        ht = ht.key_by(**hl.parse_variant(ht.ID.replace("_", ":"), reference_genome="GRCh38"))
        ht = liftover(ht)
        ht = ht.rename({"ALT_FREQS": "AF", "INFO_all": "info"})
        ht = ht.select("AF", "info")
        ht = ht.checkpoint(get_ref_freq_info_path("FG"), overwrite=args.overwrite)

    if args.ukbb:
        ht = hl.import_table(
            "gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_mfi_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X}_v3.txt",
            impute=True,
            no_header=True,
            min_partitions=100,
            source_file_field="source",
        )
        ht_qc = hl.read_table("gs://ukb31063/ukb31063.variant_qc.both_sexes_gwas_samples.autosomes.ht")
        ht = ht.key_by(locus=hl.locus(ht.source.split("_")[2].replace("chr", ""), ht.f2), alleles=[ht.f3, ht.f4])
        ht = ht.select(AF=ht_qc[ht.locus, ht.alleles].variant_qc.AF[1], info=ht.f7)
        ht = ht.checkpoint(get_ref_freq_info_path("UKBB"), overwrite=args.overwrite)

    if args.merge:
        hts = [hl.read_table(get_ref_freq_info_path(pop)).annotate(cohort=pop) for pop in POPS]
        ht = hts[0].union(*hts[1:])
        ht = ht.collect_by_key()

        def get_per_pop_values(field, pops=POPS):
            expr = hl.struct(
                **{
                    pop: hl.rbind(
                        ht.values.cohort.index(pop),
                        lambda idx: hl.or_missing(hl.is_defined(idx), ht.values[field][idx]),
                    )
                    for pop in pops
                }
            )
            return expr

        ht = ht.transmute(af_imp=get_per_pop_values("AF"), info=get_per_pop_values("info"))
        ht = ht.checkpoint(get_ref_freq_info_path("BBJ_FG_UKBB"), overwrite=args.overwrite)


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
