import argparse
import hail as hl
from fm_insights.resources import get_results_path


def main(args):
    init_kwargs = {}
    if args.n_threads:
        init_kwargs["master"] = f"local[{args.n_threads}]"
    hl.init(**init_kwargs)

    prefix = f"{args.pop}.{args.trait}"

    ht = hl.read_table(get_results_path(args.pop))
    ht = ht.filter(ht.trait == args.trait)
    ht = ht.key_by()
    ht = ht.annotate(
        rsid=hl.if_else(ht.rsid.startswith("rs"), ht.rsid, ht.variant),
        af_allele2=hl.if_else(ht.allele2 == ht.minorallele, ht.maf, 1 - ht.maf),
    )

    ht_finemap = ht.filter(hl.is_defined(ht.finemap))
    ht_susie = ht.filter(hl.is_defined(ht.susie))
    common_cols = ["chromosome", "position", "allele1", "allele2", "variant", "rsid"]
    sumstats_cols = ["af_allele2", "beta_marginal", "se_marginal", "pvalue", "region", "pip"]
    fm_cols = ["pip", "cs_id", "beta_posterior", "sd_posterior"]
    finemap_expr = {f"{k}": ht_finemap.finemap[k] for k in fm_cols}
    susie_expr = {
        **{f"{k}": ht_susie.susie[k] for k in fm_cols},
        **{f"{k}{i+1}": ht_susie.susie[k][i] for k in ["alpha", "lbf_variable"] for i in range(10)},
    }

    if args.pop == "UKBB":
        ht_hwe = hl.import_table("gs://finemapping-insights/annotations/UKBB_HWE/UKBB_HWE_ldvars.tsv.bgz")
        ht_hwe = ht_hwe.key_by(**hl.parse_variant(ht_hwe.variant))
        ht_sv = hl.import_table("gs://finemapping-insights/annotations/gnomAD_SV/gnomAD_SV_ldvars.tsv.bgz")
        ht_sv = ht_sv.key_by(**hl.parse_variant(ht_sv.variant))
        ht = ht.annotate(
            LD_HWE=hl.is_defined(ht_hwe[ht.locus, ht.alleles]), LD_SV=hl.is_defined(ht_sv[ht.locus, ht.alleles])
        )
        sumstats_cols = common_cols + sumstats_cols + ["LD_HWE", "LD_SV"]
        common_cols += ["region"]
    elif args.pop == "BBJ":
        common_cols += ["af_allele2", "beta_marginal", "se_marginal", "pvalue", "region"]
        sumstats_cols = common_cols + ["pip"]

    ht.select(*sumstats_cols).export(f"{prefix}.sumstats.tsv.bgz")
    ht_finemap.select(*common_cols, **finemap_expr).export(f"{prefix}.FINEMAP.tsv.bgz")
    ht_susie.select(*common_cols, **susie_expr).export(f"{prefix}.SuSiE.tsv.bgz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pop", type=str, required=True)
    parser.add_argument("--trait", type=str, required=True)
    parser.add_argument("--n-threads", type=int)
    args = parser.parse_args()

    main(args)
