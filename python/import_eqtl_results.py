import argparse
import os.path
import hail as hl

from fm_insights.utils import checkpoint_tmp, register_log, liftover
from fm_insights.resources import get_gtex_path, get_eqtl_catalogue_path, get_merged_eqtl_path


def import_gtex(ht_symbol, overwrite: bool = False):
    ht = hl.import_table(
        "gs://xfinemap/gtex/release/GTEx_49tissues_release1.SuSiE.tsv.bgz", impute=True, min_partitions=5000
    )

    # liftover to b37 (in UKBB/BBJ/FG)
    ht = ht.annotate(variant_hg38=ht.variant_hg38.replace("_b38", "").replace("_", ":"))
    ht = ht.annotate(**hl.parse_variant(ht.variant_hg38, reference_genome="GRCh38"))
    ht = liftover(ht)

    # flip betas so that effect allele is always ALT
    annot_expr = {
        k: hl.if_else(ht.alleles[1] == ht.minor_allele, ht[k], -ht[k]) for k in ["beta_marginal", "beta_posterior"]
    }
    ht = ht.annotate(**annot_expr)

    # annotate gene symbol
    ht = ht.annotate(symbol=ht_symbol[ht.gene.split("\\.")[0]].symbol)
    ht = ht.filter(hl.is_defined(ht.symbol) & (ht.symbol != ""))

    ht = ht.select(
        locus_hg38=ht.original_locus,
        alleles_hg38=ht.original_alleles,
        ref_alt_flip=ht.ref_alt_flip,
        eqtl=hl.struct(
            source="GTEx_v8",
            study="GTEx_v8",
            tissue=ht.tissue,
            method="SuSiE",
            symbol=ht.symbol,
            cs_id=ht.cs_id,
            pip=ht.pip,
            beta_marginal=ht.beta_marginal,
            se_marginal=ht.se_marginal,
            beta_posterior=ht.beta_posterior,
            sd_posterior=ht.sd_posterior,
        ),
    )
    ht.describe()
    ht = ht.checkpoint(get_gtex_path(), overwrite=args.overwrite)
    return ht


def import_eqtl_catalogue(ht_symbol, release_ver: str, overwrite: bool = False):
    # DS studies
    ht_ge = hl.import_table(
        "gs://eqtl-catalogue-alasoo/*/susie_full/*_ge.snp.txt.gz",
        force_bgz=True,
        impute=True,
        types={"chr": hl.tstr},
        min_partitions=2000,
        source_file_field="source",
    )

    # 8 new DS studies have a different header
    ht_ge2 = hl.import_table(
        "gs://eqtl-catalogue-alasoo/DS_8new_studies/*_ge.snp.txt.gz",
        force_bgz=True,
        impute=True,
        types={"chr": hl.tstr},
        min_partitions=2000,
        source_file_field="source",
    )
    ht_ge2 = ht_ge2.drop("alpha1:sd10")
    expr = {
        **{f"alpha{i}": hl.missing(hl.tfloat64) for i in range(1, 11)},
        **{f"mean{i}": hl.missing(hl.tfloat64) for i in range(1, 11)},
        **{f"sd{i}": hl.missing(hl.tfloat64) for i in range(1, 11)},
    }
    ht_ge2 = ht_ge2.annotate(**expr)

    # microarrays
    ht_microarray = hl.import_table(
        "gs://eqtl-catalogue-alasoo/*/susie_full/*_microarray.snp.txt.gz",
        force_bgz=True,
        impute=True,
        types={"chr": hl.tstr},
        min_partitions=1000,
        source_file_field="source",
    )
    ht_probeid = hl.import_table(
        "gs://finemapping-insights/eqtl_catalogue/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.bgz", impute=True
    )
    ht_probeid = ht_probeid.key_by("phenotype_id")
    ht_microarray = ht_microarray.annotate(phenotype_id=ht_probeid[ht_microarray.phenotype_id].gene_id)

    ht = hl.Table.union(ht_ge, ht_ge2, ht_microarray, unify=True)
    ht = checkpoint_tmp(ht)

    # filter variants
    ht = ht.filter((ht.pip > 0.001) | (hl.is_defined(ht.cs_index) & ~ht.low_purity))
    ht = ht.filter(~ht.variant_id.contains("<"))

    # annotate protein coding genes
    ht = ht.annotate(symbol=ht_symbol[ht.phenotype_id].symbol)
    ht = ht.filter(hl.is_defined(ht.symbol))

    # annotate study and tissue
    ht = ht.annotate(t=ht.source.split("/")[-1].split("\\."))
    ht = ht.annotate(study=ht.t[0], tissue=ht.t[1].replace("_ge$", "").replace("_microarray$", ""))
    ht = ht.drop("t", "source")
    ht = checkpoint_tmp(ht)
    ht.show()

    # liftover to b37 (in UKBB/BBJ/FG)
    ht = ht.annotate(variant=ht.variant_id.replace("_", ":"))
    ht = ht.annotate(**hl.parse_variant(ht.variant, reference_genome="GRCh38"))
    ht = liftover(ht, flip_rows=["posterior_mean"])

    ht = ht.select(
        locus_hg38=ht.original_locus,
        alleles_hg38=ht.original_alleles,
        ref_alt_flip=ht.ref_alt_flip,
        eqtl=hl.struct(
            source=f"eQTL_catalogue_{release_ver}",
            study=ht.study,
            tissue=ht.tissue,
            method="SuSiE",
            symbol=ht.symbol,
            cs_id=(hl.case().when(hl.is_missing(ht.cs_index) | ht.low_purity, -1).default(hl.int32(ht.cs_index[1:]))),
            pip=ht.pip,
            beta_marginal=hl.null(hl.tfloat64),
            se_marginal=hl.null(hl.tfloat64),
            beta_posterior=ht.posterior_mean,
            sd_posterior=ht.posterior_sd,
        ),
    )
    ht.describe()
    ht = ht.repartition(10000)
    ht = ht.checkpoint(get_eqtl_catalogue_path(release_ver=release_ver), overwrite=overwrite)
    return ht


def main(args):
    # ensgene <-> symbol conversion
    ht_symbol = hl.import_table(
        "gs://finemapping-insights/eqtl_catalogue/gene_counts_Ensembl_96_phenotype_metadata.tsv.bgz", impute=True
    )
    ht_symbol = ht_symbol.key_by(ensgene=ht_symbol.gene_id)
    ht_symbol = ht_symbol.select(symbol=ht_symbol.gene_name)
    ht_symbol = ht_symbol.cache()

    ht_gtex = import_gtex(ht_symbol, overwrite=args.overwrite) if args.gtex else hl.read_table(get_gtex_path())
    ht_eqtl_catalogue = (
        import_eqtl_catalogue(ht_symbol, release_ver=args.eqtl_catalogue_release_ver, overwrite=args.overwrite)
        if args.eqtl_catalogue
        else hl.read_table(get_eqtl_catalogue_path(release_ver=args.eqtl_catalogue_release_ver))
    )

    if args.merge:
        ht = ht_gtex.union(ht_eqtl_catalogue)
        ht = ht.checkpoint(get_merged_eqtl_path("merged"), overwrite=args.overwrite)

        ht_max_pip_per_gene = ht.group_by(ht.eqtl.symbol).aggregate(max_pip=hl.agg.max(ht.eqtl.pip))
        ht_max_pip_per_gene = ht_max_pip_per_gene.checkpoint(
            get_merged_eqtl_path("merged.max_pip_per_gene"), overwrite=args.overwrite
        )
        ht_max_pip_per_gene.export(get_merged_eqtl_path("merged.max_pip_per_gene", "tsv.bgz"))

        ht_max_pip_per_var = ht.group_by(ht.locus, ht.alleles).aggregate(max_pip=hl.agg.max(ht.eqtl.pip))
        ht_max_pip_per_var = ht_max_pip_per_var.checkpoint(
            get_merged_eqtl_path("merged.max_pip_per_var"), overwrite=args.overwrite
        )
        ht_max_pip_per_var.export(get_merged_eqtl_path("merged.max_pip_per_var", "tsv.bgz"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtex", action="store_true")
    parser.add_argument("--eqtl-catalogue", action="store_true")
    parser.add_argument("--eqtl-catalogue-release-ver", type=str, default="r4")
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
