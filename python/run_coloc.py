import argparse
import hail as hl
from fm_insights.utils import checkpoint_tmp, register_log
from fm_insights.resources import get_merged_results_path, get_merged_eqtl_path


def main(args):
    ht = hl.read_table(get_merged_results_path("fm_only"))
    ht_eqtl = hl.read_table(get_merged_eqtl_path("merged"))

    ht = ht.key_by("locus", "alleles")
    ht = ht.join(ht_eqtl, "left")
    ht = ht.annotate(pip_coloc=ht.pip * ht.eqtl.pip, coloc_type="CLPP")
    ht = ht.checkpoint(get_merged_results_path("coloc"), overwrite=args.overwrite)

    # CLPP (PIP_coloc > 0.001)
    ht_clpp = ht.filter(ht.pip_coloc > 0.001)
    ht_clpp = ht_clpp.repartition(1000)
    ht_clpp = ht_clpp.checkpoint(get_merged_results_path("fm_only.coloc"), overwrite=args.overwrite)

    # max pip per variant-trait-gene-cohort quadruple
    ht_vtgc = ht_clpp.group_by("locus", "alleles", "trait", "cohort", "pip", symbol=ht_clpp.eqtl.symbol).aggregate(
        max_pip_coloc=hl.agg.max(ht_clpp.pip_coloc), eqtl=hl.agg.collect_as_set(ht_clpp.eqtl),
    )
    ht_vtgc = ht_vtgc.annotate(
        tissue_max_pip_coloc=ht_vtgc.eqtl.filter(lambda x: ht_vtgc.max_pip_coloc == ht_vtgc.pip * x.pip).tissue
    )
    ht_vtgc = ht_vtgc.checkpoint(get_merged_results_path("fm_only.coloc.max_pip_per_vtgc"), overwrite=args.overwrite)

    # max pip per variant-trait-gene triple
    ht_vtg = ht_clpp.group_by("locus", "alleles", "trait", symbol=ht_clpp.eqtl.symbol).aggregate(
        max_pip=hl.agg.max(ht_clpp.pip),
        max_pip_coloc=hl.agg.max(ht_clpp.pip_coloc),
        eqtl=hl.agg.collect_as_set(ht_clpp.eqtl),
    )
    ht_vtg = ht_vtg.annotate(
        tissue_max_pip_coloc=ht_vtg.eqtl.filter(lambda x: ht_vtg.max_pip_coloc == ht_vtg.max_pip * x.pip).tissue
    )
    ht_vtg = ht_vtg.checkpoint(get_merged_results_path("fm_only.coloc.max_pip_per_vtg"), overwrite=args.overwrite)
    ht_vtg.select("max_pip_coloc").export(get_merged_results_path("fm_only.coloc.max_pip_per_vtg", "tsv.bgz"))

    ht_vtgc = ht_vtgc.annotate(
        from_gtex=ht_vtgc.eqtl.filter(lambda x: ht_vtgc.max_pip_coloc == ht_vtgc.pip * x.pip).source.contains("GTEx_v8")
    )
    ht_vtgc = ht_vtgc.filter(ht_vtgc.max_pip_coloc > 0.1)
    ht_vtgc = checkpoint_tmp(ht_vtgc)

    ht_vtgc_ukbb = ht_vtgc.filter(ht_vtgc.cohort == "UKBB")
    ht_vtg_ukbb = ht_vtgc_ukbb.group_by("locus", "alleles", "trait", symbol=ht_vtgc_ukbb.eqtl.symbol).aggregate(
        max_pip=hl.agg.max(ht_vtgc_ukbb.pip),
        max_pip_coloc=hl.agg.max(ht_vtgc_ukbb.max_pip_coloc),
        eqtl=hl.agg.collect_as_set(ht_vtgc_ukbb.eqtl),
    )
    ht_vtg_ukbb = ht_vtg_ukbb.checkpoint(
        get_merged_results_path("fm_only.coloc.max_pip_per_vtg.ukbb"), overwrite=args.overwrite
    )

    ht_vtgc_bbj_fg = ht_vtgc.filter(ht_vtgc.cohort != "UKBB")
    ht_vtg_bbj_fg = ht_vtgc_bbj_fg.group_by("locus", "alleles", "trait", symbol=ht_vtgc_bbj_fg.eqtl.symbol).aggregate(
        max_pip=hl.agg.max(ht_vtgc_bbj_fg.pip),
        max_pip_coloc=hl.agg.max(ht_vtgc_bbj_fg.max_pip_coloc),
        eqtl=hl.agg.collect_as_set(ht_vtgc_bbj_fg.eqtl),
    )
    ht_vtg_bbj_fg = ht_vtg_bbj_fg.checkpoint(
        get_merged_results_path("fm_only.coloc.max_pip_per_vtg.bbj_fg"), overwrite=args.overwrite
    )

    print(ht_vtg_ukbb.count())
    print(ht_vtg_bbj_fg.count())

    # max pip per variant-trait-cohort triple
    ht_vtc = ht_clpp.group_by("locus", "alleles", "trait", "cohort").aggregate(
        max_pip=hl.agg.max(ht_clpp.pip),
        max_pip_coloc=hl.agg.max(ht_clpp.pip_coloc),
        eqtl=hl.agg.collect_as_set(ht_clpp.eqtl),
    )
    ht_vtc = ht_vtc.annotate(
        tissue_max_pip_coloc=ht_vtc.eqtl.filter(lambda x: ht_vtc.max_pip_coloc == ht_vtc.max_pip * x.pip).tissue
    )
    ht_vtc = ht_vtc.checkpoint(get_merged_results_path("fm_only.coloc.max_pip_per_vtc"), overwrite=args.overwrite)
    ht_vtc.select("max_pip_coloc").export(get_merged_results_path("fm_only.coloc.max_pip_per_vtc", "tsv.bgz"))

    # max_pip per gene-trait pair
    ht_gt = ht_clpp.group_by("trait", symbol=ht_clpp.eqtl.symbol).aggregate(
        max_pip=hl.agg.max(ht_clpp.pip),
        max_pip_coloc=hl.agg.max(ht_clpp.pip_coloc),
        eqtl=hl.agg.collect_as_set(ht_clpp.eqtl),
    )
    ht_gt = ht_gt.checkpoint(get_merged_results_path("fm_only.coloc.max_pip_per_gt"), overwrite=args.overwrite)
    ht_gt.select("max_pip", "max_pip_coloc").export(get_merged_results_path("fm_only.coloc.max_pip_per_gt", "tsv.bgz"))

    # max_pip per variant-trait pair
    ht_vt = ht_vtc.group_by("locus", "alleles", "trait").aggregate(
        max_pip=hl.agg.max(ht_vtc.pip),
        max_pip_coloc=hl.agg.max(ht_vtc.max_pip_coloc),
        eqtl=hl.agg.explode(lambda x: hl.agg.collect_as_set(x), ht_vtc.eqtl),
    )
    ht_vt = ht_vt.annotate(max_pip_coloc_eqtl=ht_vt.eqtl.filter(lambda x: ht_vt.max_pip_coloc == ht_vt.max_pip * x.pip))
    ht_vt = ht_vt.annotate(
        tissue_max_pip_coloc=ht_vt.max_pip_coloc_eqtl.tissue,
        eqtl_study=ht_vt.max_pip_coloc_eqtl.study,
        eqtl_symbol=ht_vt.max_pip_coloc_eqtl.symbol,
        eqtl_pip=ht_vt.max_pip_coloc_eqtl.pip,
    )
    ht_vt = ht_vt.checkpoint(get_merged_results_path("fm_only.coloc.max_pip_per_vt"), overwrite=args.overwrite)
    ht_vt.select("max_pip", "max_pip_coloc").export(get_merged_results_path("fm_only.coloc.max_pip_per_vt", "tsv.bgz"))

    # max_pip per variant
    ht_var = ht_vt.group_by("locus", "alleles").aggregate(
        max_pip=hl.agg.max(ht_vt.max_pip),
        max_pip_coloc=hl.agg.max(ht_vt.max_pip_coloc),
        eqtl=hl.agg.explode(lambda x: hl.agg.collect_as_set(x), ht_vt.eqtl),
        trait_pip=hl.agg.collect_as_set(hl.struct(trait=ht_vt.trait, max_pip=ht_vt.max_pip)),
    )
    ht_var = ht_var.annotate(
        tissue_max_pip_coloc=ht_var.eqtl.filter(lambda x: ht_var.max_pip_coloc == ht_var.max_pip * x.pip).tissue,
        trait_max_pip_coloc=ht_var.trait_pip.filter(
            lambda x: ht_var.max_pip_coloc == x.max_pip * hl.max(ht_var.eqtl.pip)
        ).trait,
    )
    ht_var = ht_var.drop("trait_pip")
    ht_var = ht_var.checkpoint(get_merged_results_path("fm_only.coloc.max_pip_per_var"), overwrite=args.overwrite)

    # Regional coloc
    ht_rclpp = ht.filter(ht.susie.cs_id > 0)
    ht_rclpp = (
        ht_rclpp.group_by(
            "cohort",
            "trait",
            "region",
            ht_rclpp.susie.cs_id,
            ht_rclpp.eqtl.source,
            ht_rclpp.eqtl.study,
            ht_rclpp.eqtl.tissue,
            ht_rclpp.eqtl.symbol,
        )
        .aggregate(pip_coloc=hl.agg.sum(ht_rclpp.pip_coloc))
        .annotate(coloc_type="RCLPP")
    )
    ht_rclpp = ht_rclpp.transmute(
        eqtl=hl.struct(source=ht_rclpp.source, study=ht_rclpp.study, tissue=ht_rclpp.tissue, symbol=ht_rclpp.symbol)
    )
    ht_rclpp = ht_rclpp.repartition(1000)
    ht_rclpp = ht_rclpp.checkpoint(get_merged_results_path("fm_only.rcoloc"), overwrite=args.overwrite)

    # Regional coloc CSM
    ht_rclpp = ht.filter(ht.susie.cs_id > 0)
    ht_rclpp = ht_rclpp.key_by("cohort", "trait", "region", ht_rclpp.susie.cs_id)

    # merge CSM ID
    ht_csm = hl.import_table(get_merged_results_path("fm_only.csm_id", "tsv.bgz"), impute=True, min_partitions=200)
    ht_csm = ht_csm.transmute(susie=hl.struct(cs_id=ht_csm["susie.cs_id"]))
    ht_csm = ht_csm.key_by("cohort", "trait", "region", ht_csm.susie.cs_id)
    ht_csm = ht_csm.select("csm_id")
    ht_csm = ht_csm.distinct()

    ht_rclpp = ht_rclpp.join(ht_csm, "left")
    ht_rclpp = checkpoint_tmp(ht_rclpp)

    ht_rclpp = (
        ht_rclpp.group_by(
            "cohort",
            "trait",
            "csm_id",
            ht_rclpp.eqtl.source,
            ht_rclpp.eqtl.study,
            ht_rclpp.eqtl.tissue,
            ht_rclpp.eqtl.symbol,
        )
        .aggregate(pip_coloc=hl.agg.sum(ht_rclpp.pip_coloc))
        .annotate(coloc_type="RCLPP")
    )
    ht_rclpp = ht_rclpp.transmute(
        eqtl=hl.struct(source=ht_rclpp.source, study=ht_rclpp.study, tissue=ht_rclpp.tissue, symbol=ht_rclpp.symbol)
    )
    ht_rclpp = ht_rclpp.repartition(1000)
    ht_rclpp = ht_rclpp.checkpoint(get_merged_results_path("fm_only.rcoloc.csm"), overwrite=args.overwrite)

    # max_pip per CSM-trait pair
    ht4 = ht_rclpp.group_by("trait", "csm_id").aggregate(
        max_pip_coloc=hl.agg.max(ht_rclpp.pip_coloc), eqtl=hl.agg.collect_as_set(ht_rclpp.eqtl),
    )
    ht4 = ht4.checkpoint(get_merged_results_path("fm_only.rcoloc.csm.max_pip_per_ct"), overwrite=args.overwrite)
    ht4.select("max_pip_coloc").export(get_merged_results_path("fm_only.rcoloc.csm.max_pip_per_ct", "tsv.bgz"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
