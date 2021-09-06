from .generic import bucket


def get_raw_sumstats_path(pop, trait="*"):
    raw_sumstats_path = {
        "BBJ": f"gs://{bucket}/bbj/sumstats_formatted/BBJ.{trait}.sumstats.txt.bgz",
        "UKBB": f"gs://{bucket}/ukbb/sumstats_formatted/UKBB.{trait}.sumstats.txt.bgz",
        "FG": f"gs://{bucket}/fg_r6/sumstats_formatted/FG.{trait}.sumstats.txt.bgz",
    }
    return raw_sumstats_path[pop]


def get_raw_susie_path(pop, trait="*"):
    raw_susie_path = {
        "BBJ": f"gs://{bucket}/bbj/results/BBJ.{trait}.SuSiE.snp.bgz",
        "UKBB": f"gs://{bucket}/ukbb/results/UKBB.{trait}.SuSiE.snp.bgz",
        "FG": f"gs://{bucket}/fg_r6/results/FG.{trait}.SUSIE.snp.bgz",
    }
    return raw_susie_path[pop]


def get_raw_finemap_path(pop, trait="*"):
    raw_finemap_path = {
        "BBJ": f"gs://{bucket}/bbj/results/BBJ.{trait}.FINEMAP.snp.bgz",
        "UKBB": f"gs://{bucket}/ukbb/results/UKBB.{trait}.FINEMAP.snp.bgz",
        "FG": f"gs://{bucket}/fg_r6/results/FG.{trait}.FINEMAP.snp.bgz",
    }
    return raw_finemap_path[pop]


def get_latest_release_ver(pop):
    latest_release_ver = {"BBJ": "1.0", "UKBB": "2.0", "FG": "1.0"}
    return latest_release_ver[pop]


def get_latest_analysis_ver():
    latest_analysis_ver = "1.1"
    return latest_analysis_ver


def get_results_path(pop, release_ver=None):
    if release_ver is None:
        release_ver = get_latest_release_ver(pop)
    results_path = {
        "BBJ": f"gs://{bucket}/bbj/release/release{release_ver}/BBJ_release{release_ver}.ht",
        "UKBB": f"gs://{bucket}/ukbb/release/release{release_ver}/UKBB_release{release_ver}.ht",
        "FG": f"gs://{bucket}/fg_r6/release/release{release_ver}/FG_release{release_ver}.ht",
    }
    return results_path[pop]


def get_analysis_path(filename: str, extension: str = "ht", analysis_ver: str = None):
    if analysis_ver is None:
        analysis_ver = get_latest_analysis_ver()

    path = f"gs://{bucket}/analysis/{analysis_ver}/{filename}.{extension}"
    return path


def get_merged_results_path(suffix: str, extension: str = "ht", analysis_ver: str = None):
    return get_analysis_path(f"BBJ_FG_UKBB.{suffix}", extension, analysis_ver)


def get_gtex_path(release_ver=None):
    if release_ver is None:
        release_ver = "v8"
    return f"gs://{bucket}/gtex/gtex.{release_ver}.pip001.in_cs.ht"


def get_eqtl_catalogue_path(release_ver=None):
    if release_ver is None:
        release_ver = "r2"
    return f"gs://{bucket}/eqtl_catalogue/eqtl_catalogue.{release_ver}.pip001.in_cs.ht"


def get_merged_eqtl_path(suffix: str, extension: str = "ht", analysis_ver: str = None):
    if analysis_ver is None:
        analysis_ver = get_latest_analysis_ver()

    path = f"gs://{bucket}/analysis/{analysis_ver}/eQTL.{suffix}.{extension}"
    return path


def get_ld_txt_path(pop, lead_variant):
    lead_variant = lead_variant.replace(":", ".")
    path = {
        "BBJ": f"gs://{bucket}/bbj/ld/BBJ.{lead_variant}.ld.txt",
        "UKBB": f"gs://{bucket}/ukbb/ld/UKBB.{lead_variant}.ld.txt",
        "FG": f"gs://{bucket}/fg_r6/ld/FG.{lead_variant}.ld.txt",
    }
    return path[pop]


def get_ld_matrix_path(pop, lead_variant, window=500000):
    lead_variant = lead_variant.replace(":", ".")
    path = {
        "BBJ": f"gs://{bucket}/bbj/ld/BBJ.{lead_variant}.w{window}.ld.bgz",
        "UKBB": f"gs://{bucket}/ukbb/ld/UKBB.{lead_variant}.w{window}.ld.bgz",
        "FG": f"gs://{bucket}/fg_r6/ld/FG.{lead_variant}.w{window}.ld.bgz",
    }
    return path[pop]


def get_locus_txt_path(pop, trait, lead_variant):
    lead_variant = lead_variant.replace(":", ".")
    path = {
        "BBJ": f"gs://{bucket}/bbj/locus/BBJ.{trait}.{lead_variant}.locus.txt.bgz",
        "UKBB": f"gs://{bucket}/ukbb/locus/UKBB.{trait}.{lead_variant}.locus.txt.bgz",
        "FG": f"gs://{bucket}/fg_r6/locus/FG.{trait}.{lead_variant}.locus.txt.bgz",
    }
    return path[pop]


def get_gtex_locus_txt_path(tissue, symbol, lead_variant):
    lead_variant = lead_variant.replace(":", ".")
    return f"gs://{bucket}/gtex/locus/GTEx_v8.{tissue}.{symbol}.{lead_variant}.locus.txt.bgz"
