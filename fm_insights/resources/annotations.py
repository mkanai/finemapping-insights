import os.path
from .generic import bucket


def get_vep_annot_path(most_severe=True):
    if most_severe:
        return f"gs://{bucket}/annotations/vep/BBJ_FG_UKBB_with_LOY.vep.most_severe.ht"
    return f"gs://{bucket}/annotations/vep/BBJ_FG_UKBB_with_LOY.vep.ht"


def get_gwas_variants_path(pop: str):
    return f"gs://{bucket}/annotations/gwas_variants/{pop}.gwas_variants.ht"


def get_ref_freq_info_path(pop: str):
    return f"gs://{bucket}/annotations/ref_freq_info/{pop}.ref_freq_info.ht"


def get_gem_j_wga_path(vep: bool = False):
    if vep:
        return f"gs://{bucket}/annotations/gem_j_wga/gem_j_wga.BBJ_RIKEN_TMM_20200309.sites.vep.ht"
    return "gs://xfinemap/annotation/gem_j_wga/gem_j_wga.BBJ_RIKEN_TMM_20200309.sites.ht"


def get_clinvar_path():
    return f"gs://{bucket}/annotations/clinvar/clinvar_20210131.vep.ht"


def get_dbsnp_path():
    return f"gs://{bucket}/annotations/dbsnp/dbsnp_b154_grch37_all_GCF_000001405.25_20200514.vep.ht"


def get_human_genome_dating_path():
    return f"gs://{bucket}/annotations/human_genome_dating/human_genome_dating.ht"


def get_gtex_sqtl_signifpairs_path():
    return f"gs://{bucket}/annotations/GTEx_Analysis_v8_sQTL/GTEx.v8.sqtl_signifpairs.ht"


def get_spliceai_path():
    return f"gs://{bucket}/annotations/spliceai/spliceai_scores.masked.snv.indel.hg19.ht"


def get_ems_path(tissue: str):
    return f"gs://{bucket}/annotations/ems/max_ems_normalized.{tissue}.ht"


def get_baseline_annotations(basename=False):
    bed_files = [
        f"gs://{bucket}/annotations/baselineLD_v2.2/Ancient_Sequence_Age_Human_Enhancer.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Ancient_Sequence_Age_Human_Promoter.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/BivFlnk.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/CTCF_Hoffman.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Coding_UCSC.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Conserved_LindbladToh.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Conserved_Mammal_phastCons46way.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Conserved_Primate_phastCons46way.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Conserved_Vertebrate_phastCons46way.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/DGF_ENCODE.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/DHS_Trynka.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/DHS_peaks_Trynka.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Enhancer_Andersson.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Enhancer_Hoffman.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/FetalDHS_Trynka.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/H3K27ac_Hnisz.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/H3K27ac_PGC2.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/H3K4me1_Trynka.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/H3K4me1_peaks_Trynka.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/H3K4me3_Trynka.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/H3K4me3_peaks_Trynka.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/H3K9ac_Trynka.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/H3K9ac_peaks_Trynka.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Human_Enhancer_Villar.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Human_Promoter_Villar.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Human_Promoter_Villar_ExAC.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Intron_UCSC.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/PromoterFlanking_Hoffman.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Promoter_UCSC.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Repressed_Hoffman.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/SuperEnhancer_Hnisz.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/TFBS_ENCODE.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/TSS_Hoffman.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Transcribed_Hoffman.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/UTR_3_UCSC.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/UTR_5_UCSC.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Vahedi_Tcell_SE.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Vahedi_Tcell_TE.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/WeakEnhancer_Hoffman.bed",
        # continous
        f"gs://{bucket}/annotations/baselineLD_v2.2/ASMC.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Backgrd_Selection_Stat.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/BLUEPRINT_FE_META_TISSUE_DNAMETH_MaxCPP.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/BLUEPRINT_FE_META_TISSUE_H3K27ac_MaxCPP.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/BLUEPRINT_FE_META_TISSUE_H3K4me1_MaxCPP.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/GTEx_FE_META_TISSUE_GE_MaxCPP.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/Human_Enhancer_Villar_Species_Enhancer_Count.bed",
        f"gs://{bucket}/annotations/baselineLD_v2.2/alleleage.bed",
        # Ulrisch
        f"gs://{bucket}/annotations/Ulirsch_v1.0/DHSmerged_Ulirsch.bed",
        f"gs://{bucket}/annotations/Ulirsch_v1.0/Roadmap_H3K27ac_Ulirsch.bed",
        f"gs://{bucket}/annotations/Ulirsch_v1.0/CA_H3K27ac_Ulirsch.bed",
    ]

    if basename:
        bed_names = [os.path.basename(x).replace(".bed", "") for x in bed_files]
        return bed_names

    return bed_files
