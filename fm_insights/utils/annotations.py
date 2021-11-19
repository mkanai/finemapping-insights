import hail as hl
import os.path
from ..resources.annotations import get_vep_annot_path, get_baseline_annotations


def import_annot_bed(path, skip_invalid_intervals=True, reference_genome="GRCh37"):
    bed = hl.import_bed(
        path,
        skip_invalid_intervals=skip_invalid_intervals,
        reference_genome=reference_genome,
    )
    if "target" in bed.row:
        bed = bed.annotate(target=hl.float64(bed.target))
    else:
        bed = bed.annotate(target=hl.int32(True))
    return bed.cache()


def annotate_bed(ht, bed_files=None, reference_genome="GRCh37"):
    if bed_files is None:
        bed_files = get_baseline_annotations()

    bed_names = [os.path.basename(x).replace(".bed", "") for x in bed_files]
    annot_expr = {
        fname: hl.rbind(
            import_annot_bed(path, reference_genome=reference_genome)[ht.locus].target,
            lambda x: hl.or_else(x, hl.missing(hl.tfloat64))
            if x.dtype == hl.tfloat64
            else hl.or_else(x, 0),
        )
        for fname, path in zip(bed_names, bed_files)
    }
    ht = ht.annotate(**annot_expr)
    return ht
