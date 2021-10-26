import atexit
import datetime
import hail as hl
import os.path
import sys
import uuid
from hail.linalg.blockmatrix import BlockMatrix
from ..resources import bucket_tmp


def checkpoint_tmp(hail_obj, tmppath=f"gs://{bucket_tmp}/", tmpname=None, overwrite=True, **kwargs):
    if tmpname is None:
        tmpname = str(uuid.uuid4())
    if isinstance(hail_obj, BlockMatrix) and "force_row_major" not in kwargs:
        kwargs.update({"force_row_major": False})
    return hail_obj.checkpoint(os.path.join(tmppath, tmpname), overwrite=overwrite, **kwargs)


def register_log(name: str = None):
    if name is None and sys.argv[0] != "":
        name = os.path.basename(sys.argv[0])
        name = os.path.splitext(name)[0]
    else:
        raise ValueError("name should be specified.")

    atexit.register(
        lambda: hl.copy_log(f'gs://{bucket_tmp}/{name}_{datetime.datetime.now().strftime("%Y%m%d-%H%M%S")}.log')
    )
