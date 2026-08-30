__author__ = "Andrei Makhon"
__copyright__ = "Copyright 2026, Andrei Makhon"
__email__ = "andermachon@gmail.com"
__license__ = "MIT"


import tempfile
import shutil
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

ref = snakemake.input.ref
idx = snakemake.output.idx
mem_mb = snakemake.resources.get("mem_mb", 3072)

with tempfile.TemporaryDirectory() as temp_workdir:
    shell(
        " sortmerna --ref {ref}"
        " --workdir {temp_workdir}"
        " --index 1"
        " --threads {snakemake.threads}"
        " -m {mem_mb}"
        " {extra}"
        " {log}"
    )

    shutil.move(f"{temp_workdir}/idx", idx)
