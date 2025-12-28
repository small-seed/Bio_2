# src/pipeline.py
from __future__ import annotations

from pathlib import Path

from qc import QC_run
from align import STAR_alignment


# --- Absolute project root (Windows). utils.normalize_project_root will convert if running in WSL.
PROJECT_ROOT = r"D:\Code\Bio_2\NGS-rnaseq"

# --- Config
SAMPLES_CSV = Path("config/samples.csv")

# --- Reference paths
GENOME_FASTA = Path("reference/genome/genome.fa")
ANNOT_GTF = Path("reference/annotation/genes.gtf")
STAR_INDEX_DIR = Path("reference/star_index/BDGP6.54_115")


def run_pipeline(
    *,                             # make sure everything must be passed as a keyword argument
    qc_threads: int = 2,           # QC set up
    qc_jobs: int = 2,              # QC set up
    star_threads: int = 8,         # STAR set up
    star_jobs: int = 1,            # STAR set up
    sjdb_overhang: int = 49,       # STAR set up: read_length - 1 
    build_star_index: bool = True, # STAR set up: True first time; False after index exists
    force_index: bool = False,     # STAR set up
    force: bool = False,           # STAR set up
) -> None:
    # 1) QC
    QC_run(
        project_root=PROJECT_ROOT,
        samples_csv=SAMPLES_CSV,
        threads=qc_threads,
        jobs=qc_jobs,
        force=force,
    )

    # 2) Alignment
    STAR_alignment(
        project_root=PROJECT_ROOT,
        samples_csv=SAMPLES_CSV,
        genome_fasta=GENOME_FASTA,
        annotation_gtf=ANNOT_GTF,
        star_index_dir=STAR_INDEX_DIR,
        threads=star_threads,
        jobs=star_jobs,
        build_index=build_star_index,
        force_index=force_index,
        sjdb_overhang=sjdb_overhang,
        force=force,
    )
