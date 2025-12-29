# src/pipeline.py
from __future__ import annotations

from pathlib import Path

from qc import QC_run
from align import STAR_alignment
from count import featureCounts_run

PROJECT_ROOT = r"D:\Code\Bio_2\NGS-rnaseq"

SAMPLES_CSV = Path("config/samples.csv")

GENOME_FASTA = Path("reference/genome/genome.fa")
ANNOT_GTF = Path("reference/annotation/genes.gtf")
STAR_INDEX_DIR = Path("reference/star_index/BDGP6.54_115")


def run_pipeline(
    *,
    qc_threads: int = 2,
    qc_jobs: int = 2,
    star_threads: int = 8,
    star_jobs: int = 1,
    sjdb_overhang: int = 49,
    build_star_index: bool = True,
    force_index: bool = False,
    count_threads: int = 8,
    count_stranded: int = 0,   # 0/1/2
    force: bool = False,
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

    # 3) Counting
    featureCounts_run(
        project_root=PROJECT_ROOT,
        samples_csv=SAMPLES_CSV,
        annotation_gtf=ANNOT_GTF,
        threads=count_threads,
        stranded=count_stranded,
        force=force,
    )
