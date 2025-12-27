#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

from utils import (
    SampleRow,
    which_or_fail,
    run_cmd,
    normalize_project_root,
    parse_samples_csv,
)

# ---- Absolute main folder ----
DEFAULT_PROJECT_ROOT = r"D:\Code\Bio_2\NGS-RNASEQ"


def expected_fastqc_outputs(outdir: Path, fastq: Path) -> Tuple[Path, Path]:
    name = fastq.name
    if name.endswith(".gz"):
        name = name[:-3]
    for ext in (".fastq", ".fq"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    return outdir / f"{name}_fastqc.html", outdir / f"{name}_fastqc.zip"


def run_fastqc_one(
    fastq: Path,
    outdir: Path,
    threads: int,
    force: bool,
    logs_dir: Path,
) -> None:
    html, zipf = expected_fastqc_outputs(outdir, fastq)
    if (html.exists() and zipf.exists()) and not force:
        return

    log_file = logs_dir / f"fastqc_{fastq.stem}.log"
    cmd = ["fastqc", "--threads", str(threads), "--outdir", str(outdir), str(fastq)]
    run_cmd(cmd, log_file)


def run_multiqc(
    fastqc_dir: Path,
    multiqc_dir: Path,
    force: bool,
    logs_dir: Path,
) -> None:
    report = multiqc_dir / "multiqc_report.html"
    if report.exists() and not force:
        return

    log_file = logs_dir / "multiqc.log"
    cmd = ["multiqc", str(fastqc_dir), "--outdir", str(multiqc_dir)]
    if force:
        cmd.append("--force")
    run_cmd(cmd, log_file)


def QC_run() -> None:
    parser = argparse.ArgumentParser(
        description="QC stage using samples.csv (FastQC + MultiQC)"
    )
    parser.add_argument(
        "--project-root",
        type=str,
        default=DEFAULT_PROJECT_ROOT,
        help=r"Project root (default: D:\Code\Bio_2\NGS-rnaseq). Auto-converts to /mnt/d/... in WSL.",
    )
    parser.add_argument(
        "--samples",
        type=Path,
        default=Path("config/samples.csv"),
        help="Path to samples.csv relative to project root (default: config/samples.csv)",
    )
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    which_or_fail("fastqc")
    which_or_fail("multiqc")

    project_root = normalize_project_root(args.project_root)

    # --- Create results directories ---
    results_dir = project_root / "results"
    fastqc_dir = results_dir / "fastQC"
    multiqc_dir = results_dir / "multiQC"
    logs_dir = results_dir / "qc_logs"

    fastqc_dir.mkdir(parents=True, exist_ok=True)
    multiqc_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    samples_csv = (
        (project_root / args.samples).resolve()
        if not args.samples.is_absolute()
        else args.samples
    )
    samples = parse_samples_csv(samples_csv, project_root)

    fastqs: List[Path] = []
    for s in samples:
        fastqs.append(s.read1)
        if s.read2 is not None:
            fastqs.append(s.read2)

    errors: List[str] = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {
            ex.submit(
                run_fastqc_one,
                fq,
                fastqc_dir,
                args.threads,
                args.force,
                logs_dir,
            ): fq
            for fq in fastqs
        }
        for fut in as_completed(futs):
            fq = futs[fut]
            try:
                fut.result()
            except Exception as e:
                errors.append(f"{fq}: {e}")

    if errors:
        raise RuntimeError(
            "FastQC failures:\n" + "\n".join("- " + x for x in errors)
        )

    run_multiqc(fastqc_dir, multiqc_dir, args.force, logs_dir)

    print("QC finished.")
    print(f"Project root:   {project_root}")
    print(f"FastQC dir:     {fastqc_dir}")
    print(f"MultiQC report: {multiqc_dir / 'multiqc_report.html'}")
    print(f"Logs dir:       {logs_dir}")
