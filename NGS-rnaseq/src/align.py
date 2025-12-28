from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List

from utils import (
    SampleRow,
    which_or_fail,
    run_cmd,
    normalize_project_root,
    parse_samples_csv,
)

# Keep align.py generic: full-pipeline.py can override these via CLI args
DEFAULT_PROJECT_ROOT = r"D:\Code\Bio_2\NGS-RNASEQ"

# Default "stable" reference paths 
DEFAULT_GENOME_FASTA = Path("reference/genome.fa")
DEFAULT_ANNOT_GTF = Path("reference/genes.gtf")
DEFAULT_STAR_INDEX = Path("reference/star_index/BDGP6.54_115")


def ensure_star_index(
    star_index_dir: Path,
    genome_fasta: Path,
    annotation_gtf: Path,
    threads: int,
    force_index: bool,
    logs_dir: Path,
    sjdb_overhang: int,
) -> None:
    star_index_dir.mkdir(parents=True, exist_ok=True)

    required_any = [star_index_dir / "SA", star_index_dir / "Genome"]
    index_exists = all(p.exists() for p in required_any)
    if index_exists and not force_index:
        return

    if not genome_fasta.exists():
        raise FileNotFoundError(f"Genome FASTA not found: {genome_fasta}")
    if not annotation_gtf.exists():
        raise FileNotFoundError(f"GTF annotation not found: {annotation_gtf}")

    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = logs_dir / "star_index.log"

    cmd = [
        "STAR",
        "--runThreadN",
        str(threads),
        "--runMode",
        "genomeGenerate",
        "--genomeDir",
        str(star_index_dir),
        "--genomeFastaFiles",
        str(genome_fasta),
        "--sjdbGTFfile",
        str(annotation_gtf),
        "--sjdbOverhang",
        str(sjdb_overhang),
    ]
    run_cmd(cmd, log_file)



def star_output_bam(sample_outdir: Path) -> Path:
    return sample_outdir / "Aligned.sortedByCoord.out.bam"


def run_star_one(
    sample: SampleRow,
    star_index_dir: Path,
    sample_outdir: Path,
    threads: int,
    force: bool,
    logs_dir: Path,
) -> None:
    sample_outdir.mkdir(parents=True, exist_ok=True)

    bam = star_output_bam(sample_outdir)
    if bam.exists() and not force:
        return

    prefix = str(sample_outdir / (sample.sample + "_"))
    log_file = logs_dir / f"star_{sample.sample}.log"

    cmd = [
        "STAR",
        "--runThreadN",
        str(threads),
        "--genomeDir",
        str(star_index_dir),
        "--outFileNamePrefix",
        prefix,
        "--outSAMtype",
        "BAM",
        "SortedByCoordinate",
    ]

    read_files = [sample.read1] + ([sample.read2] if sample.read2 is not None else [])
    if any(str(p).endswith(".gz") for p in read_files):
        cmd += ["--readFilesCommand", "zcat"]

    cmd += ["--readFilesIn", str(sample.read1)]
    if sample.read2 is not None:
        cmd += [str(sample.read2)]

    run_cmd(cmd, log_file)

    produced = Path(prefix + "Aligned.sortedByCoord.out.bam")
    if produced.exists() and produced != bam:
        produced.replace(bam)

    for name in ["Log.final.out", "Log.out", "Log.progress.out", "SJ.out.tab"]:
        p = Path(prefix + name)
        if p.exists():
            p.replace(sample_outdir / name)


def STAR_alignment() -> None:
    parser = argparse.ArgumentParser(
        description="Alignment stage: STAR per sample (samples.csv), coordinate-sorted BAM."
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
        help="samples.csv path relative to project root (default: config/samples.csv)",
    )

    # Canonical reference paths (override from full-pipeline.py if needed)
    parser.add_argument("--genome-fasta", type=Path, default=DEFAULT_GENOME_FASTA)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_ANNOT_GTF)
    parser.add_argument("--star-index", type=Path, default=DEFAULT_STAR_INDEX)

    parser.add_argument("--threads", type=int, default=8, help="Threads per STAR run")
    parser.add_argument("--jobs", type=int, default=1, help="Concurrent samples to align")

    parser.add_argument("--build-index", action="store_true", help="Build STAR index if missing")
    parser.add_argument("--force-index", action="store_true", help="Force rebuild STAR index")
    parser.add_argument("--sjdb-overhang", type=int, default=149, help="Read length - 1")

    parser.add_argument("--force", action="store_true", help="Re-run even if BAM exists")
    args = parser.parse_args()

    which_or_fail("STAR")

    project_root = normalize_project_root(args.project_root)

    samples_csv = (project_root / args.samples).resolve() if not args.samples.is_absolute() else args.samples
    samples = parse_samples_csv(samples_csv, project_root)

    genome_fasta = (project_root / args.genome_fasta).resolve() if not args.genome_fasta.is_absolute() else args.genome_fasta
    annotation_gtf = (project_root / args.gtf).resolve() if not args.gtf.is_absolute() else args.gtf
    star_index_dir = (project_root / args.star_index).resolve() if not args.star_index.is_absolute() else args.star_index

    align_root = project_root / "results" / "align"
    logs_dir = align_root / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    if args.build_index or args.force_index:
        ensure_star_index(
            star_index_dir=star_index_dir,
            genome_fasta=genome_fasta,
            annotation_gtf=annotation_gtf,
            threads=args.threads,
            force_index=args.force_index,
            logs_dir=logs_dir,
            sjdb_overhang=args.sjdb_overhang,
        )

    errors: List[str] = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {
            ex.submit(run_star_one, s, star_index_dir, align_root / s.sample, args.threads, args.force, logs_dir): s.sample
            for s in samples
        }
        for fut in as_completed(futs):
            sample_name = futs[fut]
            try:
                fut.result()
            except Exception as e:
                errors.append(f"{sample_name}: {e}")

    if errors:
        raise RuntimeError("Some STAR alignments failed:\n" + "\n".join("- " + x for x in errors))

    print("Alignment finished.")
    print(f"Project root: {project_root}")
    print(f"Align root:   {align_root}")
    print(f"Logs:         {logs_dir}")
    if samples:
        print("Example BAM:  ", align_root / samples[0].sample / "Aligned.sortedByCoord.out.bam")