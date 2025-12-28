# src/count.py
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import List
from utils import normalize_project_root, parse_samples_csv, run_cmd, which_or_fail


@dataclass(frozen=True)
class CountOutputs:
    counts_tsv: Path
    summary_tsv: Path
    log_file: Path


def _detect_paired_end(samples) -> bool:
    """
    Determine whether to run featureCounts in paired-end mode (-p).

    Rules:
    - If ALL samples have read2 -> paired-end
    - If ALL samples have no read2 -> single-end
    - If mixed -> raise (pipeline consistency)
    """
    has_r2 = [(s.read2 is not None) for s in samples]
    if all(has_r2):
        return True
    if not any(has_r2):
        return False
    raise ValueError(
        "Mixed single-end and paired-end samples detected in samples.csv. "
        "Use a consistent library layout for one pipeline run."
    )


def _collect_bams(project_root: Path, samples) -> List[Path]:
    """
    Collect STAR output BAMs based on your canonical layout:
      results/align/<sample>/Aligned.sortedByCoord.out.bam
    """
    bams: List[Path] = []
    missing: List[str] = []

    for s in samples:
        bam = project_root / "results" / "align" / s.sample / "Aligned.sortedByCoord.out.bam"
        if not bam.exists():
            missing.append(f"{s.sample}: {bam}")
        else:
            bams.append(bam)

    if missing:
        raise FileNotFoundError(
            "Missing BAM(s). Run STAR_alignment first.\n" + "\n".join(missing)
        )

    return bams


def featureCounts_run(
    *,
    project_root: str,
    samples_csv: Path,
    annotation_gtf: Path,
    threads: int = 8,
    stranded: int = 0,         # 0=unstranded, 1=stranded, 2=reverse
    force: bool = False,
    feature_type: str = "exon",  # -t
    gene_attr: str = "gene_id",  # -g
) -> CountOutputs:
    """
    Run featureCounts to produce gene-level raw counts for DESeq2.

    Inputs:
      - samples.csv
      - annotation GTF
      - STAR BAMs at results/align/<sample>/Aligned.sortedByCoord.out.bam

    Outputs:
      results/counts/
        - gene_counts.raw.tsv
        - gene_counts.raw.tsv.summary
        - featureCounts.log
    """
    if stranded not in (0, 1, 2):
        raise ValueError("stranded must be 0, 1, or 2.")

    # Normalization boundary
    root = normalize_project_root(project_root)

    # Resolve project-relative paths
    samples_csv_path = samples_csv if samples_csv.is_absolute() else (root / samples_csv)
    gtf_path = annotation_gtf if annotation_gtf.is_absolute() else (root / annotation_gtf)

    if not samples_csv_path.exists():
        raise FileNotFoundError(f"samples.csv not found: {samples_csv_path}")
    if not gtf_path.exists():
        raise FileNotFoundError(f"Annotation GTF not found: {gtf_path}")

    which_or_fail("featureCounts")

    # Your utils.parse_samples_csv expects (samples_csv, project_root)
    samples = parse_samples_csv(samples_csv_path, root)

    paired_end = _detect_paired_end(samples)
    bams = _collect_bams(root, samples)

    out_dir = root / "results" / "counts"
    out_dir.mkdir(parents=True, exist_ok=True)

    counts_tsv = out_dir / "gene_counts.raw.tsv"
    summary_tsv = out_dir / "gene_counts.raw.tsv.summary"
    log_file = out_dir / "featureCounts.log"

    if counts_tsv.exists() and (not force):
        return CountOutputs(counts_tsv=counts_tsv, summary_tsv=summary_tsv, log_file=log_file)

    cmd: List[str] = [
        "featureCounts",
        "-T", str(max(1, int(threads))),
        "-a", str(gtf_path),
        "-o", str(counts_tsv),
        "-t", feature_type,
        "-g", gene_attr,
        "-s", str(stranded),
    ]

    # Paired-end fragment counting (recommended when read2 exists)
    if paired_end:
        cmd += ["-p", "-B", "-C"]

    # Add BAM inputs
    cmd += [str(b) for b in bams]

    run_cmd(cmd, log_file=log_file)

    return CountOutputs(counts_tsv=counts_tsv, summary_tsv=summary_tsv, log_file=log_file)