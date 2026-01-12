from __future__ import annotations

from pathlib import Path

from utils import normalize_project_root, run_cmd

TEMPLATE_R = Path(__file__).with_name("deseq2_template.R")
DEFAULT_COUNTS = Path("results/counts/gene_counts.raw.tsv")
DEFAULT_OUTDIR = Path("results/deseq2")


def DESeq2_run(
    *,
    project_root: str,
    samples_csv: Path,
    condition_col: str = "condition",
    counts_tsv: Path = DEFAULT_COUNTS,
    out_dir: Path = DEFAULT_OUTDIR,
    alpha: float = 0.05,
    lfc_threshold: float = 0.0,
    contrast: tuple[str, str] | None = None,  # (ref, test)
    force: bool = False,
) -> None:
    root = Path(normalize_project_root(project_root))

    counts = (root / counts_tsv).resolve() if not counts_tsv.is_absolute() else counts_tsv
    samples = (root / samples_csv).resolve() if not samples_csv.is_absolute() else samples_csv
    outdir = (root / out_dir).resolve() if not out_dir.is_absolute() else out_dir

    if not counts.exists():
        raise FileNotFoundError(f"Counts TSV not found: {counts}")
    if not samples.exists():
        raise FileNotFoundError(f"samples.csv not found: {samples}")
    if not TEMPLATE_R.exists():
        raise FileNotFoundError(f"Missing R template: {TEMPLATE_R}")

    outdir.mkdir(parents=True, exist_ok=True)

    results = outdir / "deseq_results.tsv"
    if results.exists() and not force:
        return

    template = TEMPLATE_R.read_text(encoding="utf-8")

    contrast_enabled = "1" if contrast is not None else "0"
    contrast_ref = contrast[0] if contrast else ""
    contrast_test = contrast[1] if contrast else ""

    filled = (
        template
        .replace("{{CONDITION_COL}}", condition_col)
        .replace("{{ALPHA}}", str(alpha))
        .replace("{{LFC_THRESHOLD}}", str(lfc_threshold))
        .replace("{{CONTRAST_ENABLED}}", contrast_enabled)
        .replace("{{CONTRAST_REF}}", contrast_ref)
        .replace("{{CONTRAST_TEST}}", contrast_test)
    )

    # Freeze the exact analysis script used
    r_script = outdir / "run_deseq2.R"
    r_script.write_text(filled, encoding="utf-8")

    cmd = [
        "Rscript",
        str(r_script),
        "--counts", str(counts),
        "--samples", str(samples),
        "--outdir", str(outdir),
    ]

    run_cmd(cmd, log_file=outdir / "deseq2.log")