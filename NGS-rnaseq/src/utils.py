from __future__ import annotations

import csv
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional


# -----------------------------
# Data model
# -----------------------------
@dataclass(frozen=True)
class SampleRow:
    sample: str
    read1: Path
    read2: Optional[Path] = None
    condition: Optional[str] = None


# -----------------------------
# Environment / path helpers
# -----------------------------
def is_wsl() -> bool:
    return "microsoft" in os.uname().release.lower() if hasattr(os, "uname") else False


def windows_to_wsl_path(p: str) -> str:
    """
    Convert 'D:\\Code\\Bio_2\\NGS-rnaseq' -> '/mnt/d/Code/Bio_2/NGS-rnaseq'
    Only used when running inside WSL.
    """
    p = p.strip()
    if p.startswith("/"):
        return p
    if len(p) >= 2 and p[1] == ":":
        drive = p[0].lower()
        rest = p[2:].lstrip("\\/").replace("\\", "/")
        return f"/mnt/{drive}/{rest}"
    return p.replace("\\", "/")


def normalize_project_root(raw_root: str) -> Path:
    if is_wsl() and (":" in raw_root[:3] or raw_root.startswith("\\")):
        raw_root = windows_to_wsl_path(raw_root)
    return Path(raw_root).resolve()


# -----------------------------
# Command helpers
# -----------------------------
def which_or_fail(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise FileNotFoundError(
            f"Required tool '{tool}' not found in PATH.\n"
            f"Install with: conda install -c bioconda {tool}"
        )
    return path


def run_cmd(cmd: List[str], log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("w", encoding="utf-8") as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nSee log: {log_file}"
        )


# -----------------------------
# samples.csv parsing
# -----------------------------
def parse_samples_csv(samples_csv: Path, project_root: Path) -> List[SampleRow]:
    if not samples_csv.exists():
        raise FileNotFoundError(f"samples.csv not found: {samples_csv}")

    rows: List[SampleRow] = []
    with samples_csv.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        required = {"sample", "read1"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(
                f"samples.csv must contain columns: {sorted(required)} "
                f"(found: {reader.fieldnames})"
            )

        for i, r in enumerate(reader, start=2):
            sample = (r.get("sample") or "").strip()
            read1 = (r.get("read1") or "").strip()
            read2 = (r.get("read2") or "").strip()
            condition = (r.get("condition") or "").strip() or None

            if not sample or not read1:
                raise ValueError(f"Invalid row at line {i}: sample/read1 required")

            r1 = (project_root / read1).resolve() if not Path(read1).is_absolute() else Path(read1)
            r2 = None
            if read2:
                r2 = (project_root / read2).resolve() if not Path(read2).is_absolute() else Path(read2)

            if not r1.exists():
                raise FileNotFoundError(f"read1 not found for sample '{sample}': {r1}")
            if r2 is not None and not r2.exists():
                raise FileNotFoundError(f"read2 not found for sample '{sample}': {r2}")

            rows.append(
                SampleRow(
                    sample=sample,
                    read1=r1,
                    read2=r2,
                    condition=condition,
                )
            )

    # enforce unique sample names
    seen: Dict[str, int] = {}
    for s in rows:
        seen[s.sample] = seen.get(s.sample, 0) + 1
    dupes = [k for k, v in seen.items() if v > 1]
    if dupes:
        raise ValueError(f"Duplicate sample names in samples.csv: {dupes}")

    return rows