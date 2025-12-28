# Config Setup (Windows + WSL + VS Code + Miniconda)

This document describes how to set up a reproducible environment for this reference-based RNA-seq pipeline on Windows, using WSL for Linux bioinformatics tools and Miniconda for dependency management.  
The goal is cross-user reproducibility.

---

## Why this setup exists

Most RNA-seq tools (FastQC, MultiQC, STAR, featureCounts, DESeq2) are developed and tested primarily on Linux.WSL (Windows Subsystem for Linux) provides a Linux environment on Windows, enabling consistent installation and execution. When running tools inside WSL, Linux-style paths are required.

Example:
    Windows: D:\Code\Bio_2\NGS-rnaseq
    WSL:     /mnt/d/Code/Bio_2/NGS-rnaseq

This pipeline automatically converts Windows paths to WSL paths when necessary.

---

## Prerequisites

### 1) Install WSL

Open PowerShell as Administrator and run: 
    wsl --install
Reboot if prompted.
Verify installation:
    wsl --status

---

### 2) Install VS Code + WSL extension
Install:
- Visual Studio Code
- VS Code extension: WSL (by Microsoft)
This allows opening and working with the repository directly inside the Linux (WSL) environment.

---

## Recommended repository location
Store the repository on a Windows drive but access it through WSL to avoid permission and performance issues.
Example:
    Windows: D:\Code\Bio_2\NGS-rnaseq
    WSL:     /mnt/d/Code/Bio_2/NGS-rnaseq

---

## Installing Miniconda (inside WSL)
The conda environment must be created inside WSL, not in native Windows Python.
In a WSL terminal:
    cd ~
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
Restart the shell or run:
    source ~/.bashrc
Confirm installation:
    conda --version

---

## Create the reproducible environment

### Option A: Use environment.yml (recommended)
From the project root in WSL:
    conda env create -f environment.yml
    conda activate rnaseq

### Option B: Create manually
    conda create -n rnaseq -c conda-forge -c bioconda python=3.11 \
        fastqc multiqc star subread bioconductor-deseq2 r-base
    conda activate rnaseq

---

## Verify installed tools
    python --version
    fastqc --version
    multiqc --version
    STAR --version
    featureCounts -v
    R --version

---

## Enable Miniconda in VS Code (WSL)
1. Open VS Code
2. Click the green >< icon (bottom-left) → Connect to WSL
3. Open the repository folder from within WSL (/mnt/d/...)
4. Open a new VS Code terminal (WSL terminal)
5. Activate the environment:
    conda activate rnaseq

Set Python interpreter in VS Code:
- Press Ctrl + Shift + P
- Select “Python: Select Interpreter”
- Choose the interpreter from the rnaseq conda environment (WSL)

---

## Tool summary

Quality control:
- FastQC: per-sample read quality reports
- MultiQC: aggregated QC reports

Alignment:
- STAR: spliced RNA-seq aligner producing BAM files

Counting:
- featureCounts (Subread): gene-level raw counts from BAM + GTF

Differential expression:
- DESeq2 (R / Bioconductor): statistical inference using raw counts and experimental design

---

## Cross-user reproducibility rules

- Install all dependencies via conda
- Commit and maintain environment.yml
- Avoid unpinned pip installs outside the environment
- Use project-relative paths in samples.csv
- Keep genome FASTA and GTF versions explicitly matched

---

## Notes on Windows vs WSL execution

- External tools (FastQC, MultiQC, STAR, featureCounts) must run inside WSL
- PROJECT_ROOT may be a Windows path and will be converted automatically
- Running the pipeline in native Windows Python is not recommended

---

## Quick smoke test

From the project root in WSL:

    conda activate rnaseq
    python -c "import sys; print(sys.version)"
    fastqc --version
    multiqc --version
    STAR --version
    featureCounts -v
    R -q -e "library('DESeq2'); sessionInfo()"

If all commands succeed, the environment is ready.