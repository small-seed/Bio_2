# Reference-based RNA-seq Practice Pipeline (Drosophila)

This repository contains a **practice implementation of a reference-based RNA-seq pipeline**, inspired by and following the conceptual steps of the Galaxy ARTbio tutorial:

👉 https://artbio.github.io/startbio/reference_based_RNAseq/

The goal of this project is **educational and technical practice**, not biological discovery.

---

## Overview

- Pipeline type: Reference-based RNA-seq analysis
- Organism: *Drosophila melanogaster*
- Reference genome: Ensembl BDGP6 (release 95)
- Annotation: Ensembl GTF (release 95)
- Alignment tool: STAR
- Language: Python
- Platform: Windows (WSL / local execution), **not HPC Linux**

Compared to the original Galaxy tutorial, this implementation:
- Uses a **smaller / partial subset of the original datasets**
- Replaces Galaxy tools with **command-line tools orchestrated in Python**
- Runs on **Windows (via WSL or local shell)** instead of an HPC environment

---

## Dataset Description

The RNA-seq data were generated through **deep sequencing of mRNA** from *Drosophila melanogaster* **S2-DRSC cells** that were RNAi-depleted of mRNAs encoding RNA-binding proteins.

In the original tutorial, **7 datasets** are used to study the effect of **Pasilla gene inactivation** by RNAi knock-down:

### Untreated samples
- GSM461176  
- GSM461177  
- GSM461178  
- GSM461182  

### Treated samples (Pasilla RNAi knock-down)
- GSM461179  
- GSM461180  
- GSM461181  

In this repository, **only a subset of the corresponding SRA runs** is used for practice purposes.

---

## Pipeline Steps (High-level)

The pipeline follows the same logical structure as the Galaxy tutorial:

1. Raw data quality control (FastQC)
2. Reference genome indexing (STAR)
3. Read alignment to the reference genome
4. Gene-level read counting
5. (Optional) Downstream differential expression analysis

Each step is implemented as a **Python module** that calls standard bioinformatics tools via subprocess execution.

---

## Purpose and Scope

- This repository is intended for **learning, experimentation, and pipeline development**
- Results should **not** be interpreted as biological conclusions
- The focus is on:
  - Reproducible pipeline structure
  - Correct handling of references and metadata
  - Translating Galaxy workflows into scriptable pipelines

---

## Reference

Original tutorial and data description:

> ARTbio Galaxy Training – Reference-based RNA-seq analysis  
> https://artbio.github.io/startbio/reference_based_RNAseq/

---

## Notes

- Reference genome and annotation versions must match (BDGP6, Ensembl release 95)
- Mixing reference versions may lead to incorrect results
- This project intentionally prioritizes clarity and reproducibility over performance