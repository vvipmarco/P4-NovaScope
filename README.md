# P4-NovaScope

Customized NovaScope workflow for P4 flowcell spatial transcriptomics data.

## Overview

This repository contains custom scripts and configuration files used to adapt the NovaScope pipeline to P4 flowcell spatial transcriptomics data.

The main objective is to integrate:

1. spatial barcode information from the matrix FASTQ, and
2. transcriptome FASTQ data generated from tissue placed on the spatially barcoded flowcell,

to obtain spatially resolved transcriptomic information.

## P4 Flowcell Structure

The P4 flowcell is divided into four spatial regions:

- 1000_lane1
- 1000_lane2
- 2000_lane1
- 2000_lane2

The top surface corresponds to the 1000-series regions, and the bottom surface corresponds to the 2000-series regions.

## Barcode Processing

The matrix FASTQ reads are 37 bp long.

The N32 spatial barcode is extracted from positions 5–36:

```python
barcode = seq[4:36]
