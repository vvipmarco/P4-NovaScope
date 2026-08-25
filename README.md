# P4-NovaScope

Customized NovaScope workflow for P4 flowcell spatial transcriptomics data.

## Overview

This repository contains custom scripts and configuration files used to adapt the NovaScope pipeline to P4 flowcell spatial transcriptomics data.

The main objective is to integrate:

1. spatial barcode information from the matrix FASTQ, and
2. transcriptome FASTQ data generated from tissue placed on the spatially barcoded flowcell,

to obtain spatially resolved transcriptomic information.

## P4 Flowcell Structure

For this analysis, the P4 flowcell was divided into four spatial regions:

* `1000_lane1`
* `1000_lane2`
* `2000_lane1`
* `2000_lane2`

The top surface corresponds to the 1000-series regions, and the bottom surface corresponds to the 2000-series regions.

Each region consists of an approximately 6 × 16–17 tile arrangement.

## Barcode Processing

In this dataset, the matrix FASTQ reads are 37 bp long.

The N32 spatial barcode is extracted from positions 5–36 (1-based indexing):

```python
barcode = seq[4:36]
```

Because Python uses 0-based indexing and an exclusive end position, `seq[4:36]` corresponds to positions 5–36 of the FASTQ read.

The extracted 32-bp spatial barcode is reverse-complemented before matching with the transcriptome barcode.

## Repository Structure

```text
P4-NovaScope/
├── README.md
├── .gitignore
├── scripts/
│   ├── make_whitelist_all_lanes.py
│   └── rc_barcode.py
├── configs/
│   ├── config_top1.yaml
│   └── config_bottom2.yaml
├── docker/
│   └── Dockerfile
└── docs/
    └── customization.md
```

* `scripts/`: Custom Python scripts for barcode processing and lane-specific whitelist generation.
* `configs/`: NovaScope configuration files modified for the P4 flowcell data.
* `docker/`: Docker customization used to resolve visualization-related issues.
* `docs/`: Detailed documentation of modifications made to the original NovaScope workflow.

## Main Customizations

The original NovaScope workflow was modified to accommodate the P4 flowcell dataset.

Major modifications include:

* Use of the `DraI32` barcode format for the N32 spatial barcode
* Extraction of the N32 spatial barcode from positions 5–36 of the matrix FASTQ
* Reverse-complement conversion of the extracted spatial barcode
* Generation of surface- and lane-specific whitelist files
* P4-specific tile selection and spatial layout configuration
* Modification of spatial coordinate parameters including `gap_col`, `gap_row`, and `colshift`
* Modification of spatial barcode matching parameters
* Use of paired-end transcriptome R1/R2 data for the alignment stage
* Custom Docker configuration for NovaScope visualization

Detailed parameter changes and their rationale are documented in `docs/customization.md`.

## Lane-specific Whitelist Generation

The custom script `scripts/make_whitelist_all_lanes.py` generates four surface/lane-specific whitelist files:

* `whitelist_1000_lane1_RC.txt`
* `whitelist_1000_lane2_RC.txt`
* `whitelist_2000_lane1_RC.txt`
* `whitelist_2000_lane2_RC.txt`

The script performs the following steps:

1. Loads the target tile IDs for each surface/lane.
2. Reads the matrix FASTQ.
3. Extracts the tile ID from each FASTQ header.
4. Extracts the N32 barcode using `seq[4:36]`.
5. Reverse-complements the barcode.
6. Writes the barcode to the corresponding surface/lane-specific whitelist.
7. Generates a summary TSV containing tile and read counts.

The helper script `scripts/rc_barcode.py` was additionally used for direct barcode comparison and debugging.

## Workflow

The customized analysis workflow consists of the following major steps:

1. Identify P4 flowcell tiles for each surface and lane.
2. Extract N32 spatial barcodes from the matrix FASTQ.
3. Reverse-complement the extracted spatial barcodes.
4. Generate surface- and lane-specific whitelist files.
5. Construct P4-specific spatial layouts.
6. Run the NovaScope pipeline.
7. Match transcriptome spatial barcodes to the matrix barcode-coordinate map.
8. Align matched transcriptome reads to the reference genome.
9. Perform downstream spatial transcriptomics analysis and visualization.

Within the NovaScope pipeline, the major stages are:

* `a01 (fastq2sbcd)`: Extract spatial barcode information from the matrix FASTQ.
* `a02 (sbcd2chip)`: Map spatial barcodes to flowcell coordinates.
* `a03 (smatch)`: Match transcriptome spatial barcodes to the barcode-coordinate map.
* `a04 (align)`: Align matched transcriptome reads to the reference genome.
* `a05 and downstream`: Perform downstream gene-level and spatial analyses.

The spatial barcode matching step uses transcriptome R1, whereas paired-end R1/R2 reads are used in the alignment stage.

## Usage

### 1. Generate Surface/Lane-specific Whitelists

The four whitelist files can be generated using the custom Python script.

Example with Docker:

```bash
docker run -it --rm \
  -v "/path/to/P4_test_data:/data" \
  --entrypoint /bin/bash \
  hyunminkang/novascope \
  -c "python3 /data/make_whitelist_all_lanes.py"
```

Paths can also be specified explicitly through command-line arguments.

Example:

```bash
python3 make_whitelist_all_lanes.py \
  --matrix-fastq /path/to/matrix.fastq.gz \
  --tile-1000-lane1 /path/to/1000_lane1_xynum.txt \
  --tile-1000-lane2 /path/to/1000_lane2_xynum.txt \
  --tile-2000-lane1 /path/to/2000_lane1_xynum.txt \
  --tile-2000-lane2 /path/to/2000_lane2_xynum.txt
```

### 2. Run NovaScope

A P4-specific configuration file is supplied to the NovaScope pipeline.

Example:

```bash
docker run -it --rm \
  -v "/path/to/P4_test_data:/data" \
  hyunminkang/novascope \
  -s /app/novascope/NovaScope.smk \
  --rerun-incomplete \
  -d /data/output \
  --configfile /data/config.yaml \
  -p \
  --cores 10
```

The actual configuration files used in this project are stored in the `configs/` directory.

## Major Configuration Changes

Representative changes from the original NovaScope example configuration include:

| Parameter      | Original / Example       | P4 Configuration                |
| -------------- | ------------------------ | ------------------------------- |
| Barcode format | `DraI31`                 | `DraI32`                        |
| `gap_col`      | `0.0048`                 | `0.039637`                      |
| `gap_row`      | `0.0517`                 | `0.051700`                      |
| `colshift`     | `0.1715`                 | `-0.000010`                     |
| `match_len`    | default / 32             | `27`                            |
| Whitelist      | Not externally specified | Surface/lane-specific whitelist |

The coordinate-related parameters were modified based on the P4 flowcell metadata and observed barcode-coordinate distributions.

Additional configuration changes for transcriptome samples, visualization, and downstream gene filtering are provided in the configuration files and `docs/customization.md`.

## Visualization

A custom Docker image was used because the original NovaScope Docker environment produced an ImageMagick policy error during image conversion.

The corresponding Dockerfile is provided in:

```text
docker/Dockerfile
```

The customized image enables conversion of intermediate BMP visualization files to PNG format.

## Data

Raw sequencing data are not included in this repository.

The following files are intentionally excluded because of data size and/or project confidentiality:

* Matrix FASTQ files
* Transcriptome FASTQ files
* Reference genome indexes
* Generated full whitelist files
* Large NovaScope intermediate files
* Alignment files
* NovaScope output directories

These files are stored separately and connected to the workflow through local or mounted paths.

## Current Status

Initial spatial barcode-matching analyses produced measurable match rates.

Subsequent validation showed that the same barcode-coordinate map had been reused across multiple surface/lane configurations. Therefore, the initial matching results could not be interpreted as independent surface- and lane-specific comparisons.

The spatial layouts were subsequently reconstructed so that the four P4 regions could be handled separately:

* `1000_lane1`
* `1000_lane2`
* `2000_lane1`
* `2000_lane2`

Further validation of the reconstructed surface/lane-specific layouts and downstream transcriptome mapping remains to be completed.

## Notes

This repository contains code and configuration files developed or modified specifically for the P4 flowcell analysis.

Raw sequencing data and large generated outputs are managed separately and are not tracked by Git.
