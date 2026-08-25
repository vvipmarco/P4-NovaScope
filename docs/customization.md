# NovaScope Customization for P4 Flowcell

This document summarizes the main modifications made to the original NovaScope workflow for the P4 flowcell spatial transcriptomics analysis.

## 1. P4 Flowcell Structure

For this analysis, the P4 flowcell was divided into four spatial regions:

| Surface | Lane   | Region       |
| ------- | ------ | ------------ |
| Top     | Lane 1 | `1000_lane1` |
| Top     | Lane 2 | `1000_lane2` |
| Bottom  | Lane 1 | `2000_lane1` |
| Bottom  | Lane 2 | `2000_lane2` |

Each region contains approximately a 6 × 16–17 tile arrangement.

The transcriptome samples were expected to correspond to the following surface groups:

* Top transcriptome samples → 1000-series regions
* Bottom transcriptome samples → 2000-series regions

The surface/lane-specific tile information was obtained from the corresponding `xynum` files.

Examples:

```text
1000_lane1_xynum.txt
1000_lane2_xynum.txt
2000_lane1_xynum.txt
2000_lane2_xynum.txt
```

## 2. Spatial Barcode Structure

In this dataset, the matrix FASTQ reads are 37 bp long.

The N32 spatial barcode is located at positions 5–36 using 1-based indexing.

Python slicing:

```python
barcode = seq[4:36]
```

Because Python uses 0-based indexing and an exclusive end position, `seq[4:36]` corresponds to positions 5–36 of the FASTQ read.

The extracted 32-bp barcode is reverse-complemented before transcriptome matching.

For example, the reverse-complement operation follows:

```text
A ↔ T
C ↔ G
```

The custom Python implementation is:

```python
def reverse_complement(seq: str) -> str:
    complement = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )
    return seq.translate(complement)[::-1]
```

## 3. Surface/Lane-specific Whitelist Generation

The original workflow was supplemented with a custom script for generating surface/lane-specific whitelist files directly from the matrix FASTQ.

The main script is:

```text
scripts/make_whitelist_all_lanes.py
```

The script processes the matrix FASTQ once and generates four whitelist files:

```text
whitelist_1000_lane1_RC.txt
whitelist_1000_lane2_RC.txt
whitelist_2000_lane1_RC.txt
whitelist_2000_lane2_RC.txt
```

### Processing Steps

The script performs the following operations:

1. Loads the target tile IDs from the four `xynum` files.
2. Reads the matrix FASTQ sequentially.
3. Extracts the tile ID from each FASTQ header.
4. Determines which P4 surface/lane the read belongs to.
5. Extracts the N32 barcode using `seq[4:36]`.
6. Reverse-complements the barcode.
7. Writes the barcode to the corresponding whitelist file.
8. Generates a summary TSV containing tile and read counts.

The FASTQ header is parsed using:

```python
tile_id = header.split(":")[4]
```

An example FASTQ header has the following structure:

```text
@VH01546:22:2225LYHNX:1:1102:18944:1000
```

The relevant tile ID in this example is:

```text
1102
```

The custom script also checks whether a tile appears in more than one surface/lane definition. If the tile lists overlap, the script raises an error.

## 4. Reverse-complement Helper Script

A second helper script was used during barcode debugging:

```text
scripts/rc_barcode.py
```

The script performs reverse-complement conversion for sequences provided through standard input.

```python
import sys

comp = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "N": "N",
}

for line in sys.stdin:
    seq = line.strip()

    rc_seq = "".join(
        comp.get(base, base)
        for base in reversed(seq)
    )

    print(rc_seq)
```

This helper script was mainly used for direct barcode comparison and debugging.

The main whitelist-generation script already contains its own reverse-complement function.

## 5. NovaScope Configuration Changes

The original NovaScope example configuration was modified to accommodate the P4 flowcell data.

Representative changes are summarized below.

| Parameter                        | Original / Example       | P4 Configuration                            |
| -------------------------------- | ------------------------ | ------------------------------------------- |
| `fastq2sbcd.format`              | `DraI31`                 | `DraI32`                                    |
| `sbcd_layout.tiles`              | Example tile pairs       | P4-specific tile pairs                      |
| `sbcd_layout.colshift`           | `0.1715`                 | `-0.000010`                                 |
| `sbcd2chip.gap_col`              | `0.0048`                 | `0.039637`                                  |
| `sbcd2chip.gap_row`              | `0.0517`                 | `0.051700`                                  |
| `smatch.match_len`               | earlier/default setting  | `27`                                        |
| `whitelist`                      | not externally specified | surface/lane-specific whitelist             |
| `gene_filter.min_ct_per_feature` | `50`                     | `1`                                         |
| `gene_filter.keep_gene_type`     | `protein_coding\|lncRNA` | `protein_coding\|lncRNA\|pseudogene\|miRNA` |

The exact configuration depends on the transcriptome sample and the surface/lane being analyzed.

## 6. Barcode Format

The example NovaScope configuration used:

```yaml
fastq2sbcd:
  format: DraI31
```

The P4 matrix contains an N32 spatial barcode, so the format was changed to:

```yaml
fastq2sbcd:
  format: DraI32
```

## 7. P4-specific Tile Configuration

The example NovaScope tile configuration did not correspond to the P4 flowcell layout.

The tile pairs were therefore modified using P4 flowcell information.

An example configuration was:

```yaml
sbcd_layout:
  tiles:
    - "1101,1201"
    - "1401,1501"
```

These values were selected after examining the available flowcell metadata and the actual P4 tile arrangement.

## 8. Spatial Coordinate Parameters

The spatial coordinate parameters were modified to better represent the P4 flowcell geometry.

### 8.1 `gap_col`

The original example value was:

```yaml
gap_col: 0.0048
```

The P4-specific value used in the analysis was:

```yaml
gap_col: 0.039637
```

The value was derived using information from `AutofocusReport.csv`.

The normalized horizontal coordinate was defined as:

```text
x' = (x_nm - centroid_1) / (centroid_2 - centroid_1)
```

The mean distance between horizontally adjacent tiles was calculated as:

```text
Δx' = mean(|x'_j - x'_i|)
```

The horizontal barcode spacing was then calculated as:

```text
gap_col = Δx' / 6
```

### 8.2 `gap_row`

The normalized vertical coordinate was defined as:

```text
y' = y_nm / Spot Separation
```

The mean distance between vertically adjacent tiles was calculated as:

```text
Δy' = mean(|y'_j - y'_i|)
```

The vertical barcode spacing was then calculated as:

```text
gap_row = Δy' / 16
```

The resulting value was:

```yaml
gap_row: 0.051700
```

This was effectively consistent with the original value of approximately `0.0517`.

### 8.3 `colshift`

The original value was:

```yaml
colshift: 0.1715
```

For the P4 data, the value was changed to:

```yaml
colshift: -0.000010
```

The value was estimated using barcode coordinates from `sbcds.sorted.tsv.gz`.

A subset of the data was used to calculate the residual of the x-coordinate relative to the `gap_col = 0.039637` lattice.

If the median residual is denoted by `m`, the shift was defined as:

```text
colshift = -m
```

The calculated median residual was approximately:

```text
m = 0.000010
```

resulting in:

```text
colshift = -0.000010
```

## 9. Spatial Barcode Matching Parameters

The spatial barcode matching step uses the `smatch` configuration.

The P4 configuration included:

```yaml
smatch:
  skip_sbcd: 0
  match_len: 27
```

`skip_sbcd` was set to `0` because the spatial barcode information was being regenerated and customized for the P4 dataset.

`match_len` was set to `27` during the debugging process after the earlier longer setting caused execution problems.

## 10. Transcriptome Input

The transcriptome data are paired-end FASTQ files.

Example:

```yaml
seq2nd:
  - id: "bottom_32_S4"
    fastq_R1: "/data/5th_data/transcriptome-uncut/bottom_32_S4_R1_001.fastq.gz"
    fastq_R2: "/data/5th_data/transcriptome-uncut/bottom_32_S4_R2_001.fastq.gz"
```

The roles of the two reads differ across the NovaScope workflow.

### Spatial Barcode Matching

The `a03 (smatch)` stage uses transcriptome R1 for spatial barcode matching.

### Alignment

From the `a04 (align)` stage, transcriptome R1 and R2 are used together as paired-end reads.

Therefore, the workflow should not be interpreted as using only R1 for the entire analysis.

Rather:

```text
Spatial barcode matching → primarily R1
Genome alignment → paired R1/R2
```

## 11. Alignment Configuration

An example alignment configuration used for the P4 analysis was:

```yaml
align:
  min_match_len: 30
  min_match_frac: 0.66
  len_sbcd: 32
  len_umi: 9
  len_r2: 101
  exist_action: overwrite
```

During the analysis, a STAR genome-index compatibility issue was also observed.

The existing GRCh38 index had been generated using STAR version `2.7.1a`, whereas the NovaScope Docker environment used STAR version `2.7.11b`.

This produced an incompatibility error and indicated that the genome index would need to be regenerated using a compatible STAR version before the alignment stage could be completed.

## 12. Gene Filtering

The original example configuration used:

```yaml
gene_filter:
  keep_gene_type: "protein_coding|lncRNA"
  rm_gene_regex: "^Gm\\d+|^mt-|^MT-"
  min_ct_per_feature: 50
```

For the P4 analysis, the filtering criteria were relaxed to retain a broader range of genes:

```yaml
gene_filter:
  keep_gene_type: "protein_coding|lncRNA|pseudogene|miRNA"
  rm_gene_regex: "^Gm\\d+|^mt-|^MT-"
  min_ct_per_feature: 1
```

The purpose of this modification was to initially retain low-expression features and additional noncoding gene types, followed by later filtering if necessary.

## 13. NovaScope Processing Stages

The main NovaScope workflow used in this project can be summarized as:

```text
a01 → a02 → a03 → a04 → a05 and downstream
```

### a01: `fastq2sbcd`

The matrix FASTQ is processed according to the `DraI32` barcode format.

Spatial barcode information is extracted and converted into NovaScope intermediate files.

### a02: `sbcd2chip`

Spatial barcodes are mapped to physical flowcell coordinates.

The main P4-specific parameters used in this step include:

```text
tiles
gap_col
gap_row
colshift
dup_maxnum
dup_maxdist
```

The resulting spatial barcode-coordinate file includes:

```text
1_1.sbcds.sorted.tsv.gz
```

### a03: `smatch`

Transcriptome R1 spatial barcodes are matched to the matrix barcode-coordinate map.

The summary output provides information such as:

```text
Total reads
Miss
Match
Unique
Dup(Exact)
```

### a04: `align`

Matched transcriptome reads are aligned to the reference genome using STAR.

Paired-end R1/R2 data are used at this stage.

### a05 and Downstream Analysis

After alignment, downstream processing includes gene filtering, spatial aggregation, segmentation, and visualization.

## 14. Visualization Configuration

Visualization parameters were also modified during the P4 analysis.

An example configuration included:

```yaml
visualization:
  drawxy:
    coord_per_pixel: 5000
    intensity_per_obs: 50
    icol_x: 3
    icol_y: 4

  drawsge:
    action: True
    coord_per_pixel: 5000
    auto_adjust: false
    adjust_quantile: 0.99
```

Earlier analyses used `coord_per_pixel: 1000`, while later visualization settings used `5000`.

These parameters affect rendering and visualization rather than the underlying barcode-coordinate definition.

## 15. Custom Docker Image for Visualization

During visualization, the original NovaScope Docker image produced an ImageMagick security-policy error when converting BMP files to PNG.

A custom Docker image was therefore created.

The corresponding file is:

```text
docker/Dockerfile
```

The Dockerfile contains:

```dockerfile
FROM hyunminkang/novascope

RUN sed -i \
    's/rights="none"/rights="read|write"/g' \
    /etc/ImageMagick-6/policy.xml
```

The custom image can be built using:

```bash
docker build -t novascope-img-fixed .
```

The resulting image was used for NovaScope visualization steps requiring ImageMagick.

## 16. Barcode Debugging

When spatial barcode matching produced unexpectedly high mismatch rates, the matrix and transcriptome barcodes were compared directly.

An example debugging procedure was:

### Matrix Barcode Extraction

```bash
zcat matrix.fastq.gz \
  | awk 'NR%4==2 {print substr($0,5,32)}' \
  | python3 rc_barcode.py \
  | sort -u \
  > matrix_barcodes.txt
```

### Transcriptome Barcode Extraction

```bash
zcat transcriptome_R1.fastq.gz \
  | awk 'NR%4==2 {print substr($0,2,32)}' \
  | sort -u \
  > transcriptome_barcodes.txt
```

### Barcode Intersection

```bash
comm -12 matrix_barcodes.txt transcriptome_barcodes.txt | wc -l
```

These commands were used for debugging and were not a replacement for the complete NovaScope matching workflow.

## 17. FASTQ Header Interpretation

Matrix and transcriptome FASTQ files may contain different sequencing run identifiers in their headers.

However, the FASTQ run identifier itself is not used to establish the matrix-transcriptome pairing.

The relevant relationship is the spatial barcode sequence:

```text
Matrix:
coordinate ↔ spatial barcode

Transcriptome:
spatial barcode ↔ transcript
```

The objective of the analysis is therefore:

```text
coordinate ↔ spatial barcode ↔ transcript
```

which ultimately provides:

```text
coordinate ↔ transcript
```

## 18. Data and Generated Files

Raw sequencing data are not included in this GitHub repository.

The following files are intentionally excluded:

* Matrix FASTQ files
* Transcriptome FASTQ files
* Full generated whitelist files
* STAR reference genome indexes
* BAM/SAM alignment files
* Large NovaScope intermediate files
* Full NovaScope output directories

These files are stored separately because of their large size and/or project confidentiality.

The GitHub repository is intended to contain:

* Custom scripts
* Configuration files
* Docker configuration
* Documentation
* Small reproducibility-related files when appropriate

## 19. Current Validation Status

Initial spatial barcode-matching analyses produced measurable match rates.

However, later validation showed that the same barcode-coordinate map had been reused across multiple surface/lane configurations.

Therefore, the initial matching results could not be interpreted as independent comparisons between:

```text
1000_lane1
1000_lane2
2000_lane1
2000_lane2
```

The surface/lane-specific layouts were subsequently reconstructed so that the four P4 regions could be handled separately.

Further validation of the reconstructed layouts and downstream transcriptome mapping remains to be completed.

## 20. Summary of P4-specific Modifications

The main customizations implemented for this project are:

1. Extraction of the N32 barcode from positions 5–36 of the matrix FASTQ.
2. Reverse-complement conversion of the spatial barcode.
3. Generation of four surface/lane-specific whitelist files.
4. P4-specific tile and layout construction.
5. Modification of `gap_col`, `gap_row`, and `colshift`.
6. Use of `DraI32` instead of the example `DraI31` barcode format.
7. Modification of the spatial barcode matching configuration.
8. Use of paired-end R1/R2 transcriptome reads during alignment.
9. Relaxation of downstream gene-filtering criteria.
10. Creation of a custom Docker image for visualization.
11. Direct barcode-comparison utilities for debugging.
12. Separate handling and validation of the four P4 surface/lane regions.
