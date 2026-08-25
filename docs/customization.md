# NovaScope Customization for P4 Flowcell

## 1. P4 Flowcell Structure

The P4 flowcell consists of four spatial regions:

| Surface | Lane | Region |
|---|---|---|
| Top | Lane 1 | 1000_lane1 |
| Top | Lane 2 | 1000_lane2 |
| Bottom | Lane 1 | 2000_lane1 |
| Bottom | Lane 2 | 2000_lane2 |

Each region contains approximately a 6 × 16–17 tile arrangement.

## 2. Spatial Barcode Structure

The matrix FASTQ contains 37-bp reads.

The N32 spatial barcode is located at positions 5–36.

Python slicing:

```python
barcode = seq[4:36]

