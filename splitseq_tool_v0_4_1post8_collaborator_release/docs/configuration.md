# Configuration

YAML configuration is optional. Every command can be run using command-line flags only.

## Resolution order

Values are resolved in this order:

1. explicit command-line option;
2. sample-specific YAML;
3. command-specific/project YAML;
4. `defaults:` YAML;
5. built-in default.

Boolean options support both positive and negative CLI forms, so a YAML `true` value can be overridden explicitly, for example:

```bash
--no-use-flanks
--no-round1-conversion
--no-keep-all-intermediates
```

## Example

```yaml
barcodes:
  mrna_round1: /path/to/Round1_RNA_barcodes.txt
  capture_round1: /path/to/Round1_capture_barcodes.txt
  round2: /path/to/Round2_barcodes.txt
  round3: /path/to/Round3_barcodes.txt

references:
  star_index: /path/to/STAR_index
  saf_file: /path/to/genes.saf

defaults:
  barcode_mismatches: 1
  capture_mismatches: 1
  cell_filter: keep-list
  round1_conversion: true

capture:
  single:
    target_sequence: ACGTACGTACGT
  multi:
    capture_mode: greenlist
    green_list: /path/to/approved_capture_barcodes.txt
    flank_left: AACCGGTTAACC
    flank_right: TTGGCCAATTGG
    target_length: 12

samples:
  sample_01:
    mrna_r1: /path/to/sample_01_R1.fastq.gz
    mrna_r2: /path/to/sample_01_R2.fastq.gz
    capture_merged_fastq: /path/to/capture/MergedCells_passing.fastq
    mrna_cell_keep_list: /path/to/rna/cellIDs_keep.txt
```

All nucleotide strings above are synthetic placeholders.

## Reproducibility files

When an output directory is available, each command writes:

```text
resolved_parameters.tsv
resolved_config.yaml
```

These files record the effective settings after YAML and CLI values are merged.
