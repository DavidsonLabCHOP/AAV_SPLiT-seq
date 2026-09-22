# Sequence-quality QC

The toolkit reports quality directly from FASTQ quality strings using Phred scores.

Default settings:

```text
quality threshold = 30
Phred offset      = 33
```

Q30 means Phred score **>= 30**.

## Barcode-round quality

`demux` writes:

```text
barcode_quality_qc.tsv
```

Round1, Round2, and Round3 are reported separately. For each region the table includes:

- observations evaluated;
- bases evaluated;
- percent of bases meeting the threshold;
- percent of complete barcode observations in which every base meets the threshold;
- mean Phred score.

Compact versions of the same metrics are also written to `demux_summary.tsv`.

## Capture-target quality

Capture commands write:

```text
<sample>.capture_quality_qc.tsv
```

`capture-single` evaluates the actual matched target coordinates.

`capture-multi` reports quality for valid flank-extracted targets and for targets accepted by the selected capture mode.

The quality measurement always refers to the sequenced bases in the FASTQ, even when an observed barcode is subsequently corrected to a whitelist entry.
