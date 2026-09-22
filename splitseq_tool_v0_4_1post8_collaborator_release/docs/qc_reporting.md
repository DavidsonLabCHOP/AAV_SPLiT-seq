# QC reporting

`qc-report` combines available RNA, capture, demultiplexing, and paired-cell summaries into user-facing tables and plots.

## RNA metrics

When an RNA count matrix is supplied, the report includes:

- cells in the matrix;
- total RNA counts per cell;
- genes detected per cell;
- optional minimum-count and minimum-gene flags;
- distribution summaries and histograms.

## Capture metrics

When a capture wide-count table is supplied, the report treats values generically as **capture counts**; they may represent reads or unique UMIs depending on the upstream capture command. The report includes:

- total capture counts;
- capture-positive cells;
- capture-positive fraction;
- counts per cell;
- number of positive capture features per cell;
- top capture feature;
- capture count histogram and positive-cell rank plot;
- top capture features for multiplex data.

If an RNA count matrix is also supplied, the capture-positive fraction is calculated within the RNA cell population.

## Paired-cell metrics

`paired-qc` compares RNA and capture cell-ID spaces, including optional capture Round1 -> RNA Round1 conversion, observed overlap, expected random overlap, and observed/expected enrichment.

## Outputs

```text
<sample>.qc_summary.tsv
<sample>.per_cell_qc_metrics.tsv
<sample>.qc_report.md
<sample>.qc_report_plots.pdf
<sample>_qc_plots/
```
