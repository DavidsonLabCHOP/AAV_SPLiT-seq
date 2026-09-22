# Validation notes

This repository is a validation-stage toolkit. The core Python logic is covered by lightweight regression tests, but external aligner/counting programs must still be validated on representative datasets in the environment where the toolkit will be used.

## Included regression coverage

The test suite checks:

- exact barcode matches outrank mismatch alternatives;
- ambiguous barcode ties fail under the default policy;
- variable-length green lists do not correct across lengths;
- arbitrary flank lengths and repeated/internal flank-like matches;
- conflicting conversion-table mappings;
- FASTQ truncation and sequence/quality mismatch detection;
- demux cell thresholding, Q30 output, and cleanup;
- conversion, green-list, and unrestricted capture modes;
- YAML boolean CLI overrides;
- pseudoaligner technology-string generation;
- explicit pseudoaligner cell-whitelist filtering and strict UMI-length validation;
- QC capture-count thresholds and zero-count cell retention.

Run:

```bash
PYTHONPATH=. pytest -q
```

## External-program validation still recommended

Validate at least one representative sample for:

- STAR alignment and mapping statistics;
- featureCounts assignment;
- UMI-tools gene/cell matrix dimensions and totals;
- kallisto/bustools matrix dimensions and totals if used;
- capture counts against an independently generated expected set;
- Round1 conversion and paired-cell overlap.

## Known design boundaries

- mismatch correction is substitution-only (Hamming distance);
- capture flanks are exact matches;
- mismatch tolerances greater than 1 are slower for large capture whitelists;
- `raw-kb-design` is experimental because concatenated 24-nt correction does not reproduce independent round-by-round correction exactly;
- legacy compatibility output is retained for older downstream consumers, but long/wide TSV outputs are preferred for new workflows.
