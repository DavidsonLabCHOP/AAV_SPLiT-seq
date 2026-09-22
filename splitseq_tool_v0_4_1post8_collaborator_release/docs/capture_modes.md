# Capture modes

The capture branch uses one consistent term throughout this toolkit: **capture**. It refers to any paired sequence library analyzed alongside the RNA library.

## Single-target capture

Use `capture-single` when one known sequence is the target of interest. The target is always supplied by the user with `--target-sequence`; there is no built-in target sequence.

Optional exact flanks may restrict the search region:

```bash
--use-flanks \
--flank-left AACCGGTTAACC \
--flank-right TTGGCCAATTGG
```

## Multiplex capture

`capture-multi` extracts a sequence between user-supplied exact left and right flanks. Flanks may be any non-zero length.

### Conversion mode

`--capture-mode conversion`

The conversion table contains:

```text
capture_barcode,mapped_sequence
```

The extracted capture barcode is matched/corrected against column 1. The output feature can be the corrected barcode, the mapped sequence, or an amino-acid translation of the mapped sequence.

Conflicting duplicate barcode mappings are rejected. If `--translate-aa` is used, every mapped sequence length must be divisible by three.

### Green-list mode

`--capture-mode greenlist`

Only barcodes present in the approved green list are counted. Mismatch correction may be enabled with `--capture-mismatches`. No sequence mapping or translation occurs.

### Unrestricted mode

`--capture-mode unrestricted`

Every valid A/C/G/T target extracted between the flanks is counted exactly as observed. No whitelist correction, mapping, or translation occurs.

## Target length

`--target-length N` enforces one exact length.

If `--target-length 0`, use `--min-target-length` and `--max-target-length` for a range. In whitelist-based modes, known whitelist lengths are used automatically when no explicit range is provided.

## Mismatch correction

Mismatch correction uses Hamming distance and therefore only compares sequences of equal length. Exact whitelist matches always outrank one-mismatch alternatives. Ties between equally good candidates are handled by `--ambiguous-policy`.

For large whitelists, mismatch values 0 and 1 use a precomputed lookup table. Values greater than 1 require a slower whitelist scan.

## Counting

`--count-mode reads` counts accepted reads.

`--count-mode umi` counts unique UMIs per `(cell ID, capture feature)` pair.

## Cell filtering

`--cell-filter keep-list` restricts capture analysis to a supplied RNA-derived cell list.

`--cell-filter all` keeps every parsed capture cell ID.

When RNA and capture Round1 whitelists differ but correspond row-by-row, `--round1-conversion` converts capture cell IDs into RNA cell-ID space before keep-list filtering.
