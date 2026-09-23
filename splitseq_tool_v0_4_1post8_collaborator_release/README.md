# splitseq-tool

`splitseq-tool` is a command-line toolkit for SPLiT-seq RNA processing and paired **capture** sequence analysis. The capture branch is intentionally generic: it can quantify one known target, a whitelist of approved targets, every flank-delimited target observed in the data, or targets that map through a user-supplied conversion table.

**Version:** `0.4.1.post8`  

> All experiment-specific target, flank, and example sequences shown in this repository are synthetic placeholders used only for documentation/tests. Replace them with the sequences appropriate for your experiment.

## Main commands

```text
splitseq-tool demux                 SPLiT-seq barcode correction and cell calling
splitseq-tool rna-star              STAR + featureCounts + UMI-tools RNA workflow
splitseq-tool resume-featurecounts  resume RNA processing from an aligned BAM
splitseq-tool prepare-kb-fastqs     create pseudoaligner-compatible FASTQs after demux
splitseq-tool rna-kb                kallisto/bustools RNA workflow after demux
splitseq-tool raw-kb-design         experimental direct-raw pseudoaligner design helper
splitseq-tool capture-single        quantify one capture target
splitseq-tool capture-multi         quantify flank-delimited capture targets
splitseq-tool paired-qc             compare RNA and capture cell-ID spaces
splitseq-tool qc-report             combine run summaries, matrices, and capture counts
```

Run command-specific help at any time:

```bash
splitseq-tool <command> --help
```

## Installation

The recommended installation uses conda or mamba and installs the local package without requiring PyPI during the final package-install step.

```bash
unzip splitseq_tool_v0_4_1post8_collaborator_release.zip
cd splitseq_tool_v0_4_1post8_collaborator_release
bash scripts/install_splitseq_tool.sh
conda activate splitseq-tool
splitseq-tool --help
```

For an already-created `splitseq-tool` environment:

```bash
python -m pip install -e . --no-build-isolation --no-deps
```

The full environment includes Python 3.10, STAR, featureCounts, samtools, UMI-tools, kallisto, bustools, kb-python, PyYAML, matplotlib, and pytest.

## Core barcode model

By default, the barcode read uses zero-based coordinates:

```text
UMI      0:10
Round 3  10:18
Round 2  48:56
Round 1  86:94
```

The corrected 24-nt cell ID is assembled as:

```text
Round1 + Round2 + Round3
```

The coordinate options are configurable. Barcode correction is performed separately for each 8-nt round. With `--barcode-mismatches 1`, each round can be corrected by up to one substitution independently.

RNA and capture libraries may use separate Round1 barcode lists. When `--round1-conversion` is enabled, capture Round1 barcodes are converted to RNA Round1 space by corresponding row in the two separate Round1 lists.

## Recommended paired workflow

A typical paired analysis is:

```text
RNA FASTQs
  -> demux
  -> RNA quantification (rna-star or rna-kb)
  -> final RNA cellIDs_keep.txt

capture FASTQs
  -> demux with the capture Round1 whitelist
  -> capture-single or capture-multi
       using the final RNA cell keep-list when desired
       and optional capture-Round1 -> RNA-Round1 conversion

RNA + capture
  -> paired-qc
  -> qc-report
```

The RNA branch should finish before a capture analysis that uses `--cell-filter keep-list`, because the final RNA keep-list is an input to that capture step.

## Demultiplexing

Example RNA demultiplexing:

```bash
splitseq-tool demux \
  --sample-id sample_01 \
  --r1-fastq /path/to/sample_01_R1.fastq.gz \
  --r2-fastq /path/to/sample_01_R2.fastq.gz \
  --round1 /path/to/Round1_RNA_barcodes.txt \
  --round2 /path/to/Round2_barcodes.txt \
  --round3 /path/to/Round3_barcodes.txt \
  --barcode-mismatches 1 \
  --cell-threshold-mode umis \
  --min-cell-count 100 \
  --output-dir results/sample_01/rna_demux
```

Important outputs include:

```text
MergedCells_passing.fastq
cellIDs_keep.txt
cell_barcode_umi_counts.tsv
demux_summary.tsv
barcode_quality_qc.tsv
resolved_parameters.tsv
resolved_config.yaml
cleanup_summary.tsv
```

Both `.fastq` and `.fastq.gz` inputs are supported.

## RNA quantification with STAR

```bash
splitseq-tool rna-star \
  --sample-id sample_01 \
  --merged-fastq results/sample_01/rna_demux/MergedCells_passing.fastq \
  --star-index /path/to/STAR_index \
  --saf-file /path/to/genes.saf \
  --threads 10 \
  --output-dir results/sample_01/rna
```

The RNA workflow normalizes FASTQ headers so corrected cell IDs and UMIs survive alignment, then runs STAR, featureCounts, assigned-read filtering, and UMI-tools. The final wide RNA count matrix is `counts.tsv.gz`; the cell columns are also written to `cellIDs_keep.txt`.

## Capture: single known target

`capture-single` quantifies one user-specified target sequence. There is no built-in biological target sequence.

```bash
splitseq-tool capture-single \
  --sample-id sample_01 \
  --input-fastq results/sample_01/capture_demux/MergedCells_passing.fastq \
  --target-sequence ACGTACGTACGT \
  --cell-keep-list results/sample_01/rna/sample_01/rna_star/cellIDs_keep.txt \
  --cell-filter keep-list \
  --round1-conversion \
  --mrna-round1 /path/to/Round1_RNA_barcodes.txt \
  --capture-round1 /path/to/Round1_capture_barcodes.txt \
  --capture-mismatches 1 \
  --count-mode reads \
  --output-dir results/sample_01/capture_single
```

To restrict the target search to a flank-delimited region, add explicit flanks:

```bash
--use-flanks \
--flank-left AACCGGTTAACC \
--flank-right TTGGCCAATTGG
```

Flank sequences are never assumed by the tool.

## Capture: multiplex targets

`capture-multi` always extracts a target between an explicitly supplied left and right flank. Flanks may be any non-zero length and may have different lengths.

### Mode 1: conversion

Use a two-column mapping table in which column 1 is the approved capture barcode and column 2 is the mapped nucleotide sequence.

```bash
splitseq-tool capture-multi \
  --sample-id sample_01 \
  --input-fastq /path/to/MergedCells_passing.fastq \
  --output-dir results/sample_01/capture_multi \
  --cell-filter all \
  --capture-mode conversion \
  --barcode-conversion-table /path/to/capture_barcode_map.csv \
  --flank-left AACCGGTTAACC \
  --flank-right TTGGCCAATTGG \
  --target-length 12 \
  --capture-mismatches 0 \
  --count-mode reads
```

Without translation, `--output-sequence barcode` reports the corrected capture barcode and `--output-sequence mapped` reports the mapped nucleotide sequence. The legacy aliases `rna` and `dna` remain accepted for backward compatibility.

If the mapped sequence is coding and every mapped sequence length is divisible by three, add:

```bash
--translate-aa
```

### Mode 2: green list

Count only approved capture barcodes, with no mapping or translation:

```bash
splitseq-tool capture-multi \
  --capture-mode greenlist \
  --green-list /path/to/approved_capture_barcodes.txt \
  --flank-left AACCGGTTAACC \
  --flank-right TTGGCCAATTGG \
  --target-length 12 \
  --capture-mismatches 1 \
  ...
```

`--green-list-only` is a shortcut for `--capture-mode greenlist`.

### Mode 3: unrestricted

Count every valid A/C/G/T target observed between the flanks exactly as sequenced:

```bash
splitseq-tool capture-multi \
  --capture-mode unrestricted \
  --flank-left AACCGGTTAACC \
  --flank-right TTGGCCAATTGG \
  --target-length 12 \
  ...
```

`--green-list-free` is a shortcut for `--capture-mode unrestricted`.

Unrestricted mode can produce a large number of unique features when sequencing errors are present; the long-format output is usually the easiest representation to inspect.

### Capture target length

Use an exact target length:

```bash
--target-length 12
```

or a range:

```bash
--target-length 0 \
--min-target-length 8 \
--max-target-length 20
```

If all three length controls are zero, any non-empty sequence between the flanks is accepted in unrestricted mode. In whitelist-based modes, known whitelist lengths are used automatically when no explicit length range is supplied.

## Read counts versus UMI counts

Capture modules support:

```bash
--count-mode reads
```

or:

```bash
--count-mode umi
```

Read counting is the default. UMI counting deduplicates unique UMIs separately for each `(cell ID, capture feature)` pair and requires a valid UMI in the demultiplexed FASTQ header.

## Quality-control outputs

Demultiplexing automatically writes `barcode_quality_qc.tsv` with separate Round1, Round2, and Round3 metrics. For each round it reports:

- percent of individual bases with Phred score >= the threshold;
- percent of complete barcode observations for which every base meets the threshold;
- mean Phred score.

The default threshold is Q30 (`--quality-threshold 30`, Phred+33).

Capture commands write `<sample>.capture_quality_qc.tsv`. Single-target mode evaluates the matched target bases. Multiplex mode reports quality for extracted targets and for targets accepted by the selected capture mode.

See `docs/quality_qc.md`.

## Automatic intermediate cleanup

Large staging files are removed by default only after the command completes successfully and expected final outputs are verified. A `cleanup_summary.tsv` records what was deleted.

To retain intermediates:

```bash
--keep-all-intermediates
```

or:

```bash
--keep_all_intermediates
```

Boolean options also support `--no-...` forms, which is useful for overriding YAML values from the command line.

See `docs/intermediate_cleanup.md`.

## Optional YAML configuration

Every major command can use flags only, YAML only, or YAML plus explicit overrides. CLI values take precedence over YAML values.

```bash
splitseq-tool capture-single \
  --config configs/example_project.yaml \
  --sample-id sample_01 \
  --target-sequence ACGTACGTACGT
```

Example configurations are in `configs/`. See `docs/configuration.md`.

Each run writes resolved settings to:

```text
resolved_parameters.tsv
resolved_config.yaml
```

## QC report

`qc-report` can combine demux summaries, RNA count matrices, capture count tables, and paired cell-ID summaries.

```bash
splitseq-tool qc-report \
  --sample-id sample_01 \
  --output-dir results/sample_01/qc \
  --demux-summary results/sample_01/rna_demux/demux_summary.tsv \
  --rna-counts-matrix results/sample_01/rna/sample_01/rna_star/counts.tsv.gz \
  --capture-counts-wide results/sample_01/capture_single/sample_01.single_capture_counts.wide.tsv
```

Outputs include:

```text
<sample>.qc_summary.tsv
<sample>.per_cell_qc_metrics.tsv
<sample>.qc_report.md
<sample>.qc_report_plots.pdf
<sample>_qc_plots/
```

See `docs/qc_reporting.md`.

## Legacy compatibility output

Capture commands can optionally write the historical per-cell compatibility format with:

```bash
--legacy-output
```

This preserves the historical directory/file structure for downstream tools that still require it. New analyses should preferentially use the long and wide TSV outputs.

## Pseudoalignment

The preferred pseudoalignment route is:

```text
demux -> rna-kb
```

The toolkit first performs its own round-specific SPLiT-seq barcode correction, then creates synthetic barcode/UMI and cDNA FASTQs for kb/kallisto-bustools. This keeps barcode correction semantics consistent with the STAR workflow.

`raw-kb-design` is experimental because correcting one concatenated 24-nt cell barcode is not equivalent to correcting the three 8-nt barcode rounds independently.

## Testing the installation

The repository includes regression tests for core Python logic:

```bash
PYTHONPATH=. pytest -q
```

External STAR/featureCounts/UMI-tools and kallisto/bustools executions still require validation on representative real datasets and are not fully exercised by the lightweight unit tests.

## Documentation

- `docs/capture_modes.md`
- `docs/configuration.md`
- `docs/quality_qc.md`
- `docs/intermediate_cleanup.md`
- `docs/qc_reporting.md`
- `docs/validation.md`
- `CHANGELOG.md`

## Running the small benchmarking dataset

Small example datasets are included with splitseq-tool (in mini-data directory) so that users can test the complete pipeline before running their own sequencing data. These files contain only small number of synthetic reads and are are intended for verification that the software, configuration files, barcode handling, and alignment steps are working correctly.

1. Set up and activate the environment, and install the tool as described above.
2. download the mouse reference genome and annotation
   - The example datasets were generated using the GRCm38 (mm10) mouse reference genome. Both files readily available on GENCODE, Ensembl and other platforms.
   - If the files are compressed, decompress them before building the STAR index:
   ```bash
   gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
   gunzip Mus_musculus.GRCm38.*.gtf.gz
   ```
3. build the STAR genome index
   - By default, splitseq-tool uses STAR for alignment. The STAR genome index only needs to be generated once and can then be reused for all datasets analyzed with the same genome/annotation
   - Then simply run this:
   ```bash
   STAR \
    --runMode genomeGenerate \
    --runThreadN 8 \  #adjust number of threads as necessary
    --genomeDir GRCm38_STAR_index \  
    --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa \  # this name might change depending of platform where you obtained the genome
    --sjdbGTFfile Mus_musculus.GRCm38.<release>.gtf \ # this name might change depending of platform where you obtained the genome
    --sjdbOverhang 99
   ```
   - --runThreadN can be adjusted for the number of CPU cores available on your system. The --sjdbOverhang value should ideally be read length − 1; 99 is appropriate for 100-nt reads. If your reads have a different length, this value can be changed accordingly.
4. Locate example project YAML (e.g. project_mini_single_capture_example.yaml) file provided and replace any generic paths with your actual paths as necessary.
5. run the supplied mini dataset using the provided bash script and edited YAML:
   ```bash
   bash run_project_mini_single_capture_from_yaml.sh
   ```
   or
   ```bash
   ./run_project_mini_single_capture_from_yaml.sh project_mini_single_capture_example.yaml
   ```
6. The analysis should process very quickly and you should see following message:
```bash
======================================================================
SPLiT-seq Tool YAML-driven mini single-capture benchmark
======================================================================
Config: project_mini_single_capture_example.yaml

=== 1/6 RNA demux ===

=== 2/6 RNA STAR -> featureCounts -> UMI-tools ===

=== 3/6 Capture demux ===

=== 4/6 Single-capture counting ===

=== 5/6 Paired cell-ID QC ===

=== 6/6 Combined QC report ===

======================================================================
Done.
======================================================================
 
   
