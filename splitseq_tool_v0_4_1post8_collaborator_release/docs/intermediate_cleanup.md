# Intermediate-file cleanup

Large staging files are removed by default only after the current command completes successfully and expected final outputs have been checked.

To retain intermediates:

```bash
--keep-all-intermediates
```

or:

```bash
--keep_all_intermediates
```

Each cleanup-enabled command writes:

```text
cleanup_summary.tsv
```

The summary lists each path, whether it was deleted or retained, its size, and the reason.

Current cleanup behavior:

- `demux`: removes `MergedCells_1.fastq` after final demux outputs are created;
- `rna-star`: removes the rewritten STAR-input FASTQ and large BAM staging files after the final count matrix and cell list exist;
- `resume-featurecounts`: removes only files created by that resume command and never deletes the user-supplied aligned BAM;
- `rna-kb`: removes synthetic pseudoaligner FASTQs after a successful pseudoalignment run;
- capture/QC commands do not currently create large disposable staging files.

If a command fails before completion, automatic cleanup is not performed.
