#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:-project_mini_single_capture_example.yaml}"

echo "======================================================================"
echo "SPLiT-seq Tool YAML-driven mini single-capture benchmark"
echo "======================================================================"
echo "Config: $CONFIG"
echo

echo "=== 1/6 RNA demux ==="
splitseq-tool demux \
  --config "$CONFIG" \
  --sample-id 2_2_mini_rna

echo
echo "=== 2/6 RNA STAR -> featureCounts -> UMI-tools ==="
splitseq-tool rna-star \
  --config "$CONFIG" \
  --sample-id 2_2_mini_rna

echo
echo "=== 3/6 Capture demux ==="
splitseq-tool demux \
  --config "$CONFIG" \
  --sample-id 2_2_mini_capture_demux

echo
echo "=== 4/6 Single-capture counting ==="
splitseq-tool capture-single \
  --config "$CONFIG" \
  --sample-id 2_2_mini

echo
echo "=== 5/6 Paired cell-ID QC ==="
splitseq-tool paired-qc \
  --config "$CONFIG" \
  --sample-id 2_2_mini

echo
echo "=== 6/6 Combined QC report ==="
splitseq-tool qc-report \
  --config "$CONFIG" \
  --sample-id 2_2_mini

echo
echo "======================================================================"
echo "Done."
echo "======================================================================"
