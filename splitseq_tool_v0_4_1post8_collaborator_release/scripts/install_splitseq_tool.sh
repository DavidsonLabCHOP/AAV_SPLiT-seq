#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="splitseq-tool"
WITH_PSEUDOALIGNER="yes"
WITH_STAR="yes"
PYTHON_VERSION="3.10"

usage() {
  cat <<'EOF'
Install SPLiT-seq Tool into a fresh conda/mamba environment.

Usage:
  bash scripts/install_splitseq_tool.sh [options]

Options:
  --env-name NAME           Conda environment name [splitseq-tool]
  --python VERSION          Python version [3.10]
  --no-pseudoaligner        Skip kallisto/bustools/kb-python
  --no-star                 Skip STAR/subread/samtools/umi_tools
  -h, --help                Show this help

Examples:
  bash scripts/install_splitseq_tool.sh
  bash scripts/install_splitseq_tool.sh --env-name splitseq-tool-test
  bash scripts/install_splitseq_tool.sh --no-pseudoaligner
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-name) ENV_NAME="$2"; shift 2 ;;
    --python) PYTHON_VERSION="$2"; shift 2 ;;
    --no-pseudoaligner) WITH_PSEUDOALIGNER="no"; shift ;;
    --no-star) WITH_STAR="no"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

if command -v mamba >/dev/null 2>&1; then
  CONDA_EXE="mamba"
elif command -v conda >/dev/null 2>&1; then
  CONDA_EXE="conda"
else
  echo "ERROR: conda or mamba was not found on PATH."
  echo "Install Miniforge/Mambaforge or load your HPC conda module, then rerun."
  exit 1
fi

TMP_ENV="${REPO_DIR}/.environment.generated.yml"
trap 'rm -f "$TMP_ENV"' EXIT
cat > "$TMP_ENV" <<EOF
name: ${ENV_NAME}
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=${PYTHON_VERSION}
  - pyyaml
  - matplotlib
  - pytest
EOF

if [[ "$WITH_STAR" == "yes" ]]; then
  cat >> "$TMP_ENV" <<'EOF'
  - star
  - subread
  - samtools
  - umi_tools
EOF
fi

if [[ "$WITH_PSEUDOALIGNER" == "yes" ]]; then
  cat >> "$TMP_ENV" <<'EOF'
  - kallisto
  - bustools
  - kb-python
EOF
fi

cat >> "$TMP_ENV" <<'EOF'
  - pip
  - setuptools>=61
EOF

echo "Creating/updating environment: ${ENV_NAME}"
echo "Using: ${CONDA_EXE}"
echo "Environment file: ${TMP_ENV}"

if ${CONDA_EXE} env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  ${CONDA_EXE} env update -n "${ENV_NAME}" -f "$TMP_ENV"
else
  ${CONDA_EXE} env create -f "$TMP_ENV"
fi

echo "Installing splitseq-tool itself without contacting PyPI for build isolation..."
${CONDA_EXE} run -n "${ENV_NAME}" python -m pip install -e "$REPO_DIR" --no-build-isolation --no-deps

cat <<EOF

Done.
Activate with:
  conda activate ${ENV_NAME}

Check installation:
  splitseq-tool --help

Tool versions to record:
  python --version
  splitseq-tool --help
  STAR --version 2>/dev/null || true
  featureCounts -v 2>/dev/null || true
  samtools --version 2>/dev/null | head -2 || true
  umi_tools --version 2>/dev/null || true
  kallisto version 2>/dev/null || true
  bustools version 2>/dev/null || true
EOF
