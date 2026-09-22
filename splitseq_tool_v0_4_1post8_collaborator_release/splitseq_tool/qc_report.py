from __future__ import annotations

import csv
import gzip
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from .utils import ensure_dir, write_kv


def _open_text(path):
    path = str(path)
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def _to_number(x):
    if x is None or x == '':
        return None
    try:
        if isinstance(x, (int, float)):
            return x
        s = str(x).strip()
        if s.lower() in {'na', 'nan', 'none'}:
            return None
        if any(c in s for c in ['.', 'e', 'E']):
            return float(s)
        return int(s)
    except Exception:
        return x


def read_kv_table(path: str | Path | None) -> Dict[str, object]:
    d: Dict[str, object] = {}
    if not path:
        return d
    p = Path(path)
    if not p.exists():
        return d
    with open(p) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or '\t' not in line:
                continue
            k, v = line.split('\t', 1)
            d[k] = _to_number(v)
    return d


def read_capture_wide(path: str | Path | None) -> Tuple[Dict[str, Dict[str, int]], List[str]]:
    """Read cell x capture-barcode count table.

    Expected columns: cellID, <sequence>_count, ...
    Returns: cell -> {sequence_label: count}
    """
    data: Dict[str, Dict[str, int]] = {}
    labels: List[str] = []
    if not path:
        return data, labels
    p = Path(path)
    if not p.exists():
        return data, labels
    with open(p) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None)
        if not header:
            return data, labels
        labels = [h[:-6] if h.endswith('_count') else h for h in header[1:]]
        for row in reader:
            if not row:
                continue
            cell = row[0]
            vals: Dict[str, int] = {}
            for lab, val in zip(labels, row[1:]):
                try:
                    vals[lab] = int(float(val))
                except Exception:
                    vals[lab] = 0
            data[cell] = vals
    return data, labels


def summarize_capture_counts(data: Dict[str, Dict[str, int]], min_capture_counts: int = 1) -> Tuple[List[Tuple[str, object]], Dict[str, Dict[str, object]]]:
    per_cell: Dict[str, Dict[str, object]] = {}
    total_counts = 0
    positive_cells = 0
    cells_ge_threshold = 0
    barcode_positive_counts = Counter()
    feature_count_totals = Counter()

    for cell, vals in data.items():
        total = sum(vals.values())
        n_pos = sum(1 for v in vals.values() if v > 0)
        top_label = ''
        top_count = 0
        if vals:
            top_label, top_count = max(vals.items(), key=lambda kv: kv[1])
        total_counts += total
        if total > 0:
            positive_cells += 1
        if total >= min_capture_counts:
            cells_ge_threshold += 1
        for label, c in vals.items():
            if c > 0:
                barcode_positive_counts[label] += 1
                feature_count_totals[label] += c
        per_cell[cell] = {
            'cellID': cell,
            'capture_total_counts': total,
            'capture_positive_features': n_pos,
            'capture_top_feature': top_label,
            'capture_top_feature_count': top_count,
            'capture_pass_min_counts': total >= min_capture_counts,
        }

    n_cells = len(data)
    rows = [
        ('capture_cells_in_table', n_cells),
        ('capture_total_counts', total_counts),
        ('capture_positive_cells', positive_cells),
        ('capture_positive_cells_pct', round(100 * positive_cells / n_cells, 4) if n_cells else 0),
        ('capture_min_counts_filter', min_capture_counts),
        ('capture_cells_passing_min_counts', cells_ge_threshold),
        ('capture_cells_passing_min_counts_pct', round(100 * cells_ge_threshold / n_cells, 4) if n_cells else 0),
        ('capture_mean_counts_per_cell_all', round(total_counts / n_cells, 4) if n_cells else 0),
        ('capture_mean_counts_per_positive_cell', round(total_counts / positive_cells, 4) if positive_cells else 0),
        ('capture_features_detected', len([k for k, v in feature_count_totals.items() if v > 0])),
    ]
    if feature_count_totals:
        top = feature_count_totals.most_common(1)[0]
        rows += [('capture_top_feature_by_counts', top[0]), ('capture_top_feature_count', top[1])]
    return rows, per_cell


def read_gene_by_cell_matrix(path: str | Path | None, max_genes: Optional[int] = None) -> Dict[str, Dict[str, int]]:
    """Read a simple tab-delimited gene x cell count matrix.

    Works for wide umi_tools-style matrices when first column is gene and the remaining columns are cells.
    It intentionally avoids scipy/anndata dependencies in this first implementation.
    """
    metrics: Dict[str, Dict[str, int]] = {}
    if not path:
        return metrics
    p = Path(path)
    if not p.exists():
        return metrics
    with _open_text(p) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None)
        if not header or len(header) < 2:
            return metrics
        cells = header[1:]
        for cell in cells:
            metrics[cell] = {'rna_total_counts': 0, 'rna_genes_detected': 0}
        for gi, row in enumerate(reader):
            if max_genes is not None and gi >= max_genes:
                break
            if len(row) < 2:
                continue
            for cell, val in zip(cells, row[1:]):
                try:
                    c = int(float(val))
                except Exception:
                    c = 0
                if c > 0:
                    metrics[cell]['rna_total_counts'] += c
                    metrics[cell]['rna_genes_detected'] += 1
    return metrics


def percentile(vals: List[float], q: float):
    if not vals:
        return 0
    vals = sorted(vals)
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * q / 100
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    return vals[lo] * (hi - pos) + vals[hi] * (pos - lo)


def numeric_summary(prefix: str, vals: Iterable[float]) -> List[Tuple[str, object]]:
    arr = [float(v) for v in vals if v is not None]
    if not arr:
        return [(f'{prefix}_n', 0)]
    return [
        (f'{prefix}_n', len(arr)),
        (f'{prefix}_min', min(arr)),
        (f'{prefix}_p25', round(percentile(arr, 25), 4)),
        (f'{prefix}_median', round(percentile(arr, 50), 4)),
        (f'{prefix}_mean', round(sum(arr) / len(arr), 4)),
        (f'{prefix}_p75', round(percentile(arr, 75), 4)),
        (f'{prefix}_max', max(arr)),
    ]


def write_per_cell(per_cell: Dict[str, Dict[str, object]], path: str | Path) -> None:
    fields = [
        'cellID',
        'rna_total_counts', 'rna_genes_detected', 'rna_pass_min_counts', 'rna_pass_min_genes',
        'capture_total_counts', 'capture_positive_features', 'capture_top_feature', 'capture_top_feature_count', 'capture_pass_min_counts',
        'capture_positive',
    ]
    with open(path, 'w') as out:
        out.write('\t'.join(fields) + '\n')
        for cell in sorted(per_cell):
            row = per_cell[cell]
            row.setdefault('cellID', cell)
            row.setdefault('rna_total_counts', '')
            row.setdefault('rna_genes_detected', '')
            row.setdefault('rna_pass_min_counts', '')
            row.setdefault('rna_pass_min_genes', '')
            row.setdefault('capture_total_counts', 0)
            row.setdefault('capture_positive_features', 0)
            row.setdefault('capture_top_feature', '')
            row.setdefault('capture_top_feature_count', 0)
            row.setdefault('capture_pass_min_counts', False)
            row['capture_positive'] = int(row.get('capture_total_counts') or 0) > 0
            out.write('\t'.join(str(row.get(f, '')) for f in fields) + '\n')


def make_plots(out_dir: Path, sample_id: str, per_cell: Dict[str, Dict[str, object]], capture_labels: List[str], capture_data: Dict[str, Dict[str, int]]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except Exception as e:
        with open(out_dir / f'{sample_id}.qc_plot_error.txt', 'w') as out:
            out.write(f'Could not import matplotlib: {e}\n')
        return

    plot_dir = ensure_dir(out_dir / f'{sample_id}_qc_plots')
    pdf_path = out_dir / f'{sample_id}.qc_report_plots.pdf'

    def savefig(name):
        plt.tight_layout()
        plt.savefig(plot_dir / f'{name}.png', dpi=200)
        pdf.savefig()
        plt.close()

    with PdfPages(pdf_path) as pdf:
        vals = [int(r.get('rna_total_counts') or 0) for r in per_cell.values() if r.get('rna_total_counts') != '']
        if vals:
            plt.figure()
            plt.hist(vals, bins=60)
            plt.xlabel('RNA counts / UMIs per cell')
            plt.ylabel('Cells')
            plt.title(f'{sample_id}: RNA counts per cell')
            savefig('rna_counts_per_cell_hist')

        vals = [int(r.get('rna_genes_detected') or 0) for r in per_cell.values() if r.get('rna_genes_detected') != '']
        if vals:
            plt.figure()
            plt.hist(vals, bins=60)
            plt.xlabel('Genes detected per cell')
            plt.ylabel('Cells')
            plt.title(f'{sample_id}: genes per cell')
            savefig('genes_per_cell_hist')

        cvals = [int(r.get('capture_total_counts') or 0) for r in per_cell.values()]
        if cvals:
            plt.figure()
            plt.hist(cvals, bins=60)
            plt.xlabel('Capture counts per cell')
            plt.ylabel('Cells')
            plt.title(f'{sample_id}: capture counts per cell')
            savefig('capture_counts_per_cell_hist')

            # rank plot for positive cells
            pos = sorted([v for v in cvals if v > 0], reverse=True)
            if pos:
                plt.figure()
                plt.plot(range(1, len(pos) + 1), pos, marker='o', linewidth=1)
                plt.xlabel('Positive-cell rank')
                plt.ylabel('Capture counts')
                plt.yscale('log')
                plt.title(f'{sample_id}: capture-positive cell rank plot')
                savefig('capture_positive_cell_rank')

        if capture_labels and capture_data and len(capture_labels) > 1:
            feature_reads = []
            feature_cells = []
            for lab in capture_labels:
                reads = sum(vals.get(lab, 0) for vals in capture_data.values())
                cells = sum(1 for vals in capture_data.values() if vals.get(lab, 0) > 0)
                feature_reads.append((lab, reads))
                feature_cells.append((lab, cells))
            feature_reads = sorted(feature_reads, key=lambda x: x[1], reverse=True)[:30]
            if feature_reads:
                labs = [x[0] for x in feature_reads]
                vals = [x[1] for x in feature_reads]
                plt.figure(figsize=(max(6, len(labs) * 0.25), 4))
                plt.bar(range(len(labs)), vals)
                plt.xticks(range(len(labs)), labs, rotation=90)
                plt.ylabel('Total capture counts')
                plt.title(f'{sample_id}: top capture features by counts')
                savefig('top_capture_features_by_counts')


def run_qc_report(args) -> None:
    out_dir = ensure_dir(args.output_dir)
    sample = args.sample_id or 'sample'

    summary: List[Tuple[str, object]] = [('sample_id', sample)]

    # Existing stage summaries, if supplied.
    for label, path in [
        ('demux', getattr(args, 'demux_summary', None)),
        ('rna', getattr(args, 'rna_summary', None)),
        ('capture', getattr(args, 'capture_summary', None)),
        ('paired', getattr(args, 'paired_summary', None)),
    ]:
        d = read_kv_table(path)
        if d:
            for k, v in d.items():
                summary.append((f'{label}_{k}', v))

    capture_data, capture_labels = read_capture_wide(getattr(args, 'capture_counts_wide', None))
    capture_rows, capture_per_cell = summarize_capture_counts(capture_data, getattr(args, 'min_capture_count', 1))
    summary.extend(capture_rows)

    rna_per_cell = read_gene_by_cell_matrix(getattr(args, 'rna_counts_matrix', None))
    per_cell: Dict[str, Dict[str, object]] = {}
    for cell, vals in rna_per_cell.items():
        per_cell.setdefault(cell, {'cellID': cell}).update(vals)
    for cell, vals in capture_per_cell.items():
        per_cell.setdefault(cell, {'cellID': cell}).update(vals)

    # RNA filters / summaries.
    if rna_per_cell:
        min_counts = getattr(args, 'min_rna_counts', 0)
        min_genes = getattr(args, 'min_genes', 0)
        pass_counts = 0
        pass_genes = 0
        for cell, row in per_cell.items():
            if cell in rna_per_cell:
                rc = int(row.get('rna_total_counts') or 0)
                ng = int(row.get('rna_genes_detected') or 0)
                row['rna_pass_min_counts'] = rc >= min_counts
                row['rna_pass_min_genes'] = ng >= min_genes
                if rc >= min_counts:
                    pass_counts += 1
                if ng >= min_genes:
                    pass_genes += 1
        n = len(rna_per_cell)
        summary += [
            ('rna_cells_in_count_matrix', n),
            ('rna_min_counts_filter', min_counts),
            ('rna_cells_passing_min_counts', pass_counts),
            ('rna_cells_passing_min_counts_pct', round(100 * pass_counts / n, 4) if n else 0),
            ('rna_min_genes_filter', min_genes),
            ('rna_cells_passing_min_genes', pass_genes),
            ('rna_cells_passing_min_genes_pct', round(100 * pass_genes / n, 4) if n else 0),
        ]
        summary.extend(numeric_summary('rna_counts_per_cell', [v['rna_total_counts'] for v in rna_per_cell.values()]))
        summary.extend(numeric_summary('rna_genes_per_cell', [v['rna_genes_detected'] for v in rna_per_cell.values()]))

    if capture_data:
        summary.extend(numeric_summary('capture_counts_per_cell', [sum(v.values()) for v in capture_data.values()]))
        # Capture-positive fraction. When an RNA matrix is available, both the
        # numerator and denominator are restricted to RNA cells so the percentage
        # cannot be inflated by capture-only rows.
        if rna_per_cell:
            denom = len(rna_per_cell)
            pos = sum(
                1 for cell in rna_per_cell
                if cell in capture_data and sum(capture_data[cell].values()) > 0
            )
            denom_label = 'rna_count_matrix_cells'
        else:
            denom = len(capture_data)
            pos = sum(1 for v in capture_data.values() if sum(v.values()) > 0)
            denom_label = 'capture_count_table_cells'
        summary.append(('capture_enrichment_positive_cells_over_denominator_pct', round(100 * pos / denom, 4) if denom else 0))
        summary.append(('capture_enrichment_denominator', denom_label))

    write_kv(out_dir / f'{sample}.qc_summary.tsv', summary)
    write_per_cell(per_cell, out_dir / f'{sample}.per_cell_qc_metrics.tsv')
    make_plots(out_dir, sample, per_cell, capture_labels, capture_data)

    # Short text report.
    with open(out_dir / f'{sample}.qc_report.md', 'w') as out:
        out.write(f'# SPLiT-seq QC report: {sample}\n\n')
        out.write('## Key summary metrics\n\n')
        for k, v in summary:
            out.write(f'- **{k}**: {v}\n')
        out.write('\n## Files written\n\n')
        out.write(f'- `{sample}.qc_summary.tsv`\n')
        out.write(f'- `{sample}.per_cell_qc_metrics.tsv`\n')
        out.write(f'- `{sample}.qc_report_plots.pdf` if matplotlib was available\n')
