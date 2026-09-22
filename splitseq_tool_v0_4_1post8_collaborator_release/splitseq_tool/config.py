from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional

import yaml

from .utils import ensure_dir

MISSING = object()


def load_yaml_config(path: str | Path | None) -> Dict[str, Any]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Config file not found: {p}")
    with open(p, "r") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config must contain a top-level mapping/dictionary: {p}")
    data["_config_path"] = str(p)
    return data


def _deep_get(d: Mapping[str, Any], path: str, default: Any = MISSING) -> Any:
    cur: Any = d
    for part in path.split("."):
        if not isinstance(cur, Mapping) or part not in cur:
            return default
        cur = cur[part]
    return cur


def _sample_block(config: Mapping[str, Any], sample_id: str | None) -> Dict[str, Any]:
    if not sample_id:
        return {}
    samples = config.get("samples", {})
    if isinstance(samples, Mapping):
        # Match by string value as well as exact key so sample identifiers remain user-friendly.
        block = samples.get(sample_id, MISSING)
        if block is MISSING:
            for k, v in samples.items():
                if str(k) == str(sample_id):
                    block = v
                    break
        if block is MISSING or block is None:
            return {}
        if not isinstance(block, Mapping):
            raise ValueError(f"YAML samples.{sample_id} must be a mapping")
        return dict(block)
    if isinstance(samples, list):
        for item in samples:
            if isinstance(item, Mapping) and str(item.get("id", "")) == str(sample_id):
                return dict(item)
    return {}


def _first_present(*vals: Any) -> Any:
    for v in vals:
        if v is not MISSING and v is not None:
            return v
    return MISSING


def _lookup_config_value(config: Mapping[str, Any], sample: Mapping[str, Any], cmd: str, arg_name: str) -> Any:
    """Lookup one argparse destination in sample YAML, command YAML, or project defaults.

    Precedence handled by caller:
      CLI explicit > sample-specific > command-specific > project defaults > built-in default
    """
    dashed = arg_name.replace("_", "-")
    candidates = []

    # Sample-specific values can be flat or nested by command / library type.
    candidates += [
        _deep_get(sample, arg_name),
        _deep_get(sample, dashed),
        _deep_get(sample, f"{cmd}.{arg_name}"),
        _deep_get(sample, f"{cmd}.{dashed}"),
    ]

    # Convenient aliases for FASTQs / inputs in samples.
    alias_paths = {
        "r1_fastq": ["mrna_r1", "rna_r1", "r1", "fastqs.r1", "fastqs.mrna_r1"],
        "r2_fastq": ["mrna_r2", "rna_r2", "r2", "fastqs.r2", "fastqs.mrna_r2"],
        "mrna_fastq": ["mrna_merged_fastq", "rna_merged_fastq", "rna.merged_fastq", "mrna.merged_fastq"],
        "capture_fastq": ["capture_merged_fastq", "capture.merged_fastq"],
        "input_fastq": ["capture_merged_fastq", "merged_fastq", "capture.merged_fastq"],
        "merged_fastq": ["merged_fastq", "mrna_merged_fastq", "rna_merged_fastq", "rna.merged_fastq", "mrna.merged_fastq"],
        "cell_whitelist": ["cell_whitelist", "cellIDs_keep", "mrna_cell_keep_list", "rna_cell_keep_list", "mrna.cell_keep_list", "rna.cell_keep_list"],
        "cell_keep_list": ["cell_keep_list", "cellIDs_keep", "mrna_cell_keep_list", "rna_cell_keep_list", "mrna.cell_keep_list", "rna.cell_keep_list"],
        "aligned_bam": ["aligned_bam", "rna.aligned_bam", "mrna.aligned_bam"],
        "output_dir": ["output_dir", f"outputs.{cmd}", "results_dir"],
    }
    for a in alias_paths.get(arg_name, []):
        candidates.append(_deep_get(sample, a))

    # Command-specific config blocks.
    command_aliases = {
        "capture-single": ["capture.single", "capture_single", "single_capture"],
        "capture-multi": ["capture.multi", "capture_multi", "multi_capture"],
        "paired-qc": ["paired_qc", "qc.paired", "paired"],
        "qc-report": ["qc_report", "qc.report", "report", "summary"],
        "rna-star": ["rna.star", "rna_star", "star"],
        "rna-kb": ["rna.kb", "rna_kb", "kb", "pseudoalign", "pseudoaligner"],
        "prepare-kb-fastqs": ["rna.prepare_kb_fastqs", "prepare_kb_fastqs", "kb.prepare_fastqs"],
        "raw-kb-design": ["rna.raw_kb_design", "raw_kb_design", "kb.raw_design"],
        "resume-featurecounts": ["rna.resume_featurecounts", "resume_featurecounts"],
        "demux": ["demux"],
    }
    for base in command_aliases.get(cmd, [cmd]):
        candidates += [_deep_get(config, f"{base}.{arg_name}"), _deep_get(config, f"{base}.{dashed}")]

    # Common project-level blocks.
    project_aliases = {
        "round1": ["barcodes.round1", "barcodes.mrna_round1"],
        "round2": ["barcodes.round2"],
        "round3": ["barcodes.round3"],
        "mrna_round1": ["barcodes.mrna_round1", "barcodes.round1_mrna"],
        "capture_round1": ["barcodes.capture_round1", "barcodes.round1_capture"],
        "star_index": ["references.star_index", "rna.star_index"],
        "saf_file": ["references.saf", "references.saf_file", "rna.saf_file"],
        "kallisto_index": ["references.kallisto_index", "references.kb_index", "rna.kallisto_index", "rna.kb.index"],
        "t2g_file": ["references.t2g", "references.t2g_file", "references.transcript_to_gene", "rna.t2g_file", "rna.kb.t2g_file"],
        "barcode_conversion_table": ["capture.barcode_conversion_table", "capture.multi.barcode_conversion_table"],
        "green_list": ["capture.green_list", "capture.multi.green_list", "capture.greenlist"],
        "capture_mode": ["capture.mode", "capture.multi.mode", "capture.capture_mode"],
        "target_length": ["capture.target_length", "capture.multi.target_length"],
        "flank_left": ["capture.flank_left", "capture.single.flank_left", "capture.multi.flank_left"],
        "flank_right": ["capture.flank_right", "capture.single.flank_right", "capture.multi.flank_right"],
        "target_sequence": ["capture.single.target_sequence", "capture.target_sequence"],
        "capture_counts_wide": ["qc.capture_counts_wide", "capture.counts_wide", "capture.wide_counts"],
        "min_capture_count": ["qc.min_capture_count", "qc.min_capture_reads", "min_capture_reads"],
        "rna_counts_matrix": ["qc.rna_counts_matrix", "rna.counts_matrix", "rna.matrix"],
        "demux_summary": ["qc.demux_summary", "demux.summary"],
        "rna_summary": ["qc.rna_summary", "rna.summary"],
        "capture_summary": ["qc.capture_summary", "capture.summary"],
        "paired_summary": ["qc.paired_summary", "paired_qc.summary"],
    }
    for a in project_aliases.get(arg_name, []):
        candidates.append(_deep_get(config, a))

    # Defaults block can be flat or command-specific.
    candidates += [
        _deep_get(config, f"defaults.{cmd}.{arg_name}"),
        _deep_get(config, f"defaults.{cmd}.{dashed}"),
        _deep_get(config, f"defaults.{arg_name}"),
        _deep_get(config, f"defaults.{dashed}"),
    ]

    return _first_present(*candidates)


def _arg_was_explicit(argv: Optional[list[str]], opt_strings: Iterable[str]) -> bool:
    if argv is None:
        return False
    opts = set(opt_strings)
    for token in argv:
        if token in opts:
            return True
        # Handles --option=value
        if token.startswith("--") and "=" in token and token.split("=", 1)[0] in opts:
            return True
    return False


def explicit_cli_dests(parser: argparse.ArgumentParser, argv: Optional[list[str]]) -> set[str]:
    """Return argparse dest names explicitly supplied on the command line."""
    explicit: set[str] = set()
    if argv is None:
        return explicit
    for action in parser._actions:
        if not getattr(action, "option_strings", None):
            continue
        if _arg_was_explicit(argv, action.option_strings):
            explicit.add(action.dest)
    # Also descend into selected subparser to inspect subcommand-specific actions.
    subparsers = [a for a in parser._actions if isinstance(a, argparse._SubParsersAction)]
    if subparsers and argv:
        sp = subparsers[0]
        cmd = None
        for token in argv:
            if token in sp.choices:
                cmd = token
                break
        if cmd:
            for action in sp.choices[cmd]._actions:
                if getattr(action, "option_strings", None) and _arg_was_explicit(argv, action.option_strings):
                    explicit.add(action.dest)
    return explicit


def apply_config(args: argparse.Namespace, parser: argparse.ArgumentParser, argv: Optional[list[str]]) -> argparse.Namespace:
    """Merge optional YAML config into parsed argparse args.

    CLI flags override YAML. YAML can provide sample-specific, command-specific,
    project-level, and defaults values.
    """
    config = load_yaml_config(getattr(args, "config", None))
    if not config:
        return args
    explicit = explicit_cli_dests(parser, argv or [])
    sample_id = getattr(args, "sample_id", None)
    sample = _sample_block(config, sample_id)
    cmd = getattr(args, "cmd", "")

    # Resolve sample_id itself if omitted but YAML has exactly one sample.
    if not sample_id and "sample_id" not in explicit:
        samples = config.get("samples", {})
        if isinstance(samples, Mapping) and len(samples) == 1:
            sid = next(iter(samples.keys()))
            setattr(args, "sample_id", sid)
            sample_id = sid
            sample = _sample_block(config, sid)

    for dest, current in vars(args).copy().items():
        if dest in {"config", "cmd", "func"}:
            continue
        if dest in explicit:
            continue
        val = _lookup_config_value(config, sample, cmd, dest)
        if val is not MISSING:
            setattr(args, dest, val)

    setattr(args, "_raw_config", config)
    setattr(args, "_sample_config", sample)
    return args


def require_args(args: argparse.Namespace, parser: argparse.ArgumentParser, names: Iterable[str]) -> None:
    missing = [n for n in names if getattr(args, n, None) in (None, "")]
    if missing:
        pretty = ", ".join("--" + n.replace("_", "-") for n in missing)
        parser.error(f"missing required arguments after merging YAML config: {pretty}")


def _serialize_value(v: Any) -> Any:
    if isinstance(v, Path):
        return str(v)
    if isinstance(v, set):
        return sorted(v)
    if callable(v):
        return getattr(v, "__name__", str(v))
    try:
        json.dumps(v)
        return v
    except Exception:
        return str(v)


def write_resolved_run_files(args: argparse.Namespace, command_name: str, output_dir: str | Path | None) -> None:
    """Write resolved parameters/config into output_dir when available."""
    if not output_dir:
        return
    out_dir = ensure_dir(output_dir)
    params = {}
    for k, v in sorted(vars(args).items()):
        if k in {"func", "_raw_config", "_sample_config"}:
            continue
        params[k] = _serialize_value(v)
    params["command"] = command_name

    with open(out_dir / "resolved_parameters.tsv", "w") as out:
        out.write("parameter\tvalue\n")
        for k, v in sorted(params.items()):
            out.write(f"{k}\t{v}\n")

    resolved = {
        "command": command_name,
        "parameters": params,
    }
    raw_config = getattr(args, "_raw_config", None)
    if raw_config:
        rc = copy.deepcopy(raw_config)
        rc.pop("_config_path", None)
        resolved["source_config_path"] = raw_config.get("_config_path")
        resolved["source_config"] = rc
    with open(out_dir / "resolved_config.yaml", "w") as out:
        yaml.safe_dump(resolved, out, sort_keys=False)
