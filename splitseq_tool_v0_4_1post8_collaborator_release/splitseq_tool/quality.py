from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass
class QualityAccumulator:
    """Accumulate Phred quality metrics for a fixed/variable sequence region.

    Q30 follows the standard convention of Phred >= 30.  Two complementary
    metrics are retained:
      * percent of individual bases at/above the threshold
      * percent of complete observations for which every base is at/above it
    """

    threshold: int = 30
    phred_offset: int = 33
    observations: int = 0
    bases: int = 0
    bases_passing: int = 0
    observations_all_passing: int = 0
    phred_sum: int = 0
    min_length: int | None = None
    max_length: int | None = None

    def add(self, qual: str) -> None:
        if qual is None or len(qual) == 0:
            return
        scores = [ord(ch) - self.phred_offset for ch in qual]
        self.observations += 1
        self.bases += len(scores)
        n_passing = sum(q >= self.threshold for q in scores)
        self.bases_passing += n_passing
        self.observations_all_passing += int(n_passing == len(scores))
        self.phred_sum += sum(scores)
        qlen = len(scores)
        self.min_length = qlen if self.min_length is None else min(self.min_length, qlen)
        self.max_length = qlen if self.max_length is None else max(self.max_length, qlen)

    @property
    def percent_bases_passing(self) -> float:
        return round(100.0 * self.bases_passing / self.bases, 4) if self.bases else 0.0

    @property
    def percent_observations_all_passing(self) -> float:
        return round(100.0 * self.observations_all_passing / self.observations, 4) if self.observations else 0.0

    @property
    def mean_phred(self) -> float:
        return round(self.phred_sum / self.bases, 4) if self.bases else 0.0

    def row(self, region: str, scope: str, start: int | str = "", expected_length: int | str = "") -> dict[str, object]:
        return {
            "region": region,
            "scope": scope,
            "start_0based": start,
            "expected_length": expected_length,
            "observations_evaluated": self.observations,
            "bases_evaluated": self.bases,
            "bases_Q_ge_threshold": self.bases_passing,
            "percent_bases_Q_ge_threshold": self.percent_bases_passing,
            "observations_all_bases_Q_ge_threshold": self.observations_all_passing,
            "percent_observations_all_bases_Q_ge_threshold": self.percent_observations_all_passing,
            "mean_phred": self.mean_phred,
            "observed_min_length": self.min_length if self.min_length is not None else "",
            "observed_max_length": self.max_length if self.max_length is not None else "",
            "quality_threshold": self.threshold,
            "phred_offset": self.phred_offset,
        }


QUALITY_FIELDS = [
    "region",
    "scope",
    "start_0based",
    "expected_length",
    "observations_evaluated",
    "bases_evaluated",
    "bases_Q_ge_threshold",
    "percent_bases_Q_ge_threshold",
    "observations_all_bases_Q_ge_threshold",
    "percent_observations_all_bases_Q_ge_threshold",
    "mean_phred",
    "observed_min_length",
    "observed_max_length",
    "quality_threshold",
    "phred_offset",
]


def write_quality_qc(path: str | Path, rows: Iterable[dict[str, object]]) -> None:
    with open(path, "w") as out:
        out.write("\t".join(QUALITY_FIELDS) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(field, "")) for field in QUALITY_FIELDS) + "\n")


def summary_rows(prefix: str, acc: QualityAccumulator) -> list[tuple[str, object]]:
    """Return compact key/value metrics suitable for an existing stage summary."""
    return [
        (f"{prefix}_quality_observations", acc.observations),
        (f"{prefix}_Q{acc.threshold}_base_pct", acc.percent_bases_passing),
        (f"{prefix}_all_bases_Q{acc.threshold}_pct", acc.percent_observations_all_passing),
        (f"{prefix}_mean_phred", acc.mean_phred),
    ]
