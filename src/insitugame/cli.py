"""Command-line entry point for the nf-core InSituGame module."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd

from . import __version__, InSituGameAssay


def _numeric_columns(frame: pd.DataFrame) -> list[str]:
    return [
        column
        for column in frame.columns
        if pd.api.types.is_numeric_dtype(frame[column])
    ]


def _coordinates(frame: pd.DataFrame) -> np.ndarray:
    for x_name, y_name in (("x", "y"), ("X", "Y"), ("row", "col")):
        if x_name in frame.columns and y_name in frame.columns:
            return frame[[x_name, y_name]].to_numpy(dtype=float)

    indices = np.arange(len(frame))
    width = max(2, int(np.ceil(np.sqrt(len(frame)))))
    return np.column_stack((indices % width, indices // width)).astype(float)


def _gene_spec(columns: list[str]) -> dict[str, list[str]]:
    aliases = {
        "T_S": ("T_S", "Sensitive"),
        "T_R": ("T_R", "Resistant"),
        "S_CAF": ("S_CAF", "Stroma"),
        "FITNESS": ("FITNESS", "Fitness"),
        "HYPOXIA": ("HYPOXIA", "Hypoxia"),
        "VASCULAR": ("VASCULAR", "Vascular"),
    }
    selected: dict[str, str] = {}
    remaining = [column for column in columns if column not in {"x", "y", "X", "Y", "row", "col"}]
    for key, names in aliases.items():
        match = next((name for name in names if name in columns), None)
        if match is None:
            if not remaining:
                raise ValueError(f"Input table is missing a column for {key}")
            match = remaining[min(len(selected), len(remaining) - 1)]
        selected[key] = match
    return {key: [column] for key, column in selected.items()}


def run_table(input_path: Path, output_prefix: str) -> tuple[Path, Path]:
    """Fit the assay and write its JSON summary and payoff matrix table."""
    if not input_path.is_file():
        raise FileNotFoundError(f"Input table does not exist: {input_path}")

    frame = pd.read_csv(input_path, sep=None, engine="python")
    numeric_columns = _numeric_columns(frame)
    if len(frame) < 3:
        raise ValueError("Input table must contain at least three spatial spots")
    if not numeric_columns:
        raise ValueError("Input table must contain numeric expression columns")

    prefix = Path(output_prefix)
    results_path = prefix.with_name(f"{prefix.name}_results.json")
    tables_path = prefix.with_name(f"{prefix.name}_tables.csv")
    gene_spec = _gene_spec(numeric_columns)
    counts = {column: frame[column].to_numpy(dtype=float) for column in numeric_columns}
    total_umi = np.sum(np.vstack(list(counts.values())), axis=0)
    assay = InSituGameAssay(gene_spec).fit(_coordinates(frame), counts, total_umi)
    results: dict[str, Any] = {
        "sample_id": prefix.name,
        "version": __version__,
        "status": "completed",
        "n_rows": len(frame),
        "gene_spec": gene_spec,
        **assay.summary(),
    }
    results_path.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    payoff = pd.DataFrame(
        assay.payoff_matrix,
        index=("T_S", "T_R", "S_CAF"),
        columns=("T_S", "T_R", "S_CAF"),
    )
    payoff.to_csv(tables_path, index_label="focal_phenotype")
    return results_path, tables_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="insitugame")
    parser.add_argument("--input", type=Path, required=True, help="Delimited input table")
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix for the generated results and tables",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    run_table(args.input, args.output_prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
