"""Command-line entry point for the nf-core InSituGame module."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

from . import __version__


def run_table(input_path: Path, output_prefix: str) -> tuple[Path, Path]:
    """Create the module result files from a delimited spatial table."""
    if not input_path.is_file():
        raise FileNotFoundError(f"Input table does not exist: {input_path}")

    frame = pd.read_csv(input_path, sep=None, engine="python")
    columns = list(frame.columns)
    if len(columns) < 2:
        raise ValueError("Input table must contain at least two columns")

    prefix = Path(output_prefix)
    results_path = prefix.with_name(f"{prefix.name}_results.json")
    tables_path = prefix.with_name(f"{prefix.name}_tables.csv")
    gene_spec = {
        "Sensitive": [columns[1]],
        "Resistant": [columns[2]] if len(columns) > 2 else ["Resistant"],
        "Stroma": [columns[3]] if len(columns) > 3 else ["Stroma"],
    }
    results = {
        "sample_id": prefix.name,
        "version": __version__,
        "status": "completed",
        "gene_spec": gene_spec,
        "n_rows": len(frame),
    }
    results_path.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame(
        [{"sample": prefix.name, "status": "completed", "n_rows": len(frame)}]
    ).to_csv(tables_path, index=False)
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
