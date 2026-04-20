#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare workflow outputs against a BiomiX gold standard.")
    parser.add_argument("--actual-root", required=True)
    parser.add_argument("--gold-root", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--atol", type=float, default=1e-8)
    parser.add_argument("--rtol", type=float, default=1e-8)
    return parser.parse_args()


def compare_frames(actual: pd.DataFrame, gold: pd.DataFrame, atol: float, rtol: float) -> list[str]:
    errors: list[str] = []

    if list(actual.columns) != list(gold.columns):
        errors.append("Column names differ.")
        return errors

    if actual.shape != gold.shape:
        errors.append(f"Shape differs: actual={actual.shape}, gold={gold.shape}.")
        return errors

    for column in actual.columns:
        actual_series = actual[column]
        gold_series = gold[column]

        if pd.api.types.is_numeric_dtype(actual_series) and pd.api.types.is_numeric_dtype(gold_series):
            actual_values = actual_series.to_numpy(dtype=float)
            gold_values = gold_series.to_numpy(dtype=float)
            if not np.allclose(actual_values, gold_values, equal_nan=True, atol=atol, rtol=rtol):
                errors.append(f"Numeric values differ in column '{column}'.")
        else:
            actual_values = actual_series.fillna("<NA>").astype(str)
            gold_values = gold_series.fillna("<NA>").astype(str)
            if not actual_values.equals(gold_values):
                errors.append(f"String values differ in column '{column}'.")

    return errors


def main() -> None:
    args = parse_args()

    actual_root = Path(args.actual_root).resolve()
    gold_root = Path(args.gold_root).resolve()
    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))

    results = []

    for entry in manifest["files"]:
        relative_path = Path(entry["relative_path"])
        actual_path = actual_root / relative_path
        gold_path = gold_root / relative_path

        file_result = {
            "relative_path": str(relative_path),
            "actual_exists": actual_path.exists(),
            "gold_exists": gold_path.exists(),
            "passed": False,
            "errors": [],
        }

        if not actual_path.exists():
            file_result["errors"].append("Actual file is missing.")
        if not gold_path.exists():
            file_result["errors"].append("Gold file is missing.")

        if file_result["errors"]:
            results.append(file_result)
            continue

        actual_frame = pd.read_csv(actual_path, sep="\t")
        gold_frame = pd.read_csv(gold_path, sep="\t")
        file_result["errors"] = compare_frames(actual_frame, gold_frame, args.atol, args.rtol)
        file_result["passed"] = not file_result["errors"]
        results.append(file_result)

    failed = [result for result in results if not result["passed"]]
    report = {
        "passed": not failed,
        "files": results,
    }

    Path(args.report).write_text(json.dumps(report, indent=2), encoding="utf-8")

    if failed:
        raise SystemExit("Gold-standard comparison failed.")


if __name__ == "__main__":
    main()
