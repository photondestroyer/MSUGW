#!/usr/bin/env python3
"""Check integrity and plot a random date from a merged GRIDMET NetCDF file.

The default integrity check reads every time slice of every (time, y, x)
variable in bounded time chunks. It validates the file structure, coordinates,
daily time spacing, fill values, finite data, readable shapes, and raw-data
SHA-256 fingerprints. The plotting step reads only one randomly selected date
and creates one raster panel for each data variable.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any
import hashlib

import numpy as np

try:
    import netCDF4 as nc4
except ImportError as exc:
    raise SystemExit(
        "netCDF4 is required. Install it with: python -m pip install netCDF4"
    ) from exc

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit(
        "matplotlib is required. Install it with: python -m pip install matplotlib"
    ) from exc

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None


DEFAULT_FILE = Path(r"D:\GRIDMET_Merged_Ogallala.nc")
DEFAULT_OUTPUT_DIR = DEFAULT_FILE.with_name("GRIDMET_Merged_Ogallala_plots")
DEFAULT_RANDOM_SEED = 42
DEFAULT_TIME_CHUNK = 8
EXPECTED_TIME_STEP_HOURS = 24.0


def _format_value(value: Any) -> str:
    if isinstance(value, np.generic):
        value = value.item()
    return repr(value)


def _attribute_values(variable: nc4.Variable, attribute_name: str) -> list[Any]:
    if attribute_name not in variable.ncattrs():
        return []
    value = variable.getncattr(attribute_name)
    return np.atleast_1d(np.asarray(value)).tolist()


def _invalid_mask(values: np.ndarray, variable: nc4.Variable) -> np.ndarray:
    invalid = np.zeros(values.shape, dtype=bool)
    for attribute_name in ("_FillValue", "missing_value"):
        for marker in _attribute_values(variable, attribute_name):
            try:
                if np.issubdtype(values.dtype, np.inexact) and np.isnan(marker):
                    invalid |= np.isnan(values)
                else:
                    invalid |= values == marker
            except (TypeError, ValueError):
                continue
    return invalid


def _decode_times(time_variable: nc4.Variable) -> list[datetime]:
    if "units" not in time_variable.ncattrs():
        raise ValueError("time is missing its units attribute")

    calendar = getattr(time_variable, "calendar", "standard")
    decoded = nc4.num2date(
        np.asarray(time_variable[:]),
        units=time_variable.getncattr("units"),
        calendar=calendar,
        only_use_cftime_datetimes=True,
    )
    return [
        datetime(
            value.year,
            value.month,
            value.day,
            value.hour,
            value.minute,
            value.second,
            value.microsecond,
        )
        for value in decoded
    ]


def _validate_coordinate(
    dataset: nc4.Dataset,
    name: str,
    errors: list[str],
    warnings: list[str],
) -> None:
    variable = dataset.variables.get(name)
    if variable is None:
        errors.append(f"Missing coordinate variable: {name}")
        return
    if variable.dimensions != (name,):
        errors.append(
            f"Coordinate {name} has dimensions {variable.dimensions}; "
            f"expected ({name!r},)"
        )
        return

    values = np.asarray(variable[:])
    if values.size == 0:
        errors.append(f"Coordinate {name} is empty")
        return
    if not np.issubdtype(values.dtype, np.number):
        errors.append(f"Coordinate {name} is not numeric: {values.dtype}")
        return
    if not np.all(np.isfinite(values)):
        errors.append(f"Coordinate {name} contains non-finite values")
        return

    differences = np.diff(values)
    if differences.size and not (
        np.all(differences > 0) or np.all(differences < 0)
    ):
        errors.append(f"Coordinate {name} is not strictly monotonic")
    if differences.size > 1 and not np.allclose(
        differences, differences[0], rtol=1e-6, atol=1e-9
    ):
        warnings.append(f"Coordinate {name} spacing is not perfectly regular")


def _validate_time_axis(
    time_variable: nc4.Variable,
    errors: list[str],
    warnings: list[str],
) -> list[datetime]:
    if time_variable.dimensions != ("time",):
        errors.append(
            f"time has dimensions {time_variable.dimensions}; expected ('time',)"
        )
        return []

    raw_times = np.asarray(time_variable[:])
    if raw_times.ndim != 1 or raw_times.size == 0:
        errors.append("time is empty or not one-dimensional")
        return []
    if not np.issubdtype(raw_times.dtype, np.number):
        errors.append(f"time is not numeric: {raw_times.dtype}")
        return []
    if not np.all(np.isfinite(raw_times)):
        errors.append("time contains non-finite values")
        return []
    if np.any(np.diff(raw_times) <= 0):
        errors.append("time is not strictly increasing")

    try:
        dates = _decode_times(time_variable)
    except Exception as exc:
        errors.append(f"Could not decode time: {exc}")
        return []

    if len(dates) > 1:
        decoded_deltas = np.asarray(
            [(later - earlier).total_seconds() / 3600 for earlier, later in zip(dates, dates[1:])],
            dtype=float,
        )
        if not np.allclose(decoded_deltas, EXPECTED_TIME_STEP_HOURS):
            errors.append(
                "time is not a daily series; unique hour steps are "
                f"{np.unique(decoded_deltas)}"
            )

    clock_times = {(date.hour, date.minute, date.second, date.microsecond) for date in dates}
    if len(clock_times) > 1:
        warnings.append(f"time contains more than one clock time: {sorted(clock_times)}")

    return dates


def _data_variables(dataset: nc4.Dataset) -> list[nc4.Variable]:
    return [
        variable
        for variable in dataset.variables.values()
        if variable.dimensions == ("time", "y", "x")
    ]


def _scan_variable(
    variable: nc4.Variable,
    time_length: int,
    y_length: int,
    x_length: int,
    time_chunk: int,
    errors: list[str],
) -> dict[str, Any]:
    variable.set_auto_mask(False)
    variable.set_auto_scale(False)
    digest = hashlib.sha256()
    total_cells = time_length * y_length * x_length
    valid_count = 0
    invalid_count = 0
    nonfinite_valid_count = 0
    all_fill_slices = 0
    minimum: float | None = None
    maximum: float | None = None

    starts = range(0, time_length, time_chunk)
    iterator = (
        tqdm(starts, total=(time_length + time_chunk - 1) // time_chunk, desc=f"Scanning {variable.name}", unit="chunk")
        if tqdm is not None
        else starts
    )

    for start in iterator:
        stop = min(start + time_chunk, time_length)
        try:
            raw_chunk = variable[start:stop, :, :]
        except Exception as exc:
            errors.append(f"{variable.name}: read failed for time {start}:{stop}: {exc}")
            continue

        returned_mask = np.ma.getmaskarray(raw_chunk) if np.ma.isMaskedArray(raw_chunk) else None
        values = np.asarray(raw_chunk.data if np.ma.isMaskedArray(raw_chunk) else raw_chunk)
        expected_shape = (stop - start, y_length, x_length)
        if values.shape != expected_shape:
            errors.append(
                f"{variable.name}: read shape {values.shape}; expected {expected_shape}"
            )
            continue

        values = np.ascontiguousarray(values)
        digest.update(values.tobytes(order="C"))
        invalid = _invalid_mask(values, variable)
        if returned_mask is not None:
            invalid |= returned_mask

        for offset in range(stop - start):
            frame = values[offset]
            frame_invalid = invalid[offset]
            valid = frame[~frame_invalid]
            invalid_count += int(np.count_nonzero(frame_invalid))
            valid_count += int(valid.size)
            if valid.size == 0:
                all_fill_slices += 1
                continue

            finite = np.isfinite(valid)
            nonfinite_here = int(np.count_nonzero(~finite))
            nonfinite_valid_count += nonfinite_here
            finite_values = valid[finite]
            if finite_values.size:
                local_min = float(np.min(finite_values))
                local_max = float(np.max(finite_values))
                minimum = local_min if minimum is None else min(minimum, local_min)
                maximum = local_max if maximum is None else max(maximum, local_max)

    if valid_count + invalid_count != total_cells:
        errors.append(
            f"{variable.name}: scanned {valid_count + invalid_count:,} cells; "
            f"expected {total_cells:,}"
        )
    if all_fill_slices:
        errors.append(
            f"{variable.name}: {all_fill_slices} time slice(s) contain only fill values"
        )
    if nonfinite_valid_count:
        errors.append(
            f"{variable.name}: {nonfinite_valid_count:,} non-finite valid values found"
        )

    return {
        "name": variable.name,
        "dtype": str(variable.dtype),
        "fill_value": _format_value(
            _attribute_values(variable, "_FillValue")[0]
            if _attribute_values(variable, "_FillValue")
            else None
        ),
        "cells": total_cells,
        "valid_values": valid_count,
        "fill_values": invalid_count,
        "all_fill_slices": all_fill_slices,
        "nonfinite_valid_values": nonfinite_valid_count,
        "min": minimum,
        "max": maximum,
        "raw_data_sha256": digest.hexdigest(),
    }


def validate_file(
    path: Path,
    time_chunk: int = DEFAULT_TIME_CHUNK,
    full_scan: bool = True,
) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"NetCDF file does not exist: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"NetCDF file is empty: {path}")
    if time_chunk < 1:
        raise ValueError("time_chunk must be at least 1")

    errors: list[str] = []
    warnings: list[str] = []
    reports: list[dict[str, Any]] = []
    dates: list[datetime] = []

    try:
        with nc4.Dataset(str(path), "r") as dataset:
            dataset.set_auto_mask(False)
            dataset.set_auto_scale(False)
            required_dimensions = {"time", "y", "x"}
            missing_dimensions = required_dimensions - set(dataset.dimensions)
            if missing_dimensions:
                errors.append(f"Missing dimensions: {sorted(missing_dimensions)}")
                return {
                    "file": str(path),
                    "size_bytes": path.stat().st_size,
                    "full_scan": full_scan,
                    "errors": errors,
                    "warnings": warnings,
                    "variables": reports,
                    "dates": dates,
                }

            time_length = len(dataset.dimensions["time"])
            y_length = len(dataset.dimensions["y"])
            x_length = len(dataset.dimensions["x"])
            if not time_length or not y_length or not x_length:
                errors.append(
                    f"Empty dimension: time={time_length}, y={y_length}, x={x_length}"
                )

            time_variable = dataset.variables.get("time")
            if time_variable is None:
                errors.append("Missing time variable")
            else:
                dates = _validate_time_axis(time_variable, errors, warnings)

            _validate_coordinate(dataset, "x", errors, warnings)
            _validate_coordinate(dataset, "y", errors, warnings)

            spatial_ref = dataset.variables.get("spatial_ref")
            if spatial_ref is None:
                errors.append("Missing spatial_ref variable")
            elif spatial_ref.dimensions:
                errors.append(
                    f"spatial_ref has dimensions {spatial_ref.dimensions}; expected scalar"
                )
            elif not np.isfinite(np.asarray(spatial_ref[...], dtype=float)).all():
                errors.append("spatial_ref contains non-finite values")

            data_variables = _data_variables(dataset)
            if not data_variables:
                errors.append("No variables with dimensions (time, y, x) found")

            malformed_time_variables = sorted(
                name
                for name, variable in dataset.variables.items()
                if "time" in variable.dimensions
                and variable.dimensions != ("time", "y", "x")
                and name != "time"
            )
            if malformed_time_variables:
                errors.append(
                    "Variables with unexpected time dimensions: "
                    f"{malformed_time_variables}"
                )

            for variable in data_variables:
                if variable.shape != (time_length, y_length, x_length):
                    errors.append(
                        f"{variable.name}: shape {variable.shape}; expected "
                        f"{(time_length, y_length, x_length)}"
                    )
                if not np.issubdtype(variable.dtype, np.number):
                    errors.append(f"{variable.name}: non-numeric dtype {variable.dtype}")
                if "_FillValue" not in variable.ncattrs():
                    warnings.append(f"{variable.name}: no _FillValue attribute")
                if getattr(variable, "grid_mapping", None) != "spatial_ref":
                    warnings.append(
                        f"{variable.name}: grid_mapping does not reference spatial_ref"
                    )

                if full_scan:
                    reports.append(
                        _scan_variable(
                            variable,
                            time_length,
                            y_length,
                            x_length,
                            time_chunk,
                            errors,
                        )
                    )
                else:
                    reports.append(
                        {
                            "name": variable.name,
                            "dtype": str(variable.dtype),
                            "shape": list(variable.shape),
                            "scanned": False,
                        }
                    )

    except Exception as exc:
        errors.append(f"Could not open or validate NetCDF file: {exc}")

    return {
        "file": str(path),
        "size_bytes": path.stat().st_size,
        "full_scan": full_scan,
        "dimensions": {"time": time_length, "y": y_length, "x": x_length}
        if "time_length" in locals()
        else {},
        "date_start": dates[0].isoformat() if dates else None,
        "date_end": dates[-1].isoformat() if dates else None,
        "variables": reports,
        "warnings": warnings,
        "errors": errors,
    }


def _to_float_frame(values: Any) -> np.ndarray:
    masked = np.ma.asarray(values, dtype=np.float64)
    return np.asarray(masked.filled(np.nan), dtype=np.float64)


def _oriented_frame(
    frame: np.ndarray,
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    display = frame
    x_left, x_right = float(x_coordinates[0]), float(x_coordinates[-1])
    y_bottom, y_top = float(y_coordinates[-1]), float(y_coordinates[0])
    if x_left > x_right:
        display = display[:, ::-1]
        x_left, x_right = x_right, x_left
    if y_bottom > y_top:
        display = display[::-1, :]
        y_bottom, y_top = y_top, y_bottom
    return display, (x_left, x_right, y_bottom, y_top)


def _color_limits(frame: np.ndarray) -> tuple[float, float]:
    finite = frame[np.isfinite(frame)]
    if not finite.size:
        raise ValueError("Random date contains no finite values")
    lower, upper = np.nanpercentile(finite, [2, 98])
    if lower == upper:
        padding = max(abs(float(lower)) * 0.05, 1.0)
        lower -= padding
        upper += padding
    return float(lower), float(upper)


def plot_random_date(
    path: Path,
    output_dir: Path,
    seed: int = DEFAULT_RANDOM_SEED,
) -> tuple[Path, int, datetime, list[str]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    with nc4.Dataset(str(path), "r") as dataset:
        dataset.set_auto_maskandscale(True)
        time_variable = dataset.variables["time"]
        dates = _decode_times(time_variable)
        random_index = int(rng.integers(0, len(dates)))
        random_date = dates[random_index]
        x_coordinates = np.asarray(dataset.variables["x"][:])
        y_coordinates = np.asarray(dataset.variables["y"][:])
        data_variables = _data_variables(dataset)
        frames = [
            (variable.name, _to_float_frame(variable[random_index, :, :]))
            for variable in data_variables
        ]

    column_count = 4
    row_count = (len(frames) + column_count - 1) // column_count
    figure, axes = plt.subplots(
        row_count,
        column_count,
        figsize=(20, max(5, row_count * 4.5)),
        squeeze=False,
        constrained_layout=True,
    )
    figure.suptitle(
        f"GRIDMET variables on {random_date:%Y-%m-%d %H:%M} "
        f"(time index {random_index}, seed={seed})",
        fontsize=16,
    )

    for axis, (name, frame) in zip(axes.flat, frames):
        display, extent = _oriented_frame(frame, x_coordinates, y_coordinates)
        vmin, vmax = _color_limits(frame)
        image = axis.imshow(
            display,
            extent=extent,
            origin="lower",
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
            aspect="auto",
        )
        axis.set_title(name)
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04, label=name)

    for axis in axes.flat[len(frames) :]:
        axis.set_visible(False)

    figure_path = output_dir / (
        f"{path.stem}_random_{random_date:%Y%m%dT%H%M}_seed{seed}.png"
    )
    figure.savefig(figure_path, dpi=160, bbox_inches="tight")
    plt.close(figure)
    return figure_path, random_index, random_date, [name for name, _ in frames]


def _write_report(output_dir: Path, path: Path, report: dict[str, Any]) -> Path:
    report_path = output_dir / f"{path.stem}_integrity.json"
    with report_path.open("w", encoding="utf-8") as output:
        json.dump(report, output, indent=2, default=str)
        output.write("\n")
    return report_path


def _print_report(report: dict[str, Any]) -> None:
    print(f"File: {report['file']}")
    print(f"Size: {report['size_bytes'] / 1024**3:.2f} GiB")
    print(f"Dimensions: {report.get('dimensions', {})}")
    print(f"Date range: {report.get('date_start')} through {report.get('date_end')}")
    print(f"Full scan: {report['full_scan']}")
    for variable in report["variables"]:
        if variable.get("scanned") is False:
            print(f"  {variable['name']}: metadata checked only")
            continue
        print(
            f"  {variable['name']}: valid={variable['valid_values']:,}, "
            f"fill={variable['fill_values']:,}, "
            f"range={variable['min']}..{variable['max']}, "
            f"raw_sha256={variable['raw_data_sha256']}"
        )
    if report["warnings"]:
        print("Warnings:")
        for warning in report["warnings"]:
            print(f"  - {warning}")
    if report["errors"]:
        print("Errors:")
        for error in report["errors"]:
            print(f"  - {error}")
        print("RESULT: FAIL")
    else:
        print("RESULT: PASS")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check a merged GRIDMET NetCDF and plot every variable at a random date."
    )
    parser.add_argument(
        "path",
        nargs="?",
        type=Path,
        default=DEFAULT_FILE,
        help=f"NetCDF file (default: {DEFAULT_FILE})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_RANDOM_SEED,
        help="Random seed for date selection (default: 42)",
    )
    parser.add_argument(
        "--time-chunk",
        type=int,
        default=DEFAULT_TIME_CHUNK,
        help="Number of daily slices read per integrity chunk (default: 8)",
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Skip the full data read and run structural checks only",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Run integrity checks without creating the random-date plot",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    report = validate_file(
        args.path,
        time_chunk=args.time_chunk,
        full_scan=not args.metadata_only,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = _write_report(args.output_dir, args.path, report)
    _print_report(report)
    print(f"Integrity report: {report_path}")

    if report["errors"] or args.skip_plots:
        if args.skip_plots:
            print("Plots skipped by request")
        return 1 if report["errors"] else 0

    figure_path, random_index, random_date, variable_names = plot_random_date(
        args.path,
        args.output_dir,
        seed=args.seed,
    )
    print(f"Plotted variables: {', '.join(variable_names)}")
    print(f"Random date: {random_date:%Y-%m-%d %H:%M} (time index {random_index})")
    print(f"Plot: {figure_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())