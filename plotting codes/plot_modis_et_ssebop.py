#!/usr/bin/env python3
"""Plot all time-dependent raster variables in the MODIS ET NetCDF file.

For every variable with dimensions (time, y, x), the script creates one row
containing a temporal-mean map for the complete dataset, a map for one randomly
selected date, and the spatial-mean time series with its spatial min/max band.
The raster data is read one time slice at a time so the complete 699-step file
does not need to be loaded into memory.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
from typing import Any

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
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit(
        "matplotlib is required. Install it with: python -m pip install matplotlib"
    ) from exc


DEFAULT_FILE = Path(__file__).resolve().with_name(
    "MODIS_ET_SSEBop_Merged_Ogallala.nc"
)
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().with_name("modis_et_ssebop_plots")
DEFAULT_SEED = 42
DEFAULT_TIME_CHUNK = 1


def _decode_dates(time_variable: nc4.Variable) -> list[datetime]:
    if "units" not in time_variable.ncattrs():
        raise ValueError("The time variable is missing its units attribute")

    calendar = getattr(time_variable, "calendar", "standard")
    decoded = nc4.num2date(
        np.asarray(time_variable[:]),
        units=time_variable.getncattr("units"),
        calendar=calendar,
        only_use_cftime_datetimes=True,
    )
    return [datetime(value.year, value.month, value.day) for value in decoded]


def _as_float_frame(values: Any) -> np.ndarray:
    """Convert a masked NetCDF slice to float values with missing pixels as NaN."""
    masked = np.ma.asarray(values, dtype=np.float64)
    return np.asarray(masked.filled(np.nan), dtype=np.float64)


def _display_orientation(
    frame: np.ndarray, x_coordinates: np.ndarray, y_coordinates: np.ndarray
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    display = frame
    x_left = float(x_coordinates[0])
    x_right = float(x_coordinates[-1])
    y_bottom = float(y_coordinates[-1])
    y_top = float(y_coordinates[0])

    if x_left > x_right:
        display = display[:, ::-1]
        x_left, x_right = x_right, x_left
    if y_bottom > y_top:
        display = display[::-1, :]
        y_bottom, y_top = y_top, y_bottom

    return display, (x_left, x_right, y_bottom, y_top)


def _plot_limits(*frames: np.ndarray) -> tuple[float, float]:
    finite_values = [frame[np.isfinite(frame)] for frame in frames]
    finite_values = [values for values in finite_values if values.size]
    if not finite_values:
        raise ValueError("No finite values are available for color scaling")

    combined = np.concatenate(finite_values)
    lower, upper = np.nanpercentile(combined, [2, 98])
    if not np.isfinite(lower) or not np.isfinite(upper):
        lower, upper = np.nanmin(combined), np.nanmax(combined)
    if lower == upper:
        padding = max(abs(float(lower)) * 0.05, 1.0)
        lower -= padding
        upper += padding
    return float(lower), float(upper)


def _scan_variable(
    variable: nc4.Variable,
    dates: list[datetime],
    random_index: int,
    y_length: int,
    x_length: int,
    time_chunk: int,
) -> dict[str, Any]:
    time_means = np.full(len(dates), np.nan, dtype=np.float64)
    time_mins = np.full(len(dates), np.nan, dtype=np.float64)
    time_maxs = np.full(len(dates), np.nan, dtype=np.float64)
    temporal_sum = np.zeros((y_length, x_length), dtype=np.float64)
    temporal_count = np.zeros((y_length, x_length), dtype=np.uint32)
    random_frame: np.ndarray | None = None

    variable.set_auto_maskandscale(True)
    for start in range(0, len(dates), time_chunk):
        stop = min(start + time_chunk, len(dates))
        chunk = variable[start:stop, :, :]
        if np.ma.isMaskedArray(chunk):
            chunk_values = chunk
        else:
            chunk_values = np.ma.asarray(chunk)

        for offset in range(stop - start):
            time_index = start + offset
            frame = _as_float_frame(chunk_values[offset])
            finite = np.isfinite(frame)
            if not np.any(finite):
                continue

            valid_values = frame[finite]
            time_means[time_index] = float(np.mean(valid_values))
            time_mins[time_index] = float(np.min(valid_values))
            time_maxs[time_index] = float(np.max(valid_values))
            temporal_sum[finite] += frame[finite]
            temporal_count[finite] += 1

            if time_index == random_index:
                random_frame = frame.copy()

    if random_frame is None:
        raise ValueError(
            f"Variable {variable.name} has no finite values at random index "
            f"{random_index}"
        )

    temporal_mean = np.full((y_length, x_length), np.nan, dtype=np.float64)
    valid_pixels = temporal_count > 0
    temporal_mean[valid_pixels] = (
        temporal_sum[valid_pixels] / temporal_count[valid_pixels]
    )

    if not np.any(np.isfinite(temporal_mean)):
        raise ValueError(f"Variable {variable.name} contains no finite data")

    return {
        "name": variable.name,
        "temporal_mean": temporal_mean,
        "random_frame": random_frame,
        "time_means": time_means,
        "time_mins": time_mins,
        "time_maxs": time_maxs,
        "valid_pixel_count": int(np.count_nonzero(valid_pixels)),
    }


def _plot_raster(
    axis: Any,
    frame: np.ndarray,
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
    vmin: float,
    vmax: float,
    title: str,
) -> Any:
    display, extent = _display_orientation(frame, x_coordinates, y_coordinates)
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
    axis.set_title(title)
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    return image


def _plot_time_series(
    axis: Any,
    dates: list[datetime],
    statistics: dict[str, Any],
    title: str,
) -> None:
    axis.plot(dates, statistics["time_means"], color="tab:blue", linewidth=1.2)
    axis.fill_between(
        dates,
        statistics["time_mins"],
        statistics["time_maxs"],
        color="tab:blue",
        alpha=0.18,
        linewidth=0,
        label="spatial min-max",
    )
    axis.set_title(title)
    axis.set_xlabel("Date")
    axis.set_ylabel("Spatial mean")
    axis.grid(True, alpha=0.25)
    axis.legend(loc="best")
    axis.xaxis.set_major_locator(mdates.AutoDateLocator())
    axis.xaxis.set_major_formatter(mdates.ConciseDateFormatter(axis.xaxis.get_major_locator()))


def create_plots(
    path: Path,
    output_dir: Path,
    seed: int | None = DEFAULT_SEED,
    time_chunk: int = DEFAULT_TIME_CHUNK,
) -> tuple[Path, int, datetime, list[str]]:
    """Create the complete-dataset and random-date plots."""
    if not path.is_file():
        raise FileNotFoundError(f"NetCDF file does not exist: {path}")
    if time_chunk < 1:
        raise ValueError("time_chunk must be at least 1")

    rng = np.random.default_rng(seed)
    output_dir.mkdir(parents=True, exist_ok=True)

    with nc4.Dataset(str(path), "r") as dataset:
        dataset.set_auto_maskandscale(True)
        required_dimensions = {"time", "y", "x"}
        missing_dimensions = required_dimensions - set(dataset.dimensions)
        if missing_dimensions:
            raise ValueError(f"Missing required dimensions: {sorted(missing_dimensions)}")

        time_variable = dataset.variables.get("time")
        if time_variable is None:
            raise ValueError("The dataset has no time variable")
        dates = _decode_dates(time_variable)
        if not dates:
            raise ValueError("The time coordinate is empty")

        y_length = len(dataset.dimensions["y"])
        x_length = len(dataset.dimensions["x"])
        x_coordinates = np.asarray(dataset.variables["x"][:])
        y_coordinates = np.asarray(dataset.variables["y"][:])
        random_index = int(rng.integers(0, len(dates)))
        random_date = dates[random_index]

        data_variables = [
            variable
            for variable in dataset.variables.values()
            if variable.dimensions == ("time", "y", "x")
        ]
        if not data_variables:
            raise ValueError("No variables with dimensions (time, y, x) were found")

        statistics = [
            _scan_variable(
                variable,
                dates,
                random_index,
                y_length,
                x_length,
                time_chunk,
            )
            for variable in data_variables
        ]

    row_count = len(statistics)
    figure, axes = plt.subplots(
        row_count,
        3,
        figsize=(18, max(5.5, 4.8 * row_count)),
        squeeze=False,
        constrained_layout=True,
    )
    figure.suptitle(
        f"{path.name}\nRandom date: {random_date:%Y-%m-%d} "
        f"(index {random_index}, seed={seed})",
        fontsize=14,
    )

    for row, variable_statistics in enumerate(statistics):
        name = variable_statistics["name"]
        vmin, vmax = _plot_limits(
            variable_statistics["temporal_mean"],
            variable_statistics["random_frame"],
        )
        mean_image = _plot_raster(
            axes[row, 0],
            variable_statistics["temporal_mean"],
            x_coordinates,
            y_coordinates,
            vmin,
            vmax,
            f"{name}: mean over all {len(dates)} dates",
        )
        random_image = _plot_raster(
            axes[row, 1],
            variable_statistics["random_frame"],
            x_coordinates,
            y_coordinates,
            vmin,
            vmax,
            f"{name}: {random_date:%Y-%m-%d}",
        )
        _plot_time_series(
            axes[row, 2],
            dates,
            variable_statistics,
            f"{name}: complete time series",
        )
        figure.colorbar(
            mean_image,
            ax=[axes[row, 0], axes[row, 1]],
            orientation="vertical",
            fraction=0.046,
            pad=0.04,
            label=name,
        )
        random_image.set_clim(vmin, vmax)

    figure_path = output_dir / f"{path.stem}_all_variables.png"
    figure.savefig(figure_path, dpi=160, bbox_inches="tight")
    plt.close(figure)

    return figure_path, random_index, random_date, [item["name"] for item in statistics]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot all MODIS ET raster variables for the full period and a random date."
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
        help=f"Directory for PNG output (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help="Random seed used to select the date (default: 42)",
    )
    parser.add_argument(
        "--time-chunk",
        type=int,
        default=DEFAULT_TIME_CHUNK,
        help="Number of time slices read at once (default: 1)",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    figure_path, random_index, random_date, variable_names = create_plots(
        args.path,
        args.output_dir,
        seed=args.seed,
        time_chunk=args.time_chunk,
    )
    print(f"Plotted variables: {', '.join(variable_names)}")
    print(f"Random date: {random_date:%Y-%m-%d} (time index {random_index})")
    print(f"Saved figure: {figure_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())