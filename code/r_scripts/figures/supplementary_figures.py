#!/usr/bin/env python3
"""Plot the two NECOS seeds after the matching *_lmm.R script has completed.

CONTENTS
  1. Imports and editable plotting settings
  2. Resolve paths and optional command-line arguments
  3. Load the cortical surfaces and vertex-to-parcel mapping
  4. Load the network and parcelwise model estimates
  5. Repeated surface-drawing and figure-saving helpers
  6. Draw the seven-network panels
  7. Draw the parcelwise effect and FDR panels
  8. Assemble and save the supplementary figure

Run the numbered sections in order, or use the # %% IDE cells. Edit paths and
plot settings in section 1; command-line options still override those values.
The plotting steps are visible below and execute directly, without a main().

Creates Yeo-7 network panels, unthresholded parcelwise panels, optional parcelwise
FDR maps, and the complete Supplementaryfigure_posthoc_seedbased.png/.svg.

Example (repository root):
    python code/r_scripts/figures/supplementary_figures.py

Default results: <project>/results/posthoc_seedbased. Surface geometry defaults
to the installed BrainSpace datasets. Alternatively use --brainspace-data
/path/to/brainspace/datasets, or supply BOTH --surface-dir and --parcellation.
Required files: surfaces/conte69_32k_{lh,rh}.gii and
parcellations/schaefer_400_conte69.csv (left vertices followed by right vertices).

All panels use the same symmetric color range, configurable with --vmax, so the
assembled figure's single color bar is valid for both network and parcel maps.
Unthresholded maps show coefficients, not significance. Gray marks missing or
excluded values; FDR maps additionally mask nonsignificant targets.
Dependencies: numpy, pandas, matplotlib, nibabel; BrainSpace is optional when
the surface files are provided explicitly.
"""

# %% 1. Imports and editable plotting settings
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Safe for batch runs on machines without a display.
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

SEED_NAMES = {377: "Right vPFC", 100: "Left insula / frontal operculum"}
VIEWS = (
    ("lh", 180, "Left lateral"),
    ("lh", 0, "Left medial"),
    ("rh", 180, "Right medial"),
    ("rh", 0, "Right lateral"),
)

# None uses repository-relative results or the installed BrainSpace datasets.
# Example for external geometry: BRAINSPACE_DATA = Path("/path/to/brainspace/datasets").
PROJECT_ROOT = None
SEEDBASED_ROOT = None
LABEL_FILE = None
OUTPUT_ROOT = None
BRAINSPACE_DATA = None
SURFACE_DIR = None
PARCELLATION_FILE = None
SEEDS_TO_PLOT = [377, 100]
MAP_TYPES = ["interaction_estimates", "fdr_significant"]
VMAX = 0.25  # Same symmetric range for network and parcel maps.
DPI = 300


# %% 2. Resolve paths and optional command-line arguments
# Locate the repository so the plotting script works from the project root or
# from an IDE cell.
def find_project_root() -> Path:
    """Locate the project for a full-file run or an IDE cell session."""
    start = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
    for candidate in (start, *start.parents):
        if (candidate / "code").is_dir() and (candidate / "data").is_dir():
            return candidate
    # A standalone copy can still receive --project-root on the command line.
    return Path.cwd()


# Read plotting options and resolve all input/output locations in one place.
def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT or find_project_root())
    parser.add_argument("--seedbased-root", type=Path, default=SEEDBASED_ROOT)
    parser.add_argument("--label-file", type=Path, default=LABEL_FILE)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=OUTPUT_ROOT,
        help="Default: seedbased-root; preserves original figure subfolders",
    )
    parser.add_argument(
        "--brainspace-data",
        type=Path,
        default=BRAINSPACE_DATA,
        help="BrainSpace datasets directory",
    )
    parser.add_argument(
        "--surface-dir",
        type=Path,
        default=SURFACE_DIR,
        help="Directory with conte69_32k_lh.gii and conte69_32k_rh.gii",
    )
    parser.add_argument(
        "--parcellation",
        type=Path,
        default=PARCELLATION_FILE,
        help="Schaefer-400 vertex labels, left then right hemisphere",
    )
    parser.add_argument("--seeds", type=int, nargs="+", choices=(377, 100), default=SEEDS_TO_PLOT)
    parser.add_argument(
        "--map-types",
        nargs="+",
        choices=("interaction_estimates", "fdr_significant"),
        default=MAP_TYPES,
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=VMAX,
        help="Common absolute color limit for ALL panels (default: 0.25)",
    )
    parser.add_argument("--dpi", type=int, default=DPI)
    args = parser.parse_args(argv)
    if not np.isfinite(args.vmax) or args.vmax <= 0 or args.dpi < 1:
        parser.error("--vmax and --dpi must be positive")
    args.seeds = list(dict.fromkeys(args.seeds))
    root = args.project_root.resolve()
    args.seedbased_root = args.seedbased_root or root / "results" / "posthoc_seedbased"
    args.output_root = args.output_root or args.seedbased_root
    args.label_file = (
        args.label_file or root / "data" / "yeo_networks_17_reordered_400_complete_labels.csv"
    )
    if args.surface_dir is None or args.parcellation is None:
        datasets = args.brainspace_data
        if datasets is None:
            spec = importlib.util.find_spec("brainspace")
            if spec is None or spec.origin is None:
                parser.error(
                    "Supply --brainspace-data or both --surface-dir and --parcellation, or install brainspace."
                )
            datasets = Path(spec.origin).parent / "datasets"
        args.surface_dir = args.surface_dir or datasets / "surfaces"
        args.parcellation = (
            args.parcellation or datasets / "parcellations" / "schaefer_400_conte69.csv"
        )
    return args


args = parse_args([] if "get_ipython" in globals() or "__file__" not in globals() else None)


# %% 3. Load the cortical surfaces and vertex-to-parcel mapping
# Read Conte69 geometry and the vertex-to-Schaefer parcel assignment. The zero
# label is the medial wall; parcel labels 1..400 are the cortical data values.
surfaces = {}
for hemisphere in ("lh", "rh"):
    image = nib.load(args.surface_dir / f"conte69_32k_{hemisphere}.gii")
    surfaces[hemisphere] = (image.darrays[0].data, image.darrays[1].data)
labels = np.loadtxt(args.parcellation)
if (
    labels.ndim != 1
    or not np.isfinite(labels).all()
    or not np.equal(labels, labels.astype(int)).all()
):
    raise ValueError("Surface parcellation must be a 1D array of integer parcel labels.")
labels = labels.astype(int)
if len(labels) != sum(len(surface[0]) for surface in surfaces.values()):
    raise ValueError("Parcellation vertex count does not match the two surfaces.")
if labels.min() < 0 or labels.max() > 400 or not set(range(1, 401)).issubset(set(labels)):
    raise ValueError(
        "Surface parcellation must contain all parcel IDs 1..400; 0 is the medial wall."
    )

vertex_labels = labels


# %% 4. Load the network and parcelwise model estimates
# Read the R interaction tables and expand the seven network coefficients to a
# 400-parcel vector so network and parcel maps use the same surface renderer.
root = args.seedbased_root / "LMM"
networks = pd.read_csv(root / "Yeo7" / "all_seed_yeo7_network_interactions.tsv", sep="\t")
labels = pd.read_csv(args.label_file)
if sorted(labels["7networks_order"].tolist()) != list(range(1, 401)):
    raise ValueError("The label table must contain each parcel 1..400 once.")
values = {}
for seed in args.seeds:
    subset = networks.loc[networks.seed_parcel == seed]
    if (
        len(subset) != 7
        or subset.network.duplicated().any()
        or not np.isfinite(subset.estimate).all()
    ):
        raise ValueError(f"Expected seven finite network estimates for seed {seed}.")
    network_effects = subset.set_index("network").estimate
    mapped = labels["7networks"].map(network_effects)
    if mapped.isna().any():
        raise ValueError("A parcel network has no fitted coefficient.")
    network_values = np.full(401, np.nan)
    network_values[labels["7networks_order"].to_numpy(int)] = mapped
    parcel_results = pd.read_csv(root / f"seed_{seed:03d}_interaction_results.tsv", sep="\t")
    if sorted(parcel_results.target_parcel.tolist()) != [p for p in range(1, 401) if p != seed]:
        raise ValueError(f"Expected exactly 399 distinct non-self targets for seed {seed}.")
    if not np.isfinite(parcel_results[["estimate", "p_adjusted_fdr_bh"]]).all().all():
        raise ValueError("Parcel estimates or FDR values are nonfinite.")
    parcel_values = np.full(401, np.nan)
    parcel_values[parcel_results.target_parcel.to_numpy(int)] = parcel_results.estimate
    significant = np.full(401, np.nan)
    selected = parcel_results.loc[parcel_results.p_adjusted_fdr_bh < 0.05]
    significant[selected.target_parcel.to_numpy(int)] = selected.estimate
    values[seed] = {
        "network": network_values,
        "interaction_estimates": parcel_values,
        "fdr_significant": significant,
    }


# %% 5. Repeated surface-drawing and figure-saving helpers
# The actual figures are created in sections 6-8.
# Draw one hemisphere using the mean coefficient of each triangular face.
def draw_surface(
    ax,
    coordinates: np.ndarray,
    faces: np.ndarray,
    values: np.ndarray,
    azimuth: int,
    title: str,
    vmax: float,
) -> None:
    # Match the original renderer: color each mesh face by its vertex mean.
    face_values = np.ma.masked_invalid(values[faces].mean(axis=1))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#d9d9d9")
    surface = ax.plot_trisurf(
        coordinates[:, 0],
        coordinates[:, 1],
        coordinates[:, 2],
        triangles=faces,
        cmap=cmap,
        linewidth=0,
        antialiased=False,
        shade=False,
    )
    surface.set_array(face_values)
    surface.set_clim(-vmax, vmax)
    surface.set_rasterized(True)
    ax.view_init(elev=0, azim=azimuth)
    ax.set_box_aspect(np.ptp(coordinates, axis=0))
    ax.set_axis_off()
    ax.set_title(title, fontsize=10, pad=1)


# Put the four requested views of one map into a single figure row.
def draw_row(
    fig,
    gridspec,
    row: int,
    values: np.ndarray,
    vertex_labels: np.ndarray,
    surfaces: dict,
    vmax: float,
) -> None:
    vertex_values = values[vertex_labels]
    left_count = len(surfaces["lh"][0])
    hemisphere_values = {"lh": vertex_values[:left_count], "rh": vertex_values[left_count:]}
    for column, (hemisphere, azimuth, title) in enumerate(VIEWS):
        ax = fig.add_subplot(gridspec[row, column], projection="3d")
        draw_surface(ax, *surfaces[hemisphere], hemisphere_values[hemisphere], azimuth, title, vmax)


# Add the common symmetric Fisher-z color scale used by every panel.
def add_colorbar(fig, vmax: float, bounds: list[float]) -> None:
    axis = fig.add_axes(bounds)
    colorbar = fig.colorbar(
        ScalarMappable(norm=Normalize(-vmax, vmax), cmap="RdBu_r"),
        cax=axis,
        orientation="horizontal",
    )
    ticks = np.linspace(-vmax, vmax, 5)
    colorbar.set_ticks(ticks)
    colorbar.set_ticklabels([f"{value:.3g}" for value in ticks])
    colorbar.set_label("Stress × scan estimate (Fisher-z connectivity)", fontsize=10)


# Save both publication formats using the same canvas and close the figure.
def save_figure(fig, stem: Path, dpi: int) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        fig.savefig(
            stem.with_suffix(f".{extension}"), dpi=dpi, bbox_inches="tight", facecolor="white"
        )
    plt.close(fig)


# Create one four-view surface figure for a network or parcel map.
def plot_panel(
    args: argparse.Namespace,
    seed: int,
    values: np.ndarray,
    subtitle: str,
    stem: Path,
    surfaces: dict,
    vertex_labels: np.ndarray,
) -> None:
    fig = plt.figure(figsize=(13, 6.5), facecolor="white")
    grid = fig.add_gridspec(1, 4, left=0.02, right=0.98, bottom=0.18, top=0.80, wspace=0.02)
    draw_row(fig, grid, 0, values, vertex_labels, surfaces, args.vmax)
    fig.suptitle(f"Seed parcel {seed}: {SEED_NAMES[seed]}\n{subtitle}", fontsize=14, y=0.92)
    add_colorbar(fig, args.vmax, [0.26, 0.08, 0.48, 0.028])
    save_figure(fig, stem, args.dpi)


# %% 6. Draw the seven-network panels
# Save one network-level surface map per requested seed.
for seed in args.seeds:
    plot_panel(
        args,
        seed,
        values[seed]["network"],
        "Yeo-7 network stress × scan interaction estimates",
        args.output_root / "LMM" / "Yeo7" / f"seed_{seed:03d}_yeo7_interaction_surface",
        surfaces,
        vertex_labels,
    )


# %% 7. Draw the parcelwise effect and FDR panels
# Save the unthresholded parcel map and, when requested, its FDR-masked version.
for seed in args.seeds:
    for kind in args.map_types:
        count = np.isfinite(values[seed][kind]).sum()
        subtitle = "Parcelwise stress × scan interaction estimates across 399 targets"
        if kind == "fdr_significant":
            subtitle = f"FDR-BH significant parcelwise interactions (q < 0.05; n = {count})"
        plot_panel(
            args,
            seed,
            values[seed][kind],
            subtitle,
            args.output_root / "LMM" / "parcelwise_surfaces" / f"seed_{seed:03d}_parcelwise_{kind}",
            surfaces,
            vertex_labels,
        )


# %% 8. Assemble and save the supplementary figure
# All four rows share the same coefficient color scale.
# The composite contains network and unthresholded parcel rows in seed order.
if "interaction_estimates" in args.map_types:
    rows = [(seed, kind) for seed in args.seeds for kind in ("network", "interaction_estimates")]
    fig = plt.figure(figsize=(13, 3.5 * len(rows) + 1), facecolor="white")
    grid = fig.add_gridspec(
        len(rows), 4, left=0.02, right=0.98, bottom=0.09, top=0.95, hspace=0.55, wspace=0.02
    )
    for row, (seed, kind) in enumerate(rows):
        draw_row(fig, grid, row, values[seed][kind], vertex_labels, surfaces, args.vmax)
        subtitle = (
            "Yeo-7 network interaction estimates"
            if kind == "network"
            else "Parcelwise interaction estimates (399 targets)"
        )
        fig.text(
            0.04,
            grid[row, 0].get_position(fig).y1 + 0.015,
            f"{SEED_NAMES[seed]} — {subtitle}",
            fontsize=14,
            weight="bold",
        )
    add_colorbar(fig, args.vmax, [0.26, 0.045, 0.48, 0.012])
    save_figure(fig, args.output_root / "Supplementaryfigure_posthoc_seedbased", args.dpi)

print(f"Wrote PNG/SVG surface figures to {args.output_root} (common scale ±{args.vmax}).")
