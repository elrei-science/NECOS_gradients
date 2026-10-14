#!/usr/bin/env python3
r"""Prepare the NECOS rest2/rest3 post-hoc FC analyses.

CONTENTS
  A. Set-up
  B. Seed-based connectivity
  C. Mean-edge average connectivity
  D. Subject-label permutation test

Run this file before the matching R mixed-model script and the plotting script.
The paths below are placeholders and can be replaced in an IDE, or supplied as
command-line options. The script contains two seed maps, seven within-network means, three
planned between-network means, and the three edge-heterogeneity permutations.

Example:
    python NECOS_poshoc_seedbased_meannetworkconnectivity_rs2_r3.py \
        --cache-root /path/to/timeseries --stage all

Each time-series file must be a NumPy array with shape time points x 400
parcels, in Schaefer-400/Yeo-7 order.
"""

# %% A. Set-up
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


# Editable analysis settings.
SEEDS = (100, 377)  # left insula/frontal operculum; right vPFC
RESTS = (2, 3)
N_PARCELS = 400
NETWORK_ORDER = (
    "Visual",
    "Somatomotor",
    "Dorsal Attn",
    "Salience / Ventral Attn",
    "Limbic",
    "Control",
    "Default",
)
NETWORK_LABELS = {
    "SN": "Salience / Ventral Attn",
    "DMN": "Default",
    "FPN": "Control",
}
PAIRS = (("SN", "DMN", "SN_DMN"), ("SN", "FPN", "SN_FPN"), ("DMN", "FPN", "DMN_FPN"))
TIMESERIES_FILENAME = "timeseries_hp_filt_2mni_6fwhm_timeseries.npy"
PARCEL_COLUMNS = [f"parcel_{parcel:03d}" for parcel in range(1, N_PARCELS + 1)]

# Set these in an IDE if desired. None means use the repository-relative path.
PROJECT_ROOT = None
CACHE_ROOT = None
RSFC_ROOT = None
MOTION_FILE = None
LABEL_FILE = None
SEEDBASED_ROOT = None
NETWORK_ROOT = None
RUN_STAGE = "all"  # "prepare" writes connectivity; "followups" runs permutations.
N_PERMUTATIONS = 5000
RANDOM_SEED = 20260903


# Find the project directory so defaults work both from the repository and from
# an IDE cell. A command-line --project-root can override this discovery.
def find_project_root() -> Path:
    """Find the repository when this file is run as a script or from an IDE."""
    start = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
    for candidate in (start, *start.parents):
        if (candidate / "code").is_dir() and (candidate / "data").is_dir():
            return candidate
    return Path.cwd()


# Collect paths and run settings. Keeping this separate makes the editable
# settings above easy to inspect before any data are read.
def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT or find_project_root())
    parser.add_argument("--cache-root", type=Path, default=CACHE_ROOT)
    parser.add_argument("--rsfc-root", type=Path, default=RSFC_ROOT)
    parser.add_argument("--timeseries-filename", default=TIMESERIES_FILENAME)
    parser.add_argument("--motion-file", type=Path, default=MOTION_FILE)
    parser.add_argument("--label-file", type=Path, default=LABEL_FILE)
    parser.add_argument("--seedbased-root", type=Path, default=SEEDBASED_ROOT)
    parser.add_argument("--network-root", type=Path, default=NETWORK_ROOT)
    parser.add_argument("--stage", choices=("all", "prepare", "followups"), default=RUN_STAGE)
    parser.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS)
    parser.add_argument("--random-seed", type=int, default=RANDOM_SEED)
    args = parser.parse_args(argv)
    if args.n_permutations < 1:
        parser.error("--n-permutations must be positive")

    root = args.project_root.resolve()
    args.cache_root = args.cache_root or root / "data" / "timeseries"
    args.motion_file = args.motion_file or root / "data" / "mean_FD.csv"
    args.label_file = args.label_file or root / "data" / "yeo_networks_17_reordered_400_complete_labels.csv"
    args.seedbased_root = args.seedbased_root or root / "results" / "posthoc_seedbased"
    args.network_root = args.network_root or root / "results" / "posthoc_withinbetween_fc"
    return args


# Ignore notebook/kernel arguments when this file is run cell by cell.
args = parse_args([] if "get_ipython" in globals() or "__file__" not in globals() else None)


# Use one writer for every tab-separated result so all tables have the same
# delimiter, missing-value spelling, and parent-directory handling.
def write_tsv(data: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(path, sep="\t", index=False, na_rep="NaN")


# Read the participant covariates once. These same stress and motion values are
# used by the R models and by the permutation test below.
motion = pd.read_csv(args.motion_file)
required_motion = {"ID", "stress", "rest2", "rest3"}
if not required_motion.issubset(motion.columns):
    raise ValueError(f"Motion table requires columns: {sorted(required_motion)}")
if motion.ID.isna().any() or motion.ID.duplicated().any():
    raise ValueError("Participant IDs must be present and unique.")
if not motion.stress.isin((0, 1)).all() or set(motion.stress) != {0, 1}:
    raise ValueError("stress must contain both groups: 0=control and 1=stress.")
if not np.isfinite(motion[["stress", "rest2", "rest3"]].to_numpy(float)).all():
    raise ValueError("stress, rest2, and rest3 must be finite for every participant.")
motion = motion.sort_values("ID").set_index("ID")

# Read the parcel labels and convert each Yeo network into zero-based indices
# for NumPy indexing.
labels = pd.read_csv(args.label_file)
if not {"7networks_order", "7networks"}.issubset(labels.columns):
    raise ValueError("Label table requires 7networks_order and 7networks.")
if sorted(labels["7networks_order"].tolist()) != list(range(1, N_PARCELS + 1)):
    raise ValueError("Label table must contain each parcel ID 1..400 exactly once.")
if set(labels["7networks"]) != set(NETWORK_ORDER):
    raise ValueError("Label table must contain the expected seven Yeo networks.")
membership = {
    name: labels.loc[labels["7networks"] == name, "7networks_order"].to_numpy(int) - 1
    for name in NETWORK_ORDER
}


# Load one participant's scan, checking the array shape and values before it is
# used in a correlation. The optional rsfc-root supports an alternative folder layout.
def load_timeseries(subject: str, rest: int) -> tuple[np.ndarray, str, Path]:
    """Load one scan from the cache, with an optional alternative folder layout."""
    candidates = [("cached_timeseries", args.cache_root / f"rs_{rest}" / subject / args.timeseries_filename)]
    if args.rsfc_root is not None:
        candidates.append(
            ("subject_functional", args.rsfc_root / subject / "functional" / f"rest{rest}" / args.timeseries_filename)
        )
    for source, path in candidates:
        if path.is_file():
            values = np.asarray(np.load(path, allow_pickle=False), dtype=float)
            if values.ndim != 2 or values.shape[1] != N_PARCELS or values.shape[0] < 2:
                raise ValueError(f"Expected timepoints x 400: {path}")
            if not np.isfinite(values).all() or np.any(values.std(axis=0) == 0):
                raise ValueError(f"Nonfinite data or a zero-variance parcel: {path}")
            return values, source, path
    checked = ", ".join(str(path) for _, path in candidates)
    raise FileNotFoundError(f"Missing rest{rest} time series for {subject}; checked {checked}")


# Fail early if any paired scan is missing, rather than silently changing the
# participant sample for one of the multiple-testing families.
for rest in RESTS:
    for subject in motion.index:
        load_timeseries(subject, rest)


# Calculate the two values needed for seed analysis: ordinary Pearson r for the
# optional diagnostic table and Fisher-z for the models and network averages.
def seed_connectivity(timeseries: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Return Pearson r and Fisher-z seed-to-all-parcel maps."""
    centered = timeseries - timeseries.mean(axis=0, keepdims=True)
    standardized = centered / centered.std(axis=0, ddof=1)
    correlations = (standardized[:, seed - 1, None] * standardized).sum(axis=0) / (len(timeseries) - 1)
    correlations = np.clip(correlations, -0.999999, 0.999999)
    correlations[seed - 1] = np.nan
    return correlations, np.arctanh(correlations)


# Calculate the complete 400-by-400 Fisher-z matrix used by the edge-averaged
# within- and between-network analysis. The diagonal is excluded from averaging.
def fisher_z_matrix(timeseries: np.ndarray) -> np.ndarray:
    """Calculate parcel-to-parcel correlations and Fisher-transform them."""
    correlation = np.corrcoef(timeseries, rowvar=False)
    fisher_z = np.arctanh(np.clip(correlation, -0.999999, 0.999999))
    np.fill_diagonal(fisher_z, np.nan)
    return fisher_z


# %% B. Seed-based connectivity
# For each seed, participant, and scan we retain 399 target values after
# excluding the seed parcel. The R script averages these into the seven Yeo
# networks and fits the stress x scan mixed models.
connectivity_rows = []
network_edges = {code: {rest: [] for rest in RESTS} for _, _, code in PAIRS}
availability_rows = []
run_summary = []

if args.stage in ("all", "prepare"):
    # Build one participant-by-target table for each seed and resting state.
    for rest in RESTS:
        seed_rows = {seed: {"r": [], "fisher_z": []} for seed in SEEDS}
        scan_availability = []
        for subject in motion.index:
            timeseries, source, path = load_timeseries(subject, rest)
            scan_availability.append([subject, rest, "available", source, str(path)])
            for seed in SEEDS:
                correlations, fisher_z = seed_connectivity(timeseries, seed)
                metadata = [subject, seed, source, str(path)]
                seed_rows[seed]["r"].append(metadata + correlations.tolist())
                seed_rows[seed]["fisher_z"].append(metadata + fisher_z.tolist())

        rs_output = args.seedbased_root / f"rs_{rest}"
        write_tsv(pd.DataFrame(scan_availability,
                               columns=["subject", "resting_state", "status", "source", "path_checked"]),
                  rs_output / "timeseries_availability.tsv")
        group_tables = []
        for seed in SEEDS:
            seed_output = rs_output / f"seed_{seed:03d}"
            columns = ["subject", "seed_parcel", "source", "timeseries_path", *PARCEL_COLUMNS]
            for scale in ("r", "fisher_z"):
                write_tsv(pd.DataFrame(seed_rows[seed][scale], columns=columns),
                          seed_output / f"subject_seed_maps_{scale}.tsv")
            values = np.asarray([row[4:] for row in seed_rows[seed]["fisher_z"]], float)
            counts = np.isfinite(values).sum(axis=0)
            means = np.divide(np.nansum(values, axis=0), counts,
                              out=np.full(N_PARCELS, np.nan), where=counts > 0)
            group = pd.DataFrame({
                "seed_parcel": seed,
                "seed_label": f"Schaefer400 parcel {seed}",
                "target_parcel": np.arange(1, N_PARCELS + 1),
                "mean_fisher_z": means,
                "mean_r": np.tanh(means),
                "n_subjects": counts,
            })
            write_tsv(group, seed_output / "group_seed_map.tsv")
            group_tables.append(group)
        write_tsv(pd.concat(group_tables, ignore_index=True), rs_output / "group_seed_maps.tsv")
        run_summary.append([rest, len(motion), 0, ",".join(map(str, SEEDS))])
        print(f"Prepared rest{rest}: {len(motion)} participants.", flush=True)


# %% C. Mean-edge average connectivity
# These are Fisher-z averages over unique parcel edges: seven within-network
# measures plus SN-DMN, SN-FPN, and DMN-FPN. They are the ten outcomes used by
# the R mixed-model script; no network-averaged time-series correlations are used.
if args.stage in ("all", "prepare"):
    # Re-read each scan for the parcel-edge summaries. Keeping this loop separate
    # makes the distinction between seed-based and edge-averaged FC explicit.
    for rest in RESTS:
        for subject in motion.index:
            timeseries, _, path = load_timeseries(subject, rest)
            availability_rows.append([subject, rest, "available", str(path)])
            z_fc = fisher_z_matrix(timeseries)
            for network, indices in membership.items():
                edge_values = z_fc[np.ix_(indices, indices)][np.triu_indices(len(indices), k=1)]
                measure = "within_" + network.replace(" / ", "_").replace(" ", "_")
                connectivity_rows.append(
                    [subject, rest, f"rest{rest}", "within", measure, network,
                     edge_values.size, edge_values.mean(), np.tanh(edge_values.mean())]
                )
            for first, second, code in PAIRS:
                edge_values = z_fc[np.ix_(membership[NETWORK_LABELS[first]], membership[NETWORK_LABELS[second]])]
                network_edges[code][rest].append(edge_values.ravel())
                connectivity_rows.append(
                    [subject, rest, f"rest{rest}", "between", code,
                     f"{NETWORK_LABELS[first]} - {NETWORK_LABELS[second]}",
                     edge_values.size, edge_values.mean(), np.tanh(edge_values.mean())]
                )

    write_tsv(
        pd.DataFrame(run_summary, columns=["resting_state", "n_available_subjects", "n_missing_subjects", "seed_parcels"]),
        args.seedbased_root / "run_summary.tsv",
    )
    table_dir = args.network_root / "tables"
    write_tsv(
        pd.DataFrame(connectivity_rows, columns=[
            "subject", "resting_state", "scan", "measure_type", "measure",
            "network_or_pair", "n_parcel_pairs", "mean_fisher_z", "mean_r"
        ]),
        table_dir / "subject_network_connectivity_fisher_z.tsv",
    )
    write_tsv(pd.DataFrame(availability_rows,
                           columns=["subject", "resting_state", "status", "timeseries_path"]),
              table_dir / "timeseries_availability.tsv")
    write_tsv(pd.DataFrame([[name, len(indices)] for name, indices in membership.items()],
                           columns=["network", "n_parcels"]),
              table_dir / "network_definition.tsv")


# %% D. Subject-label permutation test
# The reported follow-up asks whether edgewise stress x scan effects are more
# variable than expected under a subject-label null. For each planned pair we
# regress every edge's rest3-minus-rest2 change on stress and delta FD, then
# compare the observed SD of stress slopes with 5,000 shared label shuffles.
if args.stage in ("all", "followups"):
    # A followups-only run reconstructs the pairwise edge arrays; an all-stage
    # run already collected them while writing the mean-edge table.
    # Follow-up-only runs still need the edge arrays, so load them here.
    if not any(network_edges[code][rest] for _, _, code in PAIRS for rest in RESTS):
        for rest in RESTS:
            for subject in motion.index:
                timeseries, _, _ = load_timeseries(subject, rest)
                z_fc = fisher_z_matrix(timeseries)
                for first, second, code in PAIRS:
                    network_edges[code][rest].append(
                        z_fc[np.ix_(membership[NETWORK_LABELS[first]], membership[NETWORK_LABELS[second]])].ravel()
                    )

    # Use the rest3-minus-rest2 change for each edge and adjust every slope for
    # the participant's change in mean framewise displacement.
    stress = motion.stress.to_numpy(float)
    delta_fd = (motion.rest3 - motion.rest2).to_numpy(float)
    design = np.column_stack((np.ones(len(stress)), stress, delta_fd))
    if len(stress) <= 3 or np.linalg.matrix_rank(design) != 3:
        raise ValueError("Permutation test requires >3 participants and a full-rank design.")
    inverse = np.linalg.inv(design.T @ design)
    rng = np.random.default_rng(args.random_seed)
    result_rows = []
    # The same shuffled subject labels are applied jointly to all edges within a
    # network pair, preserving the dependence among those edge estimates.
    for _, _, code in PAIRS:
        changes = np.asarray(network_edges[code][3]) - np.asarray(network_edges[code][2])
        observed_slopes = (inverse @ design.T @ changes)[1]
        observed_sd = np.std(observed_slopes, ddof=1)
        null_sd = np.empty(args.n_permutations)
        for permutation in range(args.n_permutations):
            shuffled_stress = rng.permutation(stress)
            shuffled_design = np.column_stack((np.ones(len(stress)), shuffled_stress, delta_fd))
            shuffled_inverse = np.linalg.inv(shuffled_design.T @ shuffled_design)
            null_sd[permutation] = np.std((shuffled_inverse @ shuffled_design.T @ changes)[1], ddof=1)
        p_value = (1 + np.count_nonzero(null_sd >= observed_sd)) / (args.n_permutations + 1)
        result_rows.append([code, changes.shape[1], observed_sd, p_value])

    permutation_results = pd.DataFrame(
        result_rows,
        columns=["pair", "n_edges", "edge_slope_sd", "p_edge_heterogeneity_permutation"],
    )
    p_values = permutation_results.p_edge_heterogeneity_permutation.to_numpy()
    order = np.argsort(p_values)
    adjusted = np.minimum.accumulate((p_values[order] * len(order)
                                      / np.arange(1, len(order) + 1))[::-1])[::-1]
    fdr = np.empty(len(order))
    fdr[order] = np.minimum(adjusted, 1.0)
    permutation_results["p_edge_heterogeneity_fdr_bh"] = fdr
    followup_dir = args.network_root / "network_reconfiguration_followups"
    write_tsv(permutation_results, followup_dir / "edge_heterogeneity_permutation_tests.tsv")
    (followup_dir / "README.txt").write_text(
        "Subject-label permutations test the SD of edgewise stress x scan slopes "
        "for SN-DMN, SN-FPN, and DMN-FPN. Each model adjusts for delta mean FD. "
        f"Permutations: {args.n_permutations}; random seed: {args.random_seed}. "
        "Benjamini-Hochberg correction is across the three planned pairs.\n"
    )
    print(permutation_results.to_string(index=False), flush=True)


# Save a small record of the paths and settings used for this run so the exact
# input locations, stage, and permutation seed can be recovered later.
args.network_root.mkdir(parents=True, exist_ok=True)
metadata = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
metadata["participants"] = motion.index.tolist()
(args.network_root / f"run_metadata_{args.stage}.json").write_text(json.dumps(metadata, indent=2) + "\n")
print("Finished. Next: run code/r_scripts/functional_connectivity/*_lmm.R, then code/r_scripts/figures/supplementary_figures.py.")
