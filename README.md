# NECOS gradients: Central and peripheral dynamics of acute stress: evidence from functional cortical gradients

This repository contains the Python and R workflows used to regenerate the NECOS gradient analyses and manuscript figures from the files shipped in `data/`.

## Licensing and data reuse

- **Code license**: all repository code is licensed under `GPL-3.0`. See `LICENSE`.
- **Derivative data license**: the derivative data shipped in `data/` are licensed under `CC BY 4.0`. See `DATA_LICENSE.md`.
- **Excluded source data**: raw and preprocessed MRI data from the NECOS study are not included and cannot be redistributed because of the original consent restrictions.
- **Participant identifiers**: `NECOS_ID` is an approved pseudonymous release identifier and is retained in the released derivative datasets.

## What is in this repository?

- `code/python/`: gradient derivation and post-hoc seed-based and network connectivity preparation.
- `code/r_scripts/eccentricity/`: whole-brain and parcel eccentricity analyses.
- `code/r_scripts/dispersion/`: between-network and within-network dispersion analyses.
- `code/r_scripts/peripheral/`: peripheral physiology/questionnaire analyses linking stress markers (peak change and AUCi) to brain outcomes.
- `code/r_scripts/functional_connectivity/`: mixed models for the post-hoc seed-based and network connectivity analyses.
- `code/r_scripts/figures/`: manuscript and supplementary figure generation.
- `code/r_scripts/shared/`: shared R helper code (project paths and common input/output functions).
- `data/`: published inputs required by the R workflows.
- `results/`: generated outputs created by the scripts at runtime and excluded from version control by `.gitignore`.

## Project layout

- `code/python/NECOS_rs1_rs6_gradients.py`: gradient derivatives for all six resting states (`data/necos_gradients/`).
- `code/python/NECOS_poshoc_seedbased_meannetworkconnectivity_rs2_r3.py`: seed-based and mean network connectivity (rest2, rest3) and the edge-heterogeneity permutation test.
- `code/r_scripts/functional_connectivity/NECOS_poshoc_seedbased_meannetworkconnectivity_rs2_r3_lmm.R`: mixed models for the seed-based, seed-to-network and within/between-network connectivity outcomes.
- `code/r_scripts/shared/project_io.R`: project paths and common input/output functions used by the R workflows.
- `code/r_scripts/eccentricity/eccentricity_primary_models.Rmd`: main eccentricity analysis.
- `code/r_scripts/eccentricity/eccentricity_posthoc_tests.R`: eccentricity posthoc summaries.
- `code/r_scripts/dispersion/dispersion_between_networks_primary.Rmd`: main between-network dispersion analysis.
- `code/r_scripts/dispersion/dispersion_within_networks_primary.Rmd`: within-network dispersion analysis.
- `code/r_scripts/dispersion/dispersion_posthoc_tests.Rmd`: dispersion posthoc analysis.
- `code/r_scripts/peripheral/peripheral_brain_association_models.Rmd`: peripheral marker analysis.
- `code/r_scripts/peripheral/peripheral_auci_brain_associations.R`: AUCi of the six peripheral measures and their Group x AUCi association with the five brain outcomes (four gradient outcomes and vPFC seed-based connectivity change).
- `code/r_scripts/figures/publication_figures.R`: manuscript figure generation.
- `code/r_scripts/figures/supplementary_figures.py`: supplementary figure of the post-hoc seed-based results.

## Environment setup

### Python

The Python workflow is managed through `code/python/pyproject.toml`.

```bash
cd code/python
python3 -m venv .venv-necos-rs
source .venv-necos-rs/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Optional Jupyter support:

```bash
python -m pip install -e .[dev]
```

### R

The R workflows are managed with `renv`.

From the project root, restore the recorded package environment with:

```bash
Rscript -e "renv::restore()"
```

After the first restore, opening R in the project root automatically activates the project because `.Rprofile` sources `renv/activate.R`.

## Paths to set

All scripts find the repository root themselves and read from `data/` and write to `results/`; run them from the project root. Paths outside the repository must be set in these places:

| Script | Setting | Points to |
|---|---|---|
| `code/python/NECOS_rs1_rs6_gradients.py` | `INPUT_DATA_ROOT`, `OUTPUT_ROOT` (top of the script) | Preprocessed resting-state data, `<root>/NECOS001/functional/rest<n>/rest_preprocessed_hp_filt_2mni_6fwhm.nii.gz`, and the output folder. `--dry-run` checks both. |
| `code/python/NECOS_poshoc_seedbased_meannetworkconnectivity_rs2_r3.py` | `--cache-root` (or `CACHE_ROOT` in the script) | Cached parcel time series, `<cache-root>/rs_<n>/NECOS001/<file>.npy`. |
| `code/r_scripts/peripheral/peripheral_auci_brain_associations.R` | `seed_dir` (top of the script) | Seed-based results folder written by the connectivity workflow (`results/posthoc_seedbased`). |
| `code/r_scripts/figures/supplementary_figures.py` | `--brainspace-data` (only if BrainSpace is not installed) | BrainSpace `datasets` folder with the surface and parcellation files. |

No paths need to be set for the eccentricity, dispersion, peripheral association and publication figure scripts, or for the connectivity mixed models.

## Recommended run order

Run the scripts from the project root so the shared helper resolves `data/` and `results/` correctly.

1. Create the gradients with `code/python/NECOS_rs1_rs6_gradients.py`. This requires the preprocessed resting-state data, which are not distributed; its outputs are shipped in `data/necos_gradients/`, so all later steps can run without it.
2. Run `code/r_scripts/eccentricity/eccentricity_primary_models.Rmd`.
3. Run `code/r_scripts/eccentricity/eccentricity_posthoc_tests.R` for the eccentricity posthoc summaries.
4. Run `code/r_scripts/dispersion/dispersion_between_networks_primary.Rmd`.
5. Run `code/r_scripts/dispersion/dispersion_within_networks_primary.Rmd` for the within-network outputs.
6. Run `code/r_scripts/dispersion/dispersion_posthoc_tests.Rmd` for the dispersion posthoc summaries.
7. Run `code/r_scripts/peripheral/peripheral_brain_association_models.Rmd`.
8. Run the post-hoc connectivity workflow: `NECOS_poshoc_seedbased_meannetworkconnectivity_rs2_r3.py`, then `NECOS_poshoc_seedbased_meannetworkconnectivity_rs2_r3_lmm.R`. This requires cached parcel time series, which are not distributed (see [Post-hoc connectivity](#post-hoc-seed-based-and-network-connectivity)).
9. Run `code/r_scripts/peripheral/peripheral_auci_brain_associations.R`. It reads the outputs of steps 7 and 8; set `seed_dir` at the top of the script to the seed-based results folder from step 8 (`results/posthoc_seedbased`).
10. Run `code/r_scripts/figures/publication_figures.R`.
11. Run `code/r_scripts/figures/supplementary_figures.py`.

## How to run the scripts

### Python gradients

```bash
cd code/python
source .venv-necos-rs/bin/activate
python NECOS_rs1_rs6_gradients.py --dry-run
python NECOS_rs1_rs6_gradients.py
```

### R scripts

Example commands from the project root:

```bash
Rscript -e "rmarkdown::render('code/r_scripts/eccentricity/eccentricity_primary_models.Rmd')"
Rscript -e "rmarkdown::render('code/r_scripts/dispersion/dispersion_between_networks_primary.Rmd')"
Rscript -e "rmarkdown::render('code/r_scripts/peripheral/peripheral_brain_association_models.Rmd')"
Rscript code/r_scripts/functional_connectivity/NECOS_poshoc_seedbased_meannetworkconnectivity_rs2_r3_lmm.R
Rscript code/r_scripts/peripheral/peripheral_auci_brain_associations.R
Rscript code/r_scripts/figures/publication_figures.R
```

### Post-hoc connectivity and supplementary figure

From the repository root, with the Python environment active:

```bash
python code/python/NECOS_poshoc_seedbased_meannetworkconnectivity_rs2_r3.py --cache-root /path/to/time_series
python code/r_scripts/figures/supplementary_figures.py
```

The Python and R connectivity scripts accept `--project-root`, `--seedbased-root` and `--network-root` to read and write results elsewhere. `--stage` selects part of a script: `prepare` or `followups` for the Python script, and `parcel`, `seed-network` or `mean-network` for the R script (`seed-network` needs the datasets written by `parcel`). The supplementary figure uses the BrainSpace surface files; if BrainSpace is not installed, pass `--brainspace-data /path/to/brainspace/datasets` containing `surfaces/conte69_32k_{lh,rh}.gii` and `parcellations/schaefer_400_conte69.csv`.

## Post-hoc seed-based and network connectivity

**Inputs.** `data/mean_FD.csv` and `data/yeo_networks_17_reordered_400_complete_labels.csv` are shipped. The cached parcel time series are not distributed because of the consent restrictions on the MRI data. They are read from `<cache-root>/rs_2/NECOS001/timeseries_hp_filt_2mni_6fwhm_timeseries.npy` (and `rs_3`, all participants) as arrays of time points x 400 parcels in Schaefer-400 Yeo-7 order. `--cache-root` sets the location and `--timeseries-filename` the file name.

**Analyses.**

- Seeds: parcel 377 (right ventral PFC) and parcel 100 (left insula/frontal operculum).
- Seed-to-parcel connectivity is the Pearson correlation followed by Fisher z, excluding the seed. 399 models per seed, BH-FDR within seed.
- Seed-to-network connectivity is the mean Fisher-z edge within each Yeo-7 network. Seven models per seed, BH-FDR within seed.
- Within-network connectivity (seven networks) and between-network connectivity (SN-DMN, SN-FPN, DMN-FPN) are means of the corresponding edges. Ten models, one BH-FDR family.
- All models: `connectivity ~ stress * scan + headmotion + (1 | ID)`, REML, Satterthwaite p-values (lmerTest), reference scan rest2. A positive `stress:scanrest3` coefficient is a larger change in the stress group.
- Edge heterogeneity: the SD across edgewise stress coefficients from `delta_FC ~ stress + delta_FD`, tested with 5,000 subject-label permutations (seed 20260903), BH-FDR across the three network pairs.

## Outputs

- `results/eccentricity/`: eccentricity tables, model summaries, diagnostics.
- `results/dispersion/between_networks/`: between-network dispersion tables, model summaries, posthoc outputs, and diagnostics.
- `results/dispersion/within_networks/`: within-network dispersion tables, model summaries, and diagnostics.
- `results/peripheral/`: peripheral tables, model summaries, and diagnostics.
- `results/peripheral_auci/`: AUCi values, the vPFC connectivity change, and the AUCi x brain association table.
- `results/posthoc_seedbased/`: participant seed maps, parcelwise and seed-to-network mixed-model results (`LMM/`), and the supplementary figure.
- `results/posthoc_withinbetween_fc/`: edge-averaged connectivity inputs, mixed-model results, and the edge-heterogeneity permutation tests (`network_reconfiguration_followups/edge_heterogeneity_permutation_tests.tsv`).
- `results/figures/`: exported manuscript figure files.


### Peripheral timeline data

The file `data/peripheral_timeline_long.csv` is the input for the timeline panel in `code/r_scripts/figures/publication_figures.R` and for the AUCi calculation in `code/r_scripts/peripheral/peripheral_auci_brain_associations.R`.

It is not a duplicate of the main baseline and delta tables used by the peripheral association models. Those analysis tables are scan-level summary inputs. By contrast, `peripheral_timeline_long.csv` preserves the finer-grained repeated-measures structure needed to recreate the stacked time-course figure and to compute the AUCi.

The file combines the previously separate endocrine, autonomic, STAI, and VAS timeline tables into one tidy long-format dataset with these columns:

- `NECOS_ID`: participant identifier,
- `Group`: control vs stress group,
- `modality`: source family of the measurement,
- `measure`: specific measure name,
- `Timepoint`: repeated-measures timepoint label,
- `value`: observed value at that timepoint.

This single-file layout is simpler to publish and easier to validate downstream, while still retaining the original epoch-wise timeline information needed for the figure workflow.


## Shared helper script

The file `code/r_scripts/shared/project_io.R` holds the path and input/output functions shared by the R analysis scripts.

It provides five groups of functions:

1. **Locate the project root**
  - `pc_locate_root()` walks upward from the current working directory until it finds a folder containing `code/` and `data/`.
  - `pc_source_context_root()` tries a few possible starting points, such as the current working directory and the currently rendered R Markdown file.

2. **Create consistent output paths**
  - `pc_analysis_paths("eccentricity")` returns a list of directories like:
    - `project_root`
    - `code_root`
    - `data_root`
    - `results_root`
    - `analysis_root`
    - `tables_dir`
    - `models_dir`
    - `figures_dir`
    - `diagnostics_dir`
  - It does **not** create those directories eagerly.
  - Scripts create only the specific folders they need via `pc_ensure_dirs(...)`, which keeps the published `results/` tree clean until outputs are actually written.
  - The intention is that each analysis script can declare its output contract once at the top, then write consistently to `tables`, `models`, `figures`, and `diagnostics` without duplicating path boilerplate.

3. **Load common published inputs**
  - `pc_mean_fd()` reads `data/mean_FD.csv` and reshapes it into scan-level head-motion rows.
  - `pc_parcel_labels()` reads `data/yeo_networks_17_reordered_400_complete_labels.csv` and normalizes the label columns so the scripts can use a consistent schema even if the CSV headers differ slightly.

4. **Load published gradient files**
  - `pc_list_gradient_files()` finds the BIDS-style gradient TSV files under `data/necos_gradients/derivatives/gradients/`.
  - `pc_read_gradient_file()` reads one gradient TSV and returns parcel-wise `V1`, `V2`, and `V3` values.
  - `pc_read_session_gradient_matrices()` returns a session bundle with participant IDs and 3-column gradient matrices.
  - `pc_read_all_session_gradient_matrices()` loads all requested sessions at once.

5. **Write outputs safely**
  - `pc_write_csv()` creates the target parent directory if needed and then writes the CSV.
  - This avoids repeated `dir.create()` calls scattered throughout the analysis code.


## Notes

- The R scripts create only the `results/` subfolders they actually need, when they need them.
- The figure workflow writes exported figures directly into `results/figures/`.
- The primary workflows use `diagnostics/` for QC outputs and `figures/` for substantive plots.
- If a script fails with a missing-input error, check whether the previous script in the run order produced the required intermediate file.
- Rendered `.html` reports are generated artifacts and are intentionally excluded from version control.
- The repository is set up for code-plus-data release; local environments, caches, rendered reports, and generated `results/` outputs are intentionally excluded from version control.
