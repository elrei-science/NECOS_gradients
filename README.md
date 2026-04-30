# NECOS published code

This repository contains the Python and R workflows used to regenerate the published NECOS gradient analyses and manuscript figures from the files shipped in `data/`.

## Licensing and data reuse

- **Code license**: all repository code is licensed under `GPL-3.0`. See `LICENSE`.
- **Derivative data license**: the derivative data shipped in `data/` are licensed under `CC BY 4.0`. See `DATA_LICENSE.md`.
- **Excluded source data**: raw and preprocessed MRI data from the NECOS study are not included and cannot be redistributed because of the original consent restrictions.
- **Participant identifiers**: `NECOS_ID` is an approved pseudonymous release identifier and is retained in the released derivative datasets.

## What is in this repository?

- `code/python/`: optional Python workflow for regenerating gradient derivatives.
- `code/r_scripts/eccentricity/`: whole-brain and parcel eccentricity analyses.
- `code/r_scripts/dispersion/`: between-network and within-network dispersion analyses.
- `code/r_scripts/peripheral/`: peripheral physiology/questionnaire analyses linking stress markers to brain outcomes.
- `code/r_scripts/figures/`: manuscript figure generation.
- `code/r_scripts/shared/`: shared R helper code used by the analysis scripts.
- `data/`: published inputs required by the R workflows.
- `results/`: generated outputs created by the scripts at runtime and excluded from version control by `.gitignore`.

## Project layout

- `code/python/NECOS_rs1_rs6_gradients.py`: regenerate gradient derivatives when you need to rebuild `data/necos_gradients/`.
- `code/r_scripts/shared/project_paths.R`: shared path and input helper for the R workflows.
- `code/r_scripts/eccentricity/eccentricity_primary_models.Rmd`: main eccentricity analysis.
- `code/r_scripts/eccentricity/eccentricity_posthoc_tests.R`: eccentricity posthoc summaries.
- `code/r_scripts/dispersion/dispersion_between_networks_primary.Rmd`: main between-network dispersion analysis.
- `code/r_scripts/dispersion/dispersion_within_networks_primary.Rmd`: within-network dispersion analysis.
- `code/r_scripts/dispersion/dispersion_posthoc_tests.Rmd`: dispersion posthoc analysis.
- `code/r_scripts/peripheral/peripheral_brain_association_models.Rmd`: peripheral marker analysis.
- `code/r_scripts/figures/publication_figures.R`: manuscript figure generation.

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

If `Rscript` is not available on your system but `R+` is, the equivalent form is:

```bash
R+ --version 4.5.0 R -q -e "renv::restore()"
```

After the first restore, opening R in the project root automatically activates the project because `.Rprofile` sources `renv/activate.R`.

On this machine, R was invoked through the environment wrapper:

```bash
R+ --version 4.5.0 R
```

## Recommended run order

Run the scripts from the project root so the shared helper resolves `data/` and `results/` correctly.

1. Optional: regenerate gradient derivatives with Python if you need to rebuild `data/necos_gradients/`.
2. Run `code/r_scripts/eccentricity/eccentricity_primary_models.Rmd`.
3. Run `code/r_scripts/eccentricity/eccentricity_posthoc_tests.R` if you need eccentricity posthoc summaries.
4. Run `code/r_scripts/dispersion/dispersion_between_networks_primary.Rmd`.
5. Run `code/r_scripts/dispersion/dispersion_within_networks_primary.Rmd` if you need within-network outputs.
6. Run `code/r_scripts/dispersion/dispersion_posthoc_tests.Rmd` if you need dispersion posthoc summaries.
7. Run `code/r_scripts/peripheral/peripheral_brain_association_models.Rmd`.
8. Run `code/r_scripts/figures/publication_figures.R` last.

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
Rscript code/r_scripts/figures/publication_figures.R
```

If `Rscript` is not available on your system but `R+` is, the equivalent form is:

```bash
R+ --version 4.5.0 R -q -e "rmarkdown::render('code/r_scripts/eccentricity/eccentricity_primary_models.Rmd', knit_root_dir = getwd())"
R+ --version 4.5.0 R -q -e "source('code/r_scripts/figures/publication_figures.R')"
```

## Outputs

- `results/eccentricity/`: eccentricity tables, model summaries, diagnostics.
- `results/dispersion/between_networks/`: between-network dispersion tables, model summaries, posthoc outputs, and diagnostics.
- `results/dispersion/within_networks/`: within-network dispersion tables, model summaries, and diagnostics.
- `results/peripheral/`: peripheral tables, model summaries, and diagnostics.
- `results/figures/`: exported manuscript figure files.


### Why is there a `data/peripheral_timeline_long.csv` file?

The file `data/peripheral_timeline_long.csv` is a figure-specific input used by the timeline panel in `code/r_scripts/figures/publication_figures.R`.

It is not a duplicate of the main baseline and delta tables used by the peripheral association models. Those analysis tables are scan-level summary inputs. By contrast, `peripheral_timeline_long.csv` preserves the finer-grained repeated-measures structure needed to recreate the stacked time-course figure.

The file combines the previously separate endocrine, autonomic, STAI, and VAS timeline tables into one tidy long-format dataset with these columns:

- `NECOS_ID`: participant identifier,
- `Group`: control vs stress group,
- `modality`: source family of the measurement,
- `measure`: specific measure name,
- `Timepoint`: repeated-measures timepoint label,
- `value`: observed value at that timepoint.

This single-file layout is simpler to publish and easier to validate downstream, while still retaining the original epoch-wise timeline information needed for the figure workflow.


## The shared helper script

The file `code/r_scripts/shared/project_paths.R` exists to keep repeated path and input logic out of every analysis script.

It does four main jobs:

1. **Locate the project root**
  - `pc_locate_root()` walks upward from the current working directory until it finds a folder containing `code/`, `data/`, and `results/`.
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

In practice, this means the analysis scripts can stay short and call the helper instead of hard-coding absolute paths or re-implementing the same file-reading logic repeatedly.

More concretely, `project_paths.R` acts as the contract between the published repository layout and the analysis code: it defines where inputs live, where outputs should go, and how common NECOS gradient files should be interpreted once loaded into R.


## Notes

- The R scripts create only the `results/` subfolders they actually need, when they need them.
- The figure workflow writes exported figures directly into `results/figures/`.
- The primary workflows use `diagnostics/` for QC outputs and `figures/` for substantive plots.
- If a script fails with a missing-input error, check whether the previous script in the run order produced the required intermediate file.
- Rendered `.html` reports are generated artifacts and are intentionally excluded from version control.
- The repository is set up for code-plus-data release; local environments, caches, rendered reports, and generated `results/` outputs are intentionally excluded from version control.
