# Step-by-step analysis of the swordfish paper

[中文说明](README.zh-CN.md)

These 32 steps repeat the analysis of the paper (`docs/swordfish_bycatch_manuscript.md`) and its supplement
(`docs/supplementary_tables.md`), from the public IATTC longline file to every figure, table and simulation,
except the four items listed at the end. The steps follow the order of the Methods. Each step is one short
script that is read and run from top to bottom. Every comment is given in English and then in Chinese.

Four files hold the parts that several steps share: `common.R` compares the new tables and figures with the
archived ones and draws the bootstrap samples, `model_functions.R` fits the fleet models,
`simulation_functions.R` sets up the simulation, and `plot_functions.R` holds the drawing settings.
`run_all.R` runs the 32 steps in order. All steps are R scripts.

## What to install

- R 4.6.0 with sdmTMB 1.2.0, TMB 1.9.25, fmesher 0.8.0, ggplot2 4.0.4, ggsci 5.2.0, sf 1.1.1, patchwork 1.3.2,
  rnaturalearth 1.2.0, rnaturalearthdata 1.0.0, jsonlite 2.0.0, cluster 2.1.8.2, generics 0.1.4, scales 1.4.0,
  ragg 1.5.2 and png 0.1.9. R notes that sdmTMB was built under R 4.6.1; this does not affect the results.

## Settings

Put these scripts in the folder `analysis_walkthrough/` of the project folder. The steps read the public
IATTC longline file (https://www.iattc.org/getmedia/b8f0bdbb-595d-4c16-9965-cbedcf122aaa/PublicLLTunaBillfish.zip) from `data/raw/`, the settings in `config/` and, for the checks,
the archived results in `results/` and `tmp/`.

Open a terminal in the project folder (the folder that holds `data/`, `config/` and `results/`) and set these
before running R, so that each R process uses one core and R can write to the Chinese folder names.
`run_all.R` sets them itself.

```sh
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

## Steps

Run the steps in order. Every step writes only to `analysis_walkthrough/output/`; the files in `data/`,
`results/`, `tmp/` and `config/` are only read. All steps read `config/directed_fishery.json`.

### Fisheries, observations and prediction areas (Section 2.1)

1. `01_check_public_data.R` reads the stored IATTC archive and `PublicLLTunaBillfishNum.csv` and checks that they are the same file. It writes nothing.
2. `02_draw_study_design.R` reads nothing and writes `Fig1_study_design` (Fig. 1).
3. `03_prepare_records.R` reads the public file and writes the H2 records of the six fleets, with the position and area of each cell, to `records.rds`.
4. `04_define_areas.R` reads `records.rds` and writes the years with effort in each cell (`support.csv`), the cells with positive swordfish counts (`footprint.csv`), the 31 fits (`tasks.csv`), the cells by assessment area (`area_structure.csv`, Table S1) and the model input `input.rds`.
5. `05_map_catch_cells.R` reads `footprint.csv` and writes `Fig2_coverage` (Fig. 2) and `coastal_grid_check.csv`.
6. `06_count_catch_cells.R` reads the public file and writes the cells with swordfish by period (`period_coverage.csv`) and the years with positive counts in the cells added by the tuna-targeting longline fleets (`added_cell_positive_years.csv`) (Section 3.1, Table S1 notes).

### Annual indices in the area of overlap with Spain (Section 2.2)

7. `07_fit_baseline_models.R` reads `input.rds`, fits the baseline model of each fleet to all eligible records and writes `fits_baseline.csv`, `indices_baseline.csv` and, in `fits/`, the diagnostics, residuals by cell and relative CPUE in each cell of each fit.
8. `08_spain_area_indices.R` reads `input.rds` and `indices_baseline.csv` and writes Spain's index by assessment area (`spain_area_indices.csv`, Fig. S1).
9. `09_compare_annual_indices.R` reads `indices_baseline.csv` and writes the rank correlations with Spain (`pair_comparisons.csv`, Table S1) and the precision of the annual indices (`index_precision.csv`, Table S3).
10. `10_count_annual_overlap.R` reads `input.rds` and `indices_baseline.csv` and writes the yearly cells with records from both fleets (`pair_overlap.csv`) and the China-Spain correlations over those years (`china_spain_overlap_correlations.csv`) (Section S1).
11. `11_plot_annual_indices.R` reads `indices_baseline.csv`, `support.csv` and `spain_area_indices.csv` and writes `Fig3_annual` (Fig. 3) and `FigS1_spain_areas` (Fig. S1).
12. `12_check_residuals.R` reads the files of step 07 in `fits/`, prints Table S2 and writes `spatial_residual_cell_years.csv`, `residual_pair_correlations.csv` (Section S4), `FigS2_fleet_residuals` (Fig. S2) and `FigS3_spatial_residuals` (Fig. S3).
13. `13_retrospective.R` reads `input.rds` and `indices_baseline.csv` and writes `retrospective_fits.csv`, `retrospective_indices.csv`, `retrospective_rho.csv` (Table S4) and `FigS4_retrospective` (Fig. S4).
14. `14_quarterly_change.R` reads `input.rds` and `fits_baseline.csv` and writes the China/Spain change ratios in the same cells and quarter (`quarterly_matched_cells.csv`, `quarterly_fit_checks.csv`, `quarterly_shared_change.csv`, `quarterly_between_seasons.csv`; Fig. 4b, Section S9).
15. `15_compare_cells_fished_by_both.R` reads the relative CPUE of step 07, `support.csv`, `pair_overlap.csv` and `quarterly_shared_change.csv` and writes `shared_cell_cpue_by_year.csv`, `shared_cell_direction_by_year.csv` (Table S11), `shared_cell_cpc_spain_ratios.csv`, `Fig4_local_contrasts` (Fig. 4), `FigS11_shared_cell_cpue` (Fig. S11), `FigS12_local_contrasts_2021_2023` (Fig. S12) and `shared_cell_cpue_atlas.pdf` (Supplementary maps).

### Sensitivity to the species composition of the catch (Section 2.3)

16. `16_describe_catch_composition.R` reads `input.rds` and writes the PCA and cluster selection (`composition_selection.csv`, Table S5), `china_cluster_swordfish.csv`, the loadings, cluster profiles and annual cluster shares (with and without swordfish) and `FigS5_fleet_pca`, `FigS6_fleet_clusters`, `FigS7_fleet_pca_inclusive` and `FigS8_fleet_clusters_inclusive` (Figs. S5 to S8).
17. `17_fit_composition_models.R` reads `input.rds`, `fits_baseline.csv` and `composition_selection.csv`, fits the other 25 models and writes, for all 31 fits, `fits.csv`, `indices.csv`, `model_diagnostics.csv` (Tables S2 and S6), `model_parameters.csv` (Table S2) and `residual_quantiles.csv`.
18. `18_summarize_composition.R` reads `indices.csv` and `fits.csv` and writes `composition_adjustment_summary.csv` (Table 3) and `composition_changes.csv`.
19. `19_plot_composition_indices.R` reads `indices.csv` and writes `Fig5_composition` (Fig. 5) and `FigS9_composition_inclusive` (Fig. S9).

### Spatial coverage and simulated regional abundance (Section 2.4)

20. `20_simulation_inputs.R` reads `input.rds` and the public file and writes the sampling designs (`coverage_sampling_designs.csv`, Table 2) and the simulated catch levels (`simulated_catch_levels.csv`, Table S7).
21. `21_simulate_coverage.R` reads `input.rds` and the public file, simulates scenarios A1 to A3 under the three designs (900 fits) and writes `coverage_simulation.csv`, `coverage_annual_intervals.csv` and `coverage_uncorrected_check.csv`.
22. `22_summarize_coverage.R` reads the tables of step 21 and writes `coverage_summary.csv` (Table S9), `coverage_annual_paired_summary.csv` (Fig. 6a) and `coverage_paired_errors.csv`.
23. `23_coverage_change_intervals.R` reads `coverage_simulation.csv` and `coverage_paired_errors.csv` and writes the median change with bootstrap intervals (`coverage_change_intervals.csv`, Fig. 6b).
24. `24_simulate_equal_records.R` reads `input.rds` and the public file, simulates the two designs with the same number of records (600 fits) and writes `equal_records_samples.csv`, `equal_records_fits.csv` and `equal_records_annual.csv`.
25. `25_summarize_equal_records.R` reads the tables of steps 21 and 24 and writes `equal_records_balance.csv`, `equal_records_pairs.csv` and `equal_records_summary.csv` (Table S8).
26. `26_simulate_decline_location.R` reads `input.rds`, the public file and `equal_records_samples.csv`, simulates the decline outside the core area, inside it and in both areas (A2, L1, L2; 900 fits) and writes `decline_location_fits.csv`, `decline_location_annual.csv` and `decline_location_truth.csv`.
27. `27_summarize_decline_location.R` reads the tables of steps 21, 24 and 26 and writes `decline_location_errors.csv` and `decline_location_summary.csv` (Section S11).
28. `28_plot_decline_location.R` reads `decline_location_summary.csv` and writes `FigS14_decline_location` (Fig. S14).
29. `29_simulate_catchability_trend.R` reads `input.rds` and the public file, simulates scenarios Q1 to Q6 (1,200 fits) and writes `catchability_trend_fits.csv`, `catchability_trend_annual.csv` and `catchability_trend_count_checks.csv`.
30. `30_summarize_catchability_trend.R` reads the tables of steps 21, 23 and 29 and writes `catchability_trend_summary.csv` (Table S12, Fig. 6b).
31. `31_compare_all_and_complete.R` reads the tables of steps 21, 24 and 29 and writes the medians over all replicates and over the replicates with all fits accepted (`coverage_all_vs_complete.csv`, Table S10 and the last column of Table S12).
32. `32_plot_coverage_experiment.R` reads the tables of steps 21, 22, 23 and 30 and writes `Fig6_coverage_experiment` (Fig. 6) and `FigS13_single_replicate` (Fig. S13).

Table 1 lists the scenarios that `simulation_functions.R` defines; it has no table of its own. Each figure is
written as a PDF and a 600 dpi PNG.

The simulation steps 21, 24, 26 and 29 read the number of replicates from the config (100). A smaller number
can be given for a quick run, for example `--replicates 10`. Replicates run in parallel on `--workers 8`
processes by default. Each replicate keeps the seed of the production run, so the results do not depend on the
number of workers. Each replicate is saved as one file under `output/replicates/`, so an interrupted run
resumes; delete that folder after any change to the code, the config or the input data.

## Commands

To run all 32 steps in order, each in a new R process with the settings above:

```sh
Rscript analysis_walkthrough/run_all.R --workers 8
```

`--workers` is passed to the simulation steps 21, 24, 26 and 29. The run stops at the first step that fails and
names it. To run the steps one by one:

```sh
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
Rscript analysis_walkthrough/01_check_public_data.R
Rscript analysis_walkthrough/02_draw_study_design.R
Rscript analysis_walkthrough/03_prepare_records.R
Rscript analysis_walkthrough/04_define_areas.R
Rscript analysis_walkthrough/05_map_catch_cells.R
Rscript analysis_walkthrough/06_count_catch_cells.R
Rscript analysis_walkthrough/07_fit_baseline_models.R
Rscript analysis_walkthrough/08_spain_area_indices.R
Rscript analysis_walkthrough/09_compare_annual_indices.R
Rscript analysis_walkthrough/10_count_annual_overlap.R
Rscript analysis_walkthrough/11_plot_annual_indices.R
Rscript analysis_walkthrough/12_check_residuals.R
Rscript analysis_walkthrough/13_retrospective.R --workers 5
Rscript analysis_walkthrough/14_quarterly_change.R
Rscript analysis_walkthrough/15_compare_cells_fished_by_both.R
Rscript analysis_walkthrough/16_describe_catch_composition.R
Rscript analysis_walkthrough/17_fit_composition_models.R --workers 4
Rscript analysis_walkthrough/18_summarize_composition.R
Rscript analysis_walkthrough/19_plot_composition_indices.R
Rscript analysis_walkthrough/20_simulation_inputs.R
Rscript analysis_walkthrough/21_simulate_coverage.R --workers 8
Rscript analysis_walkthrough/22_summarize_coverage.R
Rscript analysis_walkthrough/23_coverage_change_intervals.R
Rscript analysis_walkthrough/24_simulate_equal_records.R --workers 8
Rscript analysis_walkthrough/25_summarize_equal_records.R
Rscript analysis_walkthrough/26_simulate_decline_location.R --workers 8
Rscript analysis_walkthrough/27_summarize_decline_location.R
Rscript analysis_walkthrough/28_plot_decline_location.R
Rscript analysis_walkthrough/29_simulate_catchability_trend.R --workers 8
Rscript analysis_walkthrough/30_summarize_catchability_trend.R
Rscript analysis_walkthrough/31_compare_all_and_complete.R
Rscript analysis_walkthrough/32_plot_coverage_experiment.R
```

## Checks

Each step ends by reading its new tables and the archived tables of the same name in `results/directed_fishery/`
(and, for steps 03, 04, 21 and 29, the archived model input and check tables in `tmp/directed_fishery/`). It
prints "matches" or the columns that differ and the largest difference. Numbers match when their relative
difference is below about 1e-8, which absorbs the rounding to 15 significant digits when R writes a table. For
each PNG figure the step reports whether the file is identical to the archived figure, has the same pixels, or
the share of pixels that differ. A run with fewer replicates is compared with the same replicates of the
archived run, and its summaries are expected to differ.

The production scripts drew the bootstrap samples of steps 22, 25, 27 and 30 with the default random number
generator of NumPy. `numpy_rng()` in `common.R` repeats that generator in base R, so the same seeds give the
same samples and the same intervals.

Every step has been run and checked in this way: every table matched its archived table (the largest absolute
difference was 2.5e-12, in `model_diagnostics.csv`) and every PNG figure was identical to the archived file. The supplementary maps (`shared_cell_cpue_atlas.pdf`) are not compared, because a PDF records
the time it was made. Tables S2, S5 and S6 have no archived table of their own: step 12 prints Table S2 from the
parameter and diagnostic files, step 16 compares Table S5 with the candidates kept in
`results/directed_fishery/composition_figure_receipt*.json`, and Table S6 is the `delta_AIC_vs_common` column of
`model_diagnostics.csv`.

## Not covered

- The fishing-strategy simulation and Fig. S10: this archived run used sdmTMB 1.0.0 (results in `results/sdmtmb_full/` and `results/paper_figures/main_simulation_contrasts.csv`).
- The figure book and the Word files: they assemble the existing figures and text and calculate nothing.
- A new download of the public file: IATTC replaces it at each release, so step 01 checks the stored copy downloaded on 2026-11-06.
- The refit without the spatial field: no archived fit needed it, so steps 07 and 17 stop with a message if a fit is not accepted.
