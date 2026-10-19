# Theory-space neutrino calculations

Mathematica notebooks for generating benchmark data, plotting neutrino masses and mixing angles, and analysing saved results. Each notebook has a separate task.

| Folder | Contents |
| --- | --- |
| `production_parameters/` | Exact parameter associations for all 36 production benchmarks and loading instructions. |
| `generation/` | Shared definitions, production generation, optional benchmark searches, compact-data export and validation. |
| `Plotting/` | The 60-figure production exporter and standalone spectrum, geometry, clockwork and cancellation plots. |
| `Analysis/` | Mixing-angle widths and one-generation participation and cancellation diagnostics. |
| `data/` | Location for downloaded compact records and the separate diagnostic WXF file. See its README for loading instructions. |

## Data and plotting

The large compact dataset is distributed as a ZIP attachment through the repository's **Releases** section. Download the current release's compact-data ZIP and extract it locally. The dataset contains 36 cases with 10,000 trial records per case. Keep the extracted filenames unchanged. The earlier `CompactPlotData_040fd934` export is incomplete for this workflow and should not be used.

Open `Plotting/Mass_Mix_Plots.nb` and evaluate its input cell. Select the extracted compact-data folder. No production simulation or benchmark search is needed. The exporter checks the fixed paper benchmarks before writing 60 named PDFs to a new folder beside the notebook.

For saved-data analysis, open `Analysis/Analysis.nb`. Its angle-width section reads the same compact folder. The participation and cancellation section instead reads `data/TowerDiagnostics_20261019_164715.wxf`. If that file is supplied as a separate release attachment, download it into `data/` first.

## Generation and standalone calculations

Keep the repository's folder names and relative positions. Open saved notebooks in Mathematica and use a fresh kernel when switching workflows, particularly between standalone notebooks that reuse global variables.

To generate data, evaluate the initialization cell in `generation/DefinitionsOnly.nb`, followed by the desired section. Avoid evaluating the entire notebook at once, since its search and production sections launch substantial calculations. The generation-folder README gives the section map.

The standalone notebooks in `Plotting/` contain their own matrix definitions. The cancellation plot evaluates the rationalised matrices at 100-digit precision and uses 80-digit working precision for plotting.

`RELEASE_CHECKS.md` records the completed data and plotting checks and the publication steps. `WeakSiteBenchmarks.wxf` and `Production_Parameters.wl` contain inputs rather than production samples. Full original data are retained separately from the compact release.
