# Data generation and validation

| File | Purpose |
| --- | --- |
| `TheorySpaceNeutrinos.wl` | Shared matrix construction, seeded generation, mass and mixing calculations, and precision helpers. |
| `DefinitionsOnly.nb` | Production generation, optional benchmark finding, one-generation diagnostics and compact-data export. Despite its name, evaluating its run cells performs calculations. |
| `WeakSiteBenchmarks.wxf` | Saved weak-site benchmark settings from the tuning workflow. Final production uses `../production_parameters/Production_Parameters.wl`. |
| `Validation.nb` | Usage examples and checks of the definitions, reference cases, precision and parallel execution. |
| `UnifiedChecks.wl` | Validation functions used by `Validation.nb`. |
| `UnifiedReferenceData.wl` | Stored inputs and numerical references required by `UnifiedChecks.wl`. |

## Select a workflow

In `DefinitionsOnly.nb`, evaluate the first initialization cell, then the required section. Input-cell numbers below count Input cells only, starting at one.

| Task | Input cells | Output beside the notebook |
| --- | --- | --- |
| Strong-site Dirac | Set `fig15Parameters` in 3 or 4, then run 6 | `PaperPlots_<timestamp>/`, with 6 `*_data.wxf` files and a summary. |
| Strong-site Majorana | 7, then 8 | `MajoranaPlots_<timestamp>/`, with 6 case folders and a summary. |
| Find new weak-site benchmarks | 9–14 in order | Angle-search and mass-tuning folders, then `WeakSiteBenchmarks.wxf`. Optional and separate from frozen production. |
| Final weak-site production | 15 | `WeakSite_Final_<timestamp>/`, with 12 cases, settings and a summary. |
| Final weak-hopping production | 16 | `WeakHopping_Final_<timestamp>/`, with 12 cases, settings and a summary. |
| One-generation diagnostics | 17 | `TowerDiagnostics_<timestamp>.wxf`, containing configuration, parameters and samples. |
| Compact existing production data | 18 | A new `CompactPlotData_<id>/` folder and ZIP beside the selected source folder. |

Production uses 10,000 trials per case, saved before downstream selection. The diagnostic calculation is separate and includes strong hopping as well as the three production regimes.

For compact export, select a common parent containing exactly one original production version per case. Do not include search runs or duplicate generations. This cell reads existing WXF files and does not rerun simulations. It retains sample order, failures, light masses, mixing information and numerical diagnostics while removing the large matrix and full-spectrum fields. Keep the original full records privately.

For a short validation run, use `Validation.nb` input 1, then 18–20. Input 21 checks all stored references. Input 26 checks serial/parallel agreement after loading the validation functions. These checks do not launch the 36 production runs.
