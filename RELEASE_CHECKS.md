# Release checks

## Recorded checks on 19 October 2026

| Item | Result |
| --- | --- |
| Exact production parameters | All 36 parameter associations are present. |
| Compact-export report | All 36 expected cases are present, with no duplicates or missing cases and 10,000 records per case. Total: 360,000 trial records. |
| Data size | Original case files total 4.22 GB. Compact case files total 1.85 GB before ZIP compression. The resulting ZIP was reported as approximately 540 MB. |
| Production plotting | The maintainer ran `Mass_Mix_Plots.nb` on the new compact data and reported successful generation of all 60 PDFs. |
| Diagnostic data | `TowerDiagnostics_20261019_164715.wxf` was generated and placed in `data/`. The maintainer loaded it through `Analysis.nb` and supplied the participation summary, with 1,000 resolved samples in each of the 12 cases. |
| Cancellation precision | The maintainer replaced the three machine-precision matrix conversions with `N[..., 100]`, added `WorkingPrecision -> 80` to the cancellation plot and confirmed that it runs. Save these edits in `Plotting/Spectrum_Cancellation_plot.nb`. |
| Validation dependency | `generation/UnifiedReferenceData.wl` is present beside `UnifiedChecks.wl` and `TheorySpaceNeutrinos.wl`. |

These checks establish that the updated production plotting and diagnostic-loading workflows ran on the maintainer's machine. They do not establish fresh execution of every validation test, benchmark search or standalone plotting notebook, or exact agreement of every legacy figure with the manuscript.

## Files to retain

- The notebooks, shared `.wl` files, exact production parameters and folder READMEs.
- The new compact-data ZIP and its `export_report.csv`.
- The generated `plot_export_report.wxf` and `expected_filenames.txt` from the successful 60-plot export.
- `data/TowerDiagnostics_20261019_164715.wxf` for the separate diagnostic analysis.
- The original full production data and their settings, archived on the external drive. Compact files omit full heavy spectra and large matrices and are not complete backups.

## Publication

1. Upload the code, parameters and documentation to the repository, including the saved precision edits.
2. Attach the new compact-data ZIP to a GitHub Release. Distribute the diagnostic WXF as a release attachment too if it is too large for the repository.
3. Identify the compact-data archive and diagnostic file clearly in the release description. The root and data READMEs direct users to the repository's Releases section.

No release URL or tag has been recorded here because publication has not yet been confirmed. Keep any later version-specific download links consistent with the actual uploaded assets.

## Paths and loading

Generation loads shared definitions from its own folder. Final weak-site and weak-hopping production loads `../production_parameters/Production_Parameters.wl`. New simulation outputs are created in timestamped folders beside `DefinitionsOnly.nb`.

Plotting and analysis select data through dialogs and discover their current paths. Select exactly one version of the 36-case dataset. Do not select a directory containing duplicate original and compact copies. The diagnostic loader selects its individual WXF file separately.

The four retained standalone plotting notebooks define their own matrices. No external `Get`, `Needs` or `Import` dependency on the deleted legacy modules was found in the inspected input cells.
