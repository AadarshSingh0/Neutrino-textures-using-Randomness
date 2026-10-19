# Saved numerical data

Download the current compact-data ZIP from this repository's **Releases** section and extract it locally, for example inside this folder. Keep the extracted filenames unchanged. Select the extracted compact-data folder when `Mass_Mix_Plots.nb` or the angle-analysis loader asks for a directory.

The production dataset contains:

| Regime | Cases | Trial records per case |
| --- | --- | --- |
| Strong-site Dirac | 6 | 10,000 |
| Strong-site Majorana | 6 | 10,000 |
| Weak site, both neutrino constructions | 12 | 10,000 |
| Weak hopping, both neutrino constructions | 12 | 10,000 |

There are 36 `*_plotdata.wxf` files and an `export_report.csv`. Samples retain their original order, including failures and records rejected by downstream selection. The compact format preserves light masses, mixing information, numerical diagnostics and per-sample parameter differences. Large matrices and full heavy spectra are omitted. The revised dataset was used to generate all 60 production PDFs.

Do not use the older `CompactPlotData_040fd934` export. Its strong-site files contain only parameters and its weak-regime files omit the PMNS matrices required by the current workflow.

## Separate diagnostic data

`TowerDiagnostics_20261019_164715.wxf` contains configuration, parameters and samples for the one-generation participation and cancellation diagnostics. If distributed as a release attachment, download it into this folder. In `Analysis.nb`, select this individual file in the diagnostic import dialog, then evaluate the two following analysis cells. The recorded run resolved 1,000 samples in each of 12 geometry/regime combinations.

The diagnostic file is separate from the 36-case production dataset. Neither `Production_Parameters.wl` nor `WeakSiteBenchmarks.wxf` replaces these sample files.
