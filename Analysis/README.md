# Analysis of saved data

`Analysis.nb` reads saved records without rerunning production simulations.

- **Quantify:** computes mixing-angle concentration using central 16th–84th percentile widths and the selection implemented in the notebook. Evaluate the loader first, then select one complete original or revised compact dataset. The weak-site and weak-hopping loader expects 12 cases of each. The following strong-site section should find another 12 files.
- **Participation and cancellation diagnostics:** reads `../data/TowerDiagnostics_20261019_164715.wxf` and computes the diagnostic and participation summaries. Select the file in the import dialog, then evaluate the two following analysis cells. This workflow was run successfully, with 1,000 resolved samples per case in the displayed 12-case participation summary.

Download the production ZIP from the repository's Releases section and extract it before selecting the data folder. If the diagnostic WXF is distributed separately, download it into `data/` first. It cannot be reconstructed from compact plotting records alone.

The notebook uses file dialogs and discovers current file paths. It does not depend on the original author's absolute directories. Select only one production version to avoid duplicate cases.

The loader supports both original per-sample parameter associations and the revised compact format with parameter overrides. No shared definition package is required for these saved-data calculations. Outputs appear in the notebook.
