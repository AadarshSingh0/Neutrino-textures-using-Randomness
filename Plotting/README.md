# Plotting notebooks

| File | Purpose |
| --- | --- |
| `Mass_Mix_Plots.nb` | Recreates 60 production PDFs from all 36 original or revised compact datasets, using the fixed paper filenames and angle conventions. |
| `Spectrum_Cancellation_plot.nb` | Standalone three-geometry mass spectrum and shifted-inverse cancellation plots extracted from `Site.nb`, with the high-precision cancellation update. |
| `Plots.nb` | Earlier single-generation spectrum, mode-profile and disorder plots. Settings and matrix construction are inside the notebook. |
| `l_vs_nl_vs_pet.nb` | Earlier local, non-local and Petersen-type geometry comparisons, including random ensembles and mode profiles. |
| `clockwork_vs_random.nb` | Earlier clockwork-versus-randomness comparisons and related parameter examples. |

## Production figures

Download the current compact-data ZIP from the repository's Releases section and extract it. Open the saved `Mass_Mix_Plots.nb`, evaluate its input cell, and select that data folder. The notebook discovers files recursively. Do not select a parent containing both original and compact copies or several runs of the same case.

The exporter checks parameters, record counts, selection and angle-bin counts against its embedded paper reference. It is intended for these fixed benchmarks. A changed benchmark should fail those checks rather than silently overwrite a paper plot.

Successful export creates `PaperPlots_<timestamp>_<id>/figures/` beside the notebook, plus `plot_export_report.wxf` and `expected_filenames.txt` in its parent output folder. No separate exporter `.wl` is needed. The maintainer successfully generated all 60 PDFs from the revised compact dataset.

## Standalone plots

The other four notebooks define their own matrices and do not import the deleted legacy `.m` modules. Run each in a fresh kernel and evaluate the setup for the desired calculation. Some sections generate large random ensembles. They are separate from the 36-case compact-data workflow.

In `Spectrum_Cancellation_plot.nb`, the three rationalised matrices are evaluated with `N[..., 100]`, and the cancellation `Plot` uses `WorkingPrecision -> 80`. The spectral denominators retain signed eigenvalues. The precision update leaves the plot styling unchanged, and the maintainer confirmed that the updated plot runs.

Re-evaluating matrix generation draws new random realisations. Matching a particular earlier curve requires its original realisation and settings. The recent release checks do not establish exact reproduction of every manuscript figure by the standalone notebooks.
