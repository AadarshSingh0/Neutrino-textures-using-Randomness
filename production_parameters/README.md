# Exact production parameters

- `Production_Parameters.wl` returns an Association containing all 36 benchmark parameter associations, with their original numerical representations and 30 fields per benchmark.
- `Load_Production_Parameters.md` shows how to inspect the association from a notebook saved at the repository root.

The cases cover 6 strong-site Dirac benchmarks, 6 strong-site Majorana benchmarks, 12 weak-site benchmarks and 12 weak-hopping benchmarks. These files contain inputs, not simulated samples.

From a notebook in `generation/`, `Plotting/` or `Analysis/`, load them with:

```mathematica
productionParameters = Get[
  FileNameJoin[{ParentDirectory[NotebookDirectory[]],
    "production_parameters", "Production_Parameters.wl"}]
];
Length[productionParameters]  (* 36 *)
Keys[productionParameters]
```

Use `productionParameters[caseName]` to select a case. The final weak-site and weak-hopping generation sections already load this file through the same relative path. They use the frozen production inputs independently of any new benchmark search.
