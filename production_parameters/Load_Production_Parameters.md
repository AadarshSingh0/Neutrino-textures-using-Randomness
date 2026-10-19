# Load the production parameters

Place `Production_Parameters.wl` in the `production_parameters/` folder.
Run the following cells from a saved Mathematica notebook in the repository root.

## Load all benchmark parameter sets

```mathematica
productionParameters = Get[
  FileNameJoin[{NotebookDirectory[],
    "production_parameters", "Production_Parameters.wl"}]
];

Length[productionParameters]
(* 36 *)
```

## List the available benchmark names

```mathematica
Keys[productionParameters]
```

## Select one benchmark

```mathematica
p = productionParameters[
  "petersen_site_weak_hamiltonian_majorana"
];
```

## Read individual parameters

```mathematica
Keys[p]

p["Sites"]
p["HoppingMean"]
p["HamiltonianFlavor"]
```

These commands load saved parameters. They do not generate simulation data.
