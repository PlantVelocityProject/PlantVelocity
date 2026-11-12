# Installation

Get started with PlantVelocity for RNA velocity analysis in plant single-cell data.

## Requirements

PlantVelocity declares **R 4.1.0 or later**. Use a recent R release compatible with the versions of Matrix, SeuratObject, and Seurat you install; these dependencies may require a newer R version.

Run the R commands below in an R console or RStudio. The package itself does not require compilation, although some dependencies contain compiled code.

## Install from GitHub

Use `remotes` to install PlantVelocity. [GitHub repository](https://github.com/PlantVelocityProject/PlantVelocity)

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("PlantVelocityProject/PlantVelocity")
```

## Dependencies

Required dependencies are declared in the package [DESCRIPTION](https://github.com/PlantVelocityProject/PlantVelocity/blob/HEAD/DESCRIPTION) file:

- `R6` (≥ 2.5.0) — object infrastructure.
- `Matrix` (≥ 1.5.0) — sparse matrix operations.
- `SeuratObject` (≥ 5.0.0) — Seurat data structures.
- `Seurat` (≥ 5.0.0) — creating PlantVelocity objects from Seurat data.
- `hdf5r` (≥ 1.3.0) — HDF5 / loom file access.
- `ggplot2` (≥ 3.4.0), `scales` — visualization.
- `ggrastr` — rasterized points and arrows in velocity plots.
- `RSpectra` — sparse eigendecomposition for terminal-state inference.
- `patchwork` — combined diagnostic and IR plot layouts.
- `methods`, `stats`, `utils`, `grid` — included with R.

### Optional features

- `BiocNeighbors` — neighbor search with `method = "bioc"`.
- `ggrepel` — repelled labels in IR plots.
- `testthat` (≥ 3.0.0) — package development and testing.
