# Latent Time & Terminal States

Identify candidate root cells and end points from the velocity graph, then integrate gene-wise fitted times into a shared latent-time axis.

## Before you begin

Start with the `pv` object from [RNA Velocity Analysis](rna-velocity-analysis.md), after completing dynamics recovery, velocity computation and velocity graph construction.

The examples use an existing `umap` embedding. Cell colours display the inferred scores or latent time, with velocity arrows hidden.

## Infer root cells and end points

Estimate root-cell and end-point scores from the velocity transition graph:

```r
pv <- compute_terminal_states(pv)
```

The results are stored in the cell metadata:

- `root_cells` — relative support for cells being near the beginning of the inferred dynamics.
- `end_points` — relative support for cells being near terminal regions of the inferred dynamics.

Finite scores range from 0 to 1 and represent relative spectral support, not calibrated probabilities. Isolated, unselected or constant-score cells receive NA. Forward and backward kernels match those used for fate and origin. Recomputing terminal states clears old latent-time, fate/origin and related driver caches.

### Root cells

Visualize the root-cell scores:

```r
plot_velocity_embedding(pv, reduction = "umap", group_by = "root_cells", show_arrows = FALSE)
```

```{figure} figures/root_cells.png
:figclass: figure-compact
:alt: UMAP embedding coloured by inferred root-cell scores.

Root-cell scores across the UMAP embedding.
```

Higher scores highlight candidate starting regions in the velocity graph.

### End points

Visualize the end-point scores:

```r
plot_velocity_embedding(pv, reduction = "umap", group_by = "end_points", show_arrows = FALSE)
```

```{figure} figures/end_points.png
:figclass: figure-compact
:alt: UMAP embedding coloured by inferred end-point scores.

End-point scores across the UMAP embedding.
```

Higher scores highlight candidate terminal regions. Interpret these regions alongside cell annotations and the velocity field.

## Compute latent time

Combine gene-wise fitted times into a shared temporal ordering:

```r
pv <- compute_latent_time(pv, min_likelihood = 0.1)
```

The function uses fitted gene times and automatically searches for root-cell information, including the `root_cells` scores computed above. The resulting values are stored in `pv@meta.data$latent_time`.

- `min_likelihood` controls the fitted-likelihood threshold used to select genes for time integration. At least two genes must have valid fits and sufficient time coverage; there is no fallback that reinstates rejected genes. Unsupported cells remain NA.
- `min_confidence` controls the local-consistency threshold used during neighbourhood-based correction. Its default value is `0.75`.
- By default, latent time is normalized to the interval from 0 to 1.

### Visualize latent time

```r
plot_velocity_embedding(pv, reduction = "umap", group_by = "latent_time", show_arrows = FALSE)
```

```{figure} figures/latent_time.png
:figclass: figure-compact
:alt: UMAP embedding coloured by inferred latent time.

Shared latent time across the UMAP embedding.
```

Lower values represent earlier positions and higher values represent later positions along the inferred temporal ordering. Latent time is a relative model-derived coordinate, not elapsed time in hours or days.

### Specify root and end-point information

To explicitly use the inferred root scores and add end-point scores as a weak temporal prior:

```r
pv <- compute_latent_time(pv, min_likelihood = 0.1, root_key = "root_cells", end_key = "end_points")
```

- `root_key` specifies the metadata column used to identify candidate root cells.
- `end_key` specifies the metadata column used as a weak end-point prior. This constraint is optional and is not applied by the basic call above.

This call updates the latent-time values in `pv`.

## Assign fate and origin with explicit targets

Select terminal or root cells from biological annotations and supply their names or a complete logical mask as `target_cells`. If you use spectral scores to guide selection, supply the resulting targets explicitly.

```r
# terminal_ids and root_ids are user-selected cell-name vectors.
pv <- compute_cell_fate(pv, groupby = "clusters", target_cells = terminal_ids)
pv <- compute_cell_origin(pv, groupby = "clusters", target_cells = root_ids)
```

Targets become absorbing boundaries only for this calculation. Target-cell probabilities are summed within each target group before assignment. Confidence is the largest target-group absorption probability; it is not a statistical confidence interval.

Unreachable or excluded cells have NA labels and confidence. Tied groups have an NA label and retain their maximum probability. Non-target closed classes and excluded groups consume probability mass, which is not renormalized away. The backward kernel describes reverse-graph origin support, not calibrated biological ancestry.

No complete probability matrix or additional QC fields are stored. Sparse solver failures and non-convergence error without replacing existing results. Rebuilding neighbors or moments requires recomputing invalidated downstream results.

## Inspect the results

Review the inferred values for individual cells:

```r
head(pv@meta.data[, c("root_cells", "end_points", "latent_time")])
```

Compare their distributions:

```r
summary(pv@meta.data[, c("root_cells", "end_points", "latent_time")])
```

Assess whether the inferred starting regions, terminal regions and temporal ordering agree with the cell annotations, velocity directions and gene phase portraits. These results describe inferred transcriptional dynamics and do not directly establish cell ancestry.

## Continue to driver-gene ranking

Continue with [Driver Gene Ranking](driver-gene-ranking.md) to explore genes associated with the inferred dynamics.
