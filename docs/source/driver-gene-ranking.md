# Driver Gene Ranking

Prioritize candidate genes whose RNA velocities are consistent with expression changes implied by the velocity graph.

## Before you begin

Start with the `pv` object from [RNA Velocity Analysis](rna-velocity-analysis.md), after computing expression moments, fitted dynamics, RNA velocity and the velocity graph.

Driver-gene ranking uses fitted likelihoods and velocity coherence. Latent time and terminal-state inference are not required.

## Rank candidate driver genes

Compute the ranking:

```r
pv <- rank_velocity_genes(pv)
```

The function first filters genes by fitted likelihood, then ranks the retained genes by decreasing velocity coherence.

Coherence measures the correlation across cells between a gene’s estimated velocity and its expected expression change under the velocity transition matrix. Higher scores indicate stronger agreement with the inferred dynamics.

The default settings are:

- `likelihood_quantile = 0.3` — retain approximately the top 30% of eligible fitted genes by likelihood; ties at the threshold may retain additional genes.
- `min_likelihood = 0` — apply no additional positive absolute likelihood threshold.
- `coherence_mode = "s_only"` — rank genes using spliced velocity coherence.
- `n_top = NULL` — retain all genes passing the filters and coherence checks.

Genes with missing coherence scores or no non-zero spliced velocities are excluded from the resulting table.

## Inspect the ranked genes

Retrieve the results:

```r
driver_table <- get_velocity_gene_stats(pv)
head(driver_table, 10)
```

A single all-gene table is stored in `pv@misc$velocity_gene_stats`. The getter returns selected genes in rank order by default; use `get_velocity_gene_stats(pv, include_filtered = TRUE)` for every gene. The table contains:

- `gene` — the gene identifier used in the object.
- `coherence` — the score used for ranking.
- `coherence_s` — coherence calculated from spliced velocity and expression.
- `coherence_u` — coherence calculated from unspliced velocity and expression.
- `mean_abs_velocity` — mean absolute spliced velocity across cells.
- `n_active` — number of cells with non-zero, non-missing spliced velocity.
- `fit_likelihood` — fitted likelihood from dynamics recovery.
- Raw U/S detection counts and totals, `pass_fit`, `pass_velocity`, `pass_raw`, `pass_likelihood`, `rank`, and `selected`.

Ranking settings are stored in `pv@misc$velocity_params$ranking`; velocity calculation settings stay at the top level. No separate QC table, `qc_version`, or `exclusion_reason` is stored. Rerun `rank_velocity_genes()` to update older saved objects.

Gene identifiers follow the naming scheme selected during object construction: Gene_ID for `loom_gene_attr = "Accession"` or Gene_symbol for `loom_gene_attr = "Gene"`.

## Visualize the top genes

Display the ten highest-ranked genes:

```r
plot_driver_genes(pv, n_top = 10)
```

Choose the ordering before taking the displayed top genes:

```r
plot_driver_genes(pv, n_top = 10, sort_by = "coherence")
plot_driver_genes(pv, n_top = 10, sort_by = "fit_likelihood")
```

Both modes use already selected drivers, in descending score order. Ties use gene names; missing and non-finite sorting scores come last. The x axis remains coherence and point colour remains fit likelihood. Stored ranks and selection flags are unchanged.

```{figure} figures/driver_genes.png
:figclass: figure-compact
:alt: Lollipop chart of the top ten driver genes ranked by coherence, with point colour indicating fitted likelihood.

Example ranking of the top ten candidate driver genes.
```

- The horizontal position represents velocity coherence.
- Genes are ordered by decreasing coherence.
- Point colour represents fitted likelihood.

Here, `n_top` controls how many genes are displayed without changing the stored ranking.

## Adjust the ranking

To rank genes using the mean of spliced and unspliced coherence:

```r
pv <- rank_velocity_genes(pv, coherence_mode = "mean")
```

To apply an absolute likelihood threshold and retain at most 50 ranked genes:

```r
pv <- rank_velocity_genes(pv, min_likelihood = 0.1, n_top = 50)
```

The default top-30% likelihood filter still applies. Set `likelihood_quantile = 1` to disable that relative filter:

```r
pv <- rank_velocity_genes(pv, min_likelihood = 0.1, likelihood_quantile = 1)
```

Each call replaces the ranking stored in `pv`. Retrieve the table again after changing the settings.

High coherence identifies genes consistent with the inferred velocity field. It does not by itself establish that a gene causally controls a developmental transition. Examine candidate genes alongside their phase portraits, expression patterns and biological annotations.

## Continue to intron-retention analysis

Continue with [Intron Retention Analysis](intron-retention-analysis.md) to explore IR-associated signals independently of the two-state velocity model.
