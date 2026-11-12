# plantvelo

Prepare IR-aware count layers for plant RNA velocity analysis with the upstream plantvelo Python module.

[GitHub repository](https://github.com/PlantVelocityProject/plantvelo)

## How plantvelo fits into the workflow

plantvelo performs prior-guided molecular classification upstream of the PlantVelocity R package. It uses intron-retention priors to distinguish IR-associated signals when preparing spliced and unspliced count layers for downstream analysis.

The resulting loom file connects the upstream processing workflow to PlantVelocity, where count layers are combined with Seurat cell annotations and dimensional reductions.

## Prepare the inputs

Prepare the following resources for upstream processing, checking the required formats against the repository instructions:

- Read alignments and cell barcodes for the single-cell dataset.
- Gene and repeat annotations compatible with the reference genome used for alignment.
- An IR prior matching the species and reference annotation.

Keep reference versions and gene identifiers consistent across these inputs and the Seurat object used downstream.

## Use the loom output

PlantVelocity reads the following count layers from the loom file:

- `layers/spliced` — spliced counts for velocity analysis.
- `layers/unspliced` — unspliced counts for velocity analysis.
- `layers/retained` — optional IR-associated counts, imported as the separate `ir` layer.

Cell identifiers and the selected gene identifier field must match the Seurat object. Choose `loom_gene_attr = "Accession"` for Gene_ID or `loom_gene_attr = "Gene"` for Gene_symbol.

IR-associated counts remain separate from the two-state kinetic model and can be explored through [Intron Retention Analysis](intron-retention-analysis.md).

## Continue to Data Preparation

Once the loom file is ready, continue with [Data Preparation](data-preparation.md) to align the inputs and create a `plantvelo` object in R.
