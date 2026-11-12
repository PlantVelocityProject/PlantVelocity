# Data Preparation

Prepare a Seurat object and the corresponding loom file before starting RNA velocity analysis. PlantVelocity combines count layers from the loom file with cell annotations and dimensional reductions from Seurat.

## Input data

You need two inputs:

- A Seurat v5 object containing the cells to analyse, an RNA assay and a PCA reduction. Include a UMAP reduction to follow the visualization example in Quick Start.
- A corresponding loom file containing spliced and unspliced count layers.

For IR-aware analysis, generate the loom file using the upstream `plantvelo` workflow with a species-specific IR prior, read alignments, cell barcodes, and gene and repeat annotations. Prior-based count correction takes place upstream, before importing the loom file into R.

## Prepare the Seurat object

Complete quality control and appropriate preprocessing in Seurat before import. Retain the RNA assay, cell annotations and dimensional reductions needed for downstream analysis.

```r
library(Seurat)
library(PlantVelocity)

seob <- readRDS("seurat_object.rds")
```

Inspect the available assays, reductions and metadata:

```r
Assays(seob)
Reductions(seob)
head(seob[[]])
```

The examples use the `RNA` assay, `pca` and `umap` reductions, and a metadata column named `cell_type` for colouring cells. Replace these names with those used in your own object where applicable. PCA is required for neighbour graph construction.

## Prepare the loom file

The loom file must contain:

- `layers/spliced` — spliced counts.
- `layers/unspliced` — unspliced counts.
- `col_attrs/CellID` — cell identifiers.
- A gene identifier field in `row_attrs`: `Gene` for Gene_symbol or `Accession` for Gene_ID, selected using `loom_gene_attr`.

If `layers/retained` is present, PlantVelocity imports it as the `ir` layer for independent intron-retention analysis. This layer does not enter two-state kinetic fitting.

## Match cell and gene identifiers

Cell names must match exactly between the Seurat object and `col_attrs/CellID`, including sample prefixes and barcode suffixes. The current importer requires cell names in both inputs to follow `{sample_name}_barcode-1`, for example `sample1_AAACCCAAGAAACACT-1`.

```r
head(colnames(seob))
head(rownames(seob))
```

Select the gene identifier field that matches the Seurat feature names:

- `loom_gene_attr = "Accession"` — construct the object using Gene_ID. Seurat feature names must match `row_attrs/Accession`.
- `loom_gene_attr = "Gene"` — construct the object using Gene_symbol. Seurat feature names must match `row_attrs/Gene`. This is the default.

This parameter selects the matching field; it does not convert gene IDs to gene symbols or vice versa. All Seurat cell and feature names must be present in the corresponding loom fields. The loom file may contain additional cells and genes; imported count layers are selected and reordered to match the Seurat object.

## Create the PlantVelocity object

For a Seurat object whose feature names are Gene_symbol, use `"Gene"`:

```r
pv <- create_plantvelo(seurat_object = seob, assays = "RNA", loom_path = "plantvelo_output.loom", loom_gene_attr = "Gene")
```

For a Seurat object whose feature names are Gene_ID, use `"Accession"` instead:

```r
pv <- create_plantvelo(seurat_object = seob, assays = "RNA", loom_path = "plantvelo_output.loom", loom_gene_attr = "Accession")
```

To set the active cell identities from a metadata column, also supply `input_ident`. For example, when using Gene_symbol:

```r
pv <- create_plantvelo(seurat_object = seob, assays = "RNA", loom_path = "plantvelo_output.loom", loom_gene_attr = "Gene", input_ident = "cell_type")
```

Replace `"cell_type"` with an existing metadata column. If `input_ident` is omitted, the existing Seurat identities are retained.

## Continue to RNA velocity analysis

The resulting `plantvelo` object contains aligned count layers, cell metadata and dimensional reductions. Continue with [RNA Velocity Analysis](rna-velocity-analysis.md) to construct the neighbour graph and fit RNA velocity.
