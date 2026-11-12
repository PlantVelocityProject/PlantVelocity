# About PlantVelocity

PlantVelocity is a framework for inferring potential directions of cell-state change from plant single-cell transcriptomic data. It uses species-specific long-read intron-retention (IR) priors to identify and exclude potential IR-associated signals before fitting a conventional two-state RNA velocity model.

RNA velocity estimates the rate of gene expression change from the relative abundance of unspliced and spliced RNA. In plants, widespread intron retention complicates the interpretation of intronic reads: these may originate from nascent RNA undergoing splicing or from transcripts with retained introns. Short-read, 3′-end sequencing often lacks the surrounding splicing context needed to distinguish these sources.

Long-read transcript structures provide independent evidence about intron retention and its neighbouring splicing context. PlantVelocity uses this evidence to guide molecular classification and correct input counts before neighbourhood smoothing and kinetic modelling, helping align the observed signals with the RNA states represented by the model.

## How PlantVelocity works

The workflow connects prior-guided molecular classification in Python with RNA velocity modelling and downstream analysis in R:

1. **Classify molecules.** The upstream `plantvelo` Python module combines Cell Ranger alignments, cell barcodes, gene and repeat annotations, and a species-specific IR prior. Read evidence is aggregated by cell barcode, UMI and gene to assign each molecule to one of four mutually exclusive categories: spliced, unspliced, retained or ambiguous.
2. **Correct the inputs.** Retained and ambiguous molecules are excluded from kinetic inference. The remaining spliced and unspliced counts are exported in loom format. Retained counts are preserved separately for quality control and descriptive analyses.
3. **Infer cell-state dynamics.** The PlantVelocity R package integrates the corrected counts with a Seurat v5 object, retaining its cell annotations and low-dimensional embeddings. It constructs a neighbour graph in PCA space, smooths the spliced and unspliced counts, and fits a two-state model of transcription, splicing and degradation. Gene-wise velocities are then used to construct a directed velocity graph and estimate cell-to-cell transition probabilities.

Outputs include velocity embeddings and stream plots, latent time, root and endpoint probabilities, gene phase portraits, and candidate driver-gene rankings. The retained layer supports independent summaries and visualisations of IR-associated signals across cells and genes.

## IR priors across five plant species

The study integrated **683 public long-read transcriptomic samples** spanning Iso-Seq, FLEP-seq and Oxford Nanopore sequencing to construct species-specific priors containing **8,935 high-confidence IR events** in total. The resource covers *Arabidopsis thaliana*, rice (*Oryza sativa*), tomato (*Solanum lycopersicum*), maize (*Zea mays*) and soybean (*Glycine max*). Each prior was built independently using the corresponding reference genome and annotation, combining retained and spliced read evidence with support from neighbouring intron splicing patterns.

The study also outlines an AI agent-assisted workflow for incremental resource updates, linking data discovery, read acquisition and alignment, IR evidence aggregation, quality control and versioned prior releases. Human review precedes evidence aggregation and prior updates, with provenance and version records supporting reproducible analyses.

Prior coverage depends on the tissues, conditions and sequencing data represented in the resource. An event supported by external long-read data need not occur in every cell state of a target dataset. Analyses require a prior matched to the species and reference annotation.

## Biological applications

The study evaluated PlantVelocity in Arabidopsis root hormone-response and developmental datasets, and applied it to a root-tip regeneration time series:

- **Brassinolide response.** Excluding potential IR-associated signals improved agreement between predicted directions and a cell-wall remodelling response programme in root cortex cells. The differences were most pronounced during the early response and near the boundary between non-responsive and responsive states.
- **Normal root development.** Input correction strengthened predicted transitions from meristematic to elongating states in cortex and endodermal cells, improving alignment with the known developmental sequence.
- **Root-tip regeneration.** Predicted transitions from intermediate to younger root-cap states differed most between excised roots and uncut controls at 24 h after excision. This coincided with increased expression of *PHB3* and *TCTP1*, genes associated with growth and development, in the intermediate group.

In both the brassinolide-response and normal-development evaluations, directionality metrics exceeded those obtained in each of **200 random-exclusion controls** matched to the number of unspliced counts removed per cell. These comparisons support a contribution from the composition of the excluded signals to the observed changes in directionality.

## Getting started

Start with [Installation](installation.md), then follow [Data Preparation](data-preparation.md) to combine a Seurat object with count layers generated by the upstream [plantvelo Python module](plantvelo.md). The [Quick Start](quick-start.md) walks through the core R workflow, from count import and kinetic fitting to velocity visualisation.

```{toctree}
:caption: Overview
:maxdepth: 2
:hidden:

installation
quick-start
release-notes
reference
```

```{toctree}
:caption: User Guide
:maxdepth: 2
:hidden:

data-preparation
rna-velocity-analysis
latent-time-terminal-states
driver-gene-ranking
intron-retention-analysis
```

```{toctree}
:caption: Resources
:maxdepth: 2
:hidden:

plantvelo
```
