# PlantIRPrior

PlantIRPrior provides species-specific, high-confidence intron-retention (IR) priors derived from long-read transcriptomic evidence for PlantVelocity workflows. The resource covers *Arabidopsis thaliana*, rice, tomato, maize and soybean.

[Visit the PlantIRPrior website](https://plantvelocityproject.github.io/PlantIRprior/)

## Explore the resource

The website brings together:

- **IR event browsing** by species, gene or genomic interval.
- **Source samples and reference information**, including the genome assemblies and gene annotations used to construct each prior.
- **Versioned downloads** of high-confidence IR events, PlantVelo priors, sample metadata and supporting evidence.
- **Methods** describing prior construction and the agent-assisted workflow for reviewed resource updates.

## Use with PlantVelocity

Choose a release matching the species, reference genome assembly and annotation used for your single-cell analysis. Download its **PlantVelo prior** TSV and use it with the upstream [plantvelo Python module](plantvelo.md) to prepare IR-aware count layers. Keep the release version with your analysis records so that the prior can be identified and reused.

Continue with [Data Preparation](data-preparation.md) to import the resulting loom file into PlantVelocity.
