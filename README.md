# Visium Spatial Transcriptomics Workshop — Developing Human Heart

Hands-on introduction to the analysis of 10x Genomics Visium data in R, using
sections of the developing human heart.

Written for the **2026 V4SDB Student Summer School**. No prior coding experience
is assumed — all code is provided, and the aim is to read and understand each
step rather than write it from scratch.

---

## What you will do

Working from raw Space Ranger output, you will:

1. Load a Visium section and plot genes directly on the tissue
2. Decide which genes to filter out, and justify why
3. Normalize, scale, and reduce the data
4. Cluster the spots and map the clusters back onto the anatomy
5. Identify the genes that define each cluster
6. Compare your clusters against a published atlas
7. Repeat the analysis on three sections of different ages, with batch correction
8. Explore age-associated expression and spatial autocorrelation

Approximate duration: 4 hours.

---

## Getting started

Everything runs in a pre-configured environment on the workshop server —
R, all packages, and the data are already installed.

1. Open RStudio Server and navigate to the project folder
2. Open `visium-heart-walkthrough.qmd`
3. Run chunks in order with `Cmd/Ctrl + Shift + Enter`

A rendered HTML version is included for reference and works offline.

### Running it yourself

If you want to work through the material after the workshop, you will need
R (≥ 4.2) and:

```r
install.packages(c("tidyverse", "Seurat", "patchwork", "RColorBrewer", "harmony"))
remotes::install_github("ludvigla/semla")
```

Then edit the `.libPaths()` call in the setup chunk (or delete it) and place the
data as described below.

---

## Repository contents

```
.
├── visium-heart-walkthrough.qmd    # the workshop
├── visium-heart-walkthrough.html   # rendered version
├── data/
│   ├── spaceranger/                # one folder per section
│   │   └── <slide_id>/
│   │       ├── filtered_feature_bc_matrix.h5
│   │       └── spatial/
│   ├── metadata.csv                # sample_id, organ, organism, pcw_age, sex
│   ├── GRCh38-2020-A_gene_biotypes.rds
│   └── img/
└── README.md
```

---

## The data

Three Visium sections of the developing human heart at post-conceptional
weeks 6, 9, and 10 — a small teaching subset of a published atlas:

> Lázár, E.\*, Mauron, R.\*, Andrusivová, Ž.\* et al.
> **Spatiotemporal gene expression and cellular dynamics of the developing human heart.**
> *Nature Genetics* **57**, 2756–2771 (2025).
> [doi:10.1038/s41588-025-02352-6](https://doi.org/10.1038/s41588-025-02352-6)
> (\*equal contribution)

The full study combined Visium on 16 hearts (38 sections, 69,114 spots, PCW 6–12)
with single-cell RNA-seq on 15 hearts (76,991 cells, PCW 5.5–14) and in situ
sequencing of 150 transcripts, resolving 23 spatial clusters and mapping 72
cell states onto tissue positions.

**Explore the full atlas:** <https://hdcaheart.serve.scilifelab.se>
**Full analysis code:** <https://github.com/rmauron/HDCA_heart_dev>
**Processed data:** [Mendeley 1](https://doi.org/10.17632/fhtb99mdzd.1) ·
[Mendeley 2](https://doi.org/10.17632/w65jtfsvpr.1)

Tissue was donated after elective medical abortion with written informed consent,
under Swedish ethical permit 2018/769-31.

---

## Where to go next

- [semla documentation](https://spatial-research.github.io/semla/index.html)
- [Seurat spatial vignette](https://satijalab.org/seurat/articles/spatial_vignette)
- [Deconvolution with Stereoscope](https://www.nature.com/articles/s42003-020-01247-y)
- [Spatially aware clustering with BANKSY](https://github.com/prabhakarlab/Banksy)
- [Spatial interactions with SpatialDM](https://github.com/StatBiomed/SpatialDM)

---

## Contact

Raphaël Mauron — [@rmauron](https://github.com/rmauron)

Questions about the material are welcome, during the workshop or after.