# USUV_phylodynamics_Netherlands

This repository contains BEAST XML files and R scripts supporting the manuscript:

**Phylodynamic analyses reveal persistence and extensive intermixing of Usutu virus in the Netherlands**

Emmanuelle Münger, David F. Nieuwenhuijse, Nnomzie C. Atama, James Baxter, Lu Lu, Mandev S. Gill, Henk P. van der Jeugd, Judith M.A. van den Brand, Reina S. Sikkema, Marion P.G. Koopmans, Fabiana Gámbaro&dagger;, Bas B. Oude Munnink&dagger;, Simon Dellicour&dagger;

&dagger; equal contribution

## Repository contents

```
.
├── BEAST1_XML/         BEAST v1.10 XML files used for the phylodynamic analyses
├── DataVisualisation/  R scripts used to generate the manuscript figures
│   ├── DataVisualisationSharing.R   Main script generating Figs. 1-5
│   ├── ArboNL_CustomTheme.R         Shared ggplot2 themes used across figures
│   └── basemaps.R                   Shapefile loading and basemap plotting functions
├── LICENSE
└── README.md
```

Each script contains its own header with further details on dependencies, data requirements, and (where relevant) attribution for adapted code.

## Data availability

This repository contains **analysis and plotting code**. 

- **Sequences and raw reads generated in this study** have been deposited with metadata on the European Nucleotide Archive (ENA) under project **PRJEB83966** (accession numbers **ERR16767424–ERR16767471**). The full list of accession numbers is also available at [pathogensportal.nl/sequences.html](https://www.pathogensportal.nl/sequences.html).
- **Publicly available USUV genome sequences** included in our analyses, along with their accession numbers and associated metadata, are listed in **Supplementary Data 1** of the manuscript.
