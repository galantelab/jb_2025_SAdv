<br />
<p align="center">
  <h2 align="center">Impact of intermittent lead exposure on hominid brain evolution</h2>
  <p align="center">
    This repository contains the analysis pipeline (for bulk and single-cell RNA-Seq data) used for the paper titled "Impact of intermittent lead exposure on hominid brain evolution" published by JB et al in Science Advances.
    <br />
    <a href="https://github.com/galantelab/jb_2025_SAdv/issues">Issues</a>
  </p>
</p>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h3>Table of Contents</h3></summary>
  <ol>
    <li><a href="#overview">Overview</a></li>
    <li><a href="#requirements">Requirements</a></li>
    <li><a href="#repository-structure">Repository Structure</a></li>
    <li><a href="#data">Data</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#authors">Authors</a></li>
  </ol>
</details>

<!-- OVERVIEW -->
## Overview
<p align="justify"> This repository contains the scripts and data necessary for performing a quantification, integration and differential gene expression analysis using scRNA-Seq data and bulk RNA-Seq data. The analysis was conducted using different organoids (Cortical and Thalamic) from Homo Sapiens and Other hominids in different lead exposures (10uM and 30uM). </p>

<!-- REQUIREMENTS -->
## Requirements
The pipeline was developed using the following tools and packages. Ensure that these versions are installed before running the scripts.

<b>Language</b>

&rarr; R (version 4.0 or later)

<b>R Packages:</b>

&rarr; Seurat (version 5.1.0)

&rarr; tidyverse (version 2.0.0)

&rarr; scDblFinder (version 1.16.0)

&rarr; RColorBrewer (version 1.1.3)

<!-- REPOSITORY STRUCTURE -->
## Repository Structure
The pipeline was developed using the following tools and packages. Ensure that these versions are installed before running the scripts.

The repository is organized as follows:
```bash
/jb_2025_SAdv
│
├── /scRNA_cortical               # Scripts for scRNA-Seq from cortical organoids 
│   ├── /diff_exp                 # Scripts for differential expression analysis
│   │   ├── /analysis_clusterProfiler
│   │   │   ├── script_1.sh       # Create symbolic links to unfiltered and diff exp genes exclusive per cell population
│   │   │   ├── script_2.R        # Run GO-BP enrichment analyses
│   │   │   ├── script_3.sh       # Filter only brain-related GO-BP terms
│   │   │   └── script_4.R        # Make bubbleplots without REVIGO filtering
│   │   ├── script_1.sh           # Get number of diff exp genes (log2FC > 0.1 & FDR < 0.05)
│   │   ├── script_2.R            # Filter differentially expressed protein-coding genes (|log2FC| > 0.1 & adj.pvalue < 0.05)
│   │   └── script_3.R            # Plot Venn diagrams with intersections between cell populations
│   ├── script_1.R                # Run QC filtering
│   ├── script_2.R                # Normalize samples
│   ├── script_3.R                # Integrate all samples. (Based on [PMID: 36179669] publication)
│   ├── script_4.R                # Integrate reference datasets (1.5 + 2 months). (Based on [PMID: 36179669] publication)
│   ├── script_5.R                # Map and annotate integrated dataset (1.5 + 2 months) with reference cell types. (Based on [PMID: 36179669] publication)
│   ├── script_6.R                # Plots
│   └── script_7.R                # Differential Expression
│
├── /scRNA_thalamic               # Scripts for scRNA-Seq from thalamic organoids
│   ├── script_1.R                # Run QC filtering
│   ├── script_2.R                # Normalize samples
│   ├── script_3.R                # Integrate all samples. (Based on [PMID: 37824646] publication)
│   ├── script_4.R                # Map and annotate integrated dataset with reference cell types. (Based on [PMID: 37824646] publication)
│   └── script_5.R                # Plots
│
└── /bulk_organoids               # Scripts for bulk RNA-Seq
    ├── script_1.R                # Description of what this script does
    ├── script_2.R                # Description of what this script does
    └── ...
```

<!-- DATA -->
## Data
The raw data utilized in this pipeline is available in ENA (PRJEB83863). To download the files you should use the following command:
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/EXAMPLE/EXAMPLE/EXAMPLE.fastq.gz
```

<!-- Contact -->
## Contact
Filipe F. dos Santos - (fferreira@mochsl.org.br)

Gabriela D. A. Guardia - (gguardia@mochsl.org.br)

Rafael L. V. Mercuri - (rmercuri@mochsl.org.br)

Pedro A. F. Galante - (pgalante@mochsl.org.br)

<!-- Publication -->
## Publication
Joannes-Boyau R., de Souza J.S., Austin C., Westaway K., Moffat I., Wang W., Liao W., Zhang Y., Adams J.W., 
Fiorenza L., Dérognat F., Moncel M., Schwartz G.T., Bailey M., Petroski L. P., Sanchez-Sanchez S. M., Oviedo J., 
Herai R. H., Lemos B., Tonge M., Arora M., Muotri A. R. (2025). Impact of intermittent lead exposure on hominid brain evolution. [JOURNAL], [CAP], [PAGES]. [DOI]
