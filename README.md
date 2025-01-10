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

<details open="open">
    <summary><b>Language</b></summary>
	<li>&rarr; Python3</li>

	<li>&rarr; R (version 4.0 or later)</li>
<details open="open">

    <b>Python Packages:</b>

	<li>&rarr; Kallisto-Bustools</li>

    <b>R Packages:</b>

	<li>&rarr; DeSeq2</li>

	<li>&rarr; tximport</li>

	<li>&rarr; readr</li>

	<li>&rarr; icesTAF</li>

	<li>&rarr; Seurat (version 5.1.0)</li>

	<li>&rarr; tidyverse (version 2.0.0)</li>

	<li>&rarr; scDblFinder (version 1.16.0)</li>

	<li>&rarr; RColorBrewer (version 1.1.3)</li>

	<li>&rarr; ggplot2</li>

	<li>&rarr; ggtext</li>

	<li>&rarr; ggnewscale</li>

	<li>&rarr; forcats</li>

	<li>&rarr; patchwork</li>

	<li>&rarr; Cairo</li>

	<li>&rarr; extrafont</li>

	<li>&rarr; ggpubr</li>

	<li>&rarr; ggrepel</li>
   <ol>
</details>

<!-- REPOSITORY STRUCTURE -->
## Repository Structure
The pipeline was developed using the following tools and packages. Ensure that these versions are installed before running the scripts.

The repository is organized as follows:
```bash
/jb_2025_SAdv
├── bulk_RNA-seq # Scripts for bulk RNA-Seq
│   ├── degs # Differential expression analysis
│   │   ├── automate_MASTER_call.sh
│   │   ├── source_step0_MASTER.sh
│   │   ├── source_step1_TXimport.sh
│   │   ├── source_step2_DESeq2.sh
│   │   ├── st1_TXimport.R
│   │   ├── st2_DESeq2.R
│   │   ├── st3_Filters.R
│   │   └── support
│   │       ├── edit_PC-Gene_MAP.R
│   │       ├── get_DEG_table.R
│   │       ├── get_Union-DEG.R
│   │       └── merge_TXimport_PC-GeneMAP.R
│   ├── metadata
│   │   ├── map_gene_IDs.PC.tsv
│   │   ├── map_gene_IDs.raw.tsv
│   │   ├── map_gene_IDs.tsv
│   │   └── mapID2GROUP.batch2.tsv
│   ├── metrics
│   │   └── MultiQC.html
│   ├── plots
│   │   ├── automate_PCA_call.sh
│   │   ├── CMDs_chatGPT.txt
│   │   ├── CMDs_cytoscape.txt
│   │   ├── CMDs_enrichment.txt
│   │   ├── make_bubblePlot.R
│   │   ├── make_PCA_COUNTS.R
│   │   └── make_VolcanoPlot.R
│   └── pre-processing
│       ├── source_step0_MASTER.sh
│       ├── st0_Metrics.sh
│       └── st1_Kallisto.sh
├── LICENSE
├── README.md
└── scRNA-seq # Scripts for scRNA-Seq
    ├── cortical_organoids # Scripts for scRNA-Seq from cortical organoids analysis
    │   ├── 5_diff_exp_mixed_1.5_2M # Scripts for differential expression analysis
    │   │   ├── source_step1.R
    │   │   ├── source_step2.R
    │   │   ├── source_step3.R
    │   │   ├── source_step4
    │   │   └── source_step5.R
    │   ├── list_dirs.txt
    │   ├── map_genes.txt
    │   ├── sample_groups.txt
    │   ├── source_step1.R
    │   ├── source_step2.R
    │   ├── source_step3.R
    │   ├── source_step4.R
    │   ├── source_step5.R
    │   ├── source_step6.R
    │   └── source_step7.R
    ├── Organoid_Markers.tsv
    ├── source_step0_MASTER.sh
    ├── st1_KB_Ref.sh # Build Index - kb
    ├── st2_KB_Count.geneNoMM.sh # Process Data - kb
    └── thalamic_organoids # Scripts for scRNA-Seq from thalamic organoids analysis
        ├── list_dirs.txt
        ├── source_step1.R
        ├── source_step2.R
        ├── source_step3.R
        ├── source_step4.R
        ├── source_step5.R
        └── source_step6.R
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
