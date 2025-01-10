#!/usr/bin/env bash

# "Source step #1": runs the "KB (v0.28.2)" wrapper "Kallisto (v0.50.1) + BUStools (v0.43.2)" for building the reference
## scRNA-seq: "10XV3" --> check Ivan emails about the AmpliDrop tech used!
### "CMDs" --> $ kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa dna.primary_assembly.fa.gz gtf.gz
### "CMDs" --> $ kb count -i index.idx -g t2g.txt -x 10xv3 --h5ad -t 2 read_1.fastq.gz read_2.fastq.gz
##### "PS:" do not forget to run the kb wrapper with conda! (conda activate kb_env)

# Sets 'colors' --> https://www.codegrepper.com/code-examples/shell/bash+background+color
## If I change the '0' by '1' right before the ';', the text turns 'bold'!
## If I change the '0' by '2' right before the ';', the text turns 'kind of grey'!
## If I change the '0' by '3' right before the ';', the text turns 'italic'!
## If I change the '0' by '4' right before the ';', the text turns 'underlined'!
## If I change the '0' by '5 | 6' right before the ';', the text 'blinks'!
## If I change the '0' by '9' right before the ';', the text has a 'strokeline'!
GREY='\033[1;90m'  # 'Bold' Grey
GREEN='\033[0;32m'  # Green
BOLD_GREEN='\033[1;32m'  # 'Bold' Green
ITALIC_GREEN='\033[3;32m'  # 'Italic' Green
BG_GREEN='\033[1;42m'  # 'Bold' text and green 'Background'
BLUE='\033[0;34m'  # Blue
BOLD_BLUE='\033[1;34m'  # 'Bold' Blue
ITALIC_BLUE='\033[3;34m'  # 'Bold' Blue
UND_BLUE='\033[4;34m'  # 'Underlined' Blue
PURPLE='\033[0;35m'  # Purple
BOLD_PURPLE='\033[1;35m'  # 'Bold' Purple
ORANGE='\033[0;33m'  # Brown/Orange
BOLD_ORANGE='\033[1;33m'  # 'Bold' Orange
RED='\033[0;31m'  # Red
BOLD_RED='\033[1;31m'  # 'Bold' Red
NC='\033[0m'  # No Color

# Message to 'time file'
echo -e "Setting up ${GREEN}ENV${NC}..."

# Array of 'CMDs'
declare -a CMDS=("info" "ref")

# Array of 'SUFFIXES' and 'EXTENSIONS'
declare -a SUFFIXES=(".fa" ".gtf")
declare -a EXTS=(".idx" ".txt" ".log")

# 'CLI' user ARGs
GENCODE_VERSION=$1  # 'v36' (TCGA)
GENOME=$2  # 'hg38'

# My 'OBJs'
WD=$(pwd)
TOOL="kb"  # 'wrapper' of Kallisto and BUStools
VERSION="v"$(${TOOL} ${CMDS[0]} | head -n1 | cut -f2 -d" ")  # Version: 'v0.28.2'
CITATION=$(${TOOL} ${CMDS[0]} | grep -A10 "==========" | grep -A3 "\[1\]\|\[2\]" | tr '\n' ' ')  # Kallisto and Bustools 'papers'
DONE="\t${BG_GREEN}Done!${NC}\n"

# My 'PATHs'
PATH_KB=${WD}/$(echo ${TOOL} | cut -f1 -d"_")  # To standardize the 'KB OUTDIR' regardless the version
PATH_LOGS=${PATH_KB}/logs  # 'LOGDIR'
PATH_TMP=${PATH_KB}/kbTMP  # 'TMPDIR'

# Creates 'DIRs'
## 'WARN:' we cannot create beforehand the custom TMP DIR, otherwise we get an ERROR!
mkdir -p ${PATH_KB}
mkdir -p ${PATH_LOGS}
#mkdir -p ${PATH_TMP}

# My 'FILES'
## 'WARN:' kb ref 'requires' the 'reference genome' as the 'FASTA', not the one used in Kallisto directly!
FASTA=${WD}/${GENOME}".25chr"${SUFFIXES[0]}  # 'FASTA' -> 'human ref genome hg38 25chr-filtered'
#FASTA=${WD}/${GENOME}${SUFFIXES[0]}  # 'FASTA' -> 'human ref genome hg38 with no filters --> this is the one which has indeed the chrM'
GTF=${WD}/"gencode."${GENCODE_VERSION}".annotation"${SUFFIXES[1]}  # 'GTF' -> 'original'

# Message to 'time file'
echo -e ${DONE}

# Creates 'INDEX'
if [[ -s ${PATH_KB}/"human."${TOOL}"-"${CMDS[1]}"."${GENCODE_VERSION}"."${GENOME}${EXTS[0]} && -s ${PATH_KB}/"mapT2G_human."${TOOL}"-"${CMDS[1]}"."${GENCODE_VERSION}"."${GENOME}${EXTS[1]} ]];
then
        echo -e "${BOLD_ORANGE}INDEX & T2G MAP${NC} files for ${BOLD_PURPLE}'${TOOL^^}'${NC} ${ORANGE}(Human GENCODE${GENCODE_VERSION} / ${GENOME^^})${NC} already ${BOLD_BLUE}created${NC}. ${PURPLE}Skipping '${CMDS[1]}'${NC} CMD...\n"
else
	echo -e "Creating ${BOLD_ORANGE}INDEX & T2G MAP${NC} files with ${BOLD_PURPLE}'${TOOL^^}${NC} ${PURPLE}[${CMDS[1]}]'${NC} for ${ORANGE}'Human GENCODE${GENCODE_VERSION} / ${GENOME^^}'${NC}..."

	# Runs the 'kb ref' CMD
	## 'ARGs (used):' --tmp TMP [OPT] | --verbose [OPT] & -i [REQ] | -g [REQ] | -f1 [REQ] & fasta [POS] | gtf [POS]
	### 'PS:' using DEFAULT --workflow standard and not using feature POS ARG and -d (I wanna build my own INDEX)
	#### I could, however, use a 'prebuilt index (Ensembl 108 GRCh38)' from: https://github.com/pachterlab/kallisto-transcriptome-indices
	### 'PS2:' using DEFAULT -n (no splitting in FASTA) and -k (recommended k-mer) ARGs
	##### 'WARN:' kb ref uses a 'ref genome FASTA' and a 'matching GTF (by chromosome name)' to 'create its own transcriptome FASTA'. This is the one 'which gets indexed!'
	${TOOL} ${CMDS[1]} \
		-i ${PATH_KB}/"human."${TOOL}"-"${CMDS[1]}"."${GENCODE_VERSION}"."${GENOME}${EXTS[0]} \
		-g ${PATH_KB}/"mapT2G_human."${TOOL}"-"${CMDS[1]}"."${GENCODE_VERSION}"."${GENOME}${EXTS[1]} \
		-f1 ${PATH_KB}/"transcriptome_human."${TOOL}"-"${CMDS[1]}"."${GENCODE_VERSION}"."${GENOME}${SUFFIXES[0]} \
		--tmp ${PATH_TMP} \
		--verbose \
		${FASTA} \
		${GTF} \
		2> ${PATH_LOGS}/"index_human."${TOOL}"-"${CMDS[1]}"."${GENCODE_VERSION}"."${GENOME}${EXTS[2]}
	echo -e ${DONE}
fi
