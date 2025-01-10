#!/usr/bin/env bash

# "Source step #2": runs the "KB (v0.28.2)" wrapper "Kallisto (v0.50.1) + Bustools (v0.43.2)" for generating the count matrices
## scRNA-seq: "10XV3" --> check Ivan emails about the AmpliDrop tech used!
### "CMDs" --> $ kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa dna.primary_assembly.fa.gz gtf.gz
### "CMDs" --> $ kb count -i index.idx -g t2g.txt -x 10xv3 --h5ad -t 2 read_1.fastq.gz read_2.fastq.gz
#### "PAPER" for REF --> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10690192/
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
# 'PS:' adding ref CMD in list just to automatically avoid INDEX file not found error!
declare -a CMDS=("info" "count" "ref")

# Array of 'SUFFIXES' and 'EXTENSIONS'
declare -a SUFFIXES=(".fastq.gz")
declare -a EXTS=(".idx" ".txt" ".log")

# 'CLI' user ARGs
TYPE=$1  # custom ARGs [--mm | --tcc] ('geneNoMM' | geneMM | tccNoMM | tccMM)
GENCODE_VERSION=$2  # 'v36' (TCGA)
GENOME=$3  # 'hg38'

# My 'OBJs'
WD=$(pwd)
TOOL="kb"  # 'wrapper' of Kallisto and Bustools
VERSION="v"$(${TOOL} ${CMDS[0]} | head -n1 | cut -f2 -d" ")  # Version: 'v0.28.2'
CITATION=$(${TOOL} ${CMDS[0]} | grep -A10 "==========" | grep -A3 "\[1\]\|\[2\]" | tr '\n' ' ')  # Kallisto and BUStools 'papers'
DONE="\t${BG_GREEN}Done!${NC}\n"

# My 'PATHs'
PATH_KB=${WD}/$(echo ${TOOL} | cut -f1 -d"_")  # To standardize the 'KB OUTDIR' regardless the version
PATH_OUT=${PATH_KB}/${TYPE}  # 'OUTDIR'
PATH_LOGS=${PATH_KB}/logs  # 'LOGDIR'
PATH_TMP=${PATH_KB}/kbTMP  # 'TMPDIR'
PATH_FASTQ=${WD}/fastq  # PATH to 'FASTQ' depending on the dataset and group

# Creates 'DIRs'
## 'WARN:' we cannot create beforehand the custom TMP DIR, otherwise we get an ERROR!
mkdir -p ${PATH_KB}
mkdir -p ${PATH_OUT}
mkdir -p ${PATH_LOGS}
#mkdir -p ${PATH_TMP}

# Changes KB efficiency 'PARs' here!
THREADS=15  # more if in 'PND9/10'
MEMORY="10G"  # more if in 'PND9/10'

# OBJs related to the 'progress'
TOTAL=$(ls ${PATH_FASTQ}/*"_S1_L001_R1_001"${SUFFIXES[0]} | wc -l)
ITER=0
CONVERSION=1000000000
UNIT="Gb"

# Message to 'time file'
echo -e ${DONE}
echo -e "Starting ${BOLD_PURPLE}'${TOOL}'${NC} [${PURPLE}${CMDS[1]} (no '--tcc' / no '--mm' ARGs)${NC}] tool (version: ${ITALIC_GREEN}${VERSION}${NC}; DOI: ${ITALIC_GREEN}${CITATION}${NC}) with ${ORANGE}(Human GENCODE${GENCODE_VERSION} / ${GENOME^^})${NC} using ${BOLD_ORANGE}${THREADS} threads & ${MEMORY} memory${NC}..."

# 'Iterates' over pairs of 'FASTQ' files
for F in $(ls ${PATH_FASTQ}/*"_S1_L001_R1_001"${SUFFIXES[0]}); do
	# 'Updates' ITER variable
	((ITER = ${ITER} + 1))

	# Calculates 'PERC'
	PERC=$(echo "scale=2; (${ITER} / ${TOTAL}) * 100" | bc -l)

	# Checks if file 'exists'
	echo -e "\tChecking if 'FASTQ (R1)' file ${BOLD_BLUE}exists${NC}..."
	if [[ -f ${F} ]];
	then
		echo -e "\t${DONE}\n\tChecking if it is ${BOLD_BLUE}empty${NC}..."

		# Checks if file is 'empty'
		if [[ -s ${F} ]];
		then
			# Gets 'size'
			## Multiplying by 2 to account for the paired-end (PE) nature of FASTQ reads
			SIZE=$(stat -Lc %s ${F})
			SIZE_GB=$(echo "scale=1; (${SIZE} * 2) / ${CONVERSION}" | bc -l)

			# Message to 'time file'
			echo -e "\t${DONE}"

			# Gets 'filename'
			NAME=$(basename ${F} "_S1_L001_R1_001"${SUFFIXES[0]})

			# Infers 'READ pair' names
			R1=${PATH_FASTQ}/${NAME}"_S1_L001_R1_001"${SUFFIXES[0]}
			R2=${PATH_FASTQ}/${NAME}"_S1_L001_R2_001"${SUFFIXES[0]}

			# Creates 'OUTDIR' name
			DIR=${PATH_OUT}/${NAME}
			mkdir -p ${DIR}

			# Message to 'time file'
			echo -e "\tProccessing sample ${GREEN}${NAME}${NC} (${BOLD_RED}${SIZE_GB} ${UNIT}${NC}) [${BOLD_GREEN}${ITER}/${TOTAL}${NC} ${ITALIC_GREEN}(${PERC}%)${NC}]..."

			# Runs the 'kb count' CMD
			## 'ARGs [REQ]:' -i INDEX | -g T2G (MAP) | -x TECHNOLOGY
			## 'ARGs [OPT]:' --tmp TMP | --verbose (?) | -o OUT | -t THREADS | -m MEMORY | --h5ad / --loom --> another way to efficiently store scRNAseq data
			## 'ARGs [POS]:' fastqs
			### 'PS:' using DEFAULT for --> --workflow standard | -w WHITELIST | --filter bustools
			### 'PS:' we must choose either the H5AD or the LOOM format, not both!
			### 'PS:' QUESTIONS --> --mm (includes multi-mapping reads) | --tcc (TCC matrix, not gene matrix) | --dry-run (for testing)
			#### 'WARN:' run with --tcc (TCC matrix) first and then without it (GENE matrix). Same thing for the --mm (multi-mapping reads)!
			${TOOL} ${CMDS[1]} \
				-i ${PATH_KB}/"human."${TOOL}"-"${CMDS[2]}"."${GENCODE_VERSION}"."${GENOME}${EXTS[0]} \
				-g ${PATH_KB}/"mapT2G_human."${TOOL}"-"${CMDS[2]}"."${GENCODE_VERSION}"."${GENOME}${EXTS[1]} \
				-x "10XV3" \
				--tmp ${PATH_TMP} --verbose \
				-o ${DIR} \
				-m ${MEMORY} -t ${THREADS} \
				--h5ad \
				${R1} ${R2} \
				2> ${PATH_LOGS}/"count_human_"${NAME}"."${TOOL}"-"${CMDS[1]}"."${GENCODE_VERSION}"."${GENOME}"."${TYPE}${EXTS[2]}
			echo -e "\t${DONE}"
		else
			echo -e "\t${RED}ERROR:${NC}File is empty! Unable to proceed to ${BOLD_ORANGE}'${TOOL}'${NC} pseudo-alignment / count. Checking next file...\n"
		fi
	else
		echo -e "\t${RED}ERROR:${NC}File does not exist! Verify DIR content again. Checking next file...\n"
	fi
done

# Message to 'time file'
echo -e ${DONE}

# Prints final message to 'time file'
echo -e "${UND_BLUE}Finished proccessing${NC} of all ${BOLD_ORANGE}${TOTAL}${NC} ${BOLD_RED}FASTQ (PE)${NC} files!\n"
