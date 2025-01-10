#!/usr/bin/env bash

# "Source step #1": runs a recent version of "Kallisto (v0.46.2)" for pseudo-alignment of FASTQ files

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
declare -a CMDS=("version" "cite" "index" "quant")

# 'CLI' ARGs
DATASET=$1  # 'Neandertais'
GROUP=$2  # 'Nina'
GENCODE_VERSION=$3  # 'v36'

# My 'OBJs'
WD=$(pwd)
TOOL="kallisto"  # kallisto_v0.43.1 | 'kallisto_v0.46.1 (v0.46.2)' | kallisto (v0.48.0) --> just type kallisto (Herai did not install the newest version and this one without bootstrapping support)
VERSION=$(${TOOL} ${CMDS[0]} | cut -f2 -d"," | cut -f3 -d" ")  # Version: 'v0.46.2'
CITATION=$(${TOOL} ${CMDS[1]} | grep "doi" | cut -f2 -d":")  # DOI: '10.1038/nbt.3519'
SUFFIX=".fastq.gz"  # 'zipped' files
DONE="\t${BG_GREEN}Done!${NC}\n"

# My 'PATHs'
PATH_FASTQ=${WD}/bulkRNAseq  # for 'bulk RNA-seq' data
PATH_GENCODE=${WD}/../../gencode/release$(echo ${GENCODE_VERSION} | cut -f2 -d"v")  # PATH to 'Human v36 GENCODE (TCGA / GTEx compatibility)'
PATH_INDEX=${PATH_GENCODE}/to${TOOL^}  # PATH to 'Kallisto INDEX'
PATH_KALLISTO=${WD}/$(echo ${TOOL} | cut -f1 -d"_")  # To standardize the 'Kallisto OUTDIR' regardless the version
PATH_LOGS=${PATH_KALLISTO}/logs  # 'LOGDIR'

# Creates 'DIRs' if needed
mkdir -p ${PATH_INDEX}
mkdir -p ${PATH_KALLISTO}
mkdir -p ${PATH_LOGS}

# Change Kallisto efficiency 'PARs' here!
## Careful since the server has '32 cores & 378 RAM'
### 'PS:' no bootstrap support!
THREADS=10
BOOTSTRAPS=100

# OBJs related to the 'progress'
TOTAL=$(ls ${PATH_FASTQ}/*"_R1_"*${SUFFIX} | wc -l)
ITER=0
CONVERSION=1000000000
UNIT="Gb"

# My 'FILES'
## 'PS:' I do not need the GTF, only the FASTA!
FASTA=${PATH_GENCODE}/"gencode."${GENCODE_VERSION}".transcripts.fa"

# Message to 'time file'
echo -e ${DONE}

# Creates 'indexed transcriptome' for Kallisto if necessary
if [[ -s ${PATH_INDEX}/"human_GENCODE"${GENCODE_VERSION}".idx" ]];
then
	echo -e "INDEX for '${TOOL^} (GENCODE${GENCODE_VERSION})' already exists. Skipping '${CMDS[2]}' CMD...\n"
else
	echo -e "Creating INDEX with '${TOOL^} [${CMDS[2]}]'..."

	# FASTA can be 'plaintext | gzipped'
	## 'PS:' Do NOT supply the genome FASTA; it must be a transcriptome one!
	${TOOL} ${CMDS[2]} \
		-i ${PATH_INDEX}/"human_GENCODE"${GENCODE_VERSION}".idx" \
		${FASTA} \
		2> ${PATH_INDEX}/"human_GENCODE"${GENCODE_VERSION}"_index.log"
	echo -e ${DONE}
fi

# 'Iterates' over pairs of 'FASTQ' files
for F in $(ls ${PATH_FASTQ}/*"_R1_"*${SUFFIX}); do
	# 'Updates' ITER variable
	((ITER = ${ITER} + 1))

	# Calculates 'PERC'
	PERC=$(echo "scale=2; (${ITER} / ${TOTAL}) * 100" | bc -l)

	# Checks if file 'exists & not empty'
	echo -e "Checking if 'FASTQ (R1)' file ${BOLD_BLUE}exists & is not empty${NC}..."
	if [[ -s ${F} ]];
	then
		# Gets 'size'
		## Multiplying by 2 to account for the paired-end (PE) nature of FASTQ reads
		SIZE=$(stat -Lc %s ${F})
		SIZE_GB=$(echo "scale=1; (${SIZE} * 2) / ${CONVERSION}" | bc -l)

		# Message to 'time file'
		echo -e "${DONE}\nStarting ${BOLD_PURPLE}'${TOOL}'${NC} [${PURPLE}${CMDS[3]}${NC}] tool (version: ${ITALIC_GREEN}${VERSION}${NC}; DOI: ${ITALIC_GREEN}${CITATION}${NC}) with ${BOLD_ORANGE}${THREADS} threads & ${BOOTSTRAPS} bootstraps${NC}..."

		# Gets 'filename'
		NAME=$(basename ${F} "_R1_001"${SUFFIX})

		# Infers 'READ pair' names
		R1=${PATH_FASTQ}/${NAME}"_R1_001"${SUFFIX}
		R2=${PATH_FASTQ}/${NAME}"_R2_001"${SUFFIX}

		# Creates 'OUTDIR' name
		DIR=${PATH_KALLISTO}/${NAME}
		mkdir -p ${DIR}

		# Message to 'time file'
		echo -e "\tProccessing sample ${GREEN}${NAME}${NC} (${BOLD_RED}${SIZE_GB} ${UNIT}${NC}) [${BOLD_GREEN}${ITER}/${TOTAL}${NC} ${ITALIC_GREEN}(${PERC}%)${NC}]..."

		# Runs 'Kallisto'
		${TOOL} ${CMDS[3]} \
			-i ${PATH_INDEX}/"human_GENCODE"${GENCODE_VERSION}".idx" \
			-o ${DIR} \
			-b ${BOOTSTRAPS} -t ${THREADS} \
			${R1} ${R2} \
			2> ${PATH_LOGS}/${NAME}".log"
		echo -e ${DONE}
	else
		echo -e "\t${RED}ERROR:${NC}File does not exist | is empty! Unable to proceed to ${BOLD_ORANGE}'${TOOL^}'${NC} pseudo-alignment. Verify DIR content again. Checking next file...\n"
	fi
done

################################################## 'PS:' If running this script for the first time, create the MAPS. Otherwise, it is not necessary! ##################################################
# Creates 'MAPS' between 'transcripts and genes IDs' - same logic for mouse and humans!
grep ">" ${FASTA} | sed "s/>//" | awk -F "|" '{print$0"\t"$2}' > ${WD}/"ID2Gene.txt"  # Entire FASTA header to gene '(ID --> ENSG)'
grep ">" ${FASTA} | sed "s/>//" | awk -F "|" '{print$1"\t"$2}' > ${WD}/"ENST-2-ENSG.txt"  # Just ID conversions '(ENST --> ENSG)'

# Creates 'MAPS' between 'gene / isoform IDs and symbols'
grep ">" ${FASTA} | sed "s/>//" | awk -F "|" '{print$2"\t"$6"\t"$8}' | sort | uniq > ${WD}/"map_gene_IDs.txt"  # 'Gene MAP' infos (ENSG, Symbol and Category)
grep ">" ${FASTA} | sed "s/>//" | awk -F "|" '{print$1"\t"$5"\t"$8}' | sort > ${WD}/"map_isoform_IDs.txt"  # 'Transcript MAP' infos (ENST, Symbol and Category)
grep ">" ${FASTA} | sed "s/>//" | awk -F "|" '{print$1"\t"$2"\t"$5"\t"$6"\t"$7"\t"$8}' | sort > ${WD}/"map_complete.txt"  # 'Complete MAP (all)' infos (ENST, ENSG, Symbols, Length and Categories)
#######################################################################################################################################################################################################

# Prints final message to 'time file'
echo -e "${UND_BLUE}Finished proccessing${NC} of all ${BOLD_ORANGE}${TOTAL}${NC} ${BOLD_RED}${DATASET}-${GROUP} FASTQ (PE)${NC} files!\n"
