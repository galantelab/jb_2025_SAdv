#!/usr/bin/env bash

# "Source step #2:" runs "DESeq2" differential expression gene (DEG) analyses with "R4 - Rscript (4.3.1 / 2023-06-16 / Beagle Scouts)"

# 'CLI' ARGs
COMPARISON=$1  # 'Comparison' to be done
BATCH=$2  # 'Batch' to be used

# My 'OBJs'
## It is recommended (but not required) to use only 'letters', 'numbers', and 'delimiters' ('_' or '.') for the GROUPS as these are 'safe characters' for 'R COL names'
WD=$(pwd)
DONE="\tDone!\n"
GROUP1=$(echo ${COMPARISON} | cut -f1 -d"_")  # 'Condition'
GROUP2=$(echo ${COMPARISON} | cut -f3 -d"_")  # 'Control'
NCTRL=$(ls ${WD}/../kallisto/${BATCH}/${COMPARISON}/${GROUP2} | wc -l)  # 'Sample size' for the 'control' group

# Array of 'scripts'
declare -a SCRIPTS=("st2_DESeq2.R" "st3_Filters.R")

# My 'PATHs'
PATH_INPUT=${WD}/"tximport"
PATH_OUTPUT=${WD}/"deseq2"
PATH_LOG=${WD}/"logs"

# Creates 'DIRs' if needed
mkdir -p ${PATH_OUTPUT}
mkdir -p ${PATH_LOG}

# Lists 'CONTROL file IDs' --> 'IMPORTANT CHECKPOINT (MATCH WITH THE ONES FROM R ITSELF)!'
echo -e "File IDs from 'CONTROL (${GROUP2})' [CLI]:\n$(ls ${WD}/../kallisto/${BATCH}/${COMPARISON}/${GROUP2})\n"

# 'TESTS'
echo -e "\tBatch: ${BATCH}\n"
echo -e "\tGroup 1: ${GROUP1}\n"
echo -e "\tGroup 2: ${GROUP2}\n"
echo -e "\tSample Size CONTROL: ${NCTRL}\n"

## Uses 'DESeq2 (v1.36.0)' R library
### ARG1 = 'Condition' the group to be tested;
### ARG2 = 'Control' to be the control group;
### ARG3 = 'INDIR' after TXimport;
### ARG4 = 'OUTDIR' of DESeq2;
### ARG5 = 'nCtrl' (sample size) for Control group;
### ARG6 = 'myBATCH' batch from which comparison will be made.
echo -e "Running R script '${SCRIPTS[0]}' comparing '${GROUP1} (condition)' vs. '${GROUP2} (control)' from '${BATCH}'..."
Rscript --vanilla ${WD}/${SCRIPTS[0]} \
	${GROUP1} \
	${GROUP2} \
	${PATH_INPUT} \
	${PATH_OUTPUT} \
	${NCTRL} \
	${BATCH} \
	2> ${PATH_LOG}/$(echo ${SCRIPTS[0]} | cut -f1 -d".")"_"${BATCH}"-"${COMPARISON}".log"
echo -e ${DONE}

## Filters: '|L2FC| > 1' & 'FDR <= .05'
### ARG1 = 'MAP' file (ENSG --> Symbol & Type);
### ARG2 = 'IODIR';
### ARG3 = 'COMPARISON';
### ARG4 = 'BATCH'.
echo -e "Running R script '${SCRIPTS[1]}' comparing '${GROUP1} (condition)' vs. '${GROUP2} (control)' from '${BATCH}'..."
Rscript --vanilla ${WD}/${SCRIPTS[1]} \
	${WD}/"map_pcgene_IDs.txt" \
	${PATH_OUTPUT} \
	${COMPARISON} \
	${BATCH} \
	2> ${PATH_LOG}/$(echo ${SCRIPTS[1]} | cut -f1 -d".")"_"${BATCH}"-"${COMPARISON}".log"
echo -e ${DONE}
