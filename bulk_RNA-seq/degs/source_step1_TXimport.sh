#!/usr/bin/env bash

# "Source step #1:" runs "TXImport" in order to prepare data for the next step using "R4 - Rscript (4.3.1 / 2023-06-16 / Beagle Scouts)"

# 'CLI' ARGs
COMPARISON=$1  # 'Comparison' to be done
BATCH=$2  # 'Batch' to be used

# My 'OBJs'
WD=$(pwd)
DONE="\tDone!\n"

# My 'PATHs'
# "PS:" Change the INPUT file DIR here!
PATH_INPUT=${WD}/../"kallisto"/${BATCH}/${COMPARISON}
PATH_OUTPUT=${WD}/"tximport"
PATH_LOG=${WD}/"logs"

# Creates 'DIRs' if needed
mkdir -p ${PATH_OUTPUT}
mkdir -p ${PATH_LOG}

# Uses 'tximport (v1.24.0)' R library
## ARG1 = 'Mapping IDs file (ENST --> ENSG)' as in raw Kallisto abundance.tsv (as a consequence, h5) file;
## ARG2 = 'INDIR' where Kallisto H5 (binary) OUT files are located;
## ARG3 = 'OUTDIR' where you want the outputs to be placed;
## ARG4 = 'COMPARISON' of condition vs control.
echo -e "Running R4 script: 'st1_TXimport.R' for comparison '${COMPARISON}' (${BATCH})..."
Rscript --vanilla ${WD}/st1_TXimport.R \
	${WD}/../"ID2Gene.txt" \
	${PATH_INPUT} \
	${PATH_OUTPUT} \
	${COMPARISON} \
	${BATCH} \
	2> ${PATH_LOG}/"st1_TXimport_${BATCH}-${COMPARISON}.log"
echo -e ${DONE}
