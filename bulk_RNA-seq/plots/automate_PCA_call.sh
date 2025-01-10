#!/usr/bin/env bash

# "Automates" calling of the "PCA R script" by looping through the "comparisons"

# My 'OBJs'
WD=$(pwd)
SCRIPT="make_PCA_COUNTS.R"
DONE="\tDone!\n"

# Declares 'array' of comparisons
declare -a COMPARISONS=("ASC-10_vs_ASC-30" "ASC-10_vs_ASC-CTRL" "ASC-30_vs_ASC-CTRL" \
                       "C136-10_vs_C136-30" "C136-10_vs_C136-CTRL" "C136-30_vs_C136-CTRL" \
                       "NOVA1-C15-10_vs_NOVA1-C15-30" "NOVA1-C15-10_vs_NOVA1-C15-CTRL" "NOVA1-C15-30_vs_NOVA1-C15-CTRL" \
                       "NOVA1-C27-10_vs_NOVA1-C27-30" "NOVA1-C27-10_vs_NOVA1-C27-CTRL" "NOVA1-C27-30_vs_NOVA1-C27-CTRL" \
		       "ASC-CTRL_vs_C136-CTRL" "NOVA1-C27-CTRL_vs_NOVA1-C15-CTRL")

# 'Iterates' over comparisons
echo -e "Running '${SCRIPT}' automatically for all ${#COMPARISONS[@]} comparisons from the Neanderthal Project (Nina colaboration)...\n"
for MYCOMP in ${COMPARISONS[@]}; do
	# Gets 'groups' and 'size' for the control one
	GROUP1=$(echo ${MYCOMP} | cut -f1 -d"_")  # 'Condition'
	GROUP2=$(echo ${MYCOMP} | cut -f3 -d"_")  # 'Control'
	NCTRL=2  # 'Sample size' for the 'control' group

	# Runs 'PCA' figure '(only protein-coding genes)' and saves 'LOG'
	## COUNTS = raw counts data from 'TXimport'
	## PASS = '|L2FC| > 2'
	echo -e "\tComparison: ${MYCOMP}\n"
	Rscript --vanilla ${WD}/${SCRIPT} ${WD}/TXimport ${WD}/Plots/PCAs ${GROUP1} ${GROUP2} ${NCTRL} "COUNTS" "PASS" "protein_coding" 2> ${WD}/"collab_Nina_Neanderthal_PCA_"${MYCOMP}".log"
	echo -e "\t${DONE}"
done
echo -e ${DONE}
