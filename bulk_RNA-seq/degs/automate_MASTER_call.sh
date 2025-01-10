#!/usr/bin/env bash

# "Automates" calling of the "MASTER shell script" by looping through the "batches & comparisons"

# My 'OBJs'
WD=$(pwd)
SCRIPT="source_step0_MASTER.sh"
DONE="\tDone!\n"

# Declares 'arrays' of batches & comparisons
#declare -a BATCHES=("batch1-2" "batch1" "batch2")
declare -a BATCHES=("batch2")
#declare -a COMPARISONS=("ASC-10_vs_ASC-30" "ASC-10_vs_ASC-CTRL" "ASC-30_vs_ASC-CTRL" \
#			"C136-10_vs_C136-30" "C136-10_vs_C136-CTRL" "C136-30_vs_C136-CTRL" \
#			"NOVA1-C15-10_vs_NOVA1-C15-30" "NOVA1-C15-10_vs_NOVA1-C15-CTRL" "NOVA1-C15-30_vs_NOVA1-C15-CTRL" \
#			"NOVA1-C27-10_vs_NOVA1-C27-30" "NOVA1-C27-10_vs_NOVA1-C27-CTRL" "NOVA1-C27-30_vs_NOVA1-C27-CTRL")
declare -a COMPARISONS=("NOVA1-ArAr-10uM_vs_NOVA1-ArAr-CTRL" "NOVA1-ArAr-30uM_vs_NOVA1-ArAr-CTRL" "NOVA1-HuHu-10uM_vs_NOVA1-HuHu-CTRL" "NOVA1-HuHu-30uM_vs_NOVA1-HuHu-CTRL")

# 'Iterates' over comparisons
echo -e "Running '${SCRIPT}' automatically for all ${#COMPARISONS[@]} comparisons from all ${#BATCHES[@]} batches of the Neanderthals project (Nina colaboration)...\n"
for BATCH in ${BATCHES[@]}; do
	echo -e "\tBatch: ${BATCH}\n"
	for MYCOMP in ${COMPARISONS[@]}; do
		echo -e "\t\tComparison: ${MYCOMP}\n"
		#bash ${WD}/${SCRIPT} ${MYCOMP} ${BATCH} |& tee -a ${WD}/"collab_Nina-Neanderthals_"${BATCH}"-"${MYCOMP}".log"  # saving the 'master LOGs' to separate files
		bash ${WD}/${SCRIPT} ${MYCOMP} ${BATCH}
		echo -e "\t\t${DONE}"
	done
	echo -e "\t${DONE}"
done
echo -e ${DONE}
