#!/usr/bin/env bash

###########################################################################################################################################################
# Author: Filipe Ferreira dos Santos
# Date: 08/18/2024

# This shell script aims to execute all "stXXX.sh" scripts for executing the scRNA-seq pipeline

# "PS:" in order to run the KB wrapper in the server (any machine), Daniel told me to first run:
## (i) conda activate kb_env --> in TMUX
## (ii) conda execute kb --list --> or any other kb CMD in script (actually we can just run directly the CMD!)
## (iii) conda deactivate --> in script / TMUX

# https://pt.stackoverflow.com/questions/143822/qual-a-diferen%C3%A7a-de-bin-bash-e-usr-bin-env-bash --> explica a diferenca entre linhas acima shebang etc
###########################################################################################################################################################

####################################
##### "0) - USAGE AND OPTIONS" #####
####################################

# Gets current 'Date & Time'
TIMELOG=$(date +'%d/%m/%Y %T')
DATE=$(echo ${TIMELOG} | cut -f1 -d" " | sed 's/\//-/g')

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

# 'CLI' ARGs
TRIM="no"  # TRIM reads? (yes | 'no')
GENCODE="v36"  # GENCODE transcriptome version ('v36 - TCGA' | v37 - PhD | v44 - LATEST)
GENOME="hg38"  # human reference genome version (hg19 | 'hg38' | t2t)
TYPE="geneNoMM"  # type of kb count analysis ('geneNoMM' | geneMM | tccNoMM | tccMM --> not allowed, which makes sense, actually)

# My 'OBJs'
WD=$(pwd)
LOGS=${WD}/logs
TMPDIR=${WD}/myTMPdir
MACHINE=${HOSTNAME}
MY_PID=$$
DONE="\t${BG_GREEN}Done!${NC}\n"

# Sets MYUSER as 'root' until Dani changes it to 'fferreira' (USER) inside Docker
MYUSER=${USER}

# Script 'VERSION'
VERSION="v1.0"

# Displays 'HELP'
display_usage()
{
        USAGE="${BOLD_BLUE}Usage:${NC} ${BOLD_GREEN}$0${NC} ${BOLD_ORANGE}[OPTIONS]${NC} ${GREY}-->${NC} Pipeline for scRNA-seq analyses (pre-processing)"
        echo -e "${GREY}===========================================================================================================================================================${NC}\n${USAGE}\n"
        echo -e "This shell script (${ITALIC_GREEN}${VERSION}${NC}) runs automatically the KB (Kallisto & BUStools) wrapper tool and others --> user ${ITALIC_GREEN}'${MYUSER}'${NC}.\n"
        echo -e "${BOLD_ORANGE}Options:${NC}\n\t${ORANGE}-v, --version${NC}\tdisplays this script version\n\t${ORANGE}-h, --help${NC}\tdisplays this help message"
        echo -e "\t${ORANGE}-l, --license${NC}\tdisplays the license for this script\n"
        echo -e "${BOLD_PURPLE}Contact${NC}\n${GREEN}Email:${NC} fferreira@mochsl.org.br\n${GREEN}Phone:${NC} +55 (11) 94541-9392\n"
        echo -e "${BOLD_RED}License${NC}\nGNU General Public License (GPLv3). For more detailed information, please check out this link: ${UND_BLUE}https://www.gnu.org/licenses/gpl-3.0.html${NC}"
        echo -e "or run this script with the '-l, --license' ARG option."
        echo -e "${GREY}===========================================================================================================================================================${NC}"
}

# Checks 'ARGs' provided by user
parse_args()
{
	# Checks if user used the 'help option'
	if [[ ${#} -eq 1 && ( ${1} == "--help" || ${1} == "-h" ) ]];
	then
		display_usage
		exit 0
	else
		# Checks if user passed the 'version' option
		if [[ ${#} -eq 1 && ( ${1} == "--version" || ${1} == "-v" ) ]];
		then
			echo -e "Script ${BOLD_GREEN}'${0}'${NC} version: ${ITALIC_GREEN}${VERSION}${NC}"
			exit 0
		else
			# Checks if user passed the 'license' option
			if [[ ${#} -eq 1 && ( ${1} == "--license" || ${1} == "-l" ) ]];
			then
				cat ${WD}/license/"GNU_GPLv3.txt"
				exit 0
			else
				:
			fi
		fi
	fi
}

# Calls the 'parse_args' FUN
case ${#} in
	1)
		parse_args ${1};;
	0)
		parse_args ${@}
esac

######################################
##### "1) - SETUP MASTER SCRIPT" #####
######################################

# Starts 'general time'
STARTTIME_GENERAL=$(date +%s)

# My 'FUNs'
# How to https://stackoverflow.com/questions/10822790/can-i-call-a-function-of-a-shell-script-from-another-shell-script
check_machines()
{
	if [[ $1 -eq 1 ]];
	then
		echo -e "\n============================== ${BOLD_BLUE}${TIMELOG}${NC} ==============================\n" |& tee -a ${MACHINE_FILE}
		top -b -p $2 -n1 | tail -n2 |& tee -a ${MACHINE_FILE}
	else
		top -b -p $2 -n1 | tail -n1 |& tee -a ${MACHINE_FILE}
	fi
}

# Creates specific 'time files'
log_script()
{
	LOGFILE="et_st$1.txt"
	echo -e "Creating ${BOLD_ORANGE}'${LOGFILE}'${NC} file..."
	ST_TIME_FILE=${LOGS}/${LOGFILE}
	if [[ -f ${ST_TIME_FILE} ]];
	then
		rm ${ST_TIME_FILE}
		touch ${ST_TIME_FILE}
	else
		touch ${ST_TIME_FILE}
	fi
	echo -e ${DONE}
}

# Creates 'DIRs' if needed
mkdir -p ${LOGS}
mkdir -p ${TMPDIR}

#####################################
##### "2) - CREATE TRACK FILES" #####
#####################################

# Creates 'time file'
MASTERLOG="elapsed_time_MASTER.txt"
echo -e "Creating ${BOLD_ORANGE}'${MASTERLOG}'${NC} file..."
TIME_FILE=${LOGS}/${MASTERLOG}
if [[ -f ${TIME_FILE} ]];
then
	echo -e "\tTime file already exists! Appending to it...\n"
else
        touch ${TIME_FILE}
fi
echo -e ${DONE}

# 'Prints' to time file
echo -e "${GREY}==================================================${NC} ${BLUE}${TIMELOG} (${MYUSER})${NC} ${GREY}==================================================${NC}" |& tee -a ${TIME_FILE}
echo -e "EXEC of ${BOLD_GREEN}'$0'${NC} (PID: ${ITALIC_GREEN}${MY_PID}${NC}) master script to execute pipeline in machine ${BOLD_BLUE}'${MACHINE}'${NC}\n" |& tee -a ${TIME_FILE}

# Creates 'machine file'
MASTERMACHINE="machine_usage_MASTER.txt"
echo -e "Creating ${BOLD_ORANGE}'${MASTERMACHINE}'${NC} file..." |& tee -a ${TIME_FILE}
MACHINE_FILE=${LOGS}/${MASTERMACHINE}
if [[ -f ${MACHINE_FILE} ]];
then
	echo -e "\tMachine file already exists! Appending to it..." |& tee -a ${TIME_FILE}
else
        touch ${MACHINE_FILE}
fi
echo -e ${DONE} |& tee -a ${TIME_FILE}

# Prints 'PARs' used
echo -e "${BOLD_GREEN}PARS:${NC}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}WD:${NC} ${WD}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}TRIMMING:${NC} ${TRIM^^}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}GENCODE:${NC} ${GENCODE^^}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}GENOME:${NC} ${GENOME^^}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}TYPE:${NC} ${TYPE}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}LOGDIR:${NC} ${LOGS}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}TMPDIR:${NC} ${TMPDIR}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}MACHINE:${NC} ${MACHINE^^}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}PID (VERSION):${NC} ${MY_PID} (${VERSION})" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}USER:${NC} ${MYUSER^^}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}LOGS:${NC} $(basename ${TIME_FILE}) & $(basename ${MACHINE_FILE})\n" |& tee -a ${TIME_FILE}

# 'First check' of machine usage
echo -e "Performing ${GREY}first${NC} machine check..."
check_machines 1 ${MY_PID}
echo -e ${DONE}

#####################################
##### "3) - EXECUTE SUBSCRIPTS" #####
#####################################

# 'Master FUN'
select_script()
{
	# 'Calls core function'
	# ARGS
	#   1) 'SCRIPT' to invoke
	#   2) name (part) of respective 'LOG' file
	#   3) 'SLEEP' time
	#   Others) 'EXTRA' ARGs required for certain scripts

	#run_script "st0_QC-Trim_MULTIQC.sh" "0" "5s" ${TRIM}  # Runs QC with 'FastQC (v0.11.7) / MultiQC (v1.18)' bef & aft TRIM (if needed) with 'Cutadapt (v1.16)'
	#run_script "st1_KB_Ref.sh" "1" "5s" ${GENCODE} ${GENOME}  # Builds the REF (index) with 'KB wrapper (v0.28.2)', 'GENCODE v36 (TCGA)', and 'HG38'
	run_script "st2_KB_Count."${TYPE}".sh" "2.${TYPE}" "15s" ${TYPE} ${GENCODE} ${GENOME}  # Counts the exp. values with 'KB', 'GENCODE', and 'HG38 as above'
}

# 'Core FUN'
run_script()
{
	# Starts 'local time'
	STARTTIME_LOCAL=$(date +%s)

	# Creates specific 'log file'
	echo -e "EXEC source step ${ITALIC_BLUE}'$1'${NC}...\n" |& tee -a ${TIME_FILE}
	log_script $2

	# Source step 'execution'
	## "$@" passes to script all the parameters passed - MUST use quotes!
	## "${@:4}" passes to script all the parameters passed starting from the FOURTH one - MUST use quotes!
	### 'This is what I want' because the first ARG is the script itself, the second one is the LOG and the third one is the sleep time
	## "${@:4:4}" passes to script up to 4 parameters passed starting from the FOURTH one - MUST use quotes!
	./$1 "${@:4}" |& tee -a ${ST_TIME_FILE} & echo $! > ${TMPDIR}/"myTMP_PID_"$2".txt"  # script gets 'executed' in background 'AND' I get its exact 'PID'!

	# Gets 'PID'
	ST_PID=$(cat ${TMPDIR}/"myTMP_PID_"$2".txt")

	# 'Sleeping' to give time for machine check
	sleep $3

	# Internal check of 'machine usage'
	echo -e "Performing ${GREY}internal${NC} machine check..."
	check_machines 0 ${ST_PID}
	echo -e ${DONE}

	# 'Waits' job to finish
	echo -e "\tWaiting the PID ${ITALIC_GREEN}'${ST_PID}'${NC} to finish..." |& tee -a ${TIME_FILE}
	wait ${ST_PID}  # OK if the 'sleep' time is 'less' than the 'actual' time!
	echo -e ${DONE} |& tee -a ${TIME_FILE}

	# Ends 'local time'
	ENDTIME_LOCAL=$(date +%s)

	# Calculates 'elapsed time'
	LOCAL_ELAPSED_SEC=$((${ENDTIME_LOCAL} - ${STARTTIME_LOCAL}))
	LOCAL_ELAPSED_MIN=$(echo "scale=2; ${LOCAL_ELAPSED_SEC} / 60" | bc -l)
	LOCAL_ELAPSED_HOUR=$(echo "scale=2; ${LOCAL_ELAPSED_MIN} / 60" | bc -l)

	# 'Prints' to time files
	echo -e "EXEC time (LOCAL):\n\t${PURPLE}${LOCAL_ELAPSED_HOUR} hours;\n\t${LOCAL_ELAPSED_MIN} minutes;\n\t${LOCAL_ELAPSED_SEC} seconds${NC}\n" |& tee -a ${TIME_FILE}
	echo -e "EXEC time (LOCAL):\n\t${PURPLE}${LOCAL_ELAPSED_HOUR} hours;\n\t${LOCAL_ELAPSED_MIN} minutes;\n\t${LOCAL_ELAPSED_SEC} seconds${NC}\n" |& tee -a ${ST_TIME_FILE}
}

# Calls 'master script'
select_script

##############################
##### "4) - GRAN FINALE" #####
##############################

# Last check of 'machine usage'
echo -e "Performing ${GREY}last${NC} machine check..."
check_machines 0 ${MY_PID}
echo -e ${DONE}

# Ends 'general time'
ENDTIME_GENERAL=$(date +%s)

# Calculates 'elapsed time'
GENERAL_ELAPSED_SEC=$((${ENDTIME_GENERAL} - ${STARTTIME_GENERAL}))
GENERAL_ELAPSED_MIN=$(echo "scale=2; ${GENERAL_ELAPSED_SEC} / 60" | bc -l)
GENERAL_ELAPSED_HOUR=$(echo "scale=2; ${GENERAL_ELAPSED_MIN} / 60" | bc -l)

# 'Prints' to time file
echo -e "EXEC time (GENERAL):\n\t${BOLD_PURPLE}${GENERAL_ELAPSED_HOUR} hours;\n\t${GENERAL_ELAPSED_MIN} minutes;\n\t${GENERAL_ELAPSED_SEC} seconds${NC}" |& tee -a ${TIME_FILE}
echo -e "${GREY}==================================================${NC} ${BLUE}LOG (${MACHINE})${NC} ${GREY}==================================================${NC}\n" |& tee -a ${TIME_FILE}
