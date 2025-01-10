#!/usr/bin/env bash

###########################################################################################################################################################
# Author: Filipe Ferreira dos Santos
# Date: 03/22/2024

# This shell script aims to execute all "source_stepXXX.sh" scripts for the differential expression pipeline
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
COMPARISON=$1  # comparison from SHELL SCRIPT (PVT_vs_Shaker | Others) --> 'USED TO CHANGE HERE!'
BATCH=$2  # batch from SHELL SCRIPT (1 | 2 | both) --> 'USED TO CHANGE HERE!'

# My 'OBJs'
WD=$(pwd)
LOGS=${WD}/logs
TMPDIR=${WD}/myTMPdir
MACHINE=${HOSTNAME}
MY_PID=$$
DONE="\t${BG_GREEN}Done!${NC}\n"

# Sets 'MYUSER'
MYUSER=${USER}

# Script 'VERSION'
VERSION="v3.0"

# Displays 'how to use' the script '(HELP)'
display_usage()
{
        USAGE="${BOLD_BLUE}Usage:${NC} $0 ${BOLD_ORANGE}[OPTIONS]${NC} --> program to execute internal HSL/Pedro A.F. Galante differential expression (DEG) pipeline"
        echo -e "${GREY}===========================================================================================================================================================${NC}\n${USAGE}\n"
        echo -e "This shell script (${GREEN}${VERSION}${NC}) performs automatically our internal differential expression (DEG) pipeline (TXimport + DESeq2) --> user ${BOLD_GREEN}'${MYUSER}'${NC}."
        echo -e "${ORANGE}Options:${NC}\n\t${GREEN}-v, --version${NC}\tdisplays this script version\n\t${GREEN}-h, --help${NC}\tdisplays this help message"
        echo -e "\t${GREEN}-l, --license${NC}\tdisplays the license for this script\n"
        echo -e "${BOLD_PURPLE}Contact${NC}\nEmail: fidossantos@health.ucsd.edu\nPhone: +1 (858) 933-6465\n"
        echo -e "${BOLD_RED}License${NC}\nGNU General Public License (GPLv3). For more detailed information, please check out this link: ${BOLD_BLUE}https://www.gnu.org/licenses/gpl-3.0.html${NC}"
        echo -e "or run this script with the '-l, --license' ARG option."
        echo -e "${GREY}===========================================================================================================================================================${NC}\n"
}

# Checks 'ARGs' provided by USER
parse_args()
{
	# Checks if user used the 'help' option
	if [[ ${#} -eq 1 && ( ${1} == "--help" || ${1} == "-h" ) ]];
	then
		display_usage
		exit 0
	else
		# Checks if user passed the 'version' option
		if [[ ${#} -eq 1 && ( ${1} == "--version" || ${1} == "-v" ) ]];
		then
			echo -e "Script ${GREEN}'${0}'${NC} version: ${BLUE}${VERSION}${NC}"
			exit 0
		else
			# Checks if user passed the 'license' option
			if [[ ${#} -eq 1 && ( ${1} == "--license" || ${1} == "-l" ) ]];
			then
				cat ${WD}/../../license/"GNU_GPLv3.txt"
				exit 0
			else
				:
			fi
		fi
	fi
}

# Calls the 'parse_args' function
case ${#} in
	1)
		parse_args ${1};;
	0)
		parse_args ${@}
esac

######################################
##### "1) - SETUP MASTER SCRIPT" #####
######################################

# Starts 'GENERAL time'
STARTTIME_GENERAL=$(date +%s)

# My 'FUNs'
# How to https://stackoverflow.com/questions/10822790/can-i-call-a-function-of-a-shell-script-from-another-shell-script
check_machines()
{
	if [[ $1 -eq 1 ]];
	then
		echo -e "\n============================== ${BLUE}${TIMELOG}${NC} ==============================\n" |& tee -a ${MACHINE_FILE}
		top -b -p $2 -n1 | tail -n2 |& tee -a ${MACHINE_FILE}
	else
		top -b -p $2 -n1 | tail -n1 |& tee -a ${MACHINE_FILE}
	fi
}

# Creates specific 'time files'
log_script()
{
	LOGFILE="et_st$1-DEG_${BATCH}-${COMPARISON}.txt"
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
MASTERLOG="elapsed_time_sourceSteps-MASTER.txt"
echo -e "Creating ${BOLD_ORANGE}'${MASTERLOG}'${NC} file..."
TIME_FILE=${LOGS}/${MASTERLOG}
if [[ -f ${TIME_FILE} ]];
then
	echo -e "\tTime file already exists! Appending to it...\n"
else
        touch ${TIME_FILE}
	echo -e ${DONE}
fi

# 'Prints' to time file
echo -e "\n${GREY}==================================================${NC} ${BLUE}${TIMELOG}${NC} ${GREY}==================================================${NC}" |& tee -a ${TIME_FILE}
echo -e "Execution of ${BOLD_GREEN}'$0'${NC} (PID: ${ITALIC_GREEN}${MY_PID}${NC}) master script to execute pipeline in machine ${BOLD_BLUE}'${MACHINE^^}'${NC}:\n" |& tee -a ${TIME_FILE}

# Creates 'machine file'
MASTERMACHINE="machine_usage_sourceSteps-MASTER.txt"
echo -e "Creating ${BOLD_ORANGE}'${MASTERMACHINE}'${NC} file..." |& tee -a ${TIME_FILE}
MACHINE_FILE=${LOGS}/${MASTERMACHINE}
if [[ -f ${MACHINE_FILE} ]];
then
	echo -e "\tMachine file already exists! Appending to it...\n" |& tee -a ${TIME_FILE}
else
        touch ${MACHINE_FILE}
	echo -e ${DONE} |& tee -a ${TIME_FILE}
fi

# Prints 'PARs' used
echo -e "${BOLD_GREEN}PARS:${NC}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}WD:${NC} ${WD}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}TEST:${NC} ${UND_BLUE}${BATCH}${NC} (${ITALIC_BLUE}${COMPARISON}${NC})${NC}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}LOGDIR:${NC} ${LOGS}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}TMPDIR:${NC} ${TMPDIR}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}MACHINE:${NC} ${MACHINE^^}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}PID:${NC} ${MY_PID}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}USER:${NC} ${MYUSER^^}" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}LOG:${NC} $(basename ${TIME_FILE})" |& tee -a ${TIME_FILE}
echo -e "\t${ORANGE}MACHINE LOG:${NC} $(basename ${MACHINE_FILE})\n" |& tee -a ${TIME_FILE}

# 'First check' of machine usage
echo -e "Performing ${GREY}first${NC} machine check..."
check_machines 1 ${MY_PID}
echo -e ${DONE}

#####################################
##### "3) - EXECUTE SUBSCRIPTS" #####
#####################################

# 'Master' function
select_script()
{
	# Calls core 'run_script' function
	# 'ARGS'
	#   1) 'SCRIPT' to invoke
	#   2) name (part) of respective 'LOG' file
	#   3) 'SLEEP' time
	#   4-5) 'ARGs' to script (comparison & batch)

	run_script "source_step1_TXimport.sh" "1" "15s" ${COMPARISON} ${BATCH}  # Runs 'TXimport' to get DEG values (COUNTS / TPM) by Gene (ENST2ENSG)
	run_script "source_step2_DESeq2.sh" "2" "15s" ${COMPARISON} ${BATCH}  # Runs 'DESeq2' to get DEG genes across conditions
}

# 'Core' function
run_script()
{
	# Starts 'LOCAL time'
	STARTTIME_LOCAL=$(date +%s)

	# Creates specific 'LOG file'
	echo -e "Executing source step ${ITALIC_BLUE}'$1'${NC}..." |& tee -a ${TIME_FILE}
	log_script $2

	# Source step 'execution'
	## Consider removing this '&' to avoid the necessity of using the sleep CMD and not worrying about one CMD running at the same time of others... The thing is the machine usage
	## The 'wait' CMD should do the trick, but maybe I am not getting the PID correctly... Perhaps I should use the '$!'
	### 'BEST solution:' add 'echo $!' right after the CMD and 'redirects' the output (the PID of that specific process) to 'file' (redirecting to object causes the script to freeze)
	./$1 $4 $5 |& tee -a ${ST_TIME_FILE} & echo $! > ${TMPDIR}/"myTMP_PID_"$2".txt"  # script gets 'executed' in background 'AND' I get its exact 'PID'!

	# Gets 'PID'
	ST_PID=$(cat ${TMPDIR}/"myTMP_PID_"$2".txt")

	# 'Sleeping' to give time for machine check
	sleep $3

	# Internal check of 'machine usage'
	echo -e "Performing ${GREY}internal${NC} machine check..."
	check_machines 0 ${ST_PID}
	echo -e ${DONE}

	# 'Waits' job to finish
	## So, it is OK if the 'sleep' time is 'less' than the 'actual' time!
	echo -e "Waiting the PID ${BOLD_GREEN}'${ST_PID}'${NC} to finish..." |& tee -a ${TIME_FILE}
	wait ${ST_PID}
	echo -e ${DONE} |& tee -a ${TIME_FILE}

	# Ends 'LOCAL time'
	ENDTIME_LOCAL=$(date +%s)

	# Calculates 'ELAPSED time'
	LOCAL_ELAPSED_SEC=$((${ENDTIME_LOCAL} - ${STARTTIME_LOCAL}))
	LOCAL_ELAPSED_MIN=$(echo "scale=2; ${LOCAL_ELAPSED_SEC} / 60" | bc -l)
	LOCAL_ELAPSED_HOUR=$(echo "scale=2; ${LOCAL_ELAPSED_MIN} / 60" | bc -l)

	# 'Prints' to time files
	echo -e "\tEXEC time (LOCAL):\n\t${PURPLE}${LOCAL_ELAPSED_HOUR} hours;\n\t${LOCAL_ELAPSED_MIN} minutes;\n\t${LOCAL_ELAPSED_SEC} seconds${NC}\n" |& tee -a ${TIME_FILE}
	echo -e "\tEXEC time (LOCAL):\n\t${PURPLE}${LOCAL_ELAPSED_HOUR} hours;\n\t${LOCAL_ELAPSED_MIN} minutes;\n\t${LOCAL_ELAPSED_SEC} seconds${NC}\n" |& tee -a ${ST_TIME_FILE}
}

# Calls master 'select_script' function
select_script

##############################
##### "4) - GRAN FINALE" #####
##############################

# Last check of 'machine usage'
echo -e "Performing ${GREY}last${NC} machine check..."
check_machines 0 ${MY_PID}
echo -e ${DONE}

# Ends 'GENERAL time'
ENDTIME_GENERAL=$(date +%s)

# Calculates 'ELAPSED time'
GENERAL_ELAPSED_SEC=$((${ENDTIME_GENERAL} - ${STARTTIME_GENERAL}))
GENERAL_ELAPSED_MIN=$(echo "scale=2; ${GENERAL_ELAPSED_SEC} / 60" | bc -l)
GENERAL_ELAPSED_HOUR=$(echo "scale=2; ${GENERAL_ELAPSED_MIN} / 60" | bc -l)

# 'Prints' to time file
echo -e "EXEC time (GENERAL):\n${BOLD_PURPLE}${GENERAL_ELAPSED_HOUR} hours;\n${GENERAL_ELAPSED_MIN} minutes;\n${GENERAL_ELAPSED_SEC} seconds${NC}" |& tee -a ${TIME_FILE}
echo -e "${GREY}==================================================${NC} ${BLUE}LOG (${MACHINE})${NC} ${GREY}==================================================${NC}\n" |& tee -a ${TIME_FILE}
