#!/usr/bin/env bash

# "Source step #0": Quality Check (QC) Metrics of SE|PE-FASTQ reads with "FastQC" and, sometimes, adapter SEQ trimming with 'cutadapt'
## BONUS: summarizes the results with a new tool called "MultiQC"!

# 'Prints' message
echo -e "Setting up ENV..."

# Creates 'FUN'
## 'NICE LINK:' https://phoenixnap.com/kb/bash-function

### 'FastQC'
function call_fastqc
{
	# 'ARG 1': R1
	# 'ARG 2': R2
	# 'ARG 3': name
	${TOOLS[0]} \
		-o ${PATH_QC} \
		-f $(echo ${SUFFIX} | cut -f2 -d".") \
		-t ${THREADS} \
		-d ${PATH_TMP} \
		$1 $2 \
		|& tee -a ${PATH_QC_LOGS}/$3"_"${TOOLS[0]^^}".log"
}

### 'MultiQC'
function call_multiqc
{
	# 'ARG 1': before or after trimming
	# 'ARG 2': PATH to DIR with ZIP files from FastQC
	# 'ARG 3': SUFFIX to indicate PATTERN of seeked files
	${TOOLS[2]} \
		-n ${TOOLS[2]}"_Collab-Nina_"$1".txt" \
	        -o ${PATH_MQC} \
	        --ignore-symlinks \
	        -i "MultiQC / Ancient Hominids (Nina Colaboration)" \
	        -b "Summary of all '${TOOLS[0]}' reports $1" \
	        -p \
	        -k "tsv" \
	        -z \
	        --profile-runtime \
	        $2/*$3 \
	        |& tee -a ${PATH_MQC_LOGS}/"Collab-Nina_"$1"_"${TOOLS[2]^^}".log"
}

# My 'OBJs'
WD=$(pwd)
SUFFIX=".fastq.gz"
EXT=".zip"
DONE="\tDone!\n"

# Array of 'tools'
## FastQC --> Metrics
## Cutadapt --> TRIMming
## MultiQC --> Summary
declare -a TOOLS=("fastqc" "cutadapt" "multiqc")

# 'CLI' ARGs
#DATASET=$1  # 'organoides' | Pinar
#GROUP=$2  # 'GSE181405_bRNA-seq_normal' (OrganoidDB / bulk / normal) | GSE141945_bRNA-seq_GBM (OrganoidDB / bulk / GBM) | GSE130238_scRNA-seq (Trujillo-Alysson / single-cell / normal)
#TRIM=$3  # TRIM reads? 'yes' | no

# 'Manual' ARGs
TRIM="no"  # TRIM reads? yes | 'no'

# My 'PATHs'
PATH_FASTQ=${WD}/bulkRNAseq
PATH_QC=${WD}/QC
#PATH_TRIM=${WD}/trim
PATH_MQC=${WD}/MQC
PATH_QC_LOGS=${PATH_QC}/logs
#PATH_TRIM_LOGS=${PATH_TRIM}/logs
PATH_MQC_LOGS=${PATH_MQC}/logs
PATH_TMP=${PATH_QC}/myTMPdir

# Creates 'DIRs' if needed
mkdir -p ${PATH_QC}
mkdir -p ${PATH_MQC}
mkdir -p ${PATH_QC_LOGS}
mkdir -p ${PATH_MQC_LOGS}
mkdir -p ${PATH_TMP}

# 'Fine tune' TOOLs efficiency
THREADS=8  # BEPE server has '32 CPUs'
#THREE_END_TRIM=9  # After analyzing the report of 'MultiQC' (Julia found this cutoff with Pedro)
#FIVE_END_TRIM=-5  # After analyzing the report of 'MultiQC' --> not using it right now

# Gets 'METRICs'
TOTAL=$(ls ${PATH_FASTQ}/*"_R1_001"${SUFFIX} | wc -l)
ITER=0

# 'Prints' message
echo -e ${DONE}

# 'Loops' over samples
## Using 'R1' as REF
for SAMPLE in $(ls ${PATH_FASTQ}/*"_R1_001"${SUFFIX}); do
	# 'Updates' ITER VAR
	((ITER = ${ITER} + 1))

	# Calculates 'PERC'
        PERC=$(echo "scale=2; (${ITER} / ${TOTAL}) * 100" | bc -l)

	# Gets 'filename'
	NAME=$(basename ${SAMPLE} "_R1_001"${SUFFIX})
	echo -e "Sample: ${NAME} [${ITER}/${TOTAL} (${PERC}%)]\n"

	# Obtains 'read names'
	R1=${PATH_FASTQ}/${NAME}"_R1_001"${SUFFIX}
	R2=${PATH_FASTQ}/${NAME}"_R2_001"${SUFFIX}

	# 'Checks' if analyses were already made
	## Using 'R1' as REF again
	if [[ -f ${PATH_QC}/${NAME}"_1_"${TOOLS[0]}".html" && -f ${PATH_QC}/${NAME}"_1_"${TOOLS[0]}${EXT} ]];
	then
		echo -e "\t'${TOOLS[0]^^}' files (HTML & ZIP) already exist. Analysis are not necessary (remove files and run code again for new results)!\n"
	else
		# 'FASTQC': Quality Check (QC) for FASTQ files
		echo -e "\tPerforming the 'QC' analyses on 'PE-FASTQ' files with '${TOOLS[0]}'..."
		call_fastqc ${R1} ${R2} ${NAME}
		echo -e "\t${DONE}"
	fi

	# 'Checks' if we want to trim reads or not
	if [[ ${TRIM} == "yes" ]];
	then
		echo -e "\t'Trimming' SEQs from 'PE-FASTQ' files with '${TOOLS[1]}'..."

		# 'CUTADAPT': adapter trimming
		## 'Using it', despite STAR possibly handling this
		### '-p' to redirect R2 trimmed reads to specific PE-FASTQ-trimmed file
		### PS: '-u' to trim SEQ from R1 FASTQ and '-U' from R2. '+' numbers trim from 3 end and '-' ones from 5 end
		### I might check out later the -z PAR for 'compression level'
		#### 'PS:' Not using also the -u ${FIVE_END_TRIM} and the -U ${FIVE_END_TRIM}
		#### 'PS:' not using the -j ${THREADS} ARG because of ERROR: Running in parallel is not supported on Python 2. 'Waiting Dani to fix it!'
		${TOOLS[1]} \
			-u ${THREE_END_TRIM} -U ${THREE_END_TRIM} \
			-o ${PATH_TRIM}/${NAME}"_R1-trimmed"${SUFFIX} \
			-p ${PATH_TRIM}/${NAME}"_R2-trimmed"${SUFFIX} \
			${R1} ${R2} \
			|& tee -a ${PATH_TRIM_LOGS}/${NAME}"_"${TOOLS[1]^^}".log"
		echo -e "\t${DONE}"

		# Another 'FASTQC' analysis after trimming!
		## Obtains new 'read names'
		R1_TRIM=${PATH_TRIM}/${NAME}"_R1-trimmed"${SUFFIX}
		R2_TRIM=${PATH_TRIM}/${NAME}"_R2-trimmed"${SUFFIX}

		# 'FASTQC': Quality Check (QC) for FASTQ files --> forcefully 'run again' because we are trimming again, so the 'results are necessarily different!'
		echo -e "\tPerforming 'QC' analyses again, now on 'TRIMMED PE-FASTQ' files with '${TOOLS[0]}'..."
		call_fastqc ${R1_TRIM} ${R2_TRIM} ${NAME}"_TRIM"
		echo -e "\t${DONE}"
	else
		echo -e "\tNo trimming needed!\n"
	fi

	# 'Prints' message
	echo -e ${DONE}
done

# 'Checks' if FastQC DIR 'exists' before running MultiQC
echo -e "Checking if '${TOOLS[0]^^}' DIR exists..."
if [[ -d ${PATH_QC} ]];
then  # DIR 'exists'!
	echo -e "${DONE}\nChecking if it is empty..."

	# 'Checks' if FastQC DIR is 'not empty' before running MultiQC
	if [[ -z $(ls -A ${PATH_QC}) ]];
	then  # 'Empty'
		echo -e "\tERROR: Empty! Re-run '${TOOLS[0]}' again...\n"
	else  # 'Not empty'!
		echo -e "${DONE}\nSummarizing '${TOOLS[0]}' results with '${TOOLS[2]}'...\n"

		# 'MULTIQC': aggregates results from bioinfo analyses across many samples into a single report!
		## Take care regarding 'symlinks'!
		### -l  File with list of PATHs (1/row)
		### -x  Ignore files (we can use patterns, but between double quotes!) --> it is better to specify which ones you do want to analyze (ex: .zip from FastQC)
		### 'DANGER' --pdf      PDF report (requires Pandoc to be installed)
		### Explore later the Sample handling & MultiQC behaviour OPTs
		call_multiqc "befTRIM" ${PATH_QC} "_"${TOOLS[0]}${EXT}

		# Runs it again only if 'cutadapt OPT is ON'
		if [[ ${TRIM} == "yes" ]];
		then
			call_multiqc "aftTRIM" ${PATH_QC} "-trimmed_"${TOOLS[0]}${EXT}
		else
			:
		fi

		# 'Prints' message
		echo -e ${DONE}
	fi
else  # DIR 'does not exist'
	echo -e "\tERROR: Non-existent DIR! Check script again...\n"
fi
