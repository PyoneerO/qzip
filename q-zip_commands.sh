#!/bin/bash

#WORKAROUND: DEFINE A FUNCTION FOR THE WORKFLOW TO HAVE EASY TERMINAL LOGGING
function q-zip_run {

############
### TIME ###
############

ST_DATE=`date +%Y-%m-%d`
ST_TIME=`date +%H-%M-%S`
echo "Analysis started at "${ST_DATE} ${ST_TIME}
DUR_SEC_ST=${SECONDS}

##################
### Parameters ###
##################

### GET PARAMETERS
echo "Load options from q-zip_parameters.txt"
sleep 1
. ./q-zip_parameters.txt
cat q-zip_parameters.txt

### PROJECT INFO SUMMARY
echo "Project name: "${PROJECT_NAME}
echo "Project ID: "${PROJECT_ID}
echo "PIs: "${PI_NAMES}
echo "Target: "${TARGET_MOLECULE}" "${TARGET_REGION}
echo "Reference(s): "${REF_DBS}

#########################
### Check environment ###
#########################

echo "Check availability of needed tools"
### EXIT IF IMPORTANT TOOLS NOT AVAILABLE
which ${FASTQC} > /dev/null; if [ $? -ne 0 ]; then echo "FASTQC not found; exit"; exit 1; else echo "FASTQC found"; fi
which ${TRIMMOMATIC}  > /dev/null; if [ $? -ne 0 ]; then echo "TRIMMOMATIC not found; exit"; exit 1; else echo "TRIMMOMATIC found"; fi
which ${CUTADAPT} > /dev/null; if [ $? -ne 0 ]; then echo "CUTADAPT not found; exit"; exit 1; else echo "CUTADAPT found"; fi
which ${SWARM} > /dev/null; if [ $? -ne 0 ]; then echo "SWARM not found; exit"; exit 1; else echo "SWARM found"; fi
which ${VSEARCH} > /dev/null; if [ $? -ne 0 ]; then echo "VSEARCH not found; exit"; exit 1; else echo "VSEARCH found"; fi
which ${MOTHUR} > /dev/null; if [ $? -ne 0 ]; then echo "MOTHUR not found; exit"; exit 1; else echo "MOTHUR found"; fi
which ${BIOM} > /dev/null; if [ $? -ne 0 ]; then echo "BIOM-format not found; exit"; exit 1; else echo "BIOM-format found"; fi
which parallel > /dev/null; if [ $? -ne 0 ]; then echo "GNU PARALLEL not found; exit"; exit 1; else echo "GNU PARALLEL found"; fi

echo "Check max number of open file desriptors"
### EXIT IF NUMBER OF FILES EXCEED MAX NUMBER OF OPEN FILE DESCRIPTORS
if [[ `ulimit -n` -lt `ls -l1 | grep -c _R1_001.fastq*` ]]; then
	echo "ERROR: max number of open file descriptors smaller than number of samples"
	echo `ulimit -n` "<" `ls -l1 | grep -c _R1_001.fastq*`
	echo "Please ask your sys admin to increase the number to be at least "`ls -l1 | grep -c _R1_001.fastq*`
	exit 1
else
	echo "OK: max number of open file descriptors greater than number of samples"
	echo `ulimit -n` ">" `ls -l1 | grep -c _R1_001.fastq*`
fi

################
### WORKFLOW ###
################


## CLEAN AND CREATE NEEDED DIRECTORIES
rm -rf ${LOG_FP} ${RESULTS_DIR} ${S_STRUCT} ${S_STATS} ${S_SEEDS}* ${S_SWARM} ${PROJECT_NAME}*".map" \
*"trimmed"* *OTU_table* "full_set_dereplicated.fasta" "seq_number_stats.txt" *"wang"* q-zip_seq_of_coms.txt  \
`basename ${OTU_TABLE} .csv`* ${AMPLICON_TABLE} ${REF_USED_FP} ${REF_USED_FP}
mkdir ${LOG_FP} ${RESULTS_DIR} ${REF_USED_FP}
#if [ ! -e "${REF_USED_FP}" ]; then mkdir ${REF_USED_FP}; fi
#if [ "${REUSE_REF_SEQS}" == "NO" ]; then rm -r ${REF_USED_FP}; mkdir ${REF_USED_FP}; fi
touch q-zip_seq_of_coms.txt

## PREPARE REPRODUCIBILITY VIA TRACE FILE
echo "#PROJECT ID: "${PROJECT_ID} >> q-zip_seq_of_coms.txt
echo -e "\n#source parameter" >> q-zip_seq_of_coms.txt
echo ". ./q-zip_parameters.txt" >> q-zip_seq_of_coms.txt
echo -e "\n#Sequence of commands for current project:" >> q-zip_seq_of_coms.txt
echo "PROJECT_ID="$PROJECT_ID >> q-zip_seq_of_coms.txt
echo -e "\n#create needed directory for reference sequences" >> q-zip_seq_of_coms.txt
echo "mkdir "${REF_USED_FP} >> q-zip_seq_of_coms.txt


## CREATE PRIMER-TRIMMED REFDB
#for i in ${REF_DBS}; do
#searching and trimming of primer of reference sequences
#create command
cmd='cat ${REF_RAW_PATH}"/"${REF_DBS}".fasta" |
${CUTADAPT} -g ${FORWARDPRIMER} --discard-untrimmed --minimum-length ${MIN_LEN_REF} -e ${PRIMER_MISMATCH_REF} -O ${lenFP_CUT_REF} - |
${CUTADAPT} -a ${REVERSEPRIMER_RC} --discard-untrimmed --minimum-length ${MIN_LEN_REF} --maximum-length ${MAX_LEN_REF} -e ${PRIMER_MISMATCH_REF} -O ${lenRP_CUT_REF} -
> ${REF_USED_FP}"/"${REF_DBS}"_"${FORWARDPRIMER}"_"${REVERSEPRIMER_RC}"_"${PRIMER_MISMATCH_REF}"_"${lenFP_CUT_REF}"_"${lenRP_CUT_REF}".fasta"'
#execute command
eval $cmd
#write command to record file
echo -e "\n#Search for primers in reference and trimming of primers:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt
#create corresponding taxonomy map
#create command
cmd='grep ">" ${REF_USED_FP}"/"${REF_DBS}"_"${FORWARDPRIMER}"_"${REVERSEPRIMER_RC}"_"${PRIMER_MISMATCH_REF}"_"${lenFP_CUT_REF}"_"${lenRP_CUT_REF}".fasta" |
tr -d "^>" | awk '\''FNR==NR{a[$1]=$1; next}($1 in a){print $0}'\'' - ${REF_RAW_PATH}"/"${REF_DBS}".tax"
> ${REF_USED_FP}"/"${REF_DBS}"_"${FORWARDPRIMER}"_"${REVERSEPRIMER_RC}"_"${PRIMER_MISMATCH_REF}"_"${lenFP_CUT_REF}"_"${lenRP_CUT_REF}".tax"'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#prepare corresponding taxonomy map:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt
#done

## CREATE RAW FILES TO SAMPLE MAPPING
# create command
cmd='ls -1 | grep "001.fastq.gz$\|001.fastq$" |
awk '\''{if ($0 ~ /_R1_/) {printf($0); gsub(/_L00[1-4]_R1_001.fastq.gz/,"",$0); gsub(/_L00[1-4]_R1_001.fastq/,"",$0) ; printf("\t"$0"\t")} else {print $0}}'\'' |
awk '\''{print($2"\t"$1"\t"$3)}'\'' > ${PROJECT_ID}".map"'
# execute command
eval $cmd
#write command to command record file
echo -e "\n#prepare sample name mapping file:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

## GET SAMPLE NAMES FOR STATISTICS
cut -f1 ${PROJECT_ID}*".map" > seq_number_stats.interm
## COUNT SEQUENCES IN RAW FILES FOR STATISTICS
cat ${PROJECT_ID}*".map" | cut -f2 | parallel -j ${THREADS} -k "if [ -s {} ]; then zgrep -c '^+$' {}; else echo '-'; fi " >> seq_number_stats.interm

## QUAL TRIM SEQUENCES
# create command
cmd='cut -f 1 ${PROJECT_ID}".map" |  grep -v "^#" |tr "." "_" |
parallel -j ${THREADS} java -jar ${TRIMMOMATIC}" PE -phred33 {}_L001_R1_001.fastq* {}_L001_R2_001.fastq* {.}.trimmed.R1 /dev/null {.}.trimmed.R2 /dev/null CROP:"${CROP}" SLIDINGWINDOW:"${SLIDINGWINDOW}'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#3'-end trimming:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt
## GET SEQ NUMBER AFTER TRIMMING
cat ${PROJECT_ID}*".map" | while read a b c ; do if [ -s $a.trimmed.R1 ]; then grep -c "^+$" $a.trimmed.R1; else echo "0"; fi;done >> seq_number_stats.interm

## MERGE PAIRED ENDS
#create command
cmd='ls -S1 . | grep "trimmed.R1" |
parallel -j ${THREADS} ${VSEARCH}" --fastq_mergepairs {} --reverse {.}'.R2' --log - --fastq_allowmergestagger --threads 1 --fasta_width 0 --fastq_maxdiffs "${FASTQ_MAXDIFFS}" --fastq_minovlen "${FASTQ_MINOVLEN}" --fastqout {.}.assembled.fastq"'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#paired-end merging:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

## REMOVE NOT NEEDED INTERMEDIATE FILES
if [ ${DEBUG} == "NO" ]; then rm *trimmed.R1 *trimmed.R2; fi
## GET SEQ NUMBER AFTER MERGING/ASSEMBLING
cat ${PROJECT_ID}*".map" | while read a b c ; do if [ -s $a.trimmed.assembled.fastq ]; then grep -c "^+$" $a.trimmed.assembled.fastq; else echo "0"; fi;done >> seq_number_stats.interm

## CREATE REVERSE COMPLEMENTED SEQUENCES
#create command
cmd='ls -S1 . | grep ".assembled.fastq" |
parallel -j ${THREADS} ${VSEARCH}" --fastx_revcomp {} --threads 1 --fastq_ascii 33 --fasta_width 0 --log - --fastqout {.}.revcomp.fastq"'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#create reverse complement:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

## MERGE BOTH DIRECTION
#create command
cmd='ls -1 . | grep ".assembled.fastq" |
awk '\''{gsub("fastq","",$1); printf($1"bothdir_concat.fastq\t"$1"fastq\t"$1"revcomp.fastq\n")}'\'' |
while read a b c ; do cat $b $c > $a ; done'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#concatenate both directions:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

## REMOVE NOT NEEDED INTERMEDIATE FILES
if [ ${DEBUG} == "NO" ]; then rm *assembled.fastq *assembled.revcomp.fastq; fi

## FILTER BY EXISTENENCE OF BOTH PRIMERS AND TRUNCATE PRIMER SEQUENCES
#create command
cmd='ls -S1 | grep "bothdir_concat.fastq" |
parallel -j ${THREADS} "cat {} | "${CUTADAPT}" -g "${FORWARDPRIMER}" -e "${PRIMER_MISMATCH}" -O "${lenFP_CUT}" --discard-untrimmed - | "${CUTADAPT}" -a "${REVERSEPRIMER_RC}" -e "${PRIMER_MISMATCH}" -O "${lenRP_CUT}" --discard-untrimmed -
> {.}.primer_cut.fastq"'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#filter and trim primer:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


## REMOVE NOT NEEDED INTERMEDIATE FILES
if [ ${DEBUG} == "NO" ]; then rm *bothdir_concat.fastq; fi
## GET SEQ NUMBER AFTER PRIMER FILTERING
cat ${PROJECT_ID}*".map" | while read a b c ; do if [ -s $a.trimmed.assembled.bothdir_concat.primer_cut.fastq ]; then grep -c "^+$" $a.trimmed.assembled.bothdir_concat.primer_cut.fastq; else echo "0"; fi;done >> seq_number_stats.interm


## FEATURE FILTERING AND SHA1-RELABELing
#create command
cmd='ls -S1 | grep "primer_cut.fastq" |
parallel -j ${THREADS} ${VSEARCH}" --threads 1 --fastq_ascii 33 --log - --fasta_width 0 --fastx_filter {} --fastq_maxee "${FASTQ_MAXEE}" --fastq_maxlen "${FASTQ_MAXLEN}" --fastq_minlen "${FASTQ_MINLEN}" --fastq_maxns "${FASTQ_MAXNS}" --relabel_sha1 --relabel_keep --fastqout {.}.feature_filtered.fastq"'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#Feature filtering and sha1 relabeling:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


## REMOVE NOT NEEDED INTERMEDIATE FILES
if [ ${DEBUG} == "NO" ]; then rm *primer_cut.fastq; fi
## GET SEQ NUMBER AFTER FEATURE FILTERING
cat ${PROJECT_ID}*".map" | while read a b c ; do if [ -s $a.trimmed.assembled.bothdir_concat.primer_cut.feature_filtered.fastq ]; then grep -c "^+$" $a.trimmed.assembled.bothdir_concat.primer_cut.feature_filtered.fastq; else echo "0"; fi;done >> seq_number_stats.interm

### MAKE FASTQC REPORTS
## subsample feature_filtered seqs (FASTQC_SUB)
# ls -1S | grep feature_filtered.fastq | parallel -j ${THREADS} "awk -v fqs="${FASTQC_SUB}" 'NR%(fqs*4)==1{print;getline;print;getline;print;getline;print}' {} > {.}.sub_${FASTQC_SUB}.fastq"
# ${FASTQC} -t ${THREADS} *sub_${FASTQC_SUB}.fastq
# only subsampled
# Basic Statistics => PASS;
# Per base sequence quality => PASS;
# Per sequence quality scores =>;
# Per base N content => Pass;
# Sequence length distribution => [PASS | WARN]

## DEREPLICATION ON SAMPLE LEVEL
#create command
cmd='ls -S1 | grep "feature_filtered.fastq" |
parallel -j ${THREADS} ${VSEARCH}" --derep_fulllength {} --threads 1 --log - --fasta_width 0 --relabel_sha1 --sizeout --output {.}.derep.fasta"'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#dereplication on sample level:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


## GET SEQ NUMBER AFTER DEREPLICATION
cat ${PROJECT_ID}*".map" | while read a b c ; do if [ -s $a.trimmed.assembled.bothdir_concat.primer_cut.feature_filtered.derep.fasta ]; then grep -c "^>" $a.trimmed.assembled.bothdir_concat.primer_cut.feature_filtered.derep.fasta; else echo "0"; fi;done >> seq_number_stats.interm
## REMOVE NOT NEEDED INTERMEDIATE FILES
if [ ${DEBUG} == "NO" ]; then rm *feature_filtered.fastq; fi

## DETECT AND REMOVE CHIMERAS DENOVO
#Detect chimera
# create command
cmd='ls -1S | grep "derep.fasta" |
parallel -j ${THREADS} ${VSEARCH}" --fasta_width 0 --threads 1 --log - --log - --uchime_denovo {} --fasta_score --chimeras {.}.chimeras_denovo.fasta"'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#detect chimera:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt
#consider potential chimera only if they only occur in one sample (blacklist)
# create command
cmd='awk -F";" '\''NR%2==1 {gsub(">","");print $1}'\'' *chimeras_denovo.fasta |
sort -n | uniq -c | sort -nr | awk '\''{if ($1==1) print $2}'\'' > chimera.list'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#create chimera black list:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt
#remove chimeric sequences from fastas
#create command
cmd='ls -1S | grep "derep.fasta" |
parallel -j ${THREADS} "awk '\''FNR==NR{a[\$1]=\$1; next}(!(substr(\$0,2,index(\$0,\";\")-2) in a) && (\$0 ~ />/)){if (\$0 ~ /^>/) print \$0; getline; print}'\'' chimera.list {}
> {.}.non_chimeras_denovo.fasta"'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#filter chimera from fastas files:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


## REMOVE NOT NEEDED INTERMEDIATE FILES
if [ ${DEBUG} == "NO" ]; then rm *derep.fasta; rm *derep.chimeras_denovo.fasta; fi
# GET SEQ NUMBER AFTER CHIMERA REMOVAL
cat ${PROJECT_ID}*".map" | while read a b c ; do if [ -s $a.trimmed.assembled.bothdir_concat.primer_cut.feature_filtered.derep.non_chimeras_denovo.fasta ]; then grep -c "^>" $a.trimmed.assembled.bothdir_concat.primer_cut.feature_filtered.derep.non_chimeras_denovo.fasta; else echo "0"; fi;done >> seq_number_stats.interm
# GET SEQ NUMBER AFTER CHIMERA REMOVAL, IF REREPLICATED (added abundance of each unified/dereplicated sequence)
cat ${PROJECT_ID}*".map" | while read a b c ; do if [ -s $a.trimmed.assembled.bothdir_concat.primer_cut.feature_filtered.derep.non_chimeras_denovo.fasta ]; \
then awk -F "[;=]" 'BEGIN{cnt=0}{cnt+=$3}END{print cnt}' $a.trimmed.assembled.bothdir_concat.primer_cut.feature_filtered.derep.non_chimeras_denovo.fasta ; else echo "0"; fi;done >> seq_number_stats.interm
# GET AVERAGE LENGTH OF USED SEQUENCES (...non_chimeras_denovo.fasta)
cat ${PROJECT_ID}*".map" | while read a b c ; do if [ -s $a.trimmed.assembled.bothdir_concat.primer_cut.feature_filtered.derep.non_chimeras_denovo.fasta ]; then \
awk '{gsub(";size=","_"); gsub(/;$/,""); print}' $a.trimmed.assembled.bothdir_concat.primer_cut.feature_filtered.derep.non_chimeras_denovo.fasta | \
awk -F"_" 'BEGIN{seqnr=0}NR%2==1{seqnr=$2;getline; print(seqnr*length($0), seqnr)}' | \
awk 'BEGIN{leng=0;seqnr=0}{leng+=$1;seqnr+=$2}END{print(leng/seqnr)}'; else echo "-"; fi;done >> seq_number_stats.interm

## REMOVE FASTAS WITH TOO FEW SEQS
# create command
cmd='for i in *"${WORKFLOW_SUFFIX}.fasta"; do cnt=`awk -F"[;=]" '\''BEGIN{cnt=0}{cnt+=$3}END{print cnt}'\'' $i`; 
if [[ $cnt -lt ${MIN_SAMPLE_SIZE} ]]; then echo "remove "$i", too few sequences ("$cnt")"; rm $i; fi; done'
#execute command
eval $cmd
#write sequence of commands to record file
echo -e "\n#Min sample size filter" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


## ADD SAMPLE ABUNDANCE FILTER PASSED TO SEQ NUMBER STATS
cat ${PROJECT_ID}*".map" | while read a b c ; do if [ -s $a*${WORKFLOW_SUFFIX}.fasta ]; then \
echo "yes"; else echo "no"; fi; done >> seq_number_stats.interm

## DEREPLICATE FULL STUDY + ADJUST ABUNDANCE INFORMATION FORMAT
#create command
#cmd='cat *.non_chimeras_denovo.fasta | ${VSEARCH} --threads ${THREADS} --derep_fulllength - --sizein --sizeout --fasta_width 0 --output full_set_dereplicated.fasta'
cmd='cat *.non_chimeras_denovo.fasta |
${VSEARCH} --threads ${THREADS} --derep_fulllength - --sizein --sizeout --log full_derep.log --fasta_width 0 --output - |
sed '\''s/;size=/_/; s/;$//'\'' > full_set_dereplicated.fasta'
#execute command
eval $cmd; cat full_derep.log
#write command to command trace file
echo -e "\n#Dereplicate full study + adjust abundance information:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


## ADD HEADER TO STATS FILE
echo -e "sample_ID\traw\ttrimmed\tassembled\tprimer_filtered\tfeature_filtered\tsample_derep\tchimera_filtered\tfinal_rerep\tavg_length\tsample_passed_abundance_filter" > seq_number_stats.txt
## FORMAT INTERMEDIATE STATS FILE INTO COLUMNS
lns=`awk 'END{print NR}' ${PROJECT_ID}*".map"`; cat seq_number_stats.interm | pr -ts'' --columns 11 -l $lns >>seq_number_stats.txt
## PRINT STATS FILE
column -t seq_number_stats.txt

## SWARM OTU CLUSTERING
#create command
cmd='cat full_set_dereplicated.fasta | ${SWARM} -t ${THREADS} -d ${DISTANCE} -f -b ${F_BOUNDARY} -i ${S_STRUCT} -s ${S_STATS} -w ${S_SEEDS} > ${S_SWARM}'
#execute command
eval $cmd
#write command to command trace file
echo -e "\n#Do swarming:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


## FILTER SINGLETONS (seeds, swarm and stats file; create amplicon names file)
#create command
cmd='awk -F"_" '\''{if (($0 ~ /^>/) && ($2 >= 2)) {print;getline;print} }'\'' ${S_SEEDS} > ${S_SEEDS}".no.singletons";
awk '\''{if ( NF>1 || (match($0,"_1$") == zero)) print}'\'' ${S_SWARM} > ${S_SWARM}.no.singletons;
awk '\''{if ($2>1) print}'\'' ${S_STATS} > ${S_STATS}.no.singletons;
awk '\''{gsub(" ","\n");print}'\'' ${S_SWARM}.no.singletons | awk -F"_" '\''{print $1}'\''  > swarm.no.singletons'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#filter singletons from seeds, swarm and stats file and create an amplicon names file:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


## CREATE AMPLICON TABLE
# create amplicon tables sample-wise
# create command
cmd='ls -1S | grep "${WORKFLOW_SUFFIX}.fasta" | parallel -j ${THREADS} "echo {} |
awk '\''{sub(\"${WORKFLOW_SUFFIX}.fasta\",\"\");print}'\'' > {.}.amplicon;
awk '\''/^>/{gsub(/;size=/,\"\t\");gsub(\";\",\"\");gsub(\"^>\",\"\"); print}'\'' {} |
awk '\''FNR==NR{a[\$1]=\$2;next} {if(a[\$1]!=\"\") print(a[\$1]); if (a[\$1]==\"\") print(\"0\") }'\'' - swarm.no.singletons
>> {.}.amplicon"'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#Create amplicon tables for each sample/file" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

# add header to amplicon list
#create command
cmd='sed  -i '\''1i amplicon'\'' swarm.no.singletons'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#add header to amplicon names list:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

# create list of files to be merged column-wise (amplicon names file + amplicon abundance files per sample)
#create command
cmd='slist="swarm.no.singletons "`ls -1 | grep ${WORKFLOW_SUFFIX}".amplicon"`'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#create list of files to be merged:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

# merge amplicon names and all amplicon abundances to form an amplicon contingency table by the tool 'paste'; sumup rows
# create command
cmd='paste $slist |
parallel -j ${THREADS} -k -q --pipe awk '\''{if ($0 ~ /^amplicon/) {print} else {cnt=0; printf($1); for (i=2;i<=NF;i++) {cnt+=$i; printf("\t"$i)};printf("\t"cnt);printf("\n")}}'\''
> ${AMPLICON_TABLE}"_unsorted"'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#merge amplicon names and all amplicon abundances to form an amplicon contingency table and sum-up rows (totalamplicon numbers):" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

# complete header ("total" for sum-up column)
#create command
cmd='sed -i '\''1s/$/\ttotal/'\'' ${AMPLICON_TABLE}"_unsorted"'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#complete header of amplicon table:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

# Sorting of amplicon table
#create command
cmd='head -1 ${AMPLICON_TABLE}"_unsorted" > ${AMPLICON_TABLE};
NUM_FIELDS=`head -2 ${AMPLICON_TABLE}"_unsorted" | tail -1 | awk '\''{print NF}'\''`;
LC_ALL=C; awk '\''NR > 1'\'' ${AMPLICON_TABLE}"_unsorted" | sort -T . -k${NUM_FIELDS},${NUM_FIELDS}nr -k1,1d >> ${AMPLICON_TABLE}'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#Create sorted version of amplicon table:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


## REMOVE NOT NEEDED INTERMEDIATE FILES
if [ ${DEBUG} == "NO" ]; then rm *non_chimeras_denovo.fasta *non_chimeras_denovo.amplicon; fi

## CREATE OTU TABLE
#header of OTU table
# create command
cmd='echo -e "#OTU\t$(head -n 1 "${AMPLICON_TABLE}")" > ${OTU_TABLE}'
# execute command
eval $cmd
#print command to command trace file
echo -e "\n#Create header for OTU table:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt
# OTU table
# create command
cmd='
cat ${S_STATS}.no.singletons | parallel -j ${THREADS} --pipe -l --block-size 100 --round-robin -q 
awk -v SWARM="${S_SWARM}.no.singletons"
    -v TABLE="${AMPLICON_TABLE}"
    '\''BEGIN {FS = " ";
            while ((getline < SWARM) > 0) {
                swarms[$1] = $0
            }
            FS = "\t";
            while ((getline < TABLE) > 0) {
                table[$1] = $0
            }
           }
     {
      seed = $3 "_" $4;
      n = split(swarms[seed], OTU, "[ _]");
      for (i = 1; i < n; i = i + 2) {
          s = split(table[OTU[i]], abundances, "\t");
          for (j = 1; j < s; j++) {
              samples[j] += abundances[j+1]
          }
      }
      printf "%s\t", $3;
      for (j = 1; j < s; j++) {
          printf "\t%s", samples[j]
      }
     printf "\n";
     delete samples
     }'\'' | awk '\''{print $NF,$0}'\'' | sort -nr | cut -f2- -d" " | awk '\''{gsub("\t\t","\t");print NR"\t"$0}'\'' >> ${OTU_TABLE}'
#execute command
eval $cmd
#write command to command trace file
echo -e "\n#Create OTU table:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

## TAXONOMIC ASSIGNMENT AND MERGE WITH OTU TABLE
# taxonomic assignment
# create command
cmd='
${MOTHUR} "#set.dir(output=.);classify.seqs(fasta="${S_SEEDS}".no.singletons,
reference="${REF_USED_FP}"/"${REF_DBS}"_"${FORWARDPRIMER}"_"${REVERSEPRIMER_RC}"_"${PRIMER_MISMATCH_REF}"_"${lenFP_CUT_REF}"_"${lenRP_CUT_REF}".fasta,
taxonomy="${REF_USED_FP}"/"${REF_DBS}"_"${FORWARDPRIMER}"_"${REVERSEPRIMER_RC}"_"${PRIMER_MISMATCH_REF}"_"${lenFP_CUT_REF}"_"${lenRP_CUT_REF}".tax,
processors="${THREADS}",cutoff="${RDP_CUTOFF}", probs=F);get.current();
rename.file(taxonomy=current,new=swarm."${REF_DBS}"_${FORWARDPRIMER}_${REVERSEPRIMER_RC}_${PRIMER_MISMATCH_REF}_${lenFP_CUT_REF}_${lenRP_CUT_REF}.wang.taxonomy,shorten=false);"'
# execute command
eval $cmd
#print command to command trace file
echo -e "\n#Taxonomic assignment of OTU representatives:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

mv "mothur"*"logfile" ${LOG_FP}

# merge taxonomy with OTU table
# create command
cmd='
echo -e "$(head -n 1 ${OTU_TABLE})\t${REF_DBS}_taxonomy" > ${PROJECT_ID}"_OTU_table.csv";
awk -F"[_ \t]" '\''FNR==NR{a[$2]=$0; next}($1 in a){printf a[$1]"\t"; for (i=3;i<NF;i++) {printf($i"_")}; printf($NF"\n") }'\'' ${OTU_TABLE} swarm."${REF_DBS}"_${FORWARDPRIMER}_${REVERSEPRIMER_RC}_${PRIMER_MISMATCH_REF}_${lenFP_CUT_REF}_${lenRP_CUT_REF}.wang.taxonomy 
>> ${PROJECT_ID}"_OTU_table.csv"'
# execute command
eval $cmd
#print command to command trace file
echo -e "\n#Adding taxonomy to OTU table:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


## MAKE BIOM FILE; FUTURE: also external other meta data

#make sample metadata map from raw file to sample map
# create command
cmd='cat ${PROJECT_ID}".map" | awk '\''BEGIN{print "#Sample\tR1_raw_file\tR2_raw_file"}{print}'\'' > ${PROJECT_ID}"_sample_metadata.map"'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#Make sample metadata map:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

# get pure OTU table + taxonomy
#create command
cmd='
cat ${PROJECT_ID}"_OTU_table.csv" |
awk '\''{$2="";$(NF-1)="";print}'\'' |
tr -s " " | tr " " "\t"
> ${PROJECT_ID}"_OTU_table_pure.csv"'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#Get pure OTU table + taxonomy:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

# get observation meta data "total abundance" + "representative amplicon id"
# create command
cmd='
cat ${PROJECT_ID}"_OTU_table.csv" |
awk '\''{print $1"\t"$(NF-1)"\t"$2}'\''
> ${PROJECT_ID}"_observation_metadata.map"'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#get observation meta data (total abundance + representative amplicon id):" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

# add representative amplicon sequence to observation metadata
#create command
cmd='
awk '\''BEGIN{print("#OTU\ttotal\tamplicon\tsequence")}
NR==FNR{a[$3]=$1;b[$3]=$2;c[$3]=$3;next}
(substr($0,2,index($0,"_")-2) in c)
{printf(a[substr($0,2,index($0,"_")-2)]"\t"b[substr($0,2,index($0,"_")-2)]"\t"c[substr($0,2,index($0,"_")-2)]"\t");
getline; printf $1"\n"}'\'' ${PROJECT_ID}"_observation_metadata.map" swarm.seeds.no.singletons
> ${PROJECT_ID}"_observation_metadata_full.map"'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#add representative amplicon sequence to observation metadata:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

# convert pure OTU table + taxonomy to biom format
#create command
cmd='
biom convert -i ${PROJECT_ID}"_OTU_table_pure.csv" -o ${PROJECT_ID}"_OTU_table_pure.biom" --to-hdf5'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#convert pure OTU table + taxonomy to biom format:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


# add external sample meta data to biom if provided
#create command
cmd='
if [ -s "${ADD_METADATA_FILE}" ]; 
then biom add-metadata -i ${PROJECT_ID}"_OTU_table_pure.biom" -o ${PROJECT_ID}"_OTU_table_add_hdf5.biom" --sample-metadata-fp ${ADD_METADATA_FILE}; fi;
if [ ! -s ${PROJECT_ID}"_OTU_table_add_hdf5.biom" ];
then mv ${PROJECT_ID}"_OTU_table_pure.biom" ${PROJECT_ID}"_OTU_table_add_hdf5.biom"; fi'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#add external sample meta data to biom if provided:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

#add sample meta data from raw file map and observation meta data
#create command
cmd='
biom add-metadata -i ${PROJECT_ID}"_OTU_table_add_hdf5.biom" -o ${PROJECT_ID}"_OTU_table_full_hdf5.biom" --observation-metadata-fp  ${PROJECT_ID}"_observation_metadata_full.map" --sample-metadata-fp ${PROJECT_ID}"_sample_metadata.map"'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#add sample meta data from raw file map and observation meta data:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


# convert from hdf5 to json format
#create command
cmd='
biom convert -i ${PROJECT_ID}"_OTU_table_full_hdf5.biom" -o ${PROJECT_ID}"_OTU_table_full_json.biom" --to-json'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#convert from hdf5 to json format:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


# clean from SAMPLE_PREFIX and SAMPLE_SUFFIX
# create command
cmd='
if [ "$SAMPLE_PREFIX_REGEXP" == "" ]; then SAMPLE_PREFIX_REGEXP=`date +%s%N | md5sum | cut -f1 -d " "`; fi;
if [ "$SAMPLE_SUFFIX_REGEXP" == "" ]; then SAMPLE_SUFFIX_REGEXP=`date +%s%N | md5sum | cut -f1 -d " "`; fi;
cat ${PROJECT_ID}"_OTU_table_full_json.biom" |
awk '\''{print gensub(/"id": "'${SAMPLE_PREFIX_REGEXP}'/,"\"id\": \"","g",$0)}'\'' |
awk '\''{print gensub(/'${SAMPLE_SUFFIX_REGEXP}'"/,"\"","g",$0)}'\''
> ${PROJECT_ID}"_OTU_table_full_json_clean.biom";
cat ${PROJECT_ID}_OTU_table_pure.csv |
awk '\''{print gensub(/\t'${SAMPLE_PREFIX_REGEXP}'/,"\t","g",$0)}'\'' |
awk '\''{print gensub(/'${SAMPLE_SUFFIX_REGEXP}'\t/,"\t","g",$0)}'\''
> ${PROJECT_ID}_OTU_table_pure_clean.csv'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#clean from SAMPLE_PREFIX and SAMPLE_SUFFIX:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt


# convert back from json to hdf5 format
#create command
cmd='
biom convert -i ${PROJECT_ID}"_OTU_table_full_json_clean.biom" -o ${PROJECT_ID}"_OTU_table_full_hdf5_clean.biom" --to-hdf5'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#convert back from json to hdf5 format:" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt

# add observation meta data to OTU table (csv)
#create command
cmd='
echo -e `head -1 ${PROJECT_ID}"_OTU_table_pure_clean.csv"`"\t"`head -1 ${PROJECT_ID}"_observation_metadata_full.map" | awk '\''{printf($2"\t"$3"\t"$4)}'\''`
> ${PROJECT_ID}_OTU_table_pure_clean.csv;
awk '\''FNR==NR{a[$1]=$2"\t"$3"\t"$4; next}($1 in a){printf $0"\t"a[$1]"\n"}'\''
${PROJECT_ID}_observation_metadata_full.map 
${PROJECT_ID}_OTU_table_pure_clean.csv 
>> ${PROJECT_ID}_OTU_table_pure_clean_full.csv
'
#execute command
eval $cmd
#print command to command trace file
echo -e "\n#add observation meta data to OTU table (csv):" >> q-zip_seq_of_coms.txt
echo $cmd >> q-zip_seq_of_coms.txt



# analysis done


## TIMES
END_DATE=`date +%Y-%m-%d`
END_TIME=`date +%H-%M-%S`
echo "Analysis starte at "${ST_DATE} ${ST_TIME}
echo "Analysis ended at "${END_DATE} ${END_TIME}
DUR_SEC_END=${SECONDS}
DURATION=`echo ${DUR_SEC_END}-${DUR_SEC_ST} | bc`
echo "$((${DURATION} / 3600)) hours, $(((${DURATION} / 60) % 60)) minutes and $((${DURATION} % 60)) seconds elapsed."


sleep 2

echo "Write html index file for navigation"

echo '
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html><head>
<style><endnote><head>
<style></style>
<meta content="text/html; charset=ISO-8859-1" http-equiv="content-type"><title>Qzip: a rapid meta-barcoding pipeline</title>
<style type="text/css">
body { font-family:Verdana,Arial,Helvetica,sans-serif; background-color: #ffffff; }
A:link {color: rgb(83, 141, 189) ; text-decoration:none</style></head>
<body>
<br><table style="width: 1031px; height: 101px;" border="0"><tbody><tr valign="top"><td><img style="width: 311px; height: 119px;" src="https://www.awi.de/typo3conf/ext/sms_boilerplate/Resources/Public/Images/AWI/awi_logo.svg"></td><td><div style="">
<big><br>Alfred Wegener Institute for Polar and Marine Research
<br><br></big><a href="https://www.awi.de/en/about-us/organisation/staff/stefan-neuhaus.html">StefanNeuhaus | Scientific Computing | Bioinformatics <br>
phone: +49(471)4831-2329 | e-mail: Stefan.Neuhaus@awi.de <br>
</a>
</div></td></tr></tbody></table><br>
<hr><br>
<br><h1 style="margin-left: 200px;"><small>QZIP - rapid metabarcoding pipeline</small></h1>
<br><br><table style="width: 891px; height: 215px;" border="0">
<tbody>
<tr><td style="width: 213px;"><h3>Project</h3></td><td style="width: 221px;"></td><td style="width: 203px;"><h3>Target</h3></td><td style="width: 224px;"></td></tr><tr><td style="width: 213px;">Project name</td><td style="width: 221px;">'${PROJECT_NAME}'</td><td style="width: 203px;">Molecule<br></td><td style="width: 224px;">'${TARGET_MOLECULE}'</td></tr><tr><td style="width: 213px;">
Name(s) of PI</td><td style="width: 221px;">
'${PI_NAMES}'<br></td><td style="width: 203px;">
Region</td><td style="width: 224px;">
'${TARGET_REGION}'</td></tr><tr><td style="width: 213px;">
Date of analysis</td><td style="width: 221px;">
'${ST_DATE}'</td><td style="width: 203px;">
Fwd primer</td><td style="width: 224px;">
'${FORWARDPRIMER}'</td></tr><tr><td style="width: 213px;">
Start time</td><td style="width: 221px;">
'${ST_TIME}'</td><td style="width: 203px;">
Rev primer</td><td style="width: 224px;">
'${REVERSEPRIMER}'</td></tr><tr><td style="width: 213px;">
End time</td><td style="width: 221px;">'${END_TIME}'</td><td style="width: 203px;">
Reference datasets</td><td style="width: 224px;">'${REF_DBS}'</td></tr><tr><td style="width: 213px;">
Sample meta data</td><td style="width: 221px;"><a href='${PROJECT_ID}.map'>'${PROJECT_ID}.map'</a></td><td style="width: 203px;"></td><td style="width: 224px;"></td></tr>
</tbody>
</table>
<span style="font-weight: bold;"></span><br>
<table style="width: 893px; height: 143px;" border="0">
<tbody>
<tr>
<td colspan="2" rowspan="1" style="width: 192px;"><h3>Input and settings</h3></td>
<td style="width: 183px;"><h3>Pipeline</h3></td><td style="width: 208px;"></td>
</tr>
<tr valign="top">
<td style="width: 192px;">Raw files sample map</td>
<td style="width: 200px;">&lt;Link to sample map&gt;&nbsp;</td><td style="width: 183px;">Description</td><td style="width: 208px;">&lt;Link to description vX&gt;</td>
</tr><tr><td>Reference sets used</td><td>&lt;Link to reference sets&gt;</td><td>Source code</td><td>&lt;link to script&gt;</td></tr><tr><td>Configuration file</td><td>&lt;Link to config file&gt;</td><td></td><td></td></tr>
</tbody>
</table><span style="font-weight: bold;"><br></span>
<table style="width: 893px; height: 227px;" border="0">
<tbody>
<tr>
<td style="width: 192px;"><h3>Results</h3></td>
<td style="width: 200px;"></td><td colspan="2" rowspan="1" style="width: 183px;"><h3>Sequence quality assessment</h3></td>
</tr>
<tr valign="top">
<td style="width: 192px;">OTU tables</td>
<td style="width: 200px;">&lt;Link to OTU tables&gt;&nbsp;</td><td style="width: 183px;">Workflow statistics</td><td style="width: 208px;">&lt;Link to stat file&gt;</td>
</tr><tr><td>OTU tables subsampled<br></td><td>&lt;Link to subsampled OTU tables&gt;</td><td>FastQC report<br></td><td>&lt;Link to condensed report&gt;</td></tr><tr><td>OTU sequences (seeds)</td><td>&lt;Link to seeds file&gt;<br></td><td></td><td></td></tr><tr><td>Amplicon to OTU mapping (swarm)</td><td>&lt;Link to swarm file&gt;<br></td><td></td><td></td></tr><tr><td>Taxonomic overview</td><td>&lt;Link to crona plot&gt;</td><td></td><td></td></tr>
</tbody>
</table><span style="font-weight: bold;"></span><br>
<table style="width: 450px; height: 97px;" border="0">
<tbody>
<tr>
<td colspan="2" rowspan="1" style="height: 48px; width: 200px;"><h3>Workflow documentation</h3></td>
</tr>
<tr valign="top">
<td style="height: 24px; width: 192px;">Log file</td>
<td style="height: 24px; width: 200px;">&lt;Link to log file&gt;&nbsp;</td>
</tr>
</tbody>
</table><span style="font-weight: bold;"></span><br><br><hr>Project ID: '${PROJECT_ID}'&nbsp;<br>
</body></html>
' > q-zip_index.html



ln -s ../seq_number_stats.txt ../q-zip_workflow.log ../q-zip_commands.sh ../q-zip_parameters.txt ../q-zip_index.html ../q-zip_seq_of_coms.txt \
../${AMPLICON_TABLE} ../${S_SEEDS}".no.singletons" ../${S_SWARM}".no.singletons" ../${LOG_FP} \
../${PROJECT_ID}"_sample_metadata.map" ../${ADD_METADATA_FILE} ../${PROJECT_ID}"_observation_metadata_full.map" \
../${PROJECT_ID}"_OTU_table_pure_clean.csv" ../${PROJECT_ID}"_OTU_table_full_json_clean.biom" ../${PROJECT_ID}"_OTU_table_full_hdf5_clean.biom" \
./${RESULTS_DIR}

ln -s ../${REF_USED_FP}"/"${REF_DBS}"_"${FORWARDPRIMER}"_"${REVERSEPRIMER_RC}"_"${PRIMER_MISMATCH_REF}"_"${lenFP_CUT_REF}"_"${lenRP_CUT_REF}".fasta" ./${RESULTS_DIR}
ln -s ../${REF_USED_FP}"/"${REF_DBS}"_"${FORWARDPRIMER}"_"${REVERSEPRIMER_RC}"_"${PRIMER_MISMATCH_REF}"_"${lenFP_CUT_REF}"_"${lenRP_CUT_REF}".tax" ./${RESULTS_DIR}
#ln -s ../../${REF_USED_FP}"/"${REF_DBS}"_"${FORWARDPRIMER}"_"${REVERSEPRIMER_RC}"_"${PRIMER_MISMATCH_REF}"_"${lenFP_CUT_REF}"_"${lenRP_CUT_REF}".fasta" ./${RESULTS_DIR}/${REF_USED_FP}
#ln -s ../../${REF_USED_FP}"/"${REF_DBS}"_"${FORWARDPRIMER}"_"${REVERSEPRIMER_RC}"_"${PRIMER_MISMATCH_REF}"_"${lenFP_CUT_REF}"_"${lenRP_CUT_REF}".tax" ./${RESULTS_DIR}/${REF_USED_FP}

zip -r -j -m ${RESULTS_DIR}"/"`basename ${REF_USED_FP}`".zip" ${RESULTS_DIR}"/"${REF_DBS}*
zip -r -m ${RESULTS_DIR}"/"${AMPLICON_TABLE}".zip" ${RESULTS_DIR}"/"${AMPLICON_TABLE}
zip -r -m ${RESULTS_DIR}"/"${S_SWARM}".no.singletons.zip" ${RESULTS_DIR}"/"${S_SWARM}".no.singletons"

zip -r ${RESULTS_DIR}"_"${PROJECT_ID}.zip ${RESULTS_DIR}


}

#WORKAROUND: DEFINE A FUNCTION FOR THE WORKFLOW TO HAVE EASY TERMINAL LOGGING
q-zip_run 2>&1 | tee q-zip_workflow.log

