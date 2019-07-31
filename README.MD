# Qzip - a rapid meta-barcoding pipeline

Qzip is a bash script which wraps and connects common bioinformatics tools to produce sample by observation contingency tables from demultiplexed Illumina paired-end amplicon sequences.
The usage of GNU parallel for most parallelizable processing steps and the usage of Swarm as fast clustering method makes it extremely speedy. On a 20 core linux system (physical cores) the processing of about 300 million sequences from a large scale project only takes ~ 6 h and the processing of 10 million sequences from a handy-sized project only takes a few minutes.

## features
- easy to install and use
- based on open source tools
- rigorous sequence pre-filtering
- Swarm-based clustering
- reference-independence until annotation
- contingency tables as biom file and as csv tables
- sequence number and length statistics for each step
- command-line logging and command recording

## outline
1. filter and trim of the reference sequences (cutadapt)\
Only sequences which contain both primers are passing the filter. Flanking primers and segments beyond are removed to obtain a homolog region block.
2. 3'-quality trimming of paired end query reads (trimmomatic)
3. Overlap-merge of paired-ends (vsearch)
4. Adjustment of orientation 5'-> 3' and trimming of primer segements (cutadapt)
5. Feature filter (vsearch)\
A Min and max length filter and a "max number of expected errors" filter are implemented.
6. Chimera filter (vsearch)\
Prediction of putative chimeric sequences sample-wise and removal of them from the sample if they do not occur in any other sample.
7. Swarm clustering (swarm)
8. Annotation by naive bayesian classifier (mothur)

## prerequisites:
Please have available the following tools. The brackets indicate the used and tested versions. Higher versions will probably work as well, lower versions might not.
- biom (2.1.6)
- trimmomatic (0.38)
- cutadapt (1.17)
- swarm (2.1.8)
- vsearch (2.8.2)
- mothur (1.40.5)
- GNU parallel (20141022)

## installing:
copy q-zip_commands.sh and q-zip_parameters.txt into your working directory:\
wget https://raw.githubusercontent.com/PyoneerO/qzip/master/q-zip_parameters.txt \
wget https://raw.githubusercontent.com/PyoneerO/qzip/master/q-zip_commands.sh

## input:
- demultiplexed Miseq paired-end raw data (fastq/fastq.gz) with 
  1. sufficient overlap of paired-end reads (80-100 bp overlap of raw reads recommended as \
    a. read trimming reduces length of overlap for paired-end merging (recommended min 50 bp) and \
    b. a long overlap improves amplicon quality)
  2. and retained PCR primers
- sample metadata file (optional)

## naming scheme for raw sequence data files:
SampleID_WellID/AlternativeSampleID_LaneNo_ReadNo_001.fastq.gz\
exemplary raw file name:\
95-HE492_S47_L001_R1_001.fastq.gz

or

MachineID_RunIndex_FlowcellID_SampleID_WellID/AlternativeSampleID_LaneNo_ReadNo_001.fastq.gz\
(prefix "MachineID_RunIndex_FlowcellID_" added to improve traceability)\
exemplary raw file name:\
M03457_0023_000000000-C94HP_95-HE492_S47_L001_R1_001.fastq.gz

If you provide regular expressions for Prefix (MachineID_RunIndex_FlowcellID_) and Suffix (_WellID/AlternativeSampleID) in the parameter file, the sample names in the OTU tables will be cleaned from them.

## parameter settings:
open the file q-zip_parameters.txt with an editor and set parameters especially for:
- project name
- primer sequences
- reference set (formatted according to mothur)
- paths to the needed executables

## reference sets
the pipeline accepts reference sets which consist of two files:
1. a fasta file containing the reference sequencing having unique IDs
2. a text file containing ID to taxonomy mapping. The taxonomic items need to be separated by ";" and the taxonomic tree needs to be "artifically balanced", which means that the taxonomy strings in the map are required to have equal taxonomic level depth. 

PR2 datasets and datasets distributed by the mothur webiste are already in the correct format.
Please use the same name for taxonomy map and fasta file, except of the filename suffixes which need to be "tax" for the mapping file and "fasta" for the fasta file. Please avoid any special character in the name strings, except of underscore "\_". You then need to specify the used reference set in the parameter file without filename suffix.

## workflow execution:
simply source the script:
. ./q-zip_commands.sh

## output:
a zip file ProjectName_Date_Time.zip containing the following:
- q-zip_commands.sh - used workflow file
- q-zip_parameters.txt - used parameter file
- q-zip_seq_of_coms.txt - sequence of used commands, sourcing of this file should reproduce the results
(except of the annotation as the mothur naive bayesian classifier has a random moment)
- q-zip_workflow.log - log file
- refDB.zip - primer filtered and trimmed reference set used for annotation
- seq_number_stats.txt - sequence number and length statistics for the preprocessing steps
- amplicon_table.csv - distribution table of amplicons per sample
- swarm.seeds.no.singletons - swarm/OTU representing sequences (singleton removal applied by default!)
- swarm.swarm.no.singletons - list of amplicons per swarm/OTU (singleton removal applied by default!)
- ProjectID_observation_metadata_full.map - OTU/observation metadata
- ProjectID_sample_metadata.map - sample metadata
- ProjectID_OTU_table_full_json_clean.biom - OTU table as biom, json format, all provided metadata included
- ProjectID_OTU_table_full_hdf5_clean.biom - OTU table as biom, hdf5 format, all provided metadata included
- ProjectID_OTU_table_pure_clean.csv - OTU table + taxonomie as csv
- ProjectID_OTU_table_pure_clean_full.csv - OTU table + taxonomie + sequence as csv

