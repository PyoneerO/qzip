# Qzip

Qzip is a bash script which wraps and connects common bioinformatics tools to produce OTU contingency tables from demultiplexed Illumina amplicon data.

## prerequisites (tested versions):
- biom (2.1.6)
- trimmomatic (0.38)
- cutadapt (1.17)
- swarm (2.1.8)
- vsearch (2.8.2)
- mothur (1.40.5)
- GNU parallel (20141022)

## installing:
git clone https://github.com/PyoneerO/qzip
copy q-zip_commands.sh and q-zip_parameters.txt into your working directory.

## input:
- demultiplexed Miseq paired-end raw data (fastq/fastq.gz)
- sample metadata file (optional)

## naming scheme for raw sequence data files:
<SampleID>_<WellID/AlternativeSampleID>_<LaneNo>_<ReadNo>_001.fastq.gz
exemplary raw file name:
95-HE492_S47_L001_R1_001.fastq.gz

or

<MachineID>_<RunIndex>_<FlowcellID>_<SampleID>_<WellID\/AlternativeSampleID>_<LaneNo>_<ReadNo>_001.fastq.gz
(prefix "<MachineID>_<RunIndex>_<FlowcellID>_" added for traceability reasons)
exemplary raw file name:
M03457_0023_000000000-C94HP_95-HE492_S47_L001_R1_001.fastq.gz

If you provide regular expressions for Prefix (<MachineID>_<RunIndex>_<FlowcellID>_) and Suffix 
(_<well id/alternative sample id>) in the parameter file, the sample names in the OTU tables will be cleaned from them.

## parameter settings:
open the file q-zip_parameters.txt with an editor and set parameters especially for:
- project name
- primer sequences
- reference set (formatted according to mothur)
- paths to the needed executables

## workflow execution:
simply source the script by . ./q-zip_commands.sh

## output:
a zip file <ProjectName>_<Date>_<Time>.zip contsining the following:
- q-zip_commands.sh - used workflow file
- q-zip_parameters.txt - used parameter file
- q-zip_seq_of_coms.txt - sequence of used commands, sourcing of this file should result reproduce the results
(except of the annotation as the mothur naive bayesian classifier has a random moment)
- q-zip_workflow.log - log file
- refDB.zip - primer filtered and trimmed reference set used for annotation
- seq_number_stats.txt - sequence number and length statistics for the preprocessing steps
- amplicon_table.csv - distribution table of amplicons per sample
- swarm.seeds.no.singletons - swarm/OTU representing sequences (singleton removal applied by default!)
- swarm.swarm.no.singletons - list of amplicons per swarm/OTU (singleton removal applied by default!)
- <project id>_observation_metadata_full.map - OTU/observation metadata
- <project id>_sample_metadata.map - sample metadata
- <project id>_OTU_table_full_json_clean.biom - OTU table as biom, json format, all provided metadata included
- <project id>_OTU_table_full_hdf5_clean.biom - OTU table as biom, hdf5 format, all provided metadata included
- <project id>_9_OTU_table_pure_clean.csv - OTU table + taxonomie as csv
