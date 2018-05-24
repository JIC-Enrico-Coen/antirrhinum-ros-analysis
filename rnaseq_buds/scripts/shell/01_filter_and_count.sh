#!/bin/bash
#SBATCH --job_name=fastqc
#SBATCH -o $HOME/jic/rnaseq_buds/scripts/shell/01_filter_and_count.stdout
#SBATCH -e $HOME/jic/rnaseq_buds/scripts/shell/01_filter_and_count.stderr
#SBATCH -n 48


#######################
## Setting things up ##
#######################

# Working directory
cd $HOME/jic/rnaseq_buds/data/reads

# Make directories for output
mkdir -p qc/raw   # output for QC reports of raw reads
mkdir filtered    # output for filtered reads
mkdir qc/filtered # output for QC reports of filtered reads
mkdir ../counts   # output for transcript counts

# Read CSV file and extract sample names (skip header line of file)
SAMPLES=$(cat ../sample_info.csv | cut -d "," -f 1 | tail -n +2)


##################
## QC raw reads ##
##################

# Concatenate all fastq files
for SAMPLE in $SAMPLES
do
  # Find the sample files and concatenate them
  find ./raw -type f -name "*$SAMPLE*" -exec cat {} \; > ./raw/$SAMPLE.fq.gz
done

# Run fastqc on the set of raw data
fastqc -t $SLURM_NTASKS -o ./qc/raw ./raw/*.fq.gz

# Compile fastqc reports with multiqc
multiqc --outdir ./qc/ --filename fastqc_reports_raw.html ./qc/raw/


####################
## Quality filter ##
####################

# Quality-filter reads with cutadapt
for SAMPLE in $SAMPLES
do
  cutadapt \
    -j $SLURM_NTASKS \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    --trim-n --quality-cutoff 20 --minimum-length 30 --max-n 1 \
    -o ./filtered/$SAMPLE.fq.gz \
    ./raw/$SAMPLE.fq.gz > \
    ./filtered/$SAMPLE.cutadapt.log
done

# Run fastqc on the set of filtered data
fastqc -t $SLURM_NTASKS -o ./qc/filtered ./filtered/*.fq.gz

# Compile fastqc reports with multiqc
multiqc --outdir ./qc/ --filename fastqc_reports_filtered.html ./qc/filtered/


# Compile filtering stats from cutadapt
for LOG in $(ls ./filtered/*cutadapt.log)
do
  # get the sample name
  NAME=$(basename $LOG | sed 's/\..*//')

  # Pipe to parse output to csv
  grep -A 4 "Total reads processed" $LOG | \
  sed 's/[:,]//g' | \
  sed 's/^  //g' | \
  sed 's/([^)]*)//g' | \
  sed 's/ \+ /,/g' | \
  sed 's/ $//g' | \
  sed "s/^/$NAME,/g" >> ./filtered/cutadapt_read_stats.csv

  grep "Total basepairs processed" $LOG | \
  sed 's/,//g' | \
  sed 's/ bp//g' | \
  sed 's/: \+/,/g' | \
  sed "s/^/$NAME,/" >> ./filtered/cutadapt_basepair_stats.csv

  grep "Total written" $LOG | \
  sed 's/,//g' | \
  sed 's/ bp//g' | \
  sed 's/([^)]*)//g' | \
  sed 's/ : \+/,/g' | \
  sed "s/^/$NAME,/" >> ./filtered/cutadapt_basepair_stats.csv
done



##########################
## Quantify transcripts ##
##########################

# Count reads aligned to each transcript
for SAMPLE in $SAMPLES
do
  # quantify transcripts with salmon
  salmon quant \
    -i $HOME/reference/antirrhinum/annotation/salmon_index/fmd/ \
    -l  U \
    -r ./filtered/${SAMPLE}.fq.gz \
    --output ../counts/${SAMPLE} \
    --seqBias --posBias \
    --threads ${SLURM_NTASKS}
done
