#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=4G

# Runtime
#$ -l h_rt=00:45:00

#$ -j y
#$ -cwd
##################
### Run script ###
##################

# $1 file 

DIR="$( dirname "$1" )"
base="$( basename "$1" .nodup.bam )"

source $HOME/.bashrc

cd $DIR
echo $DIR
echo $base

## Check if index is present. If not, create it:
if [[ ! -e ${base}.nodup.bam.bai ]];
  then
  echo '[INFO]: File does not seem to be indexed. Indexing now:'
  samtools index ${base}.nodup.bam
  fi

## Calculate %mtDNA:
mtReads=$(samtools idxstats "${base}.nodup.bam" | grep '^chrM\|^MT\|^M' | cut -f 3)
totalReads=$(samtools idxstats "${base}.nodup.bam" | awk '{SUM += $3} END {print SUM}')

echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'

## Remove contigs, mitochondria
samtools view -bh "${base}.nodup.bam" $(cat $HOME/scripts/hg19_chr.tsv) > "${base}.chr1-22xy.bam" 
samtools index "${base}.chr1-22xy.bam"
samtools flagstat "${base}.chr1-22xy.bam" > "${base}.chr1-22xy.bam.flagstat.qc"
echo "Done."
 
