#!/bin/bash
#$ -cwd
#$ -notify
#$ -P gentech-rqc.p
#$ -N AC-417
#$ -pe pe_slots 2
#$ -l h_rt=12:00:00
#$ -l ram.c=5G
#$ -t 1-9997
#$ -o logs
#$ -e logs

if [ ! -e logs ]
then
mkdir logs
fi

./run.sh
