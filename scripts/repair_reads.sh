#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 40G
#SBATCH -c 1
#SBATCH --time 5:00:00



REPAIR="/home/sersan/miniconda3/opt/bbmap-38.22-0/repair.sh"
FORWARD=$1
REVERSE=$2

bash $REPAIR in=$FORWARD in2=$REVERSE  out=$FORWARD\_repaired out2=$REVERSE\_repaired repair
