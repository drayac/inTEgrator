#!/usr/bin/env bash
### Get the fasta from a bed file

# source shflags
. ~/scripts/bash/common/shflags
#define args (DEFINE_class varName defaultValue explanation4help shortVarname)
FLAGS_HELP="
Convert bed to fasta using a fasta of the genome. Can work recursively.

    Usage: bed2fasta [flags] bedfile.bed genome.fa

    args:  bedfile.bed: the bed file you want to convert
           genome.fa: fasta file of the genome

  "

DEFINE_boolean 'recursive' false 'If you would like to have a recursive behaviour in folder' 'r'
DEFINE_string 'outputdir' '.' 'Specify the outputdir' 'o'
DEFINE_boolean 'useName' false 'Use the bed names as header' 'n'

# parse the command-line
FLAGS "$@" || exit $?
eval set -- "${FLAGS_ARGV}"

outputdir=${FLAGS_outputdir}
bed=$1
genome=$2

if [ ! -f "$1" ] && [ ! -d "$1" ]; then
  echo "No argumens ? or file $1" && exit 1
fi

[ ! -f "$2" ] && echo "You need the genome.fasta file to do this..." && exit 1

set -e

mkdir -p $outputdir

if [ ${FLAGS_useName} -eq ${FLAGS_TRUE} ]; then
  flags="-name"
else
  flags="-fullHeader"
fi

if [ ${FLAGS_recursive} -eq ${FLAGS_TRUE} ]; then
  for f in $1/*.bed; do
    if [ -f $f ]; then
      fname=$(basename $f) 
      echo "$fname"
      bedtools getfasta $flags -s -fi $genome -bed $f -fo $outputdir/${fname%bed}fa
    else
      echo "Skipped $f"
    fi
  done
else 
  fname=$(basename $bed) 
  bedtools getfasta $flags -s -fi $genome -bed $bed -fo $outputdir/${fname%bed}fa
fi
