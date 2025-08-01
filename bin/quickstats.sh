#!/usr/bin/env sh

# Get counts of A->G and C->T edits per contig/chromosome in a Reditools or Jacusa2 (extended BED format) output file
# This is a useful way to make quick checks of the more complex output of pluviometer.py
# Usage:
# sh quickstats.sh -i <INPUT FILE> -f <FORMAT>
# FORMAT can be either "jacusa2" or "reditools" (works for Reditools2 and Reditools3)

input=""
format=""

while getopts "i:f:" opt; do
    case $opt in
        i)
            input="$OPTARG"
            ;;
        f)
            format="$OPTARG"
            ;;
    esac
done

if [ "$format" == "reditools" ]; then
    grep -P ".+\t\d+\tA" $input | cut -f1,7 | tr -d "[]" | sed "s/, /\t/g" | awk '{sum[$1] += $4} END {for (key in sum) print "A->G", key, sum[key]}'
    grep -P ".+\t\d+\tC" $input | cut -f1,7 | tr -d "[]" | sed "s/, /\t/g" | awk '{sum[$1] += $5} END {for (key in sum) print "C->T", key, sum[key]}'
elif [ "$format" == "jacusa2" ]; then
    grep -P "\tA$" $input | cut -f1,7 | sed "s/,/\t/g" | awk '{sum[$1] += $4} END {for (key in sum) print "A->G", key, sum[key]}'
    grep -P "\tC$" $input | cut -f1,7 | sed "s/,/\t/g" | awk '{sum[$1] += $5} END {for (key in sum) print "C->T", key, sum[key]}'
fi
