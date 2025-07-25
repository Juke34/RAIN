#!/usr/bin/env sh

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
