#!/bin/bash

#set the interal field seperato

run_file="$1"
CORES=$2

#function to make gnu parallel run the commands that i want
run_pred () {
    IFS=$'\t'
    a=($1)
    eval ${a[0]} >> ${a[1]}
}
export -f run_pred
cat $run_file | parallel -j $CORES run_pred {}

