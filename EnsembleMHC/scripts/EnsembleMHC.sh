#!/bin/bash

export PATH_FILE=/Users/eawilso6/Covid-19/EnsembleMHC-Covid19/EnsembleMHC/ALGORITHM_PATHS.sh

source $PATH_FILE



protein=''
HLA=''
CORES=''

print_usage() {
    printf '%s\n' "Usage: bash run_EnsembleMHC.sh -p protein -a HLA -t "
}

while getopts 'p:ha:c' flag; do
  case "${flag}" in
    p) protein="${OPTARG}" ;;
    a) HLA="${OPTARG}" ;;
    c) CORES="${OPTARG}" ;;
    h) print_usage
       exit 1 ;;
  esac
done

if [[ -z ${HLA+x} ]]
then
    echo "specify HLA"
    exit 0
fi

if grep -Fxq "$HLA" "$ENSEMBLEMHC_PATH"/scripts/HLA_list.txt
then
    echo "HLA: "$HLA
else
    echo -e  "HLA: $HLA is not supported by EnsembleMHC. Please select one of the following:"
    cat $ENSEMBLEMHC_PATH/scripts/HLA_list.txt
    exit 0
fi



if [[ -z ${protein+x} ]]
then
    echo "specify protein"
    exit 0
fi


if [[ -z ${CORES+x} ]]
then
    echo "number of threads not set. Defaulting to all available cores"
    CORES=$(nproc)
fi


#split the protein file
Rscript $ENSEMBLEMHC_PATH/scripts/split_prot.R $protein 


##prep for binding predictions and generate run file
##the run file coordinates the parallelization of the predictions 
$ENSEMBLEMHC_PATH/scripts/prep_for_binding_pred.sh $HLA EnsembleMHC_pep_pred.tmp $PATH_FILE

cd $HLA"_predictions"
##excute prediction for the nets and MixMHC
$ENSEMBLEMHC_PATH/scripts/par_run_predict.sh "$HLA"_run_list.txt $CORES

#add HLA column
Rscript $ENSEMBLEMHC_PATH/scripts/add_HLA_column.R *_pred.csv $HLA
conda activate mhcflurry-env
##excute MHCflurry
$ENSEMBLEMHC_PATH/scripts/run_mhcflurry.sh *_pred.csv

#clean up dir
$ENSEMBLEMHC_PATH/scripts/processing_and_clean_up.sh

ls -l *out | awk '$5==0{print $9}' > emptyfiles.txt
if [[ $(cat emptyfiles.txt |wc -l) -gt 0 ]]
then
    rm $(cat emptyfiles.txt)
fi


#process files and generate whats need
Rscript $ENSEMBLEMHC_PATH/scripts/combine_files.R $HLA

cp *scored_peptides.csv ..
cd ../
mkdir tmp_out
mv $HLA"_predictions"/ tmp_out/
#mv all_*.out "$HLA"_all_raw






