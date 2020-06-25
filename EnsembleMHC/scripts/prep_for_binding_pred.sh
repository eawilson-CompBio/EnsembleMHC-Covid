#!/bin/bash

HLA=$1
file=$2
MHC=$3
#HLA_link=$(realpath $HLA_list)
#HLA=${HLA_list%%_*}
len=(8 9 10 11 12 13 14)

mkdir "$HLA"_predictions

for i in "${len[@]}"
do
   
    cat $file | awk -F ',' 'NR>1{print $1,$3}' | uniq |tr -d '"' | grep $i | grep -v "X" | awk '{print $1}' > "$HLA"_predictions/"$HLA"_"$i"mer.txt
done

cd "$HLA"_predictions



for i in $(ls)
do
    len_suffix=${i##*_}
    file=$(realpath $i)
    net_HLA=${HLA:0:7}${HLA:8}
    echo  -e "$MHC/netMHCpan-4.0/netMHCpan -f $i -p -a $HLA\tnetMHCpan-EL_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "$MHC/netMHCstabpan-1.0/netMHCstabpan -f $file -p -a $HLA\tnetMHCstab_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "$MHC/netMHC-4.0/netMHC -f $file -p -a $net_HLA\tnetMHC_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "$MHC/pickpocket-1.1/PickPocket -f $file -p -a $HLA\tpickpocket_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "$MHC/MixMHCpred-2.0.2/MixMHCpred -i $file -a $HLA -o MixMHCpred_$len_suffix\ttmp_MIX" >> "$HLA"_run_list.txt
    
done


cp ../$2 ./EnsembleMHC_pep_pred.csv
