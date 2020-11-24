#!/bin/bash

HLA=$1
file=$2
source $3
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
    echo  -e "$NETMHCPAN -f $i -p -a $HLA\tnetMHCpan-EL_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "$NETMHCSTABPAN -f $file -p -a $HLA\tnetMHCstab_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "$NETMHC -f $file -p -a $net_HLA\tnetMHC_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "$PICKPOCKET -f $file -p -a $HLA\tpickpocket_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "$MIXMHCPRED -i $file -a $HLA -o MixMHCpred_$len_suffix\ttmp_MIX" >> "$HLA"_run_list.txt
    
done


cp ../$2 ./EnsembleMHC_pep_pred.csv
