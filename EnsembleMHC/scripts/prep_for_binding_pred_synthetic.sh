#!/bin/bash

HLA_list=$1
HLA_link=$(readlink -f $HLA_list)
HLA=${HLA_list%%_*}
len=(8 9 10 11 12 13 14)

mkdir $HLA

for i in "${len[@]}"
do
   
    cat $HLA_list | awk -F ',' 'NR>1{print $1,$3}' | uniq |tr -d '"' | grep -v "X" | awk '{print $1}' > "$HLA"/"$HLA"_"$i"mer.txt
done

cd $HLA



for i in $(ls)
do
    len_suffix=${i##*_}
    file=$(readlink -f $i)
    net_HLA=${HLA:0:7}${HLA:8}
    echo  -e "/home/eawilso6/apps/nCov_design/MHC-I_software/netMHCpan-4.0/netMHCpan -f $i -p -a $HLA\tnetMHCpan-EL_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "/home/eawilso6/apps/nCov_design/MHC-I_software/netMHCpan-4.0/netMHCpan -BA -f $file -p -a $HLA\tnetMHC-BA_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "/home/eawilso6/apps/nCov_design/MHC-I_software/netMHCstabpan-1.0/netMHCstabpan -f $file -p -a $HLA\tnetMHCstab_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "/home/eawilso6/apps/nCov_design/MHC-I_software/netMHC-4.0/netMHC -f $file -p -a $net_HLA\tnetMHC_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "/home/eawilso6/apps/nCov_design/MHC-I_software/pickpocket-1.1/PickPocket -f $file -p -a $HLA\tpickpocket_$len_suffix" >> "$HLA"_run_list.txt
    echo  -e "/home/eawilso6/apps/nCov_design/MHC-I_software/MixMHCpred/MixMHCpred -i $file -a $HLA -o MixMHCpred_$len_suffix\ttmp_MIX" >> "$HLA"_run_list.txt
    
done

cp $HLA_link .

