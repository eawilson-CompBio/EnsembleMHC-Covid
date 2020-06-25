#!/bin/bash



    cat MixMHCpred* | grep -E "^[A-Z]" | awk 'NR>1{print $1,$6}' > all_MixMHCpred.out

    cat netMHC_* | grep PEPLIST | grep -v Protein  | awk '{print $3,$14}' > all_netMHC.out
    cat netMHCpan-EL_* | grep -E "^\s+[0-9]+\s+HLA" | awk '{print $3,$13}' > all_netMHCpan-EL.out
    cat netMHCstab_* | grep -E "^\s+[0-9]+\s+HLA" | awk '{print $3,$7}' > all_netstab.out
    cat pickpocket_* | grep -E "^\s+[0-9]+\s+HLA" | awk '{print $3,$5}' > all_pickpocket.out

    if [ "$(uname)" == "Darwin" ]
    then
	gsed -i '1s/^/peptide MixMHCpred\n/' all_MixMHCpred.out 
	gsed -i '1s/^/peptide netMHC_affinity\n/' all_netMHC.out
	gsed -i '1s/^/peptide netMHCpan_EL_affinity\n/' all_netMHCpan-EL.out
	gsed -i '1s/^/peptide netstab_affinity\n/' all_netstab.out
	gsed -i '1s/^/peptide pickpocket_affinity\n/' all_pickpocket.out
    else
	sed -i '1s/^/peptide MixMHCpred\n/' all_MixMHCpred.out 
	sed -i '1s/^/peptide netMHC_affinity\n/' all_netMHC.out
	sed -i '1s/^/peptide netMHCpan_EL_affinity\n/' all_netMHCpan-EL.out
	sed -i '1s/^/peptide netstab_affinity\n/' all_netstab.out
	sed -i '1s/^/peptide pickpocket_affinity\n/' all_pickpocket.out
    fi
    


    
#mkdir del_files
#mv $(ls | grep -v "all\|csv\|pred.out\|del") del_files




