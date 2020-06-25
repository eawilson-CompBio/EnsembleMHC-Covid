library(dplyr)
library(data.table)
library(stringr)
library(parallel)

##path to benmarcking data
rf_files<-list.files(path="XX",full.names = T)
##load the 52 select alleles
load("alleles_with_all_algos.R")
#read rf files
rf_files<-rf_files[which(str_extract(rf_files,"[A-C][0-9]+\\:[0-9]+")%in%sel_alleles)]
#set algorithm names
algos<-c("mhcflurry_affinity_percentile","mhcflurry_presentation_score","MixMHCpred","netMHC_affinity","netMHCpan_EL_affinity","netstab_affinity","pickpocket_affinity")

##PPV calculation function
## this will be later converted to FDr throught the following relation
## FDR = 1 - PPV
## x = the score for a given peptide, y = the idnetify of that peptide. thres = target repetorie sze 
PPV_calc<-function(x,y,thres){
    ##make data frame 
    tmp<-data.frame(as.numeric(x),y,stringsAsFactors = F)
    ##arrange in ascending order
    tmp<-tmp[order(tmp[,1]),]
    #get indicies of target peptides
    tar_ind<-which(tmp[,2]=="target")
    ##set the ind variable to be the number of target peptides x threshold
    ind<-round(length(tar_ind)*thres)
    ##calculate frequency table of peptide identities up to the ind-th peptide
    tab_PPV<-table(tmp$y[1:ind])
    ##calculate PPV and score of the ind-th peptide
    c(PPV=tab_PPV["target"]/sum(tab_PPV),value=tmp[ind,1])
    
}

## because the scores for pickpocket and MHCflurry presentation stronger binders associtated with high scores
## These scores are rescaled so that the best scoreing peptide will be 0 and the worst scoring peptide will be one.
##this is done so the inequality sign does not need to change for these algorithms
min_norm<-function(x){(max(x)-x)/(max(x)-min(x))}

##set the benchmarking threshold 
recall_threshold<-.5

##calculate bootstrap PPV and score threshold calculations
##this is written to be parallelized over 52 cores. so one core dedicated to each allele.
PPV_values<-mclapply(rf_files,function(file){
    ## read in file
    mat<-fread(file,stringsAsFactors = F)
    ## make sure there are no duplicate vakues
    mat<-mat%>%unique()
    ## min-norm the score scores for pickpocket and MHCflurry presentation
    mat$mhcflurry_presentation_score<-min_norm(mat$mhcflurry_presentation_score)
    mat$pickpocket_affinity<-min_norm(mat$pickpocket_affinity)
    ## perform 1000 boostrap iterations per allele
    bootstrap_sample<-do.call(rbind,lapply(1:1000,function(i){
        ## split the data set into the target data set and a decoy data set
        ## the decoy dataset contains a 100 fold excess of peptides
        tar<-mat%>%slice(which(ident=="target"))
        dec<-mat%>%slice(which(ident=="decoy"))
        ## ggenerate testing set by sampleing half od the target peptides and half of the decoy peptides
        set<-rbind(tar[sample(1:nrow(tar),round(nrow(tar)*.5),replace = F),],dec[sample(1:nrow(dec),nrow(dec)*.5),])

        ## calcualte the PPV for each algorithm     
        PPV<-t(sapply(colnames(set)[which(colnames(set)%in%algos)],function(p){
            PPV_calc(set[,p],set$ident,thres = recall_threshold)
        }))
        ## return a dataframe with the resulting PPV and the iteration number
        data.frame(PPV,iter=i)
    }))
},mc.cores = 52)

## name the resulting data frames
names(PPV_values)<-str_extract(rf_files,"[A-C][0-9]{2}\\:[0-9]{2}")

##find the median value for PPV and score threshold for each algorithm and allele combination
##save the results into a matrix
P_sum<-do.call(rbind,lapply(1:length(PPV_values),function(i){
  tmp<-PPV_values[[i]]%>%mutate(algo=str_remove(row.names(.),"[0-9]+"))%>%group_by(algo)%>%summarise(PPV=median(PPV.target,na.rm = T),value=median(value,na.rm = T))%>%mutate(HLA=names(PPV_values)[i])%>%slice(which(algo%in%algos))
}))

##save the results
save(P_sum,file="P_sum_median_1000_boot.R")

