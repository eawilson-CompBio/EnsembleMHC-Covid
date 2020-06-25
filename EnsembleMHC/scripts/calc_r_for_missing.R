a<-commandArgs(trailingOnly = T)
library(data.table)
library(dplyr)

tab<-lapply(list.files(pattern="out"),function(w){
    fread(w,header=T)%>%unique()
})

combo<-Reduce(merge,tab)
combo<-combo%>%unique()
ref<-fread(a)


mhcflurry<-fread(list.files(pattern="mhcflurry.pred"))
mhcflurry_pres<-mhcflurry%>%select(peptide,mhcflurry_presentation_score)
mhcflurry_aff<-mhcflurry%>%select(peptide,mhcflurry_affinity_percentile)
combo<-combo%>%merge(mhcflurry_pres)%>%merge(mhcflurry_aff)%>%merge(ref)
HLA=ref$HLA[1]
combo$HLA<-HLA
write.csv(combo,file=paste(HLA,"ready_for_rf.csv",sep="_"),row.names=F)


norm<-function(x){(max(x)-x)/(max(x)-min(x))}

sfp<-data.frame(combo)[,c(2:(1+length(tab)+2),ncol(combo)-1)]%>%unique()
sfp$pickpocket_affinity<-norm(sfp$pickpocket_affinity)
sfp$mhcflurry_presentation_score<-norm(sfp$mhcflurry_presentation_score)
if(sum(colnames(sfp)=="MixMHCpred")>0){
sfp$MixMHCpred<-as.numeric(sfp$MixMHCpred)
}

PPV_calc<-function(x,y){
    tmp<-data.frame(as.numeric(x),y,stringsAsFactors = F)
    tmp<-tmp[order(tmp[,1]),]
    tar_ind<-which(tmp[,2]=="target")
    ind<-round(length(tar_ind)*.4)
    tab_PPV<-table(tmp$y[1:ind])
    c(PPV=tab_PPV["target"]/sum(tab_PPV),value=tmp[ind,1])
    
}


output1000<-do.call(rbind,lapply(1:10,function(w){
    pos<-sfp%>%slice(which(ident=="target"))
    neg<-sfp%>%slice(which(ident=="decoy"))

    test<-rbind(pos[sample(1:nrow(pos),round(.1*nrow(pos)),replace = F),],
          neg)
    

    
    out<-t(apply(test%>%select(-ident),2,function(x){
    PPV_calc(x,test$ident)
    }))%>%data.frame()
    out$algo<-row.names(out)
    out$iter<-w
    colnames(out)[1:2]<-c("PPV","score_thres")
    out
    

}))

output100<-do.call(rbind,lapply(1:10,function(w){
    pos<-sfp%>%slice(which(ident=="target"))
    neg<-sfp%>%slice(which(ident=="decoy"))

    test<-rbind(pos[sample(1:nrow(pos),round(.5*nrow(pos)),replace = F),],
          neg[sample(1:nrow(neg),round(.5*nrow(neg)),replace = F),])
    

    out<-t(apply(test%>%select(-ident),2,function(x){
    PPV_calc(x,test$ident)
    }))%>%data.frame()
    out$algo<-row.names(out)
    out$iter<-w
    colnames(out)[1:2]<-c("PPV","score_thres")
    out
    

}))

write.csv(output1000,file=paste(HLA,"1-1000_40_PPV.csv",sep="_"),row.names=F)
write.csv(output100,file=paste(HLA,"1-100_40_PPV.csv",sep="_"),row.names=F)
