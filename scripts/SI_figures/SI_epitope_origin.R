data_path="~/Covid-19/EnsembleMHC-Covid19/datasets/"
Ensemble_PATH="~/Covid-19/EnsembleMHC-Covid19"

library(Biostrings)
library(ggthemes)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)

#read in the reference data
fasta<-readAAStringSet(paste0(data_path,"MN908947.3.fasta"))

#read in protein summary info
protein_summary<-read.csv(paste0(data_path,"MN908947.3_Summary.csv"),stringsAsFactors = F)
#name the proteins in the fasta file
names(fasta)<-protein_summary$gene
#make data frame with protein lengths and names
proteome_summary<-data.frame(protein_names=names(fasta),length=nchar(fasta))

#read in all SARS-CoV-2 sequences
SARS<-readAAStringSet(paste0(data_path,"SARS-CoV-2_strucutral_proteins.fasta"))

#select proteins 
SARS_sel<-SARS[which(nchar(SARS)%in%proteome_summary$length[str_which(proteome_summary$protein_names,"E|N|M|S")])]

#get the sequences attributed to each protein
#this is done by matching proteins by length to the reference proteins
seqs<-list(
  E=sapply(SARS_sel[which(nchar(SARS_sel)==proteome_summary$length[which(proteome_summary$protein_names=="E")])],function(i){
    str_split(i,"",simplify = T)
  }),
  
  N=sapply(SARS_sel[which(nchar(SARS_sel)==proteome_summary$length[which(proteome_summary$protein_names=="N")])],function(i){
    str_split(i,"",simplify = T)
  }),
  
  M=sapply(SARS_sel[which(nchar(SARS_sel)==proteome_summary$length[which(proteome_summary$protein_names=="M")])],function(i){
    str_split(i,"",simplify = T)
  }),
  
  S=sapply(SARS_sel[which(nchar(SARS_sel)==proteome_summary$length[which(proteome_summary$protein_names=="S")])],function(i){
    str_split(i,"",simplify = T)
  }))

#get references for structural protein sequences
ref_sel<-fasta[str_which(names(fasta),"E|N|M|S")]

#split each protein into individual amino acids
ref_seqs<-lapply(ref_sel,function(i){
  data.frame(seq=str_split(i,"",simplify = T),pos=1:nchar(i))
})

#name the list of reference sequences
names(ref_seqs)<-names(ref_sel)


all_counts<-fread(paste0(data_path,"all_peptides_prefilter.csv"))
filtered_peps<-all_counts%>%slice(which(prob<=.05))
all_filtered<-all_counts%>%slice(which(prob<=.05))%>%slice(-which(duplicated(peptide)))
all_raw<-fread(paste0(data_path,"all_peptides_with_pos.csv"))
all_raw$HLA<-str_remove_all(all_raw$HLA,"HLA-|\\*")
#get peptide pos
pep_pos<-all_raw%>%slice(which(peptide%in%all_filtered$peptide))%>%select(peptide,pos,gene)
#all repeated peptides
prot_peps<-filtered_peps%>%merge(pep_pos)%>%select(peptide,gene,pos,HLA)%>%unique()


prot_peps$peptide<-as.character(prot_peps$peptide)
cols_prot<-c("#F8766D","#00BFC4","#7CAE00","#C77CFF")
names(cols_prot)<-c("E","N","M","S")
prot_aligned<-lapply(c("E","N","M","S"),function(i){
  sel<-prot_peps[which(prot_peps$gene==i),]
  out<-sapply(1:nrow(sel),function(j){
    tmp<-rep(0,proteome_summary$length[which(proteome_summary$protein_names==i)]) 
    #tmp[sel[j,"pos"]:(nchar(sel$peptide[j])+sel[j,"pos"]-1)]<-str_split(sel$peptide[j],"",simplify = T)
    tmp[sel[j,"pos"]:(nchar(sel$peptide[j])+sel[j,"pos"]-1)]<-1
    tmp
    
  })
} )
names(prot_aligned)<-c("E","N","M","S")

al_2<-lapply(1:length(prot_aligned),function(j){
  i=prot_aligned[[j]]
  df<-data.frame(pos=1:nrow(i),count=apply(i,1,sum))
  df$prot<-names(prot_aligned)[j]
  df
})

al_2<-do.call(rbind,al_2)



align_plots<-lapply(1:length(prot_aligned),function(j){
  i=prot_aligned[[j]]
  df<-data.frame(pos=1:nrow(i),count=apply(i,1,sum))
  df$title<-names(seqs)[j]
  
  ggplot(df,aes(x=pos,y=count))+geom_line(color=cols_prot[j])+theme_calc()+theme(axis.title.x = element_blank())+ylab("# unique peptides")+scale_y_continuous(breaks = c(0:10),limits = c(0,10))+theme_pubclean()+facet_grid(.~title)
})


poly_count<-0
poly_plots<-lapply(1:length(seqs),function(j){
  i=seqs[[j]]
  df<-data.frame(pos=1:nrow(i),poly=apply(i,1,function(w){length(table(w))}))
  #print(sum(df$poly>1)/nrow(df))
  df<-data.frame(pos=1:nrow(i),poly=t(apply(i,1,function(w){
    tmp<-table(w)
    if(sum(names(tmp)=="X")>0){
    tmp<-tmp[-str_which(names(tmp),"X")]
    }
    c(length(tmp),max(tmp)/sum(tmp),sum(tmp[which(tmp!=max(tmp))]))
    })))
  #print(sum(df$poly.1>1)/nrow(df))
  print(sum(df$poly.3)/length(i))
  ggplot(df,aes(x=pos,y=poly.1-1))+geom_line(stat="identity",color=cols_prot[j])+ylab("polymorphisms")+xlab("positions")+theme_pubclean()+scale_y_continuous(breaks = c(0:3),limits = c(0,3))
  # ggarrange(
  # ggplot(df,aes(x=poly.2))+geom_density()+ylab("density")+xlab("% conservation")+theme_pubclean(),ncol = 1)
})




a<-poly_plots[[1]]+theme(text = element_text(size=5))
b<-poly_plots[[2]]+ylab("")+theme(text = element_text(size=7))
c<-poly_plots[[3]]+ylab("")+theme(text = element_text(size=7))
d<-poly_plots[[4]]+ylab("")+theme(text = element_text(size=7))

poly<-ggarrange(a,b,c,d,nrow = 1)



a<-align_plots[[1]]+theme(plot.title = element_text(hjust = 0.5))+xlab("")+theme(text = element_text(size=7))
b<-align_plots[[2]]+ylab("")+theme(plot.title = element_text(hjust = 0.5))+xlab("")+theme(text = element_text(size=7))
c<-align_plots[[3]]+ylab("")+theme(plot.title = element_text(hjust = 0.5))+xlab("")+theme(text = element_text(size=7))
d<-align_plots[[4]]+ylab("")+theme(plot.title = element_text(hjust = 0.5))+xlab("")+theme(text = element_text(size=7))

align<-ggarrange(a,b,c,d,nrow = 1)

poly_count<-ggarrange(align,poly,nrow = 2)
ggsave(poly_count,filename = "~/Covid-19/plots/main_figures/Figure_4_polymorphic_positions.pdf",width = 12,height = 6)

