setwd("~/Covid-19/new_mutations/")
library(dplyr)
library(stringr)
library(data.table)
library(parallel)
library(ggplot2)
pdf_dir <- "~/Covid-19/EnsembleMHC-Covid19/revision_requests/revision_figures/pdf/"
# functions ---------------------------------------------------------------
pt2in <- function(x) {x*0.0138889}

min_norm <- function(x) {
  (max(x) - x) / (max(x) - min(x))
}
std <- function(x) {
  (x - mean(x)) / sd(x)
}



# Read in mutation data ---------------------------------------------------

mutations_COVID_ep <- do.call(rbind, lapply(list.files(pattern = "mutation"), function(i) {
  data.frame(read.csv(i), gene = str_remove(i, "_mutation.csv"))
}))

# out of 102148 seqs 

# rename and select columns
mutations_COVID_ep <- mutations_COVID_ep %>% select(
  pos = Position,
  ref_AA = Consensus.residue,
  mut_num = Number.of.sequences.with.mutation..out.of.102148.,
  mutation = Residue.s..after.mutation, mut_freq = Number.of.sequences.corresponding.to.each.residue.after.mutation,
  gene
) 


# make heat map plot ------------------------------------------------------

mutations_COVID_ep$log<-log10(mutations_COVID_ep$mut_num)
mutations_COVID_ep$log[which(mutations_COVID_ep$log== -Inf)]<- 0
mutations_COVID_ep$group<-"0"
mutations_COVID_ep$group[which(mutations_COVID_ep$log<=1&&mutations_COVID_ep$log>0)]<-"1-10"
mutations_COVID_ep$group[which(mutations_COVID_ep$log<=2&mutations_COVID_ep$log>1)]<-"10-100"
mutations_COVID_ep$group[which(mutations_COVID_ep$log<=3&mutations_COVID_ep$log>2)]<-"100-1,000"
mutations_COVID_ep$group[which(mutations_COVID_ep$log<=4&mutations_COVID_ep$log>3)]<-"1,000-10,000"
mutations_COVID_ep$group[which(mutations_COVID_ep$log>4)]<-">10,000"
mutations_COVID_ep$group<-factor(mutations_COVID_ep$group,levels = c("0","1-10","10-100","100-1,000","1,000-10,000",">10,000")) 

ggthemr::ggthemr("fresh")

grobs<-lapply(c("E","N","M","S"),function(s){
mutations_COVID_ep %>%
  slice(which(gene == s)) %>%
  ggplot(aes(x = pos, y = 1, fill = group)) +
  geom_tile()+
  scale_fill_manual(values = wesanderson::wes_palette(name = "Zissou1", type = "continuous"))+
    ylab("")+
    theme(legend.position = "none",axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y =  element_blank())+
    xlab("position in sequence")
  
})

#just to get the legend
S<-mutations_COVID_ep %>%
  slice(which(gene == "S")) %>%
  ggplot(aes(x = pos, y = 1, fill = group)) +
  geom_tile()+
  scale_fill_manual(values = wesanderson::wes_palette(name = "Zissou1", type = "continuous"))+
  ylab("")+
  theme(legend.position = "bottom",axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y =  element_blank())+
  xlab("")



# replacements <- read.csv("replacements.csv") %>% slice(which(genomeRegion %in% c("S", "M", "N", "E")))
# insertions <- read.csv("insertions.csv") %>% slice(which(genomeRegion %in% c("S", "M", "N", "E")))
# deletions <- read.csv("deletions.csv") %>% slice(which(genomeRegion %in% c("S", "M", "N", "E")))


# map back to predicted peptides ------------------------------------------


#read in the positional data for all peptides
all_pos <- fread("~/Covid-19/EnsembleMHC-Covid19/datasets/all_peptides_with_pos.csv")

# make the mutation peptides -----

# get all structural peptides identified by ensembleMHC
wilson_peptides <- fread("~/Covid-19/EnsembleMHC-Covid19/datasets/data_not_transfered/all_peptides_prefilter.csv") %>%
  slice(which(gene %in% c("S", "M", "N", "E"))) %>%
  slice(which(prob <= 0.05)) %>% 
  select(peptide,gene) %>%
  unique()

# get the indexes of all the peptide
my_peps <- all_pos %>%
  slice(which(peptide %in% wilson_peptides$peptide)) %>%
  select(peptide, gene, pos) %>%
  unique() %>% mutate(end=pos+(nchar(peptide)-1))

# split each of the identified peptides and pair with its index 
peptide_indexs <- do.call(rbind, lapply(1:nrow(my_peps), function(i) {
  tmp <- my_peps[i, ]
  data.frame(pos = tmp$pos:tmp$end, 
             ref_AA = as.vector(str_split(tmp$peptide, pattern = "", simplify = T)), 
             gene = tmp$gene,
             peptide=tmp$peptide,
             start = tmp$pos,
             end = tmp$end)
})) %>% unique()

#get only mutations that corespond with the 108 peptide
peptides_with_mutations <- mutations_COVID_ep[-which(is.na(prodlim::row.match(mutations_COVID_ep[, c(1, 6)], peptide_indexs[, c(1, 3)]))), ]

#merge back with  mutation data 
peptides_with_mutations <- merge(peptide_indexs, peptides_with_mutations)


peps <- unique(peptides_with_mutations$peptide)

#this will create all mutated peptides
all_mutated_peptides <- do.call(rbind, lapply(peps, function(i) {
  tmp <- peptides_with_mutations %>%
    slice(which(peptide == i)) %>%
    arrange(pos)

  do.call(rbind, lapply(1:length(tmp$mutation), function(j) {
    mut_AA <- trimws(as.vector(str_split(tmp$mutation[j], pattern = ",", simplify = T)))
    mut_freq <- as.numeric(trimws(as.vector(str_split(tmp$mut_freq[j], pattern = ",", simplify = T))))/102148
    if (any(mut_AA != "")) {
      t(sapply(1:length(mut_AA), function(k) {
        pep_mut <- tmp
        pep_mut$ref_AA[j] <- mut_AA[k]
        data.frame(mutate_pep = paste(pep_mut$ref_AA, collapse = ""), mut_freq = mut_freq[k], parent_pep = unique(pep_mut$peptide), gene = unique(pep_mut$gene), mut_pos = j)
      })) %>% data.frame()
    }
  }))
}))

#fixes data class issues 
all_mutated_peptides <- apply(all_mutated_peptides, 2, unlist) %>% data.frame(., stringsAsFactors = F)
all_mutated_peptides$mut_freq <- as.numeric(all_mutated_peptides$mut_freq)
all_mutated_peptides$mut_pos <- as.numeric(all_mutated_peptides$mut_pos)

for_pred <- all_mutated_peptides %>% select(peptide=mutate_pep,gene) %>% mutate(length=nchar(peptide))
load("~/Covid-19/EnsembleMHC-Covid19/datasets/52_HLA_list.R")


#write data for predictions and use EnsembleMHC 
# lapply(sel_HLA, function(i) {
#   write.csv(for_pred, file = paste0("pred_ref/","HLA-", i, "_ready_for_predictions.csv"),row.names = F)
# })


# plot the position of EnsembleMHC predicted peptides ---------------------

endpoint<-mutations_COVID_ep %>% group_by(gene) %>% summarise(end=max(pos))
peptide_indexs

grobs_lines<-lapply(c("E","N","M","S"),function(s){
  peptide_indexs %>%
    slice(which(gene == s)) %>%
    group_by(pos) %>% summarise(count=length(pos))%>%
    ggplot(aes(x = pos, y = as.integer(count))) +
    geom_bar(stat="identity",fill="black")+
    theme_classic()+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())+
    ylab("count")+
    xlab("")+coord_cartesian(xlim =c(1,endpoint$end[which(endpoint$gene==s)]))
  
})

p1<-grobs_lines[[1]]/grobs[[1]] +plot_layout(heights = c(.9,.1))&theme(plot.margin = margin(0,0,0,0))
p2<-grobs_lines[[2]]/grobs[[2]] +plot_layout(heights = c(.9,.1))&theme(plot.margin = margin(0,0,0,0))
p3<-grobs_lines[[3]]/grobs[[3]] +plot_layout(heights = c(.9,.1))&theme(plot.margin = margin(0,0,0,0))
p4<-grobs_lines[[4]]/grobs[[4]] +plot_layout(heights = c(.9,.1))&theme(plot.margin = margin(0,0,0,0))

ggsave(p1,filename = paste0(pdf_dir,"Envelope_protein_mutations.pdf"),width = pt2in(900),height = pt2in(150))
ggsave(p2,filename = paste0(pdf_dir,"nuc_protein_mutations.pdf"),width = pt2in(900),height = pt2in(150))
ggsave(p3,filename = paste0(pdf_dir,"mem_protein_mutations.pdf"),width = pt2in(900),height = pt2in(150))
ggsave(p4,filename = paste0(pdf_dir,"spike_protein_mutations.pdf"),width = pt2in(900),height = pt2in(150))
# analysis of the results from predictions ----------

#read in teh prediction data 
pred_peptides <- do.call(rbind, lapply(list.files("transfer_dir", full.names = T), function(i) {
  read.csv(i)
}))

load("~/Covid-19/EnsembleMHC-Covid19/datasets/P_sum_median_1000_boot.R")

wilson_peptides <- fread("~/Covid-19/EnsembleMHC-Covid19/datasets/data_not_transfered/all_peptides_prefilter.csv") %>%
  slice(which(gene %in% c("S", "M", "N", "E"))) %>%
  slice(which(prob <= 0.05)) %>% 
  select(peptide,gene,HLA) %>%
  unique()



#calculate peptide probs
pep_probs <- lapply(unique(pred_peptides$HLA), function(w) {
  # create a tmp variable with that consists of the corona virus predictions for one allele
  tmp <- pred_peptides %>% dplyr::slice(which(pred_peptides$HLA == w))
  # normalize the scores for the presentation score and pickpocket
  # both of these scores are not percentiles and the identifed thresholds were based on the similarly normalized scores
  tmp$mhcflurry_presentation_score <- min_norm(tmp$mhcflurry_presentation_score)
  tmp$pickpocket_affinity <- min_norm(tmp$pickpocket_affinity)
  # coverent the gene name into a factor
  tmp$gene <- factor(tmp$gene)
  # print current HLA being processed
  print(w)
  # create the prob list which is the selection of peptides that fall within the score filter for one algorithm
  # assign the algorithm FDR to each peptide based on the benchmarking calculations
  # loop through all 7 algorithms
  prob_list <- lapply(colnames(tmp)[colnames(tmp) %in% unique(P_sum$algo)], function(q) {
    # Originally, the algorithms were assigned PPVs. These PPVs are converted to FDR through the relation FDR = 1 - PPV
    # therefore, the following line finds the PPV for the allele specified by w and algorithm q
    neg <- 1 - P_sum$PPV[which(P_sum$HLA == str_remove(w, pattern = "HLA-") & P_sum$algo == q)]
    # same as above but with the score threshold for teh w and q combo
    thres <- P_sum$value[which(P_sum$HLA == str_remove(w, pattern = "HLA-") & P_sum$algo == q)]
    # select all peptides that fall within the scoring threshold for that algorithm at that allele
    peptides <- tmp$peptide[which(tmp[, q] <= thres)]
    if(length(peptides)>0){
    # return a data frame consisting of selected peptides, the algorithm FDR, and algorithm name
    data.frame(peptide = peptides, prob = neg, algo = q)
    }
  })
  # combine all of the prob_list elements in one dataframe, calculate the products of the FDRs associatied with detecting algorithms
  # merge with information regarding the gene and HLA
  prob_combo <- do.call(rbind, prob_list) %>%
    data.frame() %>%
    group_by(peptide) %>%
    summarise(prob = prod(prob)) %>%
    merge(tmp[, c("peptide", "gene", "HLA")])
})

all_new_preds <- do.call(rbind, pep_probs)

colnames(all_new_preds)[1] <- "mutate_pep"

#merge scored peptides with the mutate peptide data
merged_all <- merge(all_new_preds, all_mutated_peptides,all = T)

#score each peptide to see if the mutant would pass the score filter
all_peptides_with_score <- merged_all %>%
  arrange(prob) %>%
  slice(-which(duplicated(mutate_pep))) %>%
  mutate(FILTER_PASS = prob <= .05)

all_peptides_with_score$FILTER_PASS[which(is.na(all_peptides_with_score$FILTER_PASS))]<-FALSE

all_peptides_with_score$HLA<-str_remove(all_peptides_with_score$HLA,"HLA-")

matched_peptides<- all_peptides_with_score[!is.na(prodlim::row.match(all_peptides_with_score%>% select(HLA,peptide=parent_pep) %>% unique(),wilson_peptides%>%select(HLA,peptide)%>%data.frame())),]

# plot the results of the mutation analysis -------------------------------
p5<-table(matched_peptides$FILTER_PASS) %>%data.frame() %>% ggplot(aes(x=Var1,y=Freq,fill=Var1))+geom_bar(stat="identity")+theme(legend.position = "none")
ggsave(p5,filename = paste0(pdf_dir,"binding_of_mutant_peptide.pdf"),width = pt2in(300),height = pt2in(150))
p6<- matched_peptides%>% ggplot(aes(mut_freq,fill=FILTER_PASS))+geom_histogram()
ggsave(p6,filename = paste0(pdf_dir,"freq_of_mutants.pdf"),width = pt2in(600),height = pt2in(150))
all_peptides_with_score %>% ggplot(aes(mut_freq,fill=FILTER_PASS))+geom_histogram(binwidth = .00001)


