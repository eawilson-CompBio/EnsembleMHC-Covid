if(!require("dplyr")){install.packages("dplyr")}
if(!require("stringr")){install.packages("stringr")}
cl_inputs <- commandArgs(trailingOnly = T)
pred_file <- read.csv(cl_inputs[1])

P_sum <- read.csv(cl_inputs[2])


#scale values with 0 being the max value and 1 being the min
min_norm <- function(x) {
  (max(x) - x) / (max(x) - min(x))
}

pep_probs <- do.call(rbind,lapply(unique(pred_file$HLA), function(w) {
  tmp <- pred_file %>% dplyr::slice(which(pred_file$HLA == w))

  # normalize the scores for the presentation score and pickpocket
  # both of these scores are not percentiles and the identifed thresholds were based on the similarly normalized scores
  tmp$mhcflurry_presentation_score <- min_norm(tmp$mhcflurry_presentation_score)
  tmp$pickpocket_affinity <- min_norm(tmp$pickpocket_affinity)
  # coverent the gene name into a factor
  tmp$gene <- factor(tmp$protein)
  # print current HLA being processed
  print(w)
  # create the prob list which is the selection of peptides that fall within the score filter for one algorithm
  # assign the algorithm FDR to each peptide based on the benchmarking calculations
  # loop through all 7 algorithms
  prob_list <- lapply(colnames(tmp)[colnames(tmp) %in% unique(P_sum$algo)], function(q) {
    # Originally, the algorithms were assigned PPVs. These PPVs are converted to FDR through the relation FDR = 1 - PPV
    # therefore, the following line finds the PPV for the allele specified by w and algorithm q
    neg <- 1 - P_sum$PPV[which(P_sum$HLA == str_remove_all(w, pattern = "\\*|HLA-") & P_sum$algo == q)]
    # same as above but with the score threshold for teh w and q combo
    thres <- P_sum$value[which(P_sum$HLA == str_remove_all(w, pattern = "\\*|HLA-") & P_sum$algo == q)]
    # select all peptides that fall within the scoring threshold for that algorithm at that allele
    peptides <- tmp$peptide[which(tmp[, q] <= thres)]
    # return a data frame consisting of selected peptides, the algorithm FDR, and algorithm name
    if(length(peptides)){ data.frame(peptide = peptides, prob = neg, algo = q) }
  })

  # combine all of the prob_list elements in one dataframe, calculate the products of the FDRs associatied with detecting algorithms
  # merge with information regarding the gene and HLA
  prob_combo <- do.call(rbind, prob_list) %>%
    data.frame() %>%
    group_by(peptide) %>%
    summarise(peptideFDR = prod(prob)) %>%
    merge(tmp[, c("peptide", "protein", "HLA")])
}))

pred_peps <- pred_file %>% left_join(pep_probs,by = c("peptide","protein","HLA"))

pred_peps$peptideFDR[which(is.na(pred_peps$peptideFDR))]<-1

HLA<-str_remove(cl_inputs[1],"_scored_peptides.csv")

write.csv(file=paste0(HLA,"_peptideFDR_pred.csv"),pred_peps,row.names = F)
