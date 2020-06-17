data_path="~/Covid-19/EnsembleMHC-Covid19/datasets/"
Ensemble_PATH="~/Covid-19/EnsembleMHC-Covid19"
library(dplyr)
library(data.table)
library(stringr)
min_norm<-function(x){(max(x)-x)/(max(x)-min(x))}
library(latex2exp)
library(lubridate)
library(viridis)
library(wpp2019)
library(parallel)
library(wpp2019)
library(ggplot2)

data("pop")
#functions

#EnsembleEMP score function
EMP_score<-function(population){
  #return dataframe
  data.frame(
    #bind dataframes by row
    do.call(rbind,
            #iterate through population list
            lapply(population,function(df){
              #iterate through select alleles in that population
              #note: because of the 95% allele threshold, there will be 50-52 alleles fpr each country
              #this is accounted for by dividing by the number of alleles
              apply(t(sapply(as.character(str_remove(df$HLA,"\\*")),function(i){
                
                #calculate the EMP score based on all proteins or on structural proteins 
                sapply(c("All proteins","Structural proteins"),function(j){
                  #pf is the peptide fraction 
                  pf=HLA_count$pep_frac[which(HLA_count$HLA==i&HLA_count$source==j)]
                  #score is the peptide fraction weighted by the normalized expression of that allele
                  score=(df$pop_freq[which(df$HLA==i)]*pf)
                })
              })),2,sum)/length(df$HLA)
            })),stringsAsFactors = F)
  
}

#calculate correlation as a function of time
death_threshold_specific_corr<-function(death_threshold,MIN_countries){
  #filter at a set number of deaths
  #set days for time series analysis
  #add days since it reached that point
  day_data<-CoV_data_sel_countries%>%slice(which(Deaths>=death_threshold))%>%
    group_by(country)%>%mutate(days=as.numeric((today()-day_correction)-mdy(date)))%>%mutate(days=max(days)-days)
  
  
  #set the threshold of days so you know when to stop iterating once the number of countries decreases below threshold
  thres<-table(day_data$days)[which(table(day_data$days)>=MIN_countries)]
  thres<-as.numeric(names(thres)[length(thres)])
  
  #mere the day data with the population data
  score_and_death_pop<-day_data%>%merge(select_population)
  
  #merge EMP scores and calculate deaths per million
  score_and_death_pop<-score_and_death_pop%>%group_by(country,days)%>%mutate(pop=`2020`/1000)%>%
    mutate(death_per_pop=Deaths/pop)%>%merge(countryEMPscore)%>%
    select(country,death_per_pop,Deaths,days,All.proteins,Structural.proteins)%>%
    slice(which(Deaths>=death_threshold))
  
  cor_data<-lapply(1:thres,function(i){
    #create subset of data at day i
    tmp<-score_and_death_pop%>%slice(which(days==i))
    list(cor_struct=cor(tmp$Structural.proteins,tmp$death_per_pop,method = "spearman"),
         cor_all_prot=cor(tmp$All.proteins,tmp$death_per_pop,method = "spearman"),
         tab=tmp,
         day=i,
         p_struct=cor.test(tmp$Structural.proteins,tmp$death_per_pop,method = "spearman")[[3]],
         p_all_prot=cor.test(tmp$All.proteins,tmp$death_per_pop,method = "spearman")[[3]])
  })
  
  p_val_struct<-sapply(1:thres,function(w){cor_data[[w]]$p_struct})
  p_val_all_prot<-sapply(1:thres,function(w){cor_data[[w]]$p_all_prot})
  #make a matrix of protein group specific correlations
  df_corr<-data.frame(days=1:thres,correlation_all_prot=sapply(cor_data,function(w){w$cor_all_prot}),correlation_struct=sapply(cor_data,function(w){w$cor_struct}))
  
  #create the p value matrix
  pvl_mat<-rbind(data.frame(days=df_corr$days,source="Structural_proteins",p_value=p_val_struct),
                 data.frame(days=df_corr$days,source="All_proteins",p_value=p_val_all_prot))
  
  #renames columns and melt matrix
  colnames(df_corr)[2:3]<-c("All_proteins","Structural_proteins")
  melted_df<-melt(df_corr,id.vars = c("days"))
  colnames(melted_df)[2:3]<-c("source","correlation")
  
  #merge with p value matrix and identify significant correlations
  melted_df<-melted_df%>%merge(pvl_mat)
  melted_df$sig<-0
  melted_df$sig[which(melted_df$p_value<=.05)]<-1
  
  #find the number of countries at each time point regradless of the minimum country threshold
  tab_melt<-table(day_data$days)
  #merge so only considered days meeting the min country remain
  melted_df<-melted_df%>%merge(data.frame(days=as.numeric(names(tab_melt)),num_country=as.numeric(tab_melt)))
  
  
  #  melt_df_pre$confirmed<-death_threshold
  melted_df$confirmed<-death_threshold
  #melt_df_pre
  #return matrix
  melted_df
}

#set algorithms vector
algos<-c("mhcflurry_affinity_percentile","mhcflurry_presentation_score","MixMHCpred","netMHC_affinity","netMHCpan_EL_affinity","netstab_affinity","pickpocket_affinity")

#load the P_sum matrix. This is teh stored algorithm and allele specific score and FDR
load(paste0(data_path,"P_sum_median_1000_boot.R"))

#day correction for when the data ends
day_correction<-as.numeric(today()-mdy("4/10/2020"))

#read in the coronavirus predictions 
corona_data<-fread(paste0(data_path,"predicted_corona_peptides.csv"))


#calculate peptide probabilities 
#start by iterating through every unique HLA
pep_probs<-lapply(unique(corona_data$HLA),function(w){
  #create a tmp variable with that consists of the corona virus predictions for one allele
  tmp <- corona_data%>%dplyr::slice(which(corona_data$HLA==w))
  #normalize the scores for the presentation score and pickpocket
  #both of these scores are not percentiles and the thresholds were bsaed on the normalized scores
  tmp$mhcflurry_presentation_score<-min_norm(tmp$mhcflurry_presentation_score)
  tmp$pickpocket_affinity<-min_norm(tmp$pickpocket_affinity)
  #coverent the gene name into a factor
  tmp$gene<-factor(tmp$gene)
  #print current HLA being processed
  print(w)
  #create the prob list which is the selection of peptides that fall within the score filter for one algorithm
  #assign the algorithm FDR to each peptide based on the benchmarking calculations
  prob_list<-lapply(colnames(tmp)[colnames(tmp)%in%unique(P_sum$algo)],function(q){
    #Originally, the algorithms were assigned PPVs. These PPVs are converted to FDR through the relation PPV = 1 - FDR
    neg=1-P_sum$PPV[which(P_sum$HLA==str_remove(w,pattern = "\\*")&P_sum$algo==q)]
    #the score threshold required for that algorithm at the HLA is recovered from P_sum matrix
    thres=P_sum$value[which(P_sum$HLA==str_remove(w,pattern = "\\*")&P_sum$algo==q)]
    #select all peptides that fall within the scoring threshold for that algorithm at that allele
    peptides<-tmp$peptide[which(tmp[,q]<=thres)]
    #return a data frame consisting of selected peptides, teh algorithm FDR, and algorithm name
    data.frame(peptide=peptides,prob=neg,algo=q)
  })
  #combine all of the prob_lists in one dataframe, calculate the products of the FDRs associatied with detecting algorithms
  #merge with information regarding the gene and HLA
  prob_combo<-do.call(rbind,prob_list)%>%data.frame()%>%
    group_by(peptide)%>%summarise(prob=prod(prob))%>%merge(tmp[,c("peptide","gene","HLA")])
  
  
  
})

#bind all the pep_prob list generated in the previous step
all_counts<-do.call(rbind,pep_probs)

#make vector of all unique alleles
sel_alleles<-unique(all_counts$HLA)

#probabiliyty threshold
prob_threshold=.05

#scramble HLA
all_counts$HLA<-as.character(all_counts$HLA)
all_counts$HLA<-sample(all_counts$HLA)

pass_filter<-all_counts%>%slice(which(prob<=.05))

pass_filter$source<-"All proteins"
pass_filter$source[which(pass_filter$gene%in%c("N","M","E","S"))]<-"Structural proteins"

HLA_count<-pass_filter%>%group_by(HLA,source)%>%summarise(count=length(peptide))%>%group_by(source)%>%mutate(pep_frac=count/sum(count))

HLA_count<-do.call(rbind,lapply(c("Structural proteins","All proteins"),function(j){
  tmp <- HLA_count[which(HLA_count$source==j),]%>%data.frame()
  if(nrow(tmp)!=52){
  rbind(tmp,data.frame(HLA=sel_alleles[-which(sel_alleles%in%tmp$HLA)],count=0,pep_frac=0,source=j))
  }else{tmp}
}))


#this is taken from the JHU github page: https://github.com/CSSEGISandData/COVID-19
#load JHU data 1/22 - 4/10
all_data_clean<-fread(paste0(data_path,"JHU_data_clean.csv" ))


#read HLA data
all_HLA_data<-fread(paste0(data_path,"HLA_data_manuscript.csv"),header = T)
#countries with at least 1 case reported case by 04/10
test_set<-unique(all_data_clean$Country.Region[which(all_data_clean$Deaths>=1)])
#save country names
country_list<-test_set

#convert to list and fix country names or add alternative names
test_set<-lapply(test_set,function(w){w})
test_set[[which(country_list=="Czechia")]]<-"Czech"
test_set[[which(country_list=="UK")]]<-c("United Kingdom|England|Wales$")
test_set[[which(country_list=="US")]]<-"USA"


#set allele vector to the 52 selected alleles
sel_HLA<-unique(HLA_count$HLA)

#normalize the allele frequencies with respect to selected alleles
allele_freq_list<-lapply(test_set,function(w){
  
  name=w
  print(name)
  w<-all_HLA_data[str_which(all_HLA_data$population,w),]
  colnames(w)<-c("HLA","pop","allele_freq","n")
  w$HLA<-str_remove(w$HLA,"\\*")
  total_n<-w%>%
    slice(which(HLA%in%sel_HLA))%>%
    mutate(country=name)%>%select(pop,n)%>%unique()%>%mutate(n_total=sum(n))%>%pull(n_total)%>%unique()
  
  
  tmp<-w%>%mutate(eff_pop=allele_freq*n*2)%>%     #convert to allele count
    slice(which(HLA%in%sel_HLA))%>%               #isolate the 52 alleles
    group_by(HLA)%>%summarise(allele_count=sum(eff_pop))%>% #group and aggreate allele count by 52 alleles
    mutate(pop_freq=allele_count/sum(allele_count))%>%      # normalize allele count
    mutate(country=name)%>%mutate(total_n=total_n)          #add country name and total number of HLA typed individual 
  
  
})

#name list
names(allele_freq_list)<-country_list

#get rid of empty lists
allele_freq_list<-allele_freq_list[which(sapply(allele_freq_list,function(w){nrow(w)>0}))]

#filter for populations with at least 1000 HLA typed individuals at 4 digit resolution
allele_freq_list<-allele_freq_list[which(sapply(allele_freq_list,function(w){unique(w$total_n)>1000}))]

#filter for number of select alleles
allele_freq_list<-allele_freq_list[which(sapply(allele_freq_list,nrow)/52>.95)]
#get coroan virus data that correspondes to the selected alleles
CoV_data_sel_countries<-all_data_clean[which(all_data_clean$Country.Region%in%names(allele_freq_list)),]



#apply the peptide score function 
countryEMPscore<-EMP_score(allele_freq_list)
#add a column for the names of the countries
countryEMPscore$country<-row.names(countryEMPscore)

#fix the column name so you can merge later 
colnames(CoV_data_sel_countries)[1]<-"country"




#make dataframe of 2020 populations by country 
#noteL: population is in thousands
select_population<-pop%>%select(name,"2020")

#some countries names  differ than the ones used in the  data
#This step changes names that differ to the names used in the  data
colnames(select_population)[1]<-"country"
select_population$country<-as.character(select_population$country)
select_population$country[str_which(select_population$country,"Korea")][2]<-"South Korea"
select_population$country[str_which(select_population$country,"Iran")]<-"Iran"
select_population$country[str_which(select_population$country,"Russia")]<-"Russia"
select_population$country[str_which(select_population$country,"United States of ")]<-"US"
select_population$country[str_which(select_population$country,"United Kingdom")]<-"UK"
select_population$country[str_which(select_population$country,"Taiw")]<-"Taiwan"
select_population$country[str_which(select_population$country,"Hong")]<-"Hong Kong"



death_threshold_iter<-mclapply(1:100,function(death_thres){death_threshold_specific_corr(death_thres,8)},mc.cores = 10)


#create data frame from death threshold correlations
#combine results into a matrix
df_Death<-do.call(rbind,death_threshold_iter)
#normalize days for visualization purposes
norm_df_Death<-do.call(rbind,lapply(death_threshold_iter,function(d){
  d$days<-((d$days)-min(d$days))/((max(d$days)-min(d$days)))
  d
}))

#create matrix of only significant correlations
#rename corr types
norm_df_Death$source<-as.character(norm_df_Death$source)
norm_df_Death$source[which(norm_df_Death$source=="Structural_proteins")]<-"SARS-CoV-2 structural proteins"
norm_df_Death$source[which(norm_df_Death$source=="All_proteins")]<-"entire SARS-CoV-2 proteome"
norm_df_Death$source<-as.factor(norm_df_Death$source)
sig<-subset(norm_df_Death,sig==1)
table(sig$source)/table(norm_df_Death$source)

ggplot(norm_df_Death,aes(days,correlation,group=confirmed))+geom_line(aes(color=confirmed))+
  scale_color_viridis()+
  geom_point(data = sig,aes(days,correlation),color="red",size=.5)+facet_wrap(source~.)+
  theme_linedraw()+theme(legend.position = "top")+labs(color="number of deaths by day 0")+guides(colour=guide_colorbar(title.position = "right"))+
  ylab(TeX("EMP score-deaths per million correlation ($\\rho$)"))+xlab("normalized days")


