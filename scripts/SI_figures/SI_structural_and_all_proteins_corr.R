data_path="~/Covid-19/EnsembleMHC-Covid19/datasets/"
Ensemble_PATH="~/Covid-19/EnsembleMHC-Covid19"
library(ggthemes)
library(ggpubr)
library(ggplot2)
library(data.table)
library(stringr)
min_norm<-function(x){(max(x)-x)/(max(x)-min(x))}
std<-function(x){(x-mean(x))/sd(x)}
library(latex2exp)
library(lubridate)
library(tidyr)
library(viridis)
library(ggrepel)
library(reshape2)
library(pwr)
library(wpp2019)
library(dplyr)
library(parallel)
library(wpp2019)
library(ggrepel)

##functions

#EnsembleEMP score function
EMP_score<-function(population){
  #return dataframe
  data.frame(
    #bind dataframes by row
    do.call(rbind,
            #iterate through populations
            lapply(population,function(df){
              #iterate through select alleles in that population
              #note: because of the 95% allele threshold, there will be 50-52 alleles
              #this is accounted for by dividing by the number of alleles
               apply(t(sapply(as.character(str_remove(df$HLA,"\\*")),function(i){
      
                 #calculate the EMP score based on all proteins or on structural proteins 
      sapply(c("All proteins","Structural proteins"),function(j){
        #pf is the peptide fraction 
        pf=HLA_count$pep_frac[which(HLA_count$HLA==i&HLA_count$source==j)]
        #score is the peptide fraction weighted by the normalized expression of that allele
        score=(df$pop_freq[which(df$HLA==i)]*pf)
        #if the co
        #if(length(score)==0){0}else{score}
        
      })
    })),2,sum)/length(df$HLA)
  })),stringsAsFactors = F)
  
}

#get individual day data for box plots
get_day_data<-function(day,death_threshold,MIN_countries){
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
  
  #return day specific data
  score_and_death_pop%>%slice(which(days==day))
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




#load population data
data("pop")
#set algorithm names
algos<-c("mhcflurry_affinity_percentile","mhcflurry_presentation_score","MixMHCpred","netMHC_affinity","netMHCpan_EL_affinity","netstab_affinity","pickpocket_affinity")

#this is taken from the JHU github page: https://github.com/CSSEGISandData/COVID-19
#load JHU data 1/22 - 4/10
all_data_clean<-fread(paste0(data_path,"JHU_data_clean.csv" ))


#day correction for when the data ends.
day_correction<-as.numeric(today()-mdy("4/10/2020"))

#read in the matrix of all proteins and structural proteins that were calculated previously
df_all_proteins<-read.csv(paste0(data_path,"identified_peptides_all_proteins_summarise_HLA_protein_counts.csv"))
df_struct<-read.csv(paste0(data_path,"identified_peptides_all_structural_proteins_only_summarise_HLA_protein_counts.csv"))

#creates the count matrix for peptide scoring. This is the total number of peptides predicted for each allele with respect to protein set
#the peptide frraction is also calculated
HLA_count<-rbind(df_all_proteins%>%mutate(source="All proteins")%>%
                   select(HLA,count=total_count_all_proteins,source),df_struct%>%mutate(source="Structural proteins")%>%
                   select(HLA,count=total_count_struct,source))%>%group_by(source)%>%
  mutate(pep_frac=count/sum(count))


#read HLA data
all_HLA_data<-fread(paste0(data_path,"HLA_data_manuscript.csv"),header = T)
#countries with at least 1 case reported case by 04/10
test_set<-unique(all_data_clean$Country.Region[which(all_data_clean$Deaths>=1)])
#Duplicate variable
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
    group_by(HLA)%>%summarise(allele_count=sum(eff_pop))%>% #group aggreate by allele
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
select_population<-pop%>%select(name,"2020")

#fix the offending names so that they match with the population data 
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

#df_itter_Death<-do.call(rbind,itter)
#create data frame from death threshold correlations
#combine results into a matrix
df_Death<-do.call(rbind,death_threshold_iter)
#normalize days for visualization purposes
norm_df_Death<-do.call(rbind,lapply(death_threshold_iter,function(d){
  d$non_norm<-
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

#table(norm_df_Death$thes[which(norm_df_Death$source=="Structural Proteins")])[2]/sum(table(norm_df_Death$thes[which(norm_df_Death$source=="Structural Proteins")]))

death_plot<-ggplot(norm_df_Death,aes(days,correlation,group=confirmed))+geom_line(aes(color=confirmed))+
  scale_color_viridis()+
  geom_point(data = sig,aes(days,correlation),color="red",size=.25)+
  facet_wrap(source~.)+theme_linedraw()+theme(legend.position = "bottom")+labs(color="number of deaths at day 0")+guides(colour=guide_colorbar(title.position = "right"))+ylab(TeX("EMP score-deaths per million correlation ($\\rho$)"))+
  xlab("normalized days")+theme(legend.text = element_text(size = 10),legend.key.size = unit(1,"cm"),
                                axis.text = element_text(size=12,face = "bold"),
                                title = element_text(size = 15,face="bold"),axis.title.y = element_text(size = 15,face="bold"))



df_spearman_pwr<-norm_df_Death
df_spearman_pwr$pwr<-apply(df_spearman_pwr,1,function(w){
  
  pwr.r.test(r = as.numeric(w[3]),n = as.numeric(w[6]),sig.level = .05)[[4]]
})


R=1
df_spearman_pwr$PPV<-(df_spearman_pwr$pwr*R)/((df_spearman_pwr$pwr*R)+df_spearman_pwr$p_value)

df<-df_spearman_pwr%>%select(source,sig,PPV)%>%group_by(source)%>%
  summarise(prop_p=table(sig)[2]/sum(table(sig)),prop_ppv=table(PPV>=.95)[2]/sum(table(PPV>=.95)))%>%melt()

#plot data
plot<-ggplot(df,aes(x=variable,y=value))+geom_bar(stat="identity",aes(fill=source))+theme_classic()+
  facet_wrap(source~.)+theme(legend.position = "none")

garg <- ggarrange(death_plot,plot,widths = c(.66,.33))

ggsave(filename = paste0(Ensemble_PATH,"/plots/SI_figures/SI_both_strucutral_n_full_prot.pdf"),garg,width = 14,height = 10.667)