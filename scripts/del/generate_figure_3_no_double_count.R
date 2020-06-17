#confirmed_case_threshold<-as.numeric(commandArgs(trailingOnly=T))
#library(Biostrings)
library(ggthemes)
library(ggpubr)
library(ggplot2)
#library(see)
library(data.table)
library(stringr)
#library(GGally)
#library(wesanderson)
min_norm<-function(x){(max(x)-x)/(max(x)-min(x))}
std<-function(x){(x-mean(x))/sd(x)}
library(latex2exp)
#library(ggridges )
library(lubridate)
#library(gridExtra)
library(tidyr)
library(viridis)
library(ggrepel)
library(reshape2)
#library(patchwork)
library(pwr)
library(wpp2019)
library(dplyr)
library(parallel)
#load population data
data("pop")
#set algorithm names
algos<-c("mhcflurry_affinity_percentile","mhcflurry_presentation_score","MixMHCpred","netMHC_affinity","netMHCpan_EL_affinity","netstab_affinity","pickpocket_affinity")


#generate coronavrius epidemilocial data
#this is taken from the JHU github page: https://github.com/CSSEGISandData/COVID-19
#the data is up april 10th 
case_data<-list.files(path = "~/Covid-19/JHU_data/4.10.corona/COVID-19/csse_covid_19_data/csse_covid_19_daily_reports/",pattern = "csv",full.names = T)
#Just as a heads up, they switched the date notation scheme about 60 days into there tracking. Here I process each indepedently then combine
all_data<-do.call(rbind,lapply(case_data[1:60],function(i){
  x<-read.csv(i,stringsAsFactors = F)
  tmp<-x[,c(2,4:6)]
  tmp<-tmp%>%group_by(Country.Region)%>%replace_na(replace = list(Confirmed=0,Deaths=0,Recovered=0 ))%>%summarise_each(sum)
  tmp$date<-str_extract(i,"[0-9]+-[0-9]+-2020")
  tmp
  
  
}))
#this for the second bit of data
all_data_cont<-do.call(rbind,lapply(case_data[61:length(case_data)],function(i){
  x<-read.csv(i,stringsAsFactors = F)
  tmp<-x[,c(4,8:10)]
  tmp<-tmp%>%group_by(Country_Region)%>%replace_na(replace = list(Confirmed=0,Deaths=0,Recovered=0 ))%>%summarise_each(sum)
  tmp$date<-str_extract(i,"[0-9]+-[0-9]+-2020")
  tmp
  
}))
#make sure column names are the same
colnames(all_data_cont)<-colnames(all_data)

#bind the list together
all_d<-data.frame(rbind(all_data,all_data_cont),stringsAsFactors = F)
#convert the data to character 
all_d$date<-as.character(all_d$date)

#fix troublesome country names so that they align with the names in the HLA data
all_d$Country.Region[str_which(all_d$Country.Region,"Bahamas")]<-"Bahamas"
all_d$Country.Region[str_which(all_d$Country.Region,"Korea")]<-"South Korea"
all_d$Country.Region[str_which(all_d$Country.Region,"Iran")]<-"Iran"
all_d$Country.Region[str_which(all_d$Country.Region,"China")]<-"China"
all_d$Country.Region[str_which(all_d$Country.Region,"United Kingdom")]<-"UK"
all_d$Country.Region[str_which(all_d$Country.Region,"Taiw")]<-"Taiwan"
all_d$Country.Region[str_which(all_d$Country.Region,"Cz")]<-"Czechia"

#this will address any problems arising from multiple entries from the same day
alldd<-all_d%>%group_by(Country.Region,date)%>%summarise_each(sum)%>%as.data.frame()

#alldd$Confirmed[is.na(alldd$Confirmed)]<-0
#alldd$Deaths[is.na(alldd$Deaths)]<-0

#calculate death rate
all_data_clean<-alldd%>%group_by(Country.Region)%>%arrange(date)#%>%mutate(death_rate=Deaths/Confirmed)

#rm(all_d)
#rm(all_data)
#rm(all_data_cont)
#rm(alldd)
#rm(case_data)


#day correction for when the data ends.
day_correction<-as.numeric(today()-mdy("4/10/2020"))

#read in the matrix of all proteins and structural proteins that were calculated previously
df_all_proteins<-read.csv("~/Covid-19/improtant_intermediate_datasets/identified_peptides_all_proteins_summarise_HLA_protein_counts.csv")
df_struct<-read.csv("~/Covid-19/improtant_intermediate_datasets/identified_peptides_all_structural_proteins_only_summarise_HLA_protein_counts.csv")

#creates the count matrix for peptide scoring. This is the total number of peptides predicted for each allele with respect to protein set
#the peptide frraction is also calculated
HLA_count<-rbind(df_all_proteins%>%mutate(source="All proteins")%>%
                   select(HLA,count=total_count_all_proteins,source),df_struct%>%mutate(source="Structural proteins")%>%
                   select(HLA,count=total_count_struct,source))%>%group_by(source)%>%
  mutate(pep_frac=count/sum(count))

#delete this
HLA_count<-rbind(HLA_count[which(HLA_count$source=="All proteins"),]%>%data.frame(),
HLA_count[which(HLA_count$source=="Structural proteins"),c("HLA","source")]%>%merge(total_syn,by="HLA")%>%select(HLA,count,source,pep_frac=frac)%>%data.frame())

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
        pf=HLA_count$count[which(HLA_count$HLA==i&HLA_count$source==j)]
        #score is the peptide fraction weighted by the normalized expression of that allele
        score=(df$pop_freq[which(df$HLA==i)]*pf)
        #if the co
        #if(length(score)==0){0}else{score}
        
      })
    })),2,sum)/length(df$HLA)
  })),stringsAsFactors = F)
  
}


all_HLA_data<-fread("~/Covid-19/HLA_data/4D_HLA_globaldata.ssv",stringsAsFactors = F,header = F)
#countries with at least 1 case reported case by 04/10
#make it so ther eis no double counting beacsue some countries partion by ethinc minority
all_HLA_data$V3[str_which(all_HLA_data$V3,"Germany")]<-"Germany"
all_HLA_data$V3[str_which(all_HLA_data$V3,"USA")]<-"USA"
all_HLA_data$V3[str_which(all_HLA_data$V3,"Israel")]<-"Israel"
confirmed_case_threshold<-1
test_set<-unique(all_data_clean$Country.Region[which(all_data_clean$Confirmed>=confirmed_case_threshold)])
#Duplicate variable
country_list<-test_set

#convert to list and fix country names or add alternative names
test_set<-lapply(test_set,function(w){w})
test_set[[which(country_list=="Czechia")]]<-"Czech"
test_set[[which(country_list=="UK")]]<-c("United Kingdom","England","Wales")
test_set[[which(country_list=="US")]]<-"USA"


#set allele vector to the 52 selected alleles
sel_HLA<-unique(HLA_count$HLA)


allele_freq_list<-lapply(test_set,function(w){
  if(length(w)>1){
    name=w[1]
    w=paste(w,collapse = "|")
    
  }else{
    name=w
  }
  
  print(name)
  w<-all_HLA_data[str_which(all_HLA_data$V3,w),]
  colnames(w)<-c("Rowname","HLA","pop","allele_freq","n")
  w$HLA<-str_remove(w$HLA,"\\*")
  total_n<-w%>%slice(which(HLA%in%sel_HLA))%>%mutate(country=name)%>%select(pop,n)%>%unique()%>%mutate(n_total=sum(n))%>%pull(n_total)%>%unique()
  
  
  tmp<-w%>%mutate(eff_pop=allele_freq*n*2)%>%     #convert to allele count
    #replace_na(list(eff_pop=0))%>%
    slice(which(HLA%in%sel_HLA))%>%               #isolate the 52 alleles
    group_by(HLA)%>%summarise(allele_count=sum(eff_pop))%>% #group aggreate by allele
    mutate(pop_freq=allele_count/sum(allele_count))%>%      # normalize allele count
    mutate(country=name)%>%mutate(total_n=total_n)          #add country name and total number of HLA typed individual 
  
  
})


allele_freq_list<-s
#name list
names(allele_freq_list)<-unique(all_data_clean$Country.Region[which(all_data_clean$Confirmed>=confirmed_case_threshold)])

#get rid of empty lists
allele_freq_list<-allele_freq_list[which(sapply(allele_freq_list,function(w){nrow(w)>0}))]

#filter for populations with at least 1000 HLA typed individuals at 4 digit resolution
#allele_freq_list<-allele_freq_list[which(sapply(allele_freq_list,function(w){unique(w$total_n)>200}))]

#filter for number of select alleles
allele_freq_list<-allele_freq_list[which(sapply(allele_freq_list,nrow)/52>.9)]
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
         day=w,
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


death_threshold_iter<-mclapply(1:100,function(death_thres){death_threshold_specific_corr(death_thres,6)},mc.cores = 10)

#df_itter_Death<-do.call(rbind,itter)
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
table(sig$source)/table(df_Death$source)
#table(norm_df_Death$thes[which(norm_df_Death$source=="Structural Proteins")])[2]/sum(table(norm_df_Death$thes[which(norm_df_Death$source=="Structural Proteins")]))

death_plot<-ggplot(norm_df_Death,aes(days,correlation,group=confirmed))+geom_line(aes(color=confirmed))+
  scale_color_viridis()+
  geom_point(data = sig,aes(days,correlation),color="red",size=.25)+
  facet_wrap(source~.)+theme_linedraw()+theme(legend.position = "bottom")+labs(color="number of deaths by day 0")+guides(colour=guide_colorbar(title.position = "bottom"))+ylab(TeX("correlation ($\\rho$)"))+
  xlab("normalized days")


#day50<-subset(df_itter_Death,confirmed>=50&source=="Structural_proteins")
#day50<-get_day_specific_death_threshold_corr(50,7)
#cols_point<-wes_palette("BottleRocket2",n=4,type="discrete")
#cols<-wes_palette("Zissou1",type="continuous")

#day_50_plot<-ggplot(day50[which(day50$source=="Structural_proteins"),],aes(days,correlation,group=confirmed))+geom_line(aes(color=num_cnt))+
 # scale_color_gradientn(colours = cols)+geom_label(data = subset(day50,days%in%c(1,7,14,21)&source=="Structural_proteins"),aes(days,correlation),label=c(1,7,14,21))+theme_linedraw()+theme(legend.position = "bottom")+labs(color="number of countries")+guides(colour=guide_colorbar(title.position = "bottom"))+ggtitle("Days since 50th reported death")

#create plot of just day 50 figure 3B
death_thres<-50
day_point_plots<-lapply(1:4,function(w){
  day<-c(1,5,10,15)[w]
  day_point<-get_day_data(day,death_thres,6)
  dd<-day_point%>%arrange(Structural.proteins)%>%mutate(struct_prot_rank=1:length(Structural.proteins))%>%arrange(death_per_pop)%>%mutate(death_rank=1:length(death_per_pop))
   ggscatter(dd,x="death_rank",y="struct_prot_rank",label="country",add="reg.line",repel = T)+
     stat_cor(method = "spearman",label.x = max(dd$death_rank)*.5,label.y = max(dd$struct_prot_rank))+
     ggtitle(paste("day",day,sep=" "))+
    xlab("death rate rank")+ylab("EnsembleMHC population rank")+scale_x_continuous(breaks = seq(2,16,2))+scale_y_continuous(breaks = seq(2,16,2))
})  

corrs<-do.call(ggarrange,day_point_plots)

#create box plots figure 3C
day_point_boxes<-lapply(c(1,5,10,15),function(day){
  #get day specific data for the selected daya
  day_point<-get_day_data(day,death_thres,6)
  
  #create data for the box plots
  dd<-day_point%>%arrange(Structural.proteins)%>%mutate(struct_prot_rank=1:length(Structural.proteins))%>%
    arrange(death_per_pop)%>%mutate(death_rank=1:length(death_per_pop))%>%mutate(group=if_else(struct_prot_rank>median(struct_prot_rank),"upper half","lower half"))
  #calculate statisitcal significance of difference
   wt<-round(wilcox.test(dd$death_per_pop[which(dd$group=="upper half")],dd$death_per_pop[which(dd$group=="lower half")])[[3]],2)
  dd%>%ggplot(aes(x=group,y=death_per_pop,fill=group))+geom_boxplot()+annotate(geom = "text",x = 1.5,y=.8*max(dd$death_per_pop),label=paste("p =", wt))+
    xlab("EnsembleMHC population rank")+ylab("deaths per 1M")+theme_pubclean()+theme(legend.position = "none")+ggtitle(paste("day",day,sep=" "))
})  

boxes<-do.call(ggarrange,day_point_boxes)

figure_3BC<-ggarrange(corrs,boxes ,widths = c(.6,.4))




figure_3  <- ggarrange(death_plot+theme(legend.key.size = unit(5,"mm")),figure_3BC,ncol=1,heights = c(1.5,2))

ggsave(figure_3,filename = "~/Covid-19/plots/main_figures/Figure_3.pdf",width = 18,height = 12)

write.csv(norm_df_Death,"~/Covid-19/improtant_intermediate_datasets/iteration_of_day_death_threshold_1_100_norm.csv",row.names = F)
write.csv(df_Death,"~/Covid-19/improtant_intermediate_datasets/iteration_of_day_death_threshold_1_100.csv",row.names = F)
write.csv(all_data_clean,"~/Covid-19/improtant_intermediate_datasets/JHU_data_clean.csv",row.names = F)



#power analyis. filter by p
df_spearman_pwr<-norm_df_Death
df_spearman_pwr$pwr<-apply(df_spearman_pwr,1,function(w){
  
  pwr.r.test(r = as.numeric(w[3]),n = as.numeric(w[6]),sig.level = .05)[[4]]
  })

pwred_norm<-df_spearman_pwr[which(df_spearman_pwr$pwr>=.8),]

#df_spearman_pwr%>%ggplot(aes(x=pwr,y=p_value,color=source))+geom_density_2d()+theme_pubclean()
#df_spearman_pwr%>%ggplot(aes(x=confirmed,y=pwr,color=source))+geom_density_2d()+theme_pubclean()+geom_vline(xintercept = 100)

#data.frame(x=seq(10,1000,10),mean_pwr=sapply(seq(1,1000,10),function(n){
#mean(df_spearman_pwr$pwr[which(df_spearman_pwr$confirmed<=n)])
#}))%>%ggplot(aes(x,mean_pwr))+geom_line()+geom_vline(xintercept = 100,color="red")+
 # ylab("mean power")+xlab("Death threshold")+ggtitle("Statisitical power as a function of death threshold by day 0")+
  #theme_classic()


# data.frame(x=1:1000,mean_pwr=sapply(1:1000,function(n){
#   df_spearman_pwr$num_cnt[which(df_spearman_pwr$confirmed==n)][14]
# }))%>%ggplot(aes(x,mean_pwr))+geom_line()+geom_vline(xintercept = 100)+geom_hline(yintercept = 10.5)+ylab("number of countries at day 7")+
#   xlab("Death threshold")+
#   theme_classic()



R=1
df_spearman_pwr$PPV<-(df_spearman_pwr$pwr*R)/((df_spearman_pwr$pwr*R)+df_spearman_pwr$p_value)


PPV_with_respect_to_R<-do.call(rbind,lapply(seq(.1,1,.1),function(R){
  
  c(w=R,table(df_spearman_pwr$source[which((df_spearman_pwr$pwr*R)/((df_spearman_pwr$pwr*R)+df_spearman_pwr$p_value)>=.95)])/table(df_spearman_pwr$source))
}))%>%data.frame()%>%melt(id.vars="w")%>%ggplot(aes(x=w,y=value,color=variable))+geom_line()+xlab("R")+ylab("proportion of correlations with PPV > 95%")+theme_pubclean()+
  theme(legend.position = "none")

# data.frame(x=1:100,mean_corr=sapply(1:100,function(n){
#   #median(df_spearman_pwr$correlation[which(df_spearman_pwr$confirmed==n&df_spearman_pwr$source=="Structural Proteins")])
#    length(df_spearman_pwr$correlation[which(df_spearman_pwr$confirmed==n&df_spearman_pwr$source=="Structural Proteins"&df_spearman_pwr$PPV>=.95)])/
#      sum(df_spearman_pwr$confirmed==n)
#   # 
# }))%>%ggplot(aes(x=x,y=mean_corr))+geom_line()

ggsave(PPV_with_respect_to_R,filename = "~/Covid-19/plots/SI_figures/SI_PPV_with_respect_to_R.pdf")

list(
  mean_cor_all=mean(norm_df_Death$correlation[which(norm_df_Death$source=="entire SARS-CoV-2 proteome")]),
  mean_cor_struct=mean(norm_df_Death$correlation[which(norm_df_Death$source=="SARS-CoV-2 structural proteins")]),
  sig_overall<-table(norm_df_Death$sig)[2]/sum(table(norm_df_Death$sig)),
  stat_sig_all=table(norm_df_Death$sig[which(norm_df_Death$source=="entire SARS-CoV-2 proteome")])[2]/sum(table(norm_df_Death$sig[which(norm_df_Death$source=="entire SARS-CoV-2 proteome")])),
  stat_sig_struct=table(norm_df_Death$sig[which(norm_df_Death$source=="SARS-CoV-2 structural proteins")])[2]/sum(table(norm_df_Death$sig[which(norm_df_Death$source=="SARS-CoV-2 structural proteins")])),
  #total_with_power=table(df_spearman_pwr$pwr>=.8)[2]/sum(table(df_spearman_pwr$pwr>=.8)),
  #all_pwr=table(df_spearman_pwr$pwr[which(df_spearman_pwr$source=="entire SARS-CoV-2 proteome")]>=.8)[2]/sum(table(df_spearman_pwr$pwr[which(df_spearman_pwr$source=="entire SARS-CoV-2 proteome")]>=.8)),
  #struct_pwr=table(df_spearman_pwr$pwr[which(df_spearman_pwr$source=="SARS-CoV-2 structural proteins")]>=.8)[2]/sum(table(df_spearman_pwr$pwr[which(df_spearman_pwr$source=="SARS-CoV-2 structural proteins")]>=.8)),
  #mean_pwr_struct=mean(pwred_norm$correlation[which(pwred_norm$source=="SARS-CoV-2 structural proteins")]),
  PPV=table(df_spearman_pwr$source[which(df_spearman_pwr$PPV>=.95)])/table(df_spearman_pwr$source),
  mean_gt_half=mean(df_spearman_pwr$correlation[which(df_spearman_pwr$p_value<=.05&df_spearman_pwr$days>.5)]),
  mean_lt_half=mean(df_spearman_pwr$correlation[which(df_spearman_pwr$p_value<=.05&df_spearman_pwr$days<.5)]),
  sd_gt_half=sd(df_spearman_pwr$correlation[which(df_spearman_pwr$p_value<=.05&df_spearman_pwr$days>.5)]),
  sd_lt_half=sd(df_spearman_pwr$correlation[which(df_spearman_pwr$p_value<=.05&df_spearman_pwr$days<.5)])
  
)


