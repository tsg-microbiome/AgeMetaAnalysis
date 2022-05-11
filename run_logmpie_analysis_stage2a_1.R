load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\LogMPie\\logmpie_age_analysis.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\logmpie_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\LogMPie\\logmpie_Analysis_2021_Revision.RData")

rm(logmpie_combined_df_sum_stat_species)
rm(logmpie_combined_df_sum_stat_species_clr)
rm(logmpie_combined_df_sum_stat_genus)
rm(logmpie_combined_df_sum_stat_genus_clr)

library(robumeta)
library(metafor)
library(dplyr)
library(effsize)
library(MASS)
library(sfsmisc)
library(gplots)
library(RColorBrewer)
library(metap)

print("Computing Associations at Species level")

HighlyDetectedSpecies <- names(which(100*colSums(apply(logmpie_select_age_final_species,2,function(x)(ifelse(x!=0,1,0))))/nrow(logmpie_select_age_final_species)>=10))

logmpie_combined_df_sum_stat_species <- as.data.frame(cbind(df_logmpie_diversity_uniqueness,logmpie_select_age_final_species[rownames(df_logmpie_diversity_uniqueness),HighlyDetectedSpecies]))

logmpie_combined_df_sum_stat_species_clr <- as.data.frame(cbind(df_logmpie_diversity_uniqueness,logmpie_select_age_final_species_clr[rownames(df_logmpie_diversity_uniqueness),HighlyDetectedSpecies]))

logmpie_clr_relab_sample_correlations_species <- as.data.frame(matrix(NA,nrow(logmpie_combined_df_sum_stat_species),2))
rownames(logmpie_clr_relab_sample_correlations_species) <- rownames(logmpie_combined_df_sum_stat_species)
colnames(logmpie_clr_relab_sample_correlations_species) <- c("study_name","kendall_correlation")
for(i in 1:nrow(logmpie_combined_df_sum_stat_species))
{
	sample_name <- rownames(logmpie_combined_df_sum_stat_species)[i]
	study_name <- "LogMPie"
	df_temp <- data.frame(clr=as.numeric(logmpie_combined_df_sum_stat_species_clr[sample_name,HighlyDetectedSpecies]),relab=as.numeric(logmpie_combined_df_sum_stat_species[sample_name,HighlyDetectedSpecies]))
	logmpie_clr_relab_sample_correlations_species[sample_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	logmpie_clr_relab_sample_correlations_species[sample_name,2] <- study_name
}

logmpie_clr_relab_species_correlations <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),2))
rownames(logmpie_clr_relab_species_correlations) <- HighlyDetectedSpecies
colnames(logmpie_clr_relab_species_correlations) <- c("kendall_correlation","study_name")
for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	study_name <- "LogMPie"
	df_temp <- data.frame(clr=as.numeric(logmpie_combined_df_sum_stat_species_clr[,species_name]),relab=as.numeric(logmpie_combined_df_sum_stat_species[,species_name]))
	logmpie_clr_relab_species_correlations[species_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	logmpie_clr_relab_species_correlations[species_name,2] <- "LogMPie"
	
}

print("Computing Associations at Species level")

print("rlm bray uniqueness (clr)")
logmpie_rlm_detected_species_bray_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(logmpie_rlm_detected_species_bray_uniqueness_clr) <- HighlyDetectedSpecies
colnames(logmpie_rlm_detected_species_bray_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	logmpie_rlm_detected_species_bray_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=logmpie_combined_df_sum_stat_species_clr),var="species_bray_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_species_bray_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=logmpie_combined_df_sum_stat_species_clr),var="species_bray_uniqueness")$p.value
	
}

logmpie_rlm_detected_species_bray_uniqueness_clr$QValue <- p.adjust(logmpie_rlm_detected_species_bray_uniqueness_clr$PValue,method="fdr")

logmpie_rlm_detected_species_bray_uniqueness_clr$Direction <- ifelse(logmpie_rlm_detected_species_bray_uniqueness_clr$QValue<=0.05,3*sign(logmpie_rlm_detected_species_bray_uniqueness_clr$Estimate),ifelse(logmpie_rlm_detected_species_bray_uniqueness_clr$PValue<=0.05,2*sign(logmpie_rlm_detected_species_bray_uniqueness_clr$Estimate),sign(logmpie_rlm_detected_species_bray_uniqueness_clr$Estimate)))

print("rlm bray uniqueness (relab)")
logmpie_rlm_detected_species_bray_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(logmpie_rlm_detected_species_bray_uniqueness) <- HighlyDetectedSpecies
colnames(logmpie_rlm_detected_species_bray_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	logmpie_rlm_detected_species_bray_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=logmpie_combined_df_sum_stat_species),var="species_bray_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_species_bray_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=logmpie_combined_df_sum_stat_species),var="species_bray_uniqueness")$p.value
	
}

logmpie_rlm_detected_species_bray_uniqueness$QValue <- p.adjust(logmpie_rlm_detected_species_bray_uniqueness$PValue,method="fdr")

logmpie_rlm_detected_species_bray_uniqueness$Direction <- ifelse(logmpie_rlm_detected_species_bray_uniqueness$QValue<=0.05,3*sign(logmpie_rlm_detected_species_bray_uniqueness$Estimate),ifelse(logmpie_rlm_detected_species_bray_uniqueness$PValue<=0.05,2*sign(logmpie_rlm_detected_species_bray_uniqueness$Estimate),sign(logmpie_rlm_detected_species_bray_uniqueness$Estimate)))

print("rlm jaccard uniqueness (clr)")
logmpie_rlm_detected_species_jaccard_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(logmpie_rlm_detected_species_jaccard_uniqueness_clr) <- HighlyDetectedSpecies
colnames(logmpie_rlm_detected_species_jaccard_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	logmpie_rlm_detected_species_jaccard_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=logmpie_combined_df_sum_stat_species_clr),var="species_jaccard_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_species_jaccard_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=logmpie_combined_df_sum_stat_species_clr),var="species_jaccard_uniqueness")$p.value
	
}

logmpie_rlm_detected_species_jaccard_uniqueness_clr$QValue <- p.adjust(logmpie_rlm_detected_species_jaccard_uniqueness_clr$PValue,method="fdr")

logmpie_rlm_detected_species_jaccard_uniqueness_clr$Direction <- ifelse(logmpie_rlm_detected_species_jaccard_uniqueness_clr$QValue<=0.05,3*sign(logmpie_rlm_detected_species_jaccard_uniqueness_clr$Estimate),ifelse(logmpie_rlm_detected_species_jaccard_uniqueness_clr$PValue<=0.05,2*sign(logmpie_rlm_detected_species_jaccard_uniqueness_clr$Estimate),sign(logmpie_rlm_detected_species_jaccard_uniqueness_clr$Estimate)))

print("rlm jaccard uniqueness (relab)")
logmpie_rlm_detected_species_jaccard_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(logmpie_rlm_detected_species_jaccard_uniqueness) <- HighlyDetectedSpecies
colnames(logmpie_rlm_detected_species_jaccard_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	logmpie_rlm_detected_species_jaccard_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=logmpie_combined_df_sum_stat_species),var="species_jaccard_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_species_jaccard_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=logmpie_combined_df_sum_stat_species),var="species_jaccard_uniqueness")$p.value
	
}

logmpie_rlm_detected_species_jaccard_uniqueness$QValue <- p.adjust(logmpie_rlm_detected_species_jaccard_uniqueness$PValue,method="fdr")

logmpie_rlm_detected_species_jaccard_uniqueness$Direction <- ifelse(logmpie_rlm_detected_species_jaccard_uniqueness$QValue<=0.05,3*sign(logmpie_rlm_detected_species_jaccard_uniqueness$Estimate),ifelse(logmpie_rlm_detected_species_jaccard_uniqueness$PValue<=0.05,2*sign(logmpie_rlm_detected_species_jaccard_uniqueness$Estimate),sign(logmpie_rlm_detected_species_jaccard_uniqueness$Estimate)))

print("rlm aitchison uniqueness (clr)")
logmpie_rlm_detected_species_aitchison_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(logmpie_rlm_detected_species_aitchison_uniqueness_clr) <- HighlyDetectedSpecies
colnames(logmpie_rlm_detected_species_aitchison_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	logmpie_rlm_detected_species_aitchison_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=logmpie_combined_df_sum_stat_species_clr),var="species_aitchison_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_species_aitchison_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=logmpie_combined_df_sum_stat_species_clr),var="species_aitchison_uniqueness")$p.value
	
}

logmpie_rlm_detected_species_aitchison_uniqueness_clr$QValue <- p.adjust(logmpie_rlm_detected_species_aitchison_uniqueness_clr$PValue,method="fdr")

logmpie_rlm_detected_species_aitchison_uniqueness_clr$Direction <- ifelse(logmpie_rlm_detected_species_aitchison_uniqueness_clr$QValue<=0.05,3*sign(logmpie_rlm_detected_species_aitchison_uniqueness_clr$Estimate),ifelse(logmpie_rlm_detected_species_aitchison_uniqueness_clr$PValue<=0.05,2*sign(logmpie_rlm_detected_species_aitchison_uniqueness_clr$Estimate),sign(logmpie_rlm_detected_species_aitchison_uniqueness_clr$Estimate)))

print("rlm aitchison uniqueness (relab)")
logmpie_rlm_detected_species_aitchison_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(logmpie_rlm_detected_species_aitchison_uniqueness) <- HighlyDetectedSpecies
colnames(logmpie_rlm_detected_species_aitchison_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	logmpie_rlm_detected_species_aitchison_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=logmpie_combined_df_sum_stat_species),var="species_aitchison_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_species_aitchison_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=logmpie_combined_df_sum_stat_species),var="species_aitchison_uniqueness")$p.value
	
}

logmpie_rlm_detected_species_aitchison_uniqueness$QValue <- p.adjust(logmpie_rlm_detected_species_aitchison_uniqueness$PValue,method="fdr")

logmpie_rlm_detected_species_aitchison_uniqueness$Direction <- ifelse(logmpie_rlm_detected_species_aitchison_uniqueness$QValue<=0.05,3*sign(logmpie_rlm_detected_species_aitchison_uniqueness$Estimate),ifelse(logmpie_rlm_detected_species_aitchison_uniqueness$PValue<=0.05,2*sign(logmpie_rlm_detected_species_aitchison_uniqueness$Estimate),sign(logmpie_rlm_detected_species_aitchison_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
logmpie_rlm_detected_species_kendall_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(logmpie_rlm_detected_species_kendall_uniqueness_clr) <- HighlyDetectedSpecies
colnames(logmpie_rlm_detected_species_kendall_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	logmpie_rlm_detected_species_kendall_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=logmpie_combined_df_sum_stat_species_clr),var="species_kendall_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_species_kendall_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=logmpie_combined_df_sum_stat_species_clr),var="species_kendall_uniqueness")$p.value
	
}

logmpie_rlm_detected_species_kendall_uniqueness_clr$QValue <- p.adjust(logmpie_rlm_detected_species_kendall_uniqueness_clr$PValue,method="fdr")

logmpie_rlm_detected_species_kendall_uniqueness_clr$Direction <- ifelse(logmpie_rlm_detected_species_kendall_uniqueness_clr$QValue<=0.05,3*sign(logmpie_rlm_detected_species_kendall_uniqueness_clr$Estimate),ifelse(logmpie_rlm_detected_species_kendall_uniqueness_clr$PValue<=0.05,2*sign(logmpie_rlm_detected_species_kendall_uniqueness_clr$Estimate),sign(logmpie_rlm_detected_species_kendall_uniqueness_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
logmpie_rlm_detected_species_kendall_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(logmpie_rlm_detected_species_kendall_uniqueness) <- HighlyDetectedSpecies
colnames(logmpie_rlm_detected_species_kendall_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	logmpie_rlm_detected_species_kendall_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=logmpie_combined_df_sum_stat_species),var="species_kendall_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_species_kendall_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=logmpie_combined_df_sum_stat_species),var="species_kendall_uniqueness")$p.value
	
}

logmpie_rlm_detected_species_kendall_uniqueness$QValue <- p.adjust(logmpie_rlm_detected_species_kendall_uniqueness$PValue,method="fdr")

logmpie_rlm_detected_species_kendall_uniqueness$Direction <- ifelse(logmpie_rlm_detected_species_kendall_uniqueness$QValue<=0.05,3*sign(logmpie_rlm_detected_species_kendall_uniqueness$Estimate),ifelse(logmpie_rlm_detected_species_kendall_uniqueness$PValue<=0.05,2*sign(logmpie_rlm_detected_species_kendall_uniqueness$Estimate),sign(logmpie_rlm_detected_species_kendall_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
logmpie_rlm_detected_species_shannon_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(logmpie_rlm_detected_species_shannon_clr) <- HighlyDetectedSpecies
colnames(logmpie_rlm_detected_species_shannon_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	logmpie_rlm_detected_species_shannon_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_shannon")),data=logmpie_combined_df_sum_stat_species_clr),var="species_shannon")$coefficients[2,1]
	logmpie_rlm_detected_species_shannon_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_shannon")),data=logmpie_combined_df_sum_stat_species_clr),var="species_shannon")$p.value
	
}

logmpie_rlm_detected_species_shannon_clr$QValue <- p.adjust(logmpie_rlm_detected_species_shannon_clr$PValue,method="fdr")

logmpie_rlm_detected_species_shannon_clr$Direction <- ifelse(logmpie_rlm_detected_species_shannon_clr$QValue<=0.05,3*sign(logmpie_rlm_detected_species_shannon_clr$Estimate),ifelse(logmpie_rlm_detected_species_shannon_clr$PValue<=0.05,2*sign(logmpie_rlm_detected_species_shannon_clr$Estimate),sign(logmpie_rlm_detected_species_shannon_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
logmpie_rlm_detected_species_shannon <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(logmpie_rlm_detected_species_shannon) <- HighlyDetectedSpecies
colnames(logmpie_rlm_detected_species_shannon) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	logmpie_rlm_detected_species_shannon[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_shannon")),data=logmpie_combined_df_sum_stat_species),var="species_shannon")$coefficients[2,1]
	logmpie_rlm_detected_species_shannon[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_shannon")),data=logmpie_combined_df_sum_stat_species),var="species_shannon")$p.value
	
}

logmpie_rlm_detected_species_shannon$QValue <- p.adjust(logmpie_rlm_detected_species_shannon$PValue,method="fdr")

logmpie_rlm_detected_species_shannon$Direction <- ifelse(logmpie_rlm_detected_species_shannon$QValue<=0.05,3*sign(logmpie_rlm_detected_species_shannon$Estimate),ifelse(logmpie_rlm_detected_species_shannon$PValue<=0.05,2*sign(logmpie_rlm_detected_species_shannon$Estimate),sign(logmpie_rlm_detected_species_shannon$Estimate)))

print("Computing Associations at Genus level")

HighlyDetectedGenus <- names(which(100*colSums(apply(logmpie_select_age_final_genus,2,function(x)(ifelse(x!=0,1,0))))/nrow(logmpie_select_age_final_genus)>=10))

logmpie_combined_df_sum_stat_genus <- as.data.frame(cbind(df_logmpie_diversity_uniqueness,logmpie_select_age_final_genus[rownames(df_logmpie_diversity_uniqueness),HighlyDetectedGenus]))

logmpie_combined_df_sum_stat_genus_clr <- as.data.frame(cbind(df_logmpie_diversity_uniqueness,logmpie_select_age_final_genus_clr[rownames(df_logmpie_diversity_uniqueness),HighlyDetectedGenus]))

colnames(logmpie_combined_df_sum_stat_genus) <- sub("-","_",colnames(logmpie_combined_df_sum_stat_genus))
colnames(logmpie_combined_df_sum_stat_genus_clr) <- sub("-","_",colnames(logmpie_combined_df_sum_stat_genus_clr))
HighlyDetectedGenus <- sub("-","_",HighlyDetectedGenus)

logmpie_clr_relab_sample_correlations_genus <- as.data.frame(matrix(NA,nrow(logmpie_combined_df_sum_stat_genus),2))
rownames(logmpie_clr_relab_sample_correlations_genus) <- rownames(logmpie_combined_df_sum_stat_genus)
colnames(logmpie_clr_relab_sample_correlations_genus) <- c("study_name","kendall_correlation")
for(i in 1:nrow(logmpie_combined_df_sum_stat_genus))
{
	sample_name <- rownames(logmpie_combined_df_sum_stat_genus)[i]
	study_name <- "LogMPie"
	df_temp <- data.frame(clr=as.numeric(logmpie_combined_df_sum_stat_genus_clr[sample_name,HighlyDetectedGenus]),relab=as.numeric(logmpie_combined_df_sum_stat_genus[sample_name,HighlyDetectedGenus]))
	logmpie_clr_relab_sample_correlations_genus[sample_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	logmpie_clr_relab_sample_correlations_genus[sample_name,2] <- study_name
}

logmpie_clr_relab_genus_correlations <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),2))
rownames(logmpie_clr_relab_genus_correlations) <- HighlyDetectedGenus
colnames(logmpie_clr_relab_genus_correlations) <- c("kendall_correlation","study_name")
for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	study_name <- "LogMPie"
	df_temp <- data.frame(clr=as.numeric(logmpie_combined_df_sum_stat_genus_clr[,genus_name]),relab=as.numeric(logmpie_combined_df_sum_stat_genus[,genus_name]))
	logmpie_clr_relab_genus_correlations[genus_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	logmpie_clr_relab_genus_correlations[genus_name,2] <- "LogMPie"
	
}


print("Computing Associations at Genus level")

print("rlm bray uniqueness (clr)")
logmpie_rlm_detected_genus_bray_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(logmpie_rlm_detected_genus_bray_uniqueness_clr) <- HighlyDetectedGenus
colnames(logmpie_rlm_detected_genus_bray_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	logmpie_rlm_detected_genus_bray_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=logmpie_combined_df_sum_stat_genus_clr),var="genus_bray_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_genus_bray_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=logmpie_combined_df_sum_stat_genus_clr),var="genus_bray_uniqueness")$p.value
	
}

logmpie_rlm_detected_genus_bray_uniqueness_clr$QValue <- p.adjust(logmpie_rlm_detected_genus_bray_uniqueness_clr$PValue,method="fdr")

logmpie_rlm_detected_genus_bray_uniqueness_clr$Direction <- ifelse(logmpie_rlm_detected_genus_bray_uniqueness_clr$QValue<=0.05,3*sign(logmpie_rlm_detected_genus_bray_uniqueness_clr$Estimate),ifelse(logmpie_rlm_detected_genus_bray_uniqueness_clr$PValue<=0.05,2*sign(logmpie_rlm_detected_genus_bray_uniqueness_clr$Estimate),sign(logmpie_rlm_detected_genus_bray_uniqueness_clr$Estimate)))

print("rlm bray uniqueness (relab)")
logmpie_rlm_detected_genus_bray_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(logmpie_rlm_detected_genus_bray_uniqueness) <- HighlyDetectedGenus
colnames(logmpie_rlm_detected_genus_bray_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	logmpie_rlm_detected_genus_bray_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=logmpie_combined_df_sum_stat_genus),var="genus_bray_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_genus_bray_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=logmpie_combined_df_sum_stat_genus),var="genus_bray_uniqueness")$p.value
	
}

logmpie_rlm_detected_genus_bray_uniqueness$QValue <- p.adjust(logmpie_rlm_detected_genus_bray_uniqueness$PValue,method="fdr")

logmpie_rlm_detected_genus_bray_uniqueness$Direction <- ifelse(logmpie_rlm_detected_genus_bray_uniqueness$QValue<=0.05,3*sign(logmpie_rlm_detected_genus_bray_uniqueness$Estimate),ifelse(logmpie_rlm_detected_genus_bray_uniqueness$PValue<=0.05,2*sign(logmpie_rlm_detected_genus_bray_uniqueness$Estimate),sign(logmpie_rlm_detected_genus_bray_uniqueness$Estimate)))

print("rlm jaccard uniqueness (clr)")
logmpie_rlm_detected_genus_jaccard_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(logmpie_rlm_detected_genus_jaccard_uniqueness_clr) <- HighlyDetectedGenus
colnames(logmpie_rlm_detected_genus_jaccard_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	logmpie_rlm_detected_genus_jaccard_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=logmpie_combined_df_sum_stat_genus_clr),var="genus_jaccard_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_genus_jaccard_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=logmpie_combined_df_sum_stat_genus_clr),var="genus_jaccard_uniqueness")$p.value
	
}

logmpie_rlm_detected_genus_jaccard_uniqueness_clr$QValue <- p.adjust(logmpie_rlm_detected_genus_jaccard_uniqueness_clr$PValue,method="fdr")

logmpie_rlm_detected_genus_jaccard_uniqueness_clr$Direction <- ifelse(logmpie_rlm_detected_genus_jaccard_uniqueness_clr$QValue<=0.05,3*sign(logmpie_rlm_detected_genus_jaccard_uniqueness_clr$Estimate),ifelse(logmpie_rlm_detected_genus_jaccard_uniqueness_clr$PValue<=0.05,2*sign(logmpie_rlm_detected_genus_jaccard_uniqueness_clr$Estimate),sign(logmpie_rlm_detected_genus_jaccard_uniqueness_clr$Estimate)))

print("rlm jaccard uniqueness (relab)")
logmpie_rlm_detected_genus_jaccard_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(logmpie_rlm_detected_genus_jaccard_uniqueness) <- HighlyDetectedGenus
colnames(logmpie_rlm_detected_genus_jaccard_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	logmpie_rlm_detected_genus_jaccard_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=logmpie_combined_df_sum_stat_genus),var="genus_jaccard_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_genus_jaccard_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=logmpie_combined_df_sum_stat_genus),var="genus_jaccard_uniqueness")$p.value
	
}

logmpie_rlm_detected_genus_jaccard_uniqueness$QValue <- p.adjust(logmpie_rlm_detected_genus_jaccard_uniqueness$PValue,method="fdr")

logmpie_rlm_detected_genus_jaccard_uniqueness$Direction <- ifelse(logmpie_rlm_detected_genus_jaccard_uniqueness$QValue<=0.05,3*sign(logmpie_rlm_detected_genus_jaccard_uniqueness$Estimate),ifelse(logmpie_rlm_detected_genus_jaccard_uniqueness$PValue<=0.05,2*sign(logmpie_rlm_detected_genus_jaccard_uniqueness$Estimate),sign(logmpie_rlm_detected_genus_jaccard_uniqueness$Estimate)))

print("rlm aitchison uniqueness (clr)")
logmpie_rlm_detected_genus_aitchison_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(logmpie_rlm_detected_genus_aitchison_uniqueness_clr) <- HighlyDetectedGenus
colnames(logmpie_rlm_detected_genus_aitchison_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	logmpie_rlm_detected_genus_aitchison_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=logmpie_combined_df_sum_stat_genus_clr),var="genus_aitchison_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_genus_aitchison_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=logmpie_combined_df_sum_stat_genus_clr),var="genus_aitchison_uniqueness")$p.value
	
}

logmpie_rlm_detected_genus_aitchison_uniqueness_clr$QValue <- p.adjust(logmpie_rlm_detected_genus_aitchison_uniqueness_clr$PValue,method="fdr")

logmpie_rlm_detected_genus_aitchison_uniqueness_clr$Direction <- ifelse(logmpie_rlm_detected_genus_aitchison_uniqueness_clr$QValue<=0.05,3*sign(logmpie_rlm_detected_genus_aitchison_uniqueness_clr$Estimate),ifelse(logmpie_rlm_detected_genus_aitchison_uniqueness_clr$PValue<=0.05,2*sign(logmpie_rlm_detected_genus_aitchison_uniqueness_clr$Estimate),sign(logmpie_rlm_detected_genus_aitchison_uniqueness_clr$Estimate)))

print("rlm aitchison uniqueness (relab)")
logmpie_rlm_detected_genus_aitchison_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(logmpie_rlm_detected_genus_aitchison_uniqueness) <- HighlyDetectedGenus
colnames(logmpie_rlm_detected_genus_aitchison_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	logmpie_rlm_detected_genus_aitchison_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=logmpie_combined_df_sum_stat_genus),var="genus_aitchison_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_genus_aitchison_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=logmpie_combined_df_sum_stat_genus),var="genus_aitchison_uniqueness")$p.value
	
}

logmpie_rlm_detected_genus_aitchison_uniqueness$QValue <- p.adjust(logmpie_rlm_detected_genus_aitchison_uniqueness$PValue,method="fdr")

logmpie_rlm_detected_genus_aitchison_uniqueness$Direction <- ifelse(logmpie_rlm_detected_genus_aitchison_uniqueness$QValue<=0.05,3*sign(logmpie_rlm_detected_genus_aitchison_uniqueness$Estimate),ifelse(logmpie_rlm_detected_genus_aitchison_uniqueness$PValue<=0.05,2*sign(logmpie_rlm_detected_genus_aitchison_uniqueness$Estimate),sign(logmpie_rlm_detected_genus_aitchison_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
logmpie_rlm_detected_genus_kendall_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(logmpie_rlm_detected_genus_kendall_uniqueness_clr) <- HighlyDetectedGenus
colnames(logmpie_rlm_detected_genus_kendall_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	logmpie_rlm_detected_genus_kendall_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=logmpie_combined_df_sum_stat_genus_clr),var="genus_kendall_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_genus_kendall_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=logmpie_combined_df_sum_stat_genus_clr),var="genus_kendall_uniqueness")$p.value
	
}

logmpie_rlm_detected_genus_kendall_uniqueness_clr$QValue <- p.adjust(logmpie_rlm_detected_genus_kendall_uniqueness_clr$PValue,method="fdr")

logmpie_rlm_detected_genus_kendall_uniqueness_clr$Direction <- ifelse(logmpie_rlm_detected_genus_kendall_uniqueness_clr$QValue<=0.05,3*sign(logmpie_rlm_detected_genus_kendall_uniqueness_clr$Estimate),ifelse(logmpie_rlm_detected_genus_kendall_uniqueness_clr$PValue<=0.05,2*sign(logmpie_rlm_detected_genus_kendall_uniqueness_clr$Estimate),sign(logmpie_rlm_detected_genus_kendall_uniqueness_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
logmpie_rlm_detected_genus_kendall_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(logmpie_rlm_detected_genus_kendall_uniqueness) <- HighlyDetectedGenus
colnames(logmpie_rlm_detected_genus_kendall_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	logmpie_rlm_detected_genus_kendall_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=logmpie_combined_df_sum_stat_genus),var="genus_kendall_uniqueness")$coefficients[2,1]
	logmpie_rlm_detected_genus_kendall_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=logmpie_combined_df_sum_stat_genus),var="genus_kendall_uniqueness")$p.value
	
}

logmpie_rlm_detected_genus_kendall_uniqueness$QValue <- p.adjust(logmpie_rlm_detected_genus_kendall_uniqueness$PValue,method="fdr")

logmpie_rlm_detected_genus_kendall_uniqueness$Direction <- ifelse(logmpie_rlm_detected_genus_kendall_uniqueness$QValue<=0.05,3*sign(logmpie_rlm_detected_genus_kendall_uniqueness$Estimate),ifelse(logmpie_rlm_detected_genus_kendall_uniqueness$PValue<=0.05,2*sign(logmpie_rlm_detected_genus_kendall_uniqueness$Estimate),sign(logmpie_rlm_detected_genus_kendall_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
logmpie_rlm_detected_genus_shannon_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(logmpie_rlm_detected_genus_shannon_clr) <- HighlyDetectedGenus
colnames(logmpie_rlm_detected_genus_shannon_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	logmpie_rlm_detected_genus_shannon_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=logmpie_combined_df_sum_stat_genus_clr),var="genus_shannon")$coefficients[2,1]
	logmpie_rlm_detected_genus_shannon_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=logmpie_combined_df_sum_stat_genus_clr),var="genus_shannon")$p.value
	
}

logmpie_rlm_detected_genus_shannon_clr$QValue <- p.adjust(logmpie_rlm_detected_genus_shannon_clr$PValue,method="fdr")

logmpie_rlm_detected_genus_shannon_clr$Direction <- ifelse(logmpie_rlm_detected_genus_shannon_clr$QValue<=0.05,3*sign(logmpie_rlm_detected_genus_shannon_clr$Estimate),ifelse(logmpie_rlm_detected_genus_shannon_clr$PValue<=0.05,2*sign(logmpie_rlm_detected_genus_shannon_clr$Estimate),sign(logmpie_rlm_detected_genus_shannon_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
logmpie_rlm_detected_genus_shannon <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(logmpie_rlm_detected_genus_shannon) <- HighlyDetectedGenus
colnames(logmpie_rlm_detected_genus_shannon) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	logmpie_rlm_detected_genus_shannon[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=logmpie_combined_df_sum_stat_genus),var="genus_shannon")$coefficients[2,1]
	logmpie_rlm_detected_genus_shannon[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=logmpie_combined_df_sum_stat_genus),var="genus_shannon")$p.value
	
}

logmpie_rlm_detected_genus_shannon$QValue <- p.adjust(logmpie_rlm_detected_genus_shannon$PValue,method="fdr")

logmpie_rlm_detected_genus_shannon$Direction <- ifelse(logmpie_rlm_detected_genus_shannon$QValue<=0.05,3*sign(logmpie_rlm_detected_genus_shannon$Estimate),ifelse(logmpie_rlm_detected_genus_shannon$PValue<=0.05,2*sign(logmpie_rlm_detected_genus_shannon$Estimate),sign(logmpie_rlm_detected_genus_shannon$Estimate)))

logmpie_HighlyDetectedSpecies <- HighlyDetectedSpecies
logmpie_HighlyDetectedGenus <- HighlyDetectedGenus

save(logmpie_combined_df_sum_stat_species,logmpie_rlm_detected_species_shannon,logmpie_rlm_detected_species_kendall_uniqueness,logmpie_rlm_detected_species_aitchison_uniqueness,logmpie_rlm_detected_species_jaccard_uniqueness,logmpie_rlm_detected_species_bray_uniqueness,logmpie_combined_df_sum_stat_species_clr,logmpie_rlm_detected_species_shannon_clr,logmpie_rlm_detected_species_kendall_uniqueness_clr,logmpie_rlm_detected_species_aitchison_uniqueness_clr,logmpie_rlm_detected_species_jaccard_uniqueness_clr,logmpie_rlm_detected_species_bray_uniqueness_clr,logmpie_combined_df_sum_stat_genus,logmpie_rlm_detected_genus_shannon,logmpie_rlm_detected_genus_kendall_uniqueness,logmpie_rlm_detected_genus_aitchison_uniqueness,logmpie_rlm_detected_genus_jaccard_uniqueness,logmpie_rlm_detected_genus_bray_uniqueness,logmpie_combined_df_sum_stat_genus_clr,logmpie_rlm_detected_genus_shannon_clr,logmpie_rlm_detected_genus_kendall_uniqueness_clr,logmpie_rlm_detected_genus_aitchison_uniqueness_clr,logmpie_rlm_detected_genus_jaccard_uniqueness_clr,logmpie_rlm_detected_genus_bray_uniqueness_clr,logmpie_HighlyDetectedSpecies,logmpie_HighlyDetectedGenus,logmpie_clr_relab_sample_correlations_species,logmpie_clr_relab_sample_correlations_genus,logmpie_clr_relab_species_correlations,logmpie_clr_relab_genus_correlations,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\logmpie_stage2a_results.RData")

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\LogMPie\\logmpie_analysis_2021_Revision.RData")

rm(list=ls())
