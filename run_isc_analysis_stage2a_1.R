load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\isc_age_analysis.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\isc_Analysis_2021_Revision.RData")

rm(isc_combined_df_sum_stat_species)
rm(isc_combined_df_sum_stat_species_clr)
rm(isc_combined_df_sum_stat_genus)
rm(isc_combined_df_sum_stat_genus_clr)

df_isc_diversity_uniqueness <- as.data.frame(df_isc_diversity_uniqueness[,c("species_shannon","genus_shannon","age","study_name","study_condition","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")])

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

HighlyDetectedSpecies <- names(which(100*colSums(apply(isc_select_age_final_species,2,function(x)(ifelse(x!=0,1,0))))/nrow(isc_select_age_final_species)>=10))

isc_combined_df_sum_stat_species <- as.data.frame(cbind(df_isc_diversity_uniqueness,isc_select_age_final_species[rownames(df_isc_diversity_uniqueness),HighlyDetectedSpecies]))

isc_combined_df_sum_stat_species_clr <- as.data.frame(cbind(df_isc_diversity_uniqueness,isc_select_age_final_species_clr[rownames(df_isc_diversity_uniqueness),HighlyDetectedSpecies]))

isc_clr_relab_sample_correlations_species <- as.data.frame(matrix(NA,nrow(isc_combined_df_sum_stat_species),2))
rownames(isc_clr_relab_sample_correlations_species) <- rownames(isc_combined_df_sum_stat_species)
colnames(isc_clr_relab_sample_correlations_species) <- c("study_name","kendall_correlation")
for(i in 1:nrow(isc_combined_df_sum_stat_species))
{
	sample_name <- rownames(isc_combined_df_sum_stat_species)[i]
	study_name <- "ISC"
	df_temp <- data.frame(clr=as.numeric(isc_combined_df_sum_stat_species_clr[sample_name,HighlyDetectedSpecies]),relab=as.numeric(isc_combined_df_sum_stat_species[sample_name,HighlyDetectedSpecies]))
	isc_clr_relab_sample_correlations_species[sample_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	isc_clr_relab_sample_correlations_species[sample_name,2] <- study_name
}

isc_clr_relab_species_correlations <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),2))
rownames(isc_clr_relab_species_correlations) <- HighlyDetectedSpecies
colnames(isc_clr_relab_species_correlations) <- c("kendall_correlation","study_name")
for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	study_name <- "ISC"
	df_temp <- data.frame(clr=as.numeric(isc_combined_df_sum_stat_species_clr[,species_name]),relab=as.numeric(isc_combined_df_sum_stat_species[,species_name]))
	isc_clr_relab_species_correlations[species_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	isc_clr_relab_species_correlations[species_name,2] <- "ISC"
	
}

print("Computing Associations at Species level")

print("rlm bray uniqueness (clr)")
isc_rlm_detected_species_bray_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(isc_rlm_detected_species_bray_uniqueness_clr) <- HighlyDetectedSpecies
colnames(isc_rlm_detected_species_bray_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	isc_rlm_detected_species_bray_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=isc_combined_df_sum_stat_species_clr),var="species_bray_uniqueness")$coefficients[2,1]
	isc_rlm_detected_species_bray_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=isc_combined_df_sum_stat_species_clr),var="species_bray_uniqueness")$p.value
	
}

isc_rlm_detected_species_bray_uniqueness_clr$QValue <- p.adjust(isc_rlm_detected_species_bray_uniqueness_clr$PValue,method="fdr")

isc_rlm_detected_species_bray_uniqueness_clr$Direction <- ifelse(isc_rlm_detected_species_bray_uniqueness_clr$QValue<=0.05,3*sign(isc_rlm_detected_species_bray_uniqueness_clr$Estimate),ifelse(isc_rlm_detected_species_bray_uniqueness_clr$PValue<=0.05,2*sign(isc_rlm_detected_species_bray_uniqueness_clr$Estimate),sign(isc_rlm_detected_species_bray_uniqueness_clr$Estimate)))

print("rlm bray uniqueness (relab)")
isc_rlm_detected_species_bray_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(isc_rlm_detected_species_bray_uniqueness) <- HighlyDetectedSpecies
colnames(isc_rlm_detected_species_bray_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	isc_rlm_detected_species_bray_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=isc_combined_df_sum_stat_species),var="species_bray_uniqueness")$coefficients[2,1]
	isc_rlm_detected_species_bray_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=isc_combined_df_sum_stat_species),var="species_bray_uniqueness")$p.value
	
}

isc_rlm_detected_species_bray_uniqueness$QValue <- p.adjust(isc_rlm_detected_species_bray_uniqueness$PValue,method="fdr")

isc_rlm_detected_species_bray_uniqueness$Direction <- ifelse(isc_rlm_detected_species_bray_uniqueness$QValue<=0.05,3*sign(isc_rlm_detected_species_bray_uniqueness$Estimate),ifelse(isc_rlm_detected_species_bray_uniqueness$PValue<=0.05,2*sign(isc_rlm_detected_species_bray_uniqueness$Estimate),sign(isc_rlm_detected_species_bray_uniqueness$Estimate)))

print("rlm jaccard uniqueness (clr)")
isc_rlm_detected_species_jaccard_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(isc_rlm_detected_species_jaccard_uniqueness_clr) <- HighlyDetectedSpecies
colnames(isc_rlm_detected_species_jaccard_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	isc_rlm_detected_species_jaccard_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=isc_combined_df_sum_stat_species_clr),var="species_jaccard_uniqueness")$coefficients[2,1]
	isc_rlm_detected_species_jaccard_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=isc_combined_df_sum_stat_species_clr),var="species_jaccard_uniqueness")$p.value
	
}

isc_rlm_detected_species_jaccard_uniqueness_clr$QValue <- p.adjust(isc_rlm_detected_species_jaccard_uniqueness_clr$PValue,method="fdr")

isc_rlm_detected_species_jaccard_uniqueness_clr$Direction <- ifelse(isc_rlm_detected_species_jaccard_uniqueness_clr$QValue<=0.05,3*sign(isc_rlm_detected_species_jaccard_uniqueness_clr$Estimate),ifelse(isc_rlm_detected_species_jaccard_uniqueness_clr$PValue<=0.05,2*sign(isc_rlm_detected_species_jaccard_uniqueness_clr$Estimate),sign(isc_rlm_detected_species_jaccard_uniqueness_clr$Estimate)))

print("rlm jaccard uniqueness (relab)")
isc_rlm_detected_species_jaccard_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(isc_rlm_detected_species_jaccard_uniqueness) <- HighlyDetectedSpecies
colnames(isc_rlm_detected_species_jaccard_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	isc_rlm_detected_species_jaccard_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=isc_combined_df_sum_stat_species),var="species_jaccard_uniqueness")$coefficients[2,1]
	isc_rlm_detected_species_jaccard_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=isc_combined_df_sum_stat_species),var="species_jaccard_uniqueness")$p.value
	
}

isc_rlm_detected_species_jaccard_uniqueness$QValue <- p.adjust(isc_rlm_detected_species_jaccard_uniqueness$PValue,method="fdr")

isc_rlm_detected_species_jaccard_uniqueness$Direction <- ifelse(isc_rlm_detected_species_jaccard_uniqueness$QValue<=0.05,3*sign(isc_rlm_detected_species_jaccard_uniqueness$Estimate),ifelse(isc_rlm_detected_species_jaccard_uniqueness$PValue<=0.05,2*sign(isc_rlm_detected_species_jaccard_uniqueness$Estimate),sign(isc_rlm_detected_species_jaccard_uniqueness$Estimate)))

print("rlm aitchison uniqueness (clr)")
isc_rlm_detected_species_aitchison_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(isc_rlm_detected_species_aitchison_uniqueness_clr) <- HighlyDetectedSpecies
colnames(isc_rlm_detected_species_aitchison_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	isc_rlm_detected_species_aitchison_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=isc_combined_df_sum_stat_species_clr),var="species_aitchison_uniqueness")$coefficients[2,1]
	isc_rlm_detected_species_aitchison_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=isc_combined_df_sum_stat_species_clr),var="species_aitchison_uniqueness")$p.value
	
}

isc_rlm_detected_species_aitchison_uniqueness_clr$QValue <- p.adjust(isc_rlm_detected_species_aitchison_uniqueness_clr$PValue,method="fdr")

isc_rlm_detected_species_aitchison_uniqueness_clr$Direction <- ifelse(isc_rlm_detected_species_aitchison_uniqueness_clr$QValue<=0.05,3*sign(isc_rlm_detected_species_aitchison_uniqueness_clr$Estimate),ifelse(isc_rlm_detected_species_aitchison_uniqueness_clr$PValue<=0.05,2*sign(isc_rlm_detected_species_aitchison_uniqueness_clr$Estimate),sign(isc_rlm_detected_species_aitchison_uniqueness_clr$Estimate)))

print("rlm aitchison uniqueness (relab)")
isc_rlm_detected_species_aitchison_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(isc_rlm_detected_species_aitchison_uniqueness) <- HighlyDetectedSpecies
colnames(isc_rlm_detected_species_aitchison_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	isc_rlm_detected_species_aitchison_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=isc_combined_df_sum_stat_species),var="species_aitchison_uniqueness")$coefficients[2,1]
	isc_rlm_detected_species_aitchison_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=isc_combined_df_sum_stat_species),var="species_aitchison_uniqueness")$p.value
	
}

isc_rlm_detected_species_aitchison_uniqueness$QValue <- p.adjust(isc_rlm_detected_species_aitchison_uniqueness$PValue,method="fdr")

isc_rlm_detected_species_aitchison_uniqueness$Direction <- ifelse(isc_rlm_detected_species_aitchison_uniqueness$QValue<=0.05,3*sign(isc_rlm_detected_species_aitchison_uniqueness$Estimate),ifelse(isc_rlm_detected_species_aitchison_uniqueness$PValue<=0.05,2*sign(isc_rlm_detected_species_aitchison_uniqueness$Estimate),sign(isc_rlm_detected_species_aitchison_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
isc_rlm_detected_species_kendall_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(isc_rlm_detected_species_kendall_uniqueness_clr) <- HighlyDetectedSpecies
colnames(isc_rlm_detected_species_kendall_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	isc_rlm_detected_species_kendall_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=isc_combined_df_sum_stat_species_clr),var="species_kendall_uniqueness")$coefficients[2,1]
	isc_rlm_detected_species_kendall_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=isc_combined_df_sum_stat_species_clr),var="species_kendall_uniqueness")$p.value
	
}

isc_rlm_detected_species_kendall_uniqueness_clr$QValue <- p.adjust(isc_rlm_detected_species_kendall_uniqueness_clr$PValue,method="fdr")

isc_rlm_detected_species_kendall_uniqueness_clr$Direction <- ifelse(isc_rlm_detected_species_kendall_uniqueness_clr$QValue<=0.05,3*sign(isc_rlm_detected_species_kendall_uniqueness_clr$Estimate),ifelse(isc_rlm_detected_species_kendall_uniqueness_clr$PValue<=0.05,2*sign(isc_rlm_detected_species_kendall_uniqueness_clr$Estimate),sign(isc_rlm_detected_species_kendall_uniqueness_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
isc_rlm_detected_species_kendall_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(isc_rlm_detected_species_kendall_uniqueness) <- HighlyDetectedSpecies
colnames(isc_rlm_detected_species_kendall_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	isc_rlm_detected_species_kendall_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=isc_combined_df_sum_stat_species),var="species_kendall_uniqueness")$coefficients[2,1]
	isc_rlm_detected_species_kendall_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=isc_combined_df_sum_stat_species),var="species_kendall_uniqueness")$p.value
	
}

isc_rlm_detected_species_kendall_uniqueness$QValue <- p.adjust(isc_rlm_detected_species_kendall_uniqueness$PValue,method="fdr")

isc_rlm_detected_species_kendall_uniqueness$Direction <- ifelse(isc_rlm_detected_species_kendall_uniqueness$QValue<=0.05,3*sign(isc_rlm_detected_species_kendall_uniqueness$Estimate),ifelse(isc_rlm_detected_species_kendall_uniqueness$PValue<=0.05,2*sign(isc_rlm_detected_species_kendall_uniqueness$Estimate),sign(isc_rlm_detected_species_kendall_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
isc_rlm_detected_species_shannon_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(isc_rlm_detected_species_shannon_clr) <- HighlyDetectedSpecies
colnames(isc_rlm_detected_species_shannon_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	isc_rlm_detected_species_shannon_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_shannon")),data=isc_combined_df_sum_stat_species_clr),var="species_shannon")$coefficients[2,1]
	isc_rlm_detected_species_shannon_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_shannon")),data=isc_combined_df_sum_stat_species_clr),var="species_shannon")$p.value
	
}

isc_rlm_detected_species_shannon_clr$QValue <- p.adjust(isc_rlm_detected_species_shannon_clr$PValue,method="fdr")

isc_rlm_detected_species_shannon_clr$Direction <- ifelse(isc_rlm_detected_species_shannon_clr$QValue<=0.05,3*sign(isc_rlm_detected_species_shannon_clr$Estimate),ifelse(isc_rlm_detected_species_shannon_clr$PValue<=0.05,2*sign(isc_rlm_detected_species_shannon_clr$Estimate),sign(isc_rlm_detected_species_shannon_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
isc_rlm_detected_species_shannon <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(isc_rlm_detected_species_shannon) <- HighlyDetectedSpecies
colnames(isc_rlm_detected_species_shannon) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	isc_rlm_detected_species_shannon[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_shannon")),data=isc_combined_df_sum_stat_species),var="species_shannon")$coefficients[2,1]
	isc_rlm_detected_species_shannon[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_shannon")),data=isc_combined_df_sum_stat_species),var="species_shannon")$p.value
	
}

isc_rlm_detected_species_shannon$QValue <- p.adjust(isc_rlm_detected_species_shannon$PValue,method="fdr")

isc_rlm_detected_species_shannon$Direction <- ifelse(isc_rlm_detected_species_shannon$QValue<=0.05,3*sign(isc_rlm_detected_species_shannon$Estimate),ifelse(isc_rlm_detected_species_shannon$PValue<=0.05,2*sign(isc_rlm_detected_species_shannon$Estimate),sign(isc_rlm_detected_species_shannon$Estimate)))

print("Computing Associations at Genus level")

HighlyDetectedGenus <- names(which(100*colSums(apply(isc_select_age_final_genus,2,function(x)(ifelse(x!=0,1,0))))/nrow(isc_select_age_final_genus)>=10))

isc_combined_df_sum_stat_genus <- as.data.frame(cbind(df_isc_diversity_uniqueness,isc_select_age_final_genus[rownames(df_isc_diversity_uniqueness),HighlyDetectedGenus]))

isc_combined_df_sum_stat_genus_clr <- as.data.frame(cbind(df_isc_diversity_uniqueness,isc_select_age_final_genus_clr[rownames(df_isc_diversity_uniqueness),HighlyDetectedGenus]))

colnames(isc_combined_df_sum_stat_genus) <- sub("-","_",colnames(isc_combined_df_sum_stat_genus))
colnames(isc_combined_df_sum_stat_genus_clr) <- sub("-","_",colnames(isc_combined_df_sum_stat_genus_clr))
HighlyDetectedGenus <- sub("-","_",HighlyDetectedGenus)

isc_clr_relab_sample_correlations_genus <- as.data.frame(matrix(NA,nrow(isc_combined_df_sum_stat_genus),2))
rownames(isc_clr_relab_sample_correlations_genus) <- rownames(isc_combined_df_sum_stat_genus)
colnames(isc_clr_relab_sample_correlations_genus) <- c("study_name","kendall_correlation")
for(i in 1:nrow(isc_combined_df_sum_stat_genus))
{
	sample_name <- rownames(isc_combined_df_sum_stat_genus)[i]
	study_name <- "ISC"
	df_temp <- data.frame(clr=as.numeric(isc_combined_df_sum_stat_genus_clr[sample_name,HighlyDetectedGenus]),relab=as.numeric(isc_combined_df_sum_stat_genus[sample_name,HighlyDetectedGenus]))
	isc_clr_relab_sample_correlations_genus[sample_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	isc_clr_relab_sample_correlations_genus[sample_name,2] <- study_name
}

isc_clr_relab_genus_correlations <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),2))
rownames(isc_clr_relab_genus_correlations) <- HighlyDetectedGenus
colnames(isc_clr_relab_genus_correlations) <- c("kendall_correlation","study_name")
for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	study_name <- "ISC"
	df_temp <- data.frame(clr=as.numeric(isc_combined_df_sum_stat_genus_clr[,genus_name]),relab=as.numeric(isc_combined_df_sum_stat_genus[,genus_name]))
	isc_clr_relab_genus_correlations[genus_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	isc_clr_relab_genus_correlations[genus_name,2] <- "ISC"
	
}


print("Computing Associations at Genus level")

print("rlm bray uniqueness (clr)")
isc_rlm_detected_genus_bray_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(isc_rlm_detected_genus_bray_uniqueness_clr) <- HighlyDetectedGenus
colnames(isc_rlm_detected_genus_bray_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	isc_rlm_detected_genus_bray_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=isc_combined_df_sum_stat_genus_clr),var="genus_bray_uniqueness")$coefficients[2,1]
	isc_rlm_detected_genus_bray_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=isc_combined_df_sum_stat_genus_clr),var="genus_bray_uniqueness")$p.value
	
}

isc_rlm_detected_genus_bray_uniqueness_clr$QValue <- p.adjust(isc_rlm_detected_genus_bray_uniqueness_clr$PValue,method="fdr")

isc_rlm_detected_genus_bray_uniqueness_clr$Direction <- ifelse(isc_rlm_detected_genus_bray_uniqueness_clr$QValue<=0.05,3*sign(isc_rlm_detected_genus_bray_uniqueness_clr$Estimate),ifelse(isc_rlm_detected_genus_bray_uniqueness_clr$PValue<=0.05,2*sign(isc_rlm_detected_genus_bray_uniqueness_clr$Estimate),sign(isc_rlm_detected_genus_bray_uniqueness_clr$Estimate)))

print("rlm bray uniqueness (relab)")
isc_rlm_detected_genus_bray_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(isc_rlm_detected_genus_bray_uniqueness) <- HighlyDetectedGenus
colnames(isc_rlm_detected_genus_bray_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	isc_rlm_detected_genus_bray_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=isc_combined_df_sum_stat_genus),var="genus_bray_uniqueness")$coefficients[2,1]
	isc_rlm_detected_genus_bray_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=isc_combined_df_sum_stat_genus),var="genus_bray_uniqueness")$p.value
	
}

isc_rlm_detected_genus_bray_uniqueness$QValue <- p.adjust(isc_rlm_detected_genus_bray_uniqueness$PValue,method="fdr")

isc_rlm_detected_genus_bray_uniqueness$Direction <- ifelse(isc_rlm_detected_genus_bray_uniqueness$QValue<=0.05,3*sign(isc_rlm_detected_genus_bray_uniqueness$Estimate),ifelse(isc_rlm_detected_genus_bray_uniqueness$PValue<=0.05,2*sign(isc_rlm_detected_genus_bray_uniqueness$Estimate),sign(isc_rlm_detected_genus_bray_uniqueness$Estimate)))

print("rlm jaccard uniqueness (clr)")
isc_rlm_detected_genus_jaccard_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(isc_rlm_detected_genus_jaccard_uniqueness_clr) <- HighlyDetectedGenus
colnames(isc_rlm_detected_genus_jaccard_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	isc_rlm_detected_genus_jaccard_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=isc_combined_df_sum_stat_genus_clr),var="genus_jaccard_uniqueness")$coefficients[2,1]
	isc_rlm_detected_genus_jaccard_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=isc_combined_df_sum_stat_genus_clr),var="genus_jaccard_uniqueness")$p.value
	
}

isc_rlm_detected_genus_jaccard_uniqueness_clr$QValue <- p.adjust(isc_rlm_detected_genus_jaccard_uniqueness_clr$PValue,method="fdr")

isc_rlm_detected_genus_jaccard_uniqueness_clr$Direction <- ifelse(isc_rlm_detected_genus_jaccard_uniqueness_clr$QValue<=0.05,3*sign(isc_rlm_detected_genus_jaccard_uniqueness_clr$Estimate),ifelse(isc_rlm_detected_genus_jaccard_uniqueness_clr$PValue<=0.05,2*sign(isc_rlm_detected_genus_jaccard_uniqueness_clr$Estimate),sign(isc_rlm_detected_genus_jaccard_uniqueness_clr$Estimate)))

print("rlm jaccard uniqueness (relab)")
isc_rlm_detected_genus_jaccard_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(isc_rlm_detected_genus_jaccard_uniqueness) <- HighlyDetectedGenus
colnames(isc_rlm_detected_genus_jaccard_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	isc_rlm_detected_genus_jaccard_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=isc_combined_df_sum_stat_genus),var="genus_jaccard_uniqueness")$coefficients[2,1]
	isc_rlm_detected_genus_jaccard_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=isc_combined_df_sum_stat_genus),var="genus_jaccard_uniqueness")$p.value
	
}

isc_rlm_detected_genus_jaccard_uniqueness$QValue <- p.adjust(isc_rlm_detected_genus_jaccard_uniqueness$PValue,method="fdr")

isc_rlm_detected_genus_jaccard_uniqueness$Direction <- ifelse(isc_rlm_detected_genus_jaccard_uniqueness$QValue<=0.05,3*sign(isc_rlm_detected_genus_jaccard_uniqueness$Estimate),ifelse(isc_rlm_detected_genus_jaccard_uniqueness$PValue<=0.05,2*sign(isc_rlm_detected_genus_jaccard_uniqueness$Estimate),sign(isc_rlm_detected_genus_jaccard_uniqueness$Estimate)))

print("rlm aitchison uniqueness (clr)")
isc_rlm_detected_genus_aitchison_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(isc_rlm_detected_genus_aitchison_uniqueness_clr) <- HighlyDetectedGenus
colnames(isc_rlm_detected_genus_aitchison_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	isc_rlm_detected_genus_aitchison_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=isc_combined_df_sum_stat_genus_clr),var="genus_aitchison_uniqueness")$coefficients[2,1]
	isc_rlm_detected_genus_aitchison_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=isc_combined_df_sum_stat_genus_clr),var="genus_aitchison_uniqueness")$p.value
	
}

isc_rlm_detected_genus_aitchison_uniqueness_clr$QValue <- p.adjust(isc_rlm_detected_genus_aitchison_uniqueness_clr$PValue,method="fdr")

isc_rlm_detected_genus_aitchison_uniqueness_clr$Direction <- ifelse(isc_rlm_detected_genus_aitchison_uniqueness_clr$QValue<=0.05,3*sign(isc_rlm_detected_genus_aitchison_uniqueness_clr$Estimate),ifelse(isc_rlm_detected_genus_aitchison_uniqueness_clr$PValue<=0.05,2*sign(isc_rlm_detected_genus_aitchison_uniqueness_clr$Estimate),sign(isc_rlm_detected_genus_aitchison_uniqueness_clr$Estimate)))

print("rlm aitchison uniqueness (relab)")
isc_rlm_detected_genus_aitchison_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(isc_rlm_detected_genus_aitchison_uniqueness) <- HighlyDetectedGenus
colnames(isc_rlm_detected_genus_aitchison_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	isc_rlm_detected_genus_aitchison_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=isc_combined_df_sum_stat_genus),var="genus_aitchison_uniqueness")$coefficients[2,1]
	isc_rlm_detected_genus_aitchison_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=isc_combined_df_sum_stat_genus),var="genus_aitchison_uniqueness")$p.value
	
}

isc_rlm_detected_genus_aitchison_uniqueness$QValue <- p.adjust(isc_rlm_detected_genus_aitchison_uniqueness$PValue,method="fdr")

isc_rlm_detected_genus_aitchison_uniqueness$Direction <- ifelse(isc_rlm_detected_genus_aitchison_uniqueness$QValue<=0.05,3*sign(isc_rlm_detected_genus_aitchison_uniqueness$Estimate),ifelse(isc_rlm_detected_genus_aitchison_uniqueness$PValue<=0.05,2*sign(isc_rlm_detected_genus_aitchison_uniqueness$Estimate),sign(isc_rlm_detected_genus_aitchison_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
isc_rlm_detected_genus_kendall_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(isc_rlm_detected_genus_kendall_uniqueness_clr) <- HighlyDetectedGenus
colnames(isc_rlm_detected_genus_kendall_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	isc_rlm_detected_genus_kendall_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=isc_combined_df_sum_stat_genus_clr),var="genus_kendall_uniqueness")$coefficients[2,1]
	isc_rlm_detected_genus_kendall_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=isc_combined_df_sum_stat_genus_clr),var="genus_kendall_uniqueness")$p.value
	
}

isc_rlm_detected_genus_kendall_uniqueness_clr$QValue <- p.adjust(isc_rlm_detected_genus_kendall_uniqueness_clr$PValue,method="fdr")

isc_rlm_detected_genus_kendall_uniqueness_clr$Direction <- ifelse(isc_rlm_detected_genus_kendall_uniqueness_clr$QValue<=0.05,3*sign(isc_rlm_detected_genus_kendall_uniqueness_clr$Estimate),ifelse(isc_rlm_detected_genus_kendall_uniqueness_clr$PValue<=0.05,2*sign(isc_rlm_detected_genus_kendall_uniqueness_clr$Estimate),sign(isc_rlm_detected_genus_kendall_uniqueness_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
isc_rlm_detected_genus_kendall_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(isc_rlm_detected_genus_kendall_uniqueness) <- HighlyDetectedGenus
colnames(isc_rlm_detected_genus_kendall_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	isc_rlm_detected_genus_kendall_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=isc_combined_df_sum_stat_genus),var="genus_kendall_uniqueness")$coefficients[2,1]
	isc_rlm_detected_genus_kendall_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=isc_combined_df_sum_stat_genus),var="genus_kendall_uniqueness")$p.value
	
}

isc_rlm_detected_genus_kendall_uniqueness$QValue <- p.adjust(isc_rlm_detected_genus_kendall_uniqueness$PValue,method="fdr")

isc_rlm_detected_genus_kendall_uniqueness$Direction <- ifelse(isc_rlm_detected_genus_kendall_uniqueness$QValue<=0.05,3*sign(isc_rlm_detected_genus_kendall_uniqueness$Estimate),ifelse(isc_rlm_detected_genus_kendall_uniqueness$PValue<=0.05,2*sign(isc_rlm_detected_genus_kendall_uniqueness$Estimate),sign(isc_rlm_detected_genus_kendall_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
isc_rlm_detected_genus_shannon_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(isc_rlm_detected_genus_shannon_clr) <- HighlyDetectedGenus
colnames(isc_rlm_detected_genus_shannon_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	isc_rlm_detected_genus_shannon_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=isc_combined_df_sum_stat_genus_clr),var="genus_shannon")$coefficients[2,1]
	isc_rlm_detected_genus_shannon_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=isc_combined_df_sum_stat_genus_clr),var="genus_shannon")$p.value
	
}

isc_rlm_detected_genus_shannon_clr$QValue <- p.adjust(isc_rlm_detected_genus_shannon_clr$PValue,method="fdr")

isc_rlm_detected_genus_shannon_clr$Direction <- ifelse(isc_rlm_detected_genus_shannon_clr$QValue<=0.05,3*sign(isc_rlm_detected_genus_shannon_clr$Estimate),ifelse(isc_rlm_detected_genus_shannon_clr$PValue<=0.05,2*sign(isc_rlm_detected_genus_shannon_clr$Estimate),sign(isc_rlm_detected_genus_shannon_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
isc_rlm_detected_genus_shannon <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(isc_rlm_detected_genus_shannon) <- HighlyDetectedGenus
colnames(isc_rlm_detected_genus_shannon) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	isc_rlm_detected_genus_shannon[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=isc_combined_df_sum_stat_genus),var="genus_shannon")$coefficients[2,1]
	isc_rlm_detected_genus_shannon[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=isc_combined_df_sum_stat_genus),var="genus_shannon")$p.value
	
}

isc_rlm_detected_genus_shannon$QValue <- p.adjust(isc_rlm_detected_genus_shannon$PValue,method="fdr")

isc_rlm_detected_genus_shannon$Direction <- ifelse(isc_rlm_detected_genus_shannon$QValue<=0.05,3*sign(isc_rlm_detected_genus_shannon$Estimate),ifelse(isc_rlm_detected_genus_shannon$PValue<=0.05,2*sign(isc_rlm_detected_genus_shannon$Estimate),sign(isc_rlm_detected_genus_shannon$Estimate)))

isc_HighlyDetectedSpecies <- HighlyDetectedSpecies
isc_HighlyDetectedGenus <- HighlyDetectedGenus

save(isc_combined_df_sum_stat_species,isc_rlm_detected_species_shannon,isc_rlm_detected_species_kendall_uniqueness,isc_rlm_detected_species_aitchison_uniqueness,isc_rlm_detected_species_jaccard_uniqueness,isc_rlm_detected_species_bray_uniqueness,isc_combined_df_sum_stat_species_clr,isc_rlm_detected_species_shannon_clr,isc_rlm_detected_species_kendall_uniqueness_clr,isc_rlm_detected_species_aitchison_uniqueness_clr,isc_rlm_detected_species_jaccard_uniqueness_clr,isc_rlm_detected_species_bray_uniqueness_clr,isc_combined_df_sum_stat_genus,isc_rlm_detected_genus_shannon,isc_rlm_detected_genus_kendall_uniqueness,isc_rlm_detected_genus_aitchison_uniqueness,isc_rlm_detected_genus_jaccard_uniqueness,isc_rlm_detected_genus_bray_uniqueness,isc_combined_df_sum_stat_genus_clr,isc_rlm_detected_genus_shannon_clr,isc_rlm_detected_genus_kendall_uniqueness_clr,isc_rlm_detected_genus_aitchison_uniqueness_clr,isc_rlm_detected_genus_jaccard_uniqueness_clr,isc_rlm_detected_genus_bray_uniqueness_clr,isc_HighlyDetectedSpecies,isc_HighlyDetectedGenus,isc_clr_relab_sample_correlations_species,isc_clr_relab_sample_correlations_genus,isc_clr_relab_species_correlations,isc_clr_relab_genus_correlations,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_stage2a_results.RData")

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\isc_analysis_2021_Revision.RData")

#rm(list=ls())
