load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\NUAGE\\nuage_age_analysis.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\nuage_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\NUAGE\\nuage_Analysis_2021_Revision.RData")

rm(nuage_combined_df_sum_stat_species)
rm(nuage_combined_df_sum_stat_species_clr)
rm(nuage_combined_df_sum_stat_genus)
rm(nuage_combined_df_sum_stat_genus_clr)

df_nuage_diversity_uniqueness <- as.data.frame(df_nuage_diversity_uniqueness[,c("species_shannon","genus_shannon","study_name","age","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")])

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

HighlyDetectedSpecies <- names(which(100*colSums(apply(nuage_select_age_final_species,2,function(x)(ifelse(x!=0,1,0))))/nrow(nuage_select_age_final_species)>=10))

nuage_combined_df_sum_stat_species <- as.data.frame(cbind(df_nuage_diversity_uniqueness,nuage_select_age_final_species[rownames(df_nuage_diversity_uniqueness),HighlyDetectedSpecies]))

nuage_combined_df_sum_stat_species_clr <- as.data.frame(cbind(df_nuage_diversity_uniqueness,nuage_select_age_final_species_clr[rownames(df_nuage_diversity_uniqueness),HighlyDetectedSpecies]))

nuage_clr_relab_sample_correlations_species <- as.data.frame(matrix(NA,nrow(nuage_combined_df_sum_stat_species),2))
rownames(nuage_clr_relab_sample_correlations_species) <- rownames(nuage_combined_df_sum_stat_species)
colnames(nuage_clr_relab_sample_correlations_species) <- c("study_name","kendall_correlation")
for(i in 1:nrow(nuage_combined_df_sum_stat_species))
{
	sample_name <- rownames(nuage_combined_df_sum_stat_species)[i]
	study_name <- "NUAGE"
	df_temp <- data.frame(clr=as.numeric(nuage_combined_df_sum_stat_species_clr[sample_name,HighlyDetectedSpecies]),relab=as.numeric(nuage_combined_df_sum_stat_species[sample_name,HighlyDetectedSpecies]))
	nuage_clr_relab_sample_correlations_species[sample_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	nuage_clr_relab_sample_correlations_species[sample_name,2] <- study_name
}

nuage_clr_relab_species_correlations <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),2))
rownames(nuage_clr_relab_species_correlations) <- HighlyDetectedSpecies
colnames(nuage_clr_relab_species_correlations) <- c("kendall_correlation","study_name")
for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	study_name <- "NUAGE"
	df_temp <- data.frame(clr=as.numeric(nuage_combined_df_sum_stat_species_clr[,species_name]),relab=as.numeric(nuage_combined_df_sum_stat_species[,species_name]))
	nuage_clr_relab_species_correlations[species_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	nuage_clr_relab_species_correlations[species_name,2] <- "NUAGE"
	
}

print("Computing Associations at Species level")

print("rlm bray uniqueness (clr)")
nuage_rlm_detected_species_bray_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(nuage_rlm_detected_species_bray_uniqueness_clr) <- HighlyDetectedSpecies
colnames(nuage_rlm_detected_species_bray_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	nuage_rlm_detected_species_bray_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=nuage_combined_df_sum_stat_species_clr),var="species_bray_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_species_bray_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=nuage_combined_df_sum_stat_species_clr),var="species_bray_uniqueness")$p.value
	
}

nuage_rlm_detected_species_bray_uniqueness_clr$QValue <- p.adjust(nuage_rlm_detected_species_bray_uniqueness_clr$PValue,method="fdr")

nuage_rlm_detected_species_bray_uniqueness_clr$Direction <- ifelse(nuage_rlm_detected_species_bray_uniqueness_clr$QValue<=0.05,3*sign(nuage_rlm_detected_species_bray_uniqueness_clr$Estimate),ifelse(nuage_rlm_detected_species_bray_uniqueness_clr$PValue<=0.05,2*sign(nuage_rlm_detected_species_bray_uniqueness_clr$Estimate),sign(nuage_rlm_detected_species_bray_uniqueness_clr$Estimate)))

print("rlm bray uniqueness (relab)")
nuage_rlm_detected_species_bray_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(nuage_rlm_detected_species_bray_uniqueness) <- HighlyDetectedSpecies
colnames(nuage_rlm_detected_species_bray_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	nuage_rlm_detected_species_bray_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=nuage_combined_df_sum_stat_species),var="species_bray_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_species_bray_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=nuage_combined_df_sum_stat_species),var="species_bray_uniqueness")$p.value
	
}

nuage_rlm_detected_species_bray_uniqueness$QValue <- p.adjust(nuage_rlm_detected_species_bray_uniqueness$PValue,method="fdr")

nuage_rlm_detected_species_bray_uniqueness$Direction <- ifelse(nuage_rlm_detected_species_bray_uniqueness$QValue<=0.05,3*sign(nuage_rlm_detected_species_bray_uniqueness$Estimate),ifelse(nuage_rlm_detected_species_bray_uniqueness$PValue<=0.05,2*sign(nuage_rlm_detected_species_bray_uniqueness$Estimate),sign(nuage_rlm_detected_species_bray_uniqueness$Estimate)))

print("rlm jaccard uniqueness (clr)")
nuage_rlm_detected_species_jaccard_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(nuage_rlm_detected_species_jaccard_uniqueness_clr) <- HighlyDetectedSpecies
colnames(nuage_rlm_detected_species_jaccard_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	nuage_rlm_detected_species_jaccard_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=nuage_combined_df_sum_stat_species_clr),var="species_jaccard_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_species_jaccard_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=nuage_combined_df_sum_stat_species_clr),var="species_jaccard_uniqueness")$p.value
	
}

nuage_rlm_detected_species_jaccard_uniqueness_clr$QValue <- p.adjust(nuage_rlm_detected_species_jaccard_uniqueness_clr$PValue,method="fdr")

nuage_rlm_detected_species_jaccard_uniqueness_clr$Direction <- ifelse(nuage_rlm_detected_species_jaccard_uniqueness_clr$QValue<=0.05,3*sign(nuage_rlm_detected_species_jaccard_uniqueness_clr$Estimate),ifelse(nuage_rlm_detected_species_jaccard_uniqueness_clr$PValue<=0.05,2*sign(nuage_rlm_detected_species_jaccard_uniqueness_clr$Estimate),sign(nuage_rlm_detected_species_jaccard_uniqueness_clr$Estimate)))

print("rlm jaccard uniqueness (relab)")
nuage_rlm_detected_species_jaccard_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(nuage_rlm_detected_species_jaccard_uniqueness) <- HighlyDetectedSpecies
colnames(nuage_rlm_detected_species_jaccard_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	nuage_rlm_detected_species_jaccard_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=nuage_combined_df_sum_stat_species),var="species_jaccard_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_species_jaccard_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=nuage_combined_df_sum_stat_species),var="species_jaccard_uniqueness")$p.value
	
}

nuage_rlm_detected_species_jaccard_uniqueness$QValue <- p.adjust(nuage_rlm_detected_species_jaccard_uniqueness$PValue,method="fdr")

nuage_rlm_detected_species_jaccard_uniqueness$Direction <- ifelse(nuage_rlm_detected_species_jaccard_uniqueness$QValue<=0.05,3*sign(nuage_rlm_detected_species_jaccard_uniqueness$Estimate),ifelse(nuage_rlm_detected_species_jaccard_uniqueness$PValue<=0.05,2*sign(nuage_rlm_detected_species_jaccard_uniqueness$Estimate),sign(nuage_rlm_detected_species_jaccard_uniqueness$Estimate)))

print("rlm aitchison uniqueness (clr)")
nuage_rlm_detected_species_aitchison_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(nuage_rlm_detected_species_aitchison_uniqueness_clr) <- HighlyDetectedSpecies
colnames(nuage_rlm_detected_species_aitchison_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	nuage_rlm_detected_species_aitchison_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=nuage_combined_df_sum_stat_species_clr),var="species_aitchison_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_species_aitchison_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=nuage_combined_df_sum_stat_species_clr),var="species_aitchison_uniqueness")$p.value
	
}

nuage_rlm_detected_species_aitchison_uniqueness_clr$QValue <- p.adjust(nuage_rlm_detected_species_aitchison_uniqueness_clr$PValue,method="fdr")

nuage_rlm_detected_species_aitchison_uniqueness_clr$Direction <- ifelse(nuage_rlm_detected_species_aitchison_uniqueness_clr$QValue<=0.05,3*sign(nuage_rlm_detected_species_aitchison_uniqueness_clr$Estimate),ifelse(nuage_rlm_detected_species_aitchison_uniqueness_clr$PValue<=0.05,2*sign(nuage_rlm_detected_species_aitchison_uniqueness_clr$Estimate),sign(nuage_rlm_detected_species_aitchison_uniqueness_clr$Estimate)))

print("rlm aitchison uniqueness (relab)")
nuage_rlm_detected_species_aitchison_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(nuage_rlm_detected_species_aitchison_uniqueness) <- HighlyDetectedSpecies
colnames(nuage_rlm_detected_species_aitchison_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	nuage_rlm_detected_species_aitchison_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=nuage_combined_df_sum_stat_species),var="species_aitchison_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_species_aitchison_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=nuage_combined_df_sum_stat_species),var="species_aitchison_uniqueness")$p.value
	
}

nuage_rlm_detected_species_aitchison_uniqueness$QValue <- p.adjust(nuage_rlm_detected_species_aitchison_uniqueness$PValue,method="fdr")

nuage_rlm_detected_species_aitchison_uniqueness$Direction <- ifelse(nuage_rlm_detected_species_aitchison_uniqueness$QValue<=0.05,3*sign(nuage_rlm_detected_species_aitchison_uniqueness$Estimate),ifelse(nuage_rlm_detected_species_aitchison_uniqueness$PValue<=0.05,2*sign(nuage_rlm_detected_species_aitchison_uniqueness$Estimate),sign(nuage_rlm_detected_species_aitchison_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
nuage_rlm_detected_species_kendall_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(nuage_rlm_detected_species_kendall_uniqueness_clr) <- HighlyDetectedSpecies
colnames(nuage_rlm_detected_species_kendall_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	nuage_rlm_detected_species_kendall_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=nuage_combined_df_sum_stat_species_clr),var="species_kendall_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_species_kendall_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=nuage_combined_df_sum_stat_species_clr),var="species_kendall_uniqueness")$p.value
	
}

nuage_rlm_detected_species_kendall_uniqueness_clr$QValue <- p.adjust(nuage_rlm_detected_species_kendall_uniqueness_clr$PValue,method="fdr")

nuage_rlm_detected_species_kendall_uniqueness_clr$Direction <- ifelse(nuage_rlm_detected_species_kendall_uniqueness_clr$QValue<=0.05,3*sign(nuage_rlm_detected_species_kendall_uniqueness_clr$Estimate),ifelse(nuage_rlm_detected_species_kendall_uniqueness_clr$PValue<=0.05,2*sign(nuage_rlm_detected_species_kendall_uniqueness_clr$Estimate),sign(nuage_rlm_detected_species_kendall_uniqueness_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
nuage_rlm_detected_species_kendall_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(nuage_rlm_detected_species_kendall_uniqueness) <- HighlyDetectedSpecies
colnames(nuage_rlm_detected_species_kendall_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	nuage_rlm_detected_species_kendall_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=nuage_combined_df_sum_stat_species),var="species_kendall_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_species_kendall_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=nuage_combined_df_sum_stat_species),var="species_kendall_uniqueness")$p.value
	
}

nuage_rlm_detected_species_kendall_uniqueness$QValue <- p.adjust(nuage_rlm_detected_species_kendall_uniqueness$PValue,method="fdr")

nuage_rlm_detected_species_kendall_uniqueness$Direction <- ifelse(nuage_rlm_detected_species_kendall_uniqueness$QValue<=0.05,3*sign(nuage_rlm_detected_species_kendall_uniqueness$Estimate),ifelse(nuage_rlm_detected_species_kendall_uniqueness$PValue<=0.05,2*sign(nuage_rlm_detected_species_kendall_uniqueness$Estimate),sign(nuage_rlm_detected_species_kendall_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
nuage_rlm_detected_species_shannon_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(nuage_rlm_detected_species_shannon_clr) <- HighlyDetectedSpecies
colnames(nuage_rlm_detected_species_shannon_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	nuage_rlm_detected_species_shannon_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_shannon")),data=nuage_combined_df_sum_stat_species_clr),var="species_shannon")$coefficients[2,1]
	nuage_rlm_detected_species_shannon_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_shannon")),data=nuage_combined_df_sum_stat_species_clr),var="species_shannon")$p.value
	
}

nuage_rlm_detected_species_shannon_clr$QValue <- p.adjust(nuage_rlm_detected_species_shannon_clr$PValue,method="fdr")

nuage_rlm_detected_species_shannon_clr$Direction <- ifelse(nuage_rlm_detected_species_shannon_clr$QValue<=0.05,3*sign(nuage_rlm_detected_species_shannon_clr$Estimate),ifelse(nuage_rlm_detected_species_shannon_clr$PValue<=0.05,2*sign(nuage_rlm_detected_species_shannon_clr$Estimate),sign(nuage_rlm_detected_species_shannon_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
nuage_rlm_detected_species_shannon <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(nuage_rlm_detected_species_shannon) <- HighlyDetectedSpecies
colnames(nuage_rlm_detected_species_shannon) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	nuage_rlm_detected_species_shannon[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_shannon")),data=nuage_combined_df_sum_stat_species),var="species_shannon")$coefficients[2,1]
	nuage_rlm_detected_species_shannon[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_shannon")),data=nuage_combined_df_sum_stat_species),var="species_shannon")$p.value
	
}

nuage_rlm_detected_species_shannon$QValue <- p.adjust(nuage_rlm_detected_species_shannon$PValue,method="fdr")

nuage_rlm_detected_species_shannon$Direction <- ifelse(nuage_rlm_detected_species_shannon$QValue<=0.05,3*sign(nuage_rlm_detected_species_shannon$Estimate),ifelse(nuage_rlm_detected_species_shannon$PValue<=0.05,2*sign(nuage_rlm_detected_species_shannon$Estimate),sign(nuage_rlm_detected_species_shannon$Estimate)))

print("Computing Associations at Genus level")

HighlyDetectedGenus <- names(which(100*colSums(apply(nuage_select_age_final_genus,2,function(x)(ifelse(x!=0,1,0))))/nrow(nuage_select_age_final_genus)>=10))

nuage_combined_df_sum_stat_genus <- as.data.frame(cbind(df_nuage_diversity_uniqueness,nuage_select_age_final_genus[rownames(df_nuage_diversity_uniqueness),HighlyDetectedGenus]))

nuage_combined_df_sum_stat_genus_clr <- as.data.frame(cbind(df_nuage_diversity_uniqueness,nuage_select_age_final_genus_clr[rownames(df_nuage_diversity_uniqueness),HighlyDetectedGenus]))

colnames(nuage_combined_df_sum_stat_genus) <- sub("-","_",colnames(nuage_combined_df_sum_stat_genus))
colnames(nuage_combined_df_sum_stat_genus_clr) <- sub("-","_",colnames(nuage_combined_df_sum_stat_genus_clr))
HighlyDetectedGenus <- sub("-","_",HighlyDetectedGenus)

nuage_clr_relab_sample_correlations_genus <- as.data.frame(matrix(NA,nrow(nuage_combined_df_sum_stat_genus),2))
rownames(nuage_clr_relab_sample_correlations_genus) <- rownames(nuage_combined_df_sum_stat_genus)
colnames(nuage_clr_relab_sample_correlations_genus) <- c("study_name","kendall_correlation")
for(i in 1:nrow(nuage_combined_df_sum_stat_genus))
{
	sample_name <- rownames(nuage_combined_df_sum_stat_genus)[i]
	study_name <- "NUAGE"
	df_temp <- data.frame(clr=as.numeric(nuage_combined_df_sum_stat_genus_clr[sample_name,HighlyDetectedGenus]),relab=as.numeric(nuage_combined_df_sum_stat_genus[sample_name,HighlyDetectedGenus]))
	nuage_clr_relab_sample_correlations_genus[sample_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	nuage_clr_relab_sample_correlations_genus[sample_name,2] <- study_name
}

nuage_clr_relab_genus_correlations <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),2))
rownames(nuage_clr_relab_genus_correlations) <- HighlyDetectedGenus
colnames(nuage_clr_relab_genus_correlations) <- c("kendall_correlation","study_name")
for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	study_name <- "NUAGE"
	df_temp <- data.frame(clr=as.numeric(nuage_combined_df_sum_stat_genus_clr[,genus_name]),relab=as.numeric(nuage_combined_df_sum_stat_genus[,genus_name]))
	nuage_clr_relab_genus_correlations[genus_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	nuage_clr_relab_genus_correlations[genus_name,2] <- "NUAGE"
	
}


print("Computing Associations at Genus level")

print("rlm bray uniqueness (clr)")
nuage_rlm_detected_genus_bray_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(nuage_rlm_detected_genus_bray_uniqueness_clr) <- HighlyDetectedGenus
colnames(nuage_rlm_detected_genus_bray_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	nuage_rlm_detected_genus_bray_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=nuage_combined_df_sum_stat_genus_clr),var="genus_bray_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_genus_bray_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=nuage_combined_df_sum_stat_genus_clr),var="genus_bray_uniqueness")$p.value
	
}

nuage_rlm_detected_genus_bray_uniqueness_clr$QValue <- p.adjust(nuage_rlm_detected_genus_bray_uniqueness_clr$PValue,method="fdr")

nuage_rlm_detected_genus_bray_uniqueness_clr$Direction <- ifelse(nuage_rlm_detected_genus_bray_uniqueness_clr$QValue<=0.05,3*sign(nuage_rlm_detected_genus_bray_uniqueness_clr$Estimate),ifelse(nuage_rlm_detected_genus_bray_uniqueness_clr$PValue<=0.05,2*sign(nuage_rlm_detected_genus_bray_uniqueness_clr$Estimate),sign(nuage_rlm_detected_genus_bray_uniqueness_clr$Estimate)))

print("rlm bray uniqueness (relab)")
nuage_rlm_detected_genus_bray_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(nuage_rlm_detected_genus_bray_uniqueness) <- HighlyDetectedGenus
colnames(nuage_rlm_detected_genus_bray_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	nuage_rlm_detected_genus_bray_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=nuage_combined_df_sum_stat_genus),var="genus_bray_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_genus_bray_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=nuage_combined_df_sum_stat_genus),var="genus_bray_uniqueness")$p.value
	
}

nuage_rlm_detected_genus_bray_uniqueness$QValue <- p.adjust(nuage_rlm_detected_genus_bray_uniqueness$PValue,method="fdr")

nuage_rlm_detected_genus_bray_uniqueness$Direction <- ifelse(nuage_rlm_detected_genus_bray_uniqueness$QValue<=0.05,3*sign(nuage_rlm_detected_genus_bray_uniqueness$Estimate),ifelse(nuage_rlm_detected_genus_bray_uniqueness$PValue<=0.05,2*sign(nuage_rlm_detected_genus_bray_uniqueness$Estimate),sign(nuage_rlm_detected_genus_bray_uniqueness$Estimate)))

print("rlm jaccard uniqueness (clr)")
nuage_rlm_detected_genus_jaccard_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(nuage_rlm_detected_genus_jaccard_uniqueness_clr) <- HighlyDetectedGenus
colnames(nuage_rlm_detected_genus_jaccard_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	nuage_rlm_detected_genus_jaccard_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=nuage_combined_df_sum_stat_genus_clr),var="genus_jaccard_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_genus_jaccard_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=nuage_combined_df_sum_stat_genus_clr),var="genus_jaccard_uniqueness")$p.value
	
}

nuage_rlm_detected_genus_jaccard_uniqueness_clr$QValue <- p.adjust(nuage_rlm_detected_genus_jaccard_uniqueness_clr$PValue,method="fdr")

nuage_rlm_detected_genus_jaccard_uniqueness_clr$Direction <- ifelse(nuage_rlm_detected_genus_jaccard_uniqueness_clr$QValue<=0.05,3*sign(nuage_rlm_detected_genus_jaccard_uniqueness_clr$Estimate),ifelse(nuage_rlm_detected_genus_jaccard_uniqueness_clr$PValue<=0.05,2*sign(nuage_rlm_detected_genus_jaccard_uniqueness_clr$Estimate),sign(nuage_rlm_detected_genus_jaccard_uniqueness_clr$Estimate)))

print("rlm jaccard uniqueness (relab)")
nuage_rlm_detected_genus_jaccard_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(nuage_rlm_detected_genus_jaccard_uniqueness) <- HighlyDetectedGenus
colnames(nuage_rlm_detected_genus_jaccard_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	nuage_rlm_detected_genus_jaccard_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=nuage_combined_df_sum_stat_genus),var="genus_jaccard_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_genus_jaccard_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=nuage_combined_df_sum_stat_genus),var="genus_jaccard_uniqueness")$p.value
	
}

nuage_rlm_detected_genus_jaccard_uniqueness$QValue <- p.adjust(nuage_rlm_detected_genus_jaccard_uniqueness$PValue,method="fdr")

nuage_rlm_detected_genus_jaccard_uniqueness$Direction <- ifelse(nuage_rlm_detected_genus_jaccard_uniqueness$QValue<=0.05,3*sign(nuage_rlm_detected_genus_jaccard_uniqueness$Estimate),ifelse(nuage_rlm_detected_genus_jaccard_uniqueness$PValue<=0.05,2*sign(nuage_rlm_detected_genus_jaccard_uniqueness$Estimate),sign(nuage_rlm_detected_genus_jaccard_uniqueness$Estimate)))

print("rlm aitchison uniqueness (clr)")
nuage_rlm_detected_genus_aitchison_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(nuage_rlm_detected_genus_aitchison_uniqueness_clr) <- HighlyDetectedGenus
colnames(nuage_rlm_detected_genus_aitchison_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	nuage_rlm_detected_genus_aitchison_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=nuage_combined_df_sum_stat_genus_clr),var="genus_aitchison_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_genus_aitchison_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=nuage_combined_df_sum_stat_genus_clr),var="genus_aitchison_uniqueness")$p.value
	
}

nuage_rlm_detected_genus_aitchison_uniqueness_clr$QValue <- p.adjust(nuage_rlm_detected_genus_aitchison_uniqueness_clr$PValue,method="fdr")

nuage_rlm_detected_genus_aitchison_uniqueness_clr$Direction <- ifelse(nuage_rlm_detected_genus_aitchison_uniqueness_clr$QValue<=0.05,3*sign(nuage_rlm_detected_genus_aitchison_uniqueness_clr$Estimate),ifelse(nuage_rlm_detected_genus_aitchison_uniqueness_clr$PValue<=0.05,2*sign(nuage_rlm_detected_genus_aitchison_uniqueness_clr$Estimate),sign(nuage_rlm_detected_genus_aitchison_uniqueness_clr$Estimate)))

print("rlm aitchison uniqueness (relab)")
nuage_rlm_detected_genus_aitchison_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(nuage_rlm_detected_genus_aitchison_uniqueness) <- HighlyDetectedGenus
colnames(nuage_rlm_detected_genus_aitchison_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	nuage_rlm_detected_genus_aitchison_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=nuage_combined_df_sum_stat_genus),var="genus_aitchison_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_genus_aitchison_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=nuage_combined_df_sum_stat_genus),var="genus_aitchison_uniqueness")$p.value
	
}

nuage_rlm_detected_genus_aitchison_uniqueness$QValue <- p.adjust(nuage_rlm_detected_genus_aitchison_uniqueness$PValue,method="fdr")

nuage_rlm_detected_genus_aitchison_uniqueness$Direction <- ifelse(nuage_rlm_detected_genus_aitchison_uniqueness$QValue<=0.05,3*sign(nuage_rlm_detected_genus_aitchison_uniqueness$Estimate),ifelse(nuage_rlm_detected_genus_aitchison_uniqueness$PValue<=0.05,2*sign(nuage_rlm_detected_genus_aitchison_uniqueness$Estimate),sign(nuage_rlm_detected_genus_aitchison_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
nuage_rlm_detected_genus_kendall_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(nuage_rlm_detected_genus_kendall_uniqueness_clr) <- HighlyDetectedGenus
colnames(nuage_rlm_detected_genus_kendall_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	nuage_rlm_detected_genus_kendall_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=nuage_combined_df_sum_stat_genus_clr),var="genus_kendall_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_genus_kendall_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=nuage_combined_df_sum_stat_genus_clr),var="genus_kendall_uniqueness")$p.value
	
}

nuage_rlm_detected_genus_kendall_uniqueness_clr$QValue <- p.adjust(nuage_rlm_detected_genus_kendall_uniqueness_clr$PValue,method="fdr")

nuage_rlm_detected_genus_kendall_uniqueness_clr$Direction <- ifelse(nuage_rlm_detected_genus_kendall_uniqueness_clr$QValue<=0.05,3*sign(nuage_rlm_detected_genus_kendall_uniqueness_clr$Estimate),ifelse(nuage_rlm_detected_genus_kendall_uniqueness_clr$PValue<=0.05,2*sign(nuage_rlm_detected_genus_kendall_uniqueness_clr$Estimate),sign(nuage_rlm_detected_genus_kendall_uniqueness_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
nuage_rlm_detected_genus_kendall_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(nuage_rlm_detected_genus_kendall_uniqueness) <- HighlyDetectedGenus
colnames(nuage_rlm_detected_genus_kendall_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	nuage_rlm_detected_genus_kendall_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=nuage_combined_df_sum_stat_genus),var="genus_kendall_uniqueness")$coefficients[2,1]
	nuage_rlm_detected_genus_kendall_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=nuage_combined_df_sum_stat_genus),var="genus_kendall_uniqueness")$p.value
	
}

nuage_rlm_detected_genus_kendall_uniqueness$QValue <- p.adjust(nuage_rlm_detected_genus_kendall_uniqueness$PValue,method="fdr")

nuage_rlm_detected_genus_kendall_uniqueness$Direction <- ifelse(nuage_rlm_detected_genus_kendall_uniqueness$QValue<=0.05,3*sign(nuage_rlm_detected_genus_kendall_uniqueness$Estimate),ifelse(nuage_rlm_detected_genus_kendall_uniqueness$PValue<=0.05,2*sign(nuage_rlm_detected_genus_kendall_uniqueness$Estimate),sign(nuage_rlm_detected_genus_kendall_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
nuage_rlm_detected_genus_shannon_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(nuage_rlm_detected_genus_shannon_clr) <- HighlyDetectedGenus
colnames(nuage_rlm_detected_genus_shannon_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	nuage_rlm_detected_genus_shannon_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=nuage_combined_df_sum_stat_genus_clr),var="genus_shannon")$coefficients[2,1]
	nuage_rlm_detected_genus_shannon_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=nuage_combined_df_sum_stat_genus_clr),var="genus_shannon")$p.value
	
}

nuage_rlm_detected_genus_shannon_clr$QValue <- p.adjust(nuage_rlm_detected_genus_shannon_clr$PValue,method="fdr")

nuage_rlm_detected_genus_shannon_clr$Direction <- ifelse(nuage_rlm_detected_genus_shannon_clr$QValue<=0.05,3*sign(nuage_rlm_detected_genus_shannon_clr$Estimate),ifelse(nuage_rlm_detected_genus_shannon_clr$PValue<=0.05,2*sign(nuage_rlm_detected_genus_shannon_clr$Estimate),sign(nuage_rlm_detected_genus_shannon_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
nuage_rlm_detected_genus_shannon <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(nuage_rlm_detected_genus_shannon) <- HighlyDetectedGenus
colnames(nuage_rlm_detected_genus_shannon) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	nuage_rlm_detected_genus_shannon[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=nuage_combined_df_sum_stat_genus),var="genus_shannon")$coefficients[2,1]
	nuage_rlm_detected_genus_shannon[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=nuage_combined_df_sum_stat_genus),var="genus_shannon")$p.value
	
}

nuage_rlm_detected_genus_shannon$QValue <- p.adjust(nuage_rlm_detected_genus_shannon$PValue,method="fdr")

nuage_rlm_detected_genus_shannon$Direction <- ifelse(nuage_rlm_detected_genus_shannon$QValue<=0.05,3*sign(nuage_rlm_detected_genus_shannon$Estimate),ifelse(nuage_rlm_detected_genus_shannon$PValue<=0.05,2*sign(nuage_rlm_detected_genus_shannon$Estimate),sign(nuage_rlm_detected_genus_shannon$Estimate)))

nuage_HighlyDetectedSpecies <- HighlyDetectedSpecies
nuage_HighlyDetectedGenus <- HighlyDetectedGenus

save(nuage_combined_df_sum_stat_species,nuage_rlm_detected_species_shannon,nuage_rlm_detected_species_kendall_uniqueness,nuage_rlm_detected_species_aitchison_uniqueness,nuage_rlm_detected_species_jaccard_uniqueness,nuage_rlm_detected_species_bray_uniqueness,nuage_combined_df_sum_stat_species_clr,nuage_rlm_detected_species_shannon_clr,nuage_rlm_detected_species_kendall_uniqueness_clr,nuage_rlm_detected_species_aitchison_uniqueness_clr,nuage_rlm_detected_species_jaccard_uniqueness_clr,nuage_rlm_detected_species_bray_uniqueness_clr,nuage_combined_df_sum_stat_genus,nuage_rlm_detected_genus_shannon,nuage_rlm_detected_genus_kendall_uniqueness,nuage_rlm_detected_genus_aitchison_uniqueness,nuage_rlm_detected_genus_jaccard_uniqueness,nuage_rlm_detected_genus_bray_uniqueness,nuage_combined_df_sum_stat_genus_clr,nuage_rlm_detected_genus_shannon_clr,nuage_rlm_detected_genus_kendall_uniqueness_clr,nuage_rlm_detected_genus_aitchison_uniqueness_clr,nuage_rlm_detected_genus_jaccard_uniqueness_clr,nuage_rlm_detected_genus_bray_uniqueness_clr,nuage_HighlyDetectedSpecies,nuage_HighlyDetectedGenus,nuage_clr_relab_sample_correlations_species,nuage_clr_relab_sample_correlations_genus,nuage_clr_relab_species_correlations,nuage_clr_relab_genus_correlations,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\nuage_stage2a_results.RData")

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\NUAGE\\nuage_analysis_2021_Revision.RData")

rm(list=ls())
