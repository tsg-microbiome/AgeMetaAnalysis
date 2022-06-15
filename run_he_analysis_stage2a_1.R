# The sub-pipeline works on the He et al data repository and performs the following tasks:
# Computes the sample-to-sample and taxa-to-taxa correlations between the clr-transformed abundances and the relative abundances of the different taxa at the species and genus level for the gut microbiomes of this dataset
# Identifies a set of Highly Detected Species within this dataset.
# Computes the association between the clr transformed abundances of this set of Highly Detected Species 
# and the different microbiome summary statistic measures (four measures of uniqueness along with Shannon Diversity) using Robust Linear Regression models

# The sub-pipeline requires the workspace he_stage1_results.RData which is the output of the previous two sub-pipelines run_he_analysis_stage1a.R, followed by
# run_he_analysis_stage1b.R, along with the he_Analysis_2021_Revision.RData which is the running workspace of the he et al dataset analysis (after the execution
# of run_he_analysis_stage1a.R and run_he_analysis_stage1b.R)

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\He\\he_age_analysis.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\he_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\He\\he_Analysis_2021_Revision.RData")

rm(he_combined_df_sum_stat_species)
rm(he_combined_df_sum_stat_species_clr)
rm(he_combined_df_sum_stat_genus)
rm(he_combined_df_sum_stat_genus_clr)

df_he_diversity_uniqueness <- as.data.frame(df_he_diversity_uniqueness[,c("species_shannon","genus_shannon","age","study_name","BMI","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")])

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

HighlyDetectedSpecies <- names(which(100*colSums(apply(he_select_age_final_species,2,function(x)(ifelse(x!=0,1,0))))/nrow(he_select_age_final_species)>=10))

he_combined_df_sum_stat_species <- as.data.frame(cbind(df_he_diversity_uniqueness,he_select_age_final_species[rownames(df_he_diversity_uniqueness),HighlyDetectedSpecies]))

he_combined_df_sum_stat_species_clr <- as.data.frame(cbind(df_he_diversity_uniqueness,he_select_age_final_species_clr[rownames(df_he_diversity_uniqueness),HighlyDetectedSpecies]))

he_clr_relab_sample_correlations_species <- as.data.frame(matrix(NA,nrow(he_combined_df_sum_stat_species),2))
rownames(he_clr_relab_sample_correlations_species) <- rownames(he_combined_df_sum_stat_species)
colnames(he_clr_relab_sample_correlations_species) <- c("study_name","kendall_correlation")
for(i in 1:nrow(he_combined_df_sum_stat_species))
{
	sample_name <- rownames(he_combined_df_sum_stat_species)[i]
	study_name <- "HE"
	df_temp <- data.frame(clr=as.numeric(he_combined_df_sum_stat_species_clr[sample_name,HighlyDetectedSpecies]),relab=as.numeric(he_combined_df_sum_stat_species[sample_name,HighlyDetectedSpecies]))
	he_clr_relab_sample_correlations_species[sample_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	he_clr_relab_sample_correlations_species[sample_name,2] <- study_name
}

he_clr_relab_species_correlations <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),2))
rownames(he_clr_relab_species_correlations) <- HighlyDetectedSpecies
colnames(he_clr_relab_species_correlations) <- c("kendall_correlation","study_name")
for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	study_name <- "HE"
	df_temp <- data.frame(clr=as.numeric(he_combined_df_sum_stat_species_clr[,species_name]),relab=as.numeric(he_combined_df_sum_stat_species[,species_name]))
	he_clr_relab_species_correlations[species_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	he_clr_relab_species_correlations[species_name,2] <- "HE"
	
}

print("Computing Associations at Species level")

print("rlm bray uniqueness (clr)")
he_rlm_detected_species_bray_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(he_rlm_detected_species_bray_uniqueness_clr) <- HighlyDetectedSpecies
colnames(he_rlm_detected_species_bray_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	he_rlm_detected_species_bray_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=he_combined_df_sum_stat_species_clr),var="species_bray_uniqueness")$coefficients[2,1]
	he_rlm_detected_species_bray_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=he_combined_df_sum_stat_species_clr),var="species_bray_uniqueness")$p.value
	
}

he_rlm_detected_species_bray_uniqueness_clr$QValue <- p.adjust(he_rlm_detected_species_bray_uniqueness_clr$PValue,method="fdr")

he_rlm_detected_species_bray_uniqueness_clr$Direction <- ifelse(he_rlm_detected_species_bray_uniqueness_clr$QValue<=0.05,3*sign(he_rlm_detected_species_bray_uniqueness_clr$Estimate),ifelse(he_rlm_detected_species_bray_uniqueness_clr$PValue<=0.05,2*sign(he_rlm_detected_species_bray_uniqueness_clr$Estimate),sign(he_rlm_detected_species_bray_uniqueness_clr$Estimate)))

print("rlm bray uniqueness (relab)")
he_rlm_detected_species_bray_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(he_rlm_detected_species_bray_uniqueness) <- HighlyDetectedSpecies
colnames(he_rlm_detected_species_bray_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	he_rlm_detected_species_bray_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=he_combined_df_sum_stat_species),var="species_bray_uniqueness")$coefficients[2,1]
	he_rlm_detected_species_bray_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=he_combined_df_sum_stat_species),var="species_bray_uniqueness")$p.value
	
}

he_rlm_detected_species_bray_uniqueness$QValue <- p.adjust(he_rlm_detected_species_bray_uniqueness$PValue,method="fdr")

he_rlm_detected_species_bray_uniqueness$Direction <- ifelse(he_rlm_detected_species_bray_uniqueness$QValue<=0.05,3*sign(he_rlm_detected_species_bray_uniqueness$Estimate),ifelse(he_rlm_detected_species_bray_uniqueness$PValue<=0.05,2*sign(he_rlm_detected_species_bray_uniqueness$Estimate),sign(he_rlm_detected_species_bray_uniqueness$Estimate)))

print("rlm jaccard uniqueness (clr)")
he_rlm_detected_species_jaccard_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(he_rlm_detected_species_jaccard_uniqueness_clr) <- HighlyDetectedSpecies
colnames(he_rlm_detected_species_jaccard_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	he_rlm_detected_species_jaccard_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=he_combined_df_sum_stat_species_clr),var="species_jaccard_uniqueness")$coefficients[2,1]
	he_rlm_detected_species_jaccard_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=he_combined_df_sum_stat_species_clr),var="species_jaccard_uniqueness")$p.value
	
}

he_rlm_detected_species_jaccard_uniqueness_clr$QValue <- p.adjust(he_rlm_detected_species_jaccard_uniqueness_clr$PValue,method="fdr")

he_rlm_detected_species_jaccard_uniqueness_clr$Direction <- ifelse(he_rlm_detected_species_jaccard_uniqueness_clr$QValue<=0.05,3*sign(he_rlm_detected_species_jaccard_uniqueness_clr$Estimate),ifelse(he_rlm_detected_species_jaccard_uniqueness_clr$PValue<=0.05,2*sign(he_rlm_detected_species_jaccard_uniqueness_clr$Estimate),sign(he_rlm_detected_species_jaccard_uniqueness_clr$Estimate)))

print("rlm jaccard uniqueness (relab)")
he_rlm_detected_species_jaccard_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(he_rlm_detected_species_jaccard_uniqueness) <- HighlyDetectedSpecies
colnames(he_rlm_detected_species_jaccard_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	he_rlm_detected_species_jaccard_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=he_combined_df_sum_stat_species),var="species_jaccard_uniqueness")$coefficients[2,1]
	he_rlm_detected_species_jaccard_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=he_combined_df_sum_stat_species),var="species_jaccard_uniqueness")$p.value
	
}

he_rlm_detected_species_jaccard_uniqueness$QValue <- p.adjust(he_rlm_detected_species_jaccard_uniqueness$PValue,method="fdr")

he_rlm_detected_species_jaccard_uniqueness$Direction <- ifelse(he_rlm_detected_species_jaccard_uniqueness$QValue<=0.05,3*sign(he_rlm_detected_species_jaccard_uniqueness$Estimate),ifelse(he_rlm_detected_species_jaccard_uniqueness$PValue<=0.05,2*sign(he_rlm_detected_species_jaccard_uniqueness$Estimate),sign(he_rlm_detected_species_jaccard_uniqueness$Estimate)))

print("rlm aitchison uniqueness (clr)")
he_rlm_detected_species_aitchison_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(he_rlm_detected_species_aitchison_uniqueness_clr) <- HighlyDetectedSpecies
colnames(he_rlm_detected_species_aitchison_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	he_rlm_detected_species_aitchison_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=he_combined_df_sum_stat_species_clr),var="species_aitchison_uniqueness")$coefficients[2,1]
	he_rlm_detected_species_aitchison_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=he_combined_df_sum_stat_species_clr),var="species_aitchison_uniqueness")$p.value
	
}

he_rlm_detected_species_aitchison_uniqueness_clr$QValue <- p.adjust(he_rlm_detected_species_aitchison_uniqueness_clr$PValue,method="fdr")

he_rlm_detected_species_aitchison_uniqueness_clr$Direction <- ifelse(he_rlm_detected_species_aitchison_uniqueness_clr$QValue<=0.05,3*sign(he_rlm_detected_species_aitchison_uniqueness_clr$Estimate),ifelse(he_rlm_detected_species_aitchison_uniqueness_clr$PValue<=0.05,2*sign(he_rlm_detected_species_aitchison_uniqueness_clr$Estimate),sign(he_rlm_detected_species_aitchison_uniqueness_clr$Estimate)))

print("rlm aitchison uniqueness (relab)")
he_rlm_detected_species_aitchison_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(he_rlm_detected_species_aitchison_uniqueness) <- HighlyDetectedSpecies
colnames(he_rlm_detected_species_aitchison_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	he_rlm_detected_species_aitchison_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=he_combined_df_sum_stat_species),var="species_aitchison_uniqueness")$coefficients[2,1]
	he_rlm_detected_species_aitchison_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=he_combined_df_sum_stat_species),var="species_aitchison_uniqueness")$p.value
	
}

he_rlm_detected_species_aitchison_uniqueness$QValue <- p.adjust(he_rlm_detected_species_aitchison_uniqueness$PValue,method="fdr")

he_rlm_detected_species_aitchison_uniqueness$Direction <- ifelse(he_rlm_detected_species_aitchison_uniqueness$QValue<=0.05,3*sign(he_rlm_detected_species_aitchison_uniqueness$Estimate),ifelse(he_rlm_detected_species_aitchison_uniqueness$PValue<=0.05,2*sign(he_rlm_detected_species_aitchison_uniqueness$Estimate),sign(he_rlm_detected_species_aitchison_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
he_rlm_detected_species_kendall_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(he_rlm_detected_species_kendall_uniqueness_clr) <- HighlyDetectedSpecies
colnames(he_rlm_detected_species_kendall_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	he_rlm_detected_species_kendall_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=he_combined_df_sum_stat_species_clr),var="species_kendall_uniqueness")$coefficients[2,1]
	he_rlm_detected_species_kendall_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=he_combined_df_sum_stat_species_clr),var="species_kendall_uniqueness")$p.value
	
}

he_rlm_detected_species_kendall_uniqueness_clr$QValue <- p.adjust(he_rlm_detected_species_kendall_uniqueness_clr$PValue,method="fdr")

he_rlm_detected_species_kendall_uniqueness_clr$Direction <- ifelse(he_rlm_detected_species_kendall_uniqueness_clr$QValue<=0.05,3*sign(he_rlm_detected_species_kendall_uniqueness_clr$Estimate),ifelse(he_rlm_detected_species_kendall_uniqueness_clr$PValue<=0.05,2*sign(he_rlm_detected_species_kendall_uniqueness_clr$Estimate),sign(he_rlm_detected_species_kendall_uniqueness_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
he_rlm_detected_species_kendall_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(he_rlm_detected_species_kendall_uniqueness) <- HighlyDetectedSpecies
colnames(he_rlm_detected_species_kendall_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	he_rlm_detected_species_kendall_uniqueness[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=he_combined_df_sum_stat_species),var="species_kendall_uniqueness")$coefficients[2,1]
	he_rlm_detected_species_kendall_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=he_combined_df_sum_stat_species),var="species_kendall_uniqueness")$p.value
	
}

he_rlm_detected_species_kendall_uniqueness$QValue <- p.adjust(he_rlm_detected_species_kendall_uniqueness$PValue,method="fdr")

he_rlm_detected_species_kendall_uniqueness$Direction <- ifelse(he_rlm_detected_species_kendall_uniqueness$QValue<=0.05,3*sign(he_rlm_detected_species_kendall_uniqueness$Estimate),ifelse(he_rlm_detected_species_kendall_uniqueness$PValue<=0.05,2*sign(he_rlm_detected_species_kendall_uniqueness$Estimate),sign(he_rlm_detected_species_kendall_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
he_rlm_detected_species_shannon_clr <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(he_rlm_detected_species_shannon_clr) <- HighlyDetectedSpecies
colnames(he_rlm_detected_species_shannon_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	he_rlm_detected_species_shannon_clr[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_shannon")),data=he_combined_df_sum_stat_species_clr),var="species_shannon")$coefficients[2,1]
	he_rlm_detected_species_shannon_clr[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_shannon")),data=he_combined_df_sum_stat_species_clr),var="species_shannon")$p.value
	
}

he_rlm_detected_species_shannon_clr$QValue <- p.adjust(he_rlm_detected_species_shannon_clr$PValue,method="fdr")

he_rlm_detected_species_shannon_clr$Direction <- ifelse(he_rlm_detected_species_shannon_clr$QValue<=0.05,3*sign(he_rlm_detected_species_shannon_clr$Estimate),ifelse(he_rlm_detected_species_shannon_clr$PValue<=0.05,2*sign(he_rlm_detected_species_shannon_clr$Estimate),sign(he_rlm_detected_species_shannon_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
he_rlm_detected_species_shannon <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),4))
rownames(he_rlm_detected_species_shannon) <- HighlyDetectedSpecies
colnames(he_rlm_detected_species_shannon) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	he_rlm_detected_species_shannon[i,1] <- summary(rlm(as.formula(paste0(species_name,"~species_shannon")),data=he_combined_df_sum_stat_species),var="species_shannon")$coefficients[2,1]
	he_rlm_detected_species_shannon[i,2] <- f.robftest(rlm(as.formula(paste0(species_name,"~species_shannon")),data=he_combined_df_sum_stat_species),var="species_shannon")$p.value
	
}

he_rlm_detected_species_shannon$QValue <- p.adjust(he_rlm_detected_species_shannon$PValue,method="fdr")

he_rlm_detected_species_shannon$Direction <- ifelse(he_rlm_detected_species_shannon$QValue<=0.05,3*sign(he_rlm_detected_species_shannon$Estimate),ifelse(he_rlm_detected_species_shannon$PValue<=0.05,2*sign(he_rlm_detected_species_shannon$Estimate),sign(he_rlm_detected_species_shannon$Estimate)))

print("Computing Associations at Genus level")

HighlyDetectedGenus <- names(which(100*colSums(apply(he_select_age_final_genus,2,function(x)(ifelse(x!=0,1,0))))/nrow(he_select_age_final_genus)>=10))

he_combined_df_sum_stat_genus <- as.data.frame(cbind(df_he_diversity_uniqueness,he_select_age_final_genus[rownames(df_he_diversity_uniqueness),HighlyDetectedGenus]))

he_combined_df_sum_stat_genus_clr <- as.data.frame(cbind(df_he_diversity_uniqueness,he_select_age_final_genus_clr[rownames(df_he_diversity_uniqueness),HighlyDetectedGenus]))

colnames(he_combined_df_sum_stat_genus) <- sub("-","_",colnames(he_combined_df_sum_stat_genus))
colnames(he_combined_df_sum_stat_genus_clr) <- sub("-","_",colnames(he_combined_df_sum_stat_genus_clr))
HighlyDetectedGenus <- sub("-","_",HighlyDetectedGenus)

he_clr_relab_sample_correlations_genus <- as.data.frame(matrix(NA,nrow(he_combined_df_sum_stat_genus),2))
rownames(he_clr_relab_sample_correlations_genus) <- rownames(he_combined_df_sum_stat_genus)
colnames(he_clr_relab_sample_correlations_genus) <- c("study_name","kendall_correlation")
for(i in 1:nrow(he_combined_df_sum_stat_genus))
{
	sample_name <- rownames(he_combined_df_sum_stat_genus)[i]
	study_name <- "HE"
	df_temp <- data.frame(clr=as.numeric(he_combined_df_sum_stat_genus_clr[sample_name,HighlyDetectedGenus]),relab=as.numeric(he_combined_df_sum_stat_genus[sample_name,HighlyDetectedGenus]))
	he_clr_relab_sample_correlations_genus[sample_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	he_clr_relab_sample_correlations_genus[sample_name,2] <- study_name
}

he_clr_relab_genus_correlations <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),2))
rownames(he_clr_relab_genus_correlations) <- HighlyDetectedGenus
colnames(he_clr_relab_genus_correlations) <- c("kendall_correlation","study_name")
for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	study_name <- "HE"
	df_temp <- data.frame(clr=as.numeric(he_combined_df_sum_stat_genus_clr[,genus_name]),relab=as.numeric(he_combined_df_sum_stat_genus[,genus_name]))
	he_clr_relab_genus_correlations[genus_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	he_clr_relab_genus_correlations[genus_name,2] <- "HE"
	
}


print("Computing Associations at Genus level")

print("rlm bray uniqueness (clr)")
he_rlm_detected_genus_bray_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(he_rlm_detected_genus_bray_uniqueness_clr) <- HighlyDetectedGenus
colnames(he_rlm_detected_genus_bray_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	he_rlm_detected_genus_bray_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=he_combined_df_sum_stat_genus_clr),var="genus_bray_uniqueness")$coefficients[2,1]
	he_rlm_detected_genus_bray_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=he_combined_df_sum_stat_genus_clr),var="genus_bray_uniqueness")$p.value
	
}

he_rlm_detected_genus_bray_uniqueness_clr$QValue <- p.adjust(he_rlm_detected_genus_bray_uniqueness_clr$PValue,method="fdr")

he_rlm_detected_genus_bray_uniqueness_clr$Direction <- ifelse(he_rlm_detected_genus_bray_uniqueness_clr$QValue<=0.05,3*sign(he_rlm_detected_genus_bray_uniqueness_clr$Estimate),ifelse(he_rlm_detected_genus_bray_uniqueness_clr$PValue<=0.05,2*sign(he_rlm_detected_genus_bray_uniqueness_clr$Estimate),sign(he_rlm_detected_genus_bray_uniqueness_clr$Estimate)))

print("rlm bray uniqueness (relab)")
he_rlm_detected_genus_bray_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(he_rlm_detected_genus_bray_uniqueness) <- HighlyDetectedGenus
colnames(he_rlm_detected_genus_bray_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	he_rlm_detected_genus_bray_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=he_combined_df_sum_stat_genus),var="genus_bray_uniqueness")$coefficients[2,1]
	he_rlm_detected_genus_bray_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=he_combined_df_sum_stat_genus),var="genus_bray_uniqueness")$p.value
	
}

he_rlm_detected_genus_bray_uniqueness$QValue <- p.adjust(he_rlm_detected_genus_bray_uniqueness$PValue,method="fdr")

he_rlm_detected_genus_bray_uniqueness$Direction <- ifelse(he_rlm_detected_genus_bray_uniqueness$QValue<=0.05,3*sign(he_rlm_detected_genus_bray_uniqueness$Estimate),ifelse(he_rlm_detected_genus_bray_uniqueness$PValue<=0.05,2*sign(he_rlm_detected_genus_bray_uniqueness$Estimate),sign(he_rlm_detected_genus_bray_uniqueness$Estimate)))

print("rlm jaccard uniqueness (clr)")
he_rlm_detected_genus_jaccard_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(he_rlm_detected_genus_jaccard_uniqueness_clr) <- HighlyDetectedGenus
colnames(he_rlm_detected_genus_jaccard_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	he_rlm_detected_genus_jaccard_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=he_combined_df_sum_stat_genus_clr),var="genus_jaccard_uniqueness")$coefficients[2,1]
	he_rlm_detected_genus_jaccard_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=he_combined_df_sum_stat_genus_clr),var="genus_jaccard_uniqueness")$p.value
	
}

he_rlm_detected_genus_jaccard_uniqueness_clr$QValue <- p.adjust(he_rlm_detected_genus_jaccard_uniqueness_clr$PValue,method="fdr")

he_rlm_detected_genus_jaccard_uniqueness_clr$Direction <- ifelse(he_rlm_detected_genus_jaccard_uniqueness_clr$QValue<=0.05,3*sign(he_rlm_detected_genus_jaccard_uniqueness_clr$Estimate),ifelse(he_rlm_detected_genus_jaccard_uniqueness_clr$PValue<=0.05,2*sign(he_rlm_detected_genus_jaccard_uniqueness_clr$Estimate),sign(he_rlm_detected_genus_jaccard_uniqueness_clr$Estimate)))

print("rlm jaccard uniqueness (relab)")
he_rlm_detected_genus_jaccard_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(he_rlm_detected_genus_jaccard_uniqueness) <- HighlyDetectedGenus
colnames(he_rlm_detected_genus_jaccard_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	he_rlm_detected_genus_jaccard_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=he_combined_df_sum_stat_genus),var="genus_jaccard_uniqueness")$coefficients[2,1]
	he_rlm_detected_genus_jaccard_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=he_combined_df_sum_stat_genus),var="genus_jaccard_uniqueness")$p.value
	
}

he_rlm_detected_genus_jaccard_uniqueness$QValue <- p.adjust(he_rlm_detected_genus_jaccard_uniqueness$PValue,method="fdr")

he_rlm_detected_genus_jaccard_uniqueness$Direction <- ifelse(he_rlm_detected_genus_jaccard_uniqueness$QValue<=0.05,3*sign(he_rlm_detected_genus_jaccard_uniqueness$Estimate),ifelse(he_rlm_detected_genus_jaccard_uniqueness$PValue<=0.05,2*sign(he_rlm_detected_genus_jaccard_uniqueness$Estimate),sign(he_rlm_detected_genus_jaccard_uniqueness$Estimate)))

print("rlm aitchison uniqueness (clr)")
he_rlm_detected_genus_aitchison_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(he_rlm_detected_genus_aitchison_uniqueness_clr) <- HighlyDetectedGenus
colnames(he_rlm_detected_genus_aitchison_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	he_rlm_detected_genus_aitchison_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=he_combined_df_sum_stat_genus_clr),var="genus_aitchison_uniqueness")$coefficients[2,1]
	he_rlm_detected_genus_aitchison_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=he_combined_df_sum_stat_genus_clr),var="genus_aitchison_uniqueness")$p.value
	
}

he_rlm_detected_genus_aitchison_uniqueness_clr$QValue <- p.adjust(he_rlm_detected_genus_aitchison_uniqueness_clr$PValue,method="fdr")

he_rlm_detected_genus_aitchison_uniqueness_clr$Direction <- ifelse(he_rlm_detected_genus_aitchison_uniqueness_clr$QValue<=0.05,3*sign(he_rlm_detected_genus_aitchison_uniqueness_clr$Estimate),ifelse(he_rlm_detected_genus_aitchison_uniqueness_clr$PValue<=0.05,2*sign(he_rlm_detected_genus_aitchison_uniqueness_clr$Estimate),sign(he_rlm_detected_genus_aitchison_uniqueness_clr$Estimate)))

print("rlm aitchison uniqueness (relab)")
he_rlm_detected_genus_aitchison_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(he_rlm_detected_genus_aitchison_uniqueness) <- HighlyDetectedGenus
colnames(he_rlm_detected_genus_aitchison_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	he_rlm_detected_genus_aitchison_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=he_combined_df_sum_stat_genus),var="genus_aitchison_uniqueness")$coefficients[2,1]
	he_rlm_detected_genus_aitchison_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=he_combined_df_sum_stat_genus),var="genus_aitchison_uniqueness")$p.value
	
}

he_rlm_detected_genus_aitchison_uniqueness$QValue <- p.adjust(he_rlm_detected_genus_aitchison_uniqueness$PValue,method="fdr")

he_rlm_detected_genus_aitchison_uniqueness$Direction <- ifelse(he_rlm_detected_genus_aitchison_uniqueness$QValue<=0.05,3*sign(he_rlm_detected_genus_aitchison_uniqueness$Estimate),ifelse(he_rlm_detected_genus_aitchison_uniqueness$PValue<=0.05,2*sign(he_rlm_detected_genus_aitchison_uniqueness$Estimate),sign(he_rlm_detected_genus_aitchison_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
he_rlm_detected_genus_kendall_uniqueness_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(he_rlm_detected_genus_kendall_uniqueness_clr) <- HighlyDetectedGenus
colnames(he_rlm_detected_genus_kendall_uniqueness_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	he_rlm_detected_genus_kendall_uniqueness_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=he_combined_df_sum_stat_genus_clr),var="genus_kendall_uniqueness")$coefficients[2,1]
	he_rlm_detected_genus_kendall_uniqueness_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=he_combined_df_sum_stat_genus_clr),var="genus_kendall_uniqueness")$p.value
	
}

he_rlm_detected_genus_kendall_uniqueness_clr$QValue <- p.adjust(he_rlm_detected_genus_kendall_uniqueness_clr$PValue,method="fdr")

he_rlm_detected_genus_kendall_uniqueness_clr$Direction <- ifelse(he_rlm_detected_genus_kendall_uniqueness_clr$QValue<=0.05,3*sign(he_rlm_detected_genus_kendall_uniqueness_clr$Estimate),ifelse(he_rlm_detected_genus_kendall_uniqueness_clr$PValue<=0.05,2*sign(he_rlm_detected_genus_kendall_uniqueness_clr$Estimate),sign(he_rlm_detected_genus_kendall_uniqueness_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
he_rlm_detected_genus_kendall_uniqueness <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(he_rlm_detected_genus_kendall_uniqueness) <- HighlyDetectedGenus
colnames(he_rlm_detected_genus_kendall_uniqueness) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	he_rlm_detected_genus_kendall_uniqueness[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=he_combined_df_sum_stat_genus),var="genus_kendall_uniqueness")$coefficients[2,1]
	he_rlm_detected_genus_kendall_uniqueness[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=he_combined_df_sum_stat_genus),var="genus_kendall_uniqueness")$p.value
	
}

he_rlm_detected_genus_kendall_uniqueness$QValue <- p.adjust(he_rlm_detected_genus_kendall_uniqueness$PValue,method="fdr")

he_rlm_detected_genus_kendall_uniqueness$Direction <- ifelse(he_rlm_detected_genus_kendall_uniqueness$QValue<=0.05,3*sign(he_rlm_detected_genus_kendall_uniqueness$Estimate),ifelse(he_rlm_detected_genus_kendall_uniqueness$PValue<=0.05,2*sign(he_rlm_detected_genus_kendall_uniqueness$Estimate),sign(he_rlm_detected_genus_kendall_uniqueness$Estimate)))

print("rlm kendall uniqueness (clr)")
he_rlm_detected_genus_shannon_clr <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(he_rlm_detected_genus_shannon_clr) <- HighlyDetectedGenus
colnames(he_rlm_detected_genus_shannon_clr) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	he_rlm_detected_genus_shannon_clr[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=he_combined_df_sum_stat_genus_clr),var="genus_shannon")$coefficients[2,1]
	he_rlm_detected_genus_shannon_clr[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=he_combined_df_sum_stat_genus_clr),var="genus_shannon")$p.value
	
}

he_rlm_detected_genus_shannon_clr$QValue <- p.adjust(he_rlm_detected_genus_shannon_clr$PValue,method="fdr")

he_rlm_detected_genus_shannon_clr$Direction <- ifelse(he_rlm_detected_genus_shannon_clr$QValue<=0.05,3*sign(he_rlm_detected_genus_shannon_clr$Estimate),ifelse(he_rlm_detected_genus_shannon_clr$PValue<=0.05,2*sign(he_rlm_detected_genus_shannon_clr$Estimate),sign(he_rlm_detected_genus_shannon_clr$Estimate)))

print("rlm kendall uniqueness (relab)")
he_rlm_detected_genus_shannon <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),4))
rownames(he_rlm_detected_genus_shannon) <- HighlyDetectedGenus
colnames(he_rlm_detected_genus_shannon) <- c("Estimate","PValue","QValue","Direction")

for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	he_rlm_detected_genus_shannon[i,1] <- summary(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=he_combined_df_sum_stat_genus),var="genus_shannon")$coefficients[2,1]
	he_rlm_detected_genus_shannon[i,2] <- f.robftest(rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=he_combined_df_sum_stat_genus),var="genus_shannon")$p.value
	
}

he_rlm_detected_genus_shannon$QValue <- p.adjust(he_rlm_detected_genus_shannon$PValue,method="fdr")

he_rlm_detected_genus_shannon$Direction <- ifelse(he_rlm_detected_genus_shannon$QValue<=0.05,3*sign(he_rlm_detected_genus_shannon$Estimate),ifelse(he_rlm_detected_genus_shannon$PValue<=0.05,2*sign(he_rlm_detected_genus_shannon$Estimate),sign(he_rlm_detected_genus_shannon$Estimate)))

he_HighlyDetectedSpecies <- HighlyDetectedSpecies
he_HighlyDetectedGenus <- HighlyDetectedGenus

save(he_combined_df_sum_stat_species,he_rlm_detected_species_shannon,he_rlm_detected_species_kendall_uniqueness,he_rlm_detected_species_aitchison_uniqueness,he_rlm_detected_species_jaccard_uniqueness,he_rlm_detected_species_bray_uniqueness,he_combined_df_sum_stat_species_clr,he_rlm_detected_species_shannon_clr,he_rlm_detected_species_kendall_uniqueness_clr,he_rlm_detected_species_aitchison_uniqueness_clr,he_rlm_detected_species_jaccard_uniqueness_clr,he_rlm_detected_species_bray_uniqueness_clr,he_combined_df_sum_stat_genus,he_rlm_detected_genus_shannon,he_rlm_detected_genus_kendall_uniqueness,he_rlm_detected_genus_aitchison_uniqueness,he_rlm_detected_genus_jaccard_uniqueness,he_rlm_detected_genus_bray_uniqueness,he_combined_df_sum_stat_genus_clr,he_rlm_detected_genus_shannon_clr,he_rlm_detected_genus_kendall_uniqueness_clr,he_rlm_detected_genus_aitchison_uniqueness_clr,he_rlm_detected_genus_jaccard_uniqueness_clr,he_rlm_detected_genus_bray_uniqueness_clr,he_HighlyDetectedSpecies,he_HighlyDetectedGenus,he_clr_relab_sample_correlations_species,he_clr_relab_sample_correlations_genus,he_clr_relab_species_correlations,he_clr_relab_genus_correlations,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\he_stage2a_results.RData")

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\He\\he_analysis_2021_Revision.RData")

rm(list=ls())
