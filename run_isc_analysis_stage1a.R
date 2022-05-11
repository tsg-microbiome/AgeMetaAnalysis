rm(list=ls())

library(vegan)
library(ade4)
library(MASS)
library(sfsmisc)
library(pcaPP)
library(robumeta)
library(metafor)
library(dplyr)
library(effsize)
library(compositions)

range_scale=function(x)
{
	y <- (x-min(x))/(max(x)-min(x));
	return(y);
}

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\isc_age_analysis.RData")

df_isc_diversity_uniqueness <- data.frame(species_shannon = diversity(isc_select_age_final_species), genus_shannon = diversity(isc_select_age_final_genus), row.names=rownames(isc_select_age_final_species))
df_isc_diversity_uniqueness$age <- isc_select_age_final_metadata[rownames(df_isc_diversity_uniqueness),"age"]
df_isc_diversity_uniqueness$study_condition <- isc_select_age_final_metadata[rownames(df_isc_diversity_uniqueness),"study_condition"]
df_isc_diversity_uniqueness$study_name <- isc_select_age_final_metadata[rownames(df_isc_diversity_uniqueness),"study_name"]


print("Genus:Bray")
dist_genus_bray <- as.matrix(vegdist(isc_select_age_final_genus,method="bray"))
diag(dist_genus_bray) <- NA
print("Genus:Jaccard")
dist_genus_jaccard <- as.matrix(vegdist(isc_select_age_final_genus,method="jaccard"))
diag(dist_genus_jaccard) <- NA

isc_select_age_final_genus_clr <- as.matrix(clr(isc_select_age_final_genus+0.00001))
isc_select_age_final_genus_clr <- as.data.frame(t(apply(isc_select_age_final_genus_clr,1,function(x)(x-min(x)))))
print("Genus:Manhattan")
dist_genus_aitchison <- as.matrix(vegdist(isc_select_age_final_genus_clr,method="euclidean"))
diag(dist_genus_aitchison) <- NA
print("Genus:Kendall")
dist_genus_kendall <- as.matrix(1-cor.fk(t(isc_select_age_final_genus_clr))/2)
diag(dist_genus_kendall) <- NA

print("Species:Bray")
dist_species_bray <- as.matrix(vegdist(isc_select_age_final_species,method="bray"))
diag(dist_species_bray) <- NA
print("Species:Jaccard")
dist_species_jaccard <- as.matrix(vegdist(isc_select_age_final_species,method="jaccard"))
diag(dist_species_jaccard) <- NA

isc_select_age_final_species_clr <- as.matrix(clr(isc_select_age_final_species+0.00001))
isc_select_age_final_species_clr <- as.data.frame(t(apply(isc_select_age_final_species_clr,1,function(x)(x-min(x)))))
print("Species:Manhattan")
dist_species_aitchison <- as.matrix(vegdist(isc_select_age_final_species_clr,method="euclidean"))
diag(dist_species_aitchison) <- NA
print("Species:Kendall")
dist_species_kendall <- as.matrix(1-cor.fk(t(isc_select_age_final_species_clr))/2)
diag(dist_species_kendall) <- NA

print("Computing Uniqueness Measures")

df_isc_diversity_uniqueness$genus_bray_uniqueness <- apply(dist_genus_bray[rownames(df_isc_diversity_uniqueness),rownames(df_isc_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_isc_diversity_uniqueness$genus_jaccard_uniqueness <- apply(dist_genus_jaccard[rownames(df_isc_diversity_uniqueness),rownames(df_isc_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_isc_diversity_uniqueness$genus_aitchison_uniqueness <- apply(dist_genus_aitchison[rownames(df_isc_diversity_uniqueness),rownames(df_isc_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_isc_diversity_uniqueness$genus_kendall_uniqueness <- apply(dist_genus_kendall[rownames(df_isc_diversity_uniqueness),rownames(df_isc_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_isc_diversity_uniqueness$species_bray_uniqueness <- apply(dist_species_bray[rownames(df_isc_diversity_uniqueness),rownames(df_isc_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_isc_diversity_uniqueness$species_jaccard_uniqueness <- apply(dist_species_jaccard[rownames(df_isc_diversity_uniqueness),rownames(df_isc_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_isc_diversity_uniqueness$species_aitchison_uniqueness <- apply(dist_species_aitchison[rownames(df_isc_diversity_uniqueness),rownames(df_isc_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_isc_diversity_uniqueness$species_kendall_uniqueness <- apply(dist_species_kendall[rownames(df_isc_diversity_uniqueness),rownames(df_isc_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\isc_analysis_2021_Revision.RData")

rm(list=ls())







