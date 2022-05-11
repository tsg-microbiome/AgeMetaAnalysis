rm(list=ls())

library(vegan)
library(ade4)
library(MASS)
library(sfsmisc)
library(pcaPP)
library(compositions)

range_scale=function(x)
{
	y <- (x-min(x))/(max(x)-min(x));
	return(y);
}

library(vegan)
library(ade4)
library(MASS)
library(sfsmisc)
library(pcaPP)

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\He\\he_age_analysis.RData")

df_he_diversity_uniqueness <- data.frame(species_shannon = diversity(he_select_age_final_species), genus_shannon = diversity(he_select_age_final_genus), row.names=rownames(he_select_age_final_species))
df_he_diversity_uniqueness$age <- he_select_age_final_metadata[rownames(df_he_diversity_uniqueness),"Age"]
df_he_diversity_uniqueness$region <- he_select_age_final_metadata[rownames(df_he_diversity_uniqueness),"Districts"]
df_he_diversity_uniqueness$study_name <- "HE"
df_he_diversity_uniqueness$BMI <- he_select_age_final_metadata[rownames(df_he_diversity_uniqueness),"BMI"]

print("Genus:Bray")
dist_genus_bray <- as.matrix(vegdist(he_select_age_final_genus,method="bray"))
diag(dist_genus_bray) <- NA
print("Genus:Jaccard")
dist_genus_jaccard <- as.matrix(vegdist(he_select_age_final_genus,method="jaccard"))
diag(dist_genus_jaccard) <- NA

he_select_age_final_genus_clr <- as.matrix(clr(he_select_age_final_genus+0.00001))
he_select_age_final_genus_clr <- as.data.frame(t(apply(he_select_age_final_genus_clr,1,function(x)(x-min(x)))))
print("Genus:Aitchison")
dist_genus_aitchison <- as.matrix(vegdist(he_select_age_final_genus_clr,method="euclidean"))
diag(dist_genus_aitchison) <- NA
print("Genus:Spearman")
dist_genus_kendall <- as.matrix(1-cor.fk(t(he_select_age_final_genus_clr))/2)
diag(dist_genus_kendall) <- NA

print("Species:Bray")
dist_species_bray <- as.matrix(vegdist(he_select_age_final_species,method="bray"))
diag(dist_species_bray) <- NA
print("Species:Jaccard")
dist_species_jaccard <- as.matrix(vegdist(he_select_age_final_species,method="jaccard"))
diag(dist_species_jaccard) <- NA

he_select_age_final_species_clr <- as.matrix(clr(he_select_age_final_species+0.00001))
he_select_age_final_species_clr <- as.data.frame(t(apply(he_select_age_final_species_clr,1,function(x)(x-min(x)))))
print("Species:Aitchison")
dist_species_aitchison <- as.matrix(vegdist(he_select_age_final_species_clr,method="euclidean"))
diag(dist_species_aitchison) <- NA
print("Species:Kendall")
dist_species_kendall <- as.matrix(1-cor.fk(t(he_select_age_final_species_clr))/2)
diag(dist_species_kendall) <- NA

print("Computing Uniqueness Measures")

df_he_diversity_uniqueness$genus_bray_uniqueness <- apply(dist_genus_bray[rownames(df_he_diversity_uniqueness),rownames(df_he_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_he_diversity_uniqueness$genus_jaccard_uniqueness <- apply(dist_genus_jaccard[rownames(df_he_diversity_uniqueness),rownames(df_he_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_he_diversity_uniqueness$genus_aitchison_uniqueness <- apply(dist_genus_aitchison[rownames(df_he_diversity_uniqueness),rownames(df_he_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_he_diversity_uniqueness$genus_kendall_uniqueness <- apply(dist_genus_kendall[rownames(df_he_diversity_uniqueness),rownames(df_he_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_he_diversity_uniqueness$species_bray_uniqueness <- apply(dist_species_bray[rownames(df_he_diversity_uniqueness),rownames(df_he_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_he_diversity_uniqueness$species_jaccard_uniqueness <- apply(dist_species_jaccard[rownames(df_he_diversity_uniqueness),rownames(df_he_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_he_diversity_uniqueness$species_aitchison_uniqueness <- apply(dist_species_aitchison[rownames(df_he_diversity_uniqueness),rownames(df_he_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_he_diversity_uniqueness$species_kendall_uniqueness <- apply(dist_species_kendall[rownames(df_he_diversity_uniqueness),rownames(df_he_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\He\\he_analysis_2021_Revision.RData")

rm(list=ls())








