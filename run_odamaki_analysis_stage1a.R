# The sub-pipeline scans/inputs the species, genus and metadata tables for the Odamaki et al data repository.
# Computes Shannon Diversity and the four different taxonomic uniqueness at the Species and Genus levels across all gut microbiomes
# The output is a single data frame with the select metadata and the diversity and uniqueness measures for all gut microbiomes

rm(list=ls())
library(vegan)
library(ade4)
library(MASS)
library(sfsmisc)
library(compositions)

range_scale=function(x)
{
	y <- (x-min(x))/(max(x)-min(x));
	return(y);
}

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\Odamaki\\odamaki_age_analysis.RData")

odamaki_select_age_final_species <- t(apply(odamaki_select_age_final_species,1,function(x)(x/sum(x))))
odamaki_select_age_final_genus <- t(apply(odamaki_select_age_final_genus,1,function(x)(x/sum(x))))

df_odamaki_diversity_uniqueness <- data.frame(species_shannon = diversity(odamaki_select_age_final_species), genus_shannon = diversity(odamaki_select_age_final_genus), row.names=rownames(odamaki_select_age_final_species))
df_odamaki_diversity_uniqueness$age <- odamaki_select_age_final_metadata[rownames(df_odamaki_diversity_uniqueness),"Age"]
df_odamaki_diversity_uniqueness$study_name <- "Odamaki"

print("Genus:Bray")
dist_genus_bray <- as.matrix(vegdist(odamaki_select_age_final_genus,method="bray"))
diag(dist_genus_bray) <- NA
print("Genus:Jaccard")
dist_genus_jaccard <- as.matrix(vegdist(odamaki_select_age_final_genus,method="jaccard"))
diag(dist_genus_jaccard) <- NA

odamaki_select_age_final_genus_clr <- as.matrix(clr(odamaki_select_age_final_genus+0.00001))
odamaki_select_age_final_genus_clr <- as.data.frame(t(apply(odamaki_select_age_final_genus_clr,1,function(x)(x-min(x)))))
print("Genus:Aitchison")
dist_genus_aitchison <- as.matrix(vegdist(odamaki_select_age_final_genus_clr,method="euclidean"))
diag(dist_genus_aitchison) <- NA
print("Genus:Kendall")
dist_genus_kendall <- as.matrix(1-cor.fk(t(odamaki_select_age_final_genus_clr))/2)
diag(dist_genus_kendall) <- NA

print("Species:Bray")
dist_species_bray <- as.matrix(vegdist(odamaki_select_age_final_species,method="bray"))
diag(dist_species_bray) <- NA
print("Species:Jaccard")
dist_species_jaccard <- as.matrix(vegdist(odamaki_select_age_final_species,method="jaccard"))
diag(dist_species_jaccard) <- NA

odamaki_select_age_final_species_clr <- as.matrix(clr(odamaki_select_age_final_species+0.00001))
odamaki_select_age_final_species_clr <- as.data.frame(t(apply(odamaki_select_age_final_species_clr,1,function(x)(x-min(x)))))
print("Species:Aitchison")
dist_species_aitchison <- as.matrix(vegdist(odamaki_select_age_final_species_clr,method="euclidean"))
diag(dist_species_aitchison) <- NA
print("Species:Spearman")
dist_species_kendall <- as.matrix(1-cor.fk(t(odamaki_select_age_final_species_clr))/2)
diag(dist_species_kendall) <- NA

print("Computing Uniqueness Measures")

df_odamaki_diversity_uniqueness$genus_bray_uniqueness <- apply(dist_genus_bray[rownames(df_odamaki_diversity_uniqueness),rownames(df_odamaki_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_odamaki_diversity_uniqueness$genus_jaccard_uniqueness <- apply(dist_genus_jaccard[rownames(df_odamaki_diversity_uniqueness),rownames(df_odamaki_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_odamaki_diversity_uniqueness$genus_aitchison_uniqueness <- apply(dist_genus_aitchison[rownames(df_odamaki_diversity_uniqueness),rownames(df_odamaki_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_odamaki_diversity_uniqueness$genus_kendall_uniqueness <- apply(dist_genus_kendall[rownames(df_odamaki_diversity_uniqueness),rownames(df_odamaki_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_odamaki_diversity_uniqueness$species_bray_uniqueness <- apply(dist_species_bray[rownames(df_odamaki_diversity_uniqueness),rownames(df_odamaki_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_odamaki_diversity_uniqueness$species_jaccard_uniqueness <- apply(dist_species_jaccard[rownames(df_odamaki_diversity_uniqueness),rownames(df_odamaki_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_odamaki_diversity_uniqueness$species_aitchison_uniqueness <- apply(dist_species_aitchison[rownames(df_odamaki_diversity_uniqueness),rownames(df_odamaki_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_odamaki_diversity_uniqueness$species_kendall_uniqueness <- apply(dist_species_kendall[rownames(df_odamaki_diversity_uniqueness),rownames(df_odamaki_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\Odamaki\\odamaki_analysis_2021_Revision.RData")
rm(list=ls())







