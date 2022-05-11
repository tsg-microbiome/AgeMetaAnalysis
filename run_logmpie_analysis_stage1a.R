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

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\LogMPie\\logmpie_age_analysis.RData")

filter_rows <- names(which(!is.nan(rowSums(logmpie_select_age_final_genus))))

logmpie_select_age_final_species <- logmpie_select_age_final_species[filter_rows,]
logmpie_select_age_final_genus <- logmpie_select_age_final_genus[filter_rows,]
logmpie_select_age_final_metadata <- logmpie_select_age_final_metadata[filter_rows,]


df_logmpie_diversity_uniqueness <- data.frame(species_shannon = diversity(logmpie_select_age_final_species), genus_shannon = diversity(logmpie_select_age_final_genus), row.names=rownames(logmpie_select_age_final_species))
df_logmpie_diversity_uniqueness$age <- logmpie_select_age_final_metadata[rownames(df_logmpie_diversity_uniqueness),"Age"]
df_logmpie_diversity_uniqueness$country <- logmpie_select_age_final_metadata[rownames(df_logmpie_diversity_uniqueness),"Districts"]
df_logmpie_diversity_uniqueness$study_name <- "LogMPie"
df_logmpie_diversity_uniqueness$BMI <- logmpie_select_age_final_metadata[rownames(df_logmpie_diversity_uniqueness),"BMI"]

print("Genus:Bray")
dist_genus_bray <- as.matrix(vegdist(logmpie_select_age_final_genus,method="bray"))
diag(dist_genus_bray) <- NA
print("Genus:Jaccard")
dist_genus_jaccard <- as.matrix(vegdist(logmpie_select_age_final_genus,method="jaccard"))
diag(dist_genus_jaccard) <- NA

logmpie_select_age_final_genus_clr <- as.matrix(clr(logmpie_select_age_final_genus+0.00001))
logmpie_select_age_final_genus_clr <- as.data.frame(t(apply(logmpie_select_age_final_genus_clr,1,function(x)(x-min(x)))))
print("Genus:Aichison")
dist_genus_aitchison <- as.matrix(vegdist(logmpie_select_age_final_genus_clr,method="euclidean"))
diag(dist_genus_aitchison) <- NA
print("Genus:Spearman")
dist_genus_kendall <- as.matrix(1-cor.fk(t(logmpie_select_age_final_genus_clr))/2)
diag(dist_genus_kendall) <- NA

print("Species:Bray")
dist_species_bray <- as.matrix(vegdist(logmpie_select_age_final_species,method="bray"))
diag(dist_species_bray) <- NA
print("Species:Jaccard")
dist_species_jaccard <- as.matrix(vegdist(logmpie_select_age_final_species,method="jaccard"))
diag(dist_species_jaccard) <- NA

logmpie_select_age_final_species_clr <- as.matrix(clr(logmpie_select_age_final_species+0.00001))
logmpie_select_age_final_species_clr <- as.data.frame(t(apply(logmpie_select_age_final_species_clr,1,function(x)(x-min(x)))))
print("Species:Manhattan")
dist_species_aitchison <- as.matrix(vegdist(logmpie_select_age_final_species_clr,method="euclidean"))
diag(dist_species_aitchison) <- NA
print("Species:Kendall")
dist_species_kendall <- as.matrix(1-cor.fk(t(logmpie_select_age_final_species_clr))/2)
diag(dist_species_kendall) <- NA

print("Computing Uniqueness Measures")

df_logmpie_diversity_uniqueness$genus_bray_uniqueness <- apply(dist_genus_bray[rownames(df_logmpie_diversity_uniqueness),rownames(df_logmpie_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_logmpie_diversity_uniqueness$genus_jaccard_uniqueness <- apply(dist_genus_jaccard[rownames(df_logmpie_diversity_uniqueness),rownames(df_logmpie_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_logmpie_diversity_uniqueness$genus_aitchison_uniqueness <- apply(dist_genus_aitchison[rownames(df_logmpie_diversity_uniqueness),rownames(df_logmpie_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_logmpie_diversity_uniqueness$genus_kendall_uniqueness <- apply(dist_genus_kendall[rownames(df_logmpie_diversity_uniqueness),rownames(df_logmpie_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_logmpie_diversity_uniqueness$species_bray_uniqueness <- apply(dist_species_bray[rownames(df_logmpie_diversity_uniqueness),rownames(df_logmpie_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_logmpie_diversity_uniqueness$species_jaccard_uniqueness <- apply(dist_species_jaccard[rownames(df_logmpie_diversity_uniqueness),rownames(df_logmpie_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_logmpie_diversity_uniqueness$species_aitchison_uniqueness <- apply(dist_species_aitchison[rownames(df_logmpie_diversity_uniqueness),rownames(df_logmpie_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))
df_logmpie_diversity_uniqueness$species_kendall_uniqueness <- apply(dist_species_kendall[rownames(df_logmpie_diversity_uniqueness),rownames(df_logmpie_diversity_uniqueness)],1,function(x)(min(x[!is.na(x)])))

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\LogMPie\\logmpie_analysis_2021_Revision.RData")

rm(list=ls())







