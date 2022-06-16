# The sub-pipeline scans/inputs the species, genus and metadata tables for the studies included in the age-microbiome analyses for the curatedMetagenomicData3 data repository.
# It computes Shannon Diversity and the four different taxonomic uniqueness at the Species and Genus levels across all gut microbiomes
# The output is a single data frame with the select metadata and the diversity and uniqueness measures for all gut microbiomes
# The data frame covers only those studies in the curatedMetagenomicData3 repository in Table 1A

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

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_age_analysis.RData")

cmd3_unique_study_names <- cmd3_final_set_sample_details[,1]
df_cmd3_diversity_uniqueness <- data.frame(species_shannon = diversity(cmd3_select_age_final_species), genus_shannon = diversity(cmd3_select_age_final_genus), row.names=rownames(cmd3_select_age_final_species))

df_cmd3_diversity_uniqueness$age <- cmd3_select_age_final_metadata[rownames(df_cmd3_diversity_uniqueness),"age"]
df_cmd3_diversity_uniqueness$country <- cmd3_select_age_final_metadata[rownames(df_cmd3_diversity_uniqueness),"country"]
df_cmd3_diversity_uniqueness$study_name <- cmd3_select_age_final_metadata[rownames(df_cmd3_diversity_uniqueness),"study_name"]
df_cmd3_diversity_uniqueness$study_condition <- cmd3_select_age_final_metadata[rownames(df_cmd3_diversity_uniqueness),"study_condition"]
df_cmd3_diversity_uniqueness$gender <- cmd3_select_age_final_metadata[rownames(df_cmd3_diversity_uniqueness),"gender"]
df_cmd3_diversity_uniqueness$BMI <- cmd3_select_age_final_metadata[rownames(df_cmd3_diversity_uniqueness),"BMI"]

select_species <- names(which(apply(cmd3_select_age_final_species,2,function(x)(length(x[x==0])))<=7500))
select_genus <- names(which(apply(cmd3_select_age_final_genus,2,function(x)(length(x[x==0])))<=7500))

#cmd3_select_age_final_genus <- cmd3_select_age_final_genus[,select_genus]
#cmd3_select_age_final_species <- cmd3_select_age_final_species[,select_species]
#cmd3_select_age_final_genus <- as.matrix(clr(cmd3_select_age_final_genus))
#save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")

cmd3_select_age_final_genus_clr <- as.matrix(clr(cmd3_select_age_final_genus+0.00001))
cmd3_select_age_final_genus_clr <- as.data.frame(t(apply(cmd3_select_age_final_genus_clr,1,function(x)(x-min(x)))))
cmd3_select_age_final_species_clr <- as.matrix(clr(cmd3_select_age_final_species+0.00001))
cmd3_select_age_final_species_clr <- as.data.frame(t(apply(cmd3_select_age_final_species_clr,1,function(x)(x-min(x)))))

df_cmd3_diversity_uniqueness$genus_bray_uniqueness <- NA
df_cmd3_diversity_uniqueness$genus_jaccard_uniqueness <- NA
df_cmd3_diversity_uniqueness$genus_aitchison_uniqueness <- NA
df_cmd3_diversity_uniqueness$genus_kendall_uniqueness <- NA
df_cmd3_diversity_uniqueness$species_bray_uniqueness <- NA
df_cmd3_diversity_uniqueness$species_jaccard_uniqueness <- NA
df_cmd3_diversity_uniqueness$species_aitchison_uniqueness <- NA
df_cmd3_diversity_uniqueness$species_kendall_uniqueness <- NA

df_cmd3_diversity_uniqueness$ac_genus_bray_uniqueness <- NA
df_cmd3_diversity_uniqueness$ac_genus_jaccard_uniqueness <- NA
df_cmd3_diversity_uniqueness$ac_genus_aitchison_uniqueness <- NA
df_cmd3_diversity_uniqueness$ac_genus_kendall_uniqueness <- NA
df_cmd3_diversity_uniqueness$ac_species_bray_uniqueness <- NA
df_cmd3_diversity_uniqueness$ac_species_jaccard_uniqueness <- NA
df_cmd3_diversity_uniqueness$ac_species_aitchison_uniqueness <- NA
df_cmd3_diversity_uniqueness$ac_species_kendall_uniqueness <- NA

print("Computing Distance Matrices and Uniqueness Measures")
for(i in 1:length(cmd3_unique_study_names))
{
	study_name <- cmd3_unique_study_names[i]
	study_samples <- rownames(cmd3_select_age_final_metadata[cmd3_select_age_final_metadata$study_name == study_name,])
	
	print(paste0(study_name,",Genus:Bray"))
	dist_genus_bray <- as.matrix(vegdist(cmd3_select_age_final_genus[study_samples,],method="bray"))
	diag(dist_genus_bray) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"genus_bray_uniqueness"] <- apply(dist_genus_bray,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"ac_genus_bray_uniqueness"] <- rlm(genus_bray_uniqueness~genus_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Genus:Jaccard"))
	dist_genus_jaccard <- as.matrix(vegdist(cmd3_select_age_final_genus[study_samples,],method="jaccard"))
	diag(dist_genus_jaccard) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"genus_jaccard_uniqueness"] <- apply(dist_genus_jaccard,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"ac_genus_jaccard_uniqueness"] <- rlm(genus_jaccard_uniqueness~genus_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Genus:Aitchison"))
	dist_genus_aitchison <- as.matrix(vegdist(cmd3_select_age_final_genus_clr[study_samples,],method="euclidean"))
	diag(dist_genus_aitchison) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"genus_aitchison_uniqueness"] <- apply(dist_genus_aitchison,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"ac_genus_aitchison_uniqueness"] <- rlm(genus_aitchison_uniqueness~genus_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Genus:Kendall"))
	dist_genus_kendall <- as.matrix(1-cor.fk(t(cmd3_select_age_final_genus_clr[study_samples,]))/2)
	diag(dist_genus_kendall) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"genus_kendall_uniqueness"] <- apply(dist_genus_kendall,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"genus_kendall_uniqueness"] <- ifelse(is.infinite(df_cmd3_diversity_uniqueness[study_samples,"genus_kendall_uniqueness"]),1,df_cmd3_diversity_uniqueness[study_samples,"genus_kendall_uniqueness"])
	df_cmd3_diversity_uniqueness[study_samples,"ac_genus_kendall_uniqueness"] <- rlm(genus_kendall_uniqueness~genus_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Species:Bray"))
	dist_species_bray <- as.matrix(vegdist(cmd3_select_age_final_species[study_samples,],method="bray"))
	diag(dist_species_bray) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"species_bray_uniqueness"] <- apply(dist_species_bray,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"ac_species_bray_uniqueness"] <- rlm(species_bray_uniqueness~species_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Species:Jaccard"))
	dist_species_jaccard <- as.matrix(vegdist(cmd3_select_age_final_species[study_samples,],method="jaccard"))
	diag(dist_species_jaccard) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"species_jaccard_uniqueness"] <- apply(dist_species_jaccard,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"ac_species_jaccard_uniqueness"] <- rlm(species_jaccard_uniqueness~species_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Species:Aitchison"))
	dist_species_aitchison <- as.matrix(vegdist(cmd3_select_age_final_species_clr[study_samples,],method="euclidean"))
	diag(dist_species_aitchison) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"species_aitchison_uniqueness"] <- apply(dist_species_aitchison,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"ac_species_aitchison_uniqueness"] <- rlm(species_aitchison_uniqueness~species_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Species:Kendall"))
	dist_species_kendall <- as.matrix(1-cor.fk(t(cmd3_select_age_final_species_clr[study_samples,]))/2)
	diag(dist_species_kendall) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"species_kendall_uniqueness"] <- apply(dist_species_kendall,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"species_kendall_uniqueness"] <- ifelse(is.infinite(df_cmd3_diversity_uniqueness[study_samples,"species_kendall_uniqueness"]),1,df_cmd3_diversity_uniqueness[study_samples,"species_kendall_uniqueness"])
	df_cmd3_diversity_uniqueness[study_samples,"ac_species_kendall_uniqueness"] <- rlm(species_kendall_uniqueness~species_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals
	
}
print("Saving Workspace")
#source("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\run_cmd_analysis_stage1a_pathway.R")
save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_pathway_analysis.RData")
cmd3_select_age_final_pathway <- apply(cmd3_stool_pathway_profile[rownames(cmd3_select_age_final_metadata),],2,as.numeric)
rownames(cmd3_select_age_final_pathway) <- rownames(cmd3_select_age_final_metadata)

cmd3_select_age_final_pathway_clr <- as.matrix(clr(cmd3_select_age_final_pathway+0.00001))
cmd3_select_age_final_pathway_clr <- as.data.frame(t(apply(cmd3_select_age_final_pathway_clr,1,function(x)(x-min(x)))))

df_cmd3_diversity_uniqueness$pathway_shannon <- diversity(cmd3_select_age_final_pathway)[rownames(df_cmd3_diversity_uniqueness)]

df_cmd3_diversity_uniqueness$pathway_bray_uniqueness <- NA
df_cmd3_diversity_uniqueness$pathway_jaccard_uniqueness <- NA
df_cmd3_diversity_uniqueness$pathway_aitchison_uniqueness <- NA
df_cmd3_diversity_uniqueness$pathway_kendall_uniqueness <- NA

df_cmd3_diversity_uniqueness$ac_pathway_bray_uniqueness <- NA
df_cmd3_diversity_uniqueness$ac_pathway_jaccard_uniqueness <- NA
df_cmd3_diversity_uniqueness$ac_pathway_aitchison_uniqueness <- NA
df_cmd3_diversity_uniqueness$ac_pathway_kendall_uniqueness <- NA

print("Computing Pathway Uniqueness Measures")

for(i in 1:length(cmd3_unique_study_names))
{
	study_name <- cmd3_unique_study_names[i]
	study_samples <- rownames(cmd3_select_age_final_metadata[cmd3_select_age_final_metadata$study_name == study_name,])
	
	print(paste0(study_name,",Pathway:Bray"))
	dist_pathway_bray <- as.matrix(vegdist(cmd3_select_age_final_pathway[study_samples,],method="bray"))
	diag(dist_pathway_bray) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"pathway_bray_uniqueness"] <- apply(dist_pathway_bray,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"ac_pathway_bray_uniqueness"] <- rlm(pathway_bray_uniqueness~pathway_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Pathway:Jaccard"))
	dist_pathway_jaccard <- as.matrix(vegdist(cmd3_select_age_final_pathway[study_samples,],method="jaccard"))
	diag(dist_pathway_jaccard) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"pathway_jaccard_uniqueness"] <- apply(dist_pathway_jaccard,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"ac_pathway_jaccard_uniqueness"] <- rlm(pathway_jaccard_uniqueness~pathway_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Pathway:Aitchison"))
	dist_pathway_aitchison <- as.matrix(vegdist(cmd3_select_age_final_pathway_clr[study_samples,],method="euclidean"))
	diag(dist_pathway_aitchison) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"pathway_aitchison_uniqueness"] <- apply(dist_pathway_aitchison,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"ac_pathway_aitchison_uniqueness"] <- rlm(pathway_aitchison_uniqueness~pathway_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Pathway:Kendall"))
	dist_pathway_kendall <- as.matrix(1-cor.fk(t(cmd3_select_age_final_pathway[study_samples,]))/2)
	diag(dist_pathway_kendall) <- NA
	df_cmd3_diversity_uniqueness[study_samples,"pathway_kendall_uniqueness"] <- apply(dist_pathway_kendall,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_diversity_uniqueness[study_samples,"pathway_kendall_uniqueness"] <- ifelse(is.infinite(df_cmd3_diversity_uniqueness[study_samples,"pathway_kendall_uniqueness"]),1,df_cmd3_diversity_uniqueness[study_samples,"pathway_kendall_uniqueness"])
	df_cmd3_diversity_uniqueness[study_samples,"ac_pathway_kendall_uniqueness"] <- rlm(pathway_kendall_uniqueness~pathway_shannon,data=df_cmd3_diversity_uniqueness[study_samples,])$residuals

}

print("Saving Workspace")
#source("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\run_cmd_analysis_stage1a_pathway.R")
save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")
rm(list=ls())


