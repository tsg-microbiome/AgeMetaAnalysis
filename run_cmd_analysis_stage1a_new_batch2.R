# The sub-pipeline works on the curatedMetagenomicData3 data repository and performs two major tasks:
# First, it computes Shannon Diversity and the four different taxonomic uniqueness at the Species and Genus levels across all gut microbiomes in curatedMetagenomicData3 that are used for the disease-microbiome studies. The output is a single data frame with the select metadata and the diversity and uniqueness measures for all gut microbiomes. The data frame covers those studies in the curatedMetagenomicData3 repository both in Table 1A and Table 1B. The uniqueness measures are computed separately for each study cohort within curatedMetagenomicData3.
# Second, it collates the Metacyc-annotated pathway and metadata tables for the studies in the curatedMetagenomicData3 repository. It also omputes Shannon Diversity and the four different taxonomic uniqueness at the Pathway levels across all gut microbiomes. The output is a single data frame with the select metadata and the diversity and uniqueness measures for all gut microbiomes.

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

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_disease_analysis.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")

cmd3_unique_disease_study_names <- rownames(cmd3_select_disease_studies_details)
cmd3_unique_disease_study_names <- setdiff(cmd3_unique_disease_study_names,c("HansenLBS_2018"))
df_cmd3_disease_diversity_uniqueness <- data.frame(species_shannon = diversity(cmd3_select_disease_final_species), genus_shannon = diversity(cmd3_select_disease_final_species), row.names=rownames(cmd3_select_disease_final_metadata))

df_cmd3_disease_diversity_uniqueness$age <- cmd3_select_disease_final_metadata[rownames(df_cmd3_disease_diversity_uniqueness),"age"]
df_cmd3_disease_diversity_uniqueness$country <- cmd3_select_disease_final_metadata[rownames(df_cmd3_disease_diversity_uniqueness),"country"]
df_cmd3_disease_diversity_uniqueness$study_name <- cmd3_select_disease_final_metadata[rownames(df_cmd3_disease_diversity_uniqueness),"study_name"]
df_cmd3_disease_diversity_uniqueness$study_condition <- cmd3_select_disease_final_metadata[rownames(df_cmd3_disease_diversity_uniqueness),"study_condition"]
df_cmd3_disease_diversity_uniqueness$gender <- cmd3_select_disease_final_metadata[rownames(df_cmd3_disease_diversity_uniqueness),"gender"]
df_cmd3_disease_diversity_uniqueness$BMI <- cmd3_select_disease_final_metadata[rownames(df_cmd3_disease_diversity_uniqueness),"BMI"]

select_species <- names(which(apply(cmd3_select_disease_final_species,2,function(x)(length(x[x==0])))<=7500))
select_genus <- names(which(apply(cmd3_select_disease_final_genus,2,function(x)(length(x[x==0])))<=7500))

cmd3_select_disease_final_genus_clr <- as.matrix(clr(cmd3_select_disease_final_genus+0.00001))
cmd3_select_disease_final_genus_clr <- as.data.frame(t(apply(cmd3_select_disease_final_genus_clr,1,function(x)(x-min(x)))))
cmd3_select_disease_final_species_clr <- as.matrix(clr(cmd3_select_disease_final_species+0.00001))
cmd3_select_disease_final_species_clr <- as.data.frame(t(apply(cmd3_select_disease_final_species_clr,1,function(x)(x-min(x)))))

#cmd3_select_disease_final_genus <- cmd3_select_disease_final_genus[,select_genus]
#cmd3_select_disease_final_species <- cmd3_select_disease_final_species[,select_species]
#cmd3_select_disease_final_genus <- as.matrix(clr(cmd3_select_disease_final_genus))
#save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")

df_cmd3_disease_diversity_uniqueness$genus_bray_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$genus_jaccard_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$genus_aitchison_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$genus_kendall_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$species_bray_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$species_jaccard_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$species_aitchison_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$species_kendall_uniqueness <- NA

df_cmd3_disease_diversity_uniqueness$ac_genus_bray_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$ac_genus_jaccard_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$ac_genus_aitchison_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$ac_genus_kendall_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$ac_species_bray_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$ac_species_jaccard_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$ac_species_aitchison_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$ac_species_kendall_uniqueness <- NA

print("Computing Distance Matrices and Uniqueness Measures")
for(i in 1:length(cmd3_unique_disease_study_names))
{
	study_name <- cmd3_unique_disease_study_names[i]
	study_samples <- rownames(cmd3_select_disease_final_metadata[cmd3_select_disease_final_metadata$study_name == study_name,])
	
	print(paste0(study_name,",Genus:Bray"))
	dist_genus_bray <- as.matrix(vegdist(cmd3_select_disease_final_genus[study_samples,],method="bray"))
	diag(dist_genus_bray) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"genus_bray_uniqueness"] <- apply(dist_genus_bray,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_genus_bray_uniqueness"] <- rlm(genus_bray_uniqueness~genus_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Genus:Jaccard"))
	dist_genus_jaccard <- as.matrix(vegdist(cmd3_select_disease_final_genus[study_samples,],method="jaccard"))
	diag(dist_genus_jaccard) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"genus_jaccard_uniqueness"] <- apply(dist_genus_jaccard,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_genus_jaccard_uniqueness"] <- rlm(genus_jaccard_uniqueness~genus_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Genus:Aitchison"))
	dist_genus_aitchison <- as.matrix(vegdist(cmd3_select_disease_final_genus_clr[study_samples,],method="euclidean"))
	diag(dist_genus_aitchison) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"genus_aitchison_uniqueness"] <- apply(dist_genus_aitchison,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_genus_aitchison_uniqueness"] <- rlm(genus_aitchison_uniqueness~genus_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Genus:Kendall"))
	dist_genus_kendall <- as.matrix(1-cor.fk(t(cmd3_select_disease_final_genus_clr[study_samples,]))/2)
	print("M here")
	diag(dist_genus_kendall) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"genus_kendall_uniqueness"] <- apply(dist_genus_kendall,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"genus_kendall_uniqueness"] <- ifelse(is.infinite(df_cmd3_disease_diversity_uniqueness[study_samples,"genus_kendall_uniqueness"]),1,df_cmd3_disease_diversity_uniqueness[study_samples,"genus_kendall_uniqueness"])
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_genus_kendall_uniqueness"] <- rlm(genus_kendall_uniqueness~genus_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Species:Bray"))
	dist_species_bray <- as.matrix(vegdist(cmd3_select_disease_final_species[study_samples,],method="bray"))
	diag(dist_species_bray) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"species_bray_uniqueness"] <- apply(dist_species_bray,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_species_bray_uniqueness"] <- rlm(species_bray_uniqueness~species_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Species:Jaccard"))
	dist_species_jaccard <- as.matrix(vegdist(cmd3_select_disease_final_species[study_samples,],method="jaccard"))
	diag(dist_species_jaccard) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"species_jaccard_uniqueness"] <- apply(dist_species_jaccard,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_species_jaccard_uniqueness"] <- rlm(species_jaccard_uniqueness~species_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Species:Aitchison"))
	dist_species_aitchison <- as.matrix(vegdist(cmd3_select_disease_final_species_clr[study_samples,],method="euclidean"))
	diag(dist_species_aitchison) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"species_aitchison_uniqueness"] <- apply(dist_species_aitchison,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_species_aitchison_uniqueness"] <- rlm(species_aitchison_uniqueness~species_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Species:Kendall"))
	dist_species_kendall <- as.matrix(1-cor.fk(t(cmd3_select_disease_final_species_clr[study_samples,]))/2)
	diag(dist_species_kendall) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"species_kendall_uniqueness"] <- apply(dist_species_kendall,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"species_kendall_uniqueness"] <- ifelse(is.infinite(df_cmd3_disease_diversity_uniqueness[study_samples,"species_kendall_uniqueness"]),1,df_cmd3_disease_diversity_uniqueness[study_samples,"species_kendall_uniqueness"])
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_species_kendall_uniqueness"] <- rlm(species_kendall_uniqueness~species_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals
	
}
print("Saving Workspace")

df_cmd3_disease_diversity_uniqueness_old <- df_cmd3_disease_diversity_uniqueness

#source("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\run_cmd_analysis_stage1a_pathway.R")
save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_pathway_analysis.RData")
cmd3_select_disease_final_pathway <- apply(cmd3_stool_disease_pathway_profile[rownames(cmd3_select_disease_final_metadata),],2,as.numeric)
rownames(cmd3_select_disease_final_pathway) <- rownames(cmd3_select_disease_final_metadata)

cmd3_select_disease_final_pathway_clr <- as.matrix(clr(cmd3_select_disease_final_pathway+0.00001))
cmd3_select_disease_final_pathway_clr <- as.data.frame(t(apply(cmd3_select_disease_final_pathway_clr,1,function(x)(x-min(x)))))


df_cmd3_disease_diversity_uniqueness$pathway_shannon <- diversity(cmd3_select_disease_final_pathway)[rownames(df_cmd3_disease_diversity_uniqueness)]

df_cmd3_disease_diversity_uniqueness$pathway_bray_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$pathway_jaccard_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$pathway_aitchison_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$pathway_kendall_uniqueness <- NA

df_cmd3_disease_diversity_uniqueness$ac_pathway_bray_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$ac_pathway_jaccard_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$ac_pathway_aitchison_uniqueness <- NA
df_cmd3_disease_diversity_uniqueness$ac_pathway_kendall_uniqueness <- NA

print("Computing Pathway Uniqueness Measures")

for(i in 1:length(cmd3_unique_disease_study_names))
{
	study_name <- cmd3_unique_disease_study_names[i]
	study_samples <- rownames(cmd3_select_disease_final_metadata[cmd3_select_disease_final_metadata$study_name == study_name,])
	
	print(paste0(study_name,",Pathway:Bray"))
	dist_pathway_bray <- as.matrix(vegdist(cmd3_select_disease_final_pathway[study_samples,],method="bray"))
	diag(dist_pathway_bray) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"pathway_bray_uniqueness"] <- apply(dist_pathway_bray,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_pathway_bray_uniqueness"] <- rlm(pathway_bray_uniqueness~pathway_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Pathway:Jaccard"))
	dist_pathway_jaccard <- as.matrix(vegdist(cmd3_select_disease_final_pathway[study_samples,],method="jaccard"))
	diag(dist_pathway_jaccard) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"pathway_jaccard_uniqueness"] <- apply(dist_pathway_jaccard,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_pathway_jaccard_uniqueness"] <- rlm(pathway_jaccard_uniqueness~pathway_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Pathway:Aitchison"))
	dist_pathway_aitchison <- as.matrix(vegdist(cmd3_select_disease_final_pathway_clr[study_samples,],method="euclidean"))
	diag(dist_pathway_aitchison) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"pathway_aitchison_uniqueness"] <- apply(dist_pathway_aitchison,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_pathway_aitchison_uniqueness"] <- rlm(pathway_aitchison_uniqueness~pathway_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals
	
	print(paste0(study_name,",Pathway:Kendall"))
	dist_pathway_kendall <- as.matrix(1-cor.fk(t(cmd3_select_disease_final_pathway_clr[study_samples,]))/2)
	diag(dist_pathway_kendall) <- NA
	df_cmd3_disease_diversity_uniqueness[study_samples,"pathway_kendall_uniqueness"] <- apply(dist_pathway_kendall,1,function(x)(min(x[!is.na(x)])))
	df_cmd3_disease_diversity_uniqueness[study_samples,"pathway_kendall_uniqueness"] <- ifelse(is.infinite(df_cmd3_disease_diversity_uniqueness[study_samples,"pathway_kendall_uniqueness"]),1,df_cmd3_disease_diversity_uniqueness[study_samples,"pathway_kendall_uniqueness"])
	df_cmd3_disease_diversity_uniqueness[study_samples,"ac_pathway_kendall_uniqueness"] <- rlm(pathway_kendall_uniqueness~pathway_shannon,data=df_cmd3_disease_diversity_uniqueness[study_samples,])$residuals

}

df_cmd3_diversity_uniqueness <- as.data.frame(rbind(df_cmd3_diversity_uniqueness,df_cmd3_disease_diversity_uniqueness[df_cmd3_disease_diversity_uniqueness$study_name == "QinJ_2012",colnames(df_cmd3_diversity_uniqueness)]))

QinJ_2012_Rows <- rownames(cmd3_select_disease_final_metadata[cmd3_select_disease_final_metadata$study_name == "QinJ_2012",])

#if(length(intersect(rownames(cmd3_select_age_final_species),QinJ_2012_Rows))>0)
#{
	temp0 <- as.data.frame(merge(t(cmd3_select_age_final_species),t(cmd3_select_disease_final_species[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,-1])
	rownames(temp0) <- merge(t(cmd3_select_age_final_species),t(cmd3_select_disease_final_species[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,1]
	temp0 <- apply(temp0,1,function(x)(ifelse(is.na(x),0,x)))
	cmd3_select_age_final_species <- as.data.frame(temp0)
	
	temp0 <- as.data.frame(merge(t(cmd3_select_age_final_species_clr),t(cmd3_select_disease_final_species_clr[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,-1])
	rownames(temp0) <- merge(t(cmd3_select_age_final_species_clr),t(cmd3_select_disease_final_species_clr[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,1]
	temp0 <- apply(temp0,1,function(x)(ifelse(is.na(x),0,x)))
	cmd3_select_age_final_species_clr <- as.data.frame(temp0)

	temp0 <- as.data.frame(merge(t(cmd3_select_age_final_genus),t(cmd3_select_disease_final_genus[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,-1])
	rownames(temp0) <- merge(t(cmd3_select_age_final_genus),t(cmd3_select_disease_final_genus[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,1]
	temp0 <- apply(temp0,1,function(x)(ifelse(is.na(x),0,x)))
	cmd3_select_age_final_genus <- as.data.frame(temp0)
	
	temp0 <- as.data.frame(merge(t(cmd3_select_age_final_genus_clr),t(cmd3_select_disease_final_genus_clr[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,-1])
	rownames(temp0) <- merge(t(cmd3_select_age_final_genus_clr),t(cmd3_select_disease_final_genus_clr[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,1]
	temp0 <- apply(temp0,1,function(x)(ifelse(is.na(x),0,x)))
	cmd3_select_age_final_genus_clr <- as.data.frame(temp0)

	temp0 <- as.data.frame(merge(t(cmd3_select_age_final_pathway),t(cmd3_select_disease_final_pathway[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,-1])
	rownames(temp0) <- merge(t(cmd3_select_age_final_pathway),t(cmd3_select_disease_final_pathway[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,1]
	temp0 <- apply(temp0,1,function(x)(ifelse(is.na(x),0,x)))
	cmd3_select_age_final_pathway <- as.data.frame(temp0)
	
	temp0 <- as.data.frame(merge(t(cmd3_select_age_final_pathway_clr),t(cmd3_select_disease_final_pathway_clr[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,-1])
	rownames(temp0) <- merge(t(cmd3_select_age_final_pathway_clr),t(cmd3_select_disease_final_pathway_clr[QinJ_2012_Rows,]),by="row.names",all=TRUE)[,1]
	temp0 <- apply(temp0,1,function(x)(ifelse(is.na(x),0,x)))
	cmd3_select_age_final_pathway_clr <- as.data.frame(temp0)

	cmd3_select_age_final_metadata <- as.data.frame(rbind(cmd3_select_age_final_metadata,cmd3_select_disease_final_metadata[QinJ_2012_Rows,]))
#}

df_cmd3_disease_diversity_uniqueness_old <- df_cmd3_disease_diversity_uniqueness

print("Saving Workspace")
#source("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\run_cmd_analysis_stage1a_pathway.R")
save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")
rm(list=ls())


