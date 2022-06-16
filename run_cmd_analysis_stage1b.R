# This sub-pipeline works on the curatedMetagenomicData3 and computes the following:
# Associations between Shannon Diversity and the four different measures of Uniqueness at the Species and Genus level.
# Associations between gut microbiome beta-diversity (computed using the four different distance measures; the same ones corresponding to the four measures of uniqueness) and age.
# Associations between the different measures of uniqueness (four measures of uniqueness each corresponding to the species and genus levels) and Shannon Diversity with age.
# All associations are computed on a per-study basis.
library(MASS)
library(sfsmisc)
library(gplots)
library(RColorBrewer)
library(pcaPP)
library(vegan)
library(compositions)

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")
#cmd3_select_age_final_genus_clr <- as.matrix(clr(cmd3_select_age_final_genus+0.00001))
#cmd3_select_age_final_species_clr <- as.matrix(clr(cmd3_select_age_final_species+0.00001))
#cmd3_select_age_final_pathway_clr <- as.matrix(clr(cmd3_select_age_final_pathway+0.00001))

cmd3_all_controls <- rownames(cmd3_select_age_final_metadata[cmd3_select_age_final_metadata$study_condition=="control",])
cmd3_study_with_controls <- unique(cmd3_select_age_final_metadata[cmd3_all_controls,"study_name"])
df_cmd3_controls_diversity_uniqueness <- df_cmd3_diversity_uniqueness[df_cmd3_diversity_uniqueness$study_condition=="control",]

cmd3_select_studies_controls <- as.data.frame(matrix(NA,length(cmd3_study_with_controls),6))
rownames(cmd3_select_studies_controls) <- cmd3_study_with_controls
colnames(cmd3_select_studies_controls) <- c("Number","Country","Minimum_Age","Maximum_Age","Percentage_60","Study_Code")
for(i in 1:length(cmd3_study_with_controls))
{
	study_name <- cmd3_study_with_controls[i]
	cmd3_select_studies_controls[i,1] <- nrow(df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name == study_name,])
	cmd3_select_studies_controls[i,2] <- paste0(unique(df_cmd3_controls_diversity_uniqueness[(df_cmd3_controls_diversity_uniqueness$study_name == study_name),"country"]),collapse=",")
	cmd3_select_studies_controls[i,3] <- min(df_cmd3_controls_diversity_uniqueness[(df_cmd3_controls_diversity_uniqueness$study_name == study_name),"age"])
	cmd3_select_studies_controls[i,4] <- max(df_cmd3_controls_diversity_uniqueness[(df_cmd3_controls_diversity_uniqueness$study_name == study_name),"age"])
	cmd3_select_studies_controls[i,5] <- round(nrow(df_cmd3_controls_diversity_uniqueness[(df_cmd3_controls_diversity_uniqueness$study_name == study_name)&(df_cmd3_controls_diversity_uniqueness$age >= 60),])/nrow(df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name == study_name,]),2)
	cmd3_select_studies_controls[i,6] <- as.integer(runif(1, min=0, max=50))
} 

dummy_map <- cmd3_select_studies_controls[,6]
names(dummy_map) <- rownames(cmd3_select_studies_controls)
df_cmd3_controls_diversity_uniqueness$study_dummy <- dummy_map[df_cmd3_controls_diversity_uniqueness$study_name]

cmd3_sorted_study_list_controls <- rownames(cmd3_select_studies_controls[rev(order(cmd3_select_studies_controls[,5])),])

cmd3_age_landscape_studies_controls <- rownames(cmd3_select_studies_controls[(cmd3_select_studies_controls$Minimum_Age<=30)&(cmd3_select_studies_controls$Maximum_Age>=60),])

EU_NA_Studies <- c("AsnicarF_2021","CosteaPI_2017","HansenLBS_2018","HMP_2019_ibdmdb","KeohaneDM_2020","NielsenHB_2014","SankaranarayananK_2015","SchirmerM_2016","WirbelJ_2018","ZellerG_2014")
Other_Studies <- c("DhakanDB_2019","GuptaA_2019","BritoIL_2016","LokmerA_2019","PasolliE_2019","PehrssonE_2016","RubelMA_2020","RampelliS_2015")
EA_Studies <- c("QinN_2014","QinJ_2012","YeZ_2018","YachidaS_2019")

sorted_cmd3_age_landscape_studies_controls <- c(EU_NA_Studies,EA_Studies,Other_Studies)

print(dim(cmd3_select_age_final_genus[QinJ_2012_Rows,]))

### START OF Step A ###
print("Step A")
cmd3_rlm_est_diversity_uniqueness <- as.data.frame(matrix(NA,length(cmd3_sorted_study_list_controls),12))
rownames(cmd3_rlm_est_diversity_uniqueness) <- cmd3_sorted_study_list_controls
colnames(cmd3_rlm_est_diversity_uniqueness) <- c("genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

cmd3_rlm_p_val_diversity_uniqueness <- as.data.frame(matrix(NA,length(cmd3_sorted_study_list_controls),12))
rownames(cmd3_rlm_p_val_diversity_uniqueness) <- cmd3_sorted_study_list_controls
colnames(cmd3_rlm_p_val_diversity_uniqueness) <- c("genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")


for(i in 1:length(cmd3_sorted_study_list_controls))
{
	study_name = cmd3_sorted_study_list_controls[i]
	
	temp_rlm <- rlm(genus_bray_uniqueness~genus_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,1] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,1] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(genus_jaccard_uniqueness~genus_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,2] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,2] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(genus_aitchison_uniqueness~genus_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,3] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,3] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(genus_kendall_uniqueness~genus_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,4] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,4] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(species_bray_uniqueness~species_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,5] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,5] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(species_jaccard_uniqueness~species_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,6] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,6] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(species_aitchison_uniqueness~species_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,7] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,7] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(species_kendall_uniqueness~species_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,8] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,8] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(pathway_bray_uniqueness~pathway_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,9] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,9] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(pathway_jaccard_uniqueness~pathway_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,10] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,10] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(pathway_aitchison_uniqueness~pathway_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,11] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,11] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(pathway_kendall_uniqueness~pathway_shannon,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_diversity_uniqueness[i,12] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_diversity_uniqueness[i,12] <- f.robftest(temp_rlm)$p.value
}

cmd3_rlm_q_val_diversity_uniqueness <- t(apply(cmd3_rlm_p_val_diversity_uniqueness,1,p.adjust))

cmd3_rlm_dir_diversity_uniqueness <- as.data.frame(matrix(NA,length(cmd3_sorted_study_list_controls),12))
rownames(cmd3_rlm_dir_diversity_uniqueness) <- cmd3_sorted_study_list_controls
colnames(cmd3_rlm_dir_diversity_uniqueness) <- c("genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

for(i in 1:nrow(cmd3_rlm_dir_diversity_uniqueness))
{
	for(j in 1:ncol(cmd3_rlm_dir_diversity_uniqueness))
	{
		cmd3_rlm_dir_diversity_uniqueness[i,j] <- ifelse(cmd3_rlm_q_val_diversity_uniqueness[i,j]<=0.05,3*sign(cmd3_rlm_est_diversity_uniqueness[i,j]),ifelse(cmd3_rlm_q_val_diversity_uniqueness[i,j]<=0.1,2*sign(cmd3_rlm_est_diversity_uniqueness[i,j]),sign(cmd3_rlm_est_diversity_uniqueness[i,j])))
	}
}

hmp_cmd3_diversity_uniqueness <- heatmap.2(t(cmd3_rlm_dir_diversity_uniqueness),density="none",trace="none",col=c("skyblue4","skyblue3","skyblue","orange","orangered","orangered4"))

#heatmap.2(t(cmd3_rlm_dir_diversity_uniqueness),density="none",trace="none",col=c("skyblue4","skyblue3","skyblue","orange","orangered","orangered4"),margins=c(15,10),lwid=c(0.5,5),lhei=c(1,5),Rowv=FALSE)

print(dim(cmd3_select_age_final_genus[QinJ_2012_Rows,]))

### END OF Step A ###

### START OF Step B ###
print("Step B")
cmd3_est_age_beta <- as.data.frame(matrix(0,length(sorted_cmd3_age_landscape_studies_controls),15))
rownames(cmd3_est_age_beta) <- sorted_cmd3_age_landscape_studies_controls
colnames(cmd3_est_age_beta) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

cmd3_p_val_age_beta <- as.data.frame(matrix(1,length(sorted_cmd3_age_landscape_studies_controls),15))
rownames(cmd3_p_val_age_beta) <- sorted_cmd3_age_landscape_studies_controls
colnames(cmd3_p_val_age_beta) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

for(i in 1:length(sorted_cmd3_age_landscape_studies_controls))
{
	study_name = sorted_cmd3_age_landscape_studies_controls[i]
	study_samples <- rownames(df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name == study_name,])
	print(paste0(study_name,",",length(study_samples)))
	
	cmd3_est_age_beta[i,1] <- 0
	cmd3_p_val_age_beta[i,1] <- 1
	
	cmd3_est_age_beta[i,2] <- 0
	cmd3_p_val_age_beta[i,2] <- 1
	
	print(paste0(study_name,",Genus:Bray"))
	temp_adonis <- adonis(vegdist(cmd3_select_age_final_genus[study_samples,],method="bray")~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,3] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,3] <- temp_adonis$aov.tab[1,6]
	
	print(paste0(study_name,",Genus:Jaccard"))
	temp_adonis <- adonis(vegdist(cmd3_select_age_final_genus[study_samples,],method="jaccard")~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,4] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,4] <- temp_adonis$aov.tab[1,6]
	
	print(paste0(study_name,",Genus:Aitchinson"))
	temp_adonis <- adonis(vegdist(cmd3_select_age_final_genus_clr[study_samples,],method="euclidean")~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,5] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,5] <- temp_adonis$aov.tab[1,6]
	
	print(paste0(study_name,",Genus:Kendall"))
	temp_adonis <- adonis(as.dist(1-cor.fk(t(cmd3_select_age_final_genus_clr[study_samples,]))/2)~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,6] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,6] <- temp_adonis$aov.tab[1,6]
	
	print(paste0(study_name,",Species:Bray"))
	temp_adonis <- adonis(vegdist(cmd3_select_age_final_species[study_samples,],method="bray")~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,7] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,7] <- temp_adonis$aov.tab[1,6]
	
	print(paste0(study_name,",Species:Jaccard"))
	temp_adonis <- adonis(vegdist(cmd3_select_age_final_species[study_samples,],method="jaccard")~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,8] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,8] <- temp_adonis$aov.tab[1,6]
	
	print(paste0(study_name,",Species:Aitchinson"))
	temp_adonis <- adonis(vegdist(cmd3_select_age_final_species_clr[study_samples,],method="euclidean")~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,9] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,9] <- temp_adonis$aov.tab[1,6]
	
	print(paste0(study_name,",Species:Kendall"))
	temp_adonis <- adonis(as.dist(1-cor.fk(t(cmd3_select_age_final_species_clr[study_samples,]))/2)~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,10] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,10] <- temp_adonis$aov.tab[1,6]
	
	##
	
	cmd3_est_age_beta[i,11] <- 0
	cmd3_p_val_age_beta[i,11] <- 1
	
	print(paste0(study_name,",Pathway:Bray"))
	temp_adonis <- adonis(vegdist(cmd3_select_age_final_pathway[study_samples,],method="bray")~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,12] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,12] <- temp_adonis$aov.tab[1,6]
	
	print(paste0(study_name,",Pathway:Jaccard"))
	temp_adonis <- adonis(vegdist(cmd3_select_age_final_pathway[study_samples,],method="jaccard")~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,13] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,13] <- temp_adonis$aov.tab[1,6]
	
	print(paste0(study_name,",Pathway:Manhattan"))
	temp_adonis <- adonis(vegdist(cmd3_select_age_final_pathway_clr[study_samples,],method="euclidean")~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,14] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,14] <- temp_adonis$aov.tab[1,6]
	
	print(paste0(study_name,",Pathway:Kendall"))
	t <- as.matrix(1-cor.fk(t(cmd3_select_age_final_pathway_clr[study_samples,]))/2)
	t[is.na(t)] <- 1
	temp_adonis <- adonis(as.dist(t)~df_cmd3_controls_diversity_uniqueness[study_samples,"age"])
	cmd3_est_age_beta[i,15] <- temp_adonis$aov.tab[1,5]
	cmd3_p_val_age_beta[i,15] <- temp_adonis$aov.tab[1,6]
	
}

### END OF Step B ###

### START OF Step C ###
print("Step C")
cmd3_rlm_est_age_sum_stat <- as.data.frame(matrix(NA,length(sorted_cmd3_age_landscape_studies_controls),15))
rownames(cmd3_rlm_est_age_sum_stat) <- sorted_cmd3_age_landscape_studies_controls
colnames(cmd3_rlm_est_age_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

cmd3_rlm_p_val_age_sum_stat <- as.data.frame(matrix(NA,length(sorted_cmd3_age_landscape_studies_controls),15))
rownames(cmd3_rlm_p_val_age_sum_stat) <- sorted_cmd3_age_landscape_studies_controls
colnames(cmd3_rlm_p_val_age_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

for(i in 1:length(sorted_cmd3_age_landscape_studies_controls))
{
	study_name = sorted_cmd3_age_landscape_studies_controls[i]
	
	temp_rlm <- rlm(genus_shannon~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,],psi = psi.bisquare)
	cmd3_rlm_est_age_sum_stat[i,1] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,1] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(species_shannon~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,],psi = psi.bisquare)
	cmd3_rlm_est_age_sum_stat[i,2] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,2] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(genus_bray_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,3] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,3] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(genus_jaccard_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,4] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,4] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(genus_aitchison_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,5] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,5] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(genus_kendall_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,6] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,6] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(species_bray_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,7] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,7] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(species_jaccard_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,8] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,8] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(species_aitchison_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,9] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,9] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(species_kendall_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,10] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,10] <- f.robftest(temp_rlm)$p.value
	
	#
	temp_rlm <- rlm(pathway_shannon~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,],psi = psi.bisquare)
	cmd3_rlm_est_age_sum_stat[i,11] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,11] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(pathway_bray_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,12] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,12] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(pathway_jaccard_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,13] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,13] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(pathway_aitchison_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,14] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,14] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(pathway_kendall_uniqueness~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_sum_stat[i,15] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_sum_stat[i,15] <- f.robftest(temp_rlm)$p.value
	
	
}

cmd3_rlm_q_val_age_sum_stat <- t(apply(cmd3_rlm_p_val_age_sum_stat,1,function(x)(p.adjust(x,method="fdr"))))

cmd3_rlm_dir_age_sum_stat <- as.data.frame(matrix(NA,length(sorted_cmd3_age_landscape_studies_controls),15))
rownames(cmd3_rlm_dir_age_sum_stat) <- sorted_cmd3_age_landscape_studies_controls
colnames(cmd3_rlm_dir_age_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

for(i in 1:nrow(cmd3_rlm_dir_age_sum_stat))
{
	for(j in 1:ncol(cmd3_rlm_dir_age_sum_stat))
	{
		cmd3_rlm_dir_age_sum_stat[i,j] <- ifelse(cmd3_rlm_q_val_age_sum_stat[i,j]<=0.10,3*sign(cmd3_rlm_est_age_sum_stat[i,j]),ifelse(cmd3_rlm_p_val_age_sum_stat[i,j]<=0.05,2*sign(cmd3_rlm_est_age_sum_stat[i,j]),1*sign(cmd3_rlm_est_age_sum_stat[i,j])))
	}
}
### END OF Step B ###

### START OF Step C ###
cmd3_rlm_est_age_ac_sum_stat <- as.data.frame(matrix(NA,length(sorted_cmd3_age_landscape_studies_controls),15))
rownames(cmd3_rlm_est_age_ac_sum_stat) <- sorted_cmd3_age_landscape_studies_controls
colnames(cmd3_rlm_est_age_ac_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

cmd3_rlm_p_val_age_ac_sum_stat <- as.data.frame(matrix(NA,length(sorted_cmd3_age_landscape_studies_controls),15))
rownames(cmd3_rlm_p_val_age_ac_sum_stat) <- sorted_cmd3_age_landscape_studies_controls
colnames(cmd3_rlm_p_val_age_ac_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

for(i in 1:length(sorted_cmd3_age_landscape_studies_controls))
{
	study_name = sorted_cmd3_age_landscape_studies_controls[i]
	
	temp_rlm <- rlm(genus_shannon~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,],psi = psi.bisquare)
	cmd3_rlm_est_age_ac_sum_stat[i,1] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_ac_sum_stat[i,1] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(species_shannon~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,],psi = psi.bisquare)
	cmd3_rlm_est_age_ac_sum_stat[i,2] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_ac_sum_stat[i,2] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(genus_bray_uniqueness~genus_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,3] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,3] <- (f.robftest(temp_rlm,var="age"))$p.value
	
	temp_rlm <- rlm(genus_jaccard_uniqueness~genus_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,4] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,4] <- (f.robftest(temp_rlm,var="age"))$p.value
	
	temp_rlm <- rlm(genus_aitchison_uniqueness~genus_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,5] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,5] <- (f.robftest(temp_rlm,var="age"))$p.value
	
	temp_rlm <- rlm(genus_kendall_uniqueness~genus_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,6] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,6] <- (f.robftest(temp_rlm,var="age"))$p.value
	
	temp_rlm <- rlm(species_bray_uniqueness~species_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,7] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,7] <- (f.robftest(temp_rlm,var="age"))$p.value
	
	temp_rlm <- rlm(species_jaccard_uniqueness~species_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,8] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,8] <- (f.robftest(temp_rlm,var="age"))$p.value
	
	temp_rlm <- rlm(species_aitchison_uniqueness~species_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,9] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,9] <- (f.robftest(temp_rlm,var="age"))$p.value
	
	temp_rlm <- rlm(species_kendall_uniqueness~species_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,10] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,10] <- (f.robftest(temp_rlm,var="age"))$p.value
	#
	temp_rlm <- rlm(pathway_shannon~age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,],psi = psi.bisquare)
	cmd3_rlm_est_age_ac_sum_stat[i,11] <- as.numeric(temp_rlm$coefficients[2])
	cmd3_rlm_p_val_age_ac_sum_stat[i,11] <- f.robftest(temp_rlm)$p.value
	
	temp_rlm <- rlm(pathway_bray_uniqueness~species_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,12] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,12] <- (f.robftest(temp_rlm,var="age"))$p.value
	
	temp_rlm <- rlm(pathway_jaccard_uniqueness~species_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,13] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,13] <- (f.robftest(temp_rlm,var="age"))$p.value
	
	temp_rlm <- rlm(pathway_aitchison_uniqueness~species_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,14] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,14] <- (f.robftest(temp_rlm,var="age"))$p.value
	
	temp_rlm <- rlm(pathway_kendall_uniqueness~species_shannon+age,df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name==study_name,])
	cmd3_rlm_est_age_ac_sum_stat[i,15] <- as.numeric(temp_rlm$coefficients[3])
	cmd3_rlm_p_val_age_ac_sum_stat[i,15] <- (f.robftest(temp_rlm,var="age"))$p.value
}

cmd3_rlm_q_val_age_ac_sum_stat <- t(apply(cmd3_rlm_p_val_age_ac_sum_stat,1,function(x)(p.adjust(x,method="fdr"))))

cmd3_rlm_dir_age_ac_sum_stat <- as.data.frame(matrix(NA,length(sorted_cmd3_age_landscape_studies_controls),15))
rownames(cmd3_rlm_dir_age_ac_sum_stat) <- sorted_cmd3_age_landscape_studies_controls
colnames(cmd3_rlm_dir_age_ac_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

for(i in 1:nrow(cmd3_rlm_dir_age_ac_sum_stat))
{
	for(j in 1:ncol(cmd3_rlm_dir_age_ac_sum_stat))
	{
		cmd3_rlm_dir_age_ac_sum_stat[i,j] <- ifelse(cmd3_rlm_q_val_age_ac_sum_stat[i,j]<=0.10,3*sign(cmd3_rlm_est_age_ac_sum_stat[i,j]),ifelse(cmd3_rlm_p_val_age_ac_sum_stat[i,j]<=0.05,2*sign(cmd3_rlm_est_age_ac_sum_stat[i,j]),1*sign(cmd3_rlm_est_age_ac_sum_stat[i,j])))
	}
}

### END OF Step C ###

### START OF Step D ###
library(robumeta)
library(metafor)
library(dplyr)
library(effsize)

compute_meta_corr <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(NA,length(grouping_list),3))
	colnames(temp_meta) <- c("dataset","ri","ni")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		temp_meta[i,2] <- cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2])
		temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	temp_meta <- escalc(measure="ZCOR",ri=ri,ni=ni,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	return(res)
}

compute_detection <- function(data,var1_list,grouping_variable,grouping_list)
{
	detection_matrix <- data.frame(matrix(0,length(var1_list),length(grouping_list)))
	rownames(detection_matrix) <- var1_list
	colnames(detection_matrix) <- grouping_list
	for(i in 1:length(var1_list))
	{
		var1 <- var1_list[i]
		for(j in 1:length(grouping_list))
		{
			group <- grouping_list[j]
			detection_matrix[i,j] <- length(which(data[data[,grouping_variable]==group,var1]>0))/length(data[data[,grouping_variable]==group,var1])
		}
	}
	return(detection_matrix)
}

compute_hedges_g <- function(data,var1_list,var2,grouping_variable,grouping_list)
{
	hedges_matrix <- data.frame(matrix(0,length(var1_list),length(grouping_list)))
	colnames(hedges_matrix) <- grouping_list
	rownames(hedges_matrix) <- var1_list
	p_value_matrix <- data.frame(matrix(1,length(var1_list),length(grouping_list)))
	colnames(p_value_matrix) <- grouping_list
	rownames(p_value_matrix) <- var1_list
	for(i in 1:length(var1_list))
	{
		var1 <- var1_list[i]
		for(j in 1:length(grouping_list))
		{
			group <- grouping_list[j]
			data_group_Case <- data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1]
			data_group_Control <- data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1]
			hedges_matrix[i,j] <- as.numeric(effsize::cohen.d(data_group_Case,data_group_Control,hedges.correction=TRUE)$estimate)
			p_value_matrix[i,j] <- as.numeric(wilcox.test(data_group_Case,data_group_Control)$p.value)
		}
	}
	return_list = list("hedges"=hedges_matrix,"p_value"=p_value_matrix)
	return(return_list)
}

compute_meta_effsize <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),7))
	colnames(temp_meta) <- c("dataset","m1i","m2i","sd1i","sd2i","n1i","n2i")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		print(group)
		data_group_Case <- data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1]
		data_group_Control <- data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1]
		temp_meta[i,2] <- mean(data_group_Case)
		temp_meta[i,3] <- mean(data_group_Control)
		temp_meta[i,4] <- sd(data_group_Case)
		temp_meta[i,5] <- sd(data_group_Control)
		temp_meta[i,6] <- length(data_group_Case)
		temp_meta[i,7] <- length(data_group_Control)
		#levels <- unique(data[data[,grouping_variable]==group,var2])
		#temp <- effsize::cohen.d(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1],data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1])$estimate
		#temp <- ifelse(is.nan(temp),0,temp)
		#temp <- ifelse(abs(temp)>1,0.99*sign(temp),temp)
		#print(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1])
		#temp_meta[i,2] <- ifelse(is.nan(temp),0,temp)
		#temp_meta[i,2] <- #cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2])
		#temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	#temp_meta <- temp_meta %>% select(study_id, ri:ni)
	temp_meta <- escalc(measure="SMD",m1i=m1i,m2i=m2i,sd1i=sd1i,sd2i=sd2i,n1i=n1i,n2i=n2i,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	#res$slabs <- rownames(temp_meta)
	return(res)
}

compute_meta_lm <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),6))
	colnames(temp_meta) <- c("dataset","ti","ni","mi","pi","di")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		print(group)
		f <- as.formula(paste0(var1,"~",var2))
		temp_rlm <- rlm(f,data=data[data[,grouping_variable]==group,])
		summary_temp_rlm <- summary(temp_rlm)
		temp_meta[i,2] <- summary_temp_rlm$coefficients[2,3]
		temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
		temp_meta[i,4] <- 1
		temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
		temp_meta[i,6] <- sign(temp_meta[i,2])
		#levels <- unique(data[data[,grouping_variable]==group,var2])
		#temp <- effsize::cohen.d(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1],data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1])$estimate
		#temp <- ifelse(is.nan(temp),0,temp)
		#temp <- ifelse(abs(temp)>1,0.99*sign(temp),temp)
		#print(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1])
		#temp_meta[i,2] <- ifelse(is.nan(temp),0,temp)
		#temp_meta[i,2] <- #cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2])
		#temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	#temp_meta <- temp_meta %>% select(study_id, ri:ni)
	temp_meta <- escalc(measure="ZPCOR",mi=mi,ni=ni,ti=ti,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	res$slabs <- rownames(temp_meta)
	return_list <- list("df_studies"=temp_meta,"model"=res)
	return(return_list)
}

cmd3_rem_genus_shannon_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_shannon","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_genus_shannon_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_shannon","age","study_name",EU_NA_Studies)
cmd3_rem_genus_shannon_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_shannon","age","study_name",EA_Studies)
cmd3_rem_genus_shannon_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_shannon","age","study_name",Other_Studies)

cmd3_rem_species_shannon_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_shannon","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_species_shannon_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_shannon","age","study_name",EU_NA_Studies)
cmd3_rem_species_shannon_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_shannon","age","study_name",EA_Studies)
cmd3_rem_species_shannon_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_shannon","age","study_name",Other_Studies)

cmd3_rem_genus_bray_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_bray_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_genus_bray_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_bray_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_genus_bray_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_bray_uniqueness","age","study_name",EA_Studies)
cmd3_rem_genus_bray_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_bray_uniqueness","age","study_name",Other_Studies)

cmd3_rem_genus_jaccard_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_genus_jaccard_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_genus_jaccard_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","study_name",EA_Studies)
cmd3_rem_genus_jaccard_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","study_name",Other_Studies)

cmd3_rem_genus_aitchison_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_genus_aitchison_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_genus_aitchison_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","study_name",EA_Studies)
cmd3_rem_genus_aitchison_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","study_name",Other_Studies)

cmd3_rem_genus_kendall_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_genus_kendall_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_genus_kendall_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","study_name",EA_Studies)
cmd3_rem_genus_kendall_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","study_name",Other_Studies)

cmd3_rem_species_bray_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_bray_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_species_bray_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_bray_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_species_bray_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_bray_uniqueness","age","study_name",EA_Studies)
cmd3_rem_species_bray_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_bray_uniqueness","age","study_name",Other_Studies)

cmd3_rem_species_jaccard_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_species_jaccard_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_species_jaccard_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","study_name",EA_Studies)
cmd3_rem_species_jaccard_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","study_name",Other_Studies)

cmd3_rem_species_aitchison_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_species_aitchison_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_species_aitchison_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","study_name",EA_Studies)
cmd3_rem_species_aitchison_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","study_name",Other_Studies)

cmd3_rem_species_kendall_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_kendall_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_species_kendall_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_kendall_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_species_kendall_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_kendall_uniqueness","age","study_name",EA_Studies)
cmd3_rem_species_kendall_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"species_kendall_uniqueness","age","study_name",Other_Studies)

cmd3_rem_pathway_shannon_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_shannon","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_pathway_shannon_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_shannon","age","study_name",EU_NA_Studies)
cmd3_rem_pathway_shannon_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_shannon","age","study_name",EA_Studies)
cmd3_rem_pathway_shannon_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_shannon","age","study_name",Other_Studies)

cmd3_rem_pathway_bray_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_bray_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_pathway_bray_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_bray_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_pathway_bray_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_bray_uniqueness","age","study_name",EA_Studies)
cmd3_rem_pathway_bray_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_bray_uniqueness","age","study_name",Other_Studies)

cmd3_rem_pathway_jaccard_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_jaccard_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_pathway_jaccard_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_jaccard_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_pathway_jaccard_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_jaccard_uniqueness","age","study_name",EA_Studies)
cmd3_rem_pathway_jaccard_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_jaccard_uniqueness","age","study_name",Other_Studies)

cmd3_rem_pathway_aitchison_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_aitchison_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_pathway_aitchison_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_aitchison_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_pathway_aitchison_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_aitchison_uniqueness","age","study_name",EA_Studies)
cmd3_rem_pathway_aitchison_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_aitchison_uniqueness","age","study_name",Other_Studies)

cmd3_rem_pathway_kendall_uniqueness_age <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_kendall_uniqueness","age","study_name",sorted_cmd3_age_landscape_studies_controls)
cmd3_rem_pathway_kendall_uniqueness_age_EU_NA <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_kendall_uniqueness","age","study_name",EU_NA_Studies)
cmd3_rem_pathway_kendall_uniqueness_age_EA_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_kendall_uniqueness","age","study_name",EA_Studies)
cmd3_rem_pathway_kendall_uniqueness_age_Other_Studies <- compute_meta_lm(df_cmd3_controls_diversity_uniqueness,"pathway_kendall_uniqueness","age","study_name",Other_Studies)


### End of Step D ###

save(list=(c(ls(pattern="df_cmd3_"),ls(pattern="est"),ls(pattern="p_val"),ls(pattern="q_val"),ls(pattern="dir"))),file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_stage1_results.RData")
save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")
rm(list=ls())
