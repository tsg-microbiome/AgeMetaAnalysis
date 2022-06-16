# This sub-pipeline works on the Odamaki et al dataset and computes the following:
# Associations between Shannon Diversity and the four different measures of Uniqueness at the Species and Genus level.
# Associations between gut microbiome beta-diversity (computed using the four different distance measures; the same ones corresponding to the four measures of uniqueness) and age.
# Associations between the different measures of uniqueness (four measures of uniqueness each corresponding to the species and genus levels) and Shannon Diversity with age.

library(MASS)
library(sfsmisc)
library(gplots)
library(RColorBrewer)
library(vegan)
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\Odamaki\\odamaki_analysis_2021_Revision.RData")


#disease_list <- c("Atherosclerosis","Cholecystitis","Colitis","Constipation","Diarrhea","Fatty_liver","Gastritis","IBS","Kidneystone","Rheumatoid_arthritis","Metabolic_syndrome","T2DM")

odamaki_all_controls <- rownames(df_odamaki_diversity_uniqueness)

df_odamaki_controls_diversity_uniqueness <- df_odamaki_diversity_uniqueness[odamaki_all_controls,]

### START OF Step A ###
print("StepA: Diversity v/s Uniqueness")
odamaki_rlm_est_diversity_uniqueness <- as.data.frame(matrix(NA,1,8))
rownames(odamaki_rlm_est_diversity_uniqueness) <- "Odamaki"
colnames(odamaki_rlm_est_diversity_uniqueness) <- c("genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

odamaki_rlm_p_val_diversity_uniqueness <- as.data.frame(matrix(NA,1,8))
rownames(odamaki_rlm_p_val_diversity_uniqueness) <- "Odamaki"
colnames(odamaki_rlm_p_val_diversity_uniqueness) <- c("genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

temp_rlm <- rlm(genus_bray_uniqueness~genus_shannon,df_odamaki_diversity_uniqueness)
odamaki_rlm_est_diversity_uniqueness[1,1] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_diversity_uniqueness[1,1] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(genus_jaccard_uniqueness~genus_shannon,df_odamaki_diversity_uniqueness)
odamaki_rlm_est_diversity_uniqueness[1,2] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_diversity_uniqueness[1,2] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(genus_aitchison_uniqueness~genus_shannon,df_odamaki_diversity_uniqueness)
odamaki_rlm_est_diversity_uniqueness[1,3] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_diversity_uniqueness[1,3] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(genus_kendall_uniqueness~genus_shannon,df_odamaki_diversity_uniqueness)
odamaki_rlm_est_diversity_uniqueness[1,4] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_diversity_uniqueness[1,4] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(species_bray_uniqueness~species_shannon,df_odamaki_diversity_uniqueness)
odamaki_rlm_est_diversity_uniqueness[1,5] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_diversity_uniqueness[1,5] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(species_jaccard_uniqueness~species_shannon,df_odamaki_diversity_uniqueness)
odamaki_rlm_est_diversity_uniqueness[1,6] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_diversity_uniqueness[1,6] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(species_aitchison_uniqueness~species_shannon,df_odamaki_diversity_uniqueness)
odamaki_rlm_est_diversity_uniqueness[1,7] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_diversity_uniqueness[1,7] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(species_kendall_uniqueness~species_shannon,df_odamaki_diversity_uniqueness)
odamaki_rlm_est_diversity_uniqueness[1,8] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_diversity_uniqueness[1,8] <- f.robftest(temp_rlm)$p.value

odamaki_rlm_q_val_diversity_uniqueness <- t(apply(odamaki_rlm_p_val_diversity_uniqueness,1,p.adjust))

odamaki_rlm_dir_diversity_uniqueness <- as.data.frame(matrix(NA,1,8))
rownames(odamaki_rlm_dir_diversity_uniqueness) <- "Odamaki"
colnames(odamaki_rlm_dir_diversity_uniqueness) <- c("genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

for(j in 1:ncol(odamaki_rlm_dir_diversity_uniqueness))
{
	odamaki_rlm_dir_diversity_uniqueness[1,j] <- ifelse(odamaki_rlm_q_val_diversity_uniqueness[1,j]<=0.05,3*sign(odamaki_rlm_est_diversity_uniqueness[1,j]),ifelse(odamaki_rlm_q_val_diversity_uniqueness[1,j]<=0.1,2*sign(odamaki_rlm_est_diversity_uniqueness[1,j]),sign(odamaki_rlm_est_diversity_uniqueness[1,j])))
}

### END OF Step A ###

### START OF Step B ###
print("StepB: Beta Diversity v/s Age")
odamaki_est_age_beta <- as.data.frame(matrix(0,1,10))
rownames(odamaki_est_age_beta) <- "Odamaki"
colnames(odamaki_est_age_beta) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")
odamaki_p_val_age_beta <- as.data.frame(matrix(0,1,10))
rownames(odamaki_p_val_age_beta) <- "Odamaki"
colnames(odamaki_p_val_age_beta) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

odamaki_est_age_beta[1,1] <- 0
odamaki_p_val_age_beta[1,1] <- 1
	
odamaki_est_age_beta[1,2] <- 0
odamaki_p_val_age_beta[1,2] <- 1

study_samples <- odamaki_all_controls

print("Genus:Bray")
temp_adonis <- adonis(as.dist(dist_genus_bray[study_samples,study_samples])~df_odamaki_controls_diversity_uniqueness[,"age"])
odamaki_est_age_beta[1,3] <- temp_adonis$aov.tab[1,5]
odamaki_p_val_age_beta[1,3] <- temp_adonis$aov.tab[1,6]
print("Genus:Jaccard")
temp_adonis <- adonis(as.dist(dist_genus_jaccard[study_samples,study_samples])~df_odamaki_controls_diversity_uniqueness[,"age"])
odamaki_est_age_beta[1,4] <- temp_adonis$aov.tab[1,5]
odamaki_p_val_age_beta[1,4] <- temp_adonis$aov.tab[1,6]
print("Genus:Manhattan")
temp_adonis <- adonis(as.dist(dist_genus_aitchison[study_samples,study_samples])~df_odamaki_controls_diversity_uniqueness[,"age"])
odamaki_est_age_beta[1,5] <- temp_adonis$aov.tab[1,5]
odamaki_p_val_age_beta[1,5] <- temp_adonis$aov.tab[1,6]
print("Genus:Kendall")
temp_adonis <- adonis(as.dist(dist_genus_kendall[study_samples,study_samples])~df_odamaki_controls_diversity_uniqueness[,"age"])
odamaki_est_age_beta[1,6] <- temp_adonis$aov.tab[1,5]
odamaki_p_val_age_beta[1,6] <- temp_adonis$aov.tab[1,6]
print("Species:Bray")
temp_adonis <- adonis(as.dist(dist_species_bray[study_samples,study_samples])~df_odamaki_controls_diversity_uniqueness[,"age"])
odamaki_est_age_beta[1,7] <- temp_adonis$aov.tab[1,5]
odamaki_p_val_age_beta[1,7] <- temp_adonis$aov.tab[1,6]
print("Species:Jaccard")
temp_adonis <- adonis(as.dist(dist_species_jaccard[study_samples,study_samples])~df_odamaki_controls_diversity_uniqueness[,"age"])
odamaki_est_age_beta[1,8] <- temp_adonis$aov.tab[1,5]
odamaki_p_val_age_beta[1,8] <- temp_adonis$aov.tab[1,6]
print("Species:Manhattan")
temp_adonis <- adonis(as.dist(dist_species_aitchison[study_samples,study_samples])~df_odamaki_controls_diversity_uniqueness[,"age"])
odamaki_est_age_beta[1,9] <- temp_adonis$aov.tab[1,5]
odamaki_p_val_age_beta[1,9] <- temp_adonis$aov.tab[1,6]
print("Species:Kendall")
temp_adonis <- adonis(as.dist(dist_species_kendall[study_samples,study_samples])~df_odamaki_controls_diversity_uniqueness[,"age"])
odamaki_est_age_beta[1,10] <- temp_adonis$aov.tab[1,5]
odamaki_p_val_age_beta[1,10] <- temp_adonis$aov.tab[1,6]
### END OF Step B ###

### START OF Step C ###
print("StepC: Uniqueness/Diversity v/s Age")
odamaki_rlm_est_age_sum_stat <- as.data.frame(matrix(NA,1,10))
rownames(odamaki_rlm_est_age_sum_stat) <- "Odamaki"
colnames(odamaki_rlm_est_age_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

odamaki_rlm_p_val_age_sum_stat <- as.data.frame(matrix(NA,1,10))
rownames(odamaki_rlm_p_val_age_sum_stat) <- "Odamaki"
colnames(odamaki_rlm_p_val_age_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

temp_rlm <- rlm(genus_shannon~age,data=df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_sum_stat[1,1] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_sum_stat[1,1] <- f.robftest(temp_rlm)$p.value
	
temp_rlm <- rlm(species_shannon~age,data=df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_sum_stat[1,2] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_sum_stat[1,2] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(genus_bray_uniqueness~age,data=df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_sum_stat[1,3] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_sum_stat[1,3] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(genus_jaccard_uniqueness~age,data=df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_sum_stat[1,4] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_sum_stat[1,4] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(genus_aitchison_uniqueness~age,data=df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_sum_stat[1,5] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_sum_stat[1,5] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(genus_kendall_uniqueness~age,data=df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_sum_stat[1,6] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_sum_stat[1,6] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(species_bray_uniqueness~age,data=df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_sum_stat[1,7] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_sum_stat[1,7] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(species_jaccard_uniqueness~age,data=df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_sum_stat[1,8] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_sum_stat[1,8] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(species_aitchison_uniqueness~age,data=df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_sum_stat[1,9] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_sum_stat[1,9] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(species_kendall_uniqueness~age,data=df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_sum_stat[1,10] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_sum_stat[1,10] <- f.robftest(temp_rlm)$p.value

odamaki_rlm_q_val_age_sum_stat <- t(apply(odamaki_rlm_p_val_age_sum_stat,1,function(x)(p.adjust(x,method="fdr"))))

odamaki_rlm_dir_age_sum_stat <- as.data.frame(matrix(NA,1,10))
rownames(odamaki_rlm_dir_age_sum_stat) <- "Odamaki"
colnames(odamaki_rlm_dir_age_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

for(j in 1:ncol(odamaki_rlm_dir_age_sum_stat))
{
	odamaki_rlm_dir_age_sum_stat[1,j] <- ifelse(odamaki_rlm_q_val_age_sum_stat[1,j]<=0.10,3*sign(odamaki_rlm_est_age_sum_stat[1,j]),ifelse(odamaki_rlm_p_val_age_sum_stat[1,j]<=0.05,2*sign(odamaki_rlm_est_age_sum_stat[1,j]),1*sign(odamaki_rlm_est_age_sum_stat[1,j])))
}

### END OF Step C ###

### START OF Step D ###
print("StepD: Alpha-Corrected Uniqueness v/s Age")
odamaki_rlm_est_age_ac_sum_stat <- as.data.frame(matrix(NA,1,10))
rownames(odamaki_rlm_est_age_ac_sum_stat) <- "Odamaki"
colnames(odamaki_rlm_est_age_ac_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

odamaki_rlm_p_val_age_ac_sum_stat <- as.data.frame(matrix(NA,1,10))
rownames(odamaki_rlm_p_val_age_ac_sum_stat) <- "Odamaki"
colnames(odamaki_rlm_p_val_age_ac_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

temp_rlm <- rlm(genus_shannon~age,data=df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_ac_sum_stat[1,1] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_ac_sum_stat[1,1] <- f.robftest(temp_rlm)$p.value
	
temp_rlm <- rlm(species_shannon~age,df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_ac_sum_stat[1,2] <- as.numeric(temp_rlm$coefficients[2])
odamaki_rlm_p_val_age_ac_sum_stat[1,2] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(genus_bray_uniqueness~genus_shannon+age,df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_ac_sum_stat[1,3] <- as.numeric(temp_rlm$coefficients[3])
odamaki_rlm_p_val_age_ac_sum_stat[1,3] <- (f.robftest(temp_rlm,var="age"))$p.value
	
temp_rlm <- rlm(genus_jaccard_uniqueness~genus_shannon+age,df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_ac_sum_stat[1,4] <- as.numeric(temp_rlm$coefficients[3])
odamaki_rlm_p_val_age_ac_sum_stat[1,4] <- (f.robftest(temp_rlm,var="age"))$p.value
	
temp_rlm <- rlm(genus_aitchison_uniqueness~genus_shannon+age,df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_ac_sum_stat[1,5] <- as.numeric(temp_rlm$coefficients[3])
odamaki_rlm_p_val_age_ac_sum_stat[1,5] <- (f.robftest(temp_rlm,var="age"))$p.value
	
temp_rlm <- rlm(genus_kendall_uniqueness~genus_shannon+age,df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_ac_sum_stat[1,6] <- as.numeric(temp_rlm$coefficients[3])
odamaki_rlm_p_val_age_ac_sum_stat[1,6] <- (f.robftest(temp_rlm,var="age"))$p.value
	
temp_rlm <- rlm(species_bray_uniqueness~species_shannon+age,df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_ac_sum_stat[1,7] <- as.numeric(temp_rlm$coefficients[3])
odamaki_rlm_p_val_age_ac_sum_stat[1,7] <- (f.robftest(temp_rlm,var="age"))$p.value
	
temp_rlm <- rlm(species_jaccard_uniqueness~species_shannon+age,df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_ac_sum_stat[1,8] <- as.numeric(temp_rlm$coefficients[3])
odamaki_rlm_p_val_age_ac_sum_stat[1,8] <- (f.robftest(temp_rlm,var="age"))$p.value
	
temp_rlm <- rlm(species_aitchison_uniqueness~species_shannon+age,df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_ac_sum_stat[1,9] <- as.numeric(temp_rlm$coefficients[3])
odamaki_rlm_p_val_age_ac_sum_stat[1,9] <- (f.robftest(temp_rlm,var="age"))$p.value
	
temp_rlm <- rlm(species_kendall_uniqueness~species_shannon+age,df_odamaki_controls_diversity_uniqueness,psi = psi.bisquare)
odamaki_rlm_est_age_ac_sum_stat[1,10] <- as.numeric(temp_rlm$coefficients[3])
odamaki_rlm_p_val_age_ac_sum_stat[1,10] <- (f.robftest(temp_rlm,var="age"))$p.value
	
	
odamaki_rlm_q_val_age_ac_sum_stat <- t(apply(odamaki_rlm_p_val_age_ac_sum_stat,1,function(x)(p.adjust(x,method="fdr"))))

odamaki_rlm_dir_age_ac_sum_stat <- as.data.frame(matrix(NA,1,10))
rownames(odamaki_rlm_dir_age_ac_sum_stat) <- "Odamaki"
colnames(odamaki_rlm_dir_age_ac_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

for(j in 1:ncol(odamaki_rlm_dir_age_ac_sum_stat))
{
	odamaki_rlm_dir_age_ac_sum_stat[1,j] <- ifelse(odamaki_rlm_q_val_age_ac_sum_stat[1,j]<=0.10,3*sign(odamaki_rlm_est_age_ac_sum_stat[1,j]),ifelse(odamaki_rlm_p_val_age_ac_sum_stat[1,j]<=0.05,2*sign(odamaki_rlm_est_age_ac_sum_stat[1,j]),1*sign(odamaki_rlm_est_age_ac_sum_stat[1,j])))
}

### END OF Step D###
save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\Odamaki\\odamaki_analysis_2021_Revision.RData")

save(list=(c(ls(pattern="df_odamaki_"),ls(pattern="est"),ls(pattern="p_val"),ls(pattern="q_val"),ls(pattern="dir"))),file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\odamaki_stage1_results.RData")

rm(list=ls())
