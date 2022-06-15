# The sub-pipeline works on the disease information corresponding to the subjects in the He et al data repository and performs the following tasks:
# Computes the number of patient gut microbiomes, the maximum and minimum age of the subjects corresponding to each disease
# Computes the association of the relative abundances of the different Highly Detected Species and the different diseases
# Subsequently specifically focusses on the disease associations of the Multiple Disease Enriched and Multiple Disease Depleted (previously identified
# in Ghosh et al eLife, 2020)

# The sub-pipeline utilizes the workspaces he_age_analysis.RData already provided in this github and he_stage1_results.RData which is generated during the execution
# of run_he_analysis_stage1b.R code (which in turn is executed after the run_he_stage1a.R code) and the he_stage2a_results.RData generated during the execution of 
# run_he_analysis_stage2a_1.R
# The workspace he_Analysis_2021_Revision.RData is a running workspace that is iteratively loaded and saved during the running of each of the previous 
# He et al analysis code (run_he_et_al_analysis_stage1a.R followed by run_he_et_al_analysis_stage1b.R followed by run_he_et_al_stage2a_1.R followed by 
# The detailed order of the running of each individual sub-pipeline is provided in the Readme file in the github corresponding to this study.

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\He\\he_age_analysis.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\he_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\he_stage2a_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\He\\he_Analysis_2021_Revision.RData")

rank_scale=function(x)
{
	x <- rank(x);
	y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
	y <- ifelse(is.nan(y),0,y)
	return(y);
}

rank_scale1=function(x)
{
	x <- ifelse(is.na(x),NA,rank(x))
	min_x <- min(x[!is.na(x)])
	max_x <- max(x[!is.na(x)])
	y <- ifelse(is.na(x),NA,(x-min_x)/(max_x-min_x))
	return(y);
}

rank_scale2=function(x,range_min,range_max)
{
	x <- rank(x);
	y <- range_min + (range_max - range_min)*(rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
	return(y);
}

range_scale=function(x)
{
	y <- (x-min(x))/(max(x)-min(x));
	return(y);
}

range_scale2=function(x,range_min,range_max)
{
	y <- range_min + (range_max - range_min)*(x-min(x))/(max(x)-min(x));
	return(y);
}


library(robumeta)
library(metafor)
library(dplyr)
library(effsize)
library(MASS)
library(sfsmisc)
library(gplots)
library(RColorBrewer)
library(metap)

disease_list <- c("Atherosclerosis","Cholecystitis","Colitis","Constipation","Diarrhea","Fatty_liver","Gastritis","IBS","Kidneystone","Rheumatoid_arthritis","Metabolic_syndrome","T2DM")

he_all_controls <- names(which(apply(he_select_age_final_metadata[,disease_list],1,function(x)(length(x[(x=="y")&(!is.na(x))])))==0))

HighlyDetectedSpecies <- names(which(100*colSums(apply(he_select_age_final_species,2,function(x)(ifelse(x>0,1,0))))/nrow(he_select_age_final_species)>=5))

he_combined_df_sum_stat_species <- as.data.frame(cbind(df_he_diversity_uniqueness,he_select_age_final_species[rownames(df_he_diversity_uniqueness),HighlyDetectedSpecies]))

disease_patient_numbers <- matrix(NA,length(disease_list),3)
rownames(disease_patient_numbers) <- c(disease_list)
colnames(disease_patient_numbers) <- c("number","minimum_age","maximum_age")

for(i in 1:length(disease_list))
{
	disease_name <- disease_list[i]
	assign(paste0("patient_list_",disease_name),rownames(he_select_age_final_metadata)[which(!is.na(he_select_age_final_metadata[,disease_name])&(he_select_age_final_metadata[,disease_name]=="y"))])
	disease_patient_numbers[disease_name,1] <- length(rownames(he_select_age_final_metadata)[which(!is.na(he_select_age_final_metadata[,disease_name])&(he_select_age_final_metadata[,disease_name]=="y"))])
	disease_patient_numbers[disease_name,2] <- min(he_select_age_final_metadata[which(!is.na(he_select_age_final_metadata[,disease_name])&(he_select_age_final_metadata[,disease_name]=="y")),"Age"])
	disease_patient_numbers[disease_name,3] <- max(he_select_age_final_metadata[which(!is.na(he_select_age_final_metadata[,disease_name])&(he_select_age_final_metadata[,disease_name]=="y")),"Age"])
}

G_Markers <- intersect(c("Bacteroides_fragilis","Actinomyces_odontolyticus","Bifidobacterium_dentium","Clostridium_clostridioforme","Clostridium_nexile","Clostridium_ramosum","Dialister_invisus","Eggerthella_lenta","Escherichia_coli","Fusobacterium_nucleatum","Granulicatella_adiacens","Lactococcus_lactis","Rothia_mucilaginosa","Ruminococcus_gnavus","Streptococcus_anginosus","Streptococcus_infantis","Streptococcus_salivarius","Streptococcus_sanguinis","Streptococcus_vestibularis","Subdoligranulum_variabile","Veillonella_atypica","Clostridium_asparagiforme","Clostridium_bolteae","Clostridium_citroniae","Clostridium_hathewayi","Clostridium_symbiosum","Streptococcus_parasanguinis","Solobacterium_moorei","Ruminococcus_torques","Streptococcus_mitis","Klebsiella_pneumoniae","Streptococcus_australis","Streptococcus_gordonii","Enterobacter_cloacae"),colnames(he_combined_df_sum_stat_species))

L_Markers <- intersect(c("Eubacterium_hallii","Dorea_longicatena","Coprococcus_comes","Coprococcus_catus","Butyrivibrio_crossotus","Bacteroides_uniformis","Alistipes_shahii","Alistipes_indistinctus","Roseburia_hominis","Pseudoflavonifractor_capillosus","Eubacterium_siraeum","Eubacterium_rectale","Coprobacter_fastidiosus","Bifidobacterium_longum","Bifidobacterium_animalis","Barnesiella_intestinihominis","Bacteroides_xylanisolvens","Alistipes_senegalensis","Alistipes_putredinis","Alistipes_onderdonkii","Akkermansia_muciniphila","Faecalibacterium_prausnitzii","Roseburia_inulinivorans","Roseburia_intestinalis"),colnames(he_combined_df_sum_stat_species))

he_combined_df_sum_stat_species$G_Markers_combined <- rowMeans(apply(he_combined_df_sum_stat_species[,G_Markers],2,rank_scale))
he_combined_df_sum_stat_species$L_Markers_combined <- rowMeans(apply(he_combined_df_sum_stat_species[,L_Markers],2,rank_scale))

HighlyDetectedSpecies <- c(HighlyDetectedSpecies,"G_Markers_combined","L_Markers_combined")

G_Markers <- c(G_Markers,"G_Markers_combined")
L_Markers <- c(L_Markers,"L_Markers_combined")

he_rlm_est_species_disease <- as.data.frame(matrix(0,length(HighlyDetectedSpecies),length(disease_list)))
rownames(he_rlm_est_species_disease) <- HighlyDetectedSpecies
colnames(he_rlm_est_species_disease) <- disease_list

he_rlm_pval_species_disease <- as.data.frame(matrix(1,length(HighlyDetectedSpecies),length(disease_list)))
rownames(he_rlm_pval_species_disease) <- HighlyDetectedSpecies
colnames(he_rlm_pval_species_disease) <- disease_list

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	for(j in 1:length(disease_list))
	{
		disease_name <- disease_list[j]
		control_list <- he_all_controls
		diseased_list <- get(paste0("patient_list_",disease_name))
		vec_control <- he_combined_df_sum_stat_species[control_list,species_name]
		vec_disease <- he_combined_df_sum_stat_species[diseased_list,species_name]
		total_vec <- c(vec_disease,vec_control)
		if(length(total_vec[total_vec > 0])>0)
		{
			temp_rlm <- rlm(he_combined_df_sum_stat_species[c(diseased_list,control_list),species_name]~as.factor(ifelse(c(diseased_list,control_list) %in% he_all_controls,0,1)))
			summary_temp_rlm <- summary(temp_rlm)
			he_rlm_est_species_disease[species_name,disease_name] <- summary_temp_rlm$coefficients[2,3]
			he_rlm_pval_species_disease[species_name,disease_name] <- f.robftest(temp_rlm)$p.value
		}
	}
}

he_rlm_qval_species_disease <- apply(he_rlm_pval_species_disease,2,function(x)(p.adjust(x,method="fdr")))

he_rlm_dir_species_disease <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),length(disease_list)))
rownames(he_rlm_dir_species_disease) <- HighlyDetectedSpecies
colnames(he_rlm_dir_species_disease) <- disease_list

he_rlm_est_select_species_disease <- he_rlm_est_species_disease[c(G_Markers,L_Markers),]
he_rlm_pval_select_species_disease <- he_rlm_pval_species_disease[c(G_Markers,L_Markers),]
he_rlm_qval_select_species_disease <- t(apply(he_rlm_pval_select_species_disease,1,function(x)(p.adjust(x,method="fdr"))))

select_markers <- c(G_Markers,L_Markers)

he_rlm_dir_select_species_disease <- as.data.frame(matrix(NA,length(select_markers),length(disease_list)))
rownames(he_rlm_dir_select_species_disease) <- select_markers
colnames(he_rlm_dir_select_species_disease) <- disease_list

for(i in 1:length(select_markers))
{
	for(j in 1:length(disease_list))
	{
		he_rlm_dir_select_species_disease[i,j] <- ifelse(he_rlm_qval_select_species_disease[i,j] <= 0.1,3*sign(he_rlm_est_select_species_disease[i,j]),ifelse(he_rlm_pval_select_species_disease[i,j] <= 0.05,2*sign(he_rlm_est_select_species_disease[i,j]),sign(he_rlm_est_select_species_disease[i,j])))
	}
}

AssociationDirection <- as.data.frame(cbind(apply(he_rlm_dir_select_species_disease,1,function(x)(length(x[x<= -2]))),apply(he_rlm_dir_select_species_disease,1,function(x)(length(x[x>=2])))))

G_Markers_ordered <- intersect(G_Markers,rownames(AssociationDirection[order(AssociationDirection[,2]-AssociationDirection[,1]),]))
L_Markers_ordered <- intersect(L_Markers,rownames(AssociationDirection[order(AssociationDirection[,2]-AssociationDirection[,1]),]))

Marker_List_Sorted <- c(setdiff(c(G_Markers_ordered,L_Markers_ordered),c("G_Markers_combined","L_Markers_combined")),c("G_Markers_combined","L_Markers_combined"))

mat_disease_association <- apply(he_rlm_dir_select_species_disease[Marker_List_Sorted,],2,function(x)(ifelse(x<0,x-0.1,x)))

#heatmap.2(t(mat_disease_association),density="none",trace="none",col=c("skyblue4","skyblue3","powderblue","white","mistyrose1","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(t(mat_disease_association)),rowsep=1:nrow(t(mat_disease_association)),key=FALSE,margins=c(10,5),lhei=c(0.1,5),lwid=c(0.1,5))
for(i in 1:length(HighlyDetectedSpecies))
{
	for(j in 1:length(disease_list))
	{
		he_rlm_dir_species_disease[i,j] <- ifelse(he_rlm_qval_species_disease[i,j] <= 0.1,3*sign(he_rlm_est_species_disease[i,j]),ifelse(he_rlm_pval_species_disease[i,j] <= 0.05,2*sign(he_rlm_est_species_disease[i,j]),sign(he_rlm_est_species_disease[i,j])))
	}
}

he_mat_disease_association <- mat_disease_association
save(df_he_diversity_uniqueness,df_he_controls_diversity_uniqueness,he_select_age_final_species,he_combined_df_sum_stat_species,he_rlm_est_select_species_disease,he_rlm_pval_select_species_disease,he_rlm_qval_select_species_disease,he_mat_disease_association,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\he_stage2b_results.RData") 

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\He\\he_analysis_2021_Revision.RData")

rm(list=ls())
