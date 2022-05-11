load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\AGP\\ag_age_analysis.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ag_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ag_stage2a_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\AGP\\ag_Analysis_2021_Revision.RData")

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

disease_list <- c("alzheimers","asd","cancer","cardiovascular_disease","cdiff","diabetes","ibd","ibs","kidney_disease","liver_disease","lung_disease","migraine","sibo")

ag_all_controls <- names(which(apply(ag_select_age_final_metadata[,disease_list],1,function(x)(length(x[(x=="Yes")&(!is.na(x))])))==0))

HighlyDetectedSpecies <- names(which(100*colSums(apply(ag_select_age_final_species,2,function(x)(ifelse(x>0,1,0))))/nrow(ag_select_age_final_species)>=5))

ag_combined_df_sum_stat_species <- as.data.frame(cbind(df_ag_diversity_uniqueness,ag_select_age_final_species[rownames(df_ag_diversity_uniqueness),HighlyDetectedSpecies]))

disease_patient_numbers <- matrix(NA,length(disease_list),3)
rownames(disease_patient_numbers) <- c(disease_list)
colnames(disease_patient_numbers) <- c("number","minimum_age","maximum_age")

for(i in 1:length(disease_list))
{
	disease_name <- disease_list[i]
	assign(paste0("patient_list_",disease_name),rownames(ag_select_age_final_metadata[!is.na(ag_select_age_final_metadata[,disease_name])&(ag_select_age_final_metadata[,disease_name]=="Yes"),]))
	disease_patient_numbers[disease_name,1] <- length(rownames(ag_select_age_final_metadata[!is.na(ag_select_age_final_metadata[,disease_name])&(ag_select_age_final_metadata[,disease_name]=="Yes"),]))
	disease_patient_numbers[disease_name,2] <- min(ag_select_age_final_metadata[!is.na(ag_select_age_final_metadata[,disease_name])&(ag_select_age_final_metadata[,disease_name]=="Yes"),"age_years"])
	disease_patient_numbers[disease_name,3] <- max(ag_select_age_final_metadata[!is.na(ag_select_age_final_metadata[,disease_name])&(ag_select_age_final_metadata[,disease_name]=="Yes"),"age_years"])
}

all_diseased_list <- NULL
for(i in 1:length(disease_list))
{
	disease_name <- disease_list[i]
	all_diseased_list <- c(all_diseased_list,get(paste0("patient_list_",disease_name)))
}
#patient_list_all_diseases <- unique(all_diseased_list)
#disease_list <- c(disease_list,"all_diseases")

G_Markers <- intersect(c("Bacteroides_fragilis","Actinomyces_odontolyticus","Bifidobacterium_dentium","Clostridium_clostridioforme","Clostridium_nexile","Clostridium_ramosum","Dialister_invisus","Eggerthella_lenta","Escherichia_coli","Fusobacterium_nucleatum","Granulicatella_adiacens","Lactococcus_lactis","Rothia_mucilaginosa","Ruminococcus_gnavus","Streptococcus_anginosus","Streptococcus_infantis","Streptococcus_salivarius","Streptococcus_sanguinis","Streptococcus_vestibularis","Subdoligranulum_variabile","Veillonella_atypica","Clostridium_asparagiforme","Clostridium_bolteae","Clostridium_citroniae","Clostridium_hathewayi","Clostridium_symbiosum","Streptococcus_parasanguinis","Solobacterium_moorei","Ruminococcus_torques","Streptococcus_mitis","Klebsiella_pneumoniae","Streptococcus_australis","Streptococcus_gordonii","Enterobacter_cloacae"),colnames(ag_combined_df_sum_stat_species))

L_Markers <- intersect(c("Eubacterium_hallii","Dorea_longicatena","Coprococcus_comes","Coprococcus_catus","Butyrivibrio_crossotus","Bacteroides_uniformis","Alistipes_shahii","Alistipes_indistinctus","Roseburia_hominis","Pseudoflavonifractor_capillosus","Eubacterium_siraeum","Eubacterium_rectale","Coprobacter_fastidiosus","Bifidobacterium_longum","Bifidobacterium_animalis","Barnesiella_intestinihominis","Bacteroides_xylanisolvens","Alistipes_senegalensis","Alistipes_putredinis","Alistipes_onderdonkii","Akkermansia_muciniphila","Faecalibacterium_prausnitzii","Roseburia_inulinivorans","Roseburia_intestinalis"),colnames(ag_combined_df_sum_stat_species))

ag_combined_df_sum_stat_species$G_Markers_combined <- rowMeans(apply(ag_combined_df_sum_stat_species[,G_Markers],2,rank_scale))
ag_combined_df_sum_stat_species$L_Markers_combined <- rowMeans(apply(ag_combined_df_sum_stat_species[,L_Markers],2,rank_scale))

HighlyDetectedSpecies <- c(HighlyDetectedSpecies,"G_Markers_combined","L_Markers_combined")

G_Markers <- c(G_Markers,"G_Markers_combined")
L_Markers <- c(L_Markers,"L_Markers_combined")

ag_rlm_est_species_disease <- as.data.frame(matrix(0,length(HighlyDetectedSpecies),length(disease_list)))
rownames(ag_rlm_est_species_disease) <- HighlyDetectedSpecies
colnames(ag_rlm_est_species_disease) <- disease_list

ag_rlm_pval_species_disease <- as.data.frame(matrix(1,length(HighlyDetectedSpecies),length(disease_list)))
rownames(ag_rlm_pval_species_disease) <- HighlyDetectedSpecies
colnames(ag_rlm_pval_species_disease) <- disease_list

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	for(j in 1:length(disease_list))
	{
		disease_name <- disease_list[j]
		control_list <- ag_all_controls
		diseased_list <- get(paste0("patient_list_",disease_name))
		vec_control <- ag_combined_df_sum_stat_species[control_list,species_name]
		vec_disease <- ag_combined_df_sum_stat_species[diseased_list,species_name]
		total_vec <- c(vec_disease,vec_control)
		if(length(total_vec[total_vec > 0])>0)
		{
			temp_rlm <- rlm(ag_combined_df_sum_stat_species[c(diseased_list,control_list),species_name]~as.factor(ifelse(c(diseased_list,control_list) %in% ag_all_controls,0,1)))
			summary_temp_rlm <- summary(temp_rlm)
			ag_rlm_est_species_disease[species_name,disease_name] <- summary_temp_rlm$coefficients[2,3]
			ag_rlm_pval_species_disease[species_name,disease_name] <- f.robftest(temp_rlm)$p.value
		}
	}
}

ag_rlm_qval_species_disease <- apply(ag_rlm_pval_species_disease,2,function(x)(p.adjust(x,method="fdr")))

ag_rlm_dir_species_disease <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),length(disease_list)))
rownames(ag_rlm_dir_species_disease) <- HighlyDetectedSpecies
colnames(ag_rlm_dir_species_disease) <- disease_list

for(i in 1:length(HighlyDetectedSpecies))
{
	for(j in 1:length(disease_list))
	{
		ag_rlm_dir_species_disease[i,j] <- ifelse(ag_rlm_qval_species_disease[i,j] <= 0.1,3*sign(ag_rlm_est_species_disease[i,j]),ifelse(ag_rlm_pval_species_disease[i,j] <= 0.05,2*sign(ag_rlm_est_species_disease[i,j]),sign(ag_rlm_est_species_disease[i,j])))
	}
}

ag_rlm_est_select_species_disease <- ag_rlm_est_species_disease[c(G_Markers,L_Markers),]
ag_rlm_pval_select_species_disease <- ag_rlm_pval_species_disease[c(G_Markers,L_Markers),]
ag_rlm_qval_select_species_disease <- t(apply(ag_rlm_pval_select_species_disease,1,function(x)(p.adjust(x,method="fdr"))))

select_markers <- c(G_Markers,L_Markers)

ag_rlm_dir_select_species_disease <- as.data.frame(matrix(NA,length(select_markers),length(disease_list)))
rownames(ag_rlm_dir_select_species_disease) <- select_markers
colnames(ag_rlm_dir_select_species_disease) <- disease_list

for(i in 1:length(select_markers))
{
	for(j in 1:length(disease_list))
	{
		ag_rlm_dir_select_species_disease[i,j] <- ifelse(ag_rlm_qval_select_species_disease[i,j] <= 0.1,3*sign(ag_rlm_est_select_species_disease[i,j]),ifelse(ag_rlm_pval_select_species_disease[i,j] <= 0.05,2*sign(ag_rlm_est_select_species_disease[i,j]),sign(ag_rlm_est_select_species_disease[i,j])))
	}
}

AssociationDirection <- as.data.frame(cbind(apply(ag_rlm_dir_select_species_disease,1,function(x)(length(x[x<= -2]))),apply(ag_rlm_dir_select_species_disease,1,function(x)(length(x[x>=2])))))

G_Markers_ordered <- intersect(G_Markers,rownames(AssociationDirection[order(AssociationDirection[,2]-AssociationDirection[,1]),]))
L_Markers_ordered <- intersect(L_Markers,rownames(AssociationDirection[order(AssociationDirection[,2]-AssociationDirection[,1]),]))

Marker_List_Sorted <- c(setdiff(c(G_Markers_ordered,L_Markers_ordered),c("G_Markers_combined","L_Markers_combined")),c("G_Markers_combined","L_Markers_combined"))

mat_disease_association <- apply(ag_rlm_dir_select_species_disease[Marker_List_Sorted,],2,function(x)(ifelse(x<0,x-0.1,x)))

#heatmap.2(t(mat_disease_association),density="none",trace="none",col=c("skyblue4","skyblue3","powderblue","white","mistyrose1","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(t(mat_disease_association)),rowsep=1:nrow(t(mat_disease_association)),key=FALSE,margins=c(10,5),lhei=c(0.1,5),lwid=c(0.1,5))

ag_mat_disease_association <- mat_disease_association
save(df_ag_controls_diversity_uniqueness,df_ag_diversity_uniqueness,ag_select_age_final_species,ag_combined_df_sum_stat_species,ag_rlm_est_select_species_disease,ag_rlm_pval_select_species_disease,ag_rlm_qval_select_species_disease,ag_mat_disease_association,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\ag_stage2b_results.RData") 

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\AGP\\ag_analysis_2021_Revision.RData")

#rm(list=ls())
