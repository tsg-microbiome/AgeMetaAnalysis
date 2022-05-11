load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\isc_age_analysis.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_stage2a_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\isc_Analysis_2021_Revision.RData")

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

ibs_samples <- rownames(isc_select_age_final_metadata[(isc_select_age_final_metadata$study_name == "IBS")&(isc_select_age_final_metadata$study_condition == "IBS"),])
control_samples <- rownames(isc_select_age_final_metadata[(isc_select_age_final_metadata$study_name == "IBS")&(isc_select_age_final_metadata$study_condition == "control"),])


HighlyDetectedSpecies <- names(which(100*apply(isc_select_age_final_species[c(ibs_samples,control_samples),],2,function(x)(length(x[x>0])))/length(c(ibs_samples,control_samples))>=5))

isc_combined_df_sum_stat_species <- as.data.frame(cbind(df_isc_diversity_uniqueness,isc_select_age_final_species[rownames(df_isc_diversity_uniqueness),HighlyDetectedSpecies]))

G_Markers <- intersect(c("Bacteroides_fragilis","Actinomyces_odontolyticus","Bifidobacterium_dentium","Clostridium_clostridioforme","Clostridium_nexile","Clostridium_ramosum","Dialister_invisus","Eggerthella_lenta","Escherichia_coli","Fusobacterium_nucleatum","Granulicatella_adiacens","Lactococcus_lactis","Prevotella_copri","Rothia_mucilaginosa","Ruminococcus_gnavus","Streptococcus_anginosus","Streptococcus_infantis","Streptococcus_salivarius","Streptococcus_sanguinis","Streptococcus_vestibularis","Veillonella_atypica","Clostridium_asparagiforme","Clostridiales_bacterium_1_7_47FAA","Clostridium_bolteae","Clostridium_citroniae","Clostridium_hathewayi","Clostridium_symbiosum","Streptococcus_parasanguinis","Solobacterium_moorei","Ruminococcus_torques","Streptococcus_mitis","Klebsiella_pneumoniae","Streptococcus_australis","Streptococcus_gordonii","Enterobacter_cloacae"),colnames(isc_combined_df_sum_stat_species))

L_Markers <- intersect(c("Eubacterium_hallii","Dorea_longicatena","Coprococcus_comes","Coprococcus_catus","Butyrivibrio_crossotus","Bacteroides_uniformis","Alistipes_shahii","Alistipes_indistinctus","Roseburia_hominis","Pseudoflavonifractor_capillosus","Eubacterium_siraeum","Eubacterium_rectale","Coprobacter_fastidiosus","Bifidobacterium_longum","Bifidobacterium_animalis","Barnesiella_intestinihominis","Bacteroides_xylanisolvens","Alistipes_senegalensis","Alistipes_putredinis","Alistipes_onderdonkii","Akkermansia_muciniphila","Faecalibacterium_prausnitzii","Roseburia_inulinivorans","Roseburia_intestinalis"),colnames(isc_combined_df_sum_stat_species))

ibs_samples <- rownames(df_isc_diversity_uniqueness[(df_isc_diversity_uniqueness$study_condition=="IBS")&(df_isc_diversity_uniqueness$study_name=="IBS"),])
control_samples <- rownames(df_isc_diversity_uniqueness[(df_isc_diversity_uniqueness$study_condition=="control")&(df_isc_diversity_uniqueness$study_name=="IBS"),])

disease_list <- "IBS"

select_markers <- c(G_Markers,L_Markers)

isc_combined_df_sum_stat_species$G_Markers_combined <- rowMeans(apply(isc_combined_df_sum_stat_species[,G_Markers],2,rank_scale))
isc_combined_df_sum_stat_species$L_Markers_combined <- rowMeans(apply(isc_combined_df_sum_stat_species[,L_Markers],2,rank_scale))

isc_rlm_est_species_disease <- as.data.frame(matrix(0,length(HighlyDetectedSpecies),length(disease_list)))
rownames(isc_rlm_est_species_disease) <- HighlyDetectedSpecies
colnames(isc_rlm_est_species_disease) <- disease_list

isc_rlm_pval_species_disease <- as.data.frame(matrix(1,length(HighlyDetectedSpecies),length(disease_list)))
rownames(isc_rlm_pval_species_disease) <- HighlyDetectedSpecies
colnames(isc_rlm_pval_species_disease) <- disease_list

for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	for(j in 1:length(disease_list))
	{
		disease_name <- disease_list[j]
		control_list <- control_samples
		print(paste0(species_name,",",disease_name))
		diseased_list <- ibs_samples
		vec_control <- isc_combined_df_sum_stat_species[control_list,species_name]
		vec_disease <- isc_combined_df_sum_stat_species[diseased_list,species_name]
		total_vec <- c(vec_disease,vec_control)
		if(length(total_vec[total_vec > 0])>0)
		{
			temp_rlm <- rlm(isc_combined_df_sum_stat_species[c(diseased_list,control_list),species_name]~as.factor(ifelse(c(diseased_list,control_list) %in% control_list,0,1)))
			summary_temp_rlm <- summary(temp_rlm)
			isc_rlm_est_species_disease[species_name,disease_name] <- summary_temp_rlm$coefficients[2,3]
			isc_rlm_pval_species_disease[species_name,disease_name] <- wilcox.test(vec_disease,vec_control)$p.value
		}
	}
}

isc_rlm_qval_species_disease <- apply(isc_rlm_pval_species_disease,2,function(x)(p.adjust(x,method="fdr")))

isc_rlm_dir_species_disease <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),length(disease_list)))
rownames(isc_rlm_dir_species_disease) <- HighlyDetectedSpecies
colnames(isc_rlm_dir_species_disease) <- disease_list

for(i in 1:length(HighlyDetectedSpecies))
{
	for(j in 1:length(disease_list))
	{
		isc_rlm_dir_species_disease[i,j] <- ifelse(isc_rlm_qval_species_disease[i,j] <= 0.1,3*sign(isc_rlm_est_species_disease[i,j]),ifelse(isc_rlm_pval_species_disease[i,j] <= 0.05,2*sign(isc_rlm_est_species_disease[i,j]),sign(isc_rlm_est_species_disease[i,j])))
	}
}

isc_rlm_est_select_species_disease <- data.frame(ibs=isc_rlm_est_species_disease[c(G_Markers,L_Markers),],row.names=c(G_Markers,L_Markers))
isc_rlm_pval_select_species_disease <- data.frame(ibs=isc_rlm_pval_species_disease[c(G_Markers,L_Markers),],row.names=c(G_Markers,L_Markers))
isc_rlm_qval_select_species_disease <- t(t(apply(isc_rlm_pval_select_species_disease,1,function(x)(p.adjust(x,method="fdr")))))
isc_rlm_dir_select_species_disease <- data.frame(ibs=ifelse(isc_rlm_qval_select_species_disease[,1] <= 0.1,3*sign(isc_rlm_est_select_species_disease[,1]),ifelse(isc_rlm_pval_select_species_disease[,1] <= 0.05,2*sign(isc_rlm_est_select_species_disease[,1]),sign(isc_rlm_est_select_species_disease[,1]))),row.names=select_markers)

mat_disease_association <- apply(isc_rlm_dir_select_species_disease,2,function(x)(ifelse(x<0,x-0.1,x)))

isc_mat_disease_association <- mat_disease_association
save(isc_select_age_final_species,df_isc_diversity_uniqueness,df_isc_controls_diversity_uniqueness,isc_combined_df_sum_stat_species,isc_rlm_est_select_species_disease,isc_rlm_pval_select_species_disease,isc_rlm_qval_select_species_disease,isc_mat_disease_association,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_stage2b_results.RData")
save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\isc_analysis_2021_Revision.RData")

rm(list=ls()) 


