# The sub-pipeline works on the disease information corresponding to the subjects in the NUAGE data repository and performs the following tasks:
# Computes the association of the relative abundances of the different Highly Detected Species and the different diseases
# Subsequently specifically focusses on the disease associations of the Multiple Disease Enriched and Multiple Disease Depleted (previously identified
# in Ghosh et al eLife, 2020)

# The sub-pipeline utilizes the workspaces nuage_age_analysis.RData already provided in this github and nuage_stage1_results.RData which is generated during the execution
# of run_nuage_analysis_stage1b.R code (which in turn is executed after the run_nuage_stage1a.R code) and the nuage_stage2a_results.RData generated during the execution of 
# run_nuage_analysis_stage2a_1.R
# The workspace nuage_analysis_2021_Revision.RData is a running workspace that is iteratively loaded and saved during the running of each of the previous 
# He et al analysis code (run_nuage_analysis_stage1a.R followed by run_nuage_analysis_stage1b.R followed by run_nuage_stage2a_1.R followed by 
# The detailed order of the running of each individual sub-pipeline is provided in the Readme file in the github corresponding to this study.

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\NUAGE\\nuage_age_analysis.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\nuage_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\nuage_stage2a_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\NUAGE\\nuage_Analysis_2021_Revision.RData")

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


nuage_all_controls <- rownames(df_nuage_diversity_uniqueness)

HighlyDetectedSpecies <- names(which(100*colSums(apply(nuage_select_age_final_species,2,function(x)(ifelse(x>0,1,0))))/nrow(nuage_select_age_final_species)>=5))

nuage_combined_df_sum_stat_species <- as.data.frame(cbind(df_nuage_diversity_uniqueness,nuage_select_age_final_species[rownames(df_nuage_diversity_uniqueness),HighlyDetectedSpecies]))


G_Markers <- intersect(c("Bacteroides_fragilis","Actinomyces_odontolyticus","Bifidobacterium_dentium","Clostridium_clostridioforme","Clostridium_nexile","Clostridium_ramosum","Dialister_invisus","Eggerthella_lenta","Escherichia_coli","Fusobacterium_nucleatum","Granulicatella_adiacens","Lactococcus_lactis","Rothia_mucilaginosa","Ruminococcus_gnavus","Streptococcus_anginosus","Streptococcus_infantis","Streptococcus_salivarius","Streptococcus_sanguinis","Streptococcus_vestibularis","Subdoligranulum_variabile","Veillonella_atypica","Clostridium_asparagiforme","Clostridium_bolteae","Clostridium_citroniae","Clostridium_hathewayi","Clostridium_symbiosum","Streptococcus_parasanguinis","Solobacterium_moorei","Ruminococcus_torques","Streptococcus_mitis","Klebsiella_pneumoniae","Streptococcus_australis","Streptococcus_gordonii","Enterobacter_cloacae"),colnames(nuage_combined_df_sum_stat_species))

L_Markers <- intersect(c("Eubacterium_hallii","Dorea_longicatena","Coprococcus_comes","Coprococcus_catus","Butyrivibrio_crossotus","Bacteroides_uniformis","Alistipes_shahii","Alistipes_indistinctus","Roseburia_hominis","Pseudoflavonifractor_capillosus","Eubacterium_siraeum","Eubacterium_rectale","Coprobacter_fastidiosus","Bifidobacterium_longum","Bifidobacterium_animalis","Barnesiella_intestinihominis","Bacteroides_xylanisolvens","Alistipes_senegalensis","Alistipes_putredinis","Alistipes_onderdonkii","Akkermansia_muciniphila","Faecalibacterium_prausnitzii","Roseburia_inulinivorans","Roseburia_intestinalis"),colnames(nuage_combined_df_sum_stat_species))

nuage_combined_df_sum_stat_species$G_Markers_combined <- rowMeans(apply(nuage_combined_df_sum_stat_species[,G_Markers],2,rank_scale))
nuage_combined_df_sum_stat_species$L_Markers_combined <- rowMeans(apply(nuage_combined_df_sum_stat_species[,L_Markers],2,rank_scale))

HighlyDetectedSpecies <- c(HighlyDetectedSpecies,"G_Markers_combined","L_Markers_combined")

G_Markers <- c(G_Markers,"G_Markers_combined")
L_Markers <- c(L_Markers,"L_Markers_combined")

select_markers <- c(G_Markers,L_Markers)

save(df_nuage_diversity_uniqueness,df_nuage_controls_diversity_uniqueness,nuage_select_age_final_species,nuage_combined_df_sum_stat_species,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\nuage_stage2b_results.RData") 

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\NUAGE\\nuage_analysis_2021_Revision.RData")

rm(list=ls())
