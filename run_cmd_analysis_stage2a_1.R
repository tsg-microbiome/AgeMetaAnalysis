load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_disease_analysis.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")
rm(cmd3_combined_df_sum_stat_genus)
rm(cmd3_combined_df_sum_stat_species)
df_cmd3_diversity_uniqueness <- as.data.frame(df_cmd3_diversity_uniqueness[,c("species_shannon","genus_shannon","age","country","study_name","study_condition","gender","BMI","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","ac_genus_bray_uniqueness","ac_genus_jaccard_uniqueness","ac_genus_aitchison_uniqueness","ac_genus_kendall_uniqueness","ac_species_bray_uniqueness","ac_species_jaccard_uniqueness","ac_species_aitchison_uniqueness","ac_species_kendall_uniqueness","pathway_shannon","pathway_bray_uniqueness","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness","ac_pathway_bray_uniqueness","ac_pathway_jaccard_uniqueness","ac_pathway_aitchison_uniqueness","ac_pathway_kendall_uniqueness")])
df_cmd3_disease_diversity_uniqueness <- as.data.frame(df_cmd3_disease_diversity_uniqueness[,c("species_shannon","genus_shannon","age","country","study_name","study_condition","gender","BMI","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","ac_genus_bray_uniqueness","ac_genus_jaccard_uniqueness","ac_genus_aitchison_uniqueness","ac_genus_kendall_uniqueness","ac_species_bray_uniqueness","ac_species_jaccard_uniqueness","ac_species_aitchison_uniqueness","ac_species_kendall_uniqueness","pathway_shannon","pathway_bray_uniqueness","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness","ac_pathway_bray_uniqueness","ac_pathway_jaccard_uniqueness","ac_pathway_aitchison_uniqueness","ac_pathway_kendall_uniqueness")])


library(robumeta)
library(metafor)
library(dplyr)
library(effsize)
library(MASS)
library(sfsmisc)
library(gplots)
library(RColorBrewer)
library(metap)
library(compositions)

#cmd3_select_age_final_species_clr <- as.matrix(clr(cmd3_select_age_final_species+0.00001))
#cmd3_select_age_final_genus_clr <- as.matrix(clr(cmd3_select_age_final_genus+0.00001))

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
		vec1 <- data[data[,grouping_variable]==group,var1]
		vec2 <- data[data[,grouping_variable]==group,var2]
		#print(paste0(group,",",length(vec1[vec1>0])))
		if((length(vec1[vec1>0]) > 0)&&(length(vec2[vec2>0]) > 0))
		{
			#print(data[data[,grouping_variable]==group,c(var1,var2)])
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
		else
		{
			temp_meta[i,2] <- 0
			temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
			temp_meta[i,4] <- 1
			temp_meta[i,5] <- 1
			temp_meta[i,6] <- 1
		}
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

batch_rlm_grouped <- function(data,metadata,variable_group,metadata_feature,grouping_feature,grouping_list)
{
	df_est <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
	rownames(df_est) <- grouping_list
	colnames(df_est) <- variable_group
	df_p_val <- as.data.frame(matrix(1,length(grouping_list),length(variable_group)))
	rownames(df_p_val) <- grouping_list
	colnames(df_p_val) <- variable_group
	for(i in 1:length(grouping_list))
	{
		study_name <- grouping_list[i]
		#print(study_name)
		study_samples <- rownames(metadata[metadata$study_name == study_name,])
		for(j in 1:length(variable_group))
		{
			species_name <- variable_group[j]
			species_val <- data[study_samples,species_name]
			#print(length(species_val[species_val>0]))
			if(length(species_val[species_val>0])>0)
			{
				temp_rlm <- rlm(data[study_samples,species_name]~metadata[study_samples,metadata_feature])
				summary_temp_rlm <- summary(temp_rlm)
				df_est[i,j] <- summary_temp_rlm$coefficients[2,3]
				df_p_val[i,j] <- f.robftest(temp_rlm)$p.value
			}
		}
		
	}
	df_q_val <- apply(df_p_val,2,p.adjust)
	l_fisher <- p.adjust(apply(df_q_val,2,function(x)(sumlog(x)$p)),method="fdr")
	print("l_fisher generated")
	df_dir <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
	rownames(df_dir) <- grouping_list
	colnames(df_dir) <- variable_group
	for(i in 1:length(grouping_list))
	{
		for(j in 1:length(variable_group))
		{
			df_dir[i,j] <- ifelse(df_q_val[i,j]<=0.10,3*sign(df_est[i,j]),ifelse(df_p_val[i,j]<=0.05,2*sign(df_est[i,j]),1*sign(df_est[i,j])))
		}
	}
	return_list <- list("est"=df_est,"p.value"=df_p_val,"q.value"=df_q_val,"fisher"=l_fisher,"dir"=df_dir)
	return(return_list)
}

batch_rem_grouped <- function(data,metadata,variable_group,metadata_feature,grouping_feature,grouping_list)
{
	common_rows <- intersect(rownames(data),rownames(metadata))
	df_est <- as.data.frame(matrix(0,length(variable_group),5))
	rownames(df_est) <- variable_group
	colnames(df_est) <- c("est","pval","qval","dir","metadata_feature")
	df_est[,1] <- 0
	df_est[,2] <- 1
	df_est[,3] <- 1
	df_est[,4] <- 0
	for(j in 1:length(variable_group))
	{
			species_name <- variable_group[j]
			#print(species_name)
			species_val <- data[common_rows,species_name]
			metadata_val <- metadata[common_rows,metadata_feature]
			metadata_grouping <- metadata[common_rows,grouping_feature]
			print(length(species_val[species_val>0]))
			if(length(species_val[species_val>0])>0)
			{
				tryCatch(               
								expr = {                     
										df_temp <- data.frame(sp=species_val,meta_val=metadata_val,meta_grp=metadata_grouping,row.names=common_rows)
										temp_res <- compute_meta_lm(df_temp,"sp","meta_val","meta_grp",grouping_list)
										#print(temp_res$model)
										df_est[j,1] <- as.numeric(temp_res$model$beta)
										df_est[j,2] <- temp_res$model$pval
		
									},
									error = function(e){         
										print("Error observed. Moving to next")
									},
									finally = {            
										print("finally Executed")
									}
								)
				
			}
		
	}
	df_est[,3] <- p.adjust(df_est[,2],method="fdr")
	df_est[,4] <- ifelse(df_est[,3]<=0.10,3*sign(df_est[,1]),ifelse(df_est[,2]<=0.05,2*sign(df_est[,1]),1*sign(df_est[,1])))
	df_est[,5] <- metadata_feature
	return(df_est)
}

rank_scale=function(x)
{
	x <- rank(x);
	y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
	y <- ifelse(is.nan(y),0,y)
	return(y);
}

SortedStudyList <- c("HMP_2019_ibdmdb","SankaranarayananK_2015","CosteaPI_2017","AsnicarF_2021","HansenLBS_2018","NielsenHB_2014","SchirmerM_2016","WirbelJ_2018","ZellerG_2014","KeohaneDM_2020","QinN_2014","YeZ_2018","QinJ_2012","YachidaS_2019","DhakanDB_2019","GuptaA_2019","BritoIL_2016","PehrssonE_2016","LokmerA_2019","PasolliE_2019","RubelMA_2020","RampelliS_2015")

print("Computing Associations at Species Level")

colnames(cmd3_select_age_final_species_clr)[grep("Hungatella_hathewayi",colnames(cmd3_select_age_final_species_clr))] <- "Clostridium_hathewayi"
colnames(cmd3_select_age_final_species_clr)[grep("Tyzzerella_nexilis",colnames(cmd3_select_age_final_species_clr))]  <- "Clostridium_nexile"
colnames(cmd3_select_age_final_species_clr)[grep("Streptococcus_anginosus_group",colnames(cmd3_select_age_final_species_clr))]  <- "Streptococcus_anginosus"
colnames(cmd3_select_age_final_species_clr)[grep("Enterobacter_cloacae_group",colnames(cmd3_select_age_final_species_clr))]  <- "Enterobacter_cloacae"
colnames(cmd3_select_age_final_species_clr)[grep("Erysipelatoclostridium_ramosum",colnames(cmd3_select_age_final_species_clr))]  <- "Clostridium_ramosum"

colnames(cmd3_select_age_final_species)[grep("Hungatella_hathewayi",colnames(cmd3_select_age_final_species))] <- "Clostridium_hathewayi"
colnames(cmd3_select_age_final_species)[grep("Tyzzerella_nexilis",colnames(cmd3_select_age_final_species))]  <- "Clostridium_nexile"
colnames(cmd3_select_age_final_species)[grep("Streptococcus_anginosus_group",colnames(cmd3_select_age_final_species))]  <- "Streptococcus_anginosus"
colnames(cmd3_select_age_final_species)[grep("Enterobacter_cloacae_group",colnames(cmd3_select_age_final_species))]  <- "Enterobacter_cloacae"
colnames(cmd3_select_age_final_species)[grep("Erysipelatoclostridium_ramosum",colnames(cmd3_select_age_final_species))]  <- "Clostridium_ramosum"

rm(temp_data)
rm(temp_data_clr)

temp_data_clr <- as.data.frame(cbind(cmd3_select_age_final_species_clr[rownames(df_cmd3_diversity_uniqueness),],cmd3_select_age_final_metadata[rownames(df_cmd3_diversity_uniqueness),c("study_condition","study_name")]))
colnames(temp_data_clr)[1638] <- "study_condition"
colnames(temp_data_clr)[1639] <- "study_name"

temp_data <- as.data.frame(cbind(cmd3_select_age_final_species[rownames(df_cmd3_diversity_uniqueness),],cmd3_select_age_final_metadata[rownames(df_cmd3_diversity_uniqueness),c("study_condition","study_name")]))
colnames(temp_data)[1638] <- "study_condition"
colnames(temp_data)[1639] <- "study_name"

DetectionRatesSpecies <- compute_detection(temp_data,colnames(temp_data)[1:1637],"study_name",SortedStudyList)
HighlyDetectedSpecies <- grep("_CAG_",names(which((apply(DetectionRatesSpecies,1,function(x)(length(x[x>=0.05])))/ncol(DetectionRatesSpecies))>=0.5)),value=TRUE,invert=TRUE)

cmd3_clr_relab_sample_correlations_species <- as.data.frame(matrix(NA,nrow(temp_data),2))
rownames(cmd3_clr_relab_sample_correlations_species) <- rownames(temp_data)
colnames(cmd3_clr_relab_sample_correlations_species) <- c("study_name","kendall_correlation")
for(i in 1:nrow(temp_data))
{
	sample_name <- rownames(temp_data)[i]
	study_name <- temp_data[sample_name,"study_name"]
	df_temp <- data.frame(clr=as.numeric(temp_data_clr[sample_name,HighlyDetectedSpecies]),relab=as.numeric(temp_data[sample_name,HighlyDetectedSpecies]))
	cmd3_clr_relab_sample_correlations_species[sample_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	cmd3_clr_relab_sample_correlations_species[sample_name,2] <- temp_data[sample_name,"study_name"]
}

cmd3_clr_relab_species_correlations <- as.data.frame(matrix(NA,length(HighlyDetectedSpecies),length(SortedStudyList)))
rownames(cmd3_clr_relab_species_correlations) <- HighlyDetectedSpecies
colnames(cmd3_clr_relab_species_correlations) <- SortedStudyList
for(i in 1:length(HighlyDetectedSpecies))
{
	species_name <- HighlyDetectedSpecies[i]
	for(j in 1:length(SortedStudyList))
	{
		study_name <- SortedStudyList[j]
		study_samples <- rownames(temp_data[temp_data$study_name==study_name,])
		df_temp <- data.frame(clr=as.numeric(temp_data_clr[study_samples,species_name]),relab=as.numeric(temp_data[study_samples,species_name]))
		cmd3_clr_relab_species_correlations[species_name,study_name] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	}
}


print("performing batch lm (relab)")
print("species bray uniqueness")
cmd3_rlm_detected_species_shannon <- batch_rlm_grouped(temp_data,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_shannon","study_name",SortedStudyList)
print("species bray uniqueness")
cmd3_rlm_detected_species_bray_uniqueness <- batch_rlm_grouped(temp_data,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_bray_uniqueness","study_name",SortedStudyList)
print("species jaccard uniqueness")
cmd3_rlm_detected_species_jaccard_uniqueness <- batch_rlm_grouped(temp_data,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_jaccard_uniqueness","study_name",SortedStudyList)
print("species aitchison uniqueness")
cmd3_rlm_detected_species_aitchison_uniqueness <- batch_rlm_grouped(temp_data,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_aitchison_uniqueness","study_name",SortedStudyList)
print("species kendall uniqueness")
cmd3_rlm_detected_species_kendall_uniqueness <- batch_rlm_grouped(temp_data,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_kendall_uniqueness","study_name",SortedStudyList)

print("performing batch lm (clr)")
print("species bray uniqueness")
cmd3_rlm_detected_species_shannon_clr <- batch_rlm_grouped(temp_data_clr,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_shannon","study_name",SortedStudyList)
print("species bray uniqueness")
cmd3_rlm_detected_species_bray_uniqueness_clr <- batch_rlm_grouped(temp_data_clr,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_bray_uniqueness","study_name",SortedStudyList)
print("species jaccard uniqueness")
cmd3_rlm_detected_species_jaccard_uniqueness_clr <- batch_rlm_grouped(temp_data_clr,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_jaccard_uniqueness","study_name",SortedStudyList)
print("species aitchison uniqueness")
cmd3_rlm_detected_species_aitchison_uniqueness_clr <- batch_rlm_grouped(temp_data_clr,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_aitchison_uniqueness","study_name",SortedStudyList)
print("species kendall uniqueness")
cmd3_rlm_detected_species_kendall_uniqueness_clr <- batch_rlm_grouped(temp_data_clr,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_kendall_uniqueness","study_name",SortedStudyList)

rm(temp_data)
rm(temp_data_clr)

print("performing batch rem (clr)")

cmd3_rem_species_versus_shannon_clr_Overall <- batch_rem_grouped(cmd3_select_age_final_species_clr,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_shannon","study_name",SortedStudyList)
cmd3_rem_species_versus_bray_uniqueness_clr_Overall <- batch_rem_grouped(cmd3_select_age_final_species_clr,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_bray_uniqueness","study_name",SortedStudyList)
cmd3_rem_species_versus_jaccard_uniqueness_clr_Overall <- batch_rem_grouped(cmd3_select_age_final_species_clr,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_jaccard_uniqueness","study_name",SortedStudyList)
cmd3_rem_species_versus_aitchison_uniqueness_clr_Overall <- batch_rem_grouped(cmd3_select_age_final_species_clr,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_aitchison_uniqueness","study_name",SortedStudyList)
cmd3_rem_species_versus_kendall_uniqueness_clr_Overall <- batch_rem_grouped(cmd3_select_age_final_species_clr,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_kendall_uniqueness","study_name",SortedStudyList)

print("performing batch rem (relab)")

cmd3_rem_species_versus_shannon_Overall <- batch_rem_grouped(cmd3_select_age_final_species,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_shannon","study_name",SortedStudyList)
cmd3_rem_species_versus_bray_uniqueness_Overall <- batch_rem_grouped(cmd3_select_age_final_species,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_bray_uniqueness","study_name",SortedStudyList)
cmd3_rem_species_versus_jaccard_uniqueness_Overall <- batch_rem_grouped(cmd3_select_age_final_species,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_jaccard_uniqueness","study_name",SortedStudyList)
cmd3_rem_species_versus_aitchison_uniqueness_Overall <- batch_rem_grouped(cmd3_select_age_final_species,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_aitchison_uniqueness","study_name",SortedStudyList)
cmd3_rem_species_versus_kendall_uniqueness_Overall <- batch_rem_grouped(cmd3_select_age_final_species,df_cmd3_diversity_uniqueness,HighlyDetectedSpecies,"species_kendall_uniqueness","study_name",SortedStudyList)

#cmd3_rem_all_versus_shannon_Overall <- batch_rem_grouped(df_cmd3_diversity_uniqueness,df_cmd3_diversity_uniqueness,c("species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness"),"species_shannon","study_name",SortedStudyList)
#cmd3_rem_all_rest_uniqueness_versus_kendall_uniqueness_Overall <- batch_rem_grouped(df_cmd3_diversity_uniqueness,df_cmd3_diversity_uniqueness,c("species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness"),"species_kendall_uniqueness","study_name",SortedStudyList)

cmd3_df_species_versus_sum_stat_all <- data.frame(shannon=cmd3_rem_species_versus_shannon_Overall[,4],bray_uniqueness=cmd3_rem_species_versus_bray_uniqueness_Overall[,4],jaccard_uniqueness=cmd3_rem_species_versus_jaccard_uniqueness_Overall[,4],aitchison_uniqueness=cmd3_rem_species_versus_aitchison_uniqueness_Overall[,4],kendall_uniqueness=cmd3_rem_species_versus_kendall_uniqueness_Overall[,4],row.names=HighlyDetectedSpecies)

#hmp_species_versus_sum_stat_all <- heatmap.2(apply(cmd3_df_species_versus_sum_stat_all,1,function(x)(ifelse(abs(x)==3,sign(x),0))),density="none",trace="none",col=c("blue3","white","red"),key=FALSE,lwid=c(0.5,5),cexCol=0.5,margins=c(20,10))


#cmd3_taxa_groups <- read.table("C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_taxa_properties.txt",sep="\t",row.names=1,header=TRUE)

cmd3_combined_df_sum_stat_species <- as.data.frame(cbind(df_cmd3_diversity_uniqueness,cmd3_select_age_final_species[rownames(df_cmd3_diversity_uniqueness),HighlyDetectedSpecies]))

cmd3_combined_df_sum_stat_species_clr <- as.data.frame(cbind(df_cmd3_diversity_uniqueness,cmd3_select_age_final_species_clr[rownames(df_cmd3_diversity_uniqueness),HighlyDetectedSpecies]))

#cmd3_direction_detection_species_bray_uniqueness <- as.data.frame(cbind(apply(cmd3_rlm_detected_species_bray_uniqueness$dir,2,function(x)(length(x[x>=2]))),apply(cmd3_rlm_detected_species_bray_uniqueness$dir,2,function(x)(length(x[x<= -2]))),cmd3_rlm_detected_species_bray_uniqueness$fisher))

cmd3_HighlyDetectedSpecies <- HighlyDetectedSpecies

print("Computing Associations at Genus Level")

genus_new_names <- data.frame(new_name = colnames(cmd3_select_age_final_genus),row.names=colnames(cmd3_select_age_final_genus))

genus_new_names[grep("Hungatella",colnames(cmd3_select_age_final_genus)),] <- "Clostridium"
genus_new_names[grep("Tyzzerella",colnames(cmd3_select_age_final_genus)),] <- "Clostridium"
genus_new_names[grep("Erysipelatoclostridium",colnames(cmd3_select_age_final_genus)),] <- "Clostridium"

temp0 <- aggregate(t(cmd3_select_age_final_genus),by=list(genus_new_names[,1]),FUN=sum)[,-1]
rownames(temp0) <- aggregate(t(cmd3_select_age_final_genus),by=list(genus_new_names[,1]),FUN=sum)[,1]
cmd3_select_age_final_genus <- t(temp0)

cmd3_select_age_final_genus_clr <- as.matrix(clr(cmd3_select_age_final_genus+0.00001))

temp_data_clr <- as.data.frame(cbind(cmd3_select_age_final_genus_clr[rownames(df_cmd3_diversity_uniqueness),],cmd3_select_age_final_metadata[rownames(df_cmd3_diversity_uniqueness),c("study_condition","study_name")]))
colnames(temp_data_clr)[429] <- "study_condition"
colnames(temp_data_clr)[430] <- "study_name"

temp_data <- as.data.frame(cbind(cmd3_select_age_final_genus[rownames(df_cmd3_diversity_uniqueness),],cmd3_select_age_final_metadata[rownames(df_cmd3_diversity_uniqueness),c("study_condition","study_name")]))
colnames(temp_data)[429] <- "study_condition"
colnames(temp_data)[430] <- "study_name"

DetectionRatesGenus <- compute_detection(temp_data,colnames(temp_data)[1:428],"study_name",SortedStudyList)
HighlyDetectedGenus <- grep("_CAG_",names(which((apply(DetectionRatesGenus,1,function(x)(length(x[x>=0.05])))/ncol(DetectionRatesGenus))>=0.5)),value=TRUE,invert=TRUE)

cmd3_clr_relab_sample_correlations_genus <- as.data.frame(matrix(NA,nrow(temp_data),2))
rownames(cmd3_clr_relab_sample_correlations_genus) <- rownames(temp_data)
colnames(cmd3_clr_relab_sample_correlations_genus) <- c("study_name","kendall_correlation")
for(i in 1:nrow(temp_data))
{
	sample_name <- rownames(temp_data)[i]
	study_name <- temp_data[sample_name,"study_name"]
	df_temp <- data.frame(clr=as.numeric(temp_data_clr[sample_name,HighlyDetectedGenus]),relab=as.numeric(temp_data[sample_name,HighlyDetectedGenus]))
	cmd3_clr_relab_sample_correlations_genus[sample_name,1] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	cmd3_clr_relab_sample_correlations_genus[sample_name,2] <- temp_data[sample_name,"study_name"]
}

cmd3_clr_relab_genus_correlations <- as.data.frame(matrix(NA,length(HighlyDetectedGenus),length(SortedStudyList)))
rownames(cmd3_clr_relab_genus_correlations) <- HighlyDetectedGenus
colnames(cmd3_clr_relab_genus_correlations) <- SortedStudyList
for(i in 1:length(HighlyDetectedGenus))
{
	genus_name <- HighlyDetectedGenus[i]
	for(j in 1:length(SortedStudyList))
	{
		study_name <- SortedStudyList[j]
		study_samples <- rownames(temp_data[temp_data$study_name==study_name,])
		df_temp <- data.frame(clr=as.numeric(temp_data_clr[study_samples,genus_name]),relab=as.numeric(temp_data[study_samples,genus_name]))
		cmd3_clr_relab_genus_correlations[genus_name,study_name] <- as.numeric(cor(df_temp[(df_temp[,1] != 0)&(df_temp[,2]!=0),])[1,2])
	}
}


print("performing batch lm (relab)")
print("genus bray uniqueness")
cmd3_rlm_detected_genus_shannon <- batch_rlm_grouped(temp_data,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_shannon","study_name",SortedStudyList)
print("genus bray uniqueness")
cmd3_rlm_detected_genus_bray_uniqueness <- batch_rlm_grouped(temp_data,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_bray_uniqueness","study_name",SortedStudyList)
print("genus jaccard uniqueness")
cmd3_rlm_detected_genus_jaccard_uniqueness <- batch_rlm_grouped(temp_data,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_jaccard_uniqueness","study_name",SortedStudyList)
print("genus aitchison uniqueness")
cmd3_rlm_detected_genus_aitchison_uniqueness <- batch_rlm_grouped(temp_data,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_aitchison_uniqueness","study_name",SortedStudyList)
print("genus kendall uniqueness")
cmd3_rlm_detected_genus_kendall_uniqueness <- batch_rlm_grouped(temp_data,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_kendall_uniqueness","study_name",SortedStudyList)

print("performing batch lm (clr)")
print("genus bray uniqueness")
cmd3_rlm_detected_genus_shannon_clr <- batch_rlm_grouped(temp_data_clr,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_shannon","study_name",SortedStudyList)
print("genus bray uniqueness")
cmd3_rlm_detected_genus_bray_uniqueness_clr <- batch_rlm_grouped(temp_data_clr,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_bray_uniqueness","study_name",SortedStudyList)
print("genus jaccard uniqueness")
cmd3_rlm_detected_genus_jaccard_uniqueness_clr <- batch_rlm_grouped(temp_data_clr,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_jaccard_uniqueness","study_name",SortedStudyList)
print("genus aitchison uniqueness")
cmd3_rlm_detected_genus_aitchison_uniqueness_clr <- batch_rlm_grouped(temp_data_clr,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_aitchison_uniqueness","study_name",SortedStudyList)
print("genus kendall uniqueness")
cmd3_rlm_detected_genus_kendall_uniqueness_clr <- batch_rlm_grouped(temp_data_clr,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_kendall_uniqueness","study_name",SortedStudyList)

rm(temp_data)
rm(temp_data_clr)

print("performing batch rem (clr)")

cmd3_rem_genus_versus_shannon_clr_Overall <- batch_rem_grouped(cmd3_select_age_final_genus_clr,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_shannon","study_name",SortedStudyList)
cmd3_rem_genus_versus_bray_uniqueness_clr_Overall <- batch_rem_grouped(cmd3_select_age_final_genus_clr,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_bray_uniqueness","study_name",SortedStudyList)
cmd3_rem_genus_versus_jaccard_uniqueness_clr_Overall <- batch_rem_grouped(cmd3_select_age_final_genus_clr,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_jaccard_uniqueness","study_name",SortedStudyList)
cmd3_rem_genus_versus_aitchison_uniqueness_clr_Overall <- batch_rem_grouped(cmd3_select_age_final_genus_clr,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_aitchison_uniqueness","study_name",SortedStudyList)
cmd3_rem_genus_versus_kendall_uniqueness_clr_Overall <- batch_rem_grouped(cmd3_select_age_final_genus_clr,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_kendall_uniqueness","study_name",SortedStudyList)

print("performing batch rem (relab)")

cmd3_rem_genus_versus_shannon_Overall <- batch_rem_grouped(cmd3_select_age_final_genus,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_shannon","study_name",SortedStudyList)
cmd3_rem_genus_versus_bray_uniqueness_Overall <- batch_rem_grouped(cmd3_select_age_final_genus,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_bray_uniqueness","study_name",SortedStudyList)
cmd3_rem_genus_versus_jaccard_uniqueness_Overall <- batch_rem_grouped(cmd3_select_age_final_genus,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_jaccard_uniqueness","study_name",SortedStudyList)
cmd3_rem_genus_versus_aitchison_uniqueness_Overall <- batch_rem_grouped(cmd3_select_age_final_genus,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_aitchison_uniqueness","study_name",SortedStudyList)
cmd3_rem_genus_versus_kendall_uniqueness_Overall <- batch_rem_grouped(cmd3_select_age_final_genus,df_cmd3_diversity_uniqueness,HighlyDetectedGenus,"genus_kendall_uniqueness","study_name",SortedStudyList)

#cmd3_rem_all_versus_shannon_Overall <- batch_rem_grouped(df_cmd3_diversity_uniqueness,df_cmd3_diversity_uniqueness,c("genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness"),"genus_shannon","study_name",SortedStudyList)
#cmd3_rem_all_rest_uniqueness_versus_kendall_uniqueness_Overall <- batch_rem_grouped(df_cmd3_diversity_uniqueness,df_cmd3_diversity_uniqueness,c("genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness"),"genus_kendall_uniqueness","study_name",SortedStudyList)

cmd3_df_genus_versus_sum_stat_all <- data.frame(shannon=cmd3_rem_genus_versus_shannon_Overall[,4],bray_uniqueness=cmd3_rem_genus_versus_bray_uniqueness_Overall[,4],jaccard_uniqueness=cmd3_rem_genus_versus_jaccard_uniqueness_Overall[,4],aitchison_uniqueness=cmd3_rem_genus_versus_aitchison_uniqueness_Overall[,4],kendall_uniqueness=cmd3_rem_genus_versus_kendall_uniqueness_Overall[,4],row.names=HighlyDetectedGenus)

#hmp_genus_versus_sum_stat_all <- heatmap.2(apply(cmd3_df_genus_versus_sum_stat_all,1,function(x)(ifelse(abs(x)==3,sign(x),0))),density="none",trace="none",col=c("blue3","white","red"),key=FALSE,lwid=c(0.5,5),cexCol=0.5,margins=c(20,10))


#cmd3_taxa_groups <- read.table("C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_taxa_properties.txt",sep="\t",row.names=1,header=TRUE)

cmd3_combined_df_sum_stat_genus <- as.data.frame(cbind(df_cmd3_diversity_uniqueness,cmd3_select_age_final_genus[rownames(df_cmd3_diversity_uniqueness),HighlyDetectedGenus]))

cmd3_combined_df_sum_stat_genus_clr <- as.data.frame(cbind(df_cmd3_diversity_uniqueness,cmd3_select_age_final_genus_clr[rownames(df_cmd3_diversity_uniqueness),HighlyDetectedGenus]))

#cmd3_direction_detection_genus_bray_uniqueness <- as.data.frame(cbind(apply(cmd3_rlm_detected_genus_bray_uniqueness$dir,2,function(x)(length(x[x>=2]))),apply(cmd3_rlm_detected_genus_bray_uniqueness$dir,2,function(x)(length(x[x<= -2]))),cmd3_rlm_detected_genus_bray_uniqueness$fisher))

cmd3_HighlyDetectedGenus <- HighlyDetectedGenus

save(df_cmd3_diversity_uniqueness,cmd3_rlm_detected_species_bray_uniqueness,cmd3_rlm_detected_species_jaccard_uniqueness,cmd3_rlm_detected_species_aitchison_uniqueness,cmd3_rlm_detected_species_kendall_uniqueness,cmd3_rlm_detected_species_shannon,cmd3_combined_df_sum_stat_species,cmd3_combined_df_sum_stat_species_clr,cmd3_HighlyDetectedSpecies,cmd3_clr_relab_sample_correlations_species,cmd3_clr_relab_species_correlations,cmd3_rlm_detected_genus_bray_uniqueness,cmd3_rlm_detected_genus_jaccard_uniqueness,cmd3_rlm_detected_genus_aitchison_uniqueness,cmd3_rlm_detected_genus_kendall_uniqueness,cmd3_rlm_detected_genus_shannon,cmd3_combined_df_sum_stat_genus,cmd3_combined_df_sum_stat_genus_clr,cmd3_HighlyDetectedGenus,cmd3_clr_relab_sample_correlations_genus,cmd3_clr_relab_genus_correlations,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_stage2a_results.RData")

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")

rm(list=ls())







