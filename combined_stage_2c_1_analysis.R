print("Loading Libraries")
library(compositions)
library(igraph)
library(metap)
library(robumeta)
library(metafor)
library(dplyr)
library(effsize)
library(MASS)
library(gplots)
library(RColorBrewer)
library(sfsmisc)
library(pcaPP)
library(psych)
library(dendextend)
library(gplots)
library(RColorBrewer)

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_pathway_profiles.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_pathway_profiles.RData")

compute_meta_lm_pathway <- function(data,var1,var2,grouping_variable,grouping_list)
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
			df_temp <- data.frame(var1=vec1,var2=vec2)
			temp_rlm <- rlm(var1~var2,data=df_temp)
			summary_temp_rlm <- summary(temp_rlm)
			temp_meta[i,2] <- summary_temp_rlm$coefficients[2,3]
			temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
			temp_meta[i,4] <- 1
			temp_meta[i,5] <- f.robftest(temp_rlm,var="var2")$p.value
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

compute_meta_lm_group <- function(data,feature_list,metadata_var,grouping_var,grouping_list)
{
	return_out <- as.data.frame(matrix(NA,length(feature_list),10))
	rownames(return_out) <- feature_list
	colnames(return_out) <- c("beta","pval","ci.ub","ci.lb","tau2","QE","QEp","qval","dir","consistency")
	return_out[,1] <- 0
	return_out[,2] <- 1
	return_out[,3] <- 0
	return_out[,4] <- 0
	return_out[,5] <- 0
	return_out[,6] <- 0
	return_out[,7] <- 1
	return_out[,10] <- 0
	
	for(i in 1:length(feature_list))
	{
		species_name <- feature_list[i]
		#print(species_name)
		tryCatch(               
					expr = {                     
						temp_res <- compute_meta_lm_pathway(data,species_name,metadata_var,grouping_var,grouping_list)
						print(species_name)
						return_out[i,"beta"] <- temp_res$model$beta
						return_out[i,"pval"] <- temp_res$model$pval
						return_out[i,"ci.ub"] <- temp_res$model$ci.ub
						return_out[i,"ci.lb"] <- temp_res$model$ci.lb
						return_out[i,"tau2"] <- temp_res$model$tau2
						return_out[i,"QE"] <- temp_res$model$QE
						return_out[i,"QEp"] <- temp_res$model$QEp
						return_out[i,"consistency"] <- length(which(temp_res$df_studies[temp_res$df_studies$ti!=0,"di"]==as.numeric(sign(temp_res$model$b))))/length(temp_res$df_studies[temp_res$df_studies$ti!=0,"di"])
					},
					error = function(e){    
						print(e)
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
	}
	return_out$qval <- p.adjust(return_out$pval,method="fdr")
	return_out$dir <- ifelse(return_out$qval <= 0.1,3*sign(return_out$beta),ifelse(return_out$pval <= 0.08,2*sign(return_out$beta),sign(return_out$beta)))
	return(return_out)
}

combined_shotgun_pathway_profile_clr <- merge(t(cmd3_select_age_final_pathway_clr),t(isc_select_age_final_pathway_clr),by="row.names",all=TRUE)[,-1]
rownames(combined_shotgun_pathway_profile_clr) <- merge(t(cmd3_select_age_final_pathway_clr),t(isc_select_age_final_pathway_clr),by="row.names",all=TRUE)[,1]
combined_shotgun_pathway_profile_clr <- apply(combined_shotgun_pathway_profile_clr,1,function(x)(ifelse(is.na(x),0,x)))

combined_shotgun_pathway_profile_clr <- as.data.frame(combined_shotgun_pathway_profile_clr)
combined_shotgun_pathway_profile_clr$study_name <- NA
combined_shotgun_pathway_profile_clr[intersect(rownames(combined_shotgun_pathway_profile_clr),rownames(combined_df_sum_stat_species)),"study_name"] <- combined_df_sum_stat_species[intersect(rownames(combined_shotgun_pathway_profile_clr),rownames(combined_df_sum_stat_species)),"study_name"]

combined_shotgun_pathway_profile_clr$G1 <- NA
combined_shotgun_pathway_profile_clr[intersect(rownames(combined_shotgun_pathway_profile_clr),rownames(combined_df_sum_stat_species_clr)),"G1"] <- combined_df_sum_stat_species_clr[intersect(rownames(combined_shotgun_pathway_profile_clr),rownames(combined_df_sum_stat_species_clr)),"G1"]

combined_shotgun_pathway_profile_clr$G2 <- NA
combined_shotgun_pathway_profile_clr[intersect(rownames(combined_shotgun_pathway_profile_clr),rownames(combined_df_sum_stat_species_clr)),"G2"] <- combined_df_sum_stat_species_clr[intersect(rownames(combined_shotgun_pathway_profile_clr),rownames(combined_df_sum_stat_species_clr)),"G2"]

select_pathway_study_names <- unique(combined_shotgun_pathway_profile_clr$study_name)[!is.na(unique(combined_shotgun_pathway_profile_clr$study_name))]

select_pathway_study_names <- intersect(sorted_study_rows_shotgun,select_pathway_study_names)

select_pathways <- names(which(apply(compute_detection(combined_shotgun_pathway_profile_clr,colnames(combined_shotgun_pathway_profile_clr)[1:925],"study_name",select_pathway_study_names),1,function(x)(length(x[x>0.05])))>=16))

combined_rem_pathways_G1 <- compute_meta_lm_group(combined_shotgun_pathway_profile_clr,select_pathways,"G1","study_name",select_pathway_study_names)

combined_rem_pathways_G2 <- compute_meta_lm_group(combined_shotgun_pathway_profile_clr,select_pathways,"G2","study_name",select_pathway_study_names)

G2_significant_positive_pathways <- rownames(combined_rem_pathways_G2[(combined_rem_pathways_G2$dir == 3)&(combined_rem_pathways_G2$consistency>=0.65),])
G2_significant_negative_pathways <- rownames(combined_rem_pathways_G2[(combined_rem_pathways_G2$dir == -3)&(combined_rem_pathways_G2$consistency>=0.65),])
G2_positive_pathways <- rownames(combined_rem_pathways_G2[(combined_rem_pathways_G2$dir >= 2)&(combined_rem_pathways_G2$consistency>=0.65),])
G2_negative_pathways <- rownames(combined_rem_pathways_G2[(combined_rem_pathways_G2$dir <= -2)&(combined_rem_pathways_G2$consistency>=0.65),])

G1_significant_positive_pathways <- rownames(combined_rem_pathways_G1[(combined_rem_pathways_G1$dir == 3)&(combined_rem_pathways_G1$consistency>=0.65),])
G1_significant_negative_pathways <- rownames(combined_rem_pathways_G1[(combined_rem_pathways_G1$dir == -3)&(combined_rem_pathways_G1$consistency>=0.65),])
G1_positive_pathways <- rownames(combined_rem_pathways_G1[(combined_rem_pathways_G1$dir >= 2)&(combined_rem_pathways_G1$consistency>=0.65),])
G1_negative_pathways <- rownames(combined_rem_pathways_G1[(combined_rem_pathways_G1$dir <= -2)&(combined_rem_pathways_G1$consistency>=0.65),])
