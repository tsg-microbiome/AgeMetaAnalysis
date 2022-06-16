library(robumeta)
library(metafor)
library(dplyr)
library(effsize)
library(MASS)
library(sfsmisc)
library(compositions)
library(igraph)
library(metap)
library(ggplot2)
library(reshape)
library(reshape2)
library(gplots)
library(RColorBrewer)

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
			detection_matrix[i,j] <- length(which(data[data[,grouping_variable]==group,var1]!=0))/length(data[data[,grouping_variable]==group,var1])
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
		#print(group)
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
compute_meta_lm_single_adjust <- function(data,var1,var2,var3,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),6))
	colnames(temp_meta) <- c("dataset","ti","ni","mi","pi","di")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		print(group)
		f <- as.formula(paste0(var1,"~",var3,"+",var2))
		temp_rlm <- rlm(f,data=data[data[,grouping_variable]==group,])
		summary_temp_rlm <- summary(temp_rlm)
		temp_meta[i,2] <- summary_temp_rlm$coefficients[3,3]
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
						temp_res <- compute_meta_lm(data,species_name,metadata_var,grouping_var,grouping_list)
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
	#print("l_fisher generated")
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
	colnames(df_est) <- c("est","pval","qval","dir",metadata_feature)
	for(j in 1:length(variable_group))
	{
			species_name <- variable_group[j]
			#print(species_name)
			species_val <- data[common_rows,species_name]
			metadata_val <- metadata[common_rows,metadata_feature]
			metadata_grouping <- metadata[common_rows,grouping_feature]
			#print(length(species_val[species_val>0]))
			if(length(species_val[species_val>0])>0)
			{
				df_temp <- data.frame(sp=species_val,meta_val=metadata_val,meta_grp=metadata_grouping,row.names=common_rows)
				temp_res <- compute_meta_lm(df_temp,"sp","meta_val","meta_grp",grouping_list)
				print(temp_res$model)
				df_est[j,1] <- as.numeric(temp_res$model$beta)
				df_est[j,2] <- temp_res$model$pval
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

range_scale=function(x)
{
	y <- (x-min(x))/(max(x)-min(x));
	return(y);
}

rem_network=function(data,species_list,group_name,study_list)
{
	species_data <- as.data.frame(t(apply(data[,species_list],1,clr)))
	species_data$group <- data[,group_name]
	source_array <- NULL
	target_array <- NULL
	direction_array <- NULL
	k <- 0
	for(i in 1:length(species_list))
	{
		species1 <- species_list[i]
		#print(species1)
		for(j in 1:length(species_list))
		{
			species2 <- species_list[j]
			#print(paste0(species1,",",species2))
			if(species1 != species2)
			{
				temp_rem <- compute_meta_lm(data,species1,species2,group_name,study_list)
				#print(temp_rem)
				if(temp_rem$model$pval <= 0.05)
				{
					source_array[k] <- species1
					target_array[k] <- species2
					direction_array[k] <- sign(as.numeric(temp_rem$model$beta))
					k=k+1
				}
			}
		}
	}
	return(data.frame(source=source_array,target=target_array,direction=direction_array))
}

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_stage2a_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ag_stage2a_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\he_stage2a_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\odamaki_stage2a_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_stage2a_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\nuage_stage2a_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\logmpie_stage2a_results.RData")

print("Creating Study Groups")

isc_combined_df_sum_stat_species$study_name <- "ISC"
isc_combined_df_sum_stat_species_clr$study_name <- "ISC"

study_rows_shotgun <- c("HMP_2019_ibdmdb","SankaranarayananK_2015","CosteaPI_2017","AsnicarF_2021","HansenLBS_2018","NielsenHB_2014","SchirmerM_2016","WirbelJ_2018","ZellerG_2014","KeohaneDM_2020","QinN_2014","YeZ_2018","QinJ_2012","YachidaS_2019","DhakanDB_2019","GuptaA_2019","BritoIL_2016","PehrssonE_2016","LokmerA_2019","PasolliE_2019","RubelMA_2020","RampelliS_2015","ISC")

shotgun_common_species <- intersect(rownames(isc_rlm_detected_species_bray_uniqueness),colnames(cmd3_rlm_detected_species_bray_uniqueness$est))

sorted_study_rows_shotgun <- c("HMP_2019_ibdmdb","SankaranarayananK_2015","CosteaPI_2017","AsnicarF_2021","HansenLBS_2018","NielsenHB_2014","SchirmerM_2016","WirbelJ_2018","ZellerG_2014","KeohaneDM_2020","ISC","QinN_2014","YeZ_2018","QinJ_2012","YachidaS_2019","DhakanDB_2019","GuptaA_2019","BritoIL_2016","PehrssonE_2016","LokmerA_2019","PasolliE_2019","RubelMA_2020","RampelliS_2015")

EU_NA_studies_shotgun <- c("HMP_2019_ibdmdb","SankaranarayananK_2015","CosteaPI_2017","AsnicarF_2021","HansenLBS_2018","NielsenHB_2014","SchirmerM_2016","WirbelJ_2018","ZellerG_2014","KeohaneDM_2020","ISC")
Est_Asia_studies_shotgun <- c("QinN_2014","YeZ_2018","QinJ_2012","YachidaS_2019")
Other_studies_shotgun <- c("DhakanDB_2019","GuptaA_2019","BritoIL_2016","PehrssonE_2016","LokmerA_2019","PasolliE_2019","RubelMA_2020","RampelliS_2015")

study_rows_full <- c("HMP_2019_ibdmdb","SankaranarayananK_2015","CosteaPI_2017","AsnicarF_2021","HansenLBS_2018","NielsenHB_2014","SchirmerM_2016","WirbelJ_2018","ZellerG_2014","KeohaneDM_2020","ISC","QinN_2014","YeZ_2018","QinJ_2012","YachidaS_2019","DhakanDB_2019","GuptaA_2019","BritoIL_2016","PehrssonE_2016","LokmerA_2019","PasolliE_2019","RubelMA_2020","RampelliS_2015","AG","NUAGE","Odamaki","HE","LogMPie")

sorted_study_rows <- study_rows_full

study_rows_16S <- c("AG","NUAGE","Odamaki","He","LogMPie")

print("Computing CLR v/s RELAB Correlations")

cmd3_temp <- as.data.frame(melt(cmd3_clr_relab_species_correlations))
ag_temp <- ag_clr_relab_species_correlations[,c(2,1)]
colnames(ag_temp) <- c("variable","value")
he_temp <- he_clr_relab_species_correlations[,c(2,1)]
colnames(he_temp) <- c("variable","value")
odamaki_temp <- odamaki_clr_relab_species_correlations[,c(2,1)]
colnames(odamaki_temp) <- c("variable","value")
isc_temp <- isc_clr_relab_species_correlations[,c(2,1)]
colnames(isc_temp) <- c("variable","value")
nuage_temp <- nuage_clr_relab_species_correlations[,c(2,1)]
colnames(nuage_temp) <- c("variable","value")
logmpie_temp <- logmpie_clr_relab_species_correlations[,c(2,1)]
colnames(logmpie_temp) <- c("variable","value")

combined_clr_relab_species_correlation <- as.data.frame(rbind(cmd3_temp,ag_temp,isc_temp,nuage_temp,odamaki_temp,he_temp,logmpie_temp))

#ggplot(combined_clr_relab_genus_correlation[combined_clr_relab_genus_correlation$variable %in% sorted_study_rows,],aes(y=value,x=variable,group=variable))+geom_boxplot()+ylim(0,1)+xlab("")+ylab("")+theme_bw()+theme(axis.text.y=element_text(size=15))

cmd3_temp <- as.data.frame(melt(cmd3_clr_relab_genus_correlations))
ag_temp <- ag_clr_relab_genus_correlations[,c(2,1)]
colnames(ag_temp) <- c("variable","value")
he_temp <- he_clr_relab_genus_correlations[,c(2,1)]
colnames(he_temp) <- c("variable","value")
odamaki_temp <- odamaki_clr_relab_genus_correlations[,c(2,1)]
colnames(odamaki_temp) <- c("variable","value")
isc_temp <- isc_clr_relab_genus_correlations[,c(2,1)]
colnames(isc_temp) <- c("variable","value")
nuage_temp <- nuage_clr_relab_genus_correlations[,c(2,1)]
colnames(nuage_temp) <- c("variable","value")
logmpie_temp <- logmpie_clr_relab_genus_correlations[,c(2,1)]
colnames(logmpie_temp) <- c("variable","value")

#ggplot(combined_clr_relab_species_correlation[combined_clr_relab_species_correlation$variable %in% sorted_study_rows,],aes(y=value,x=variable,group=variable))+geom_boxplot()+ylim(0,1)+xlab("")+ylab("")+theme_bw()+theme(axis.text.y=element_text(size=15))

combined_clr_relab_genus_correlation <- as.data.frame(rbind(cmd3_temp,ag_temp,isc_temp,nuage_temp,odamaki_temp,he_temp,logmpie_temp))

colnames(cmd3_clr_relab_sample_correlations_species) <- c("correlation","study_name")
colnames(ag_clr_relab_sample_correlations_species) <- c("correlation","study_name")
colnames(isc_clr_relab_sample_correlations_species) <- c("correlation","study_name")
colnames(nuage_clr_relab_sample_correlations_species) <- c("correlation","study_name")
colnames(odamaki_clr_relab_sample_correlations_species) <- c("correlation","study_name")
colnames(he_clr_relab_sample_correlations_species) <- c("correlation","study_name")
colnames(logmpie_clr_relab_sample_correlations_species) <- c("correlation","study_name")

combined_clr_relab_sample_correlations_species <- as.data.frame(rbind(cmd3_clr_relab_sample_correlations_species,ag_clr_relab_sample_correlations_species,isc_clr_relab_sample_correlations_species,nuage_clr_relab_sample_correlations_species,odamaki_clr_relab_sample_correlations_species,he_clr_relab_sample_correlations_species,logmpie_clr_relab_sample_correlations_species))

combined_clr_relab_sample_correlations_species <- combined_clr_relab_sample_correlations_species[combined_clr_relab_sample_correlations_species$study_name %in% study_rows_full,]

combined_clr_relab_sample_correlations_species$study_name <- factor(combined_clr_relab_sample_correlations_species$study_name,levels=sorted_study_rows)
ggplot(combined_clr_relab_sample_correlations_species[combined_clr_relab_sample_correlations_species$study_name %in% sorted_study_rows,],aes(y=correlation,x=study_name,group=study_name))+geom_boxplot()+ylim(0,1)+xlab("")+ylab("")+theme_bw()+theme(axis.text.y=element_text(size=15))

colnames(cmd3_clr_relab_sample_correlations_genus) <- c("correlation","study_name")
colnames(ag_clr_relab_sample_correlations_genus) <- c("correlation","study_name")
colnames(isc_clr_relab_sample_correlations_genus) <- c("correlation","study_name")
colnames(nuage_clr_relab_sample_correlations_genus) <- c("correlation","study_name")
colnames(odamaki_clr_relab_sample_correlations_genus) <- c("correlation","study_name")
colnames(he_clr_relab_sample_correlations_genus) <- c("correlation","study_name")
colnames(logmpie_clr_relab_sample_correlations_genus) <- c("correlation","study_name")

combined_clr_relab_sample_correlations_genus <- as.data.frame(rbind(cmd3_clr_relab_sample_correlations_genus,ag_clr_relab_sample_correlations_genus,isc_clr_relab_sample_correlations_genus,nuage_clr_relab_sample_correlations_genus,odamaki_clr_relab_sample_correlations_genus,he_clr_relab_sample_correlations_genus,logmpie_clr_relab_sample_correlations_genus))

combined_clr_relab_sample_correlations_genus$study_name <- factor(combined_clr_relab_sample_correlations_genus$study_name,levels=sorted_study_rows)
ggplot(combined_clr_relab_sample_correlations_genus[combined_clr_relab_sample_correlations_genus$study_name %in% sorted_study_rows,],aes(y=correlation,x=study_name,group=study_name))+geom_boxplot()+ylim(0,1)+xlab("")+ylab("")+theme_bw()+theme(axis.text.y=element_text(size=15))


print("Shotgun Analysis")

print("preparing Species Table")

cmd3_summary_featuress <- c("species_shannon","age","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

cmd3_species_features <- colnames(cmd3_combined_df_sum_stat_species)[35:225]

isc_summary_features <- c("species_shannon","age","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")

isc_species_features <- grep("unclassified",colnames(isc_combined_df_sum_stat_species)[14:173],value=TRUE,invert=TRUE)

cmd3_temp_df_species <- as.data.frame(cmd3_combined_df_sum_stat_species[,cmd3_species_features])
isc_temp_df_species <- as.data.frame(isc_combined_df_sum_stat_species[,isc_species_features])

shotgun_combined_df_sum_stat_species <- merge(t(cmd3_temp_df_species),t(isc_temp_df_species),by="row.names",all=TRUE)[,-1]
rownames(shotgun_combined_df_sum_stat_species) <- merge(t(cmd3_temp_df_species),t(isc_temp_df_species),by="row.names",all=TRUE)[,1]
shotgun_combined_df_sum_stat_species <- as.data.frame(apply(shotgun_combined_df_sum_stat_species,1,function(x)(ifelse(is.na(x),0,x))))
shotgun_combined_df_sum_stat_species_clr <- as.data.frame(as.matrix(clr(shotgun_combined_df_sum_stat_species+0.00001)))
shotgun_combined_df_sum_stat_species_clr <- as.data.frame(t(apply(shotgun_combined_df_sum_stat_species_clr,1,function(x)(x-min(x)))))

shotgun_combined_df_sum_stat_species$species_shannon <- NA 
shotgun_combined_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"species_shannon"] <- cmd3_combined_df_sum_stat_species[,"species_shannon"]
shotgun_combined_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"species_shannon"] <- isc_combined_df_sum_stat_species[,"species_shannon"]

shotgun_combined_df_sum_stat_species_clr$species_shannon <- NA 
shotgun_combined_df_sum_stat_species_clr[rownames(cmd3_combined_df_sum_stat_species_clr),"species_shannon"] <- cmd3_combined_df_sum_stat_species_clr[,"species_shannon"]
shotgun_combined_df_sum_stat_species_clr[rownames(isc_combined_df_sum_stat_species_clr),"species_shannon"] <- isc_combined_df_sum_stat_species_clr[,"species_shannon"]

shotgun_combined_df_sum_stat_species$age <- NA 
shotgun_combined_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"age"] <- cmd3_combined_df_sum_stat_species[,"age"]
shotgun_combined_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"age"] <- isc_combined_df_sum_stat_species[,"age"]

shotgun_combined_df_sum_stat_species_clr$age <- NA 
shotgun_combined_df_sum_stat_species_clr[rownames(cmd3_combined_df_sum_stat_species_clr),"age"] <- cmd3_combined_df_sum_stat_species_clr[,"age"]
shotgun_combined_df_sum_stat_species_clr[rownames(isc_combined_df_sum_stat_species_clr),"age"] <- isc_combined_df_sum_stat_species_clr[,"age"]

shotgun_combined_df_sum_stat_species$species_bray_uniqueness <- NA 
shotgun_combined_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"species_bray_uniqueness"] <- cmd3_combined_df_sum_stat_species[,"species_bray_uniqueness"]
shotgun_combined_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"species_bray_uniqueness"] <- isc_combined_df_sum_stat_species[,"species_bray_uniqueness"]

shotgun_combined_df_sum_stat_species_clr$species_bray_uniqueness <- NA 
shotgun_combined_df_sum_stat_species_clr[rownames(cmd3_combined_df_sum_stat_species_clr),"species_bray_uniqueness"] <- cmd3_combined_df_sum_stat_species_clr[,"species_bray_uniqueness"]
shotgun_combined_df_sum_stat_species_clr[rownames(isc_combined_df_sum_stat_species_clr),"species_bray_uniqueness"] <- isc_combined_df_sum_stat_species_clr[,"species_bray_uniqueness"]

shotgun_combined_df_sum_stat_species$species_jaccard_uniqueness <- NA 
shotgun_combined_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"species_jaccard_uniqueness"] <- cmd3_combined_df_sum_stat_species[,"species_jaccard_uniqueness"]
shotgun_combined_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"species_jaccard_uniqueness"] <- isc_combined_df_sum_stat_species[,"species_jaccard_uniqueness"]

shotgun_combined_df_sum_stat_species_clr$species_jaccard_uniqueness <- NA 
shotgun_combined_df_sum_stat_species_clr[rownames(cmd3_combined_df_sum_stat_species_clr),"species_jaccard_uniqueness"] <- cmd3_combined_df_sum_stat_species_clr[,"species_jaccard_uniqueness"]
shotgun_combined_df_sum_stat_species_clr[rownames(isc_combined_df_sum_stat_species_clr),"species_jaccard_uniqueness"] <- isc_combined_df_sum_stat_species_clr[,"species_jaccard_uniqueness"]

shotgun_combined_df_sum_stat_species$species_aitchison_uniqueness <- NA 
shotgun_combined_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"species_aitchison_uniqueness"] <- cmd3_combined_df_sum_stat_species[,"species_aitchison_uniqueness"]
shotgun_combined_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"species_aitchison_uniqueness"] <- isc_combined_df_sum_stat_species[,"species_aitchison_uniqueness"]

shotgun_combined_df_sum_stat_species_clr$species_aitchison_uniqueness <- NA 
shotgun_combined_df_sum_stat_species_clr[rownames(cmd3_combined_df_sum_stat_species_clr),"species_aitchison_uniqueness"] <- cmd3_combined_df_sum_stat_species_clr[,"species_aitchison_uniqueness"]
shotgun_combined_df_sum_stat_species_clr[rownames(isc_combined_df_sum_stat_species_clr),"species_aitchison_uniqueness"] <- isc_combined_df_sum_stat_species_clr[,"species_aitchison_uniqueness"]

shotgun_combined_df_sum_stat_species$species_kendall_uniqueness <- NA 
shotgun_combined_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"species_kendall_uniqueness"] <- cmd3_combined_df_sum_stat_species[,"species_kendall_uniqueness"]
shotgun_combined_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"species_kendall_uniqueness"] <- isc_combined_df_sum_stat_species[,"species_kendall_uniqueness"]

shotgun_combined_df_sum_stat_species_clr$species_kendall_uniqueness <- NA 
shotgun_combined_df_sum_stat_species_clr[rownames(cmd3_combined_df_sum_stat_species_clr),"species_kendall_uniqueness"] <- cmd3_combined_df_sum_stat_species_clr[,"species_kendall_uniqueness"]
shotgun_combined_df_sum_stat_species_clr[rownames(isc_combined_df_sum_stat_species_clr),"species_kendall_uniqueness"] <- isc_combined_df_sum_stat_species_clr[,"species_kendall_uniqueness"]

shotgun_combined_df_sum_stat_species$study_name <- NA 
shotgun_combined_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"study_name"] <- cmd3_combined_df_sum_stat_species[,"study_name"]
shotgun_combined_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"study_name"] <- isc_combined_df_sum_stat_species[,"study_name"]

shotgun_combined_df_sum_stat_species_clr$study_name <- NA 
shotgun_combined_df_sum_stat_species_clr[rownames(cmd3_combined_df_sum_stat_species_clr),"study_name"] <- cmd3_combined_df_sum_stat_species_clr[,"study_name"]
shotgun_combined_df_sum_stat_species_clr[rownames(isc_combined_df_sum_stat_species_clr),"study_name"] <- isc_combined_df_sum_stat_species_clr[,"study_name"]

shotgun_combined_df_sum_stat_species_clr$study_condition <- NA 
shotgun_combined_df_sum_stat_species_clr[rownames(cmd3_combined_df_sum_stat_species_clr),"study_condition"] <- cmd3_combined_df_sum_stat_species_clr[,"study_condition"]
shotgun_combined_df_sum_stat_species_clr[rownames(isc_combined_df_sum_stat_species_clr),"study_condition"] <- isc_combined_df_sum_stat_species_clr[,"study_condition"]

shotgun_combined_df_sum_stat_species$study_condition <- NA 
shotgun_combined_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"study_condition"] <- cmd3_combined_df_sum_stat_species[,"study_condition"]
shotgun_combined_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"study_condition"] <- isc_combined_df_sum_stat_species[,"study_condition"]

print("preparing Genus Table")

cmd3_summary_features <- c("genus_shannon","age","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness")

cmd3_genus_features <- colnames(cmd3_combined_df_sum_stat_genus)[34:129]

isc_genus_features <- c("genus_shannon","age","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness")

isc_genus_features <- grep("unclassified",colnames(isc_combined_df_sum_stat_genus)[14:68],value=TRUE,invert=TRUE)

shotgun_combined_df_sum_stat_genus <- merge(t(cmd3_combined_df_sum_stat_genus[,cmd3_genus_features]),t(isc_combined_df_sum_stat_genus[,isc_genus_features]),by="row.names",all=TRUE)[,-1]
rownames(shotgun_combined_df_sum_stat_genus) <- merge(t(cmd3_combined_df_sum_stat_genus[,cmd3_genus_features]),t(isc_combined_df_sum_stat_genus[,isc_genus_features]),by="row.names",all=TRUE)[,1]
shotgun_combined_df_sum_stat_genus <- as.data.frame(apply(shotgun_combined_df_sum_stat_genus,1,function(x)(ifelse(is.na(x),0,x))))
shotgun_combined_df_sum_stat_genus_clr <- as.data.frame(as.matrix(clr(shotgun_combined_df_sum_stat_genus+0.00001)))
shotgun_combined_df_sum_stat_genus_clr <- as.data.frame(t(apply(shotgun_combined_df_sum_stat_genus_clr,1,function(x)(x-min(x)))))

shotgun_combined_df_sum_stat_genus$genus_shannon <- NA 
shotgun_combined_df_sum_stat_genus[rownames(cmd3_combined_df_sum_stat_genus),"genus_shannon"] <- cmd3_combined_df_sum_stat_genus[,"genus_shannon"]
shotgun_combined_df_sum_stat_genus[rownames(isc_combined_df_sum_stat_genus),"genus_shannon"] <- isc_combined_df_sum_stat_genus[,"genus_shannon"]

shotgun_combined_df_sum_stat_genus_clr$genus_shannon <- NA 
shotgun_combined_df_sum_stat_genus_clr[rownames(cmd3_combined_df_sum_stat_genus_clr),"genus_shannon"] <- cmd3_combined_df_sum_stat_genus_clr[,"genus_shannon"]
shotgun_combined_df_sum_stat_genus_clr[rownames(isc_combined_df_sum_stat_genus_clr),"genus_shannon"] <- isc_combined_df_sum_stat_genus_clr[,"genus_shannon"]

shotgun_combined_df_sum_stat_genus$age <- NA 
shotgun_combined_df_sum_stat_genus[rownames(cmd3_combined_df_sum_stat_genus),"age"] <- cmd3_combined_df_sum_stat_genus[,"age"]
shotgun_combined_df_sum_stat_genus[rownames(isc_combined_df_sum_stat_genus),"age"] <- isc_combined_df_sum_stat_genus[,"age"]

shotgun_combined_df_sum_stat_genus_clr$age <- NA 
shotgun_combined_df_sum_stat_genus_clr[rownames(cmd3_combined_df_sum_stat_genus_clr),"age"] <- cmd3_combined_df_sum_stat_genus_clr[,"age"]
shotgun_combined_df_sum_stat_genus_clr[rownames(isc_combined_df_sum_stat_genus_clr),"age"] <- isc_combined_df_sum_stat_genus_clr[,"age"]

shotgun_combined_df_sum_stat_genus$genus_bray_uniqueness <- NA 
shotgun_combined_df_sum_stat_genus[rownames(cmd3_combined_df_sum_stat_genus),"genus_bray_uniqueness"] <- cmd3_combined_df_sum_stat_genus[,"genus_bray_uniqueness"]
shotgun_combined_df_sum_stat_genus[rownames(isc_combined_df_sum_stat_genus),"genus_bray_uniqueness"] <- isc_combined_df_sum_stat_genus[,"genus_bray_uniqueness"]

shotgun_combined_df_sum_stat_genus_clr$genus_bray_uniqueness <- NA 
shotgun_combined_df_sum_stat_genus_clr[rownames(cmd3_combined_df_sum_stat_genus_clr),"genus_bray_uniqueness"] <- cmd3_combined_df_sum_stat_genus_clr[,"genus_bray_uniqueness"]
shotgun_combined_df_sum_stat_genus_clr[rownames(isc_combined_df_sum_stat_genus_clr),"genus_bray_uniqueness"] <- isc_combined_df_sum_stat_genus_clr[,"genus_bray_uniqueness"]

shotgun_combined_df_sum_stat_genus$genus_jaccard_uniqueness <- NA 
shotgun_combined_df_sum_stat_genus[rownames(cmd3_combined_df_sum_stat_genus),"genus_jaccard_uniqueness"] <- cmd3_combined_df_sum_stat_genus[,"genus_jaccard_uniqueness"]
shotgun_combined_df_sum_stat_genus[rownames(isc_combined_df_sum_stat_genus),"genus_jaccard_uniqueness"] <- isc_combined_df_sum_stat_genus[,"genus_jaccard_uniqueness"]

shotgun_combined_df_sum_stat_genus_clr$genus_jaccard_uniqueness <- NA 
shotgun_combined_df_sum_stat_genus_clr[rownames(cmd3_combined_df_sum_stat_genus_clr),"genus_jaccard_uniqueness"] <- cmd3_combined_df_sum_stat_genus_clr[,"genus_jaccard_uniqueness"]
shotgun_combined_df_sum_stat_genus_clr[rownames(isc_combined_df_sum_stat_genus_clr),"genus_jaccard_uniqueness"] <- isc_combined_df_sum_stat_genus_clr[,"genus_jaccard_uniqueness"]

shotgun_combined_df_sum_stat_genus$genus_aitchison_uniqueness <- NA 
shotgun_combined_df_sum_stat_genus[rownames(cmd3_combined_df_sum_stat_genus),"genus_aitchison_uniqueness"] <- cmd3_combined_df_sum_stat_genus[,"genus_aitchison_uniqueness"]
shotgun_combined_df_sum_stat_genus[rownames(isc_combined_df_sum_stat_genus),"genus_aitchison_uniqueness"] <- isc_combined_df_sum_stat_genus[,"genus_aitchison_uniqueness"]

shotgun_combined_df_sum_stat_genus_clr$genus_aitchison_uniqueness <- NA 
shotgun_combined_df_sum_stat_genus_clr[rownames(cmd3_combined_df_sum_stat_genus_clr),"genus_aitchison_uniqueness"] <- cmd3_combined_df_sum_stat_genus_clr[,"genus_aitchison_uniqueness"]
shotgun_combined_df_sum_stat_genus_clr[rownames(isc_combined_df_sum_stat_genus_clr),"genus_aitchison_uniqueness"] <- isc_combined_df_sum_stat_genus_clr[,"genus_aitchison_uniqueness"]

shotgun_combined_df_sum_stat_genus$genus_kendall_uniqueness <- NA 
shotgun_combined_df_sum_stat_genus[rownames(cmd3_combined_df_sum_stat_genus),"genus_kendall_uniqueness"] <- cmd3_combined_df_sum_stat_genus[,"genus_kendall_uniqueness"]
shotgun_combined_df_sum_stat_genus[rownames(isc_combined_df_sum_stat_genus),"genus_kendall_uniqueness"] <- isc_combined_df_sum_stat_genus[,"genus_kendall_uniqueness"]

shotgun_combined_df_sum_stat_genus_clr$genus_kendall_uniqueness <- NA 
shotgun_combined_df_sum_stat_genus_clr[rownames(cmd3_combined_df_sum_stat_genus_clr),"genus_kendall_uniqueness"] <- cmd3_combined_df_sum_stat_genus_clr[,"genus_kendall_uniqueness"]
shotgun_combined_df_sum_stat_genus_clr[rownames(isc_combined_df_sum_stat_genus_clr),"genus_kendall_uniqueness"] <- isc_combined_df_sum_stat_genus_clr[,"genus_kendall_uniqueness"]

shotgun_combined_df_sum_stat_genus$study_name <- NA 
shotgun_combined_df_sum_stat_genus[rownames(cmd3_combined_df_sum_stat_genus),"study_name"] <- cmd3_combined_df_sum_stat_genus[,"study_name"]
shotgun_combined_df_sum_stat_genus[rownames(isc_combined_df_sum_stat_genus),"study_name"] <- isc_combined_df_sum_stat_genus[,"study_name"]

shotgun_combined_df_sum_stat_genus_clr$study_name <- NA 
shotgun_combined_df_sum_stat_genus_clr[rownames(cmd3_combined_df_sum_stat_genus_clr),"study_name"] <- cmd3_combined_df_sum_stat_genus_clr[,"study_name"]
shotgun_combined_df_sum_stat_genus_clr[rownames(isc_combined_df_sum_stat_genus_clr),"study_name"] <- isc_combined_df_sum_stat_genus_clr[,"study_name"]

shotgun_combined_df_sum_stat_genus_clr$study_condition <- NA 
shotgun_combined_df_sum_stat_genus_clr[rownames(cmd3_combined_df_sum_stat_genus_clr),"study_condition"] <- cmd3_combined_df_sum_stat_genus_clr[,"study_condition"]
shotgun_combined_df_sum_stat_genus_clr[rownames(isc_combined_df_sum_stat_genus_clr),"study_condition"] <- isc_combined_df_sum_stat_genus_clr[,"study_condition"]

shotgun_combined_df_sum_stat_genus$study_condition <- NA 
shotgun_combined_df_sum_stat_genus[rownames(cmd3_combined_df_sum_stat_genus),"study_condition"] <- cmd3_combined_df_sum_stat_genus[,"study_condition"]
shotgun_combined_df_sum_stat_genus[rownames(isc_combined_df_sum_stat_genus),"study_condition"] <- isc_combined_df_sum_stat_genus[,"study_condition"]

print("Doing RLM Shotgun Species")
select_shotgun_species <- setdiff(colnames(shotgun_combined_df_sum_stat_species_clr),c("species_shannon","age","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","study_condition","study_name"))

print("RLM: species_bray_uniqueness")
shotgun_rlm_est_species_bray_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_est_species_bray_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_est_species_bray_uniqueness_clr) <- sorted_study_rows_shotgun
shotgun_rlm_p_val_species_bray_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_p_val_species_bray_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_p_val_species_bray_uniqueness_clr) <- sorted_study_rows_shotgun

for(i in 1:length(select_shotgun_species))
{
	species_name <- select_shotgun_species[i]
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		study_name <- sorted_study_rows_shotgun[j]
		study_samples <- rownames(shotgun_combined_df_sum_stat_species_clr[shotgun_combined_df_sum_stat_species_clr$study_name == study_name,])
		tryCatch(               
					expr = {  
						print(paste0("species_bray_uniqueness, ",species_name))
						df_temp <- as.data.frame(shotgun_combined_df_sum_stat_species_clr[study_samples,c(species_name,"species_bray_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						shotgun_rlm_est_species_bray_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						shotgun_rlm_p_val_species_bray_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){  
						#print(e)
						print(paste0("species_bray_uniqueness, ",species_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}
shotgun_rlm_est_species_bray_uniqueness_clr <- apply(shotgun_rlm_est_species_bray_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
shotgun_rlm_q_val_species_bray_uniqueness_clr <- apply(shotgun_rlm_p_val_species_bray_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
shotgun_rlm_dir_species_bray_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_dir_species_bray_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_dir_species_bray_uniqueness_clr) <- sorted_study_rows_shotgun
for(i in 1:length(select_shotgun_species))
{
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		shotgun_rlm_dir_species_bray_uniqueness_clr[i,j] <- ifelse(shotgun_rlm_q_val_species_bray_uniqueness_clr[i,j]<=0.1,3*sign(shotgun_rlm_est_species_bray_uniqueness_clr[i,j]),ifelse(shotgun_rlm_p_val_species_bray_uniqueness_clr[i,j]<=0.05,2*sign(shotgun_rlm_est_species_bray_uniqueness_clr[i,j]),sign(shotgun_rlm_est_species_bray_uniqueness_clr[i,j])))
	}
}

print("RLM: species_jaccard_uniqueness")
shotgun_rlm_est_species_jaccard_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_est_species_jaccard_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_est_species_jaccard_uniqueness_clr) <- sorted_study_rows_shotgun
shotgun_rlm_p_val_species_jaccard_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_p_val_species_jaccard_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_p_val_species_jaccard_uniqueness_clr) <- sorted_study_rows_shotgun

for(i in 1:length(select_shotgun_species))
{
	species_name <- select_shotgun_species[i]
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		study_name <- sorted_study_rows_shotgun[j]
		study_samples <- rownames(shotgun_combined_df_sum_stat_species_clr[shotgun_combined_df_sum_stat_species_clr$study_name == study_name,])
		tryCatch(               
					expr = { 
						print(paste0("species_jaccard_uniqueness, ",species_name))
						df_temp <- as.data.frame(shotgun_combined_df_sum_stat_species_clr[study_samples,c(species_name,"species_jaccard_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						shotgun_rlm_est_species_jaccard_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						shotgun_rlm_p_val_species_jaccard_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){         
						print(paste0("species_jaccard_uniqueness, ",species_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}
shotgun_rlm_est_species_jaccard_uniqueness_clr <- apply(shotgun_rlm_est_species_jaccard_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
shotgun_rlm_q_val_species_jaccard_uniqueness_clr <- apply(shotgun_rlm_p_val_species_jaccard_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
shotgun_rlm_dir_species_jaccard_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_dir_species_jaccard_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_dir_species_jaccard_uniqueness_clr) <- sorted_study_rows_shotgun
for(i in 1:length(select_shotgun_species))
{
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		shotgun_rlm_dir_species_jaccard_uniqueness_clr[i,j] <- ifelse(shotgun_rlm_q_val_species_jaccard_uniqueness_clr[i,j]<=0.1,3*sign(shotgun_rlm_est_species_jaccard_uniqueness_clr[i,j]),ifelse(shotgun_rlm_p_val_species_jaccard_uniqueness_clr[i,j]<=0.05,2*sign(shotgun_rlm_est_species_jaccard_uniqueness_clr[i,j]),sign(shotgun_rlm_est_species_jaccard_uniqueness_clr[i,j])))
	}
}

print("RLM: species_aitchison_uniqueness")
shotgun_rlm_est_species_aitchison_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_est_species_aitchison_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_est_species_aitchison_uniqueness_clr) <- sorted_study_rows_shotgun
shotgun_rlm_p_val_species_aitchison_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_p_val_species_aitchison_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_p_val_species_aitchison_uniqueness_clr) <- sorted_study_rows_shotgun

for(i in 1:length(select_shotgun_species))
{
	species_name <- select_shotgun_species[i]
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		study_name <- sorted_study_rows_shotgun[j]
		study_samples <- rownames(shotgun_combined_df_sum_stat_species_clr[shotgun_combined_df_sum_stat_species_clr$study_name == study_name,])
		tryCatch(               
					expr = { 
						print(paste0("species_aitchison_uniqueness, ",species_name))
						df_temp <- as.data.frame(shotgun_combined_df_sum_stat_species_clr[study_samples,c(species_name,"species_aitchison_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						shotgun_rlm_est_species_aitchison_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						shotgun_rlm_p_val_species_aitchison_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){  
						print(paste0("species_aitchison_uniqueness, ",species_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}
shotgun_rlm_est_species_aitchison_uniqueness_clr <- apply(shotgun_rlm_est_species_aitchison_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
shotgun_rlm_q_val_species_aitchison_uniqueness_clr <- apply(shotgun_rlm_p_val_species_aitchison_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
shotgun_rlm_dir_species_aitchison_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_dir_species_aitchison_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_dir_species_aitchison_uniqueness_clr) <- sorted_study_rows_shotgun
for(i in 1:length(select_shotgun_species))
{
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		shotgun_rlm_dir_species_aitchison_uniqueness_clr[i,j] <- ifelse(shotgun_rlm_q_val_species_aitchison_uniqueness_clr[i,j]<=0.1,3*sign(shotgun_rlm_est_species_aitchison_uniqueness_clr[i,j]),ifelse(shotgun_rlm_p_val_species_aitchison_uniqueness_clr[i,j]<=0.05,2*sign(shotgun_rlm_est_species_aitchison_uniqueness_clr[i,j]),sign(shotgun_rlm_est_species_aitchison_uniqueness_clr[i,j])))
	}
}

print("RLM: species_kendall_uniqueness")
shotgun_rlm_est_species_kendall_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_est_species_kendall_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_est_species_kendall_uniqueness_clr) <- sorted_study_rows_shotgun
shotgun_rlm_p_val_species_kendall_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_p_val_species_kendall_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_p_val_species_kendall_uniqueness_clr) <- sorted_study_rows_shotgun

for(i in 1:length(select_shotgun_species))
{
	species_name <- select_shotgun_species[i]
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		study_name <- sorted_study_rows_shotgun[j]
		study_samples <- rownames(shotgun_combined_df_sum_stat_species_clr[shotgun_combined_df_sum_stat_species_clr$study_name == study_name,])
		tryCatch(               
					expr = {    
						print(paste0("species_kendall_uniqueness, ",species_name))
						df_temp <- as.data.frame(shotgun_combined_df_sum_stat_species_clr[study_samples,c(species_name,"species_kendall_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						shotgun_rlm_est_species_kendall_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						shotgun_rlm_p_val_species_kendall_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){   
						print(paste0("species_kendall_uniqueness, ",species_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}
shotgun_rlm_est_species_kendall_uniqueness_clr <- apply(shotgun_rlm_est_species_kendall_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
shotgun_rlm_q_val_species_kendall_uniqueness_clr <- apply(shotgun_rlm_p_val_species_kendall_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
shotgun_rlm_dir_species_kendall_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_dir_species_kendall_uniqueness_clr) <- select_shotgun_species
colnames(shotgun_rlm_dir_species_kendall_uniqueness_clr) <- sorted_study_rows_shotgun
for(i in 1:length(select_shotgun_species))
{
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		shotgun_rlm_dir_species_kendall_uniqueness_clr[i,j] <- ifelse(shotgun_rlm_q_val_species_kendall_uniqueness_clr[i,j]<=0.1,3*sign(shotgun_rlm_est_species_kendall_uniqueness_clr[i,j]),ifelse(shotgun_rlm_p_val_species_kendall_uniqueness_clr[i,j]<=0.05,2*sign(shotgun_rlm_est_species_kendall_uniqueness_clr[i,j]),sign(shotgun_rlm_est_species_kendall_uniqueness_clr[i,j])))
	}
}

print("RLM: species_shannon")
shotgun_rlm_est_species_shannon_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_est_species_shannon_clr) <- select_shotgun_species
colnames(shotgun_rlm_est_species_shannon_clr) <- sorted_study_rows_shotgun
shotgun_rlm_p_val_species_shannon_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_p_val_species_shannon_clr) <- select_shotgun_species
colnames(shotgun_rlm_p_val_species_shannon_clr) <- sorted_study_rows_shotgun

for(i in 1:length(select_shotgun_species))
{
	species_name <- select_shotgun_species[i]
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		study_name <- sorted_study_rows_shotgun[j]
		study_samples <- rownames(shotgun_combined_df_sum_stat_species_clr[shotgun_combined_df_sum_stat_species_clr$study_name == study_name,])
		tryCatch(               
					expr = {   
						print(paste0("species_shannon, ",species_name))
						df_temp <- as.data.frame(shotgun_combined_df_sum_stat_species_clr[study_samples,c(species_name,"species_shannon")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~species_shannon")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						shotgun_rlm_est_species_shannon_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						shotgun_rlm_p_val_species_shannon_clr[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){   
						print(paste0("species_shannon, ",species_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}
shotgun_rlm_est_species_shannon_clr <- apply(shotgun_rlm_est_species_shannon_clr,2,function(x)(ifelse(is.nan(x),0,x)))
shotgun_rlm_q_val_species_shannon_clr <- apply(shotgun_rlm_p_val_species_shannon_clr,2,function(x)(p.adjust(x,method="fdr")))
shotgun_rlm_dir_species_shannon_clr <- as.data.frame(matrix(0,length(select_shotgun_species),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_dir_species_shannon_clr) <- select_shotgun_species
colnames(shotgun_rlm_dir_species_shannon_clr) <- sorted_study_rows_shotgun
for(i in 1:length(select_shotgun_species))
{
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		shotgun_rlm_dir_species_shannon_clr[i,j] <- ifelse(shotgun_rlm_q_val_species_shannon_clr[i,j]<=0.1,3*sign(shotgun_rlm_est_species_shannon_clr[i,j]),ifelse(shotgun_rlm_p_val_species_shannon_clr[i,j]<=0.05,2*sign(shotgun_rlm_est_species_shannon_clr[i,j]),sign(shotgun_rlm_est_species_shannon_clr[i,j])))
	}
}

print("Doing RLM Shotgun Genus")
select_shotgun_genus <- setdiff(colnames(shotgun_combined_df_sum_stat_genus_clr),c("genus_shannon","age","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","study_name","study_condition"))

print("RLM: genus_bray_uniqueness")
shotgun_rlm_est_genus_bray_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_est_genus_bray_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_est_genus_bray_uniqueness_clr) <- sorted_study_rows_shotgun
shotgun_rlm_p_val_genus_bray_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_p_val_genus_bray_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_p_val_genus_bray_uniqueness_clr) <- sorted_study_rows_shotgun

for(i in 1:length(select_shotgun_genus))
{
	genus_name <- select_shotgun_genus[i]
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		study_name <- sorted_study_rows_shotgun[j]
		study_samples <- rownames(shotgun_combined_df_sum_stat_genus_clr[shotgun_combined_df_sum_stat_genus_clr$study_name == study_name,])
		tryCatch(               
					expr = {  
						print(paste0("genus_bray_uniqueness, ",genus_name))
						df_temp <- as.data.frame(shotgun_combined_df_sum_stat_genus_clr[study_samples,c(genus_name,"genus_bray_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(genus_name,"~genus_bray_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						shotgun_rlm_est_genus_bray_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						shotgun_rlm_p_val_genus_bray_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){  
						#print(e)
						print(paste0("genus_bray_uniqueness, ",genus_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}
shotgun_rlm_est_genus_bray_uniqueness_clr <- apply(shotgun_rlm_est_genus_bray_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
shotgun_rlm_q_val_genus_bray_uniqueness_clr <- apply(shotgun_rlm_p_val_genus_bray_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
shotgun_rlm_dir_genus_bray_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_dir_genus_bray_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_dir_genus_bray_uniqueness_clr) <- sorted_study_rows_shotgun
for(i in 1:length(select_shotgun_genus))
{
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		shotgun_rlm_dir_genus_bray_uniqueness_clr[i,j] <- ifelse(shotgun_rlm_q_val_genus_bray_uniqueness_clr[i,j]<=0.1,3*sign(shotgun_rlm_est_genus_bray_uniqueness_clr[i,j]),ifelse(shotgun_rlm_p_val_genus_bray_uniqueness_clr[i,j]<=0.05,2*sign(shotgun_rlm_est_genus_bray_uniqueness_clr[i,j]),sign(shotgun_rlm_est_genus_bray_uniqueness_clr[i,j])))
	}
}

print("RLM: genus_jaccard_uniqueness")
shotgun_rlm_est_genus_jaccard_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_est_genus_jaccard_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_est_genus_jaccard_uniqueness_clr) <- sorted_study_rows_shotgun
shotgun_rlm_p_val_genus_jaccard_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_p_val_genus_jaccard_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_p_val_genus_jaccard_uniqueness_clr) <- sorted_study_rows_shotgun

for(i in 1:length(select_shotgun_genus))
{
	genus_name <- select_shotgun_genus[i]
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		study_name <- sorted_study_rows_shotgun[j]
		study_samples <- rownames(shotgun_combined_df_sum_stat_genus_clr[shotgun_combined_df_sum_stat_genus_clr$study_name == study_name,])
		tryCatch(               
					expr = { 
						print(paste0("genus_jaccard_uniqueness, ",genus_name))
						df_temp <- as.data.frame(shotgun_combined_df_sum_stat_genus_clr[study_samples,c(genus_name,"genus_jaccard_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(genus_name,"~genus_jaccard_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						shotgun_rlm_est_genus_jaccard_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						shotgun_rlm_p_val_genus_jaccard_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){         
						print(paste0("genus_jaccard_uniqueness, ",genus_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}
shotgun_rlm_est_genus_jaccard_uniqueness_clr <- apply(shotgun_rlm_est_genus_jaccard_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
shotgun_rlm_q_val_genus_jaccard_uniqueness_clr <- apply(shotgun_rlm_p_val_genus_jaccard_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
shotgun_rlm_dir_genus_jaccard_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_dir_genus_jaccard_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_dir_genus_jaccard_uniqueness_clr) <- sorted_study_rows_shotgun
for(i in 1:length(select_shotgun_genus))
{
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		shotgun_rlm_dir_genus_jaccard_uniqueness_clr[i,j] <- ifelse(shotgun_rlm_q_val_genus_jaccard_uniqueness_clr[i,j]<=0.1,3*sign(shotgun_rlm_est_genus_jaccard_uniqueness_clr[i,j]),ifelse(shotgun_rlm_p_val_genus_jaccard_uniqueness_clr[i,j]<=0.05,2*sign(shotgun_rlm_est_genus_jaccard_uniqueness_clr[i,j]),sign(shotgun_rlm_est_genus_jaccard_uniqueness_clr[i,j])))
	}
}

print("RLM: genus_aitchison_uniqueness")
shotgun_rlm_est_genus_aitchison_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_est_genus_aitchison_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_est_genus_aitchison_uniqueness_clr) <- sorted_study_rows_shotgun
shotgun_rlm_p_val_genus_aitchison_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_p_val_genus_aitchison_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_p_val_genus_aitchison_uniqueness_clr) <- sorted_study_rows_shotgun

for(i in 1:length(select_shotgun_genus))
{
	genus_name <- select_shotgun_genus[i]
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		study_name <- sorted_study_rows_shotgun[j]
		study_samples <- rownames(shotgun_combined_df_sum_stat_genus_clr[shotgun_combined_df_sum_stat_genus_clr$study_name == study_name,])
		tryCatch(               
					expr = { 
						print(paste0("genus_aitchison_uniqueness, ",genus_name))
						df_temp <- as.data.frame(shotgun_combined_df_sum_stat_genus_clr[study_samples,c(genus_name,"genus_aitchison_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(genus_name,"~genus_aitchison_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						shotgun_rlm_est_genus_aitchison_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						shotgun_rlm_p_val_genus_aitchison_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){  
						print(paste0("genus_aitchison_uniqueness, ",genus_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}
shotgun_rlm_est_genus_aitchison_uniqueness_clr <- apply(shotgun_rlm_est_genus_aitchison_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
shotgun_rlm_q_val_genus_aitchison_uniqueness_clr <- apply(shotgun_rlm_p_val_genus_aitchison_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
shotgun_rlm_dir_genus_aitchison_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_dir_genus_aitchison_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_dir_genus_aitchison_uniqueness_clr) <- sorted_study_rows_shotgun
for(i in 1:length(select_shotgun_genus))
{
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		shotgun_rlm_dir_genus_aitchison_uniqueness_clr[i,j] <- ifelse(shotgun_rlm_q_val_genus_aitchison_uniqueness_clr[i,j]<=0.1,3*sign(shotgun_rlm_est_genus_aitchison_uniqueness_clr[i,j]),ifelse(shotgun_rlm_p_val_genus_aitchison_uniqueness_clr[i,j]<=0.05,2*sign(shotgun_rlm_est_genus_aitchison_uniqueness_clr[i,j]),sign(shotgun_rlm_est_genus_aitchison_uniqueness_clr[i,j])))
	}
}

print("RLM: genus_kendall_uniqueness")
shotgun_rlm_est_genus_kendall_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_est_genus_kendall_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_est_genus_kendall_uniqueness_clr) <- sorted_study_rows_shotgun
shotgun_rlm_p_val_genus_kendall_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_p_val_genus_kendall_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_p_val_genus_kendall_uniqueness_clr) <- sorted_study_rows_shotgun

for(i in 1:length(select_shotgun_genus))
{
	genus_name <- select_shotgun_genus[i]
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		study_name <- sorted_study_rows_shotgun[j]
		study_samples <- rownames(shotgun_combined_df_sum_stat_genus_clr[shotgun_combined_df_sum_stat_genus_clr$study_name == study_name,])
		tryCatch(               
					expr = {    
						print(paste0("genus_kendall_uniqueness, ",genus_name))
						df_temp <- as.data.frame(shotgun_combined_df_sum_stat_genus_clr[study_samples,c(genus_name,"genus_kendall_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(genus_name,"~genus_kendall_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						shotgun_rlm_est_genus_kendall_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						shotgun_rlm_p_val_genus_kendall_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){   
						print(paste0("genus_kendall_uniqueness, ",genus_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}
shotgun_rlm_est_genus_kendall_uniqueness_clr <- apply(shotgun_rlm_est_genus_kendall_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
shotgun_rlm_q_val_genus_kendall_uniqueness_clr <- apply(shotgun_rlm_p_val_genus_kendall_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
shotgun_rlm_dir_genus_kendall_uniqueness_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_dir_genus_kendall_uniqueness_clr) <- select_shotgun_genus
colnames(shotgun_rlm_dir_genus_kendall_uniqueness_clr) <- sorted_study_rows_shotgun
for(i in 1:length(select_shotgun_genus))
{
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		shotgun_rlm_dir_genus_kendall_uniqueness_clr[i,j] <- ifelse(shotgun_rlm_q_val_genus_kendall_uniqueness_clr[i,j]<=0.1,3*sign(shotgun_rlm_est_genus_kendall_uniqueness_clr[i,j]),ifelse(shotgun_rlm_p_val_genus_kendall_uniqueness_clr[i,j]<=0.05,2*sign(shotgun_rlm_est_genus_kendall_uniqueness_clr[i,j]),sign(shotgun_rlm_est_genus_kendall_uniqueness_clr[i,j])))
	}
}

print("RLM: genus_shannon")
shotgun_rlm_est_genus_shannon_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_est_genus_shannon_clr) <- select_shotgun_genus
colnames(shotgun_rlm_est_genus_shannon_clr) <- sorted_study_rows_shotgun
shotgun_rlm_p_val_genus_shannon_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_p_val_genus_shannon_clr) <- select_shotgun_genus
colnames(shotgun_rlm_p_val_genus_shannon_clr) <- sorted_study_rows_shotgun

for(i in 1:length(select_shotgun_genus))
{
	genus_name <- select_shotgun_genus[i]
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		study_name <- sorted_study_rows_shotgun[j]
		study_samples <- rownames(shotgun_combined_df_sum_stat_genus_clr[shotgun_combined_df_sum_stat_genus_clr$study_name == study_name,])
		tryCatch(               
					expr = {   
						print(paste0("genus_shannon, ",genus_name))
						df_temp <- as.data.frame(shotgun_combined_df_sum_stat_genus_clr[study_samples,c(genus_name,"genus_shannon")])
						temp_rlm <- rlm(as.formula(paste0(genus_name,"~genus_shannon")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						shotgun_rlm_est_genus_shannon_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						shotgun_rlm_p_val_genus_shannon_clr[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){   
						print(paste0("genus_shannon, ",genus_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}
shotgun_rlm_est_genus_shannon_clr <- apply(shotgun_rlm_est_genus_shannon_clr,2,function(x)(ifelse(is.nan(x),0,x)))
shotgun_rlm_q_val_genus_shannon_clr <- apply(shotgun_rlm_p_val_genus_shannon_clr,2,function(x)(p.adjust(x,method="fdr")))
shotgun_rlm_dir_genus_shannon_clr <- as.data.frame(matrix(0,length(select_shotgun_genus),length(sorted_study_rows_shotgun)))
rownames(shotgun_rlm_dir_genus_shannon_clr) <- select_shotgun_genus
colnames(shotgun_rlm_dir_genus_shannon_clr) <- sorted_study_rows_shotgun
for(i in 1:length(select_shotgun_genus))
{
	for(j in 1:length(sorted_study_rows_shotgun))
	{
		shotgun_rlm_dir_genus_shannon_clr[i,j] <- ifelse(shotgun_rlm_q_val_genus_shannon_clr[i,j]<=0.1,3*sign(shotgun_rlm_est_genus_shannon_clr[i,j]),ifelse(shotgun_rlm_p_val_genus_shannon_clr[i,j]<=0.05,2*sign(shotgun_rlm_est_genus_shannon_clr[i,j]),sign(shotgun_rlm_est_genus_shannon_clr[i,j])))
	}
}

print("16S Analysis")

print("preparing Species Table")

study_rows_16S <- c("AG","NUAGE","HE","Odamaki","LogMPie")

nuage_species <- grep("Unclassified",colnames(nuage_combined_df_sum_stat_species)[13:208],value=TRUE,invert=TRUE)
logmpie_species <- colnames(logmpie_combined_df_sum_stat_species)[14:382]
ag_species <- colnames(ag_combined_df_sum_stat_species)[14:142]
colnames(he_combined_df_sum_stat_species)[17] <- "Shigella_dysenteriae"
colnames(he_combined_df_sum_stat_species)[18] <- "Escherichia_fergusonii"
colnames(he_combined_df_sum_stat_species)[36] <- "Shigella_flexneri"
colnames(he_combined_df_sum_stat_species)[118] <- "Shigella_sonnei"
colnames(he_combined_df_sum_stat_species)[165] <- "Shigella_boydii"
colnames(he_combined_df_sum_stat_species)[204] <- "Escherichia_coli"
he_species <- colnames(he_combined_df_sum_stat_species)[14:283]
he_species <- he_species[!is.na(he_species)]
odamaki_species <- colnames(odamaki_combined_df_sum_stat_species)[13:255]

common_species_16S <- setdiff(names(which(table(c(nuage_species,ag_species,odamaki_species,he_species,logmpie_species))>=3)),c("L_Markers","G_Markers","G1","G2"))

temp0 <- merge(t(odamaki_combined_df_sum_stat_species[,odamaki_species]),t(he_combined_df_sum_stat_species[,he_species]),by="row.names",all=TRUE)[,-1]
rownames(temp0) <- merge(t(odamaki_combined_df_sum_stat_species[,odamaki_species]),t(he_combined_df_sum_stat_species[,he_species]),by="row.names",all=TRUE)[,1]
temp0 <- apply(temp0,1,function(x)(ifelse(is.na(x),0,x)))

temp1 <- merge(t(temp0),t(ag_combined_df_sum_stat_species[,ag_species]),by="row.names",all=TRUE)[,-1]
rownames(temp1) <- merge(t(temp0),t(ag_combined_df_sum_stat_species[,ag_species]),by="row.names",all=TRUE)[,1]
temp1 <- apply(temp1,1,function(x)(ifelse(is.na(x),0,x)))

temp2 <- merge(t(temp1),t(nuage_combined_df_sum_stat_species[,nuage_species]),by="row.names",all=TRUE)[,-1]
rownames(temp2) <- merge(t(temp1),t(nuage_combined_df_sum_stat_species[,nuage_species]),by="row.names",all=TRUE)[,1]
temp2 <- apply(temp2,1,function(x)(ifelse(is.na(x),0,x)))

temp3 <- merge(t(temp2),t(logmpie_combined_df_sum_stat_species[,logmpie_species]),by="row.names",all=TRUE)[,-1]
rownames(temp3) <- merge(t(temp2),t(logmpie_combined_df_sum_stat_species[,logmpie_species]),by="row.names",all=TRUE)[,1]
temp3 <- apply(temp3,1,function(x)(ifelse(is.na(x),0,x)))

combined_16S_df_select_age_final_species <- temp3[,common_species_16S]
combined_16S_df_select_age_final_species_clr <- as.data.frame(as.matrix(clr(combined_16S_df_select_age_final_species+0.00001)))
combined_16S_df_select_age_final_species_clr <- as.data.frame(t(apply(combined_16S_df_select_age_final_species_clr,1,function(x)(x-min(x)))))

combined_16S_df_diversity_uniqueness <- as.data.frame(rbind(odamaki_combined_df_sum_stat_species[,c("species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","study_name","age")],he_combined_df_sum_stat_species[,c("species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","study_name","age")],ag_combined_df_sum_stat_species[,c("species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","study_name","age")],nuage_combined_df_sum_stat_species[,c("species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","study_name","age")],logmpie_combined_df_sum_stat_species[,c("species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","study_name","age")]))

combined_16S_df_sum_stat_species_clr <- as.data.frame(cbind(combined_16S_df_diversity_uniqueness,combined_16S_df_select_age_final_species_clr[rownames(combined_16S_df_diversity_uniqueness),]))

combined_16S_df_sum_stat_species <- as.data.frame(cbind(combined_16S_df_diversity_uniqueness,combined_16S_df_select_age_final_species[rownames(combined_16S_df_diversity_uniqueness),]))

select_16S_species <- common_species_16S

print("RLM: species_bray_uniqueness 16S")
rlm_16S_est_species_bray_uniqueness_clr <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_est_species_bray_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_est_species_bray_uniqueness_clr) <- study_rows_16S
rlm_16S_p_val_species_bray_uniqueness_clr  <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_p_val_species_bray_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_p_val_species_bray_uniqueness_clr) <- study_rows_16S

for(i in 1:length(select_16S_species))
{
	species_name <- select_16S_species[i]
	for(j in 1:length(study_rows_16S))
	{
		study_name <- study_rows_16S[j]
		study_samples <- rownames(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$study_name == study_name,])
		tryCatch(               
					expr = {   
						print(paste0("species_bray_uniqueness 16S,",species_name))
						df_temp <- as.data.frame(combined_16S_df_sum_stat_species_clr[study_samples,c(species_name,"species_bray_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~species_bray_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						rlm_16S_est_species_bray_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						rlm_16S_p_val_species_bray_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
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
}
rlm_16S_est_species_bray_uniqueness_clr <- apply(rlm_16S_est_species_bray_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
rlm_16S_q_val_species_bray_uniqueness_clr <- apply(rlm_16S_p_val_species_bray_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
rlm_16S_dir_species_bray_uniqueness_clr <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_dir_species_bray_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_dir_species_bray_uniqueness_clr) <- study_rows_16S
for(i in 1:length(select_16S_species))
{
	for(j in 1:length(study_rows_16S))
	{
		rlm_16S_dir_species_bray_uniqueness_clr[i,j] <- ifelse(rlm_16S_q_val_species_bray_uniqueness_clr[i,j]<=0.1,3*sign(rlm_16S_est_species_bray_uniqueness_clr[i,j]),ifelse(rlm_16S_p_val_species_bray_uniqueness_clr[i,j]<=0.05,2*sign(rlm_16S_est_species_bray_uniqueness_clr[i,j]),sign(rlm_16S_est_species_bray_uniqueness_clr[i,j])))
	}
}

print("RLM: species_jaccard_uniqueness 16S")
rlm_16S_est_species_jaccard_uniqueness_clr <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_est_species_jaccard_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_est_species_jaccard_uniqueness_clr) <- study_rows_16S
rlm_16S_p_val_species_jaccard_uniqueness_clr  <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_p_val_species_jaccard_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_p_val_species_jaccard_uniqueness_clr) <- study_rows_16S

for(i in 1:length(select_16S_species))
{
	species_name <- select_16S_species[i]
	for(j in 1:length(study_rows_16S))
	{
		study_name <- study_rows_16S[j]
		study_samples <- rownames(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$study_name == study_name,])
		tryCatch(               
					expr = {   
						print(paste0("species_jaccard_uniqueness 16S,",species_name))
						df_temp <- as.data.frame(combined_16S_df_sum_stat_species_clr[study_samples,c(species_name,"species_jaccard_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~species_jaccard_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						rlm_16S_est_species_jaccard_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						rlm_16S_p_val_species_jaccard_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
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
rlm_16S_est_species_jaccard_uniqueness_clr <- apply(rlm_16S_est_species_jaccard_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
rlm_16S_q_val_species_jaccard_uniqueness_clr <- apply(rlm_16S_p_val_species_jaccard_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
rlm_16S_dir_species_jaccard_uniqueness_clr <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_dir_species_jaccard_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_dir_species_jaccard_uniqueness_clr) <- study_rows_16S
for(i in 1:length(select_16S_species))
{
	for(j in 1:length(study_rows_16S))
	{
		rlm_16S_dir_species_jaccard_uniqueness_clr[i,j] <- ifelse(rlm_16S_q_val_species_jaccard_uniqueness_clr[i,j]<=0.1,3*sign(rlm_16S_est_species_jaccard_uniqueness_clr[i,j]),ifelse(rlm_16S_p_val_species_jaccard_uniqueness_clr[i,j]<=0.05,2*sign(rlm_16S_est_species_jaccard_uniqueness_clr[i,j]),sign(rlm_16S_est_species_jaccard_uniqueness_clr[i,j])))
	}
}

print("RLM: species_aitchison_uniqueness 16S")
rlm_16S_est_species_aitchison_uniqueness_clr <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_est_species_aitchison_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_est_species_aitchison_uniqueness_clr) <- study_rows_16S
rlm_16S_p_val_species_aitchison_uniqueness_clr  <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_p_val_species_aitchison_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_p_val_species_aitchison_uniqueness_clr) <- study_rows_16S

for(i in 1:length(select_16S_species))
{
	species_name <- select_16S_species[i]
	for(j in 1:length(study_rows_16S))
	{
		study_name <- study_rows_16S[j]
		study_samples <- rownames(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$study_name == study_name,])
		tryCatch(               
					expr = {   
						print(paste0("species_aitchison_uniqueness 16S,",species_name))
						df_temp <- as.data.frame(combined_16S_df_sum_stat_species_clr[study_samples,c(species_name,"species_aitchison_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~species_aitchison_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						rlm_16S_est_species_aitchison_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						rlm_16S_p_val_species_aitchison_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
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
rlm_16S_est_species_aitchison_uniqueness_clr <- apply(rlm_16S_est_species_aitchison_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
rlm_16S_q_val_species_aitchison_uniqueness_clr <- apply(rlm_16S_p_val_species_aitchison_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
rlm_16S_dir_species_aitchison_uniqueness_clr <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_dir_species_aitchison_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_dir_species_aitchison_uniqueness_clr) <- study_rows_16S
for(i in 1:length(select_16S_species))
{
	for(j in 1:length(study_rows_16S))
	{
		rlm_16S_dir_species_aitchison_uniqueness_clr[i,j] <- ifelse(rlm_16S_q_val_species_aitchison_uniqueness_clr[i,j]<=0.1,3*sign(rlm_16S_est_species_aitchison_uniqueness_clr[i,j]),ifelse(rlm_16S_p_val_species_aitchison_uniqueness_clr[i,j]<=0.05,2*sign(rlm_16S_est_species_aitchison_uniqueness_clr[i,j]),sign(rlm_16S_est_species_aitchison_uniqueness_clr[i,j])))
	}
}

print("RLM: species_kendall_uniqueness 16S")
rlm_16S_est_species_kendall_uniqueness_clr <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_est_species_kendall_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_est_species_kendall_uniqueness_clr) <- study_rows_16S
rlm_16S_p_val_species_kendall_uniqueness_clr  <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_p_val_species_kendall_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_p_val_species_kendall_uniqueness_clr) <- study_rows_16S

for(i in 1:length(select_16S_species))
{
	species_name <- select_16S_species[i]
	for(j in 1:length(study_rows_16S))
	{
		study_name <- study_rows_16S[j]
		study_samples <- rownames(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$study_name == study_name,])
		tryCatch(               
					expr = {   
						print(paste0("species_kendall_uniqueness 16S,",species_name))
						df_temp <- as.data.frame(combined_16S_df_sum_stat_species_clr[study_samples,c(species_name,"species_kendall_uniqueness")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~species_kendall_uniqueness")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						rlm_16S_est_species_kendall_uniqueness_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						rlm_16S_p_val_species_kendall_uniqueness_clr[i,j] <- f.robftest(temp_rlm)$p.value
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
rlm_16S_est_species_kendall_uniqueness_clr <- apply(rlm_16S_est_species_kendall_uniqueness_clr,2,function(x)(ifelse(is.nan(x),0,x)))
rlm_16S_q_val_species_kendall_uniqueness_clr <- apply(rlm_16S_p_val_species_kendall_uniqueness_clr,2,function(x)(p.adjust(x,method="fdr")))
rlm_16S_dir_species_kendall_uniqueness_clr <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_dir_species_kendall_uniqueness_clr) <- select_16S_species
colnames(rlm_16S_dir_species_kendall_uniqueness_clr) <- study_rows_16S
for(i in 1:length(select_16S_species))
{
	for(j in 1:length(study_rows_16S))
	{
		rlm_16S_dir_species_kendall_uniqueness_clr[i,j] <- ifelse(rlm_16S_q_val_species_kendall_uniqueness_clr[i,j]<=0.1,3*sign(rlm_16S_est_species_kendall_uniqueness_clr[i,j]),ifelse(rlm_16S_p_val_species_kendall_uniqueness_clr[i,j]<=0.05,2*sign(rlm_16S_est_species_kendall_uniqueness_clr[i,j]),sign(rlm_16S_est_species_kendall_uniqueness_clr[i,j])))
	}
}

print("RLM: species_shannon 16S")
rlm_16S_est_species_shannon_clr <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_est_species_shannon_clr) <- select_16S_species
colnames(rlm_16S_est_species_shannon_clr) <- study_rows_16S
rlm_16S_p_val_species_shannon_clr  <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_p_val_species_shannon_clr) <- select_16S_species
colnames(rlm_16S_p_val_species_shannon_clr) <- study_rows_16S

for(i in 1:length(select_16S_species))
{
	species_name <- select_16S_species[i]
	for(j in 1:length(study_rows_16S))
	{
		study_name <- study_rows_16S[j]
		study_samples <- rownames(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$study_name == study_name,])
		tryCatch(               
					expr = {   
						print(paste0("species_shannon 16S,",species_name))
						df_temp <- as.data.frame(combined_16S_df_sum_stat_species_clr[study_samples,c(species_name,"species_shannon")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~species_shannon")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						rlm_16S_est_species_shannon_clr[i,j] <- summary_temp_rlm$coefficients[2,3]
						rlm_16S_p_val_species_shannon_clr[i,j] <- f.robftest(temp_rlm)$p.value
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
rlm_16S_est_species_shannon_clr <- apply(rlm_16S_est_species_shannon_clr,2,function(x)(ifelse(is.nan(x),0,x)))
rlm_16S_q_val_species_shannon_clr <- apply(rlm_16S_p_val_species_shannon_clr,2,function(x)(p.adjust(x,method="fdr")))
rlm_16S_dir_species_shannon_clr <- as.data.frame(matrix(0,length(select_16S_species),length(study_rows_16S)))
rownames(rlm_16S_dir_species_shannon_clr) <- select_16S_species
colnames(rlm_16S_dir_species_shannon_clr) <- study_rows_16S
for(i in 1:length(select_16S_species))
{
	for(j in 1:length(study_rows_16S))
	{
		rlm_16S_dir_species_shannon_clr[i,j] <- ifelse(rlm_16S_q_val_species_shannon_clr[i,j]<=0.1,3*sign(rlm_16S_est_species_shannon_clr[i,j]),ifelse(rlm_16S_p_val_species_shannon_clr[i,j]<=0.05,2*sign(rlm_16S_est_species_shannon_clr[i,j]),sign(rlm_16S_est_species_shannon_clr[i,j])))
	}
}

print("Merging Data sets")

sorted_study_rows <- c(sorted_study_rows_shotgun,study_rows_16S)

combined_df_final_species <- merge(t(shotgun_combined_df_sum_stat_species[,select_shotgun_species]),t(combined_16S_df_sum_stat_species[,select_16S_species]),by="row.names",all=TRUE)[,-1]
rownames(combined_df_final_species) <- merge(t(shotgun_combined_df_sum_stat_species[,select_shotgun_species]),t(combined_16S_df_sum_stat_species[,select_16S_species]),by="row.names",all=TRUE)[,1]
combined_df_final_species <- apply(combined_df_final_species,1,function(x)(ifelse(is.na(x),0,x)))
combined_df_final_species_clr <- as.data.frame(as.matrix(clr(combined_df_final_species + 0.00001)))
combined_df_final_species_clr <- as.data.frame(t(apply(combined_df_final_species,1,function(x)(x-min(x)))))

combined_df_sum_stat <- merge(t(shotgun_combined_df_sum_stat_species[,c("species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","species_shannon","age")]),t(combined_16S_df_sum_stat_species[c("species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","species_shannon","age")]),by="row.names",all=TRUE)[,-1]
rownames(combined_df_sum_stat) <- merge(t(shotgun_combined_df_sum_stat_species[,c("species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","species_shannon","age")]),t(combined_16S_df_sum_stat_species[c("species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","species_shannon","age")]),by="row.names",all=TRUE)[,1]
combined_df_sum_stat <- apply(combined_df_sum_stat,1,function(x)(ifelse(is.na(x),0,x)))

combined_df_sum_stat_species_clr <- as.data.frame(cbind(combined_df_sum_stat,combined_df_final_species_clr[rownames(combined_df_sum_stat),]))

combined_df_sum_stat_species <- as.data.frame(cbind(combined_df_sum_stat,combined_df_final_species[rownames(combined_df_sum_stat),]))

combined_df_sum_stat_species_clr$study_name <- NA
combined_df_sum_stat_species_clr[rownames(shotgun_combined_df_sum_stat_species_clr),"study_name"] <- shotgun_combined_df_sum_stat_species_clr[,"study_name"]
combined_df_sum_stat_species_clr[rownames(combined_16S_df_sum_stat_species),"study_name"] <- combined_16S_df_sum_stat_species[,"study_name"]

combined_df_sum_stat_species$study_name <- NA
combined_df_sum_stat_species[rownames(shotgun_combined_df_sum_stat_species_clr),"study_name"] <- shotgun_combined_df_sum_stat_species_clr[,"study_name"]
combined_df_sum_stat_species[rownames(combined_16S_df_sum_stat_species),"study_name"] <- combined_16S_df_sum_stat_species[,"study_name"]

species_columns <- setdiff(colnames(combined_df_sum_stat_species),c("species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","species_shannon","age","study_name"))

combined_rem_species_bray_uniqueness <- compute_meta_lm_group(combined_df_sum_stat_species,species_columns,"species_bray_uniqueness","study_name",sorted_study_rows)
combined_rem_species_jaccard_uniqueness <- compute_meta_lm_group(combined_df_sum_stat_species,species_columns,"species_jaccard_uniqueness","study_name",sorted_study_rows)
combined_rem_species_aitchison_uniqueness <- compute_meta_lm_group(combined_df_sum_stat_species,species_columns,"species_aitchison_uniqueness","study_name",sorted_study_rows)
combined_rem_species_kendall_uniqueness <- compute_meta_lm_group(combined_df_sum_stat_species,species_columns,"species_kendall_uniqueness","study_name",sorted_study_rows)
combined_rem_species_shannon <- compute_meta_lm_group(combined_df_sum_stat_species,species_columns,"species_shannon","study_name",sorted_study_rows)

combined_rem_species_bray_uniqueness_clr <- compute_meta_lm_group(combined_df_sum_stat_species_clr,species_columns,"species_bray_uniqueness","study_name",sorted_study_rows)
combined_rem_species_jaccard_uniqueness_clr <- compute_meta_lm_group(combined_df_sum_stat_species_clr,species_columns,"species_jaccard_uniqueness","study_name",sorted_study_rows)
combined_rem_species_aitchison_uniqueness_clr <- compute_meta_lm_group(combined_df_sum_stat_species_clr,species_columns,"species_aitchison_uniqueness","study_name",sorted_study_rows)
combined_rem_species_kendall_uniqueness_clr <- compute_meta_lm_group(combined_df_sum_stat_species_clr,species_columns,"species_kendall_uniqueness","study_name",sorted_study_rows)
combined_rem_species_shannon_clr <- compute_meta_lm_group(combined_df_sum_stat_species_clr,species_columns,"species_shannon","study_name",sorted_study_rows)

overall_species_association_pattern_clr <- data.frame(bray_uniqueness=ifelse((abs(combined_rem_species_bray_uniqueness_clr$dir)>=3)&(combined_rem_species_bray_uniqueness_clr$consistency>=0.65),sign(combined_rem_species_bray_uniqueness_clr$dir),0),jaccard_uniqueness=ifelse((abs(combined_rem_species_jaccard_uniqueness_clr$dir)>=3)&(combined_rem_species_jaccard_uniqueness_clr$consistency>=0.65),sign(combined_rem_species_jaccard_uniqueness_clr$dir),0),aitchison_uniqueness=ifelse((abs(combined_rem_species_aitchison_uniqueness_clr$dir)>=3)&(combined_rem_species_aitchison_uniqueness_clr$consistency>=0.65),sign(combined_rem_species_aitchison_uniqueness_clr$dir),0),kendall_uniqueness=ifelse((abs(combined_rem_species_kendall_uniqueness_clr$dir)>=3)&(combined_rem_species_kendall_uniqueness_clr$consistency>=0.65),sign(combined_rem_species_kendall_uniqueness_clr$dir),0),shannon=ifelse((abs(combined_rem_species_shannon_clr$dir)>=3)&(combined_rem_species_shannon_clr$consistency>=0.65),sign(combined_rem_species_shannon_clr$dir),0),row.names=species_columns)

select_detected_species_16S <- names(which(apply(compute_detection(combined_df_sum_stat_species,species_columns,"study_name",study_rows_16S),1,function(x)(length(x[x>0.05])))>=3))
select_detected_species_shotgun <- names(which(apply(compute_detection(combined_df_sum_stat_species,species_columns,"study_name",sorted_study_rows_shotgun),1,function(x)(length(x[x>0.05])))>=14))
common_detected_species <- intersect(select_detected_species_shotgun,select_detected_species_16S)
select_associated_species <- names(which(apply(overall_species_association_pattern_clr,1,function(x)(length(x[x==0]))!=5)))
mat <- as.matrix(overall_species_association_pattern_clr[intersect(intersect(select_shotgun_species,select_16S_species),select_associated_species),])
hmpSpeciesGroupings <- heatmap.2(mat,density.info="none",trace="none",srtCol=90)

#heatmap.2(t(mat),density.info="none",trace="none",srtCol=90,lhei=c(1,5),lwid=c(0.5,5),margins=c(20,5),Rowv=FALSE,col=c("Purple","White","Orange"))



