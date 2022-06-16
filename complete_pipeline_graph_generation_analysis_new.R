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



compute_meta_corr <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(NA,length(grouping_list),3))
	colnames(temp_meta) <- c("dataset","ri","ni")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		dat1 <- data[data[,grouping_variable]==group,var1]
		dat2 <- data[data[,grouping_variable]==group,var2]
		temp_meta[i,2] <- cor.fk(dat1,dat2)
		temp_meta[i,3] <- length(dat1)
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
	print(temp_meta)
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
	return_out <- as.data.frame(matrix(NA,length(feature_list),9))
	rownames(return_out) <- feature_list
	colnames(return_out) <- c("beta","pval","ci.ub","ci.lb","tau2","QE","QEp","qval","dir")
	for(i in 1:length(feature_list))
	{
		species_name <- feature_list[i]
		temp_res <- compute_meta_lm(data,species_name,metadata_var,grouping_var,grouping_list)
		return_out[i,"beta"] <- temp_res$model$beta
		return_out[i,"pval"] <- temp_res$model$pval
		return_out[i,"ci.ub"] <- temp_res$model$ci.ub
		return_out[i,"ci.lb"] <- temp_res$model$ci.lb
		return_out[i,"tau2"] <- temp_res$model$tau2
		return_out[i,"QE"] <- temp_res$model$QE
		return_out[i,"QEp"] <- temp_res$model$QEp
	}
	return_out$qval <- p.adjust(return_out$pval,method="fdr")
	return_out$dir <- ifelse(return_out$qval <= 0.1,3*sign(return_out$beta),ifelse(return_out$pval <= 0.05,2*sign(return_out$beta),sign(return_out$beta)))
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
	colnames(df_est) <- c("est","pval","qval","dir",metadata_feature)
	for(j in 1:length(variable_group))
	{
			species_name <- variable_group[j]
			print(species_name)
			species_val <- data[common_rows,species_name]
			metadata_val <- metadata[common_rows,metadata_feature]
			metadata_grouping <- metadata[common_rows,grouping_feature]
			print(length(species_val[species_val>0]))
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

rem_network1 <- function(data,species_list,group_name,study_list)
{
        species_data <- data[,species_list]
        species_data$group <- data[,group_name]
        est_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
		rownames(est_matrix) <- species_list
		colnames(est_matrix) <- species_list
		pval_matrix <- as.data.frame(matrix(1,length(species_list),length(species_list)))
		rownames(pval_matrix) <- species_list
		colnames(pval_matrix) <- species_list
        for(i in 1:length(species_list))
        {
                species1 <- species_list[i]
                #print(species1)
                for(j in 1:length(species_list))
                {
                        species2 <- species_list[j]
                        print(paste0(species1,",",species2))
                        if(species1 != species2)
                        {
								tryCatch(               
								expr = {                     
										temp_rem <- compute_meta_corr(data,species1,species2,group_name,study_list)
										#print(temp_rem)
										est_matrix[species1,species2] <- as.numeric(temp_rem$beta)
										pval_matrix[species1,species2] <- temp_rem$pval
		
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
        }
        qval_matrix <- apply(pval_matrix,2,function(x)(p.adjust(x,method="fdr")))
		dir_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
		rownames(dir_matrix) <- species_list
		colnames(dir_matrix) <- species_list
		for(i in 1:length(species_list))
        {
          for(j in 1:length(species_list))
          {
			dir_matrix[i,j] <- ifelse(qval_matrix[i,j]<=0.001,sign(est_matrix[i,j]),0)
		  }
		}
		return_list <- list("est"=est_matrix,"pval"=pval_matrix,"qval"=qval_matrix,"dir"=dir_matrix)
		return(return_list)
}

rem_network2 <- function(data,species_list,feature_list,group_name,study_list)
{
        species_data <- data[,species_list]
        species_data$group <- data[,group_name]
        est_matrix <- as.data.frame(matrix(0,length(species_list),length(feature_list)))
		rownames(est_matrix) <- species_list
		colnames(est_matrix) <- feature_list
		pval_matrix <- as.data.frame(matrix(1,length(species_list),length(feature_list)))
		rownames(pval_matrix) <- species_list
		colnames(pval_matrix) <- feature_list
        for(i in 1:length(species_list))
        {
                species1 <- species_list[i]
                #print(species1)
                for(j in 1:length(feature_list))
                {
                        species2 <- feature_list[j]
                        print(paste0(species1,",",species2))
                        if(species1 != species2)
                        {
								tryCatch(               
								expr = {                     
										temp_rem <- compute_meta_lm(data,species1,species2,group_name,study_list)
										#print(temp_rem)
										est_matrix[species1,species2] <- as.numeric(temp_rem$model$beta)
										pval_matrix[species1,species2] <- temp_rem$model$pval
		
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
        }
        qval_matrix <- apply(pval_matrix,2,function(x)(p.adjust(x,method="fdr")))
		dir_matrix <- as.data.frame(matrix(0,length(species_list),length(feature_list)))
		rownames(dir_matrix) <- species_list
		colnames(dir_matrix) <- feature_list
		for(i in 1:length(species_list))
        {
          for(j in 1:length(feature_list))
          {
			dir_matrix[i,j] <- ifelse(qval_matrix[i,j]<=0.15,2*sign(est_matrix[i,j]),ifelse(pval_matrix[i,j]<=0.05,sign(est_matrix[i,j]),0))
		  }
		}
		return_list <- list("est"=est_matrix,"pval"=pval_matrix,"qval"=qval_matrix,"dir"=dir_matrix)
		return(return_list)
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



load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\he_combined_df_sum_species.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ag_combined_df_sum_species.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\odamaki_combined_df_sum_species.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\logmpie_combined_df_sum_species.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\nuage_combined_df_sum_species.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\taxa_groupings_new.RData")

print("Preprocessing species profiles")
## AG All Clear
ag_df_species <- ag_combined_df_sum_stat_species[,c(13:141)]
#ag_df_species <- clr(ag_df_species)
ag_df_species <- as.matrix(clr(ag_df_species+0.00001))
ag_df_species <- as.data.frame(t(apply(ag_df_species,1,function(x)(x-min(x)))))
## He
colnames(he_combined_df_sum_stat_species)[grep("Escherichia.Shigella_dysenteriae",colnames(he_combined_df_sum_stat_species))] <- "Shigella_dysenteriae"
colnames(he_combined_df_sum_stat_species)[grep("Escherichia.Shigella_fergusonii",colnames(he_combined_df_sum_stat_species))] <- "Escherichia_fergusonii"
colnames(he_combined_df_sum_stat_species)[grep("Escherichia.Shigella_flexneri",colnames(he_combined_df_sum_stat_species))] <- "Shigella_flexneri"
colnames(he_combined_df_sum_stat_species)[grep("Escherichia.Shigella_sonnei",colnames(he_combined_df_sum_stat_species))] <- "Shigella_sonnei"
colnames(he_combined_df_sum_stat_species)[grep("Escherichia.Shigella_boydii",colnames(he_combined_df_sum_stat_species))] <- "Shigella_boydii"
colnames(he_combined_df_sum_stat_species)[grep("Escherichia.Shigella_coli",colnames(he_combined_df_sum_stat_species))] <- "Escherichia_coli"
he_df_species <- he_combined_df_sum_stat_species[,c(15:283)]
#he_df_species <- clr(he_df_species)
he_df_species <- as.matrix(clr(he_df_species+0.00001))
he_df_species <- as.data.frame(t(apply(he_df_species,1,function(x)(x-min(x)))))
## Odamaki
colnames(odamaki_combined_df_sum_stat_species)[grep("Escherichia.Shigella_coli",colnames(odamaki_combined_df_sum_stat_species))] <- "Escherichia_coli"
colnames(odamaki_combined_df_sum_stat_species)[grep("Escherichia.Shigella_sonnei",colnames(odamaki_combined_df_sum_stat_species))] <- "Shigella_sonnei"
colnames(odamaki_combined_df_sum_stat_species)[grep("Escherichia.Shigella_flexneri",colnames(odamaki_combined_df_sum_stat_species))] <- "Shigella_flexneri"
odamaki_df_species <- odamaki_combined_df_sum_stat_species[,c(13:255)]
#odamaki_df_species <- clr(odamaki_df_species)
odamaki_df_species <- as.matrix(clr(odamaki_df_species+0.00001))
odamaki_df_species <- as.data.frame(t(apply(odamaki_df_species,1,function(x)(x-min(x)))))
## NUAGE
nuage_combined_df_sum_stat_species <- nuage_combined_df_sum_stat_species[,grep("classified",colnames(nuage_combined_df_sum_stat_species),value=TRUE,invert=TRUE)]
nuage_df_species <- nuage_combined_df_sum_stat_species[,c(13:208)]
#nuage_df_species <- clr(nuage_df_species)
nuage_df_species <- as.matrix(clr(nuage_df_species+0.00001))
nuage_df_species <- as.data.frame(t(apply(nuage_df_species,1,function(x)(x-min(x)))))
## LogMPie
colnames(logmpie_combined_df_sum_stat_species)[grep("Escherichia.Shigella_boydii",colnames(logmpie_combined_df_sum_stat_species))] <- "Shigella_boydii"
colnames(logmpie_combined_df_sum_stat_species)[grep("Escherichia.Shigella_coli",colnames(logmpie_combined_df_sum_stat_species))] <- "Escherichia_coli"
colnames(logmpie_combined_df_sum_stat_species)[grep("Escherichia.Shigella_dysenteriae",colnames(logmpie_combined_df_sum_stat_species))] <- "Shigella_dysenteriae"
colnames(logmpie_combined_df_sum_stat_species)[grep("Escherichia.Shigella_sonnei",colnames(logmpie_combined_df_sum_stat_species))] <- "Shigella_sonnei"
colnames(logmpie_combined_df_sum_stat_species)[grep("Escherichia.Shigella_flexneri",colnames(logmpie_combined_df_sum_stat_species))] <- "Shigella_flexneri"
logmpie_df_species <- logmpie_combined_df_sum_stat_species[,c(14:286)]
#logmpie_df_species <- clr(logmpie_df_species)
logmpie_df_species <- as.matrix(clr(logmpie_df_species+0.00001))
logmpie_df_species <- as.data.frame(t(apply(logmpie_df_species,1,function(x)(x-min(x)))))

print("Merging Data Sets")
temp0_16S <- merge(t(ag_df_species),t(he_df_species),by="row.names",all=TRUE)[,-1]
rownames(temp0_16S) <- merge(t(ag_df_species),t(he_df_species),by="row.names",all=TRUE)[,1]
temp0_16S <- apply(temp0_16S,1,function(x)(ifelse(is.na(x),0,x)))

temp1_16S <- merge(t(temp0_16S),t(odamaki_df_species),by="row.names",all=TRUE)[,-1]
rownames(temp1_16S) <- merge(t(temp0_16S),t(odamaki_df_species),by="row.names",all=TRUE)[,1]
temp1_16S <- apply(temp1_16S,1,function(x)(ifelse(is.na(x),0,x)))

temp2_16S <- merge(t(temp1_16S),t(logmpie_df_species),by="row.names",all=TRUE)[,-1]
rownames(temp2_16S) <- merge(t(temp1_16S),t(logmpie_df_species),by="row.names",all=TRUE)[,1]
temp2_16S <- apply(temp2_16S,1,function(x)(ifelse(is.na(x),0,x)))

temp3_16S <- merge(t(temp2_16S),t(nuage_df_species),by="row.names",all=TRUE)[,-1]
rownames(temp3_16S) <- merge(t(temp2_16S),t(nuage_df_species),by="row.names",all=TRUE)[,1]
combined_16S_df_sum_stat_species <- as.data.frame(apply(temp3_16S,1,function(x)(ifelse(is.na(x),0,x))))

combined_16S_df_sum_stat_species$study_name <- NA
combined_16S_df_sum_stat_species[rownames(ag_combined_df_sum_stat_species),"study_name"] <- "AG"
combined_16S_df_sum_stat_species[rownames(he_combined_df_sum_stat_species),"study_name"] <- "He"
combined_16S_df_sum_stat_species[rownames(odamaki_combined_df_sum_stat_species),"study_name"] <- "Odamaki"
combined_16S_df_sum_stat_species[rownames(logmpie_combined_df_sum_stat_species),"study_name"] <- "LogMPie"
combined_16S_df_sum_stat_species[rownames(nuage_combined_df_sum_stat_species),"study_name"] <- "NUAGE"

combined_16S_df_sum_stat_species$age <- NA
combined_16S_df_sum_stat_species[rownames(ag_combined_df_sum_stat_species),"age"] <- ag_combined_df_sum_stat_species[,"age"]
combined_16S_df_sum_stat_species[rownames(he_combined_df_sum_stat_species),"age"] <-he_combined_df_sum_stat_species[,"age"]
combined_16S_df_sum_stat_species[rownames(odamaki_combined_df_sum_stat_species),"age"] <-odamaki_combined_df_sum_stat_species[,"age"]
combined_16S_df_sum_stat_species[rownames(logmpie_combined_df_sum_stat_species),"age"] <-logmpie_combined_df_sum_stat_species[,"age"]
combined_16S_df_sum_stat_species[rownames(nuage_combined_df_sum_stat_species),"age"] <-nuage_combined_df_sum_stat_species[,"age"]

combined_16S_df_sum_stat_species$species_shannon <- NA
combined_16S_df_sum_stat_species[rownames(ag_combined_df_sum_stat_species),"species_shannon"] <- ag_combined_df_sum_stat_species[,"species_shannon"]
combined_16S_df_sum_stat_species[rownames(he_combined_df_sum_stat_species),"species_shannon"] <-he_combined_df_sum_stat_species[,"species_shannon"]
combined_16S_df_sum_stat_species[rownames(odamaki_combined_df_sum_stat_species),"species_shannon"] <-odamaki_combined_df_sum_stat_species[,"species_shannon"]
combined_16S_df_sum_stat_species[rownames(logmpie_combined_df_sum_stat_species),"species_shannon"] <-logmpie_combined_df_sum_stat_species[,"species_shannon"]
combined_16S_df_sum_stat_species[rownames(nuage_combined_df_sum_stat_species),"species_shannon"] <-nuage_combined_df_sum_stat_species[,"species_shannon"]

combined_16S_df_sum_stat_species$species_bray_uniqueness <- NA
combined_16S_df_sum_stat_species[rownames(ag_combined_df_sum_stat_species),"species_bray_uniqueness"] <- ag_combined_df_sum_stat_species[,"species_bray_uniqueness"]
combined_16S_df_sum_stat_species[rownames(he_combined_df_sum_stat_species),"species_bray_uniqueness"] <-he_combined_df_sum_stat_species[,"species_bray_uniqueness"]
combined_16S_df_sum_stat_species[rownames(odamaki_combined_df_sum_stat_species),"species_bray_uniqueness"] <-odamaki_combined_df_sum_stat_species[,"species_bray_uniqueness"]
combined_16S_df_sum_stat_species[rownames(logmpie_combined_df_sum_stat_species),"species_bray_uniqueness"] <-logmpie_combined_df_sum_stat_species[,"species_bray_uniqueness"]
combined_16S_df_sum_stat_species[rownames(nuage_combined_df_sum_stat_species),"species_bray_uniqueness"] <-nuage_combined_df_sum_stat_species[,"species_bray_uniqueness"]

combined_16S_df_sum_stat_species$species_jaccard_uniqueness <- NA
combined_16S_df_sum_stat_species[rownames(ag_combined_df_sum_stat_species),"species_jaccard_uniqueness"] <- ag_combined_df_sum_stat_species[,"species_jaccard_uniqueness"]
combined_16S_df_sum_stat_species[rownames(he_combined_df_sum_stat_species),"species_jaccard_uniqueness"] <-he_combined_df_sum_stat_species[,"species_jaccard_uniqueness"]
combined_16S_df_sum_stat_species[rownames(odamaki_combined_df_sum_stat_species),"species_jaccard_uniqueness"] <-odamaki_combined_df_sum_stat_species[,"species_jaccard_uniqueness"]
combined_16S_df_sum_stat_species[rownames(logmpie_combined_df_sum_stat_species),"species_jaccard_uniqueness"] <-logmpie_combined_df_sum_stat_species[,"species_jaccard_uniqueness"]
combined_16S_df_sum_stat_species[rownames(nuage_combined_df_sum_stat_species),"species_jaccard_uniqueness"] <-nuage_combined_df_sum_stat_species[,"species_jaccard_uniqueness"]

combined_16S_df_sum_stat_species$species_manhattan_uniqueness <- NA
combined_16S_df_sum_stat_species[rownames(ag_combined_df_sum_stat_species),"species_manhattan_uniqueness"] <- ag_combined_df_sum_stat_species[,"species_manhattan_uniqueness"]
combined_16S_df_sum_stat_species[rownames(he_combined_df_sum_stat_species),"species_manhattan_uniqueness"] <-he_combined_df_sum_stat_species[,"species_manhattan_uniqueness"]
combined_16S_df_sum_stat_species[rownames(odamaki_combined_df_sum_stat_species),"species_manhattan_uniqueness"] <-odamaki_combined_df_sum_stat_species[,"species_manhattan_uniqueness"]
combined_16S_df_sum_stat_species[rownames(logmpie_combined_df_sum_stat_species),"species_manhattan_uniqueness"] <-logmpie_combined_df_sum_stat_species[,"species_manhattan_uniqueness"]
combined_16S_df_sum_stat_species[rownames(nuage_combined_df_sum_stat_species),"species_manhattan_uniqueness"] <-nuage_combined_df_sum_stat_species[,"species_manhattan_uniqueness"]

combined_16S_df_sum_stat_species$species_kendall_uniqueness <- NA
combined_16S_df_sum_stat_species[rownames(ag_combined_df_sum_stat_species),"species_kendall_uniqueness"] <- ag_combined_df_sum_stat_species[,"species_kendall_uniqueness"]
combined_16S_df_sum_stat_species[rownames(he_combined_df_sum_stat_species),"species_kendall_uniqueness"] <-he_combined_df_sum_stat_species[,"species_kendall_uniqueness"]
combined_16S_df_sum_stat_species[rownames(odamaki_combined_df_sum_stat_species),"species_kendall_uniqueness"] <-odamaki_combined_df_sum_stat_species[,"species_kendall_uniqueness"]
combined_16S_df_sum_stat_species[rownames(logmpie_combined_df_sum_stat_species),"species_kendall_uniqueness"] <-logmpie_combined_df_sum_stat_species[,"species_kendall_uniqueness"]
combined_16S_df_sum_stat_species[rownames(nuage_combined_df_sum_stat_species),"species_kendall_uniqueness"] <-nuage_combined_df_sum_stat_species[,"species_kendall_uniqueness"]

combined_16S_df_sum_stat_species$species_shannon <- NA
combined_16S_df_sum_stat_species[rownames(ag_combined_df_sum_stat_species),"species_shannon"] <- ag_combined_df_sum_stat_species[,"species_shannon"]
combined_16S_df_sum_stat_species[rownames(he_combined_df_sum_stat_species),"species_shannon"] <-he_combined_df_sum_stat_species[,"species_shannon"]
combined_16S_df_sum_stat_species[rownames(odamaki_combined_df_sum_stat_species),"species_shannon"] <-odamaki_combined_df_sum_stat_species[,"species_shannon"]
combined_16S_df_sum_stat_species[rownames(logmpie_combined_df_sum_stat_species),"species_shannon"] <-logmpie_combined_df_sum_stat_species[,"species_shannon"]
combined_16S_df_sum_stat_species[rownames(nuage_combined_df_sum_stat_species),"species_shannon"] <-nuage_combined_df_sum_stat_species[,"species_shannon"]

select_species_16S <- rownames(species_group_16S)

print("Association with Microbiome Properties")
rem_16S_sum_stat_species_young <- rem_network2(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$age<60,],select_species_16S,c("species_bray_uniqueness","species_jaccard_uniqueness","species_manhattan_uniqueness","species_kendall_uniqueness","species_shannon"),"study_name",unique(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$age<60,"study_name"]))

rem_16S_sum_stat_species_elderly <- rem_network2(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$age>=60,],select_species_16S,c("species_bray_uniqueness","species_jaccard_uniqueness","species_manhattan_uniqueness","species_kendall_uniqueness","species_shannon"),"study_name",unique(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$age>=60,"study_name"]))

filtered_select_species_16S <- names(which(apply(rem_16S_sum_stat_species$dir,1,function(x)(length(x[x!=0])))>0))

print("Building Networks")

rem_network_16S_elderly <- rem_network1(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$age>=60,],filtered_select_species_16S,"study_name",unique(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$age>=60,"study_name"]))

co_occurrence_graph_16S_elderly <- graph_from_adjacency_matrix(apply(rem_network_16S_elderly$dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper")

rem_network_16S_young <- rem_network1(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$age<60,],filtered_select_species_16S,"study_name",unique(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$age<60,"study_name"]))

co_occurrence_graph_16S_young <- graph_from_adjacency_matrix(apply(rem_network_16S_young$dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper")

study_names_16S <- unique(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$age>=60,"study_name"])

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_combined_df_sum_species.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_combined_df_sum_species.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\taxa_groupings.RData")

print("Preprocessing species profiles")
## CMD3 All Clear
cmd3_df_species <- cmd3_combined_df_sum_stat_species[,c(41:210)]
cmd3_df_species <- clr(cmd3_df_species)
## ISC All Clear
isc_df_species <- isc_combined_df_sum_stat_species[,c(14:124)]
isc_df_species <- clr(isc_df_species)

print("Merging Data Sets")
temp0_shotgun <- merge(t(cmd3_df_species),t(isc_df_species),by="row.names",all=TRUE)[,-1]
rownames(temp0_shotgun) <- merge(t(cmd3_df_species),t(isc_df_species),by="row.names",all=TRUE)[,1]
combined_shotgun_df_sum_stat_species <- as.data.frame(apply(temp0_shotgun,1,function(x)(ifelse(is.na(x),0,x))))
rm(temp0_shotgun)

combined_shotgun_df_sum_stat_species$study_name <- NA
combined_shotgun_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"study_name"] <- cmd3_combined_df_sum_stat_species[,"study_name"]
combined_shotgun_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"study_name"] <- "ISC"

combined_shotgun_df_sum_stat_species$age <- NA
combined_shotgun_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"age"] <- cmd3_combined_df_sum_stat_species[,"age"]
combined_shotgun_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"age"] <- isc_combined_df_sum_stat_species[,"age"]

combined_shotgun_df_sum_stat_species$species_bray_uniqueness <- NA
combined_shotgun_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"species_bray_uniqueness"] <- cmd3_combined_df_sum_stat_species[,"species_bray_uniqueness"]
combined_shotgun_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"species_bray_uniqueness"] <- isc_combined_df_sum_stat_species[,"species_bray_uniqueness"]

combined_shotgun_df_sum_stat_species$species_jaccard_uniqueness <- NA
combined_shotgun_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"species_jaccard_uniqueness"] <- cmd3_combined_df_sum_stat_species[,"species_jaccard_uniqueness"]
combined_shotgun_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"species_jaccard_uniqueness"] <- isc_combined_df_sum_stat_species[,"species_jaccard_uniqueness"]

combined_shotgun_df_sum_stat_species$species_manhattan_uniqueness <- NA
combined_shotgun_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"species_manhattan_uniqueness"] <- cmd3_combined_df_sum_stat_species[,"species_manhattan_uniqueness"]
combined_shotgun_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"species_manhattan_uniqueness"] <- isc_combined_df_sum_stat_species[,"species_manhattan_uniqueness"]

combined_shotgun_df_sum_stat_species$species_kendall_uniqueness <- NA
combined_shotgun_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"species_kendall_uniqueness"] <- cmd3_combined_df_sum_stat_species[,"species_kendall_uniqueness"]
combined_shotgun_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"species_kendall_uniqueness"] <- isc_combined_df_sum_stat_species[,"species_kendall_uniqueness"]

combined_shotgun_df_sum_stat_species$species_shannon <- NA
combined_shotgun_df_sum_stat_species[rownames(cmd3_combined_df_sum_stat_species),"species_shannon"] <- cmd3_combined_df_sum_stat_species[,"species_shannon"]
combined_shotgun_df_sum_stat_species[rownames(isc_combined_df_sum_stat_species),"species_shannon"] <- isc_combined_df_sum_stat_species[,"species_shannon"]

select_species_shotgun <- rownames(species_group_shotgun)

print("Association with Microbiome Properties")
rem_shotgun_sum_stat_species_young <- rem_network2(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$age<60,],select_species_shotgun,c("species_bray_uniqueness","species_jaccard_uniqueness","species_manhattan_uniqueness","species_kendall_uniqueness","species_shannon"),"study_name",names(which(table(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$age<60,"study_name"])>10)))

rem_shotgun_sum_stat_species_elderly <- rem_network2(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$age>=60,],select_species_shotgun,c("species_bray_uniqueness","species_jaccard_uniqueness","species_manhattan_uniqueness","species_kendall_uniqueness","species_shannon"),"study_name",names(which(table(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$age>=60,"study_name"])>10)))

filtered_select_species_shotgun <- names(which(apply(rem_shotgun_sum_stat_species$dir,1,function(x)(length(x[x!=0])))>0))

print("Building Networks")

rem_network_shotgun_elderly <- rem_network1(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$age>=60,],filtered_select_species_shotgun,"study_name",names(which(table(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$age>=60,"study_name"])>10)))

co_occurrence_graph_shotgun_elderly <- graph_from_adjacency_matrix(apply(rem_network_shotgun_elderly$dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper")

rem_network_shotgun_young <- rem_network1(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$age<60,],filtered_select_species_shotgun,"study_name",names(which(table(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$age<60,"study_name"])>10)))

co_occurrence_graph_shotgun_young <- graph_from_adjacency_matrix(apply(rem_network_shotgun_young$dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper")

study_names_shotgun <- unique(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$age>=60,"study_name"])


sorted_study_rows <- c("HMP_2019_ibdmdb","SankaranarayananK_2015","CosteaPI_2017","AsnicarF_2021","HansenLBS_2018","NielsenHB_2014","SchirmerM_2016","WirbelJ_2018","ZellerG_2014","KeohaneDM_2020","QinN_2014","YeZ_2018","YachidaS_2019","DhakanDB_2019","GuptaA_2019","BritoIL_2016","PehrssonE_2016","LokmerA_2019","PasolliE_2019","RubelMA_2020","RampelliS_2015","AG","NUAGE","ISC","HE","Odamaki","LogMPie")

sorted_study_rows_16S <- sorted_study_rows[c(22,23,25:27)]

sorted_study_rows_shotgun <- sorted_study_rows[c(1:21,24)]

study_profiles <- as.data.frame(matrix(0,length(sorted_study_rows),4))
rownames(study_profiles) <- sorted_study_rows
colnames(study_profiles) <- c("Min_Age","Min_Age","Number_Young","Number_Elderly")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	print(study_name)
	if(study_name %in% sorted_study_rows_shotgun)
	{
		study_profiles[i,1] <- min(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$study_name == study_name,"age"])
		study_profiles[i,2] <- max(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$study_name == study_name,"age"])
		study_profiles[i,3] <- nrow(combined_shotgun_df_sum_stat_species[(combined_shotgun_df_sum_stat_species$study_name == study_name)&(combined_shotgun_df_sum_stat_species$age < 60),])
		study_profiles[i,4] <- nrow(combined_shotgun_df_sum_stat_species[(combined_shotgun_df_sum_stat_species$study_name == study_name)&(combined_shotgun_df_sum_stat_species$age >= 60),])
	}
	if(study_name %in% sorted_study_rows_16S)
	{
		study_profiles[i,1] <- min(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$study_name == study_name,"age"])
		study_profiles[i,2] <- max(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$study_name == study_name,"age"])
		study_profiles[i,3] <- nrow(combined_16S_df_sum_stat_species[(combined_16S_df_sum_stat_species$study_name == study_name)&(combined_16S_df_sum_stat_species$age < 60),])
		study_profiles[i,4] <- nrow(combined_16S_df_sum_stat_species[(combined_16S_df_sum_stat_species$study_name == study_name)&(combined_16S_df_sum_stat_species$age >= 60),])
	}
	
}

selected_network_studies <- rownames(study_profiles[(study_profiles[,3]>=50)&(study_profiles[,4]>=50),])

for(i in 1:length(selected_network_studies))
{
	study_name <- selected_network_studies[i]
	if(study_name %in% sorted_study_rows_shotgun)
	{
		temp_kendall <- cor.fk(combined_shotgun_df_sum_stat_species[(combined_shotgun_df_sum_stat_species$study_name==study_name)&(combined_shotgun_df_sum_stat_species$age>=60),select_species_shotgun])
		n_samples <- nrow(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$study_name==study_name,select_species_shotgun])
		temp_kendall <- apply(temp_kendall,2,function(x)(ifelse(is.nan(x),0,x)))
		assign(paste0("kendall_est_elderly_",study_name),temp_kendall)
		qval <- corr.p(temp_kendall,n=n_samples,adjust="fdr")
		temp_qval <- qval$p
		temp_qval <- apply(temp_qval,2,function(x)(ifelse(is.nan(x),0,x)))
		assign(paste0("kendall_qval_elderly_",study_name),temp_qval)
		temp_dir <- matrix(0,nrow(temp_kendall),ncol(temp_kendall))
		rownames(temp_dir) <- rownames(temp_kendall)
		colnames(temp_dir) <- colnames(temp_kendall)
		for(i in 1:nrow(temp_dir))
		{
			for(j in 1:ncol(temp_dir))
			{
				temp_dir[i,j] <- ifelse(temp_qval[i,j]<=0.01,sign(temp_qval[i,j]),0)
			}
		}
		assign(paste0("kendall_dir_elderly_",study_name),temp_dir)
		temp_graph1 <- graph_from_adjacency_matrix(apply(temp_dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper")
		assign(paste0("co_occurrence_graph_elderly_",study_name),graph_from_adjacency_matrix(apply(temp_dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper"))
		temp_properties <- data.frame(species_name = names(V(temp_graph1)),kendall_uniqueness = rem_shotgun_sum_stat_species_elderly$est[names(V(temp_graph1)),1],bray_uniqueness = rem_shotgun_sum_stat_species_elderly$est[names(V(temp_graph1)),2],jaccard_uniqueness = rem_shotgun_sum_stat_species_elderly$est[names(V(temp_graph1)),3],manhattan_uniqueness = rem_shotgun_sum_stat_species_elderly$est[names(V(temp_graph1)),4],shannon = rem_shotgun_sum_stat_species_elderly$est[names(V(temp_graph1)),5],degree=range_scale(degree(temp_graph1)),betweenness=range_scale(betweenness(temp_graph1)),hub_score = range_scale(hub_score(temp_graph1)$vector),group=species_group_shotgun[names(V(temp_graph1)),7])
		temp_properties$study_name <- study_name
		assign(paste0("df_species_properties_elderly_",study_name),temp_properties)
		
		temp_kendall <- cor.fk(combined_shotgun_df_sum_stat_species[(combined_shotgun_df_sum_stat_species$study_name==study_name)&(combined_shotgun_df_sum_stat_species$age<60),select_species_shotgun])
		n_samples <- nrow(combined_shotgun_df_sum_stat_species[combined_shotgun_df_sum_stat_species$study_name==study_name,select_species_shotgun])
		temp_kendall <- apply(temp_kendall,2,function(x)(ifelse(is.nan(x),0,x)))
		assign(paste0("kendall_est_young_",study_name),temp_kendall)
		qval <- corr.p(temp_kendall,n=n_samples,adjust="fdr")
		temp_qval <- qval$p
		temp_qval <- apply(temp_qval,2,function(x)(ifelse(is.nan(x),0,x)))
		assign(paste0("kendall_qval_young_",study_name),temp_qval)
		temp_dir <- matrix(0,nrow(temp_kendall),ncol(temp_kendall))
		rownames(temp_dir) <- rownames(temp_kendall)
		colnames(temp_dir) <- colnames(temp_kendall)
		for(i in 1:nrow(temp_dir))
		{
			for(j in 1:ncol(temp_dir))
			{
				temp_dir[i,j] <- ifelse(temp_qval[i,j]<=0.01,sign(temp_qval[i,j]),0)
			}
		}
		assign(paste0("kendall_dir_young_",study_name),temp_dir)
		temp_graph <- graph_from_adjacency_matrix(apply(temp_dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper")
		assign(paste0("co_occurrence_graph_young_",study_name),graph_from_adjacency_matrix(apply(temp_dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper"))
		temp_properties <- data.frame(species_name = names(V(temp_graph)),kendall_uniqueness = rem_shotgun_sum_stat_species_young$est[names(V(temp_graph)),1],bray_uniqueness = rem_shotgun_sum_stat_species_young$est[names(V(temp_graph)),2],jaccard_uniqueness = rem_shotgun_sum_stat_species_young$est[names(V(temp_graph)),3],manhattan_uniqueness = rem_shotgun_sum_stat_species_young$est[names(V(temp_graph)),4],shannon = rem_shotgun_sum_stat_species_young$est[names(V(temp_graph)),5],degree=range_scale(degree(temp_graph)),betweenness=range_scale(betweenness(temp_graph)),hub_score = range_scale(hub_score(temp_graph)$vector),group=species_group_shotgun[names(V(temp_graph)),7])
		temp_properties$study_name <- study_name
		assign(paste0("df_species_properties_young_",study_name),temp_properties)
		
	}
	if(study_name %in% sorted_study_rows_16S)
	{
		temp_kendall <- cor.fk(combined_16S_df_sum_stat_species[(combined_16S_df_sum_stat_species$study_name==study_name)&(combined_16S_df_sum_stat_species$age>=60),select_species_16S])
		n_samples <- nrow(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$study_name==study_name,select_species_16S])
		temp_kendall <- apply(temp_kendall,2,function(x)(ifelse(is.nan(x),0,x)))
		assign(paste0("kendall_est_elderly_",study_name),temp_kendall)
		qval <- corr.p(temp_kendall,n=n_samples,adjust="fdr")
		temp_qval <- qval$p
		temp_qval <- apply(temp_qval,2,function(x)(ifelse(is.nan(x),0,x)))
		assign(paste0("kendall_qval_elderly_",study_name),temp_qval)
		temp_dir <- matrix(0,nrow(temp_kendall),ncol(temp_kendall))
		rownames(temp_dir) <- rownames(temp_kendall)
		colnames(temp_dir) <- colnames(temp_kendall)
		for(i in 1:nrow(temp_dir))
		{
			for(j in 1:ncol(temp_dir))
			{
				temp_dir[i,j] <- ifelse(temp_qval[i,j]<=0.01,sign(temp_qval[i,j]),0)
			}
		}
		assign(paste0("kendall_dir_elderly_",study_name),temp_dir)
		temp_graph <- graph_from_adjacency_matrix(apply(temp_dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper")
		assign(paste0("co_occurrence_graph_elderly_",study_name),graph_from_adjacency_matrix(apply(temp_dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper"))
		temp_properties <- data.frame(species_name = names(V(temp_graph)),kendall_uniqueness = rem_16S_sum_stat_species_elderly$est[names(V(temp_graph)),1],bray_uniqueness = rem_16S_sum_stat_species_elderly$est[names(V(temp_graph)),2],jaccard_uniqueness = rem_16S_sum_stat_species_elderly$est[names(V(temp_graph)),3],manhattan_uniqueness = rem_16S_sum_stat_species_elderly$est[names(V(temp_graph)),4],shannon = rem_16S_sum_stat_species_elderly$est[names(V(temp_graph)),5],degree=range_scale(degree(temp_graph)),betweenness=range_scale(betweenness(temp_graph)),hub_score = range_scale(hub_score(temp_graph)$vector),group=species_group_16S[names(V(temp_graph)),7])
		temp_properties$study_name <- study_name
		assign(paste0("df_species_properties_elderly_",study_name),temp_properties)
		
		temp_kendall <- cor.fk(combined_16S_df_sum_stat_species[(combined_16S_df_sum_stat_species$study_name==study_name)&(combined_16S_df_sum_stat_species$age<60),select_species_16S])
		n_samples <- nrow(combined_16S_df_sum_stat_species[combined_16S_df_sum_stat_species$study_name==study_name,select_species_16S])
		temp_kendall <- apply(temp_kendall,2,function(x)(ifelse(is.nan(x),0,x)))
		assign(paste0("kendall_est_young_",study_name),temp_kendall)
		qval <- corr.p(temp_kendall,n=n_samples,adjust="fdr")
		temp_qval <- qval$p
		temp_qval <- apply(temp_qval,2,function(x)(ifelse(is.nan(x),0,x)))
		assign(paste0("kendall_qval_young_",study_name),temp_qval)
		temp_dir <- matrix(0,nrow(temp_kendall),ncol(temp_kendall))
		rownames(temp_dir) <- rownames(temp_kendall)
		colnames(temp_dir) <- colnames(temp_kendall)
		for(i in 1:nrow(temp_dir))
		{
			for(j in 1:ncol(temp_dir))
			{
				temp_dir[i,j] <- ifelse(temp_qval[i,j]<=0.01,sign(temp_qval[i,j]),0)
			}
		}
		assign(paste0("kendall_dir_young_",study_name),temp_dir)
		temp_graph <- graph_from_adjacency_matrix(apply(temp_dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper")
		assign(paste0("co_occurrence_graph_young_",study_name),graph_from_adjacency_matrix(apply(temp_dir,2,function(x)(ifelse(x>=1,1,0))),diag=FALSE,mode="upper"))
		temp_properties <- data.frame(species_name = names(V(temp_graph)),kendall_uniqueness = rem_16S_sum_stat_species_young$est[names(V(temp_graph)),1],bray_uniqueness = rem_16S_sum_stat_species_young$est[names(V(temp_graph)),2],jaccard_uniqueness = rem_16S_sum_stat_species_young$est[names(V(temp_graph)),3],manhattan_uniqueness = rem_16S_sum_stat_species_young$est[names(V(temp_graph)),4],shannon = rem_16S_sum_stat_species_elderly$est[names(V(temp_graph)),5],degree=range_scale(degree(temp_graph)),betweenness=range_scale(betweenness(temp_graph)),hub_score = range_scale(hub_score(temp_graph)$vector),group=species_group_16S[names(V(temp_graph)),7])
		temp_properties$study_name <- study_name
		assign(paste0("df_species_properties_young_",study_name),temp_properties)
	
	}
}

select_network_studies_shotgun <- intersect(selected_network_studies,sorted_study_rows_shotgun)
combined_df_species_properties_shotgun_young <- get(paste0("df_species_properties_young_",select_network_studies_shotgun[1]))
for(i in 2:length(select_network_studies_shotgun))
{
	print(select_network_studies_shotgun[i])
	temp_df <- get(paste0("df_species_properties_young_",select_network_studies_shotgun[i]))
	combined_df_species_properties_shotgun_young <- as.data.frame(rbind(combined_df_species_properties_shotgun_young,temp_df))
}

combined_df_species_properties_shotgun_elderly <- get(paste0("df_species_properties_elderly_",select_network_studies_shotgun[1]))
for(i in 2:length(select_network_studies_shotgun))
{
	print(select_network_studies_shotgun[i])
	temp_df <- get(paste0("df_species_properties_elderly_",select_network_studies_shotgun[i]))
	combined_df_species_properties_shotgun_elderly <- as.data.frame(rbind(combined_df_species_properties_shotgun_elderly,temp_df))
	temp_df <- NULL
}


select_network_studies_16S <- intersect(selected_network_studies,sorted_study_rows_16S)
combined_df_species_properties_16S_young <- get(paste0("df_species_properties_young_",select_network_studies_16S[1]))
for(i in 2:length(select_network_studies_16S))
{
	temp_df <- get(paste0("df_species_properties_young_",select_network_studies_16S[i]))
	combined_df_species_properties_16S_young <- as.data.frame(rbind(combined_df_species_properties_16S_young,temp_df))
}
combined_df_species_properties_16S_elderly <- get(paste0("df_species_properties_elderly_",select_network_studies_16S[1]))
for(i in 2:length(select_network_studies_16S))
{
	temp_df <- get(paste0("df_species_properties_elderly_",select_network_studies_16S[i]))
	combined_df_species_properties_16S_elderly <- as.data.frame(rbind(combined_df_species_properties_16S_elderly,temp_df))
}

combined_df_species_properties_young <- as.data.frame(rbind(combined_df_species_properties_shotgun_young,combined_df_species_properties_16S_young))

combined_df_species_properties_elderly <- as.data.frame(rbind(combined_df_species_properties_shotgun_elderly,combined_df_species_properties_16S_elderly))

common_species <- intersect(rownames(species_group_shotgun),rownames(species_group_16S))

degree_comp <- as.data.frame(matrix(NA,length(common_species),3))
rownames(degree_comp) <- common_species
colnames(degree_comp) <- c("direction","p_value","q_value")
for(i in 1:length(common_species))
{
	species_name <- common_species[i]
	degree_comp[species_name,1] <- mean(combined_df_species_properties_elderly[combined_df_species_properties_elderly$species_name==species_name,"degree"]-combined_df_species_properties_young[combined_df_species_properties_young$species_name==species_name,"degree"])
	temp_wilcox <- wilcox.test(combined_df_species_properties_elderly[combined_df_species_properties_elderly$species_name==species_name,"degree"],combined_df_species_properties_young[combined_df_species_properties_young$species_name==species_name,"degree"],paired=TRUE)
	degree_comp[species_name,2] <- temp_wilcox$p.value
}
degree_comp[,3] <- p.adjust(degree_comp[,2],method="fdr")

betweenness_comp <- as.data.frame(matrix(NA,length(common_species),3))
rownames(betweenness_comp) <- common_species
colnames(betweenness_comp) <- c("direction","p_value","q_value")
for(i in 1:length(common_species))
{
	species_name <- common_species[i]
	betweenness_comp[species_name,1] <- mean(combined_df_species_properties_elderly[combined_df_species_properties_elderly$species_name==species_name,"betweenness"]-combined_df_species_properties_young[combined_df_species_properties_young$species_name==species_name,"betweenness"])
	temp_wilcox <- wilcox.test(combined_df_species_properties_elderly[combined_df_species_properties_elderly$species_name==species_name,"betweenness"],combined_df_species_properties_young[combined_df_species_properties_young$species_name==species_name,"betweenness"],paired=TRUE)
	betweenness_comp[species_name,2] <- temp_wilcox$p.value
}
betweenness_comp[,3] <- p.adjust(betweenness_comp[,2],method="fdr")

hub_score_comp <- as.data.frame(matrix(NA,length(common_species),3))
rownames(hub_score_comp) <- common_species
colnames(hub_score_comp) <- c("direction","p_value","q_value")
for(i in 1:length(common_species))
{
	species_name <- common_species[i]
	hub_score_comp[species_name,1] <- mean(combined_df_species_properties_elderly[combined_df_species_properties_elderly$species_name==species_name,"hub_score"]-combined_df_species_properties_young[combined_df_species_properties_young$species_name==species_name,"hub_score"])
	temp_wilcox <- wilcox.test(combined_df_species_properties_elderly[combined_df_species_properties_elderly$species_name==species_name,"hub_score"],combined_df_species_properties_young[combined_df_species_properties_young$species_name==species_name,"hub_score"],paired=TRUE)
	hub_score_comp[species_name,2] <- temp_wilcox$p.value
}
hub_score_comp[,3] <- p.adjust(hub_score_comp[,2],method="fdr")

### G2: Elderly###
png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_shannon_betweenness_elderly.png"))
rem_G2_shannon_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"shannon","betweenness","study_name",selected_network_studies)
forest(rem_G2_shannon_betweenness_elderly$model)
title(paste0("G2: Shannon versus Betweenness (elderly): E: ",round(as.numeric(rem_G2_shannon_betweenness_elderly$model$beta),3),", P: ",format(round(rem_G2_shannon_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_kendall_uniqueness_betweenness_elderly.png"))
rem_G2_kendall_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"kendall_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G2_kendall_uniqueness_betweenness_elderly$model)
title(paste0("G2: Kendall Uniqueness v/s Betweenness (elderly): E: ",round(as.numeric(rem_G2_kendall_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_G2_kendall_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_manhattan_uniqueness_betweenness_elderly.png"))
rem_G2_manhattan_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"manhattan_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G2_manhattan_uniqueness_betweenness_elderly$model)
title(paste0("G2: Manhattan Uniqueness v/s Betweenness (elderly): E: ",round(as.numeric(rem_G2_manhattan_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_G2_manhattan_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_jaccard_uniqueness_betweenness_elderly.png"))
rem_G2_jaccard_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"jaccard_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G2_jaccard_uniqueness_betweenness_elderly$model)
title(paste0("G2: Jaccard Uniqueness v/s Betweenness (elderly): E: ",round(as.numeric(rem_G2_jaccard_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_G2_jaccard_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_bray_uniqueness_betweenness_elderly.png"))
rem_G2_bray_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"bray_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G2_bray_uniqueness_betweenness_elderly$model)
title(paste0("G2: Bray Uniqueness v/s Betweenness (elderly): E: ",round(as.numeric(rem_G2_bray_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_G2_bray_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_shannon_degree_elderly.png"))
rem_G2_shannon_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"shannon","degree","study_name",selected_network_studies)
forest(rem_G2_shannon_degree_elderly$model)
title(paste0("G2: Shannon v/s Degree (elderly): E: ",round(as.numeric(rem_G2_shannon_degree_elderly$model$beta),3),", P: ",format(round(rem_G2_shannon_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_kendall_uniqueness_degree_elderly.png"))
rem_G2_kendall_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"kendall_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G2_kendall_uniqueness_degree_elderly$model)
title(paste0("G2: Kendall Uniqueness v/s Degree (elderly): E: ",round(as.numeric(rem_G2_kendall_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_G2_kendall_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_manhattan_uniqueness_degree_elderly.png"))
rem_G2_manhattan_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"manhattan_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G2_manhattan_uniqueness_degree_elderly$model)
title(paste0("G2: Manhattan Uniqueness v/s Degree (elderly): E: ",round(as.numeric(rem_G2_manhattan_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_G2_manhattan_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_jaccard_uniqueness_degree_elderly.png"))
rem_G2_jaccard_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"jaccard_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G2_jaccard_uniqueness_degree_elderly$model)
title(paste0("G2: Jaccard Uniqueness v/s Degree (elderly): E: ",round(as.numeric(rem_G2_jaccard_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_G2_jaccard_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_bray_uniqueness_degree_elderly.png"))
rem_G2_bray_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"bray_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G2_bray_uniqueness_degree_elderly$model)
title(paste0("G2: Bray Uniqueness v/s Degree (elderly): E: ",round(as.numeric(rem_G2_bray_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_G2_bray_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_shannon_hub_score_elderly.png"))
rem_G2_shannon_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"shannon","hub_score","study_name",selected_network_studies)
forest(rem_G2_shannon_hub_score_elderly$model)
title(paste0("G2: Shannon v/s Hub-Score (elderly): E: ",round(as.numeric(rem_G2_shannon_hub_score_elderly$model$beta),3),", P: ",format(round(rem_G2_shannon_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_kendall_uniqueness_hub_score_elderly.png"))
rem_G2_kendall_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"kendall_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G2_kendall_uniqueness_hub_score_elderly$model)
title(paste0("G2: Kendall Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_G2_kendall_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_G2_kendall_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_manhattan_uniqueness_hub_score_elderly.png"))
rem_G2_manhattan_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"manhattan_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G2_manhattan_uniqueness_hub_score_elderly$model)
title(paste0("G2: Manhattan Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_G2_manhattan_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_G2_manhattan_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_jaccard_uniqueness_hub_score_elderly.png"))
rem_G2_jaccard_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"jaccard_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G2_jaccard_uniqueness_hub_score_elderly$model)
title(paste0("G2: Jaccard Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_G2_jaccard_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_G2_jaccard_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_bray_uniqueness_hub_score_elderly.png"))
rem_G2_bray_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G2",],"bray_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G2_bray_uniqueness_hub_score_elderly$model)
title(paste0("G2: Bray Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_G2_bray_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_G2_bray_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

### G4: Elderly###
png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_shannon_betweenness_elderly.png"))
rem_G4_shannon_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"shannon","betweenness","study_name",selected_network_studies)
forest(rem_G4_shannon_betweenness_elderly$model)
title(paste0("G4: Shannon versus Betweenness (elderly): E: ",round(as.numeric(rem_G4_shannon_betweenness_elderly$model$beta),3),", P: ",format(round(rem_G4_shannon_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_kendall_uniqueness_betweenness_elderly.png"))
rem_G4_kendall_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"kendall_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G4_kendall_uniqueness_betweenness_elderly$model)
title(paste0("G4: Kendall Uniqueness v/s Betweenness (elderly): E: ",round(as.numeric(rem_G4_kendall_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_G4_kendall_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_manhattan_uniqueness_betweenness_elderly.png"))
rem_G4_manhattan_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"manhattan_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G4_manhattan_uniqueness_betweenness_elderly$model)
title(paste0("G4: Manhattan Uniqueness v/s Betweenness (elderly): E: ",round(as.numeric(rem_G4_manhattan_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_G4_manhattan_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_jaccard_uniqueness_betweenness_elderly.png"))
rem_G4_jaccard_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"jaccard_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G4_jaccard_uniqueness_betweenness_elderly$model)
title(paste0("G4: Jaccard Uniqueness v/s Betweenness (elderly): E: ",round(as.numeric(rem_G4_jaccard_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_G4_jaccard_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_bray_uniqueness_betweenness_elderly.png"))
rem_G4_bray_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"bray_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G4_bray_uniqueness_betweenness_elderly$model)
title(paste0("G4: Bray Uniqueness v/s Betweenness (elderly): E: ",round(as.numeric(rem_G4_bray_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_G4_bray_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_shannon_degree_elderly.png"))
rem_G4_shannon_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"shannon","degree","study_name",selected_network_studies)
forest(rem_G4_shannon_degree_elderly$model)
title(paste0("G4: Shannon v/s Degree (elderly): E: ",round(as.numeric(rem_G4_shannon_degree_elderly$model$beta),3),", P: ",format(round(rem_G4_shannon_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_kendall_uniqueness_degree_elderly.png"))
rem_G4_kendall_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"kendall_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G4_kendall_uniqueness_degree_elderly$model)
title(paste0("G4: Kendall Uniqueness v/s Degree (elderly): E: ",round(as.numeric(rem_G4_kendall_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_G4_kendall_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_manhattan_uniqueness_degree_elderly.png"))
rem_G4_manhattan_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"manhattan_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G4_manhattan_uniqueness_degree_elderly$model)
title(paste0("G4: Manhattan Uniqueness v/s Degree (elderly): E: ",round(as.numeric(rem_G4_manhattan_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_G4_manhattan_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_jaccard_uniqueness_degree_elderly.png"))
rem_G4_jaccard_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"jaccard_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G4_jaccard_uniqueness_degree_elderly$model)
title(paste0("G4: Jaccard Uniqueness v/s Degree (elderly): E: ",round(as.numeric(rem_G4_jaccard_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_G4_jaccard_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_bray_uniqueness_degree_elderly.png"))
rem_G4_bray_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"bray_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G4_bray_uniqueness_degree_elderly$model)
title(paste0("G4: Bray Uniqueness v/s Degree (elderly): E: ",round(as.numeric(rem_G4_bray_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_G4_bray_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_shannon_hub_score_elderly.png"))
rem_G4_shannon_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"shannon","hub_score","study_name",selected_network_studies)
forest(rem_G4_shannon_hub_score_elderly$model)
title(paste0("G4: Shannon v/s Hub-Score (elderly): E: ",round(as.numeric(rem_G4_shannon_hub_score_elderly$model$beta),3),", P: ",format(round(rem_G4_shannon_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_kendall_uniqueness_hub_score_elderly.png"))
rem_G4_kendall_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"kendall_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G4_kendall_uniqueness_hub_score_elderly$model)
title(paste0("G4: Kendall Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_G4_kendall_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_G4_kendall_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_manhattan_uniqueness_hub_score_elderly.png"))
rem_G4_manhattan_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"manhattan_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G4_manhattan_uniqueness_hub_score_elderly$model)
title(paste0("G4: Manhattan Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_G4_manhattan_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_G4_manhattan_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_jaccard_uniqueness_hub_score_elderly.png"))
rem_G4_jaccard_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"jaccard_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G4_jaccard_uniqueness_hub_score_elderly$model)
title(paste0("G4: Jaccard Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_G4_jaccard_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_G4_jaccard_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_bray_uniqueness_hub_score_elderly.png"))
rem_G4_bray_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[combined_df_species_properties_elderly$group=="G4",],"bray_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G4_bray_uniqueness_hub_score_elderly$model)
title(paste0("G4: Bray Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_G4_bray_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_G4_bray_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

### Others: Elderly###
png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_shannon_betweenness_elderly.png"))
rem_Others_shannon_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"shannon","betweenness","study_name",selected_network_studies)
forest(rem_Others_shannon_betweenness_elderly$model)
title(paste0("Others: Shannon versus Betweenness (elderly): E: ",round(as.numeric(rem_Others_shannon_betweenness_elderly$model$beta),3),", P: ",format(round(rem_Others_shannon_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_kendall_uniqueness_betweenness_elderly.png"))
rem_Others_kendall_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"kendall_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_Others_kendall_uniqueness_betweenness_elderly$model)
title(paste0("Others: Kendall Uniqueness v/s Betweenness (elderly): E: ",round(as.numeric(rem_Others_kendall_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_Others_kendall_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_manhattan_uniqueness_betweenness_elderly.png"))
rem_Others_manhattan_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"manhattan_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_Others_manhattan_uniqueness_betweenness_elderly$model)
title(paste0("Others: Manhattan Uniqueness v/s Betweenness (elderly): E: ",round(as.numeric(rem_Others_manhattan_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_Others_manhattan_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_jaccard_uniqueness_betweenness_elderly.png"))
rem_Others_jaccard_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"jaccard_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_Others_jaccard_uniqueness_betweenness_elderly$model)
title(paste0("Others: Jaccard Uniqueness v/s Betweenness (elderly): E: ",round(as.numeric(rem_Others_jaccard_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_Others_jaccard_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_bray_uniqueness_betweenness_elderly.png"))
rem_Others_bray_uniqueness_betweenness_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"bray_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_Others_bray_uniqueness_betweenness_elderly$model)
title(paste0("Others: Bray Unique. v/s Betweenness (elderly): E: ",round(as.numeric(rem_Others_bray_uniqueness_betweenness_elderly$model$beta),3),", P: ",format(round(rem_Others_bray_uniqueness_betweenness_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_shannon_degree_elderly.png"))
rem_Others_shannon_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"shannon","degree","study_name",selected_network_studies)
forest(rem_Others_shannon_degree_elderly$model)
title(paste0("Others: Shannon v/s Degree (elderly): E: ",round(as.numeric(rem_Others_shannon_degree_elderly$model$beta),3),", P: ",format(round(rem_Others_shannon_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_kendall_uniqueness_degree_elderly.png"))
rem_Others_kendall_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"kendall_uniqueness","degree","study_name",selected_network_studies)
forest(rem_Others_kendall_uniqueness_degree_elderly$model)
title(paste0("Others: Kendall Uniqueness v/s Degree (elderly): E: ",round(as.numeric(rem_Others_kendall_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_Others_kendall_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_manhattan_uniqueness_degree_elderly.png"))
rem_Others_manhattan_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"manhattan_uniqueness","degree","study_name",selected_network_studies)
forest(rem_Others_manhattan_uniqueness_degree_elderly$model)
title(paste0("Others: Manhattan Uniqueness v/s Degree (elderly): E: ",round(as.numeric(rem_Others_manhattan_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_Others_manhattan_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_jaccard_uniqueness_degree_elderly.png"))
rem_Others_jaccard_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"jaccard_uniqueness","degree","study_name",selected_network_studies)
forest(rem_Others_jaccard_uniqueness_degree_elderly$model)
title(paste0("Others: Jaccard Uniqueness v/s Degree (elderly): E: ",round(as.numeric(rem_Others_jaccard_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_Others_jaccard_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_bray_uniqueness_degree_elderly.png"))
rem_Others_bray_uniqueness_degree_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"bray_uniqueness","degree","study_name",selected_network_studies)
forest(rem_Others_bray_uniqueness_degree_elderly$model)
title(paste0("Others: Bray Unique. v/s Degree (elderly): E: ",round(as.numeric(rem_Others_bray_uniqueness_degree_elderly$model$beta),3),", P: ",format(round(rem_Others_bray_uniqueness_degree_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_shannon_hub_score_elderly.png"))
rem_Others_shannon_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"shannon","hub_score","study_name",selected_network_studies)
forest(rem_Others_shannon_hub_score_elderly$model)
title(paste0("Others: Shannon v/s Hub-Score (elderly): E: ",round(as.numeric(rem_Others_shannon_hub_score_elderly$model$beta),3),", P: ",format(round(rem_Others_shannon_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_kendall_uniqueness_hub_score_elderly.png"))
rem_Others_kendall_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"kendall_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_Others_kendall_uniqueness_hub_score_elderly$model)
title(paste0("Others: Kendall Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_Others_kendall_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_Others_kendall_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_manhattan_uniqueness_hub_score_elderly.png"))
rem_Others_manhattan_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"manhattan_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_Others_manhattan_uniqueness_hub_score_elderly$model)
title(paste0("Others: Manhattan Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_Others_manhattan_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_Others_manhattan_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_jaccard_uniqueness_hub_score_elderly.png"))
rem_Others_jaccard_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"jaccard_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_Others_jaccard_uniqueness_hub_score_elderly$model)
title(paste0("Others: Jaccard Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_Others_jaccard_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_Others_jaccard_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_bray_uniqueness_hub_score_elderly.png"))
rem_Others_bray_uniqueness_hub_score_elderly <- compute_meta_lm(combined_df_species_properties_elderly[!(combined_df_species_properties_elderly$group %in% c("G2","G4")),],"bray_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_Others_bray_uniqueness_hub_score_elderly$model)
title(paste0("Others: Bray Uniqueness v/s Hub-Score (elderly): E: ",round(as.numeric(rem_Others_bray_uniqueness_hub_score_elderly$model$beta),3),", P: ",format(round(rem_Others_bray_uniqueness_hub_score_elderly$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

### G2: Young###
png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_shannon_betweenness_young.png"))
rem_G2_shannon_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"shannon","betweenness","study_name",selected_network_studies)
forest(rem_G2_shannon_betweenness_young$model)
title(paste0("G2: Shannon versus Betweenness (young): E: ",round(as.numeric(rem_G2_shannon_betweenness_young$model$beta),3),", P: ",format(round(rem_G2_shannon_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_kendall_uniqueness_betweenness_young.png"))
rem_G2_kendall_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"kendall_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G2_kendall_uniqueness_betweenness_young$model)
title(paste0("G2: Kendall Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_G2_kendall_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_G2_kendall_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_manhattan_uniqueness_betweenness_young.png"))
rem_G2_manhattan_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"manhattan_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G2_manhattan_uniqueness_betweenness_young$model)
title(paste0("G2: Manhattan Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_G2_manhattan_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_G2_manhattan_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_jaccard_uniqueness_betweenness_young.png"))
rem_G2_jaccard_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"jaccard_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G2_jaccard_uniqueness_betweenness_young$model)
title(paste0("G2: Jaccard Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_G2_jaccard_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_G2_jaccard_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_bray_uniqueness_betweenness_young.png"))
rem_G2_bray_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"bray_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G2_bray_uniqueness_betweenness_young$model)
title(paste0("G2: Bray Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_G2_bray_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_G2_bray_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_shannon_degree_young.png"))
rem_G2_shannon_degree_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"shannon","degree","study_name",selected_network_studies)
forest(rem_G2_shannon_degree_young$model)
title(paste0("G2: Shannon v/s Degree (young): E: ",round(as.numeric(rem_G2_shannon_degree_young$model$beta),3),", P: ",format(round(rem_G2_shannon_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_kendall_uniqueness_degree_young.png"))
rem_G2_kendall_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"kendall_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G2_kendall_uniqueness_degree_young$model)
title(paste0("G2: Kendall Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_G2_kendall_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_G2_kendall_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_manhattan_uniqueness_degree_young.png"))
rem_G2_manhattan_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"manhattan_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G2_manhattan_uniqueness_degree_young$model)
title(paste0("G2: Manhattan Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_G2_manhattan_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_G2_manhattan_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_jaccard_uniqueness_degree_young.png"))
rem_G2_jaccard_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"jaccard_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G2_jaccard_uniqueness_degree_young$model)
title(paste0("G2: Jaccard Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_G2_jaccard_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_G2_jaccard_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_bray_uniqueness_degree_young.png"))
rem_G2_bray_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"bray_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G2_bray_uniqueness_degree_young$model)
title(paste0("G2: Bray Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_G2_bray_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_G2_bray_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_shannon_hub_score_young.png"))
rem_G2_shannon_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"shannon","hub_score","study_name",selected_network_studies)
forest(rem_G2_shannon_hub_score_young$model)
title(paste0("G2: Shannon v/s Hub-Score (young): E: ",round(as.numeric(rem_G2_shannon_hub_score_young$model$beta),3),", P: ",format(round(rem_G2_shannon_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_kendall_uniqueness_hub_score_young.png"))
rem_G2_kendall_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"kendall_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G2_kendall_uniqueness_hub_score_young$model)
title(paste0("G2: Kendall Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_G2_kendall_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_G2_kendall_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_manhattan_uniqueness_hub_score_young.png"))
rem_G2_manhattan_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"manhattan_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G2_manhattan_uniqueness_hub_score_young$model)
title(paste0("G2: Manhattan Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_G2_manhattan_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_G2_manhattan_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_jaccard_uniqueness_hub_score_young.png"))
rem_G2_jaccard_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"jaccard_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G2_jaccard_uniqueness_hub_score_young$model)
title(paste0("G2: Jaccard Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_G2_jaccard_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_G2_jaccard_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G2_bray_uniqueness_hub_score_young.png"))
rem_G2_bray_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G2",],"bray_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G2_bray_uniqueness_hub_score_young$model)
title(paste0("G2: Bray Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_G2_bray_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_G2_bray_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

### G4: Young###
png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_shannon_betweenness_young.png"))
rem_G4_shannon_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"shannon","betweenness","study_name",selected_network_studies)
forest(rem_G4_shannon_betweenness_young$model)
title(paste0("G4: Shannon versus Betweenness (young): E: ",round(as.numeric(rem_G4_shannon_betweenness_young$model$beta),3),", P: ",format(round(rem_G4_shannon_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_kendall_uniqueness_betweenness_young.png"))
rem_G4_kendall_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"kendall_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G4_kendall_uniqueness_betweenness_young$model)
title(paste0("G4: Kendall Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_G4_kendall_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_G4_kendall_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_manhattan_uniqueness_betweenness_young.png"))
rem_G4_manhattan_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"manhattan_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G4_manhattan_uniqueness_betweenness_young$model)
title(paste0("G4: Manhattan Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_G4_manhattan_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_G4_manhattan_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_jaccard_uniqueness_betweenness_young.png"))
rem_G4_jaccard_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"jaccard_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G4_jaccard_uniqueness_betweenness_young$model)
title(paste0("G4: Jaccard Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_G4_jaccard_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_G4_jaccard_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_bray_uniqueness_betweenness_young.png"))
rem_G4_bray_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"bray_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_G4_bray_uniqueness_betweenness_young$model)
title(paste0("G4: Bray Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_G4_bray_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_G4_bray_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_shannon_degree_young.png"))
rem_G4_shannon_degree_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"shannon","degree","study_name",selected_network_studies)
forest(rem_G4_shannon_degree_young$model)
title(paste0("G4: Shannon v/s Degree (young): E: ",round(as.numeric(rem_G4_shannon_degree_young$model$beta),3),", P: ",format(round(rem_G4_shannon_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_kendall_uniqueness_degree_young.png"))
rem_G4_kendall_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"kendall_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G4_kendall_uniqueness_degree_young$model)
title(paste0("G4: Kendall Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_G4_kendall_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_G4_kendall_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_manhattan_uniqueness_degree_young.png"))
rem_G4_manhattan_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"manhattan_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G4_manhattan_uniqueness_degree_young$model)
title(paste0("G4: Manhattan Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_G4_manhattan_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_G4_manhattan_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_jaccard_uniqueness_degree_young.png"))
rem_G4_jaccard_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"jaccard_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G4_jaccard_uniqueness_degree_young$model)
title(paste0("G4: Jaccard Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_G4_jaccard_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_G4_jaccard_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_bray_uniqueness_degree_young.png"))
rem_G4_bray_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"bray_uniqueness","degree","study_name",selected_network_studies)
forest(rem_G4_bray_uniqueness_degree_young$model)
title(paste0("G4: Bray Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_G4_bray_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_G4_bray_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_shannon_hub_score_young.png"))
rem_G4_shannon_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"shannon","hub_score","study_name",selected_network_studies)
forest(rem_G4_shannon_hub_score_young$model)
title(paste0("G4: Shannon v/s Hub-Score (young): E: ",round(as.numeric(rem_G4_shannon_hub_score_young$model$beta),3),", P: ",format(round(rem_G4_shannon_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_kendall_uniqueness_hub_score_young.png"))
rem_G4_kendall_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"kendall_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G4_kendall_uniqueness_hub_score_young$model)
title(paste0("G4: Kendall Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_G4_kendall_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_G4_kendall_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_manhattan_uniqueness_hub_score_young.png"))
rem_G4_manhattan_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"manhattan_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G4_manhattan_uniqueness_hub_score_young$model)
title(paste0("G4: Manhattan Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_G4_manhattan_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_G4_manhattan_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_jaccard_uniqueness_hub_score_young.png"))
rem_G4_jaccard_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"jaccard_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G4_jaccard_uniqueness_hub_score_young$model)
title(paste0("G4: Jaccard Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_G4_jaccard_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_G4_jaccard_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_G4_bray_uniqueness_hub_score_young.png"))
rem_G4_bray_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[combined_df_species_properties_young$group=="G4",],"bray_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_G4_bray_uniqueness_hub_score_young$model)
title(paste0("G4: Bray Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_G4_bray_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_G4_bray_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

### Others: Young###
png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_shannon_betweenness_young.png"))
rem_Others_shannon_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"shannon","betweenness","study_name",selected_network_studies)
forest(rem_Others_shannon_betweenness_young$model)
title(paste0("Others: Shannon versus Betweenness (young): E: ",round(as.numeric(rem_Others_shannon_betweenness_young$model$beta),3),", P: ",format(round(rem_Others_shannon_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_kendall_uniqueness_betweenness_young.png"))
rem_Others_kendall_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"kendall_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_Others_kendall_uniqueness_betweenness_young$model)
title(paste0("Others: Kendall Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_Others_kendall_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_Others_kendall_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_manhattan_uniqueness_betweenness_young.png"))
rem_Others_manhattan_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"manhattan_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_Others_manhattan_uniqueness_betweenness_young$model)
title(paste0("Others: Manhattan Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_Others_manhattan_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_Others_manhattan_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_jaccard_uniqueness_betweenness_young.png"))
rem_Others_jaccard_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"jaccard_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_Others_jaccard_uniqueness_betweenness_young$model)
title(paste0("Others: Jaccard Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_Others_jaccard_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_Others_jaccard_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_bray_uniqueness_betweenness_young.png"))
rem_Others_bray_uniqueness_betweenness_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"bray_uniqueness","betweenness","study_name",selected_network_studies)
forest(rem_Others_bray_uniqueness_betweenness_young$model)
title(paste0("Others: Bray Uniqueness v/s Betweenness (young): E: ",round(as.numeric(rem_Others_bray_uniqueness_betweenness_young$model$beta),3),", P: ",format(round(rem_Others_bray_uniqueness_betweenness_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_shannon_degree_young.png"))
rem_Others_shannon_degree_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"shannon","degree","study_name",selected_network_studies)
forest(rem_Others_shannon_degree_young$model)
title(paste0("Others: Shannon v/s Degree (young): E: ",round(as.numeric(rem_Others_shannon_degree_young$model$beta),3),", P: ",format(round(rem_Others_shannon_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_kendall_uniqueness_degree_young.png"))
rem_Others_kendall_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"kendall_uniqueness","degree","study_name",selected_network_studies)
forest(rem_Others_kendall_uniqueness_degree_young$model)
title(paste0("Others: Kendall Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_Others_kendall_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_Others_kendall_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_manhattan_uniqueness_degree_young.png"))
rem_Others_manhattan_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"manhattan_uniqueness","degree","study_name",selected_network_studies)
forest(rem_Others_manhattan_uniqueness_degree_young$model)
title(paste0("Others: Manhattan Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_Others_manhattan_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_Others_manhattan_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_jaccard_uniqueness_degree_young.png"))
rem_Others_jaccard_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"jaccard_uniqueness","degree","study_name",selected_network_studies)
forest(rem_Others_jaccard_uniqueness_degree_young$model)
title(paste0("Others: Jaccard Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_Others_jaccard_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_Others_jaccard_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_bray_uniqueness_degree_young.png"))
rem_Others_bray_uniqueness_degree_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"bray_uniqueness","degree","study_name",selected_network_studies)
forest(rem_Others_bray_uniqueness_degree_young$model)
title(paste0("Others: Bray Uniqueness v/s Degree (young): E: ",round(as.numeric(rem_Others_bray_uniqueness_degree_young$model$beta),3),", P: ",format(round(rem_Others_bray_uniqueness_degree_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_shannon_hub_score_young.png"))
rem_Others_shannon_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"shannon","hub_score","study_name",selected_network_studies)
forest(rem_Others_shannon_hub_score_young$model)
title(paste0("Others: Shannon v/s Hub-Score (young): E: ",round(as.numeric(rem_Others_shannon_hub_score_young$model$beta),3),", P: ",format(round(rem_Others_shannon_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_kendall_uniqueness_hub_score_young.png"))
rem_Others_kendall_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"kendall_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_Others_kendall_uniqueness_hub_score_young$model)
title(paste0("Others: Kendall Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_Others_kendall_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_Others_kendall_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_manhattan_uniqueness_hub_score_young.png"))
rem_Others_manhattan_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"manhattan_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_Others_manhattan_uniqueness_hub_score_young$model)
title(paste0("Others: Manhattan Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_Others_manhattan_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_Others_manhattan_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_jaccard_uniqueness_hub_score_young.png"))
rem_Others_jaccard_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"jaccard_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_Others_jaccard_uniqueness_hub_score_young$model)
title(paste0("Others: Jaccard Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_Others_jaccard_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_Others_jaccard_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

png(paste0("C:\\Projects\\ELDERMET\\NatureAgingRevision\\GroupNetwork\\Centrality_To_Properties\\","rem_Others_bray_uniqueness_hub_score_young.png"))
rem_Others_bray_uniqueness_hub_score_young <- compute_meta_lm(combined_df_species_properties_young[!(combined_df_species_properties_young$group %in% c("G2","G4")),],"bray_uniqueness","hub_score","study_name",selected_network_studies)
forest(rem_Others_bray_uniqueness_hub_score_young$model)
title(paste0("Others: Bray Uniqueness v/s Hub-Score (young): E: ",round(as.numeric(rem_Others_bray_uniqueness_hub_score_young$model$beta),3),", P: ",format(round(rem_Others_bray_uniqueness_hub_score_young$model$pval,4),scientific=TRUE)),cex=0.5)
dev.off()

group_color <- NULL
group_color["G1"] <- "blue"
group_color["G2"] <- "green"
group_color["G3"] <- "yellow"
group_color["G4"] <- "red"

degree_to_species_elderly <- as.data.frame(matrix(0,length(common_species),length(selected_network_studies)))
rownames(degree_to_species_elderly) <- common_species
colnames(degree_to_species_elderly) <- selected_network_studies

for(i in 1:length(common_species))
{
	species_name <- common_species[i]
	for(j in 1:length(selected_network_studies))
	{
		study_name = selected_network_studies[j]
		degree_to_species_elderly[species_name,study_name] <- combined_df_species_properties_elderly[(combined_df_species_properties_elderly$species_name == species_name)&(combined_df_species_properties_elderly$study_name == study_name),"degree"]
	}
}

degree_to_species_young <- as.data.frame(matrix(0,length(common_species),length(selected_network_studies)))
rownames(degree_to_species_young) <- common_species
colnames(degree_to_species_young) <- selected_network_studies

for(i in 1:length(common_species))
{
	species_name <- common_species[i]
	for(j in 1:length(selected_network_studies))
	{
		study_name = selected_network_studies[j]
		degree_to_species_young[species_name,study_name] <- combined_df_species_properties_young[(combined_df_species_properties_young$species_name == species_name)&(combined_df_species_properties_young$study_name == study_name),"degree"]
	}
}

betweenness_to_species_elderly <- as.data.frame(matrix(0,length(common_species),length(selected_network_studies)))
rownames(betweenness_to_species_elderly) <- common_species
colnames(betweenness_to_species_elderly) <- selected_network_studies

for(i in 1:length(common_species))
{
	species_name <- common_species[i]
	for(j in 1:length(selected_network_studies))
	{
		study_name = selected_network_studies[j]
		betweenness_to_species_elderly[species_name,study_name] <- combined_df_species_properties_elderly[(combined_df_species_properties_elderly$species_name == species_name)&(combined_df_species_properties_elderly$study_name == study_name),"betweenness"]
	}
}

betweenness_to_species_young <- as.data.frame(matrix(0,length(common_species),length(selected_network_studies)))
rownames(betweenness_to_species_young) <- common_species
colnames(betweenness_to_species_young) <- selected_network_studies

for(i in 1:length(common_species))
{
	species_name <- common_species[i]
	for(j in 1:length(selected_network_studies))
	{
		study_name = selected_network_studies[j]
		betweenness_to_species_young[species_name,study_name] <- combined_df_species_properties_young[(combined_df_species_properties_young$species_name == species_name)&(combined_df_species_properties_young$study_name == study_name),"betweenness"]
	}
}

hub_score_to_species_elderly <- as.data.frame(matrix(0,length(common_species),length(selected_network_studies)))
rownames(hub_score_to_species_elderly) <- common_species
colnames(hub_score_to_species_elderly) <- selected_network_studies

for(i in 1:length(common_species))
{
	species_name <- common_species[i]
	for(j in 1:length(selected_network_studies))
	{
		study_name = selected_network_studies[j]
		hub_score_to_species_elderly[species_name,study_name] <- combined_df_species_properties_elderly[(combined_df_species_properties_elderly$species_name == species_name)&(combined_df_species_properties_elderly$study_name == study_name),"hub_score"]
	}
}

hub_score_to_species_young <- as.data.frame(matrix(0,length(common_species),length(selected_network_studies)))
rownames(hub_score_to_species_young) <- common_species
colnames(hub_score_to_species_young) <- selected_network_studies

for(i in 1:length(common_species))
{
	species_name <- common_species[i]
	for(j in 1:length(selected_network_studies))
	{
		study_name = selected_network_studies[j]
		hub_score_to_species_young[species_name,study_name] <- combined_df_species_properties_young[(combined_df_species_properties_young$species_name == species_name)&(combined_df_species_properties_young$study_name == study_name),"hub_score"]
	}
}

#heatmap.2(t(degree_to_species_young),density="none",trace="none",lhei=c(1,5),lwid=c(0.5,5),margins=c(20,7),col=rev(brewer.pal(8,"RdYlBu")))
#heatmap.2(t(degree_to_species_young),density="none",trace="none",lhei=c(1,5),lwid=c(0.5,5),margins=c(13,7),ColSideColors=group_color[species_group_shotgun[rownames(degree_to_species_young),7]],col=rev(brewer.pal(8,"RdYlBu")))

#heatmap.2(t(degree_to_species_elderly),density="none",trace="none",lhei=c(1,5),lwid=c(0.5,5),margins=c(20,7),col=rev(brewer.pal(8,"RdYlBu")))
#heatmap.2(t(degree_to_species_elderly),density="none",trace="none",lhei=c(1,5),lwid=c(0.5,5),margins=c(13,7),ColSideColors=group_color[species_group_shotgun[rownames(degree_to_species_elderly),7]],col=rev(brewer.pal(8,"RdYlBu")))
