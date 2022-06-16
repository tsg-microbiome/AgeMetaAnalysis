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
		print(group)
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




load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\nuage_stage_2c_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ag_stage_2c_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\he_stage_2c_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\odamaki_stage_2c_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\logmpie_stage_2c_results.RData")


two_cohort_association_direction <- merge(ag_association_direction,he_association_direction,by="row.names",all=TRUE)[,-1]
rownames(two_cohort_association_direction) <- merge(ag_association_direction,he_association_direction,by="row.names",all=TRUE)[,1]
two_cohort_association_direction <- apply(two_cohort_association_direction,2,function(x)(ifelse(is.na(x),0,x)))

association_direction_16S <- merge(two_cohort_association_direction,nuage_association_direction,by="row.names",all=TRUE)[,-1]
rownames(association_direction_16S) <- merge(two_cohort_association_direction,nuage_association_direction,by="row.names",all=TRUE)[,1]
association_direction_16S <- apply(association_direction_16S,2,function(x)(ifelse(is.na(x),0,x)))

PositiveWithUnhealthyAging_16S <- names(which((rowSums((apply(association_direction_16S,2,function(x)(ifelse(x>=2,1,0))))[,c(2,4,6)])>=2)&(rowSums((apply(association_direction_16S,2,function(x)(ifelse(x>=2,1,0))))[,c(1,3,5)])==0)))
NegativeWithUnhealthyAging_16S <- names(which((rowSums((apply(association_direction_16S,2,function(x)(ifelse(x>=2,1,0))))[,c(2,4,6)])==0)&(rowSums((apply(association_direction_16S,2,function(x)(ifelse(x>=2,1,0))))[,c(1,3,5)])>=2)))
ArrangedRankedMarkers_16S <- c(PositiveWithUnhealthyAging_16S,NegativeWithUnhealthyAging_16S)[order(rowSums(association_direction_16S[c(PositiveWithUnhealthyAging_16S,NegativeWithUnhealthyAging_16S),c(1,3,5)])-rowSums(association_direction_16S[c(PositiveWithUnhealthyAging_16S,NegativeWithUnhealthyAging_16S),c(2,4,6)]))]

colnames(he_rlm_dir_unhealthy_aging) <- paste0("he_et_al","_",colnames(he_rlm_dir_unhealthy_aging))
colnames(nuage_rlm_dir_unhealthy_aging) <- paste0("nuage","_",colnames(nuage_rlm_dir_unhealthy_aging))
colnames(ag_rlm_dir_unhealthy_aging) <- paste0("ag","_",colnames(ag_rlm_dir_unhealthy_aging))

combined_dir_unhealthy_aging_16S <- as.data.frame(cbind(he_rlm_dir_unhealthy_aging[ArrangedRankedMarkers_16S,],nuage_rlm_dir_unhealthy_aging[ArrangedRankedMarkers_16S,],ag_rlm_dir_unhealthy_aging[ArrangedRankedMarkers_16S,]))

#heatmap.2(apply(combined_dir_unhealthy_aging_16S,2,function(x)(ifelse(x<0,x-0.1,x))),density="none",trace="none",col=c("blue4","blue3","white","white","white","orange","red"),Colv=FALSE,Rowv=FALSE,margins=c(10,10),lhei=c(0.5,5),cexRow=0.5,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(combined_dir_unhealthy_aging_16S),rowsep=1:nrow(combined_dir_unhealthy_aging_16S))

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\em_stage_2c_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_stage_2c_results.RData")

colnames(cmd3_rlm_dir_unhealthy_aging) <- paste0("cmd3","_",colnames(cmd3_rlm_dir_unhealthy_aging))
colnames(em_rlm_dir_unhealthy_aging) <- paste0("em","_",colnames(em_rlm_dir_unhealthy_aging))

association_direction_shotgun <- as.data.frame(cbind(cmd3_association_direction[intersect(rownames(cmd3_association_direction),rownames(em_association_direction)),],em_association_direction[intersect(rownames(cmd3_association_direction),rownames(em_association_direction)),]))

colnames(association_direction_shotgun) <- c("decreased_cmd3","increased_cmd3","increased_em","decreased_em")

PositiveWithUnhealthyAging_shotgun <- names(which(rowSums(apply(association_direction_shotgun[,c(2,4)],2,function(x)(ifelse(x>=2,1,0))))==2))
NegativeWithUnhealthyAging_shotgun <- names(which(rowSums(apply(association_direction_shotgun[,c(1,3)],2,function(x)(ifelse(x>=2,1,0))))==2))

ArrangedRankedMarkers_shotgun <- c(PositiveWithUnhealthyAging_shotgun,NegativeWithUnhealthyAging_shotgun)[order(rowSums(association_direction_shotgun[c(PositiveWithUnhealthyAging_shotgun,NegativeWithUnhealthyAging_shotgun),c(1,3)])-rowSums(association_direction_shotgun[c(PositiveWithUnhealthyAging_shotgun,NegativeWithUnhealthyAging_shotgun),c(2,4)]))]

combined_dir_unhealthy_aging_shotgun <- as.data.frame(cbind(cmd3_rlm_dir_unhealthy_aging[ArrangedRankedMarkers_shotgun,],em_rlm_dir_unhealthy_aging[ArrangedRankedMarkers_shotgun,]))

#heatmap.2(apply(combined_dir_unhealthy_aging_shotgun,2,function(x)(ifelse(x<0,x-0.1,x))),density="none",trace="none",col=c("blue4","blue3","white","white","white","orange","red"),Colv=FALSE,Rowv=FALSE,margins=c(10,10),lhei=c(0.5,5),cexRow=0.5,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(combined_dir_unhealthy_aging_shotgun),rowsep=1:nrow(combined_dir_unhealthy_aging_shotgun))

shotgun_association_direction_young <- merge(cmd3_association_direction_young,isc_association_direction_young,by="row.names",all=TRUE)[,-1]
rownames(shotgun_association_direction_young) <- merge(cmd3_association_direction_young,isc_association_direction_young,by="row.names",all=TRUE)[,1]
shotgun_association_direction_young <- apply(shotgun_association_direction_young,2,function(x)(ifelse(is.na(x),0,x)))

PositiveWithUnhealthyYoung_shotgun <- rownames(shotgun_association_direction_young[(shotgun_association_direction_young[,2]>=2)&(shotgun_association_direction_young[,1]<=2)&(shotgun_association_direction_young[,4]==1),])
NegativeWithUnhealthyYoung_shotgun <- rownames(shotgun_association_direction_young[(shotgun_association_direction_young[,1]>=2)&(shotgun_association_direction_young[,2]<=2)&(shotgun_association_direction_young[,3]==1),])

ArrangedRankedMarkers_shotgun_young <- c(PositiveWithUnhealthyYoung_shotgun,NegativeWithUnhealthyYoung_shotgun)[order(rowSums(shotgun_association_direction_young[c(PositiveWithUnhealthyYoung_shotgun,NegativeWithUnhealthyYoung_shotgun),c(1,3)])-rowSums(shotgun_association_direction_young[c(PositiveWithUnhealthyYoung_shotgun,NegativeWithUnhealthyYoung_shotgun),c(2,4)]))]

combined_dir_unhealthy_young_shotgun <- as.data.frame(cbind(cmd3_rlm_dir_unhealthy_young[ArrangedRankedMarkers_shotgun_young,],isc_rlm_dir_unhealthy_young[ArrangedRankedMarkers_shotgun_young,]))
colnames(combined_dir_unhealthy_young_shotgun)[11] <- "ISC_IBS"

#heatmap.2(apply(combined_dir_unhealthy_young_shotgun,2,function(x)(ifelse(x<0,x-0.1,x))),density="none",trace="none",col=c("blue4","blue3","white","white","white","orange","red"),Colv=FALSE,Rowv=FALSE,margins=c(10,10),lhei=c(0.5,5),cexRow=0.5,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(combined_dir_unhealthy_young_shotgun),rowsep=1:nrow(combined_dir_unhealthy_young_shotgun))

association_direction_16S_young <- merge(ag_association_direction_young,he_association_direction_young,by="row.names",all=TRUE)[,-1]
rownames(association_direction_16S_young) <- merge(ag_association_direction_young,he_association_direction_young,by="row.names",all=TRUE)[,1]
association_direction_16S_young <- apply(association_direction_16S_young,2,function(x)(ifelse(is.na(x),0,x)))

PositiveWithUnhealthyYoung_16S <- rownames(association_direction_16S_young[(association_direction_16S_young[,2]>=2)&(association_direction_16S_young[,4]>=2)&(association_direction_16S_young[,1]<=1)&(association_direction_16S_young[,3]<=1),])
NegativeWithUnhealthyYoung_16S <- rownames(association_direction_16S_young[(association_direction_16S_young[,2]<=1)&(association_direction_16S_young[,4]<=1)&(association_direction_16S_young[,1]>=2)&(association_direction_16S_young[,3]>=2),])

ArrangedRankedMarkers_16S_young <- c(PositiveWithUnhealthyYoung_16S,NegativeWithUnhealthyYoung_16S)[order(rowSums(association_direction_16S_young[c(PositiveWithUnhealthyYoung_16S,NegativeWithUnhealthyYoung_16S),c(1,3)])-rowSums(association_direction_16S_young[c(PositiveWithUnhealthyYoung_16S,NegativeWithUnhealthyYoung_16S),c(2,4)]))]

combined_dir_unhealthy_young_16S <- as.data.frame(cbind(ag_rlm_dir_unhealthy_young[ArrangedRankedMarkers_16S_young,],he_rlm_dir_unhealthy_young[ArrangedRankedMarkers_16S_young,]))

#heatmap.2(apply(combined_dir_unhealthy_young_16S,2,function(x)(ifelse(x<0,x-0.1,x))),density="none",trace="none",col=c("blue4","blue3","white","white","white","orange","red"),Colv=FALSE,Rowv=FALSE,margins=c(10,10),lhei=c(0.5,5),cexRow=0.5,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(combined_dir_unhealthy_young_16S),rowsep=1:nrow(combined_dir_unhealthy_young_16S))

colnames(he_association_direction) <- paste0("he_",colnames(he_association_direction))
colnames(ag_association_direction) <- paste0("ag_",colnames(ag_association_direction))
colnames(nuage_association_direction) <- paste0("nuage_",colnames(nuage_association_direction))
colnames(cmd3_association_direction) <- paste0("cmd3_",colnames(cmd3_association_direction))
colnames(em_association_direction) <- paste0("em_",colnames(em_association_direction))

temp0 <- merge(he_association_direction,ag_association_direction,by="row.names",all=TRUE)[,-1]
rownames(temp0) <- merge(he_association_direction,ag_association_direction,by="row.names",all=TRUE)[,1]
temp0 <- apply(temp0,2,function(x)(ifelse(is.na(x),0,x)))

temp1 <- merge(temp0,nuage_association_direction,by="row.names",all=TRUE)[,-1]
rownames(temp1) <- merge(temp0,nuage_association_direction,by="row.names",all=TRUE)[,1]
temp1 <- apply(temp1,2,function(x)(ifelse(is.na(x),0,x)))

temp2 <- merge(temp1,cmd3_association_direction,by="row.names",all=TRUE)[,-1]
rownames(temp2) <- merge(temp1,cmd3_association_direction,by="row.names",all=TRUE)[,1]
temp2 <- apply(temp2,2,function(x)(ifelse(is.na(x),0,x)))

temp3 <- merge(temp2,em_association_direction,by="row.names",all=TRUE)[,-1]
rownames(temp3) <- merge(temp2,em_association_direction,by="row.names",all=TRUE)[,1]
temp3 <- apply(temp3,2,function(x)(ifelse(is.na(x),0,x)))

combined_association_direction <- temp3


CombinedPositiveWithUnhealthyAging <- names(which(apply(temp3[,c(2,4,6,8,10)],1,function(x)(length(x[x>1])))>=3))
CombinedNegativeWithUnhealthyAging <- names(which(apply(temp3[,c(1,3,5,7,9)],1,function(x)(length(x[x>1])))>=3))
OpposingPositiveWithUnhealthyAging <- names(which(apply(temp3[,c(1,3,5,7,9)],1,function(x)(length(x[x>1])))<=1))
OpposingNegativeWithUnhealthyAging <- names(which(apply(temp3[,c(2,4,6,8,10)],1,function(x)(length(x[x>1])))<=1))

FinalPositiveWithUnhealthyAging <- intersect(CombinedPositiveWithUnhealthyAging,OpposingPositiveWithUnhealthyAging)
FinalNegativeWithUnhealthyAging <- intersect(CombinedNegativeWithUnhealthyAging,OpposingNegativeWithUnhealthyAging)

SortedMarkersWithUnhealthyAging <- c(FinalPositiveWithUnhealthyAging,FinalNegativeWithUnhealthyAging)[order(rowSums(combined_association_direction[c(FinalPositiveWithUnhealthyAging,FinalNegativeWithUnhealthyAging),c(1,3,5,7,9)])-rowSums(combined_association_direction[c(FinalPositiveWithUnhealthyAging,FinalNegativeWithUnhealthyAging),c(2,4,6,8,10)]))]

colnames(cmd3_rlm_dir_unhealthy_aging) <- paste0("cmd3_",colnames(cmd3_rlm_dir_unhealthy_aging))
colnames(em_rlm_dir_unhealthy_aging) <- paste0("em_",colnames(em_rlm_dir_unhealthy_aging))
colnames(ag_rlm_dir_unhealthy_aging) <- paste0("ag_",colnames(ag_rlm_dir_unhealthy_aging))
colnames(nuage_rlm_dir_unhealthy_aging) <- paste0("nuage_",colnames(nuage_rlm_dir_unhealthy_aging))
colnames(he_rlm_dir_unhealthy_aging) <- paste0("he_",colnames(he_rlm_dir_unhealthy_aging))

temp0 <- merge(cmd3_rlm_dir_unhealthy_aging,em_rlm_dir_unhealthy_aging,by="row.names",all=TRUE)[,-1]
rownames(temp0) <- merge(cmd3_rlm_dir_unhealthy_aging,em_rlm_dir_unhealthy_aging,by="row.names",all=TRUE)[,1]
temp0 <- apply(temp0,2,function(x)(ifelse(is.na(x),0,x)))

temp1 <- merge(temp0,ag_rlm_dir_unhealthy_aging,by="row.names",all=TRUE)[,-1]
rownames(temp1) <- merge(temp0,ag_rlm_dir_unhealthy_aging,by="row.names",all=TRUE)[,1]
temp0 <- apply(temp1,2,function(x)(ifelse(is.na(x),0,x)))

temp2 <- merge(temp1,nuage_rlm_dir_unhealthy_aging,by="row.names",all=TRUE)[,-1]
rownames(temp2) <- merge(temp1,nuage_rlm_dir_unhealthy_aging,by="row.names",all=TRUE)[,1]
temp2 <- apply(temp2,2,function(x)(ifelse(is.na(x),0,x)))

temp3 <- merge(temp2,he_rlm_dir_unhealthy_aging,by="row.names",all=TRUE)[,-1]
rownames(temp3) <- merge(temp2,he_rlm_dir_unhealthy_aging,by="row.names",all=TRUE)[,1]
temp3 <- apply(temp3,2,function(x)(ifelse(is.na(x),0,x)))

combined_rlm_dir_unhealthy_aging <- temp3

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\taxa_groupings_new.RData")
marker_directions_unhealthy_aging <- data.frame(increased=apply(combined_rlm_dir_unhealthy_aging[SortedMarkersWithUnhealthyAging,],1,function(x)(length(x[x>=2]))),decreased=apply(combined_rlm_dir_unhealthy_aging[SortedMarkersWithUnhealthyAging,],1,function(x)(length(x[x<=-2]))))
FilteredMarkersWithUnhealthyAging <- rownames(marker_directions_unhealthy_aging)[apply(marker_directions_unhealthy_aging,1,min)<=2]

mat <- combined_rlm_dir_unhealthy_aging[FilteredMarkersWithUnhealthyAging,]

G1Species <- rownames(species_groupings[species_groupings[,7]=="G1",])
G2Species <- rownames(species_groupings[species_groupings[,7]=="G2",])

#heatmap.2(apply(mat,2,function(x)(ifelse(x<0,x-0.1,x))),density="none",trace="none",col=c("blue4","blue3","white","white","white","orange","red"),Colv=FALSE,Rowv=FALSE,margins=c(5,10),lhei=c(0.1,5),lwid=c(2,5),cexRow=0.5,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),RowSideColors=ifelse(rownames(mat) %in% G1Species,"Green2",ifelse(rownames(mat) %in% G2Species,"Orange2","Grey")))

colnames(he_association_direction_young) <- paste0("he_",colnames(he_association_direction_young))
colnames(ag_association_direction_young) <- paste0("ag_",colnames(ag_association_direction_young))
colnames(cmd3_association_direction_young) <- paste0("cmd3_",colnames(cmd3_association_direction_young))
colnames(isc_association_direction_young) <- paste0("isc_",colnames(isc_association_direction_young))

temp0 <- merge(he_association_direction_young,ag_association_direction_young,by="row.names",all=TRUE)[,-1]
rownames(temp0) <- merge(he_association_direction_young,ag_association_direction_young,by="row.names",all=TRUE)[,1]
temp0 <- apply(temp0,2,function(x)(ifelse(is.na(x),0,x)))

temp2 <- merge(temp0,cmd3_association_direction_young,by="row.names",all=TRUE)[,-1]
rownames(temp2) <- merge(temp0,cmd3_association_direction_young,by="row.names",all=TRUE)[,1]
temp2 <- apply(temp2,2,function(x)(ifelse(is.na(x),0,x)))

temp3 <- merge(temp2,isc_association_direction_young,by="row.names",all=TRUE)[,-1]
rownames(temp3) <- merge(temp2,isc_association_direction_young,by="row.names",all=TRUE)[,1]
temp3 <- apply(temp3,2,function(x)(ifelse(is.na(x),0,x)))

combined_association_direction_young <- temp3



CombinedPositiveWithUnhealthyYoung <- names(which(apply(combined_association_direction_young[,c(2,4,6,8)],1,function(x)(length(x[x>1])))>=2))
CombinedNegativeWithUnhealthyYoung <- names(which(apply(combined_association_direction_young[,c(1,3,5,7)],1,function(x)(length(x[x>1])))>=2))
OpposingPositiveWithUnhealthyYoung <- names(which(apply(combined_association_direction_young[,c(1,3,5,7)],1,function(x)(length(x[x>1])))<=1))
OpposingNegativeWithUnhealthyYoung <- names(which(apply(combined_association_direction_young[,c(2,4,6,8)],1,function(x)(length(x[x>1])))<=1))

FinalPositiveWithUnhealthyYoung <- intersect(CombinedPositiveWithUnhealthyYoung,OpposingPositiveWithUnhealthyYoung)
FinalNegativeWithUnhealthyYoung <- intersect(CombinedNegativeWithUnhealthyYoung,OpposingNegativeWithUnhealthyYoung)

SortedMarkersWithUnhealthyYoung <- c(FinalPositiveWithUnhealthyYoung,FinalNegativeWithUnhealthyYoung)[order(rowSums(combined_association_direction_young[c(FinalPositiveWithUnhealthyYoung,FinalNegativeWithUnhealthyYoung),c(1,3,5,7)])-rowSums(combined_association_direction_young[c(FinalPositiveWithUnhealthyYoung,FinalNegativeWithUnhealthyYoung),c(2,4,6,8)]))]


colnames(cmd3_rlm_dir_unhealthy_young) <- paste0("cmd3_",colnames(cmd3_rlm_dir_unhealthy_young))
colnames(isc_rlm_dir_unhealthy_young) <- paste0("isc_",colnames(isc_rlm_dir_unhealthy_young))
colnames(ag_rlm_dir_unhealthy_young) <- paste0("ag_",colnames(ag_rlm_dir_unhealthy_young))
colnames(he_rlm_dir_unhealthy_young) <- paste0("he_",colnames(he_rlm_dir_unhealthy_young))

temp0 <- merge(cmd3_rlm_dir_unhealthy_young,isc_rlm_dir_unhealthy_young,by="row.names",all=TRUE)[,-1]
rownames(temp0) <- merge(cmd3_rlm_dir_unhealthy_young,isc_rlm_dir_unhealthy_young,by="row.names",all=TRUE)[,1]
temp0 <- apply(temp0,2,function(x)(ifelse(is.na(x),0,x)))

temp1 <- merge(temp0,ag_rlm_dir_unhealthy_young,by="row.names",all=TRUE)[,-1]
rownames(temp1) <- merge(temp0,ag_rlm_dir_unhealthy_young,by="row.names",all=TRUE)[,1]
temp0 <- apply(temp1,2,function(x)(ifelse(is.na(x),0,x)))

temp3 <- merge(temp1,he_rlm_dir_unhealthy_young,by="row.names",all=TRUE)[,-1]
rownames(temp3) <- merge(temp1,he_rlm_dir_unhealthy_young,by="row.names",all=TRUE)[,1]
temp3 <- apply(temp3,2,function(x)(ifelse(is.na(x),0,x)))

combined_rlm_dir_unhealthy_young <- temp3

marker_directions_unhealthy_young <- data.frame(increased=apply(combined_rlm_dir_unhealthy_young[SortedMarkersWithUnhealthyYoung,],1,function(x)(length(x[x>=2]))),decreased=apply(combined_rlm_dir_unhealthy_young[SortedMarkersWithUnhealthyYoung,],1,function(x)(length(x[x<=-2]))))
FilteredMarkersWithUnhealthyYoung <- rownames(marker_directions_unhealthy_young)[apply(marker_directions_unhealthy_young,1,min)<=2]

mat <- combined_rlm_dir_unhealthy_young[FilteredMarkersWithUnhealthyYoung,]

heatmap.2(apply(mat,2,function(x)(ifelse(x<0,x-0.1,x))),density="none",trace="none",col=c("blue4","blue3","white","white","white","orange","red"),Colv=FALSE,Rowv=FALSE,margins=c(5,10),lhei=c(0.1,5),lwid=c(2,5),cexRow=0.5,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),RowSideColors=ifelse(rownames(mat) %in% G1Species,"Green2",ifelse(rownames(mat) %in% G2Species,"Orange2","Grey")))
