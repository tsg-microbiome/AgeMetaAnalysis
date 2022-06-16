# The sub-pipeline works on the curatedMetagenomicData3 data repository
# The objective of the sub-pipeline is to perform an integrated association analysis of microbiome features with the different microbiome features with the different 
# diseases in the He et al dataset
# The features are the 112 Highly Detected Species (across the seven data repositories) (listed in Figure 2) and scanned here in from the workspace 
# taxa_groupings_new.RData; the four microbiome uniqueness measures, shannon diversity, the range scaled abundances of the Kendall-Positive (Disease-Enriched) G2, 
# Kendall-Negative (Health-Enriched) G1, Kendall-Not-Associated species groups; the Multiple-Disease-Enriched and Multiple-Disease-Depleted taxa identified in Ghosh
# et al eLife 2020.
# The disease associations are computed separately for the gut microbiomes from elderly (age >= 60y) and young (age < 60y) subjects.
# The outputs of this sub-pipeline when integrated with those from the other data repositories provide the ranked order of microbiome features associated with
# the unhealthy phenotype in the elderly (Figure 5) and the young (Supplementary Figure S12).

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

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")

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

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\taxa_groupings_new.RData")

colnames(cmd3_select_disease_final_species_clr)[grep("Hungatella_hathewayi",colnames(cmd3_select_disease_final_species_clr))] <- "Clostridium_hathewayi"
colnames(cmd3_select_disease_final_species_clr)[grep("Tyzzerella_nexilis",colnames(cmd3_select_disease_final_species_clr))]  <- "Clostridium_nexile"
colnames(cmd3_select_disease_final_species_clr)[grep("Streptococcus_anginosus_group",colnames(cmd3_select_disease_final_species_clr))]  <- "Streptococcus_anginosus"
colnames(cmd3_select_disease_final_species_clr)[grep("Enterobacter_cloacae_group",colnames(cmd3_select_disease_final_species_clr))]  <- "Enterobacter_cloacae"
colnames(cmd3_select_disease_final_species_clr)[grep("Erysipelatoclostridium_ramosum",colnames(cmd3_select_disease_final_species_clr))]  <- "Clostridium_ramosum"

G_Markers <- intersect(c("Bacteroides_fragilis","Actinomyces_odontolyticus","Bifidobacterium_dentium","Clostridium_clostridioforme","Clostridium_nexile","Clostridium_ramosum","Dialister_invisus","Eggerthella_lenta","Escherichia_coli","Fusobacterium_nucleatum","Granulicatella_adiacens","Lactococcus_lactis","Rothia_mucilaginosa","Ruminococcus_gnavus","Streptococcus_anginosus","Streptococcus_infantis","Streptococcus_salivarius","Streptococcus_sanguinis","Streptococcus_vestibularis","Subdoligranulum_variabile","Veillonella_atypica","Clostridium_asparagiforme","Clostridium_bolteae","Clostridium_citroniae","Clostridium_hathewayi","Clostridium_symbiosum","Streptococcus_parasanguinis","Solobacterium_moorei","Ruminococcus_torques","Streptococcus_mitis","Klebsiella_pneumoniae","Streptococcus_australis","Streptococcus_gordonii","Enterobacter_cloacae"),colnames(cmd3_select_disease_final_species_clr))

L_Markers <- intersect(c("Eubacterium_hallii","Dorea_longicatena","Coprococcus_comes","Coprococcus_catus","Butyrivibrio_crossotus","Bacteroides_uniformis","Alistipes_shahii","Alistipes_indistinctus","Roseburia_hominis","Pseudoflavonifractor_capillosus","Eubacterium_siraeum","Eubacterium_rectale","Coprobacter_fastidiosus","Bifidobacterium_longum","Bifidobacterium_animalis","Barnesiella_intestinihominis","Bacteroides_xylanisolvens","Alistipes_senegalensis","Alistipes_putredinis","Alistipes_onderdonkii","Akkermansia_muciniphila","Faecalibacterium_prausnitzii","Roseburia_inulinivorans","Roseburia_intestinalis"),colnames(cmd3_select_disease_final_species_clr))

df_cmd3_disease_investigation <- as.data.frame(df_cmd3_disease_diversity_uniqueness[rownames(cmd3,c("study_name","age","study_condition","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")])

df_cmd3_disease_investigation$G1 <- rowMeans(apply(apply(cmd3_select_disease_final_species_clr[rownames(df_cmd3_disease_investigation),intersect(colnames(cmd3_select_disease_final_species_clr),rownames(species_groupings[species_groupings$BroadGroup=="G1",]))],2,range_scale),2,function(x)(ifelse(is.nan(x),0,x))))
df_cmd3_disease_investigation$G2 <- rowMeans(apply(apply(cmd3_select_disease_final_species_clr[rownames(df_cmd3_disease_investigation),intersect(colnames(cmd3_select_disease_final_species_clr),rownames(species_groupings[species_groupings$BroadGroup=="G2",]))],2,range_scale),2,function(x)(ifelse(is.nan(x),0,x))))
df_cmd3_disease_investigation$Other <- rowMeans(apply(apply(cmd3_select_disease_final_species_clr[rownames(df_cmd3_disease_investigation),intersect(colnames(cmd3_select_disease_final_species_clr),rownames(species_groupings[(species_groupings$BroadGroup!="G1")&(species_groupings$BroadGroup!="G2"),]))],2,range_scale),2,function(x)(ifelse(is.nan(x),0,x))))

df_cmd3_disease_investigation$L_Markers <- rowMeans(apply(apply(cmd3_select_disease_final_species_clr[rownames(df_cmd3_disease_investigation),L_Markers],2,range_scale),2,function(x)(ifelse(is.nan(x),0,x))))
df_cmd3_disease_investigation$G_Markers <- rowMeans(apply(apply(cmd3_select_disease_final_species_clr[rownames(df_cmd3_disease_investigation),G_Markers],2,range_scale),2,function(x)(ifelse(is.nan(x),0,x))))

G1Species <- intersect(rownames(species_groupings[species_groupings$BroadGroup=="G1",]),colnames(cmd3_select_disease_final_species_clr))
G2Species <- intersect(rownames(species_groupings[species_groupings$BroadGroup=="G2",]),colnames(cmd3_select_disease_final_species_clr))
OtherSpecies <- intersect(setdiff(rownames(species_groupings),c(G1Species,G2Species)),colnames(cmd3_select_disease_final_species_clr))

df_cmd3_disease_investigation <- as.data.frame(cbind(df_cmd3_disease_investigation,cmd3_select_disease_final_species_clr[rownames(df_cmd3_disease_investigation),c(G1Species,G2Species,OtherSpecies)]))

cmd3_disease_names <- setdiff(unique(df_cmd3_disease_investigation$study_condition),"control")
cmd3_elderly_ds <- rownames(df_cmd3_disease_investigation[df_cmd3_disease_investigation$age>=60,])
cmd3_young_ds <- rownames(df_cmd3_disease_investigation[df_cmd3_disease_investigation$age<60,])

study_dummy <- NULL
study_dummy[1:54] <- NA
names(study_dummy) <- unique(df_cmd3_disease_investigation$study_name)
for(i in 1:54)
{
	study_dummy[i] <- as.integer(rnorm(1,0,40))
}

df_cmd3_disease_investigation$study_dummy <- as.numeric(study_dummy[df_cmd3_disease_investigation$study_name])

cmd3_disease_unique_study_names <- setdiff(unique(cmd3_select_disease_final_metadata$study_name),c("HMP_2019_t2d","HMP_2019_ibdmdb"))

df_cmd3_disease_numbers_elderly <- as.data.frame(matrix(0,length(cmd3_disease_names),3))
rownames(df_cmd3_disease_numbers_elderly) <- cmd3_disease_names
colnames(df_cmd3_disease_numbers_elderly) <- c("Number_Patients","Number_Controls","Disease_Studies")
for(i in 1:length(cmd3_disease_names))
{
	disease_name <- cmd3_disease_names[i]
	
	disease_study_list <- intersect(unique(df_cmd3_disease_investigation[df_cmd3_disease_investigation$study_condition == disease_name,"study_name"]),unique(df_cmd3_disease_investigation[df_cmd3_disease_investigation$study_condition == "control","study_name"]))
	assign(paste0("disease_list_",disease_name),disease_study_list)
	df_cmd3_disease_numbers_elderly[disease_name,3] <- paste0(disease_study_list,collapse=", ")
	
	patient_list <- intersect(rownames(df_cmd3_disease_investigation[(df_cmd3_disease_investigation$study_name %in% disease_study_list)&(df_cmd3_disease_investigation$study_condition == disease_name),]),cmd3_elderly_ds)
	assign(paste0("patient_list_",disease_name),patient_list)
	df_cmd3_disease_numbers_elderly[disease_name,1] <- length(intersect(cmd3_elderly_ds,patient_list))
	
	control_list <- intersect(rownames(df_cmd3_disease_investigation[(df_cmd3_disease_investigation$study_name %in% disease_study_list)&(df_cmd3_disease_investigation$study_condition == "control"),]),cmd3_elderly_ds)
	assign(paste0("control_list_",disease_name),control_list)
	df_cmd3_disease_numbers_elderly[disease_name,2] <- length(intersect(cmd3_elderly_ds,control_list))

}

df_cmd3_disease_numbers_young <- as.data.frame(matrix(0,length(cmd3_disease_names),3))
rownames(df_cmd3_disease_numbers_young) <- cmd3_disease_names
colnames(df_cmd3_disease_numbers_young) <- c("Number_Patients","Number_Controls","Disease_Studies")
for(i in 1:length(cmd3_disease_names))
{
	disease_name <- cmd3_disease_names[i]
	
	disease_study_list <- intersect(intersect(unique(df_cmd3_disease_investigation[df_cmd3_disease_investigation$study_condition == disease_name,"study_name"]),unique(df_cmd3_disease_investigation[df_cmd3_disease_investigation$study_condition == "control","study_name"])),cmd3_disease_unique_study_names)
	assign(paste0("disease_list_",disease_name),disease_study_list)
	df_cmd3_disease_numbers_young[disease_name,3] <- paste0(disease_study_list,collapse=", ")
	
	patient_list_young <- intersect(rownames(df_cmd3_disease_investigation[(df_cmd3_disease_investigation$study_name %in% disease_study_list)&(df_cmd3_disease_investigation$study_condition == disease_name),]),cmd3_young_ds)
	assign(paste0("patient_list_young_",disease_name),patient_list_young)
	df_cmd3_disease_numbers_young[disease_name,1] <- length(intersect(cmd3_young_ds,patient_list_young))
	
	control_list_young <- intersect(rownames(df_cmd3_disease_investigation[(df_cmd3_disease_investigation$study_name %in% disease_study_list)&(df_cmd3_disease_investigation$study_condition == "control"),]),cmd3_young_ds)
	assign(paste0("control_list_young_",disease_name),control_list_young)
	df_cmd3_disease_numbers_young[disease_name,2] <- length(intersect(cmd3_young_ds,control_list_young))

}

FinalIndexList <- rownames(df_cmd3_disease_numbers_elderly[(df_cmd3_disease_numbers_elderly[,1]>=10)&(df_cmd3_disease_numbers_elderly[,2]>=10),])

FinalFeatureList <- setdiff(c("species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers","L_Markers","G1","G2",intersect(colnames(df_cmd3_disease_investigation),c(G1Species,G2Species,OtherSpecies))),c("L_Markers_combined","G_Markers_combined"))

cmd3_rlm_est_unhealthy_aging <- as.data.frame(matrix(0,length(FinalFeatureList),length(FinalIndexList)))
rownames(cmd3_rlm_est_unhealthy_aging) <- FinalFeatureList
colnames(cmd3_rlm_est_unhealthy_aging) <- FinalIndexList
cmd3_rlm_p_val_unhealthy_aging <- as.data.frame(matrix(0,length(FinalFeatureList),length(FinalIndexList)))
rownames(cmd3_rlm_p_val_unhealthy_aging) <- FinalFeatureList
colnames(cmd3_rlm_p_val_unhealthy_aging) <- FinalIndexList
for(i in 1:length(FinalFeatureList))
{
	feature_name <- FinalFeatureList[i]
	print(feature_name)
	for(j in 1:length(FinalIndexList))
	{
		disease_name <- FinalIndexList[j]
		study_list <- get(paste0("disease_list_",disease_name))
		patient_list <- get(paste0("patient_list_",disease_name))
		control_list <- get(paste0("control_list_",disease_name))
		vec <- df_cmd3_disease_investigation[c(patient_list,control_list),feature_name]
		if(length(vec[vec>0])>0)
		{
			temp_rlm <- rlm(df_cmd3_disease_investigation[c(patient_list,control_list),feature_name]~as.factor(ifelse(df_cmd3_disease_investigation[c(patient_list,control_list),"study_condition"]=="control",0,1)))
			summary_temp_rlm <- summary(temp_rlm)
			cmd3_rlm_est_unhealthy_aging[feature_name,disease_name] <- summary_temp_rlm$coefficients[2,3]
			cmd3_rlm_p_val_unhealthy_aging[feature_name,disease_name] <- wilcox.test(df_cmd3_disease_investigation[patient_list,feature_name],df_cmd3_disease_investigation[control_list,feature_name])$p.value
		}
	}
}

cmd3_rlm_q_val_unhealthy_aging <- apply(cmd3_rlm_p_val_unhealthy_aging,2,function(x)(p.adjust(x,method="fdr")))

cmd3_rlm_dir_unhealthy_aging <- as.data.frame(matrix(0,length(FinalFeatureList),length(FinalIndexList)))
rownames(cmd3_rlm_dir_unhealthy_aging) <- FinalFeatureList
colnames(cmd3_rlm_dir_unhealthy_aging) <- FinalIndexList
for(i in 1:length(FinalFeatureList))
{
	for(j in 1:length(FinalIndexList))
	{
		cmd3_rlm_dir_unhealthy_aging[i,j] <- ifelse(cmd3_rlm_q_val_unhealthy_aging[i,j] <= 0.1,3*sign(cmd3_rlm_est_unhealthy_aging[i,j]),ifelse(cmd3_rlm_p_val_unhealthy_aging[i,j] <= 0.05,2*sign(cmd3_rlm_est_unhealthy_aging[i,j]),sign(cmd3_rlm_est_unhealthy_aging[i,j])))
	}
}

cmd3_association_direction <- data.frame(decreased=apply(cmd3_rlm_dir_unhealthy_aging,1,function(x)(length(x[x<= -2]))),increased=apply(cmd3_rlm_dir_unhealthy_aging,1,function(x)(length(x[x>=2]))))

cmd3_rlm_est_unhealthy_young <- as.data.frame(matrix(0,length(FinalFeatureList),length(FinalIndexList)))
rownames(cmd3_rlm_est_unhealthy_young) <- FinalFeatureList
colnames(cmd3_rlm_est_unhealthy_young) <- FinalIndexList
cmd3_rlm_p_val_unhealthy_young <- as.data.frame(matrix(1,length(FinalFeatureList),length(FinalIndexList)))
rownames(cmd3_rlm_p_val_unhealthy_young) <- FinalFeatureList
colnames(cmd3_rlm_p_val_unhealthy_young) <- FinalIndexList
for(i in 1:length(FinalFeatureList))
{
	feature_name <- FinalFeatureList[i]
	print(feature_name)
	for(j in 1:length(FinalIndexList))
	{
		disease_name <- FinalIndexList[j]
		study_list <- get(paste0("disease_list_",disease_name))
		patient_list <- get(paste0("patient_list_young_",disease_name))
		control_list <- get(paste0("control_list_young_",disease_name))
		vec <- df_cmd3_disease_investigation[c(patient_list,control_list),feature_name]
		if(length(vec[vec>0])>0)
		{
			temp_rlm <- rlm(df_cmd3_disease_investigation[c(patient_list,control_list),feature_name]~as.factor(ifelse(df_cmd3_disease_investigation[c(patient_list,control_list),"study_condition"]=="control",0,1)))
			summary_temp_rlm <- summary(temp_rlm)
			cmd3_rlm_est_unhealthy_young[feature_name,disease_name] <- summary_temp_rlm$coefficients[2,3]
			cmd3_rlm_p_val_unhealthy_young[feature_name,disease_name] <- wilcox.test(df_cmd3_disease_investigation[patient_list,feature_name],df_cmd3_disease_investigation[control_list,feature_name])$p.value
		}
	}
}

cmd3_rlm_q_val_unhealthy_young <- apply(cmd3_rlm_p_val_unhealthy_young,2,function(x)(p.adjust(x,method="fdr")))

cmd3_rlm_dir_unhealthy_young <- as.data.frame(matrix(0,length(FinalFeatureList),length(FinalIndexList)))
rownames(cmd3_rlm_dir_unhealthy_young) <- FinalFeatureList
colnames(cmd3_rlm_dir_unhealthy_young) <- FinalIndexList
for(i in 1:length(FinalFeatureList))
{
	for(j in 1:length(FinalIndexList))
	{
		cmd3_rlm_dir_unhealthy_young[i,j] <- ifelse(cmd3_rlm_q_val_unhealthy_young[i,j] <= 0.1,3*sign(cmd3_rlm_est_unhealthy_young[i,j]),ifelse(cmd3_rlm_p_val_unhealthy_young[i,j] <= 0.05,2*sign(cmd3_rlm_est_unhealthy_young[i,j]),sign(cmd3_rlm_est_unhealthy_young[i,j])))
	}
}

cmd3_association_direction_young <- data.frame(decreased=apply(cmd3_rlm_dir_unhealthy_young,1,function(x)(length(x[x<= -2]))),increased=apply(cmd3_rlm_dir_unhealthy_young,1,function(x)(length(x[x>=2]))))

save(cmd3_rlm_est_unhealthy_aging,cmd3_rlm_p_val_unhealthy_aging,cmd3_rlm_q_val_unhealthy_aging,cmd3_rlm_dir_unhealthy_aging,cmd3_association_direction,cmd3_rlm_est_unhealthy_young,cmd3_rlm_p_val_unhealthy_young,cmd3_rlm_q_val_unhealthy_young,cmd3_rlm_dir_unhealthy_young,cmd3_association_direction_young,df_cmd3_diversity_uniqueness,df_cmd3_disease_investigation,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_stage_2c_results.RData")

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")

rm(list=ls())
