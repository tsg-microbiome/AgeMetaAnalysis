library(robumeta)
library(metafor)
library(dplyr)
library(effsize)
library(MASS)
library(sfsmisc)
library(compositions)
library(igraph)
library(metap)


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
		
		#print(group)
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
	res <- rma(yi, vi, data=temp_meta,control=list(maxiter=1000,stepadj=0.5))
	res$ids <- rownames(temp_meta)
	res$slabs <- rownames(temp_meta)
	return_list <- list("df_studies"=temp_meta,"model"=res)
	return(return_list)
}

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ag_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\he_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\odamaki_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_stage1_pathway_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\nuage_stage1_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\logmpie_stage1_results.RData")

print("Data Loaded")

cmd3_pathway_p_val_age_beta <- cmd3_p_val_age_beta[,c(11:15)]
cmd3_p_val_age_beta <- cmd3_p_val_age_beta[,c(1:10)]

cmd3_pathway_est_age_beta <- cmd3_est_age_beta[,c(11:15)]
cmd3_est_age_beta <- cmd3_est_age_beta[,c(1:10)]

sorted_study_rows <- c("HMP_2019_ibdmdb","SankaranarayananK_2015","CosteaPI_2017","AsnicarF_2021","HansenLBS_2018","NielsenHB_2014","SchirmerM_2016","WirbelJ_2018","ZellerG_2014","KeohaneDM_2020","QinN_2014","YeZ_2018","QinJ_2012","YachidaS_2019","DhakanDB_2019","GuptaA_2019","BritoIL_2016","PehrssonE_2016","LokmerA_2019","PasolliE_2019","RubelMA_2020","RampelliS_2015","AG","NUAGE","ISC","HE","Odamaki","LOGMPIE")

sorted_study_rows_pathways <- c("HMP_2019_ibdmdb","SankaranarayananK_2015","CosteaPI_2017","AsnicarF_2021","HansenLBS_2018","NielsenHB_2014","SchirmerM_2016","WirbelJ_2018","ZellerG_2014","KeohaneDM_2020","ISC","QinN_2014","YeZ_2018","QinJ_2012","YachidaS_2019","DhakanDB_2019","GuptaA_2019","BritoIL_2016","PehrssonE_2016","LokmerA_2019","PasolliE_2019","RubelMA_2020","RampelliS_2015")

combined_p_val_age_beta <- as.data.frame(rbind(cmd3_p_val_age_beta,ag_p_val_age_beta,he_p_val_age_beta,isc_p_val_age_beta,odamaki_p_val_age_beta,nuage_p_val_age_beta,logmpie_p_val_age_beta))

combined_pathway_p_val_age_beta <- as.data.frame(rbind(cmd3_pathway_p_val_age_beta,isc_pathway_p_val_age_beta))

combined_p_val_age_beta <- combined_p_val_age_beta[sorted_study_rows,c(3:10)]
#mat <- apply(combined_p_val_age_beta,1,function(x)(ifelse(x<=0.05,1,0)))
#heatmap.2(apply(combined_p_val_age_beta,1,function(x)(ifelse(x<=0.05,1,0))),density="none",trace="none",col=c("white","cyan"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),key=FALSE,margins=c(25,5),lhei=c(0.1,5),lwid=c(1,5),cellnote=apply(t(combined_p_val_age_beta),2,function(x)(ifelse(x<=0.001,"***",ifelse(x<=0.01,"**",ifelse(x<=0.05,"*",""))))),notecol="black",notecex=1.5)

combined_pathway_p_val_age_beta <- combined_pathway_p_val_age_beta[sorted_study_rows_pathways,c(2:5)]
#mat <- apply(combined_pathway_p_val_age_beta,1,function(x)(ifelse(x<=0.05,1,0)))
#heatmap.2(apply(combined_pathway_p_val_age_beta,1,function(x)(ifelse(x<=0.05,1,0))),density="none",trace="none",col=c("white","cyan"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),key=FALSE,margins=c(25,5),lhei=c(0.1,5),lwid=c(1,5),cellnote=apply(t(combined_pathway_p_val_age_beta),2,function(x)(ifelse(x<=0.001,"***",ifelse(x<=0.01,"**",ifelse(x<=0.05,"*",""))))),notecol="black",notecex=1.5)

cmd3_studies <- c("HMP_2019_ibdmdb","SankaranarayananK_2015","CosteaPI_2017","AsnicarF_2021","HansenLBS_2018","NielsenHB_2014","SchirmerM_2016","WirbelJ_2018","ZellerG_2014","KeohaneDM_2020","QinN_2014","YeZ_2018","QinJ_2012","YachidaS_2019","DhakanDB_2019","GuptaA_2019","BritoIL_2016","PehrssonE_2016","LokmerA_2019","PasolliE_2019","RubelMA_2020","RampelliS_2015")

combined_study_details_controls <- as.data.frame(matrix(NA,length(sorted_study_rows),8))
rownames(combined_study_details_controls) <- sorted_study_rows
colnames(combined_study_details_controls) <- c("data_repository","data_type","nationality","region","minimum_age","maximum_age","total","percentage_elderly")

for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	print(study_name)
	if(study_name %in% cmd3_studies)
	{
		combined_study_details_controls[study_name,"data_repository"] <- "CMD3"
		combined_study_details_controls[study_name,"data_type"] <- "Shotgun"
		combined_study_details_controls[study_name,"nationality"] <- paste0(unique(df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name %in% study_name,"country"]),collapse=",")
		combined_study_details_controls[study_name,"minimum_age"] <- min(df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name %in% study_name,"age"])
		combined_study_details_controls[study_name,"maximum_age"] <- max(df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name %in% study_name,"age"])
		combined_study_details_controls[study_name,"total"] <- nrow(df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name %in% study_name,])
		combined_study_details_controls[study_name,"percentage_elderly"] <- 100*(nrow(df_cmd3_controls_diversity_uniqueness[(df_cmd3_controls_diversity_uniqueness$study_name %in% study_name)&(df_cmd3_controls_diversity_uniqueness$age >= 60),]))/(nrow(df_cmd3_controls_diversity_uniqueness[df_cmd3_controls_diversity_uniqueness$study_name %in% study_name,]))
	}
	if(study_name == "AG")
	{
		combined_study_details_controls[study_name,"data_repository"] <- "AGP"
		combined_study_details_controls[study_name,"data_type"] <- "16S"
		combined_study_details_controls[study_name,"nationality"] <- paste0(c("USA","GBR"),collapse=",")
		combined_study_details_controls[study_name,"minimum_age"] <- min(df_ag_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"maximum_age"] <- max(df_ag_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"total"] <- nrow(df_ag_controls_diversity_uniqueness)
		combined_study_details_controls[study_name,"percentage_elderly"] <- 100*(nrow(df_ag_controls_diversity_uniqueness[(df_ag_controls_diversity_uniqueness$age >= 60),]))/(nrow(df_ag_controls_diversity_uniqueness))
	}
	if(study_name == "HE")
	{
		combined_study_details_controls[study_name,"data_repository"] <- "HE"
		combined_study_details_controls[study_name,"data_type"] <- "16S"
		combined_study_details_controls[study_name,"nationality"] <- "CHN"
		combined_study_details_controls[study_name,"minimum_age"] <- min(df_he_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"maximum_age"] <- max(df_he_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"total"] <- nrow(df_he_controls_diversity_uniqueness)
		combined_study_details_controls[study_name,"percentage_elderly"] <- 100*(nrow(df_he_controls_diversity_uniqueness[(df_he_controls_diversity_uniqueness$age >= 60),]))/(nrow(df_he_controls_diversity_uniqueness))
	}
	if(study_name == "Odamaki")
	{
		combined_study_details_controls[study_name,"data_repository"] <- "Odamaki"
		combined_study_details_controls[study_name,"data_type"] <- "16S"
		combined_study_details_controls[study_name,"nationality"] <- "JPN"
		combined_study_details_controls[study_name,"minimum_age"] <- min(df_odamaki_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"maximum_age"] <- max(df_odamaki_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"total"] <- nrow(df_odamaki_controls_diversity_uniqueness)
		combined_study_details_controls[study_name,"percentage_elderly"] <- 100*(nrow(df_odamaki_controls_diversity_uniqueness[(df_odamaki_controls_diversity_uniqueness$age >= 60),]))/(nrow(df_odamaki_controls_diversity_uniqueness))
	}
	if(study_name == "NUAGE")
	{
		combined_study_details_controls[study_name,"data_repository"] <- "NUAGE"
		combined_study_details_controls[study_name,"data_type"] <- "16S"
		combined_study_details_controls[study_name,"nationality"] <- "EU"
		combined_study_details_controls[study_name,"minimum_age"] <- min(df_nuage_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"maximum_age"] <- max(df_nuage_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"total"] <- nrow(df_nuage_controls_diversity_uniqueness)
		combined_study_details_controls[study_name,"percentage_elderly"] <- 100*(nrow(df_nuage_controls_diversity_uniqueness[(df_nuage_controls_diversity_uniqueness$age >= 60),]))/(nrow(df_nuage_controls_diversity_uniqueness))
	}
	if(study_name == "ISC")
	{
		combined_study_details_controls[study_name,"data_repository"] <- "ISC"
		combined_study_details_controls[study_name,"data_type"] <- "Shotgun"
		combined_study_details_controls[study_name,"nationality"] <- "IRE"
		combined_study_details_controls[study_name,"minimum_age"] <- min(df_isc_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"maximum_age"] <- max(df_isc_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"total"] <- nrow(df_isc_controls_diversity_uniqueness)
		combined_study_details_controls[study_name,"percentage_elderly"] <- 100*(nrow(df_isc_controls_diversity_uniqueness[(df_isc_controls_diversity_uniqueness$age >= 60),]))/(nrow(df_isc_controls_diversity_uniqueness))
	}
	if(study_name == "LOGMPIE")
	{
		combined_study_details_controls[study_name,"data_repository"] <- "LOGMPIE"
		combined_study_details_controls[study_name,"data_type"] <- "16S"
		combined_study_details_controls[study_name,"nationality"] <- "IND"
		combined_study_details_controls[study_name,"minimum_age"] <- min(df_logmpie_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"maximum_age"] <- max(df_logmpie_controls_diversity_uniqueness[,"age"])
		combined_study_details_controls[study_name,"total"] <- nrow(df_logmpie_controls_diversity_uniqueness)
		combined_study_details_controls[study_name,"percentage_elderly"] <- 100*(nrow(df_logmpie_controls_diversity_uniqueness[(df_logmpie_controls_diversity_uniqueness$age >= 60),]))/(nrow(df_logmpie_controls_diversity_uniqueness))
	}
	
	
}

combined_study_details_all_samples <- as.data.frame(matrix(NA,length(sorted_study_rows),8))
rownames(combined_study_details_all_samples) <- sorted_study_rows
colnames(combined_study_details_all_samples) <- c("data_repository","data_type","nationality","region","minimum_age","maximum_age","total","total_elderly")

for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	print(paste0(study_name,"1"))
	if(study_name %in% cmd3_studies)
	{
		combined_study_details_all_samples[study_name,"data_repository"] <- "CMD3"
		combined_study_details_all_samples[study_name,"data_type"] <- "Shotgun"
		combined_study_details_all_samples[study_name,"nationality"] <- paste0(unique(df_cmd3_diversity_uniqueness[df_cmd3_diversity_uniqueness$study_name %in% study_name,"country"]),collapse=",")
		combined_study_details_all_samples[study_name,"minimum_age"] <- min(df_cmd3_diversity_uniqueness[df_cmd3_diversity_uniqueness$study_name %in% study_name,"age"])
		combined_study_details_all_samples[study_name,"maximum_age"] <- max(df_cmd3_diversity_uniqueness[df_cmd3_diversity_uniqueness$study_name %in% study_name,"age"])
		combined_study_details_all_samples[study_name,"total"] <- nrow(df_cmd3_diversity_uniqueness[df_cmd3_diversity_uniqueness$study_name %in% study_name,])
		#combined_study_details_all_samples[study_name,"percentage_elderly"] <- 100*(nrow(df_cmd3_diversity_uniqueness[(df_cmd3_diversity_uniqueness$study_name %in% study_name)&(df_cmd3_diversity_uniqueness$age >= 60),]))/(nrow(df_cmd3_diversity_uniqueness[df_cmd3_diversity_uniqueness$study_name %in% study_name,]))
		combined_study_details_all_samples[study_name,"total_elderly"] <- nrow(df_cmd3_diversity_uniqueness[(df_cmd3_diversity_uniqueness$study_name %in% study_name)&(df_cmd3_diversity_uniqueness$age >= 60),])
	}
	if(study_name == "AG")
	{
		combined_study_details_all_samples[study_name,"data_repository"] <- "AGP"
		combined_study_details_all_samples[study_name,"data_type"] <- "16S"
		combined_study_details_all_samples[study_name,"nationality"] <- paste0(c("USA","GBR"),collapse=",")
		combined_study_details_all_samples[study_name,"minimum_age"] <- min(df_ag_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"maximum_age"] <- max(df_ag_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"total"] <- nrow(df_ag_diversity_uniqueness)
		#combined_study_details_all_samples[study_name,"percentage_elderly"] <- 100*(nrow(df_ag_diversity_uniqueness[(df_ag_diversity_uniqueness$age >= 60),]))/(nrow(df_ag_diversity_uniqueness))
		combined_study_details_all_samples[study_name,"total_elderly"] <- nrow(df_ag_diversity_uniqueness[(df_ag_diversity_uniqueness$age >= 60),])
	}
	if(study_name == "HE")
	{
		combined_study_details_all_samples[study_name,"data_repository"] <- "HE"
		combined_study_details_all_samples[study_name,"data_type"] <- "16S"
		combined_study_details_all_samples[study_name,"nationality"] <- "CHN"
		combined_study_details_all_samples[study_name,"minimum_age"] <- min(df_he_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"maximum_age"] <- max(df_he_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"total"] <- nrow(df_he_diversity_uniqueness)
		#combined_study_details_all_samples[study_name,"percentage_elderly"] <- 100*(nrow(df_he_diversity_uniqueness[(df_he_diversity_uniqueness$age >= 60),]))/(nrow(df_he_diversity_uniqueness))
		combined_study_details_all_samples[study_name,"total_elderly"] <- nrow(df_he_diversity_uniqueness[(df_he_diversity_uniqueness$age >= 60),])
	}
	if(study_name == "Odamaki")
	{
		combined_study_details_all_samples[study_name,"data_repository"] <- "Odamaki"
		combined_study_details_all_samples[study_name,"data_type"] <- "16S"
		combined_study_details_all_samples[study_name,"nationality"] <- "JPN"
		combined_study_details_all_samples[study_name,"minimum_age"] <- min(df_odamaki_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"maximum_age"] <- max(df_odamaki_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"total"] <- nrow(df_odamaki_diversity_uniqueness)
		#combined_study_details_all_samples[study_name,"percentage_elderly"] <- 100*(nrow(df_odamaki_diversity_uniqueness[(df_odamaki_diversity_uniqueness$age >= 60),]))/(nrow(df_odamaki_diversity_uniqueness))
		combined_study_details_all_samples[study_name,"total_elderly"] <- nrow(df_odamaki_diversity_uniqueness[(df_odamaki_diversity_uniqueness$age >= 60),])
		
	}
	if(study_name == "NUAGE")
	{
		combined_study_details_all_samples[study_name,"data_repository"] <- "NUAGE"
		combined_study_details_all_samples[study_name,"data_type"] <- "16S"
		combined_study_details_all_samples[study_name,"nationality"] <- "EU"
		combined_study_details_all_samples[study_name,"minimum_age"] <- min(df_nuage_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"maximum_age"] <- max(df_nuage_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"total"] <- nrow(df_nuage_diversity_uniqueness)
		#combined_study_details_all_samples[study_name,"percentage_elderly"] <- 100*(nrow(df_nuage_diversity_uniqueness[(df_nuage_diversity_uniqueness$age >= 60),]))/(nrow(df_nuage_diversity_uniqueness))
		combined_study_details_all_samples[study_name,"total_elderly"] <- nrow(df_nuage_diversity_uniqueness[(df_nuage_diversity_uniqueness$age >= 60),])
		
	}
	if(study_name == "ISC")
	{
		combined_study_details_all_samples[study_name,"data_repository"] <- "ISC"
		combined_study_details_all_samples[study_name,"data_type"] <- "Shotgun"
		combined_study_details_all_samples[study_name,"nationality"] <- "IRE"
		combined_study_details_all_samples[study_name,"minimum_age"] <- min(df_isc_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"maximum_age"] <- max(df_isc_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"total"] <- nrow(df_isc_diversity_uniqueness)
		#combined_study_details_all_samples[study_name,"percentage_elderly"] <- 100*(nrow(df_isc_diversity_uniqueness[(df_isc_diversity_uniqueness$age >= 60),]))/(nrow(df_isc_diversity_uniqueness))
		combined_study_details_all_samples[study_name,"total_elderly"] <- nrow(df_isc_diversity_uniqueness[(df_isc_diversity_uniqueness$age >= 60),])
	}
	if(study_name == "LOGMPIE")
	{
		combined_study_details_all_samples[study_name,"data_repository"] <- "LOGMPIE"
		combined_study_details_all_samples[study_name,"data_type"] <- "16S"
		combined_study_details_all_samples[study_name,"nationality"] <- "IND"
		combined_study_details_all_samples[study_name,"minimum_age"] <- min(df_logmpie_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"maximum_age"] <- max(df_logmpie_diversity_uniqueness[,"age"])
		combined_study_details_all_samples[study_name,"total"] <- nrow(df_logmpie_diversity_uniqueness)
		combined_study_details_all_samples[study_name,"percentage_elderly"] <- 100*(nrow(df_logmpie_diversity_uniqueness[(df_logmpie_diversity_uniqueness$age >= 60),]))/(nrow(df_logmpie_diversity_uniqueness))
		combined_study_details_all_samples[study_name,"total_elderly"] <- nrow(df_logmpie_diversity_uniqueness[(df_logmpie_diversity_uniqueness$age >= 60),])
	}
	
	
}
combined_study_details_all_samples[1:11,"region"] <- "EU_NA"
combined_study_details_all_samples[12:14,"region"] <- "EstAsia"
combined_study_details_all_samples[15:16,"region"] <- "SthAsia"
combined_study_details_all_samples[17:18,"region"] <- "Pc_SA"
combined_study_details_all_samples[19:23,"region"] <- "Afr"
combined_study_details_all_samples["AG","region"] <- "EU_NA"
combined_study_details_all_samples["NUAGE","region"] <- "EU_NA"
combined_study_details_all_samples["ISC","region"] <- "EU_NA"
combined_study_details_all_samples["HE","region"] <- "EstAsia"
combined_study_details_all_samples["Odamaki","region"] <- "EstAsia"
combined_study_details_all_samples["LOGMPIE","region"] <- "SthAsia"
combined_study_details_controls[1:11,"region"] <- "EU_NA"
combined_study_details_controls[12:14,"region"] <- "EstAsia"
combined_study_details_controls[15:16,"region"] <- "SthAsia"
combined_study_details_controls[17:18,"region"] <- "Pc_SA"
combined_study_details_controls[19:23,"region"] <- "Afr"
combined_study_details_controls["AG","region"] <- "EU_NA"
combined_study_details_controls["NUAGE","region"] <- "EU_NA"
combined_study_details_controls["ISC","region"] <- "EU_NA"
combined_study_details_controls["HE","region"] <- "EstAsia"
combined_study_details_controls["Odamaki","region"] <- "EstAsia"
combined_study_details_controls["LOGMPIE","region"] <- "SthAsia"

print("Study Profiles created")

cmd3_rlm_est_pathway_age_sum_stat <- cmd3_rlm_est_age_sum_stat[,c(11:15)]
cmd3_rlm_est_age_sum_stat <- cmd3_rlm_est_age_sum_stat[,c(1:10)]

cmd3_rlm_p_val_pathway_age_sum_stat <- cmd3_rlm_p_val_age_sum_stat[,c(11:15)]
cmd3_rlm_p_val_age_sum_stat <- cmd3_rlm_p_val_age_sum_stat[,c(1:10)]

cmd3_rlm_q_val_pathway_age_sum_stat <- cmd3_rlm_q_val_age_sum_stat[,c(11:15)]
cmd3_rlm_q_val_age_sum_stat <- cmd3_rlm_q_val_age_sum_stat[,c(1:10)]

cmd3_rlm_dir_pathway_age_sum_stat <- cmd3_rlm_dir_age_sum_stat[,c(11:15)]
cmd3_rlm_dir_age_sum_stat <- cmd3_rlm_dir_age_sum_stat[,c(1:10)]

print("Study DFs created")

combined_dir_age_sum_stat <- as.data.frame(rbind(cmd3_rlm_dir_age_sum_stat,ag_rlm_dir_age_sum_stat,he_rlm_dir_age_sum_stat,odamaki_rlm_dir_age_sum_stat,nuage_rlm_dir_age_sum_stat,isc_rlm_dir_age_sum_stat,logmpie_rlm_dir_age_sum_stat))
combined_dir_age_sum_stat <- combined_dir_age_sum_stat[sorted_study_rows,]
#mat <- t(combined_dir_age_sum_stat)
#mat <- apply(mat,2,function(x)(ifelse(x<0,x+0.1,x)))
#heatmap.2(mat,density="none",trace="none",col=c("skyblue4","skyblue3","powderblue","mistyrose1","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),key=FALSE,margins=c(25,5),lhei=c(0.1,5),lwid=c(1,5))

print("Merged Taxa Directionalities")

combined_dir_pathway_age_sum_stat <- as.data.frame(rbind(cmd3_rlm_dir_pathway_age_sum_stat,isc_rlm_dir_pathway_age_sum_stat))
combined_dir_pathway_age_sum_stat <- combined_dir_pathway_age_sum_stat[sorted_study_rows_pathways,]
#mat <- t(combined_dir_pathway_age_sum_stat)
#mat <- apply(mat,2,function(x)(ifelse(x<0,x+0.1,x)))
#heatmap.2(mat,density="none",trace="none",col=c("skyblue4","skyblue3","powderblue","mistyrose1","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),key=FALSE,margins=c(25,5),lhei=c(0.1,5),lwid=c(1,5))

print("Merged Pathway Directionalities")

df_ag_diversity_uniqueness$study_name = "AG"
df_he_diversity_uniqueness$study_name = "HE"
df_nuage_diversity_uniqueness$study_name = "NUAGE"
df_odamaki_diversity_uniqueness$study_name = "Odamaki"
df_isc_diversity_uniqueness$study_name = "ISC"
df_logmpie_diversity_uniqueness$study_name = "LOGMPIE"

df_combined_diversity_uniqueness <- as.data.frame(rbind(df_cmd3_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_ag_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_he_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_nuage_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_odamaki_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_isc_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_logmpie_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")]))

print("Merged Taxonomy Diversity Uniqueness")

df_isc_pathway_diversity_uniqueness$study_name = "ISC"

df_combined_pathway_diversity_uniqueness <- as.data.frame(rbind(df_cmd3_diversity_uniqueness[,c("age","study_name","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")],df_isc_pathway_diversity_uniqueness[,c("age","study_name","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")]))

print("Merged Pathway Diversity Uniqueness")

df_ag_controls_diversity_uniqueness$study_name = "AG"
df_he_controls_diversity_uniqueness$study_name = "HE"
df_nuage_controls_diversity_uniqueness$study_name = "NUAGE"
df_odamaki_controls_diversity_uniqueness$study_name = "Odamaki"
df_isc_controls_diversity_uniqueness$study_name = "ISC"
df_logmpie_controls_diversity_uniqueness$study_name = "LOGMPIE"

df_combined_controls_diversity_uniqueness <- as.data.frame(rbind(df_cmd3_controls_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_ag_controls_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_he_controls_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_nuage_controls_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_odamaki_controls_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_isc_controls_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")],df_logmpie_controls_diversity_uniqueness[,c("age","study_name","genus_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")]))

print("Merged Control Taxonomy Diversity Uniqueness")

df_combined_controls_pathway_diversity_uniqueness <- as.data.frame(rbind(df_cmd3_controls_diversity_uniqueness[,c("age","study_name","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")],df_isc_pathway_diversity_uniqueness[,c("age","study_name","pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")]))

print("Merged Control Pathway Diversity Uniqueness")

Other_Studies <- rownames(combined_study_details_controls[!(combined_study_details_controls$region %in% c("EU_NA","EstAsia")),])
EU_NA_Studies <- rownames(combined_study_details_controls[(combined_study_details_controls$region == "EU_NA"),])
EstAsia_Studies <- rownames(combined_study_details_controls[(combined_study_details_controls$region == "EstAsia"),])

combined_rem_genus_shannon_age_Overall <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_shannon","age","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_genus_shannon_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_shannon","age","study_name",EU_NA_Studies)
combined_rem_genus_shannon_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_shannon","age","study_name",EstAsia_Studies)
combined_rem_genus_shannon_age_Other_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_shannon","age","study_name",Other_Studies)

combined_rem_species_shannon_age_Overall <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_shannon","age","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_species_shannon_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_shannon","age","study_name",EU_NA_Studies)
combined_rem_species_shannon_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_shannon","age","study_name",EstAsia_Studies)
combined_rem_species_shannon_age_Other_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_shannon","age","study_name",Other_Studies)

combined_rem_genus_bray_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_bray_uniqueness","age","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_genus_bray_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_bray_uniqueness","age","study_name",EU_NA_Studies)
combined_rem_genus_bray_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_bray_uniqueness","age","study_name",EstAsia_Studies)
combined_rem_genus_bray_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_bray_uniqueness","age","study_name",Other_Studies)

combined_rem_genus_jaccard_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_genus_jaccard_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","study_name",EU_NA_Studies)
combined_rem_genus_jaccard_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","study_name",EstAsia_Studies)
combined_rem_genus_jaccard_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","study_name",Other_Studies)

combined_rem_genus_aitchison_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_genus_aitchison_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","study_name",EU_NA_Studies)
combined_rem_genus_aitchison_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","study_name",EstAsia_Studies)
combined_rem_genus_aitchison_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","study_name",Other_Studies)

combined_rem_genus_kendall_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_genus_kendall_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","study_name",EU_NA_Studies)
combined_rem_genus_kendall_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","study_name",EstAsia_Studies)
combined_rem_genus_kendall_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","study_name",Other_Studies)

combined_rem_species_bray_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_bray_uniqueness","age","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_species_bray_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_bray_uniqueness","age","study_name",EU_NA_Studies)
combined_rem_species_bray_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_bray_uniqueness","age","study_name",EstAsia_Studies)
combined_rem_species_bray_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_bray_uniqueness","age","study_name",Other_Studies)

combined_rem_species_jaccard_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_species_jaccard_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","study_name",EU_NA_Studies)
combined_rem_species_jaccard_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","study_name",EstAsia_Studies)
combined_rem_species_jaccard_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","study_name",Other_Studies)

combined_rem_species_aitchison_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_species_aitchison_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","study_name",EU_NA_Studies)
combined_rem_species_aitchison_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","study_name",EstAsia_Studies)
combined_rem_species_aitchison_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","study_name",Other_Studies)

combined_rem_species_kendall_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_kendall_uniqueness","age","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_species_kendall_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_kendall_uniqueness","age","study_name",EU_NA_Studies)
combined_rem_species_kendall_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_kendall_uniqueness","age","study_name",EstAsia_Studies)
combined_rem_species_kendall_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_diversity_uniqueness,"species_kendall_uniqueness","age","study_name",Other_Studies)

print("Completed Taxonomy REMs")

Other_Studies_pathway <- intersect(sorted_study_rows_pathways,rownames(combined_study_details_controls[!(combined_study_details_controls$region %in% c("EU_NA","EstAsia")),]))
EU_NA_Studies_pathway <- intersect(sorted_study_rows_pathways,rownames(combined_study_details_controls[(combined_study_details_controls$region == "EU_NA"),]))
EstAsia_Studies_pathway <- intersect(sorted_study_rows_pathways,rownames(combined_study_details_controls[(combined_study_details_controls$region == "EstAsia"),]))

combined_rem_pathway_shannon_age_Overall <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_shannon","age","study_name",c(EU_NA_Studies_pathway,EstAsia_Studies_pathway,Other_Studies_pathway))
combined_rem_pathway_shannon_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_shannon","age","study_name",EU_NA_Studies_pathway)
combined_rem_pathway_shannon_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_shannon","age","study_name",EstAsia_Studies_pathway)
combined_rem_pathway_shannon_age_Other_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_shannon","age","study_name",Other_Studies_pathway)

combined_rem_pathway_bray_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_bray_uniqueness","age","study_name",c(EU_NA_Studies_pathway,EstAsia_Studies_pathway,Other_Studies_pathway))
combined_rem_pathway_bray_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_bray_uniqueness","age","study_name",EU_NA_Studies_pathway)
combined_rem_pathway_bray_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_bray_uniqueness","age","study_name",EstAsia_Studies_pathway)
combined_rem_pathway_bray_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_bray_uniqueness","age","study_name",Other_Studies_pathway)

combined_rem_pathway_jaccard_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_jaccard_uniqueness","age","study_name",c(EU_NA_Studies_pathway,EstAsia_Studies_pathway,Other_Studies_pathway))
combined_rem_pathway_jaccard_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_jaccard_uniqueness","age","study_name",EU_NA_Studies_pathway)
combined_rem_pathway_jaccard_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_jaccard_uniqueness","age","study_name",EstAsia_Studies_pathway)
combined_rem_pathway_jaccard_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_jaccard_uniqueness","age","study_name",Other_Studies_pathway)

combined_rem_pathway_aitchison_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_aitchison_uniqueness","age","study_name",c(EU_NA_Studies_pathway,EstAsia_Studies_pathway,Other_Studies_pathway))
combined_rem_pathway_aitchison_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_aitchison_uniqueness","age","study_name",EU_NA_Studies_pathway)
combined_rem_pathway_aitchison_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_aitchison_uniqueness","age","study_name",EstAsia_Studies_pathway)
combined_rem_pathway_aitchison_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_aitchison_uniqueness","age","study_name",Other_Studies_pathway)

combined_rem_pathway_kendall_uniqueness_age_Overall <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_kendall_uniqueness","age","study_name",c(EU_NA_Studies_pathway,EstAsia_Studies_pathway,Other_Studies_pathway))
combined_rem_pathway_kendall_uniqueness_age_EU_NA_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_kendall_uniqueness","age","study_name",EU_NA_Studies_pathway)
combined_rem_pathway_kendall_uniqueness_age_EstAsia_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_kendall_uniqueness","age","study_name",EstAsia_Studies_pathway)
combined_rem_pathway_kendall_uniqueness_age_Other_Studies <- compute_meta_lm(df_combined_controls_pathway_diversity_uniqueness,"pathway_kendall_uniqueness","age","study_name",Other_Studies_pathway)

print("Completed Pathway REMs")

cmd3_rlm_dir_pathway_diversity_uniqueness <- cmd3_rlm_dir_diversity_uniqueness[,c(9:12)]
cmd3_rlm_dir_diversity_uniqueness <- cmd3_rlm_dir_diversity_uniqueness[,c(1:8)]

combined_dir_diversity_uniqueness <- as.data.frame(rbind(cmd3_rlm_dir_diversity_uniqueness,ag_rlm_dir_diversity_uniqueness,he_rlm_dir_diversity_uniqueness,isc_rlm_dir_diversity_uniqueness,odamaki_rlm_dir_diversity_uniqueness,nuage_rlm_dir_diversity_uniqueness,logmpie_rlm_dir_diversity_uniqueness))
combined_dir_diversity_uniqueness <- combined_dir_diversity_uniqueness[sorted_study_rows,]
#mat <- t(combined_dir_diversity_uniqueness)
#mat <- apply(mat,2,function(x)(ifelse(x<0,x+0.1,x)))
#heatmap.2(mat,density="none",trace="none",col=c("skyblue4","skyblue3","powderblue","mistyrose1","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),key=FALSE,margins=c(25,5),lhei=c(0.1,5),lwid=c(1,5))

combined_dir_pathway_diversity_uniqueness <- as.data.frame(rbind(cmd3_rlm_dir_pathway_diversity_uniqueness,isc_rlm_dir_pathway_diversity_uniqueness))
combined_dir_pathway_diversity_uniqueness <- combined_dir_pathway_diversity_uniqueness[sorted_study_rows_pathways,]
#mat <- t(combined_dir_pathway_diversity_uniqueness)
#mat <- apply(mat,2,function(x)(ifelse(x<0,x+0.1,x)))
#heatmap.2(mat,density="none",trace="none",col=c("skyblue4","skyblue3","powderblue","mistyrose1","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),key=FALSE,margins=c(25,5),lhei=c(0.1,5),lwid=c(1,5))

cmd3_rlm_dir_pathway_age_ac_sum_stat <- cmd3_rlm_dir_age_ac_sum_stat[,c(11:15)]
cmd3_rlm_dir_age_ac_sum_stat <- cmd3_rlm_dir_age_ac_sum_stat[,c(1:10)]

combined_dir_age_ac_sum_stat <- as.data.frame(rbind(cmd3_rlm_dir_age_ac_sum_stat,ag_rlm_dir_age_ac_sum_stat,he_rlm_dir_age_ac_sum_stat,odamaki_rlm_dir_age_ac_sum_stat,nuage_rlm_dir_age_ac_sum_stat,isc_rlm_dir_age_ac_sum_stat,logmpie_rlm_dir_age_ac_sum_stat))
combined_dir_age_ac_sum_stat <- combined_dir_age_ac_sum_stat[sorted_study_rows,]
#mat <- t(combined_dir_age_ac_sum_stat)
#mat <- apply(mat,2,function(x)(ifelse(x<0,x+0.1,x)))
#heatmap.2(mat,density="none",trace="none",col=c("skyblue4","skyblue3","powderblue","mistyrose1","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),key=FALSE,margins=c(25,5),lhei=c(0.1,5),lwid=c(1,5))

combined_dir_pathway_age_ac_sum_stat <- as.data.frame(rbind(cmd3_rlm_dir_pathway_age_ac_sum_stat,isc_rlm_dir_pathway_age_ac_sum_stat))
combined_dir_pathway_age_ac_sum_stat <- combined_dir_pathway_age_ac_sum_stat[sorted_study_rows_pathways,]
#mat <- t(combined_dir_pathway_age_ac_sum_stat)
#mat <- apply(mat,2,function(x)(ifelse(x<0,x+0.1,x)))
#heatmap.2(mat,density="none",trace="none",col=c("skyblue4","skyblue3","powderblue","mistyrose1","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),key=FALSE,margins=c(25,5),lhei=c(0.1,5),lwid=c(1,5))

print("Computing for alpha-adjusted bray uniqueness (genus)")
combined_rem_ac_genus_bray_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_bray_uniqueness","age","genus_shannon","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_ac_genus_bray_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_bray_uniqueness","age","genus_shannon","study_name",EU_NA_Studies)
combined_rem_ac_genus_bray_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_bray_uniqueness","age","genus_shannon","study_name",EstAsia_Studies)
combined_rem_ac_genus_bray_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_bray_uniqueness","age","genus_shannon","study_name",Other_Studies)

print("Computing for alpha-adjusted jaccard uniqueness (genus)")
combined_rem_ac_genus_jaccard_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","genus_shannon","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_ac_genus_jaccard_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","genus_shannon","study_name",EU_NA_Studies)
combined_rem_ac_genus_jaccard_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","genus_shannon","study_name",EstAsia_Studies)
combined_rem_ac_genus_jaccard_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_jaccard_uniqueness","age","genus_shannon","study_name",Other_Studies)

print("Computing for alpha-adjusted aitchison uniqueness (genus)")
combined_rem_ac_genus_aitchison_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","genus_shannon","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_ac_genus_aitchison_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","genus_shannon","study_name",EU_NA_Studies)
combined_rem_ac_genus_aitchison_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","genus_shannon","study_name",EstAsia_Studies)
combined_rem_ac_genus_aitchison_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_aitchison_uniqueness","age","genus_shannon","study_name",Other_Studies)

print("Computing for alpha-adjusted kendall uniqueness (genus)")
combined_rem_ac_genus_kendall_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","genus_shannon","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_ac_genus_kendall_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","genus_shannon","study_name",EU_NA_Studies)
combined_rem_ac_genus_kendall_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","genus_shannon","study_name",EstAsia_Studies)
combined_rem_ac_genus_kendall_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"genus_kendall_uniqueness","age","genus_shannon","study_name",Other_Studies)

print("Computing for alpha-adjusted bray uniqueness (species)")
combined_rem_ac_species_bray_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_bray_uniqueness","age","species_shannon","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_ac_species_bray_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_bray_uniqueness","age","species_shannon","study_name",EU_NA_Studies)
combined_rem_ac_species_bray_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_bray_uniqueness","age","species_shannon","study_name",EstAsia_Studies)
combined_rem_ac_species_bray_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_bray_uniqueness","age","species_shannon","study_name",Other_Studies)

print("Computing for alpha-adjusted jaccard uniqueness (species)")
combined_rem_ac_species_jaccard_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","species_shannon","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_ac_species_jaccard_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","species_shannon","study_name",EU_NA_Studies)
combined_rem_ac_species_jaccard_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","species_shannon","study_name",EstAsia_Studies)
combined_rem_ac_species_jaccard_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_jaccard_uniqueness","age","species_shannon","study_name",Other_Studies)

print("Computing for alpha-adjusted aitchison uniqueness (species)")
combined_rem_ac_species_aitchison_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","species_shannon","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_ac_species_aitchison_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","species_shannon","study_name",EU_NA_Studies)
combined_rem_ac_species_aitchison_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","species_shannon","study_name",EstAsia_Studies)
combined_rem_ac_species_aitchison_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_aitchison_uniqueness","age","species_shannon","study_name",Other_Studies)

print("Computing for alpha-adjusted kendall uniqueness (species)")
combined_rem_ac_species_kendall_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_kendall_uniqueness","age","species_shannon","study_name",c(EU_NA_Studies,EstAsia_Studies,Other_Studies))
combined_rem_ac_species_kendall_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_kendall_uniqueness","age","species_shannon","study_name",EU_NA_Studies)
combined_rem_ac_species_kendall_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_kendall_uniqueness","age","species_shannon","study_name",EstAsia_Studies)
combined_rem_ac_species_kendall_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_diversity_uniqueness,"species_kendall_uniqueness","age","species_shannon","study_name",Other_Studies)

print("Computing for alpha-adjusted bray uniqueness (pathways)")
combined_rem_ac_pathway_bray_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_bray_uniqueness","age","pathway_shannon","study_name",c(EU_NA_Studies_pathway,EstAsia_Studies_pathway,Other_Studies_pathway))
combined_rem_ac_pathway_bray_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_bray_uniqueness","age","pathway_shannon","study_name",EU_NA_Studies_pathway)
combined_rem_ac_pathway_bray_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_bray_uniqueness","age","pathway_shannon","study_name",EstAsia_Studies_pathway)
combined_rem_ac_pathway_bray_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_bray_uniqueness","age","pathway_shannon","study_name",Other_Studies_pathway)

print("Computing for alpha-adjusted jaccard uniqueness (pathways)")
combined_rem_ac_pathway_jaccard_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_jaccard_uniqueness","age","pathway_shannon","study_name",c(EU_NA_Studies_pathway,EstAsia_Studies_pathway,Other_Studies_pathway))
combined_rem_ac_pathway_jaccard_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_jaccard_uniqueness","age","pathway_shannon","study_name",EU_NA_Studies_pathway)
combined_rem_ac_pathway_jaccard_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_jaccard_uniqueness","age","pathway_shannon","study_name",EstAsia_Studies_pathway)
combined_rem_ac_pathway_jaccard_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_jaccard_uniqueness","age","pathway_shannon","study_name",Other_Studies_pathway)

print("Computing for alpha-adjusted aitchison uniqueness (pathways)")
combined_rem_ac_pathway_aitchison_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_aitchison_uniqueness","age","pathway_shannon","study_name",c(EU_NA_Studies_pathway,EstAsia_Studies_pathway,Other_Studies_pathway))
combined_rem_ac_pathway_aitchison_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_aitchison_uniqueness","age","pathway_shannon","study_name",EU_NA_Studies_pathway)
combined_rem_ac_pathway_aitchison_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_aitchison_uniqueness","age","pathway_shannon","study_name",EstAsia_Studies_pathway)
combined_rem_ac_pathway_aitchison_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_aitchison_uniqueness","age","pathway_shannon","study_name",Other_Studies_pathway)

print("Computing for alpha-adjusted kendall uniqueness (pathways)")
combined_rem_ac_pathway_kendall_uniqueness_age_Overall <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_kendall_uniqueness","age","pathway_shannon","study_name",c(EU_NA_Studies_pathway,EstAsia_Studies_pathway,Other_Studies_pathway))
combined_rem_ac_pathway_kendall_uniqueness_age_EU_NA_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_kendall_uniqueness","age","pathway_shannon","study_name",EU_NA_Studies_pathway)
combined_rem_ac_pathway_kendall_uniqueness_age_EstAsia_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_kendall_uniqueness","age","pathway_shannon","study_name",EstAsia_Studies_pathway)
combined_rem_ac_pathway_kendall_uniqueness_age_Other_Studies <- compute_meta_lm_single_adjust(df_combined_controls_pathway_diversity_uniqueness,"pathway_kendall_uniqueness","age","pathway_shannon","study_name",Other_Studies_pathway)

print("Completed Alpha adjusted Uniqueness REMs")

combined_rlm_est_diversity_uniqueness <- as.data.frame(rbind(cmd3_rlm_est_diversity_uniqueness[,c(1:8)],ag_rlm_est_diversity_uniqueness,nuage_rlm_est_diversity_uniqueness,isc_rlm_est_diversity_uniqueness,he_rlm_est_diversity_uniqueness,odamaki_rlm_est_diversity_uniqueness,logmpie_rlm_est_diversity_uniqueness))
write.table(combined_rlm_est_diversity_uniqueness[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_est_diversity_uniqueness.txt",row.names=T,quote=FALSE,col.names=T)

combined_rlm_p_val_diversity_uniqueness <- as.data.frame(rbind(cmd3_rlm_p_val_diversity_uniqueness[,c(1:8)],ag_rlm_p_val_diversity_uniqueness,nuage_rlm_p_val_diversity_uniqueness,isc_rlm_p_val_diversity_uniqueness,he_rlm_p_val_diversity_uniqueness,odamaki_rlm_p_val_diversity_uniqueness,logmpie_rlm_p_val_diversity_uniqueness))
write.table(combined_rlm_p_val_diversity_uniqueness[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_p_val_diversity_uniqueness.txt",row.names=T,quote=FALSE,col.names=T)

combined_rlm_q_val_diversity_uniqueness <- as.data.frame(rbind(cmd3_rlm_q_val_diversity_uniqueness[,c(1:8)],ag_rlm_q_val_diversity_uniqueness,nuage_rlm_q_val_diversity_uniqueness,isc_rlm_q_val_diversity_uniqueness,he_rlm_q_val_diversity_uniqueness,odamaki_rlm_q_val_diversity_uniqueness,logmpie_rlm_q_val_diversity_uniqueness))
write.table(combined_rlm_q_val_diversity_uniqueness[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_q_val_diversity_uniqueness.txt",row.names=T,quote=FALSE,col.names=T)

combined_rlm_est_pathways_diversity_uniqueness <- as.data.frame(rbind(cmd3_rlm_est_diversity_uniqueness[,c(9:12)],isc_rlm_est_pathway_diversity_uniqueness))
write.table(combined_rlm_est_pathways_diversity_uniqueness[sorted_study_rows_pathways,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_est_pathways_diversity_uniqueness.txt",row.names=T,quote=FALSE,col.names=T)

combined_rlm_p_val_pathways_diversity_uniqueness <- as.data.frame(rbind(cmd3_rlm_p_val_diversity_uniqueness[,c(9:12)],isc_rlm_p_val_pathway_diversity_uniqueness))
write.table(combined_rlm_p_val_pathways_diversity_uniqueness[sorted_study_rows_pathways,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_p_val_pathways_diversity_uniqueness.txt",row.names=T,quote=FALSE,col.names=T)

combined_rlm_q_val_pathways_diversity_uniqueness <- as.data.frame(rbind(cmd3_rlm_q_val_diversity_uniqueness[,c(9:12)],isc_rlm_q_val_pathway_diversity_uniqueness))
write.table(combined_rlm_q_val_pathways_diversity_uniqueness[sorted_study_rows_pathways,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_q_val_pathways_diversity_uniqueness.txt",row.names=T,quote=FALSE,col.names=T)

combined_rlm_genus_kendall_bray <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_genus_kendall_bray) <- sorted_study_rows
colnames(combined_rlm_genus_kendall_bray) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(genus_kendall_uniqueness~genus_bray_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_genus_kendall_bray[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_genus_kendall_bray[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_genus_kendall_bray[,3] <- p.adjust(combined_rlm_genus_kendall_bray[,2],method="fdr")

combined_rlm_genus_kendall_jaccard <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_genus_kendall_jaccard) <- sorted_study_rows
colnames(combined_rlm_genus_kendall_jaccard) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(genus_kendall_uniqueness~genus_jaccard_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_genus_kendall_jaccard[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_genus_kendall_jaccard[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_genus_kendall_jaccard[,3] <- p.adjust(combined_rlm_genus_kendall_jaccard[,2],method="fdr")

combined_rlm_genus_kendall_aitchison <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_genus_kendall_aitchison) <- sorted_study_rows
colnames(combined_rlm_genus_kendall_aitchison) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(genus_kendall_uniqueness~genus_aitchison_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_genus_kendall_aitchison[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_genus_kendall_aitchison[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_genus_kendall_aitchison[,3] <- p.adjust(combined_rlm_genus_kendall_aitchison[,2],method="fdr")

combined_rlm_genus_aitchison_bray <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_genus_aitchison_bray) <- sorted_study_rows
colnames(combined_rlm_genus_aitchison_bray) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(genus_aitchison_uniqueness~genus_bray_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_genus_aitchison_bray[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_genus_aitchison_bray[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_genus_aitchison_bray[,3] <- p.adjust(combined_rlm_genus_aitchison_bray[,2],method="fdr")

combined_rlm_genus_aitchison_jaccard <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_genus_aitchison_jaccard) <- sorted_study_rows
colnames(combined_rlm_genus_aitchison_jaccard) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(genus_aitchison_uniqueness~genus_jaccard_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_genus_aitchison_jaccard[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_genus_aitchison_jaccard[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_genus_aitchison_jaccard[,3] <- p.adjust(combined_rlm_genus_aitchison_jaccard[,2],method="fdr")

combined_rlm_genus_bray_jaccard <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_genus_bray_jaccard) <- sorted_study_rows
colnames(combined_rlm_genus_bray_jaccard) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(genus_bray_uniqueness~genus_jaccard_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_genus_bray_jaccard[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_genus_bray_jaccard[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_genus_bray_jaccard[,3] <- p.adjust(combined_rlm_genus_bray_jaccard[,2],method="fdr")

combined_rlm_species_kendall_bray <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_species_kendall_bray) <- sorted_study_rows
colnames(combined_rlm_species_kendall_bray) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(species_kendall_uniqueness~species_bray_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_species_kendall_bray[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_species_kendall_bray[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_species_kendall_bray[,3] <- p.adjust(combined_rlm_species_kendall_bray[,2],method="fdr")

combined_rlm_species_kendall_jaccard <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_species_kendall_jaccard) <- sorted_study_rows
colnames(combined_rlm_species_kendall_jaccard) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(species_kendall_uniqueness~species_jaccard_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_species_kendall_jaccard[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_species_kendall_jaccard[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_species_kendall_jaccard[,3] <- p.adjust(combined_rlm_species_kendall_jaccard[,2],method="fdr")

combined_rlm_species_kendall_aitchison <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_species_kendall_aitchison) <- sorted_study_rows
colnames(combined_rlm_species_kendall_aitchison) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(species_kendall_uniqueness~species_aitchison_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_species_kendall_aitchison[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_species_kendall_aitchison[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_species_kendall_aitchison[,3] <- p.adjust(combined_rlm_species_kendall_aitchison[,2],method="fdr")

combined_rlm_species_aitchison_bray <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_species_aitchison_bray) <- sorted_study_rows
colnames(combined_rlm_species_aitchison_bray) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(species_aitchison_uniqueness~species_bray_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_species_aitchison_bray[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_species_aitchison_bray[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_species_aitchison_bray[,3] <- p.adjust(combined_rlm_species_aitchison_bray[,2],method="fdr")

combined_rlm_species_aitchison_jaccard <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_species_aitchison_jaccard) <- sorted_study_rows
colnames(combined_rlm_species_aitchison_jaccard) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(species_aitchison_uniqueness~species_jaccard_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_species_aitchison_jaccard[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_species_aitchison_jaccard[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_species_aitchison_jaccard[,3] <- p.adjust(combined_rlm_species_aitchison_jaccard[,2],method="fdr")

combined_rlm_species_bray_jaccard <- as.data.frame(matrix(NA,length(sorted_study_rows),3))
rownames(combined_rlm_species_bray_jaccard) <- sorted_study_rows
colnames(combined_rlm_species_bray_jaccard) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows))
{
	study_name <- sorted_study_rows[i]
	temp_rlm <- rlm(species_bray_uniqueness~species_jaccard_uniqueness,df_combined_diversity_uniqueness[df_combined_diversity_uniqueness$study_name==study_name,])
	combined_rlm_species_bray_jaccard[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_species_bray_jaccard[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_species_bray_jaccard[,3] <- p.adjust(combined_rlm_species_bray_jaccard[,2],method="fdr")

combined_rlm_pathway_kendall_bray <- as.data.frame(matrix(NA,length(sorted_study_rows_pathways),3))
rownames(combined_rlm_pathway_kendall_bray) <- sorted_study_rows_pathways
colnames(combined_rlm_pathway_kendall_bray) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows_pathways))
{
	study_name <- sorted_study_rows_pathways[i]
	temp_rlm <- rlm(pathway_kendall_uniqueness~pathway_bray_uniqueness,df_combined_pathway_diversity_uniqueness[df_combined_pathway_diversity_uniqueness$study_name==study_name,])
	combined_rlm_pathway_kendall_bray[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_pathway_kendall_bray[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_pathway_kendall_bray[,3] <- p.adjust(combined_rlm_pathway_kendall_bray[,2],method="fdr")

combined_rlm_pathway_kendall_jaccard <- as.data.frame(matrix(NA,length(sorted_study_rows_pathways),3))
rownames(combined_rlm_pathway_kendall_jaccard) <- sorted_study_rows_pathways
colnames(combined_rlm_pathway_kendall_jaccard) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows_pathways))
{
	study_name <- sorted_study_rows_pathways[i]
	temp_rlm <- rlm(pathway_kendall_uniqueness~pathway_jaccard_uniqueness,df_combined_pathway_diversity_uniqueness[df_combined_pathway_diversity_uniqueness$study_name==study_name,])
	combined_rlm_pathway_kendall_jaccard[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_pathway_kendall_jaccard[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_pathway_kendall_jaccard[,3] <- p.adjust(combined_rlm_pathway_kendall_jaccard[,2],method="fdr")

combined_rlm_pathway_kendall_aitchison <- as.data.frame(matrix(NA,length(sorted_study_rows_pathways),3))
rownames(combined_rlm_pathway_kendall_aitchison) <- sorted_study_rows_pathways
colnames(combined_rlm_pathway_kendall_aitchison) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows_pathways))
{
	study_name <- sorted_study_rows_pathways[i]
	temp_rlm <- rlm(pathway_kendall_uniqueness~pathway_aitchison_uniqueness,df_combined_pathway_diversity_uniqueness[df_combined_pathway_diversity_uniqueness$study_name==study_name,])
	combined_rlm_pathway_kendall_aitchison[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_pathway_kendall_aitchison[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_pathway_kendall_aitchison[,3] <- p.adjust(combined_rlm_pathway_kendall_aitchison[,2],method="fdr")

combined_rlm_pathway_aitchison_bray <- as.data.frame(matrix(NA,length(sorted_study_rows_pathways),3))
rownames(combined_rlm_pathway_aitchison_bray) <- sorted_study_rows_pathways
colnames(combined_rlm_pathway_aitchison_bray) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows_pathways))
{
	study_name <- sorted_study_rows_pathways[i]
	temp_rlm <- rlm(pathway_aitchison_uniqueness~pathway_bray_uniqueness,df_combined_pathway_diversity_uniqueness[df_combined_pathway_diversity_uniqueness$study_name==study_name,])
	combined_rlm_pathway_aitchison_bray[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_pathway_aitchison_bray[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_pathway_aitchison_bray[,3] <- p.adjust(combined_rlm_pathway_aitchison_bray[,2],method="fdr")

combined_rlm_pathway_aitchison_jaccard <- as.data.frame(matrix(NA,length(sorted_study_rows_pathways),3))
rownames(combined_rlm_pathway_aitchison_jaccard) <- sorted_study_rows_pathways
colnames(combined_rlm_pathway_aitchison_jaccard) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows_pathways))
{
	study_name <- sorted_study_rows_pathways[i]
	temp_rlm <- rlm(pathway_aitchison_uniqueness~pathway_jaccard_uniqueness,df_combined_pathway_diversity_uniqueness[df_combined_pathway_diversity_uniqueness$study_name==study_name,])
	combined_rlm_pathway_aitchison_jaccard[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_pathway_aitchison_jaccard[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_pathway_aitchison_jaccard[,3] <- p.adjust(combined_rlm_pathway_aitchison_jaccard[,2],method="fdr")

combined_rlm_pathway_bray_jaccard <- as.data.frame(matrix(NA,length(sorted_study_rows_pathways),3))
rownames(combined_rlm_pathway_bray_jaccard) <- sorted_study_rows_pathways
colnames(combined_rlm_pathway_bray_jaccard) <- c("est","p_val","q_val")
for(i in 1:length(sorted_study_rows_pathways))
{
	study_name <- sorted_study_rows_pathways[i]
	temp_rlm <- rlm(pathway_bray_uniqueness~pathway_jaccard_uniqueness,df_combined_pathway_diversity_uniqueness[df_combined_pathway_diversity_uniqueness$study_name==study_name,])
	combined_rlm_pathway_bray_jaccard[i,1] <- as.numeric(temp_rlm$coefficients[2])
	combined_rlm_pathway_bray_jaccard[i,2] <- f.robftest(temp_rlm)$p.value
}
combined_rlm_pathway_bray_jaccard[,3] <- p.adjust(combined_rlm_pathway_bray_jaccard[,2],method="fdr")

combined_est_age_beta <- as.data.frame(rbind(cmd3_est_age_beta[,c(3:10)],ag_est_age_beta[,c(3:10)],nuage_est_age_beta[,c(3:10)],isc_est_age_beta[,c(3:10)],he_est_age_beta[,c(3:10)],odamaki_est_age_beta[,c(3:10)],logmpie_est_age_beta[,c(3:10)]))
write.table(combined_est_age_beta[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_est_age_beta.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)
write.table(combined_p_val_age_beta[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_p_val_age_beta.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_pathway_p_val_age_beta <- as.data.frame(rbind(cmd3_pathway_p_val_age_beta[,c(2:5)],isc_pathway_p_val_age_beta[,c(2:5)]))
write.table(combined_pathway_p_val_age_beta[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\ccombined_pathway_p_val_age_beta.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)
combined_pathway_est_age_beta <- as.data.frame(rbind(cmd3_pathway_est_age_beta[,c(2:5)],isc_pathway_est_age_beta[,c(2:5)]))
write.table(combined_pathway_est_age_beta[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\ccombined_pathway_est_age_beta.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_est_age_sum_stat <- as.data.frame(rbind(cmd3_rlm_est_age_sum_stat,ag_rlm_est_age_sum_stat,nuage_rlm_est_age_sum_stat,isc_rlm_est_age_sum_stat,he_rlm_est_age_sum_stat,odamaki_rlm_est_age_sum_stat,logmpie_rlm_est_age_sum_stat))
write.table(combined_rlm_est_age_sum_stat[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_est_age_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_p_val_age_sum_stat <- as.data.frame(rbind(cmd3_rlm_p_val_age_sum_stat,ag_rlm_p_val_age_sum_stat,nuage_rlm_p_val_age_sum_stat,isc_rlm_p_val_age_sum_stat,he_rlm_p_val_age_sum_stat,odamaki_rlm_p_val_age_sum_stat,logmpie_rlm_p_val_age_sum_stat))
write.table(combined_rlm_p_val_age_sum_stat[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_p_val_age_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_q_val_age_sum_stat <- as.data.frame(rbind(cmd3_rlm_q_val_age_sum_stat,ag_rlm_q_val_age_sum_stat,nuage_rlm_q_val_age_sum_stat,isc_rlm_q_val_age_sum_stat,he_rlm_q_val_age_sum_stat,odamaki_rlm_q_val_age_sum_stat,logmpie_rlm_q_val_age_sum_stat))
write.table(combined_rlm_q_val_age_sum_stat[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_q_val_age_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_est_pathway_age_sum_stat <- as.data.frame(rbind(cmd3_rlm_est_pathway_age_sum_stat,isc_rlm_est_pathway_age_sum_stat))
write.table(combined_rlm_est_pathway_age_sum_stat[sorted_study_rows_pathways,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_est_pathway_age_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_p_val_pathway_age_sum_stat <- as.data.frame(rbind(cmd3_rlm_p_val_pathway_age_sum_stat,isc_rlm_p_val_pathway_age_sum_stat))
write.table(combined_rlm_p_val_pathway_age_sum_stat[sorted_study_rows_pathways,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_p_val_pathway_age_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_q_val_pathway_age_sum_stat <- as.data.frame(rbind(cmd3_rlm_q_val_pathway_age_sum_stat,isc_rlm_q_val_pathway_age_sum_stat))
write.table(combined_rlm_q_val_pathway_age_sum_stat[sorted_study_rows_pathways,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_q_val_pathway_age_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_est_age_ac_sum_stat <- as.data.frame(rbind(cmd3_rlm_est_age_ac_sum_stat[,c(1:10)],ag_rlm_est_age_ac_sum_stat,nuage_rlm_est_age_ac_sum_stat,isc_rlm_est_age_ac_sum_stat,he_rlm_est_age_ac_sum_stat,odamaki_rlm_est_age_ac_sum_stat,logmpie_rlm_est_age_ac_sum_stat))
write.table(combined_rlm_est_age_ac_sum_stat[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_est_age_ac_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_p_val_age_ac_sum_stat <- as.data.frame(rbind(cmd3_rlm_p_val_age_ac_sum_stat[,c(1:10)],ag_rlm_p_val_age_ac_sum_stat,nuage_rlm_p_val_age_ac_sum_stat,isc_rlm_p_val_age_ac_sum_stat,he_rlm_p_val_age_ac_sum_stat,odamaki_rlm_p_val_age_ac_sum_stat,logmpie_rlm_p_val_age_ac_sum_stat))
write.table(combined_rlm_p_val_age_ac_sum_stat[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_p_val_age_ac_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_q_val_age_ac_sum_stat <- as.data.frame(rbind(cmd3_rlm_q_val_age_ac_sum_stat[,c(1:10)],ag_rlm_q_val_age_ac_sum_stat,nuage_rlm_q_val_age_ac_sum_stat,isc_rlm_q_val_age_ac_sum_stat,he_rlm_q_val_age_ac_sum_stat,odamaki_rlm_q_val_age_ac_sum_stat,logmpie_rlm_q_val_age_ac_sum_stat))
write.table(combined_rlm_q_val_age_ac_sum_stat[sorted_study_rows,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_q_val_age_ac_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_est_pathway_age_ac_sum_stat <- as.data.frame(rbind(cmd3_rlm_est_age_ac_sum_stat[,c(11:15)],isc_rlm_est_pathway_age_ac_sum_stat))
write.table(combined_rlm_est_pathway_age_ac_sum_stat[sorted_study_rows_pathways,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_est_pathway_age_ac_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_p_val_pathway_age_ac_sum_stat <- as.data.frame(rbind(cmd3_rlm_p_val_age_ac_sum_stat[,c(11:15)],isc_rlm_p_val_pathway_age_ac_sum_stat))
write.table(combined_rlm_p_val_pathway_age_ac_sum_stat[sorted_study_rows_pathways,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_p_val_pathway_age_ac_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

combined_rlm_q_val_pathway_age_ac_sum_stat <- as.data.frame(rbind(cmd3_rlm_q_val_age_ac_sum_stat[,c(11:15)],isc_rlm_q_val_pathway_age_ac_sum_stat))
write.table(combined_rlm_q_val_pathway_age_ac_sum_stat[sorted_study_rows_pathways,],file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\RawOutputs\\combined_rlm_q_val_pathway_age_ac_sum_stat.txt",sep="\t",row.names=T,quote=FALSE,col.names=T)

rem_est_age_sum_stat <- as.data.frame(matrix(NA,10,4))
rownames(rem_est_age_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")
colnames(rem_est_age_sum_stat) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")


rem_est_age_sum_stat["genus_shannon","Overall"] <- as.numeric(combined_rem_genus_shannon_age_Overall$model$b)
rem_est_age_sum_stat["species_shannon","Overall"] <- as.numeric(combined_rem_species_shannon_age_Overall$model$b)
rem_est_age_sum_stat["genus_bray_uniqueness","Overall"] <- as.numeric(combined_rem_genus_bray_uniqueness_age_Overall$model$b)
rem_est_age_sum_stat["genus_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_genus_jaccard_uniqueness_age_Overall$model$b)
rem_est_age_sum_stat["genus_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_genus_aitchison_uniqueness_age_Overall$model$b)
rem_est_age_sum_stat["genus_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_genus_kendall_uniqueness_age_Overall$model$b)
rem_est_age_sum_stat["species_bray_uniqueness","Overall"] <- as.numeric(combined_rem_species_bray_uniqueness_age_Overall$model$b)
rem_est_age_sum_stat["species_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_species_jaccard_uniqueness_age_Overall$model$b)
rem_est_age_sum_stat["species_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_species_aitchison_uniqueness_age_Overall$model$b)
rem_est_age_sum_stat["species_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_species_kendall_uniqueness_age_Overall$model$b)

rem_est_age_sum_stat["genus_shannon","EU_NA_Studies"] <- as.numeric(combined_rem_genus_shannon_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat["species_shannon","EU_NA_Studies"] <- as.numeric(combined_rem_species_shannon_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat["genus_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_genus_bray_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat["genus_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_genus_jaccard_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat["genus_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_genus_aitchison_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat["genus_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_genus_kendall_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat["species_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_species_bray_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat["species_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_species_jaccard_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat["species_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_species_aitchison_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat["species_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_species_kendall_uniqueness_age_EU_NA_Studies$model$b)

rem_est_age_sum_stat["genus_shannon","EstAsia_Studies"] <- as.numeric(combined_rem_genus_shannon_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat["species_shannon","EstAsia_Studies"] <- as.numeric(combined_rem_species_shannon_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat["genus_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_genus_bray_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat["genus_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_genus_jaccard_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat["genus_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_genus_aitchison_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat["genus_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_genus_kendall_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat["species_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_species_bray_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat["species_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_species_jaccard_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat["species_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_species_aitchison_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat["species_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_species_kendall_uniqueness_age_EstAsia_Studies$model$b)

rem_est_age_sum_stat["genus_shannon","Other_Studies"] <- as.numeric(combined_rem_genus_shannon_age_Other_Studies$model$b)
rem_est_age_sum_stat["species_shannon","Other_Studies"] <- as.numeric(combined_rem_species_shannon_age_Other_Studies$model$b)
rem_est_age_sum_stat["genus_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_genus_bray_uniqueness_age_Other_Studies$model$b)
rem_est_age_sum_stat["genus_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_genus_jaccard_uniqueness_age_Other_Studies$model$b)
rem_est_age_sum_stat["genus_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_genus_aitchison_uniqueness_age_Other_Studies$model$b)
rem_est_age_sum_stat["genus_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_genus_kendall_uniqueness_age_Other_Studies$model$b)
rem_est_age_sum_stat["species_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_species_bray_uniqueness_age_Other_Studies$model$b)
rem_est_age_sum_stat["species_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_species_jaccard_uniqueness_age_Other_Studies$model$b)
rem_est_age_sum_stat["species_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_species_aitchison_uniqueness_age_Other_Studies$model$b)
rem_est_age_sum_stat["species_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_species_kendall_uniqueness_age_Other_Studies$model$b)

rem_p_val_age_sum_stat <- as.data.frame(matrix(NA,10,4))
rownames(rem_p_val_age_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")
colnames(rem_p_val_age_sum_stat) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")

rem_p_val_age_sum_stat["genus_shannon","Overall"] <- as.numeric(combined_rem_genus_shannon_age_Overall$model$pval)
rem_p_val_age_sum_stat["species_shannon","Overall"] <- as.numeric(combined_rem_species_shannon_age_Overall$model$pval)
rem_p_val_age_sum_stat["genus_bray_uniqueness","Overall"] <- as.numeric(combined_rem_genus_bray_uniqueness_age_Overall$model$pval)
rem_p_val_age_sum_stat["genus_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_genus_jaccard_uniqueness_age_Overall$model$pval)
rem_p_val_age_sum_stat["genus_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_genus_aitchison_uniqueness_age_Overall$model$pval)
rem_p_val_age_sum_stat["genus_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_genus_kendall_uniqueness_age_Overall$model$pval)
rem_p_val_age_sum_stat["species_bray_uniqueness","Overall"] <- as.numeric(combined_rem_species_bray_uniqueness_age_Overall$model$pval)
rem_p_val_age_sum_stat["species_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_species_jaccard_uniqueness_age_Overall$model$pval)
rem_p_val_age_sum_stat["species_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_species_aitchison_uniqueness_age_Overall$model$pval)
rem_p_val_age_sum_stat["species_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_species_kendall_uniqueness_age_Overall$model$pval)

rem_p_val_age_sum_stat["genus_shannon","EU_NA_Studies"] <- as.numeric(combined_rem_genus_shannon_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat["species_shannon","EU_NA_Studies"] <- as.numeric(combined_rem_species_shannon_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat["genus_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_genus_bray_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat["genus_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_genus_jaccard_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat["genus_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_genus_aitchison_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat["genus_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_genus_kendall_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat["species_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_species_bray_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat["species_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_species_jaccard_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat["species_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_species_aitchison_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat["species_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_species_kendall_uniqueness_age_EU_NA_Studies$model$pval)

rem_p_val_age_sum_stat["genus_shannon","EstAsia_Studies"] <- as.numeric(combined_rem_genus_shannon_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat["species_shannon","EstAsia_Studies"] <- as.numeric(combined_rem_species_shannon_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat["genus_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_genus_bray_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat["genus_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_genus_jaccard_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat["genus_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_genus_aitchison_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat["genus_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_genus_kendall_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat["species_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_species_bray_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat["species_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_species_jaccard_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat["species_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_species_aitchison_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat["species_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_species_kendall_uniqueness_age_EstAsia_Studies$model$pval)

rem_p_val_age_sum_stat["genus_shannon","Other_Studies"] <- as.numeric(combined_rem_genus_shannon_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat["species_shannon","Other_Studies"] <- as.numeric(combined_rem_species_shannon_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat["genus_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_genus_bray_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat["genus_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_genus_jaccard_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat["genus_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_genus_aitchison_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat["genus_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_genus_kendall_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat["species_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_species_bray_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat["species_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_species_jaccard_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat["species_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_species_aitchison_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat["species_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_species_kendall_uniqueness_age_Other_Studies$model$pval)

rem_q_val_age_sum_stat <- t(apply(rem_p_val_age_sum_stat,1,function(x)(p.adjust(x,method="fdr"))))

rem_dir_age_sum_stat <- as.data.frame(matrix(NA,10,4))
rownames(rem_dir_age_sum_stat) <- c("genus_shannon","species_shannon","genus_bray_uniqueness","genus_jaccard_uniqueness","genus_aitchison_uniqueness","genus_kendall_uniqueness","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness")
colnames(rem_dir_age_sum_stat) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")
for(i in 1:nrow(rem_dir_age_sum_stat))
{
	for(j in 1:ncol(rem_dir_age_sum_stat))
	{
		rem_dir_age_sum_stat[i,j] <- ifelse(rem_q_val_age_sum_stat[i,j] <= 0.1,3*sign(rem_est_age_sum_stat[i,j]),ifelse(rem_p_val_age_sum_stat[i,j] <= 0.05,2*sign(rem_est_age_sum_stat[i,j]),sign(rem_est_age_sum_stat[i,j])))
	}
}

rem_est_age_ac_sum_stat <- as.data.frame(matrix(NA,8,4))
rownames(rem_est_age_ac_sum_stat) <- c("ac_genus_bray_uniqueness","ac_genus_jaccard_uniqueness","ac_genus_aitchison_uniqueness","ac_genus_kendall_uniqueness","ac_species_bray_uniqueness","ac_species_jaccard_uniqueness","ac_species_aitchison_uniqueness","ac_species_kendall_uniqueness")
colnames(rem_est_age_ac_sum_stat) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")

rem_est_age_ac_sum_stat["ac_genus_bray_uniqueness","Overall"] <- as.numeric(combined_rem_ac_genus_bray_uniqueness_age_Overall$model$b)
rem_est_age_ac_sum_stat["ac_genus_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_ac_genus_jaccard_uniqueness_age_Overall$model$b)
rem_est_age_ac_sum_stat["ac_genus_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_ac_genus_aitchison_uniqueness_age_Overall$model$b)
rem_est_age_ac_sum_stat["ac_genus_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_ac_genus_kendall_uniqueness_age_Overall$model$b)
rem_est_age_ac_sum_stat["ac_species_bray_uniqueness","Overall"] <- as.numeric(combined_rem_ac_species_bray_uniqueness_age_Overall$model$b)
rem_est_age_ac_sum_stat["ac_species_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_ac_species_jaccard_uniqueness_age_Overall$model$b)
rem_est_age_ac_sum_stat["ac_species_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_ac_species_aitchison_uniqueness_age_Overall$model$b)
rem_est_age_ac_sum_stat["ac_species_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_ac_species_kendall_uniqueness_age_Overall$model$b)

rem_est_age_ac_sum_stat["ac_genus_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_genus_bray_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_ac_sum_stat["ac_genus_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_genus_jaccard_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_ac_sum_stat["ac_genus_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_genus_aitchison_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_ac_sum_stat["ac_genus_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_genus_kendall_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_species_bray_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_species_jaccard_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_species_aitchison_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_species_kendall_uniqueness_age_EU_NA_Studies$model$b)

rem_est_age_ac_sum_stat["ac_genus_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_genus_bray_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_ac_sum_stat["ac_genus_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_genus_jaccard_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_ac_sum_stat["ac_genus_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_genus_aitchison_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_ac_sum_stat["ac_genus_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_genus_kendall_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_species_bray_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_species_jaccard_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_species_aitchison_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_species_kendall_uniqueness_age_EstAsia_Studies$model$b)

rem_est_age_ac_sum_stat["ac_genus_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_genus_bray_uniqueness_age_Other_Studies$model$b)
rem_est_age_ac_sum_stat["ac_genus_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_genus_jaccard_uniqueness_age_Other_Studies$model$b)
rem_est_age_ac_sum_stat["ac_genus_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_genus_aitchison_uniqueness_age_Other_Studies$model$b)
rem_est_age_ac_sum_stat["ac_genus_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_genus_kendall_uniqueness_age_Other_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_species_bray_uniqueness_age_Other_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_species_jaccard_uniqueness_age_Other_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_species_aitchison_uniqueness_age_Other_Studies$model$b)
rem_est_age_ac_sum_stat["ac_species_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_species_kendall_uniqueness_age_Other_Studies$model$b)

rem_p_val_age_ac_sum_stat <- as.data.frame(matrix(NA,8,4))
rownames(rem_p_val_age_ac_sum_stat) <- c("ac_genus_bray_uniqueness","ac_genus_jaccard_uniqueness","ac_genus_aitchison_uniqueness","ac_genus_kendall_uniqueness","ac_species_bray_uniqueness","ac_species_jaccard_uniqueness","ac_species_aitchison_uniqueness","ac_species_kendall_uniqueness")
colnames(rem_p_val_age_ac_sum_stat) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")

rem_p_val_age_ac_sum_stat["ac_genus_bray_uniqueness","Overall"] <- as.numeric(combined_rem_ac_genus_bray_uniqueness_age_Overall$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_ac_genus_jaccard_uniqueness_age_Overall$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_ac_genus_aitchison_uniqueness_age_Overall$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_ac_genus_kendall_uniqueness_age_Overall$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_bray_uniqueness","Overall"] <- as.numeric(combined_rem_ac_species_bray_uniqueness_age_Overall$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_ac_species_jaccard_uniqueness_age_Overall$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_ac_species_aitchison_uniqueness_age_Overall$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_ac_species_kendall_uniqueness_age_Overall$model$pval)

rem_p_val_age_ac_sum_stat["ac_genus_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_genus_bray_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_genus_jaccard_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_genus_aitchison_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_genus_kendall_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_species_bray_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_species_jaccard_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_species_aitchison_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_species_kendall_uniqueness_age_EU_NA_Studies$model$pval)

rem_p_val_age_ac_sum_stat["ac_genus_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_genus_bray_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_genus_jaccard_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_genus_aitchison_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_genus_kendall_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_species_bray_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_species_jaccard_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_species_aitchison_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_species_kendall_uniqueness_age_EstAsia_Studies$model$pval)

rem_p_val_age_ac_sum_stat["ac_genus_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_genus_bray_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_genus_jaccard_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_genus_aitchison_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_genus_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_genus_kendall_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_species_bray_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_species_jaccard_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_species_aitchison_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_ac_sum_stat["ac_species_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_species_kendall_uniqueness_age_Other_Studies$model$pval)

rem_q_val_age_ac_sum_stat <- t(apply(rem_p_val_age_ac_sum_stat,1,function(x)(p.adjust(x,method="fdr"))))

rem_dir_age_ac_sum_stat <- as.data.frame(matrix(NA,8,4))
rownames(rem_dir_age_ac_sum_stat) <- c("ac_genus_bray_uniqueness","ac_genus_jaccard_uniqueness","ac_genus_aitchison_uniqueness","ac_genus_kendall_uniqueness","ac_species_bray_uniqueness","ac_species_jaccard_uniqueness","ac_species_aitchison_uniqueness","ac_species_kendall_uniqueness")
colnames(rem_dir_age_ac_sum_stat) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")
for(i in 1:nrow(rem_dir_age_ac_sum_stat))
{
	for(j in 1:ncol(rem_dir_age_ac_sum_stat))
	{
		rem_dir_age_ac_sum_stat[i,j] <- ifelse(rem_q_val_age_ac_sum_stat[i,j] <= 0.1,3*sign(rem_est_age_ac_sum_stat[i,j]),ifelse(rem_p_val_age_ac_sum_stat[i,j] <= 0.05,2*sign(rem_est_age_ac_sum_stat[i,j]),sign(rem_est_age_ac_sum_stat[i,j])))
	}
}

rem_est_age_sum_stat_pathway <- as.data.frame(matrix(NA,5,4))
rownames(rem_est_age_sum_stat_pathway) <- c("pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")
colnames(rem_est_age_sum_stat_pathway) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")

rem_est_age_sum_stat_pathway["pathway_shannon","Overall"] <- as.numeric(combined_rem_pathway_shannon_age_Overall$model$b)
rem_est_age_sum_stat_pathway["pathway_bray_uniqueness","Overall"] <- as.numeric(combined_rem_pathway_bray_uniqueness_age_Overall$model$b)
rem_est_age_sum_stat_pathway["pathway_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_pathway_jaccard_uniqueness_age_Overall$model$b)
rem_est_age_sum_stat_pathway["pathway_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_pathway_aitchison_uniqueness_age_Overall$model$b)
rem_est_age_sum_stat_pathway["pathway_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_pathway_kendall_uniqueness_age_Overall$model$b)

rem_est_age_sum_stat_pathway["pathway_shannon","EU_NA_Studies"] <- as.numeric(combined_rem_pathway_shannon_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_pathway_bray_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_pathway_jaccard_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_pathway_aitchison_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_pathway_kendall_uniqueness_age_EU_NA_Studies$model$b)

rem_est_age_sum_stat_pathway["pathway_shannon","EstAsia_Studies"] <- as.numeric(combined_rem_pathway_shannon_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_pathway_bray_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_pathway_jaccard_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_pathway_aitchison_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_pathway_kendall_uniqueness_age_EstAsia_Studies$model$b)

rem_est_age_sum_stat_pathway["pathway_shannon","Other_Studies"] <- as.numeric(combined_rem_pathway_shannon_age_Other_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_pathway_bray_uniqueness_age_Other_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_pathway_jaccard_uniqueness_age_Other_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_pathway_aitchison_uniqueness_age_Other_Studies$model$b)
rem_est_age_sum_stat_pathway["pathway_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_pathway_kendall_uniqueness_age_Other_Studies$model$b)

rem_p_val_age_sum_stat_pathway <- as.data.frame(matrix(NA,5,4))
rownames(rem_p_val_age_sum_stat_pathway) <- c("pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")
colnames(rem_p_val_age_sum_stat_pathway) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")

rem_p_val_age_sum_stat_pathway["pathway_shannon","Overall"] <- as.numeric(combined_rem_pathway_shannon_age_Overall$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_bray_uniqueness","Overall"] <- as.numeric(combined_rem_pathway_bray_uniqueness_age_Overall$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_pathway_jaccard_uniqueness_age_Overall$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_pathway_aitchison_uniqueness_age_Overall$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_pathway_kendall_uniqueness_age_Overall$model$pval)

rem_p_val_age_sum_stat_pathway["pathway_shannon","EU_NA_Studies"] <- as.numeric(combined_rem_pathway_shannon_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_pathway_bray_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_pathway_jaccard_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_pathway_aitchison_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_pathway_kendall_uniqueness_age_EU_NA_Studies$model$pval)

rem_p_val_age_sum_stat_pathway["pathway_shannon","EstAsia_Studies"] <- as.numeric(combined_rem_pathway_shannon_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_pathway_bray_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_pathway_jaccard_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_pathway_aitchison_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_pathway_kendall_uniqueness_age_EstAsia_Studies$model$pval)

rem_p_val_age_sum_stat_pathway["pathway_shannon","Other_Studies"] <- as.numeric(combined_rem_pathway_shannon_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_pathway_bray_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_pathway_jaccard_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_pathway_aitchison_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_sum_stat_pathway["pathway_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_pathway_kendall_uniqueness_age_Other_Studies$model$pval)

rem_q_val_age_sum_stat_pathway <- t(apply(rem_p_val_age_sum_stat_pathway,1,function(x)(p.adjust(x,method="fdr"))))

rem_dir_age_sum_stat_pathway <- as.data.frame(matrix(NA,5,4))
rownames(rem_dir_age_sum_stat_pathway) <- c("pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")
colnames(rem_dir_age_sum_stat_pathway) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")
for(i in 1:nrow(rem_dir_age_sum_stat_pathway))
{
	for(j in 1:ncol(rem_dir_age_sum_stat_pathway))
	{
		rem_dir_age_sum_stat_pathway[i,j] <- ifelse(rem_q_val_age_sum_stat_pathway[i,j] <= 0.1,3*sign(rem_est_age_sum_stat_pathway[i,j]),ifelse(rem_p_val_age_sum_stat_pathway[i,j] <= 0.05,2*sign(rem_est_age_sum_stat_pathway[i,j]),sign(rem_est_age_sum_stat_pathway[i,j])))
	}
}

rem_est_age_ac_sum_stat_pathway <- as.data.frame(matrix(NA,4,4))
rownames(rem_est_age_ac_sum_stat_pathway) <- c("pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")
colnames(rem_est_age_ac_sum_stat_pathway) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")

rem_est_age_ac_sum_stat_pathway["pathway_bray_uniqueness","Overall"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_Overall$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_ac_pathway_jaccard_uniqueness_age_Overall$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_Overall$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_ac_pathway_kendall_uniqueness_age_Overall$model$b)

rem_est_age_ac_sum_stat_pathway["pathway_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_pathway_jaccard_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_EU_NA_Studies$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_pathway_kendall_uniqueness_age_EU_NA_Studies$model$b)

rem_est_age_ac_sum_stat_pathway["pathway_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_pathway_jaccard_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_EstAsia_Studies$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_pathway_kendall_uniqueness_age_EstAsia_Studies$model$b)

rem_est_age_ac_sum_stat_pathway["pathway_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_Other_Studies$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_pathway_jaccard_uniqueness_age_Other_Studies$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_Other_Studies$model$b)
rem_est_age_ac_sum_stat_pathway["pathway_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_pathway_kendall_uniqueness_age_Other_Studies$model$b)

rem_p_val_age_ac_sum_stat_pathway <- as.data.frame(matrix(NA,4,4))
rownames(rem_p_val_age_ac_sum_stat_pathway) <- c("pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")
colnames(rem_p_val_age_ac_sum_stat_pathway) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")

rem_p_val_age_ac_sum_stat_pathway["pathway_bray_uniqueness","Overall"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_Overall$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_jaccard_uniqueness","Overall"] <- as.numeric(combined_rem_ac_pathway_jaccard_uniqueness_age_Overall$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_aitchison_uniqueness","Overall"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_Overall$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_kendall_uniqueness","Overall"] <- as.numeric(combined_rem_ac_pathway_kendall_uniqueness_age_Overall$model$pval)

rem_p_val_age_ac_sum_stat_pathway["pathway_bray_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_jaccard_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_pathway_jaccard_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_aitchison_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_EU_NA_Studies$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_kendall_uniqueness","EU_NA_Studies"] <- as.numeric(combined_rem_ac_pathway_kendall_uniqueness_age_EU_NA_Studies$model$pval)

rem_p_val_age_ac_sum_stat_pathway["pathway_bray_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_jaccard_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_pathway_jaccard_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_aitchison_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_EstAsia_Studies$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_kendall_uniqueness","EstAsia_Studies"] <- as.numeric(combined_rem_ac_pathway_kendall_uniqueness_age_EstAsia_Studies$model$pval)

rem_p_val_age_ac_sum_stat_pathway["pathway_bray_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_jaccard_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_pathway_jaccard_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_aitchison_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_pathway_bray_uniqueness_age_Other_Studies$model$pval)
rem_p_val_age_ac_sum_stat_pathway["pathway_kendall_uniqueness","Other_Studies"] <- as.numeric(combined_rem_ac_pathway_kendall_uniqueness_age_Other_Studies$model$pval)

rem_q_val_age_ac_sum_stat_pathway <- t(apply(rem_p_val_age_ac_sum_stat_pathway,1,function(x)(p.adjust(x,method="fdr"))))

rem_dir_age_ac_sum_stat_pathway <- as.data.frame(matrix(NA,4,4))
rownames(rem_dir_age_ac_sum_stat_pathway) <- c("pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")
colnames(rem_dir_age_ac_sum_stat_pathway) <- c("Overall","EU_NA_Studies","EstAsia_Studies","Other_Studies")
for(i in 1:nrow(rem_dir_age_ac_sum_stat_pathway))
{
	for(j in 1:ncol(rem_dir_age_sum_stat_pathway))
	{
		rem_dir_age_ac_sum_stat_pathway[i,j] <- ifelse(rem_q_val_age_ac_sum_stat_pathway[i,j] <= 0.1,3*sign(rem_est_age_ac_sum_stat_pathway[i,j]),ifelse(rem_p_val_age_sum_stat_pathway[i,j] <= 0.05,2*sign(rem_est_age_ac_sum_stat_pathway[i,j]),sign(rem_est_age_ac_sum_stat_pathway[i,j])))
	}
}

#mat <- rem_dir_age_sum_stat
#mat <- apply(mat,2,function(x)(ifelse(x<0,x+0.1,x)))
#heatmap.2(mat,density="none",trace="none",col=c("skyblue4","skyblue3","powderblue","mistyrose1","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),key=FALSE,margins=c(25,5),lhei=c(0.1,5),lwid=c(1,5))

#mat <- rem_dir_age_ac_sum_stat
#mat <- apply(mat,2,function(x)(ifelse(x<0,x+0.1,x)))
#heatmap.2(mat,density="none",trace="none",col=c("skyblue4","skyblue3","powderblue","mistyrose1","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),key=FALSE,margins=c(25,5),lhei=c(0.1,5),lwid=c(1,5))

#mat <- rem_dir_age_sum_stat_pathway
#mat <- apply(mat,2,function(x)(ifelse(x<0,x+0.1,x)))
#heatmap.2(mat,density="none",trace="none",col=c("skyblue4","skyblue3","powderblue","mistyrose1","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="purple",colsep=1:ncol(mat),rowsep=1:nrow(mat),key=FALSE,margins=c(25,5),lhei=c(0.1,5),lwid=c(1,5))
