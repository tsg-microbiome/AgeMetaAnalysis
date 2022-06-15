# The sub-pipeline collates the Metacyc-annotated pathway and metadata tables for the Irish Shotgun Cohorts (ISC) data repository.
# Computes Shannon Diversity and the four different taxonomic uniqueness at the Pathway levels across all gut microbiomes
# The output is a single data frame with the select metadata and the diversity and uniqueness measures for all gut microbiomes
# It subsequently performs the following functions:
# Associations between Shannon Diversity and the four different measures of Uniqueness at the Pathway level.
# Associations between pathway-level gut microbiome beta-diversity (computed using the four different distance measures; the same ones corresponding to the four measures of uniqueness) and age.
# Associations between the different measures of uniqueness (four measures of uniqueness at the pathway level) and Shannon Diversity with age.


library(vegan)
library(ade4)
library(MASS)
library(sfsmisc)
library(pcaPP)
library(robumeta)
library(metafor)
library(dplyr)
library(effsize)

range_scale=function(x)
{
	y <- (x-min(x))/(max(x)-min(x));
	return(y);
}

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\EM_Pathabundance.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\ExM_Pathabundance.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\IBS4DPathabundance.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\ISC_Metadata.RData")

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\isc_analysis_2021_Revision.RData")

IBS4D_Pathabundance <- t(apply(IBS4D_pathabundance[,-295],1,function(x)(x/sum(x))))

temp0 <- merge(t(EM_Pathabundance),t(ExM_Pathabundance),by="row.names",all=TRUE)[,-1]
rownames(temp0) <- merge(t(EM_Pathabundance),t(ExM_Pathabundance),by="row.names",all=TRUE)[,1]
temp0 <- apply(temp0,1,function(x)(ifelse(is.na(x),0,x)))

temp1 <- merge(t(temp0),t(IBS4D_Pathabundance),by="row.names",all=TRUE)[,-1]
rownames(temp1) <- merge(t(temp0),t(IBS4D_Pathabundance),by="row.names",all=TRUE)[,1]
temp1 <- apply(temp1,1,function(x)(ifelse(is.na(x),0,x)))

isc_select_age_final_pathway <- temp1

rownames(isc_select_metadata) <- sub("_S","_S0",rownames(isc_select_metadata))
rownames(isc_select_metadata) <- sub("^S","S0",rownames(isc_select_metadata))

common_rows <- intersect(rownames(isc_select_age_final_pathway),rownames(isc_select_metadata))

control_rows <- rownames(isc_select_metadata[isc_select_metadata$study_condition == "control",])

common_rows <- intersect(common_rows,control_rows)

isc_select_age_final_pathway_metadata <- isc_select_metadata[common_rows,]
isc_select_age_final_pathway <- isc_select_age_final_pathway[common_rows,]

df_isc_pathway_diversity_uniqueness <- as.data.frame(isc_select_age_final_pathway_metadata)

df_isc_pathway_diversity_uniqueness$pathway_shannon <- diversity(isc_select_age_final_pathway)

print("Pathway:Bray")
dist_pathway_bray <- as.matrix(vegdist(isc_select_age_final_pathway,method="bray"))
diag(dist_pathway_bray) <- NA

print("Pathway:Jaccard")
dist_pathway_jaccard <- as.matrix(vegdist(isc_select_age_final_pathway,method="jaccard"))
diag(dist_pathway_jaccard) <- NA

isc_select_age_final_pathway_clr <- as.matrix(clr(isc_select_age_final_pathway+0.00001))
isc_select_age_final_pathway_clr <- as.data.frame(t(apply(isc_select_age_final_pathway_clr,1,function(x)(x-min(x)))))
print("Pathway:Manhattan")
dist_pathway_aitchison <- as.matrix(vegdist(isc_select_age_final_pathway_clr,method="euclidean"))
diag(dist_pathway_aitchison) <- NA

print("Pathway:Kendall")
dist_pathway_kendall <- as.matrix(1-cor.fk(t(isc_select_age_final_pathway_clr))/2)
diag(dist_pathway_kendall) <- NA

df_isc_pathway_diversity_uniqueness$pathway_bray_uniqueness <- NA
df_isc_pathway_diversity_uniqueness$pathway_jaccard_uniqueness <- NA
df_isc_pathway_diversity_uniqueness$pathway_aitchison_uniqueness <- NA
df_isc_pathway_diversity_uniqueness$pathway_kendall_uniqueness <- NA

study_names <- unique(isc_select_age_final_pathway_metadata$study_name)

for(i in 1:length(study_names))
{
	study_name <- study_names[i]
	study_samples <- rownames(isc_select_age_final_pathway_metadata[isc_select_age_final_pathway_metadata$study_name == study_name,])
	
	study_samples <- intersect(common_rows,study_samples)
	
	df_isc_pathway_diversity_uniqueness[study_samples,"pathway_bray_uniqueness"] <- apply(dist_pathway_bray[study_samples,study_samples],1,function(x)(min(x[!is.na(x)])))

	df_isc_pathway_diversity_uniqueness[study_samples,"pathway_jaccard_uniqueness"] <- apply(dist_pathway_jaccard[study_samples,study_samples],1,function(x)(min(x[!is.na(x)])))
	
	df_isc_pathway_diversity_uniqueness[study_samples,"pathway_aitchison_uniqueness"] <- apply(dist_pathway_aitchison[study_samples,study_samples],1,function(x)(min(x[!is.na(x)])))
	
	df_isc_pathway_diversity_uniqueness[study_samples,"pathway_kendall_uniqueness"] <- apply(dist_pathway_kendall[study_samples,study_samples],1,function(x)(min(x[!is.na(x)])))
}

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
		print(group)
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

#### START OF Step A ###
print("Step A")

isc_rlm_est_pathway_diversity_uniqueness <- as.data.frame(matrix(NA,1,4))
rownames(isc_rlm_est_pathway_diversity_uniqueness) <- "ISC"
colnames(isc_rlm_est_pathway_diversity_uniqueness) <- c("pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

isc_rlm_p_val_pathway_diversity_uniqueness <- as.data.frame(matrix(NA,1,4))
rownames(isc_rlm_p_val_pathway_diversity_uniqueness) <- "ISC"
colnames(isc_rlm_p_val_pathway_diversity_uniqueness) <- c("pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

temp_rlm <- rlm(pathway_bray_uniqueness~pathway_shannon,df_isc_pathway_diversity_uniqueness)
isc_rlm_est_pathway_diversity_uniqueness[1,1] <- as.numeric(temp_rlm$coefficients[2])
isc_rlm_p_val_pathway_diversity_uniqueness[1,1] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(pathway_jaccard_uniqueness~pathway_shannon,df_isc_pathway_diversity_uniqueness)
isc_rlm_est_pathway_diversity_uniqueness[1,2] <- as.numeric(temp_rlm$coefficients[2])
isc_rlm_p_val_pathway_diversity_uniqueness[1,2] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(pathway_aitchison_uniqueness~pathway_shannon,df_isc_pathway_diversity_uniqueness)
isc_rlm_est_pathway_diversity_uniqueness[1,3] <- as.numeric(temp_rlm$coefficients[2])
isc_rlm_p_val_pathway_diversity_uniqueness[1,3] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(pathway_kendall_uniqueness~pathway_shannon,df_isc_pathway_diversity_uniqueness)
isc_rlm_est_pathway_diversity_uniqueness[1,4] <- as.numeric(temp_rlm$coefficients[2])
isc_rlm_p_val_pathway_diversity_uniqueness[1,4] <- f.robftest(temp_rlm)$p.value

isc_rlm_q_val_pathway_diversity_uniqueness <- t(apply(isc_rlm_p_val_pathway_diversity_uniqueness,1,function(x)(p.adjust(x,method="fdr"))))

isc_rlm_dir_pathway_diversity_uniqueness <- as.data.frame(matrix(NA,1,4))
rownames(isc_rlm_dir_pathway_diversity_uniqueness) <- "ISC"
colnames(isc_rlm_dir_pathway_diversity_uniqueness) <- c("pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

for(j in 1:ncol(isc_rlm_dir_pathway_diversity_uniqueness))
{
	isc_rlm_dir_pathway_diversity_uniqueness[1,j] <- ifelse(isc_rlm_q_val_pathway_diversity_uniqueness[1,j]<=0.05,3*sign(isc_rlm_est_pathway_diversity_uniqueness[1,j]),ifelse(isc_rlm_q_val_pathway_diversity_uniqueness[1,j]<=0.1,2*sign(isc_rlm_est_pathway_diversity_uniqueness[1,j]),sign(isc_rlm_est_pathway_diversity_uniqueness[1,j])))
}
#### END OF Step A ###

#### START OF Step B ###
print("Step B")
isc_pathway_est_age_beta <- as.data.frame(matrix(NA,1,5))
rownames(isc_pathway_est_age_beta) <- "ISC"
colnames(isc_pathway_est_age_beta) <- c("pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

isc_pathway_p_val_age_beta <- as.data.frame(matrix(NA,1,5))
rownames(isc_pathway_p_val_age_beta) <- "ISC"
colnames(isc_pathway_p_val_age_beta) <- c("pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

isc_pathway_est_age_beta[1,1] <- 0
isc_pathway_p_val_age_beta[1,1] <- 1

temp_adonis <- adonis(as.dist(dist_pathway_bray[rownames(df_isc_pathway_diversity_uniqueness),rownames(df_isc_pathway_diversity_uniqueness)])~df_isc_pathway_diversity_uniqueness[rownames(df_isc_pathway_diversity_uniqueness),"age"])
isc_pathway_est_age_beta[1,2] <- temp_adonis$aov.tab[1,5]
isc_pathway_p_val_age_beta[1,2] <- temp_adonis$aov.tab[1,6]

temp_adonis <- adonis(as.dist(dist_pathway_jaccard[rownames(df_isc_pathway_diversity_uniqueness),rownames(df_isc_pathway_diversity_uniqueness)])~df_isc_pathway_diversity_uniqueness[rownames(df_isc_pathway_diversity_uniqueness),"age"])
isc_pathway_est_age_beta[1,3] <- temp_adonis$aov.tab[1,5]
isc_pathway_p_val_age_beta[1,3] <- temp_adonis$aov.tab[1,6]

temp_adonis <- adonis(as.dist(dist_pathway_aitchison[rownames(df_isc_pathway_diversity_uniqueness),rownames(df_isc_pathway_diversity_uniqueness)])~df_isc_pathway_diversity_uniqueness[rownames(df_isc_pathway_diversity_uniqueness),"age"])
isc_pathway_est_age_beta[1,4] <- temp_adonis$aov.tab[1,5]
isc_pathway_p_val_age_beta[1,4] <- temp_adonis$aov.tab[1,6]

temp_adonis <- adonis(as.dist(dist_pathway_kendall[rownames(df_isc_pathway_diversity_uniqueness),rownames(df_isc_pathway_diversity_uniqueness)])~df_isc_pathway_diversity_uniqueness[rownames(df_isc_pathway_diversity_uniqueness),"age"])
isc_pathway_est_age_beta[1,5] <- temp_adonis$aov.tab[1,5]
isc_pathway_p_val_age_beta[1,5] <- temp_adonis$aov.tab[1,6]

#### END OF Step B ###

### START OF Step C ###
print("StepC: Uniqueness/Diversity v/s Age")
isc_rlm_est_pathway_age_sum_stat <- as.data.frame(matrix(NA,1,5))
rownames(isc_rlm_est_pathway_age_sum_stat) <- "ISC"
colnames(isc_rlm_est_pathway_age_sum_stat) <- c("pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

isc_rlm_p_val_pathway_age_sum_stat <- as.data.frame(matrix(NA,1,5))
rownames(isc_rlm_p_val_pathway_age_sum_stat) <- "ISC"
colnames(isc_rlm_p_val_pathway_age_sum_stat) <- c("pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

temp_rlm <- rlm(pathway_shannon~age,data=df_isc_pathway_diversity_uniqueness,psi = psi.bisquare)
isc_rlm_est_pathway_age_sum_stat[1,1] <- as.numeric(temp_rlm$coefficients[2])
isc_rlm_p_val_pathway_age_sum_stat[1,1] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(pathway_bray_uniqueness~age,data=df_isc_pathway_diversity_uniqueness,psi = psi.bisquare)
isc_rlm_est_pathway_age_sum_stat[1,2] <- as.numeric(temp_rlm$coefficients[2])
isc_rlm_p_val_pathway_age_sum_stat[1,2] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(pathway_jaccard_uniqueness~age,data=df_isc_pathway_diversity_uniqueness,psi = psi.bisquare)
isc_rlm_est_pathway_age_sum_stat[1,3] <- as.numeric(temp_rlm$coefficients[2])
isc_rlm_p_val_pathway_age_sum_stat[1,3] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(pathway_aitchison_uniqueness~age,data=df_isc_pathway_diversity_uniqueness,psi = psi.bisquare)
isc_rlm_est_pathway_age_sum_stat[1,4] <- as.numeric(temp_rlm$coefficients[2])
isc_rlm_p_val_pathway_age_sum_stat[1,4] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(pathway_kendall_uniqueness~age,data=df_isc_pathway_diversity_uniqueness,psi = psi.bisquare)
isc_rlm_est_pathway_age_sum_stat[1,5] <- as.numeric(temp_rlm$coefficients[2])
isc_rlm_p_val_pathway_age_sum_stat[1,5] <- f.robftest(temp_rlm)$p.value

isc_rlm_q_val_pathway_age_sum_stat <- t(apply(isc_rlm_p_val_pathway_age_sum_stat,1,function(x)(p.adjust(x,method="fdr"))))

isc_rlm_dir_pathway_age_sum_stat <- as.data.frame(matrix(NA,1,5))
rownames(isc_rlm_dir_pathway_age_sum_stat) <- "ISC"
colnames(isc_rlm_dir_pathway_age_sum_stat) <- c("pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

for(j in 1:ncol(isc_rlm_dir_pathway_age_sum_stat))
{
	isc_rlm_dir_pathway_age_sum_stat[1,j] <- ifelse(isc_rlm_q_val_pathway_age_sum_stat[1,j]<=0.10,3*sign(isc_rlm_est_pathway_age_sum_stat[1,j]),ifelse(isc_rlm_p_val_pathway_age_sum_stat[1,j]<=0.05,2*sign(isc_rlm_est_pathway_age_sum_stat[1,j]),1*sign(isc_rlm_est_pathway_age_sum_stat[1,j])))
}

#### END OF Step C ###

### START OF Step D ###
print("StepD: Alpha-Corrected Uniqueness v/s Age")
isc_rlm_est_pathway_age_ac_sum_stat <- as.data.frame(matrix(NA,1,5))
rownames(isc_rlm_est_pathway_age_ac_sum_stat) <- "ISC"
colnames(isc_rlm_est_pathway_age_ac_sum_stat) <- c("pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

isc_rlm_p_val_pathway_age_ac_sum_stat <- as.data.frame(matrix(NA,1,5))
rownames(isc_rlm_p_val_pathway_age_ac_sum_stat) <- "ISC"
colnames(isc_rlm_p_val_pathway_age_ac_sum_stat) <- c("pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

temp_rlm <- rlm(pathway_shannon~age,data=df_isc_pathway_diversity_uniqueness,psi = psi.bisquare)
isc_rlm_est_pathway_age_ac_sum_stat[1,1] <- as.numeric(temp_rlm$coefficients[2])
isc_rlm_p_val_pathway_age_ac_sum_stat[1,1] <- f.robftest(temp_rlm)$p.value

temp_rlm <- rlm(pathway_bray_uniqueness~pathway_shannon+age,df_isc_pathway_diversity_uniqueness,psi = psi.bisquare)
isc_rlm_est_pathway_age_ac_sum_stat[1,2] <- as.numeric(temp_rlm$coefficients[3])
isc_rlm_p_val_pathway_age_ac_sum_stat[1,2] <- (f.robftest(temp_rlm,var="age"))$p.value

temp_rlm <- rlm(pathway_jaccard_uniqueness~pathway_shannon+age,df_isc_pathway_diversity_uniqueness,psi = psi.bisquare)
isc_rlm_est_pathway_age_ac_sum_stat[1,3] <- as.numeric(temp_rlm$coefficients[3])
isc_rlm_p_val_pathway_age_ac_sum_stat[1,3] <- (f.robftest(temp_rlm,var="age"))$p.value

temp_rlm <- rlm(pathway_aitchison_uniqueness~pathway_shannon+age,df_isc_pathway_diversity_uniqueness,psi = psi.bisquare)
isc_rlm_est_pathway_age_ac_sum_stat[1,4] <- as.numeric(temp_rlm$coefficients[3])
isc_rlm_p_val_pathway_age_ac_sum_stat[1,4] <- (f.robftest(temp_rlm,var="age"))$p.value

temp_rlm <- rlm(pathway_kendall_uniqueness~pathway_shannon+age,df_isc_pathway_diversity_uniqueness,psi = psi.bisquare)
isc_rlm_est_pathway_age_ac_sum_stat[1,5] <- as.numeric(temp_rlm$coefficients[3])
isc_rlm_p_val_pathway_age_ac_sum_stat[1,5] <- (f.robftest(temp_rlm,var="age"))$p.value

isc_rlm_q_val_pathway_age_ac_sum_stat <- t(apply(isc_rlm_p_val_pathway_age_ac_sum_stat,1,function(x)(p.adjust(x,method="fdr"))))

isc_rlm_dir_pathway_age_ac_sum_stat <- as.data.frame(matrix(NA,1,5))
rownames(isc_rlm_dir_pathway_age_ac_sum_stat) <- "ISC"
colnames(isc_rlm_dir_pathway_age_ac_sum_stat) <- c("pathway_shannon","pathway_bray_uniqueness","pathway_jaccard_uniqueness","pathway_aitchison_uniqueness","pathway_kendall_uniqueness")

for(j in 1:ncol(isc_rlm_dir_pathway_age_ac_sum_stat))
{
	isc_rlm_dir_pathway_age_ac_sum_stat[1,j] <- ifelse(isc_rlm_q_val_pathway_age_ac_sum_stat[1,j]<=0.10,3*sign(isc_rlm_est_pathway_age_ac_sum_stat[1,j]),ifelse(isc_rlm_p_val_pathway_age_ac_sum_stat[1,j]<=0.05,2*sign(isc_rlm_est_pathway_age_ac_sum_stat[1,j]),1*sign(isc_rlm_est_pathway_age_ac_sum_stat[1,j])))
}

save(list=(c(ls(pattern="df_isc_"),ls(pattern="est"),ls(pattern="p_val"),ls(pattern="q_val"),ls(pattern="dir"))),file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_stage1_pathway_results.RData")

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\isc_analysis_2021_Revision.RData")
save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ISC\\isc_analysis_2021_Revision.RData")
rm(list=ls())



