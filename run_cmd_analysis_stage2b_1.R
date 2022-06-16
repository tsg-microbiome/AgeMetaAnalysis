load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_disease_analysis.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")


library(robumeta)
library(metafor)
library(dplyr)
library(effsize)

print("Wilcox")
wilcox_batch <- function(x,y)
{
	p_array <- NULL;
	type_array <- NULL;
	mean1_array <- NULL;
	mean2_array <- NULL;
	x <- x[abs(rowSums(x,na.rm=TRUE)) > 0,];
	y <- y[abs(rowSums(y,na.rm=TRUE)) > 0,];
	z <- intersect(rownames(x),rownames(y));
	for(i in 1:length(z))
	{
		p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]))$p.value;
		type_array[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
		mean1_array[i] <- mean(as.numeric(x[z[i],]),na.rm=TRUE);
		mean2_array[i] <- mean(as.numeric(y[z[i],]),na.rm=TRUE);
		i <- i + 1;
	}
	out <- as.data.frame(cbind(p_array,type_array,p.adjust(p_array,method="fdr"),mean1_array,mean2_array));
	rownames(out) <- z;
	out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
	return(t(out));
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
		#print(group)
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

#cmd3_select_disease_final_metadata$study_condition = ifelse(cmd3_select_disease_final_metadata$study_condition %in% c("T2D","IGT"),"T2D_IGT",cmd3_select_disease_final_metadata$study_condition)

studies_previously_considered <- c("FengQ_2015","NielsenHB_2014","HMP_2019_ibdmdb","QinJ_2012","KarlssonFH_2013","QinN_2014","VogtmannE_2016","ZellerG_2014")

cmd3_select_disease_final_metadata$study_with_condition <- paste0(cmd3_select_disease_final_metadata$study_condition,":",cmd3_select_disease_final_metadata$study_name)

new_studies_to_check <- setdiff(rownames(cmd3_select_disease_studies_details),studies_previously_considered)

diseases_to_check <- setdiff(unique(cmd3_select_disease_final_metadata[cmd3_select_disease_final_metadata$study_name %in% new_studies_to_check,"study_condition"]),c("control","IGT","carcinoma_surgery_history"))

conditions_to_check <- unique(cmd3_select_disease_final_metadata[(cmd3_select_disease_final_metadata$study_condition %in% diseases_to_check)&(cmd3_select_disease_final_metadata$study_name %in% new_studies_to_check),"study_with_condition"])

conditions_to_check <- conditions_to_check[order(conditions_to_check)]

temp_data <- as.data.frame(cbind(cmd3_select_disease_final_species,cmd3_select_disease_final_metadata[,c("study_condition","study_with_condition","study_name","age")]))
colnames(temp_data)[1638] <- "study_condition"
colnames(temp_data)[1639] <- "study_with_condition"
colnames(temp_data)[1640] <- "study_name"
colnames(temp_data)[1641] <- "age"

cmd3_select_df_species <- temp_data[temp_data$study_name %in% new_studies_to_check,]
## Doing name change to revert to old taxonomy scheme
## Earlier Clostridium_nexile is now Tyzzerella_nexilis reverted back
colnames(cmd3_select_df_species)[267] <- "Clostridium_nexile"
## Earlier Clostridium_hathewayi now Hungatella_hathewayi reverted back
colnames(cmd3_select_df_species)[261] <- "Clostridium_hathewayi"
## Earlier Earlier Streptococcus_anginosus is now Streptococcus_anginosus_group reverted back
colnames(cmd3_select_df_species)[389] <- "Streptococcus_anginosus"
## Earlier Enterobacter_cloacae is now Enterobacter_cloacae_group reverted back
colnames(cmd3_select_df_species)[25] <- "Enterobacter_cloacae"
## Earlier Clostridium_ramosum is now Erysipelatoclostridium_ramosum reverted back
colnames(cmd3_select_df_species)[42] <- "Clostridium_ramosum"

G_Markers <- intersect(c("Bacteroides_fragilis","Actinomyces_odontolyticus","Bifidobacterium_dentium","Clostridium_clostridioforme","Clostridium_nexile","Clostridium_ramosum","Dialister_invisus","Eggerthella_lenta","Escherichia_coli","Fusobacterium_nucleatum","Granulicatella_adiacens","Lactococcus_lactis","Prevotella_copri","Rothia_mucilaginosa","Ruminococcus_gnavus","Streptococcus_anginosus","Streptococcus_infantis","Streptococcus_salivarius","Streptococcus_sanguinis","Streptococcus_vestibularis","Veillonella_atypica","Clostridium_asparagiforme","Clostridiales_bacterium_1_7_47FAA","Clostridium_bolteae","Clostridium_citroniae","Clostridium_hathewayi","Clostridium_symbiosum","Streptococcus_parasanguinis","Solobacterium_moorei","Ruminococcus_torques","Streptococcus_mitis","Klebsiella_pneumoniae","Streptococcus_australis","Streptococcus_gordonii","Enterobacter_cloacae"),colnames(cmd3_select_df_species))

L_Markers <- intersect(c("Eubacterium_hallii","Dorea_longicatena","Coprococcus_comes","Coprococcus_catus","Butyrivibrio_crossotus","Bacteroides_uniformis","Alistipes_shahii","Alistipes_indistinctus","Roseburia_hominis","Pseudoflavonifractor_capillosus","Eubacterium_siraeum","Eubacterium_rectale","Coprobacter_fastidiosus","Bifidobacterium_longum","Bifidobacterium_animalis","Barnesiella_intestinihominis","Bacteroides_xylanisolvens","Alistipes_senegalensis","Alistipes_putredinis","Alistipes_onderdonkii","Akkermansia_muciniphila","Faecalibacterium_prausnitzii","Roseburia_inulinivorans","Roseburia_intestinalis"),colnames(cmd3_select_df_species))

cmd3_select_df_species$G_Markers_combined <- rowMeans(apply(cmd3_select_df_species[,G_Markers],2,rank_scale))
cmd3_select_df_species$L_Markers_combined <- rowMeans(apply(cmd3_select_df_species[,L_Markers],2,rank_scale))

SpeciesDetection <- compute_detection(cmd3_select_df_species,colnames(cmd3_select_df_species)[1:1637],"study_name",new_studies_to_check)

select_species <- names(which(apply(SpeciesDetection,1,function(x)(length(x[x>0])))>=9))

disease_list <- diseases_to_check

disease_patient_numbers <- matrix(NA,length(disease_list),4)
rownames(disease_patient_numbers) <- c(disease_list)
colnames(disease_patient_numbers) <- c("number","minimum_age","maximum_age","study_names")

for(i in 1:length(disease_list))
{
	disease_name <- disease_list[i]
	study_list <- unique(cmd3_select_df_species[cmd3_select_df_species$study_condition == disease_name,"study_name"])
	assign(paste0("study_list_",disease_name),study_list)
	assign(paste0("patient_list_",disease_name),rownames(cmd3_select_df_species[(cmd3_select_df_species$study_name %in% study_list)&(cmd3_select_df_species$study_condition == disease_name),]))
	assign(paste0("control_list_",disease_name),rownames(cmd3_select_df_species[(cmd3_select_df_species$study_name %in% study_list)&(cmd3_select_df_species$study_condition == "control"),]))
	disease_patient_numbers[disease_name,1] <- nrow(cmd3_select_df_species[which(!is.na(cmd3_select_df_species[,"study_name"]))&(cmd3_select_df_species[,"study_name"] %in% study_list),])
	disease_patient_numbers[disease_name,2] <- min(cmd3_select_df_species[which(!is.na(cmd3_select_df_species[,"study_name"]))&(cmd3_select_df_species[,"study_name"] %in% study_list),"age"])
	disease_patient_numbers[disease_name,3] <- max(cmd3_select_df_species[which(!is.na(cmd3_select_df_species[,"study_name"]))&(cmd3_select_df_species[,"study_name"] %in% study_list),"age"])
	disease_patient_numbers[disease_name,4] <- paste0(unique(cmd3_select_df_species[cmd3_select_df_species$study_condition == disease_name,"study_name"]),collapse=",")
}

cmd3_rlm_est_species_disease <- as.data.frame(matrix(0,length(select_species),length(disease_list)))
rownames(cmd3_rlm_est_species_disease) <- select_species
colnames(cmd3_rlm_est_species_disease) <- disease_list

cmd3_rlm_pval_species_disease <- as.data.frame(matrix(1,length(select_species),length(disease_list)))
rownames(cmd3_rlm_pval_species_disease) <- select_species
colnames(cmd3_rlm_pval_species_disease) <- disease_list

for(i in 1:length(select_species))
{
	species_name <- select_species[i]
	for(j in 1:length(disease_list))
	{
		#print(paste0(species_name,",",disease_name))
		disease_name <- disease_list[j]
		control_list <- get(paste0("control_list_",disease_name))
		diseased_list <- get(paste0("patient_list_",disease_name))
		vec_control <- cmd3_select_df_species[control_list,species_name]
		vec_disease <- cmd3_select_df_species[diseased_list,species_name]
		total_vec <- c(vec_disease,vec_control)
		if(length(total_vec[total_vec > 0])>0)
		{
			temp_rlm <- rlm(cmd3_select_df_species[c(diseased_list,control_list),species_name]~as.factor(ifelse(c(diseased_list,control_list) %in% control_list,0,1)))
			summary_temp_rlm <- summary(temp_rlm)
			if((!is.nan(summary_temp_rlm$coefficients[2,3]))&(summary_temp_rlm$coefficients[2,1] != 0)&(abs(summary_temp_rlm$coefficients[2,3]) <= 10000))
			{
				cmd3_rlm_est_species_disease[species_name,disease_name] <- summary_temp_rlm$coefficients[2,3]
				cmd3_rlm_pval_species_disease[species_name,disease_name] <- wilcox.test(vec_disease,vec_control)$p.value
			}
			else
			{
				cmd3_rlm_est_species_disease[species_name,disease_name] <- 0
				cmd3_rlm_pval_species_disease[species_name,disease_name] <- 1
			}
		}
	}
}
cmd3_rlm_qval_species_disease <- t(apply(cmd3_rlm_pval_species_disease,1,function(x)(p.adjust(x,method="fdr"))))
cmd3_rlm_est_select_species_disease <- cmd3_rlm_est_species_disease[c(G_Markers,L_Markers),]
cmd3_rlm_pval_select_species_disease <- cmd3_rlm_pval_species_disease[c(G_Markers,L_Markers),]
cmd3_rlm_qval_select_species_disease <- apply(cmd3_rlm_pval_select_species_disease,2,function(x)(p.adjust(x,method="fdr")))

select_markers <- c(G_Markers,L_Markers)

cmd3_rlm_dir_select_species_disease <- as.data.frame(matrix(NA,length(select_markers),length(disease_list)))
rownames(cmd3_rlm_dir_select_species_disease) <- select_markers
colnames(cmd3_rlm_dir_select_species_disease) <- disease_list

for(i in 1:length(select_markers))
{
	for(j in 1:length(disease_list))
	{
		cmd3_rlm_dir_select_species_disease[i,j] <- ifelse(cmd3_rlm_qval_select_species_disease[i,j] <= 0.1,3*sign(cmd3_rlm_est_select_species_disease[i,j]),ifelse(cmd3_rlm_pval_select_species_disease[i,j] <= 0.05,2*sign(cmd3_rlm_est_select_species_disease[i,j]),sign(cmd3_rlm_est_select_species_disease[i,j])))
	}
}

AssociationDirection <- as.data.frame(cbind(apply(cmd3_rlm_dir_select_species_disease,1,function(x)(length(x[x<= -2]))),apply(cmd3_rlm_dir_select_species_disease,1,function(x)(length(x[x>=2])))))

G_Markers_ordered <- rownames(AssociationDirection[intersect(rownames(AssociationDirection[order(AssociationDirection[,1]-AssociationDirection[,2]),]),G_Markers),])
L_Markers_ordered <- rownames(AssociationDirection[intersect(rownames(AssociationDirection[order(AssociationDirection[,1]-AssociationDirection[,2]),]),L_Markers),])

Marker_List_Sorted <- c(G_Markers_ordered,L_Markers_ordered)

mat_disease_association <- apply(cmd3_rlm_dir_select_species_disease[Marker_List_Sorted,],2,function(x)(ifelse(x<0,x-0.1,x)))

cmd3_rlm_dir_species_disease <- as.data.frame(matrix(NA,length(select_species),length(disease_list)))
rownames(cmd3_rlm_dir_species_disease) <- select_species
colnames(cmd3_rlm_dir_species_disease) <- disease_list

for(i in 1:length(select_species))
{
	for(j in 1:length(disease_list))
	{
		cmd3_rlm_dir_species_disease[i,j] <- ifelse(cmd3_rlm_qval_species_disease[i,j] <= 0.1,3*sign(cmd3_rlm_est_species_disease[i,j]),ifelse(cmd3_rlm_pval_species_disease[i,j] <= 0.05,2*sign(cmd3_rlm_est_species_disease[i,j]),sign(cmd3_rlm_est_species_disease[i,j])))
	}
}

cmd3_mat_disease_association <- mat_disease_association
save(EU_NA_Studies,EA_Studies,Other_Studies,df_cmd3_diversity_uniqueness,df_cmd3_controls_diversity_uniqueness,cmd3_select_disease_final_species,cmd3_select_age_final_species,cmd3_rlm_est_select_species_disease,cmd3_rlm_pval_select_species_disease,cmd3_rlm_qval_select_species_disease,cmd3_mat_disease_association,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_stage2b_results.RData") 

save.image("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_analysis_2021_Revision.RData")
rm(list=ls())








