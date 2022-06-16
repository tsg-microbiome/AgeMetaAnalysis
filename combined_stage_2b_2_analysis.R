library(robumeta)
library(metafor)
library(dplyr)
library(effsize)
library(MASS)
library(gplots)
library(RColorBrewer)
library(sfsmisc)

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

compute_meta_corr <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(NA,length(grouping_list),3))
	colnames(temp_meta) <- c("dataset","ri","ni")
	
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		var1_length <- length(which(data[data[,grouping_variable]==group,var1]>0))
		var2_length <- length(which(data[data[,grouping_variable]==group,var2]>0))
		temp_meta[i,1] <- group
		if((var1_length > 0)&(var2_length > 0))
		{
			temp_meta[i,2] <- cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2],method="kendall")
			temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
		}
		else
		{
			temp_meta[i,2] <- 0
			temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
		}
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	temp_meta <- escalc(measure="ZCOR",ri=ri,ni=ni,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	res$slabs <- rownames(temp_meta)
	return_list <- list("df_studies"=temp_meta,"model"=res)
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
		var1_length <- length(which(data[data[,grouping_variable]==group,var1]>0))
		var2_length <- length(which(data[data[,grouping_variable]==group,var2]>0))
		#print(group)
		if((var1_length > 0)&(var2_length > 0))
		{
			f <- as.formula(paste0(var1,"~",var2))
			temp_rlm <- rlm(f,data=data[data[,grouping_variable]==group,])
			summary_temp_rlm <- summary(temp_rlm)
			temp_meta[i,2] <- summary_temp_rlm$coefficients[2,3]
			temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
			temp_meta[i,4] <- 1
			temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
			temp_meta[i,6] <- sign(temp_meta[i,2])
		}
		else
		{
			temp_meta[i,2] <- 0
			temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
			temp_meta[i,4] <- 1
			temp_meta[i,5] <- 1
			temp_meta[i,6] <- 0
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

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\he_stage2b_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\ag_stage2b_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\nuage_stage2b_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\logmpie_stage2b_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\odamaki_stage2b_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\taxa_groupings_new.RData")


species_group_16S <- species_groupings

group_colors <- NA
group_colors["G1"] <- rgb(102/255,139/255,69/255)
group_colors["G2"] <- rgb(255/255,0/255,0/255)
group_colors["Others"] <- rgb(51/255,51/255,255/255)


combined_mat_16S_disease_association <- merge(he_mat_disease_association,ag_mat_disease_association,by="row.names",all=TRUE)[,-1]
rownames(combined_mat_16S_disease_association) <- merge(he_mat_disease_association,ag_mat_disease_association,by="row.names",all=TRUE)[,1]
combined_mat_16S_disease_association <- apply(combined_mat_16S_disease_association,2,function(x)(ifelse(is.na(x),0,x)))

G_Markers_16S <- intersect(c("Bacteroides_fragilis","Actinomyces_odontolyticus","Bifidobacterium_dentium","Clostridium_clostridioforme","Clostridium_nexile","Clostridium_ramosum","Dialister_invisus","Eggerthella_lenta","Escherichia_coli","Fusobacterium_nucleatum","Granulicatella_adiacens","Lactococcus_lactis","Rothia_mucilaginosa","Ruminococcus_gnavus","Streptococcus_anginosus","Streptococcus_infantis","Streptococcus_salivarius","Streptococcus_sanguinis","Streptococcus_vestibularis","Subdoligranulum_variabile","Veillonella_atypica","Clostridium_asparagiforme","Clostridium_bolteae","Clostridium_citroniae","Clostridium_hathewayi","Clostridium_symbiosum","Streptococcus_parasanguinis","Solobacterium_moorei","Ruminococcus_torques","Streptococcus_mitis","Klebsiella_pneumoniae","Streptococcus_australis","Streptococcus_gordonii","Enterobacter_cloacae","G_Markers_combined"),rownames(combined_mat_16S_disease_association))

L_Markers_16S <- intersect(c("Eubacterium_hallii","Dorea_longicatena","Coprococcus_comes","Coprococcus_catus","Butyrivibrio_crossotus","Bacteroides_uniformis","Alistipes_shahii","Alistipes_indistinctus","Roseburia_hominis","Pseudoflavonifractor_capillosus","Eubacterium_siraeum","Eubacterium_rectale","Coprobacter_fastidiosus","Bifidobacterium_longum","Bifidobacterium_animalis","Barnesiella_intestinihominis","Bacteroides_xylanisolvens","Alistipes_senegalensis","Alistipes_putredinis","Alistipes_onderdonkii","Akkermansia_muciniphila","Faecalibacterium_prausnitzii","Roseburia_inulinivorans","Roseburia_intestinalis","L_Markers_combined"),rownames(combined_mat_16S_disease_association))

AssociationDirection_16S <- as.data.frame(cbind(apply(combined_mat_16S_disease_association[c(G_Markers_16S,L_Markers_16S),],1,function(x)(length(x[x<=-2]))),apply(combined_mat_16S_disease_association[c(G_Markers_16S,L_Markers_16S),],1,function(x)(length(x[x>=2])))))
colnames(AssociationDirection_16S) <- c("Decreased","Increased")

G_Markers_16S_ordered <- setdiff(intersect(rownames(AssociationDirection_16S[order(AssociationDirection_16S[,1]-AssociationDirection_16S[,2]),]),G_Markers_16S),"G_Markers_combined")

L_Markers_16S_ordered <- setdiff(intersect(rownames(AssociationDirection_16S[order(AssociationDirection_16S[,1]-AssociationDirection_16S[,2]),]),L_Markers_16S),"L_Markers_combined")

Markers_16S_sorted <- c(rownames(AssociationDirection_16S[c(G_Markers_16S_ordered,L_Markers_16S_ordered),]),"G_Markers_combined","L_Markers_combined")

combined_mat_16S_disease_association <- combined_mat_16S_disease_association[setdiff(Markers_16S_sorted,c("G_Markers_combined","L_Markers_combined")),]

heatmap.2(t(combined_mat_16S_disease_association),density="none",trace="none",col=c("skyblue4","skyblue3","white","white","white","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="grey",colsep=1:ncol(t(combined_mat_16S_disease_association)),rowsep=1:nrow(t(combined_mat_16S_disease_association)),key=FALSE,margins=c(10,5),lhei=c(0.1,5),lwid=c(0.1,5),ColSideColors=ifelse(is.na(species_group_16S[rownames(combined_mat_16S_disease_association),7]),"white",group_colors[species_group_16S[rownames(combined_mat_16S_disease_association),7]]))

validated_L_Markers_16S <- L_Markers_16S_ordered[9:23]
validated_G_Markers_16S <- G_Markers_16S_ordered[1:16]

df_ag_diversity_uniqueness$G_Markers_combined <- rowMeans(apply(ag_select_age_final_species[rownames(df_ag_diversity_uniqueness),intersect(validated_G_Markers_16S,colnames(ag_select_age_final_species))],2,rank_scale))

df_ag_diversity_uniqueness$L_Markers_combined <- rowMeans(apply(ag_select_age_final_species[rownames(df_ag_diversity_uniqueness),intersect(validated_L_Markers_16S,colnames(ag_select_age_final_species))],2,rank_scale))

df_ag_diversity_uniqueness[,intersect(G_Markers_16S,colnames(ag_select_age_final_species))] <- ag_select_age_final_species[rownames(df_ag_diversity_uniqueness),intersect(G_Markers_16S,colnames(ag_select_age_final_species))]

df_ag_diversity_uniqueness[,intersect(L_Markers_16S,colnames(ag_select_age_final_species))] <- ag_select_age_final_species[rownames(df_ag_diversity_uniqueness),intersect(L_Markers_16S,colnames(ag_select_age_final_species))]

df_he_diversity_uniqueness$G_Markers_combined <- rowMeans(apply(he_select_age_final_species[rownames(df_he_diversity_uniqueness),intersect(validated_G_Markers_16S,colnames(he_select_age_final_species))],2,rank_scale))

df_he_diversity_uniqueness$L_Markers_combined <- rowMeans(apply(he_select_age_final_species[rownames(df_he_diversity_uniqueness),intersect(validated_L_Markers_16S,colnames(he_select_age_final_species))],2,rank_scale))

df_he_diversity_uniqueness[,intersect(G_Markers_16S,colnames(he_select_age_final_species))] <- he_select_age_final_species[rownames(df_he_diversity_uniqueness),intersect(G_Markers_16S,colnames(he_select_age_final_species))]

df_he_diversity_uniqueness[,intersect(L_Markers_16S,colnames(he_select_age_final_species))] <- he_select_age_final_species[rownames(df_he_diversity_uniqueness),intersect(L_Markers_16S,colnames(he_select_age_final_species))]

df_nuage_diversity_uniqueness$G_Markers_combined <- rowMeans(apply(nuage_select_age_final_species[rownames(df_nuage_diversity_uniqueness),intersect(validated_G_Markers_16S,colnames(nuage_select_age_final_species))],2,rank_scale))

df_nuage_diversity_uniqueness$L_Markers_combined <- rowMeans(apply(nuage_select_age_final_species[rownames(df_nuage_diversity_uniqueness),intersect(validated_L_Markers_16S,colnames(nuage_select_age_final_species))],2,rank_scale))

df_nuage_diversity_uniqueness[,intersect(G_Markers_16S,colnames(nuage_select_age_final_species))] <- nuage_select_age_final_species[rownames(df_nuage_diversity_uniqueness),intersect(G_Markers_16S,colnames(nuage_select_age_final_species))]

df_nuage_diversity_uniqueness[,intersect(L_Markers_16S,colnames(nuage_select_age_final_species))] <- nuage_select_age_final_species[rownames(df_nuage_diversity_uniqueness),intersect(L_Markers_16S,colnames(nuage_select_age_final_species))]

df_odamaki_diversity_uniqueness$G_Markers_combined <- rowMeans(apply(odamaki_select_age_final_species[rownames(df_odamaki_diversity_uniqueness),intersect(validated_G_Markers_16S,colnames(odamaki_select_age_final_species))],2,rank_scale))

df_odamaki_diversity_uniqueness$L_Markers_combined <- rowMeans(apply(odamaki_select_age_final_species[rownames(df_odamaki_diversity_uniqueness),intersect(validated_L_Markers_16S,colnames(odamaki_select_age_final_species))],2,rank_scale))

df_odamaki_diversity_uniqueness[,intersect(G_Markers_16S,colnames(odamaki_select_age_final_species))] <- odamaki_select_age_final_species[rownames(df_odamaki_diversity_uniqueness),intersect(G_Markers_16S,colnames(odamaki_select_age_final_species))]

df_odamaki_diversity_uniqueness[,intersect(L_Markers_16S,colnames(odamaki_select_age_final_species))] <- odamaki_select_age_final_species[rownames(df_odamaki_diversity_uniqueness),intersect(L_Markers_16S,colnames(odamaki_select_age_final_species))]

df_logmpie_diversity_uniqueness$G_Markers_combined <- rowMeans(apply(logmpie_select_age_final_species[rownames(df_logmpie_diversity_uniqueness),intersect(validated_G_Markers_16S,colnames(logmpie_select_age_final_species))],2,rank_scale))

df_logmpie_diversity_uniqueness$L_Markers_combined <- rowMeans(apply(logmpie_select_age_final_species[rownames(df_logmpie_diversity_uniqueness),intersect(validated_L_Markers_16S,colnames(logmpie_select_age_final_species))],2,rank_scale))

df_logmpie_diversity_uniqueness[,intersect(G_Markers_16S,colnames(logmpie_select_age_final_species))] <- logmpie_select_age_final_species[rownames(df_logmpie_diversity_uniqueness),intersect(G_Markers_16S,colnames(logmpie_select_age_final_species))]

df_logmpie_diversity_uniqueness[,intersect(L_Markers_16S,colnames(logmpie_select_age_final_species))] <- logmpie_select_age_final_species[rownames(df_logmpie_diversity_uniqueness),intersect(L_Markers_16S,colnames(logmpie_select_age_final_species))]

temp0 <- merge(t(df_ag_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(G_Markers_16S,L_Markers_16S),colnames(df_ag_diversity_uniqueness)))]),t(df_he_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(G_Markers_16S,L_Markers_16S),colnames(df_he_diversity_uniqueness)))]),by="row.names",all=TRUE)[,-1]
rownames(temp0) <- merge(t(df_ag_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(G_Markers_16S,L_Markers_16S),colnames(df_ag_diversity_uniqueness)))]),t(df_he_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(G_Markers_16S,L_Markers_16S),colnames(df_he_diversity_uniqueness)))]),by="row.names",all=TRUE)[,1]
temp0<- as.data.frame(apply(temp0,1,function(x)(ifelse(is.na(x),0,x))))

temp1 <- merge(t(temp0),t(df_nuage_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(G_Markers_16S,L_Markers_16S),colnames(df_nuage_diversity_uniqueness)))]),by="row.names",all=TRUE)[,-1]
rownames(temp1) <- merge(t(temp0),t(df_nuage_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(G_Markers_16S,L_Markers_16S),colnames(df_nuage_diversity_uniqueness)))]),by="row.names",all=TRUE)[,1]
temp1 <- as.data.frame(apply(temp1,1,function(x)(ifelse(is.na(x),0,x))))

temp2 <- merge(t(temp1),t(df_odamaki_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(G_Markers_16S,L_Markers_16S),colnames(df_odamaki_diversity_uniqueness)))]),by="row.names",all=TRUE)[,-1]
rownames(temp2) <- merge(t(temp1),t(df_odamaki_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(G_Markers_16S,L_Markers_16S),colnames(df_odamaki_diversity_uniqueness)))]),by="row.names",all=TRUE)[,1]
temp2 <- as.data.frame(apply(temp2,1,function(x)(ifelse(is.na(x),0,x))))

temp3 <- merge(t(temp2),t(df_logmpie_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(G_Markers_16S,L_Markers_16S),colnames(df_logmpie_diversity_uniqueness)))]),by="row.names",all=TRUE)[,-1]
rownames(temp3) <- merge(t(temp2),t(df_logmpie_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(G_Markers_16S,L_Markers_16S),colnames(df_logmpie_diversity_uniqueness)))]),by="row.names",all=TRUE)[,1]
temp3 <- as.data.frame(apply(temp3,1,function(x)(ifelse(is.na(x),0,x))))

df_16S_diversity_uniqueness <- temp3
df_16S_diversity_uniqueness$study_name <- NA
df_16S_diversity_uniqueness[rownames(df_ag_diversity_uniqueness),"study_name"] <- "AG"
df_16S_diversity_uniqueness[rownames(df_he_diversity_uniqueness),"study_name"] <- "He"
df_16S_diversity_uniqueness[rownames(df_nuage_diversity_uniqueness),"study_name"] <- "NUAGE"
df_16S_diversity_uniqueness[rownames(df_odamaki_diversity_uniqueness),"study_name"] <- "Odamaki"
df_16S_diversity_uniqueness[rownames(df_logmpie_diversity_uniqueness),"study_name"] <- "LogMPie"

validated_Markers_16S <- c(validated_G_Markers_16S,validated_L_Markers_16S)
study_list_16S <- unique(df_16S_diversity_uniqueness$study_name)
markers_16S_age_est_study_level <- as.data.frame(matrix(0,length(validated_Markers_16S),length(study_list_16S)))
rownames(markers_16S_age_est_study_level) <- validated_Markers_16S
colnames(markers_16S_age_est_study_level) <- study_list_16S
markers_16S_age_pval_study_level <- as.data.frame(matrix(0,length(validated_Markers_16S),length(study_list_16S)))
rownames(markers_16S_age_pval_study_level) <- validated_Markers_16S
colnames(markers_16S_age_pval_study_level) <- study_list_16S

for(i in 1:length(validated_Markers_16S))
{
	species_name <- validated_Markers_16S[i]
	for(j in 1:length(study_list_16S))
	{
		#print(paste0(species_name,",",disease_name))
		study_name <- study_list_16S[j]
		sample_list <- rownames(df_16S_diversity_uniqueness[df_16S_diversity_uniqueness$study_name == study_name,])
		vec_species <- df_16S_diversity_uniqueness[sample_list,species_name]
		vec_age <- df_16S_diversity_uniqueness[sample_list,"age"]
		if((length(vec_species[vec_species > 0])>0)&(length(vec_age[vec_age > 0])>0))
		{
			temp_rlm <- rlm(vec_species~vec_age)
			summary_temp_rlm <- summary(temp_rlm)
			if((!is.nan(summary_temp_rlm$coefficients[2,3]))&(summary_temp_rlm$coefficients[2,1] != 0)&(abs(summary_temp_rlm$coefficients[2,3]) <= 10000))
			{
				markers_16S_age_est_study_level[species_name,study_name] <- summary_temp_rlm$coefficients[2,3]
				markers_16S_age_pval_study_level[species_name,study_name] <- f.robftest(temp_rlm)$p.value
			}
			else
			{
				markers_16S_age_est_study_level[species_name,study_name] <- 0
				markers_16S_age_pval_study_level[species_name,study_name] <- 1
			}
		}
	}
}

markers_16S_age_qval_study_level <- t(apply(markers_16S_age_pval_study_level,1,function(x)(p.adjust(x,method="fdr"))))

markers_16S_age_dir_study_level <- as.data.frame(matrix(0,length(validated_Markers_16S),length(study_list_16S)))
rownames(markers_16S_age_dir_study_level) <- validated_Markers_16S
colnames(markers_16S_age_dir_study_level) <- study_list_16S
for(i in 1:length(validated_Markers_16S))
{
	species_name <- validated_Markers_16S[i]
	for(j in 1:length(study_list_16S))
	{
		markers_16S_age_dir_study_level[i,j] <- ifelse(markers_16S_age_qval_study_level[i,j] <= 0.1,3*sign(markers_16S_age_est_study_level[i,j]),ifelse(markers_16S_age_pval_study_level[i,j] <= 0.05,2*sign(markers_16S_age_est_study_level[i,j]),sign(markers_16S_age_est_study_level[i,j])))
	}
}

markers_16S_rem_age <- as.data.frame(matrix(NA,length(validated_Markers_16S),5))
rownames(markers_16S_rem_age) <- validated_Markers_16S
colnames(markers_16S_rem_age) <- c("PVal","QVal","REM_Med","REM_Upper","REM_Lower")
for(i in 1:length(validated_Markers_16S))
{
	species_name <- validated_Markers_16S[i]
	rem_marker <- compute_meta_lm(df_16S_diversity_uniqueness,species_name,"age","study_name",setdiff(study_list_16S,"NUAGE"))
	markers_16S_rem_age[i,1] <- as.numeric(rem_marker$model$pval)
	markers_16S_rem_age[i,3] <- as.numeric(rem_marker$model$beta)
	markers_16S_rem_age[i,4] <- as.numeric(rem_marker$model$ci.ub)
	markers_16S_rem_age[i,5] <- as.numeric(rem_marker$model$ci.lb)
}
markers_16S_rem_age[,2] <- p.adjust(markers_16S_rem_age[,1],method="fdr")

###


species_group_shotgun <- species_groupings

group_colors <- NA
group_colors["G1"] <- rgb(102/255,139/255,69/255)
group_colors["G2"] <- rgb(255/255,0/255,0/255)
group_colors["Others"] <- rgb(51/255,51/255,255/255)

load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\cmd3_stage2b_results.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\isc_stage2b_results.RData")

## Doing name change to revert to old taxonomy scheme
## Earlier Clostridium_nexile is now Tyzzerella_nexilis reverted back
colnames(cmd3_select_age_final_species)[267] <- "Clostridium_nexile"
## Earlier Clostridium_hathewayi now Hungatella_hathewayi reverted back
colnames(cmd3_select_age_final_species)[261] <- "Clostridium_hathewayi"
## Earlier Earlier Streptococcus_anginosus is now Streptococcus_anginosus_group reverted back
colnames(cmd3_select_age_final_species)[389] <- "Streptococcus_anginosus"
## Earlier Enterobacter_cloacae is now Enterobacter_cloacae_group reverted back
colnames(cmd3_select_age_final_species)[25] <- "Enterobacter_cloacae"
## Earlier Clostridium_ramosum is now Erysipelatoclostridium_ramosum reverted back
colnames(cmd3_select_age_final_species)[42] <- "Clostridium_ramosum"

combined_mat_shotgun_disease_association <- merge(cmd3_mat_disease_association,isc_mat_disease_association,by="row.names",all=TRUE)[,-1]
rownames(combined_mat_shotgun_disease_association) <- merge(cmd3_mat_disease_association,isc_mat_disease_association,by="row.names",all=TRUE)[,1]
combined_mat_shotgun_disease_association <- apply(combined_mat_shotgun_disease_association,2,function(x)(ifelse(is.na(x),0,x)))

G_Markers_shotgun <- intersect(c("Bacteroides_fragilis","Actinomyces_odontolyticus","Bifidobacterium_dentium","Clostridium_clostridioforme","Clostridium_nexile","Clostridium_ramosum","Dialister_invisus","Eggerthella_lenta","Escherichia_coli","Fusobacterium_nucleatum","Granulicatella_adiacens","Lactococcus_lactis","Rothia_mucilaginosa","Ruminococcus_gnavus","Streptococcus_anginosus","Streptococcus_infantis","Streptococcus_salivarius","Streptococcus_sanguinis","Streptococcus_vestibularis","Veillonella_atypica","Clostridium_asparagiforme","Clostridiales_bacterium_1_7_47FAA","Clostridium_bolteae","Clostridium_citroniae","Clostridium_hathewayi","Clostridium_symbiosum","Streptococcus_parasanguinis","Solobacterium_moorei","Ruminococcus_torques","Streptococcus_mitis","Klebsiella_pneumoniae","Streptococcus_australis","Streptococcus_gordonii","Enterobacter_cloacae"),rownames(combined_mat_shotgun_disease_association))

L_Markers_shotgun <- intersect(c("Eubacterium_hallii","Dorea_longicatena","Coprococcus_comes","Coprococcus_catus","Butyrivibrio_crossotus","Bacteroides_uniformis","Alistipes_shahii","Alistipes_indistinctus","Roseburia_hominis","Pseudoflavonifractor_capillosus","Eubacterium_siraeum","Eubacterium_rectale","Coprobacter_fastidiosus","Bifidobacterium_longum","Bifidobacterium_animalis","Barnesiella_intestinihominis","Bacteroides_xylanisolvens","Alistipes_senegalensis","Alistipes_putredinis","Alistipes_onderdonkii","Akkermansia_muciniphila","Faecalibacterium_prausnitzii","Roseburia_inulinivorans","Roseburia_intestinalis"),rownames(combined_mat_shotgun_disease_association))

AssociationDirection_shotgun <- as.data.frame(cbind(apply(combined_mat_shotgun_disease_association[c(G_Markers_shotgun,L_Markers_shotgun),],1,function(x)(length(x[x<=-2]))),apply(combined_mat_shotgun_disease_association[c(G_Markers_shotgun,L_Markers_shotgun),],1,function(x)(length(x[x>=2])))))
colnames(AssociationDirection_shotgun) <- c("Decreased","Increased")

G_Markers_shotgun_ordered <- intersect(rownames(AssociationDirection_shotgun[order(AssociationDirection_shotgun[,1]-AssociationDirection_shotgun[,2]),]),G_Markers_shotgun)

L_Markers_shotgun_ordered <- intersect(rownames(AssociationDirection_shotgun[order(AssociationDirection_shotgun[,1]-AssociationDirection_shotgun[,2]),]),L_Markers_shotgun)

Markers_shotgun_sorted <- rownames(AssociationDirection_shotgun[c(G_Markers_shotgun_ordered,L_Markers_shotgun_ordered),])

combined_mat_shotgun_disease_association <- combined_mat_shotgun_disease_association[Markers_shotgun_sorted,]

heatmap.2(t(combined_mat_shotgun_disease_association),density="none",trace="none",col=c("skyblue4","skyblue3","white","white","white","orangered3","orangered4"),Colv=FALSE,Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="grey",colsep=1:ncol(t(combined_mat_shotgun_disease_association)),rowsep=1:nrow(t(combined_mat_shotgun_disease_association)),key=FALSE,margins=c(12,5),lhei=c(0.1,5),lwid=c(0.1,5),ColSideColors=ifelse(is.na(species_group_shotgun[rownames(combined_mat_shotgun_disease_association),7]),"white",group_colors[species_group_shotgun[rownames(combined_mat_shotgun_disease_association),7]]))

validated_L_Markers_shotgun <- L_Markers_shotgun_ordered[7:24]
validated_G_Markers_shotgun <- G_Markers_shotgun_ordered[1:24]

df_cmd3_controls_diversity_uniqueness$G_Markers_combined <- rank_scale(apply(cmd3_select_age_final_species[rownames(df_cmd3_controls_diversity_uniqueness),intersect(validated_G_Markers_shotgun,colnames(cmd3_select_age_final_species))],1,sum))

df_cmd3_controls_diversity_uniqueness$L_Markers_combined <- rank_scale(apply(cmd3_select_age_final_species[rownames(df_cmd3_controls_diversity_uniqueness),intersect(validated_L_Markers_shotgun,colnames(cmd3_select_age_final_species))],1,sum))

df_cmd3_controls_diversity_uniqueness[,intersect(validated_L_Markers_shotgun,colnames(cmd3_select_age_final_species))] <- cmd3_select_age_final_species[rownames(df_cmd3_controls_diversity_uniqueness),intersect(validated_L_Markers_shotgun,colnames(cmd3_select_age_final_species))]

df_cmd3_controls_diversity_uniqueness[,intersect(validated_G_Markers_shotgun,colnames(cmd3_select_age_final_species))] <- cmd3_select_age_final_species[rownames(df_cmd3_controls_diversity_uniqueness),intersect(validated_G_Markers_shotgun,colnames(cmd3_select_age_final_species))]

#df_cmd3_diversity_uniqueness <- df_cmd3_diversity_uniqueness[intersect()]

df_cmd3_diversity_uniqueness$G_Markers_combined <- rowMeans(apply(cmd3_select_age_final_species[rownames(df_cmd3_diversity_uniqueness),intersect(validated_G_Markers_shotgun,colnames(cmd3_select_age_final_species))],2,rank_scale))

df_cmd3_diversity_uniqueness$L_Markers_combined <- rowMeans(apply(cmd3_select_age_final_species[rownames(df_cmd3_diversity_uniqueness),intersect(validated_L_Markers_shotgun,colnames(cmd3_select_age_final_species))],2,rank_scale))

#df_cmd3_diversity_uniqueness$G_Markers_combined <- rowMeans(apply(cmd3_select_age_final_species[rownames(df_cmd3_diversity_uniqueness),intersect(validated_G_Markers_shotgun,colnames(cmd3_select_age_final_species))],2,function(x)(ifelse(x>0,1,0))))

#df_cmd3_diversity_uniqueness$L_Markers_combined <- rowMeans(apply(cmd3_select_age_final_species[rownames(df_cmd3_diversity_uniqueness),intersect(validated_L_Markers_shotgun,colnames(cmd3_select_age_final_species))],2,function(x)(ifelse(x>0,1,0))))

df_cmd3_diversity_uniqueness[,intersect(validated_G_Markers_shotgun,colnames(cmd3_select_age_final_species))] <- cmd3_select_age_final_species[rownames(df_cmd3_diversity_uniqueness),intersect(validated_G_Markers_shotgun,colnames(cmd3_select_age_final_species))]

df_cmd3_diversity_uniqueness[,intersect(validated_L_Markers_shotgun,colnames(cmd3_select_age_final_species))] <- cmd3_select_age_final_species[rownames(df_cmd3_diversity_uniqueness),intersect(validated_L_Markers_shotgun,colnames(cmd3_select_age_final_species))]

df_isc_diversity_uniqueness$G_Markers_combined <- rowMeans(apply(isc_select_age_final_species[rownames(df_isc_diversity_uniqueness),intersect(validated_G_Markers_shotgun,colnames(isc_select_age_final_species))],2,rank_scale))

df_isc_diversity_uniqueness$L_Markers_combined <- rowMeans(apply(isc_select_age_final_species[rownames(df_isc_diversity_uniqueness),intersect(validated_L_Markers_shotgun,colnames(isc_select_age_final_species))],2,rank_scale))

#df_isc_diversity_uniqueness$G_Markers_combined <- rowMeans(apply(isc_select_age_final_species[rownames(df_isc_diversity_uniqueness),intersect(validated_G_Markers_shotgun,colnames(isc_select_age_final_species))],2,function(x)(ifelse(x>0,1,0))))

#df_isc_diversity_uniqueness$L_Markers_combined <- rowMeans(apply(isc_select_age_final_species[rownames(df_isc_diversity_uniqueness),intersect(validated_L_Markers_shotgun,colnames(isc_select_age_final_species))],2,function(x)(ifelse(x>0,1,0))))

df_isc_diversity_uniqueness[,intersect(validated_G_Markers_shotgun,colnames(isc_select_age_final_species))] <- isc_select_age_final_species[rownames(df_isc_diversity_uniqueness),intersect(validated_G_Markers_shotgun,colnames(isc_select_age_final_species))]

df_isc_diversity_uniqueness[,intersect(validated_L_Markers_shotgun,colnames(isc_select_age_final_species))] <- isc_select_age_final_species[rownames(df_isc_diversity_uniqueness),intersect(validated_L_Markers_shotgun,colnames(isc_select_age_final_species))]

temp <- merge(t(df_cmd3_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(validated_G_Markers_shotgun,validated_L_Markers_shotgun),colnames(df_cmd3_diversity_uniqueness)))]),t(df_isc_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(validated_G_Markers_shotgun,validated_L_Markers_shotgun),colnames(df_isc_diversity_uniqueness)))]),by="row.names",all=TRUE)[,-1]
rownames(temp) <- merge(t(df_cmd3_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(validated_G_Markers_shotgun,validated_L_Markers_shotgun),colnames(df_cmd3_diversity_uniqueness)))]),t(df_isc_diversity_uniqueness[,c("age","species_shannon","species_bray_uniqueness","species_jaccard_uniqueness","species_aitchison_uniqueness","species_kendall_uniqueness","G_Markers_combined","L_Markers_combined",intersect(c(validated_G_Markers_shotgun,validated_L_Markers_shotgun),colnames(df_isc_diversity_uniqueness)))]),by="row.names",all=TRUE)[,1]
df_shotgun_diversity_uniqueness <- as.data.frame(apply(temp,1,function(x)(ifelse(is.na(x),0,x))))
df_shotgun_diversity_uniqueness$study_name <- NA
df_shotgun_diversity_uniqueness[rownames(df_cmd3_diversity_uniqueness),"study_name"] <- df_cmd3_diversity_uniqueness$study_name
df_shotgun_diversity_uniqueness[rownames(df_isc_diversity_uniqueness),"study_name"] <- "ISC"

EU_NA_Studies <- c(EU_NA_Studies,"ISC")

df_shotgun_diversity_uniqueness$age_cat <- ifelse(df_shotgun_diversity_uniqueness$age<40,1,ifelse(df_shotgun_diversity_uniqueness$age<60,2,3))

validated_Markers_shotgun <- c(validated_G_Markers_shotgun,validated_L_Markers_shotgun)
shotgun_study_list <- c(EU_NA_Studies,EA_Studies,Other_Studies)
shotgun_markers_age_est_study_level <- as.data.frame(matrix(0,length(validated_Markers_shotgun),length(shotgun_study_list)))
rownames(shotgun_markers_age_est_study_level) <- validated_Markers_shotgun
colnames(shotgun_markers_age_est_study_level) <- shotgun_study_list
shotgun_markers_age_pval_study_level <- as.data.frame(matrix(1,length(validated_Markers_shotgun),length(shotgun_study_list)))
rownames(shotgun_markers_age_pval_study_level) <- validated_Markers_shotgun
colnames(shotgun_markers_age_pval_study_level) <- shotgun_study_list

for(i in 1:length(validated_Markers_shotgun))
{
	species_name <- validated_Markers_shotgun[i]
	for(j in 1:length(shotgun_study_list))
	{
		#print(paste0(species_name,",",disease_name))
		study_name <- shotgun_study_list[j]
		sample_list <- rownames(df_shotgun_diversity_uniqueness[df_shotgun_diversity_uniqueness$study_name == study_name,])
		vec_species <- df_shotgun_diversity_uniqueness[sample_list,species_name]
		vec_age <- df_shotgun_diversity_uniqueness[sample_list,"age"]
		if((length(vec_species[vec_species > 0])>0)&(length(vec_age[vec_age > 0])>0))
		{
			temp_rlm <- rlm(vec_species~vec_age)
			summary_temp_rlm <- summary(temp_rlm)
			if((!is.nan(summary_temp_rlm$coefficients[2,3]))&(summary_temp_rlm$coefficients[2,1] != 0)&(abs(summary_temp_rlm$coefficients[2,3]) <= 10000))
			{
				shotgun_markers_age_est_study_level[species_name,study_name] <- summary_temp_rlm$coefficients[2,3]
				shotgun_markers_age_pval_study_level[species_name,study_name] <- f.robftest(temp_rlm)$p.value
			}
			else
			{
				shotgun_markers_age_est_study_level[species_name,study_name] <- 0
				shotgun_markers_age_pval_study_level[species_name,study_name] <- 1
			}
		}
	}
}

shotgun_markers_age_qval_study_level <- t(apply(shotgun_markers_age_pval_study_level,1,function(x)(p.adjust(x,method="fdr"))))

shotgun_markers_age_dir_study_level <- as.data.frame(matrix(0,length(validated_Markers_shotgun),length(shotgun_study_list)))
rownames(shotgun_markers_age_dir_study_level) <- validated_Markers_shotgun
colnames(shotgun_markers_age_dir_study_level) <- shotgun_study_list
for(i in 1:length(validated_Markers_shotgun))
{
	species_name <- validated_Markers_shotgun[i]
	for(j in 1:length(shotgun_study_list))
	{
		shotgun_markers_age_dir_study_level[i,j] <- ifelse(shotgun_markers_age_qval_study_level[i,j] <= 0.1,3*sign(shotgun_markers_age_est_study_level[i,j]),ifelse(shotgun_markers_age_pval_study_level[i,j] <= 0.05,2*sign(shotgun_markers_age_est_study_level[i,j]),sign(shotgun_markers_age_est_study_level[i,j])))
	}
}

shotgun_markers_age_est_study_grp_level <- as.data.frame(matrix(0,length(validated_Markers_shotgun),4))
rownames(shotgun_markers_age_est_study_grp_level) <- validated_Markers_shotgun
colnames(shotgun_markers_age_est_study_grp_level) <- c("All_Studies","EU_NA_Studies","EA_Studies","Other_Studies")
shotgun_markers_age_upper_study_grp_level <- as.data.frame(matrix(1,length(validated_Markers_shotgun),4))
rownames(shotgun_markers_age_upper_study_grp_level) <- validated_Markers_shotgun
colnames(shotgun_markers_age_upper_study_grp_level) <- c("All_Studies","EU_NA_Studies","EA_Studies","Other_Studies")
shotgun_markers_age_lower_study_grp_level <- as.data.frame(matrix(1,length(validated_Markers_shotgun),4))
rownames(shotgun_markers_age_lower_study_grp_level) <- validated_Markers_shotgun
colnames(shotgun_markers_age_lower_study_grp_level) <- c("All_Studies","EU_NA_Studies","EA_Studies","Other_Studies")
shotgun_markers_age_pval_study_grp_level <- as.data.frame(matrix(1,length(validated_Markers_shotgun),4))
rownames(shotgun_markers_age_pval_study_grp_level) <- validated_Markers_shotgun
colnames(shotgun_markers_age_pval_study_grp_level) <- c("All_Studies","EU_NA_Studies","EA_Studies","Other_Studies")
for(i in 1:length(validated_Markers_shotgun))
{
	species_name <- validated_Markers_shotgun[i]
	res1 <- compute_meta_lm(df_shotgun_diversity_uniqueness,species_name,"age","study_name",c(EU_NA_Studies,EA_Studies,Other_Studies))
	res2 <- compute_meta_lm(df_shotgun_diversity_uniqueness,species_name,"age","study_name",c(EU_NA_Studies))
	res3 <- compute_meta_lm(df_shotgun_diversity_uniqueness,species_name,"age","study_name",c(EA_Studies))
	res4 <- compute_meta_lm(df_shotgun_diversity_uniqueness,species_name,"age","study_name",c(Other_Studies))
	
	shotgun_markers_age_est_study_grp_level[i,1] <- as.numeric(res1$model$beta)
	shotgun_markers_age_pval_study_grp_level[i,1] <- as.numeric(res1$model$pval)
	shotgun_markers_age_upper_study_grp_level[i,1] <- as.numeric(res1$model$ci.ub)
	shotgun_markers_age_lower_study_grp_level[i,1] <- as.numeric(res1$model$ci.lb)
	
	shotgun_markers_age_est_study_grp_level[i,2] <- as.numeric(res2$model$beta)
	shotgun_markers_age_pval_study_grp_level[i,2] <- as.numeric(res2$model$pval)
	shotgun_markers_age_upper_study_grp_level[i,2] <- as.numeric(res2$model$ci.ub)
	shotgun_markers_age_lower_study_grp_level[i,2] <- as.numeric(res2$model$ci.lb)
	
	shotgun_markers_age_est_study_grp_level[i,3] <- as.numeric(res3$model$beta)
	shotgun_markers_age_pval_study_grp_level[i,3] <- as.numeric(res3$model$pval)
	shotgun_markers_age_upper_study_grp_level[i,3] <- as.numeric(res3$model$ci.ub)
	shotgun_markers_age_lower_study_grp_level[i,3] <- as.numeric(res3$model$ci.lb)
	
	shotgun_markers_age_est_study_grp_level[i,4] <- as.numeric(res4$model$beta)
	shotgun_markers_age_pval_study_grp_level[i,4] <- as.numeric(res4$model$pval)
	shotgun_markers_age_upper_study_grp_level[i,4] <- as.numeric(res4$model$ci.ub)
	shotgun_markers_age_lower_study_grp_level[i,4] <- as.numeric(res4$model$ci.lb)
	
}
shotgun_markers_age_qval_study_grp_level <- t(apply(shotgun_markers_age_pval_study_grp_level,1,function(x)(p.adjust(x,method="fdr"))))
shotgun_markers_age_dir_study_grp_level <- as.data.frame(matrix(0,length(validated_Markers_shotgun),4))
rownames(shotgun_markers_age_dir_study_grp_level) <- validated_Markers_shotgun
colnames(shotgun_markers_age_dir_study_grp_level) <- c("All_Studies","EU_NA_Studies","EA_Studies","Other_Studies")
for(i in 1:length(validated_Markers_shotgun))
{
	for(j in 1:4)
	{
		shotgun_markers_age_dir_study_grp_level[i,j] <- ifelse(shotgun_markers_age_qval_study_grp_level[i,j] <= 0.1,3*sign(shotgun_markers_age_est_study_grp_level[i,j]),ifelse(shotgun_markers_age_pval_study_grp_level[i,j] <= 0.05,2*sign(shotgun_markers_age_est_study_grp_level[i,j]),sign(shotgun_markers_age_est_study_grp_level[i,j])))
	}
}


study_dummy <- 1:length(shotgun_study_list)
names(study_dummy) <- shotgun_study_list
for(i in 1:length(shotgun_study_list))
{
	study_dummy[i] <- as.integer(runif(1, min=0, max=50))
}

df_shotgun_diversity_uniqueness$study_dummy <- as.numeric(study_dummy[df_shotgun_diversity_uniqueness$study_name])

AgeAssociationDirection_shotgun <- as.data.frame(cbind(apply(shotgun_markers_age_dir_study_level,1,function(x)(length(x[x<=-2]))),apply(shotgun_markers_age_dir_study_level,1,function(x)(length(x[x>=2])))))
AgeAssociationDirection_shotgun$Rank <- rank(shotgun_markers_age_est_study_grp_level[,1])+rank(AgeAssociationDirection_shotgun[,2]-AgeAssociationDirection_shotgun[,1])
shotgun_Markers_age_ordered <- rownames(AgeAssociationDirection_shotgun[order(AgeAssociationDirection_shotgun$Rank),])


shotgun_markers_age_rlm_est_study_grp_level <- as.data.frame(matrix(0,length(validated_Markers_shotgun),4))
rownames(shotgun_markers_age_rlm_est_study_grp_level) <- validated_Markers_shotgun
colnames(shotgun_markers_age_rlm_est_study_grp_level) <- c("All_Studies","EU_NA_Studies","EA_Studies","Other_Studies")
shotgun_markers_age_rlm_pval_study_grp_level <- as.data.frame(matrix(1,length(validated_Markers_shotgun),4))
rownames(shotgun_markers_age_rlm_pval_study_grp_level) <- validated_Markers_shotgun
colnames(shotgun_markers_age_rlm_pval_study_grp_level) <- c("All_Studies","EU_NA_Studies","EA_Studies","Other_Studies")
for(i in 1:length(validated_Markers_shotgun))
{
	species_name <- validated_Markers_shotgun[i]
	
	
	vec1 <- df_shotgun_diversity_uniqueness[df_shotgun_diversity_uniqueness$study_name %in% c(EU_NA_Studies,EA_Studies,Other_Studies),species_name]
	if(length(vec1[vec1>0])>0)
	{
		rlm1 <- rlm(as.formula(paste0(species_name,"~study_dummy+age")),df_shotgun_diversity_uniqueness[df_shotgun_diversity_uniqueness$study_name %in% c(EU_NA_Studies,EA_Studies,Other_Studies),])
		shotgun_markers_age_rlm_est_study_grp_level[i,1] <- summary(rlm1)$coefficients[3,3]
		shotgun_markers_age_rlm_pval_study_grp_level[i,1] <- f.robftest(rlm1,var="age")$p.value
	}
	else
	{
		shotgun_markers_age_rlm_est_study_grp_level[i,1] <- 0
		shotgun_markers_age_rlm_pval_study_grp_level[i,1] <- 1
	}
	
	vec2 <- df_shotgun_diversity_uniqueness[df_shotgun_diversity_uniqueness$study_name %in% c(EU_NA_Studies),species_name]
	if(length(vec2[vec2>0])>0)
	{
		rlm2 <- rlm(as.formula(paste0(species_name,"~study_dummy+age")),df_shotgun_diversity_uniqueness[df_shotgun_diversity_uniqueness$study_name %in% c(EU_NA_Studies),])
		shotgun_markers_age_rlm_est_study_grp_level[i,2] <- summary(rlm2)$coefficients[3,3]
		shotgun_markers_age_rlm_pval_study_grp_level[i,2] <- f.robftest(rlm2,var="age")$p.value
	}
	else
	{
		shotgun_markers_age_rlm_est_study_grp_level[i,2] <- 0
		shotgun_markers_age_rlm_pval_study_grp_level[i,2] <- 1
	}
	
	vec3 <- df_shotgun_diversity_uniqueness[df_shotgun_diversity_uniqueness$study_name %in% c(EA_Studies),species_name]
	if(length(vec3[vec3>0])>0)
	{
		rlm3 <- rlm(as.formula(paste0(species_name,"~study_dummy+age")),df_shotgun_diversity_uniqueness[df_shotgun_diversity_uniqueness$study_name %in% c(EA_Studies),])
		shotgun_markers_age_rlm_est_study_grp_level[i,3] <- summary(rlm3)$coefficients[3,3]
		shotgun_markers_age_rlm_pval_study_grp_level[i,3] <-  f.robftest(rlm3,var="age")$p.value
	}
	else
	{
		shotgun_markers_age_rlm_est_study_grp_level[i,3] <- 0
		shotgun_markers_age_rlm_pval_study_grp_level[i,3] <- 1
	}
	
	vec4 <- df_shotgun_diversity_uniqueness[df_shotgun_diversity_uniqueness$study_name %in% c(Other_Studies),species_name]
	if(length(vec4[vec4>0])>0)
	{
		rlm4 <- rlm(as.formula(paste0(species_name,"~study_dummy+age")),df_shotgun_diversity_uniqueness[df_shotgun_diversity_uniqueness$study_name %in% c(Other_Studies),])
		shotgun_markers_age_rlm_est_study_grp_level[i,4] <- summary(rlm4)$coefficients[3,3]
		shotgun_markers_age_rlm_pval_study_grp_level[i,4] <- f.robftest(rlm4,var="age")$p.value
	}
	else
	{
		shotgun_markers_age_rlm_est_study_grp_level[i,4] <- 0
		shotgun_markers_age_rlm_pval_study_grp_level[i,4] <- 1
	}
	
}
shotgun_markers_age_rlm_qval_study_grp_level <- t(apply(shotgun_markers_age_rlm_pval_study_grp_level,1,function(x)(p.adjust(x,method="fdr"))))
shotgun_markers_age_rlm_dir_study_grp_level <- as.data.frame(matrix(0,length(validated_Markers_shotgun),4))
rownames(shotgun_markers_age_rlm_dir_study_grp_level) <- validated_Markers_shotgun
colnames(shotgun_markers_age_rlm_dir_study_grp_level) <- c("All_Studies","EU_NA_Studies","EA_Studies","Other_Studies")
for(i in 1:length(validated_Markers_shotgun))
{
	for(j in 1:4)
	{
		shotgun_markers_age_rlm_dir_study_grp_level[i,j] <- ifelse(shotgun_markers_age_rlm_qval_study_grp_level[i,j] <= 0.1,3*sign(shotgun_markers_age_rlm_est_study_grp_level[i,j]),ifelse(shotgun_markers_age_rlm_pval_study_grp_level[i,j] <= 0.05,2*sign(shotgun_markers_age_rlm_est_study_grp_level[i,j]),sign(shotgun_markers_age_rlm_est_study_grp_level[i,j])))
	}
}

shotgun_markers_final_age_association <- data.frame(rem_dir=shotgun_markers_age_dir_study_grp_level[shotgun_Markers_age_ordered,1],rlm_dir=shotgun_markers_age_rlm_dir_study_grp_level[shotgun_Markers_age_ordered,1],rem_lower=shotgun_markers_age_lower_study_grp_level[shotgun_Markers_age_ordered,1],rem_med=shotgun_markers_age_est_study_grp_level[shotgun_Markers_age_ordered,1],rem_upper=shotgun_markers_age_upper_study_grp_level[shotgun_Markers_age_ordered,1],studies_decreased=AgeAssociationDirection_shotgun[shotgun_Markers_age_ordered,1]
,studies_increased=AgeAssociationDirection_shotgun[shotgun_Markers_age_ordered,2],q.value=p.adjust(as.numeric(apply(shotgun_markers_age_pval_study_level,1,function(x)(sumlog(x)$p))),method="fdr"),row.names=shotgun_Markers_age_ordered
)

filtered_shotgun_markers_final_age_association <- shotgun_markers_final_age_association[shotgun_markers_final_age_association[,8]<=0.01,]
filtered_shotgun_markers_final_age_association$name <- 1:nrow(filtered_shotgun_markers_final_age_association)
filtered_shotgun_markers_final_age_association$color <- ifelse(rownames(filtered_shotgun_markers_final_age_association) %in% validated_G_Markers_shotgun,"Red","Blue")

AgeEnriched_G_Markers <- c("Streptococcus_vestibularis","Bifidobacterium_dentium","Clostridium_nexile","Clostridium_asparagiforme","Clostridium_symbiosum","Clostridium_citroniae","Clostridium_hathewayi","Enterobacter_cloacae","Eggerthella_lenta")

AgeDepleted_L_Markers <- c("Faecalibacterium_prausnitzii","Eubacterium_hallii","Eubacterium_rectale","Bifidobacterium_longum")


df_shotgun_diversity_uniqueness$AgeEnriched_G_Markers <- apply(apply(df_shotgun_diversity_uniqueness[,AgeEnriched_G_Markers],2,rank_scale),1,sum)

df_shotgun_diversity_uniqueness$AgeDepleted_L_Markers <- apply(apply(df_shotgun_diversity_uniqueness[,AgeDepleted_L_Markers],2,rank_scale),1,sum)

all_control_rows <- c(rownames(df_cmd3_controls_diversity_uniqueness),rownames(df_ag_controls_diversity_uniqueness),rownames(df_he_controls_diversity_uniqueness),rownames(df_nuage_controls_diversity_uniqueness),rownames(df_isc_controls_diversity_uniqueness),rownames(df_odamaki_controls_diversity_uniqueness),rownames(df_logmpie_controls_diversity_uniqueness))
combined_df_sum_stat_species$study_condition <- NA
combined_df_sum_stat_species[intersect(all_control_rows,rownames(combined_df_sum_stat_species)),"study_condition"] = "control"
combined_df_sum_stat_species[setdiff(rownames(combined_df_sum_stat_species),all_control_rows),"study_condition"] = "not_control"

compute_meta_corr_fem <- function(data,var1,var2,grouping_variable,grouping_list)
{
        temp_meta <- data.frame(matrix(NA,length(grouping_list),3))
        colnames(temp_meta) <- c("dataset","ri","ni")
        
        for(i in 1:length(grouping_list))
        {
                group <- grouping_list[i]
                var1_length <- length(which(data[data[,grouping_variable]==group,var1]>0))
                var2_length <- length(which(data[data[,grouping_variable]==group,var2]>0))
                temp_meta[i,1] <- group
                if((var1_length > 0)&(var2_length > 0))
                {
                        temp_meta[i,2] <- cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2],method="kendall")
                        temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
                }
                else
                {
                        temp_meta[i,2] <- 0
                        temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
                }
        }
        temp_meta <- mutate(temp_meta,study_id=grouping_list)
        rownames(temp_meta) <- grouping_list
        temp_meta <- escalc(measure="ZCOR",ri=ri,ni=ni,data=temp_meta)
        res <- rma(yi, vi, data=temp_meta, method="FE")
        res$ids <- rownames(temp_meta)
        res$slabs <- rownames(temp_meta)
        return_list <- list("df_studies"=temp_meta,"model"=res)
        return(return_list)
        
        
}


combined_rem_G2_age <- compute_meta_corr_fem(combined_df_sum_stat_species_clr[(combined_df_sum_stat_species_clr$age>=60),],"G2","age","study_name",selected_elderly_studies)
combined_rem_G1_age <- compute_meta_corr_fem(combined_df_sum_stat_species_clr[(combined_df_sum_stat_species_clr$age>=60),],"G1","age","study_name",selected_elderly_studies)

selected_elderly_studies_wo_NUAGE <- setdiff(selected_elderly_studies,"NUAGE")

print("RLM: G1 Species versus age")
combined_rlm_est_G1_age <- as.data.frame(matrix(0,length(G1Species),length(selected_elderly_studies)))
rownames(combined_rlm_est_G1_age) <- G1Species
colnames(combined_rlm_est_G1_age) <- selected_elderly_studies
combined_rlm_p_val_G1_age <- as.data.frame(matrix(0,length(G1Species),length(selected_elderly_studies)))
rownames(combined_rlm_p_val_G1_age) <- G1Species
colnames(combined_rlm_p_val_G1_age) <- selected_elderly_studies
for(i in 1:length(G1Species))
{
	species_name <- G1Species[i]
	for(j in 1:length(selected_elderly_studies))
	{
		study_name <- selected_elderly_studies[j]
		study_samples <- rownames(combined_df_sum_stat_species_clr[(combined_df_sum_stat_species_clr$study_name == study_name)&(combined_df_sum_stat_species_clr$age>=60),])
		tryCatch(               
					expr = {   
						print(paste0("age, ",species_name))
						df_temp <- as.data.frame(combined_df_sum_stat_species_clr[study_samples,c(species_name,"age")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~age")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						combined_rlm_est_G1_age[i,j] <- summary_temp_rlm$coefficients[2,3]
						combined_rlm_p_val_G1_age[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){   
						print(paste0("age, ",species_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}

combined_rlm_q_val_G1_age <- apply(combined_rlm_p_val_G1_age,2,function(x)(p.adjust(x,method="fdr")))

combined_rlm_dir_G1_age <- as.data.frame(matrix(0,length(G1Species),length(selected_elderly_studies)))
rownames(combined_rlm_dir_G1_age) <- G1Species
colnames(combined_rlm_dir_G1_age) <- selected_elderly_studies
for(i in 1:length(G1Species))
{
	species_name <- G1Species[i]
	for(j in 1:length(selected_elderly_studies_wo_NUAGE))
	{
		combined_rlm_dir_G1_age[i,j] <- ifelse(combined_rlm_q_val_G1_age[i,j]<=0.1,3*sign(combined_rlm_est_G1_age[i,j]),ifelse(combined_rlm_p_val_G1_age[i,j]<=0.05,2*sign(combined_rlm_est_G1_age[i,j]),sign(combined_rlm_est_G1_age[i,j])))
	}
}

mat <- as.matrix(combined_rlm_dir_G1_age)


hmpG1Age <- heatmap.2(t(mat),density="none",trace="none",col=c("skyblue4","skyblue3","white","white","white","orangered3","orangered4"),Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="grey",colsep=1:ncol(t(mat)),rowsep=1:nrow(t(mat)),key=FALSE,margins=c(12,7),lhei=c(0.5,5),lwid=c(0.1,5))

G1AgeAssociationAcrossStudies <- data.frame(increased=apply(hmpG1Age$carpet,1,function(x)(length(x[x>=2]))),decreased=apply(hmpG1Age$carpet,1,function(x)(length(x[x<= -2]))))

print("RLM: G2 Species versus age")
combined_rlm_est_G2_age <- as.data.frame(matrix(0,length(G2Species),length(selected_elderly_studies)))
rownames(combined_rlm_est_G2_age) <- G2Species
colnames(combined_rlm_est_G2_age) <- selected_elderly_studies
combined_rlm_p_val_G2_age <- as.data.frame(matrix(0,length(G2Species),length(selected_elderly_studies)))
rownames(combined_rlm_p_val_G2_age) <- G2Species
colnames(combined_rlm_p_val_G2_age) <- selected_elderly_studies
for(i in 1:length(G2Species))
{
	species_name <- G2Species[i]
	for(j in 1:length(selected_elderly_studies))
	{
		study_name <- selected_elderly_studies[j]
		study_samples <- rownames(combined_df_sum_stat_species_clr[(combined_df_sum_stat_species_clr$study_name == study_name)&(combined_df_sum_stat_species_clr$age >= 60),])
		tryCatch(               
					expr = {   
						print(paste0("age, ",species_name))
						df_temp <- as.data.frame(combined_df_sum_stat_species_clr[study_samples,c(species_name,"age")])
						temp_rlm <- rlm(as.formula(paste0(species_name,"~age")),data=df_temp)
						summary_temp_rlm <- summary(temp_rlm)
						combined_rlm_est_G2_age[i,j] <- summary_temp_rlm$coefficients[2,3]
						combined_rlm_p_val_G2_age[i,j] <- f.robftest(temp_rlm)$p.value
					},
					error = function(e){   
						print(paste0("age, ",species_name))
						print("Error observed. Moving to next")
					},
					finally = {            
						print("finally Executed")
					}
				)
		
	}
}

combined_rlm_q_val_G2_age <- apply(combined_rlm_p_val_G2_age,2,function(x)(p.adjust(x,method="fdr")))

combined_rlm_dir_G2_age <- as.data.frame(matrix(0,length(G2Species),length(selected_elderly_studies)))
rownames(combined_rlm_dir_G2_age) <- G2Species
colnames(combined_rlm_dir_G2_age) <- selected_elderly_studies
for(i in 1:length(G2Species))
{
	species_name <- G2Species[i]
	for(j in 1:length(selected_elderly_studies))
	{
		combined_rlm_dir_G2_age[i,j] <- ifelse(combined_rlm_q_val_G2_age[i,j]<=0.1,3*sign(combined_rlm_est_G2_age[i,j]),ifelse(combined_rlm_p_val_G2_age[i,j]<=0.05,2*sign(combined_rlm_est_G2_age[i,j]),sign(combined_rlm_est_G2_age[i,j])))
	}
}

mat <- as.matrix(combined_rlm_dir_G2_age)

hmpG2Age <- heatmap.2(t(mat),density="none",trace="none",col=c("skyblue4","skyblue3","white","white","white","orangered3","orangered4"),Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="grey",colsep=1:ncol(t(mat)),rowsep=1:nrow(t(mat)),key=FALSE,margins=c(12,7),lhei=c(0.5,5),lwid=c(0.1,5))

G2AgeAssociationAcrossStudies <- data.frame(increased=apply(hmpG2Age$carpet,1,function(x)(length(x[x>=2]))),decreased=apply(hmpG2Age$carpet,1,function(x)(length(x[x<= -2]))))

compute_meta_corr_group <- function(data,feature_list,metadata_var,grouping_var,grouping_list)
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
						temp_res <- compute_meta_corr(data,species_name,metadata_var,grouping_var,grouping_list)
						print(species_name)
						return_out[i,"beta"] <- temp_res$model$beta
						return_out[i,"pval"] <- temp_res$model$pval
						return_out[i,"ci.ub"] <- temp_res$model$ci.ub
						return_out[i,"ci.lb"] <- temp_res$model$ci.lb
						return_out[i,"tau2"] <- temp_res$model$tau2
						return_out[i,"QE"] <- temp_res$model$QE
						return_out[i,"QEp"] <- temp_res$model$QEp
						return_out[i,"consistency"] <- length(which(sign(temp_res$df_studies[temp_res$df_studies$ri!=0,"ri"])==sign(as.numeric(temp_res$model$beta))))/length(temp_res$df_studies[temp_res$df_studies$ri!=0,"ri"])
						
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
	return_list <- list("model" = temp_res$model,"df_studies"=temp_res$df_studies)
	#return(return_list)
	return(return_out)
}


combined_rem_G1Species_age <- compute_meta_lm_group(combined_df_sum_stat_species_clr[combined_df_sum_stat_species_clr$age>=60,],G1Species,"age","study_name",selected_elderly_studies)
combined_rem_G1Species_age$species_group <- "G1"

combined_rem_G2Species_age <- compute_meta_lm_group(combined_df_sum_stat_species_clr[combined_df_sum_stat_species_clr$age>=60,],G2Species,"age","study_name",selected_elderly_studies_wo_NUAGE)
combined_rem_G2Species_age$species_group <- "G2"

OtherSpecies <- setdiff(rownames(species_groupings),c(G1Species,G2Species))
combined_df_sum_stat_species_clr$Other <- rowMeans(apply(combined_df_sum_stat_species_clr[,OtherSpecies],2,range_scale))

combined_rem_Other_age <- compute_meta_corr_fem(combined_df_sum_stat_species_clr[(combined_df_sum_stat_species_clr$age>=60),],"Other","age","study_name",setdiff(selected_elderly_studies,"NUAGE"))

combined_rem_OtherSpecies_age <- compute_meta_lm_group(combined_df_sum_stat_species_clr[combined_df_sum_stat_species_clr$age>=60,],OtherSpecies,"age","study_name",selected_elderly_studies_wo_NUAGE)
combined_rem_OtherSpecies_age$species_group <- "Other"

combined_rem_bothSpecies_age <- as.data.frame(rbind(combined_rem_G1Species_age,combined_rem_G2Species_age))

mat_rem <- combined_rem_bothSpecies_age[combined_rem_bothSpecies_age$consistency>=0.60,]
mat_rem[,"qval"] <- ifelse(mat_rem[,"qval"]<=0.01,0.01,mat_rem[,"qval"])

#ggplot(mat_rem,aes(x=beta,y=-log(qval,10)))+geom_point(col=ifelse(mat_rem$species_group %in% "G1","green","red"),size=ifelse(abs(mat_rem$dir)>=2,5,2),pch=ifelse(abs(mat_rem$dir)>2,15,19))+geom_text_repel(label=ifelse(abs(mat_rem$dir)>=2,rownames(mat_rem),""),size=ifelse(abs(mat_rem$dir)>=3,5,ifelse(abs(mat_rem$dir)>=2,3,1)),max.overlaps=30)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15))+geom_hline(yintercept=0)+geom_vline(xintercept=0)

combined_df_sum_stat_species_clr$study_condition <- NA
combined_df_sum_stat_species_clr[intersect(rownames(combined_df_sum_stat_species_clr),rownames(df_cmd3_disease_investigation)),"study_condition"] <- df_cmd3_disease_investigation[intersect(rownames(combined_df_sum_stat_species_clr),rownames(df_cmd3_disease_investigation)),"study_condition"]
combined_df_sum_stat_species_clr[intersect(rownames(combined_df_sum_stat_species_clr),rownames(df_ag_controls_diversity_uniqueness)),"study_condition"] <- "control"
combined_df_sum_stat_species_clr[intersect(rownames(combined_df_sum_stat_species_clr),rownames(df_he_controls_diversity_uniqueness)),"study_condition"] <- "control"
combined_df_sum_stat_species_clr[intersect(rownames(combined_df_sum_stat_species_clr),rownames(df_logmpie_controls_diversity_uniqueness)),"study_condition"] <- "control"
combined_df_sum_stat_species_clr[intersect(rownames(combined_df_sum_stat_species_clr),rownames(df_odamaki_controls_diversity_uniqueness)),"study_condition"] <- "control"
combined_df_sum_stat_species_clr[intersect(rownames(combined_df_sum_stat_species_clr),rownames(df_nuage_controls_diversity_uniqueness)),"study_condition"] <- "control"
combined_df_sum_stat_species_clr[intersect(rownames(combined_df_sum_stat_species_clr),rownames(df_isc_controls_diversity_uniqueness)),"study_condition"] <- "control"

combined_rem_OtherSpecies_controls_age <- compute_meta_lm_group(combined_df_sum_stat_species_clr[(combined_df_sum_stat_species_clr$age >= 60)&(combined_df_sum_stat_species_clr$study_condition=="control"),],OtherSpecies,"age","study_name",selected_elderly_studies)
combined_rem_G1Species_controls_age <- compute_meta_lm_group(combined_df_sum_stat_species_clr[(combined_df_sum_stat_species_clr$age >= 60)&(combined_df_sum_stat_species_clr$study_condition=="control"),],G1Species,"age","study_name",selected_elderly_studies)
combined_rem_G2Species_controls_age <- compute_meta_lm_group(combined_df_sum_stat_species_clr[(combined_df_sum_stat_species_clr$age >= 60)&(combined_df_sum_stat_species_clr$study_condition=="control"),],G2Species,"age","study_name",selected_elderly_studies)

