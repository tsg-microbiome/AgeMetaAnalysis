DiseaseNumbers <- table(curatedMetagenomicData_metadata[rownames(df_cmd_genus_uniqueness),"study_condition"])
ElderlyDiseaseNumbers <- table(df_cmd_genus_uniqueness[df_cmd_genus_uniqueness$age>=60,"study_condition_2"])
ElderlyDiseases <- names(which(table(df_cmd_genus_uniqueness[df_cmd_genus_uniqueness$age>=60,"study_condition_2"])>=5))
ElderlyDiseases <- setdiff(ElderlyDiseases,"control")
DiseaseAssociationMatrix <- as.data.frame(matrix(NA,length(ElderlyDiseases),33))
rownames(DiseaseAssociationMatrix) <- ElderlyDiseases
colnames(DiseaseAssociationMatrix) <- c("Uniqueness_Estimate","Uniqueness_T","Uniqueness_P","Diversity_Estimate","Diversity_T","Diversity_P","AllG_Estimate","AllG_T","AllG_P","G1_Estimate","G1_T","G1_P","G2_Estimate","G2_T","G2_P","G3_Estimate","G3_T","G3_P","L1_Estimate","L1_T","L1_P","L2_Estimate","L2_T","L2_P","L3_Estimate","L3_T","L3_P","AllL_Estimate","AllL_T","AllL_P","GtoL_Estimate","GToL_T","GToL_P")
for(i in 1:length(ElderlyDiseases))
{
	Disease <- ElderlyDiseases[i]
	print(Disease)
	DiseaseStudyNames <- unique(df_cmd_genus_uniqueness[df_cmd_genus_uniqueness$study_condition_2 == Disease,"study_name"])
	assign(paste0(Disease,"_study_names"),DiseaseStudyNames)
	controls <- rownames(df_cmd_genus_uniqueness[(df_cmd_genus_uniqueness$study_name %in% DiseaseStudyNames)&(df_cmd_genus_uniqueness$study_condition_2=="control")&(df_cmd_genus_uniqueness$age>=60),])
	diseased <- rownames(df_cmd_genus_uniqueness[(df_cmd_genus_uniqueness$study_name %in% DiseaseStudyNames)&(df_cmd_genus_uniqueness$study_condition_2==Disease)&(df_cmd_genus_uniqueness$age>=60),])
	summary_fit <- summary(lm(bray_uniqueness~diversity+as.factor(ifelse(study_condition_2=="control",0,1)),data=df_cmd_genus_uniqueness[c(controls,diseased),]))
	DiseaseAssociationMatrix[i,1] <- summary_fit$coefficients[3,1]	
	DiseaseAssociationMatrix[i,2] <- summary_fit$coefficients[3,3]
	DiseaseAssociationMatrix[i,3] <- summary_fit$coefficients[3,4]
	summary_fit <- summary(lm(diversity~as.factor(ifelse(study_condition_2=="control",0,1)),data=df_cmd_genus_uniqueness[c(controls,diseased),]))
	DiseaseAssociationMatrix[i,4] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,5] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,6] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(AllG~as.factor(ifelse(study_condition_2=="control",0,1)),data=df_cmd_genus_uniqueness[c(controls,diseased),]))
	DiseaseAssociationMatrix[i,7] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,8] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,9] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(G1~as.factor(ifelse(study_condition_2=="control",0,1)),data=df_cmd_genus_uniqueness[c(controls,diseased),]))
	DiseaseAssociationMatrix[i,10] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,11] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,12] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(G2~as.factor(ifelse(study_condition_2=="control",0,1)),data=df_cmd_genus_uniqueness[c(controls,diseased),]))
	DiseaseAssociationMatrix[i,13] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,14] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,15] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(G3~as.factor(ifelse(study_condition_2=="control",0,1)),data=df_cmd_genus_uniqueness[c(controls,diseased),]))
	DiseaseAssociationMatrix[i,16] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,17] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,18] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(L1~as.factor(ifelse(study_condition_2=="control",0,1)),data=df_cmd_genus_uniqueness[c(controls,diseased),]))
	DiseaseAssociationMatrix[i,19] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,20] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,21] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(L2~as.factor(ifelse(study_condition_2=="control",0,1)),data=df_cmd_genus_uniqueness[c(controls,diseased),]))
	DiseaseAssociationMatrix[i,22] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,23] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,24] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(L3~as.factor(ifelse(study_condition_2=="control",0,1)),data=df_cmd_genus_uniqueness[c(controls,diseased),]))
	DiseaseAssociationMatrix[i,25] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,26] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,27] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(AllL~as.factor(ifelse(study_condition_2=="control",0,1)),data=df_cmd_genus_uniqueness[c(controls,diseased),]))
	DiseaseAssociationMatrix[i,28] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,29] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,30] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(G_to_L_Ratio~as.factor(ifelse(study_condition_2=="control",0,1)),data=df_cmd_genus_uniqueness[c(controls,diseased),]))
	DiseaseAssociationMatrix[i,31] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,32] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,33] <- summary_fit$coefficients[2,4]
}

DiseaseAssociationSummary <- matrix(NA,nrow(DiseaseAssociationMatrix),5)
rownames(DiseaseAssociationSummary) <- rownames(DiseaseAssociationMatrix)
colnames(DiseaseAssociationSummary) <-c("Uniqueness","Diversity","Pathobiont","Commensal","PathobiontCommensalRatio")

DiseaseAssociationSummary[,1] <- ifelse(DiseaseAssociationMatrix[,3] <= 0.001,5*sign(DiseaseAssociationMatrix[,1]),ifelse(DiseaseAssociationMatrix[,3] <= 0.01,4*sign(DiseaseAssociationMatrix[,1]),ifelse(DiseaseAssociationMatrix[,3] <= 0.05,3*sign(DiseaseAssociationMatrix[,1]),ifelse(DiseaseAssociationMatrix[,3] <= 0.10,2*sign(DiseaseAssociationMatrix[,1]),sign(DiseaseAssociationMatrix[,1])))))
DiseaseAssociationSummary[,2] <- ifelse(DiseaseAssociationMatrix[,6] <= 0.001,5*sign(DiseaseAssociationMatrix[,4]),ifelse(DiseaseAssociationMatrix[,6] <= 0.01,4*sign(DiseaseAssociationMatrix[,4]),ifelse(DiseaseAssociationMatrix[,6] <= 0.05,3*sign(DiseaseAssociationMatrix[,4]),ifelse(DiseaseAssociationMatrix[,6] <= 0.10,2*sign(DiseaseAssociationMatrix[,4]),sign(DiseaseAssociationMatrix[,4])))))
DiseaseAssociationSummary[,3] <- ifelse(DiseaseAssociationMatrix[,9] <= 0.001,5*sign(DiseaseAssociationMatrix[,7]),ifelse(DiseaseAssociationMatrix[,9] <= 0.01,4*sign(DiseaseAssociationMatrix[,7]),ifelse(DiseaseAssociationMatrix[,9] <= 0.05,3*sign(DiseaseAssociationMatrix[,7]),ifelse(DiseaseAssociationMatrix[,9] <= 0.10,2*sign(DiseaseAssociationMatrix[,7]),sign(DiseaseAssociationMatrix[,7])))))
DiseaseAssociationSummary[,4] <- ifelse(DiseaseAssociationMatrix[,30] <= 0.001,5*sign(DiseaseAssociationMatrix[,28]),ifelse(DiseaseAssociationMatrix[,30] <= 0.01,4*sign(DiseaseAssociationMatrix[,28]),ifelse(DiseaseAssociationMatrix[,30] <= 0.05,3*sign(DiseaseAssociationMatrix[,28]),ifelse(DiseaseAssociationMatrix[,30] <= 0.10,2*sign(DiseaseAssociationMatrix[,28]),sign(DiseaseAssociationMatrix[,28])))))
DiseaseAssociationSummary[,5] <- ifelse(DiseaseAssociationMatrix[,33] <= 0.001,5*sign(DiseaseAssociationMatrix[,31]),ifelse(DiseaseAssociationMatrix[,33] <= 0.01,4*sign(DiseaseAssociationMatrix[,31]),ifelse(DiseaseAssociationMatrix[,33] <= 0.05,3*sign(DiseaseAssociationMatrix[,31]),ifelse(DiseaseAssociationMatrix[,33] <= 0.10,2*sign(DiseaseAssociationMatrix[,31]),sign(DiseaseAssociationMatrix[,31])))))




