IndexList <- c("FIM","Barthel_Score","MMSE","Charlson_Comorbidity","Geriatric_Depression")

DiseaseAssociationMatrix <- as.data.frame(matrix(NA,length(IndexList),30))
rownames(DiseaseAssociationMatrix) <- c("Frailty_FIM","Frailty_BarthelScore","CognitiveImpairment_MMSE","Charlson_Comorbidity","Geriatric_Depression")
colnames(DiseaseAssociationMatrix) <- c("Uniqueness_Estimate","Uniqueness_T","Uniqueness_P","Diversity_Estimate","Diversity_T","Diversity_P","AllG_Estimate","AllG_T","AllG_P","G1_Estimate","G1_T","G1_P","G2_Estimate","G2_T","G2_P","G3_Estimate","G3_T","G3_P","L1_Estimate","L1_T","L1_P","L2_Estimate","L2_T","L2_P","L3_Estimate","L3_T","L3_P","AllL_Estimate","AllL_T","AllL_P")

for(i in 1:length(IndexList))
{
	Index <- IndexList[i]
	summary_fit <- summary(lm(as.formula(paste0("bray_uniqueness~",Index)),data=df_em_genus_uniqueness,na.action=na.omit))
	DiseaseAssociationMatrix[i,1] <- -summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,2] <- -summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,3] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("diversity~",Index)),data=df_em_genus_uniqueness,na.action=na.omit))
	DiseaseAssociationMatrix[i,4] <- -summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,5] <- -summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,6] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("AllG~",Index)),data=df_em_genus_uniqueness,na.action=na.omit))
	DiseaseAssociationMatrix[i,7] <- -summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,8] <- -summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,9] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("G1~",Index)),data=df_em_genus_uniqueness,na.action=na.omit))
	DiseaseAssociationMatrix[i,10] <- -summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,11] <- -summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,12] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("G2~",Index)),data=df_em_genus_uniqueness,na.action=na.omit))
	DiseaseAssociationMatrix[i,13] <- -summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,14] <- -summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,15] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("G3~",Index)),data=df_em_genus_uniqueness,na.action=na.omit))
	DiseaseAssociationMatrix[i,16] <- -summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,17] <- -summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,18] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("L1~",Index)),data=df_em_genus_uniqueness,na.action=na.omit))
	DiseaseAssociationMatrix[i,19] <- -summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,20] <- -summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,21] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("L2~",Index)),data=df_em_genus_uniqueness,na.action=na.omit))
	DiseaseAssociationMatrix[i,22] <- -summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,23] <- -summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,24] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("L3~",Index)),data=df_em_genus_uniqueness,na.action=na.omit))
	DiseaseAssociationMatrix[i,25] <- -summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,26] <- -summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,27] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("AllL~",Index)),data=df_em_genus_uniqueness,na.action=na.omit))
	DiseaseAssociationMatrix[i,28] <- -summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,29] <- -summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,30] <- summary_fit$coefficients[2,4]
}

DiseaseAssociationSummary <- matrix(NA,nrow(DiseaseAssociationMatrix),4)
rownames(DiseaseAssociationSummary) <- rownames(DiseaseAssociationMatrix)
colnames(DiseaseAssociationSummary) <-c("Uniqueness","Diversity","Pathobiont","Commensal")

DiseaseAssociationSummary[,1] <- ifelse(DiseaseAssociationMatrix[,3] <= 0.001,5*sign(DiseaseAssociationMatrix[,1]),ifelse(DiseaseAssociationMatrix[,3] <= 0.01,4*sign(DiseaseAssociationMatrix[,1]),ifelse(DiseaseAssociationMatrix[,3] <= 0.05,3*sign(DiseaseAssociationMatrix[,1]),ifelse(DiseaseAssociationMatrix[,3] <= 0.10,2*sign(DiseaseAssociationMatrix[,1]),sign(DiseaseAssociationMatrix[,1])))))
DiseaseAssociationSummary[,2] <- ifelse(DiseaseAssociationMatrix[,6] <= 0.001,5*sign(DiseaseAssociationMatrix[,4]),ifelse(DiseaseAssociationMatrix[,6] <= 0.01,4*sign(DiseaseAssociationMatrix[,4]),ifelse(DiseaseAssociationMatrix[,6] <= 0.05,3*sign(DiseaseAssociationMatrix[,4]),ifelse(DiseaseAssociationMatrix[,6] <= 0.10,2*sign(DiseaseAssociationMatrix[,4]),sign(DiseaseAssociationMatrix[,4])))))
DiseaseAssociationSummary[,3] <- ifelse(DiseaseAssociationMatrix[,9] <= 0.001,5*sign(DiseaseAssociationMatrix[,7]),ifelse(DiseaseAssociationMatrix[,9] <= 0.01,4*sign(DiseaseAssociationMatrix[,7]),ifelse(DiseaseAssociationMatrix[,9] <= 0.05,3*sign(DiseaseAssociationMatrix[,7]),ifelse(DiseaseAssociationMatrix[,9] <= 0.10,2*sign(DiseaseAssociationMatrix[,7]),sign(DiseaseAssociationMatrix[,7])))))
DiseaseAssociationSummary[,4] <- ifelse(DiseaseAssociationMatrix[,30] <= 0.001,5*sign(DiseaseAssociationMatrix[,28]),ifelse(DiseaseAssociationMatrix[,30] <= 0.01,4*sign(DiseaseAssociationMatrix[,28]),ifelse(DiseaseAssociationMatrix[,30] <= 0.05,3*sign(DiseaseAssociationMatrix[,28]),ifelse(DiseaseAssociationMatrix[,30] <= 0.10,2*sign(DiseaseAssociationMatrix[,28]),sign(DiseaseAssociationMatrix[,28])))))