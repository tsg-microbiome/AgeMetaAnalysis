IndexList <- c("cspraxis","hgtdommean","gstfsttime","FriedScore","hsCRP","IL17","microbiome_scores")

df_NUAGE_genus_uniqueness_1$Inv_cspraxis <- (-1)* CombinedDataFinal[rownames(df_NUAGE_genus_uniqueness_1),"cspraxis"]
df_NUAGE_genus_uniqueness_1$Inv_hgtdommean <- (-1)* CombinedDataFinal[rownames(df_NUAGE_genus_uniqueness_1),"hgtdommean"]
#df_NUAGE_genus_uniqueness_1$Inv_microbiome_scores <- (-1)* CombinedDataFinal[rownames(df_NUAGE_genus_uniqueness_1),"microbiome_scores"]
df_NUAGE_genus_uniqueness_1$gstfsttime <- CombinedDataFinal[rownames(df_NUAGE_genus_uniqueness_1),"gstfsttime"]
df_NUAGE_genus_uniqueness_1$hsCRP <- rank_scale(CombinedDataFinal[rownames(df_NUAGE_genus_uniqueness_1),"hsCRP"])
df_NUAGE_genus_uniqueness_1$IL17 <- rank_scale(CombinedDataFinal[rownames(df_NUAGE_genus_uniqueness_1),"IL17"])
df_NUAGE_genus_uniqueness_1$FriedScore <- rank_scale(CombinedDataFinal[rownames(df_NUAGE_genus_uniqueness_1),"FriedScore"])
df_NUAGE_genus_uniqueness_1$Inv_BabcokMemory <- -CombinedData[rownames(df_NUAGE_genus_uniqueness_1),"bsrtscore"]
df_NUAGE_genus_uniqueness_1$Inv_VerbalFluency <- -CombinedData[rownames(df_NUAGE_genus_uniqueness_1),"bostontot"]
df_NUAGE_genus_uniqueness_1$GDS <- CombinedData[rownames(df_NUAGE_genus_uniqueness_1),"gdsscore"]
df_NUAGE_genus_uniqueness_1$Inv_MMSE <- -CombinedData[rownames(df_NUAGE_genus_uniqueness_1),"mmtotalsc"]

FinalIndexList <- c("Inv_hgtdommean","gstfsttime","FriedScore","Inv_cspraxis","Inv_BabcokMemory","Inv_VerbalFluency","GDS","Inv_MMSE","hsCRP","IL17")

DiseaseAssociationMatrix <- as.data.frame(matrix(NA,length(FinalIndexList),30))
rownames(DiseaseAssociationMatrix) <- FinalIndexList
colnames(DiseaseAssociationMatrix) <- c("Uniqueness_Estimate","Uniqueness_T","Uniqueness_P","Diversity_Estimate","Diversity_T","Diversity_P","AllG_Estimate","AllG_T","AllG_P","G1_Estimate","G1_T","G1_P","G2_Estimate","G2_T","G2_P","G3_Estimate","G3_T","G3_P","L1_Estimate","L1_T","L1_P","L2_Estimate","L2_T","L2_P","L3_Estimate","L3_T","L3_P","AllL_Estimate","AllL_T","AllL_P")

for(i in 1:length(FinalIndexList))
{
	Index <- FinalIndexList[i]
	summary_fit <- summary(lm(as.formula(paste0("bray_uniqueness~",Index)),data=df_NUAGE_genus_uniqueness_1,na.action=na.omit))
	DiseaseAssociationMatrix[i,1] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,2] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,3] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("diversity~",Index)),data=df_NUAGE_genus_uniqueness_1,na.action=na.omit))
	DiseaseAssociationMatrix[i,4] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,5] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,6] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("AllG~",Index)),data=df_NUAGE_genus_uniqueness_1,na.action=na.omit))
	DiseaseAssociationMatrix[i,7] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,8] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,9] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("G1~",Index)),data=df_NUAGE_genus_uniqueness_1,na.action=na.omit))
	DiseaseAssociationMatrix[i,10] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,11] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,12] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("G2~",Index)),data=df_NUAGE_genus_uniqueness_1,na.action=na.omit))
	DiseaseAssociationMatrix[i,13] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,14] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,15] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("G3~",Index)),data=df_NUAGE_genus_uniqueness_1,na.action=na.omit))
	DiseaseAssociationMatrix[i,16] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,17] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,18] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("L1~",Index)),data=df_NUAGE_genus_uniqueness_1,na.action=na.omit))
	DiseaseAssociationMatrix[i,19] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,20] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,21] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("L2~",Index)),data=df_NUAGE_genus_uniqueness_1,na.action=na.omit))
	DiseaseAssociationMatrix[i,22] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,23] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,24] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("L3~",Index)),data=df_NUAGE_genus_uniqueness_1,na.action=na.omit))
	DiseaseAssociationMatrix[i,25] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,26] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,27] <- summary_fit$coefficients[2,4]
	summary_fit <- summary(lm(as.formula(paste0("AllL~",Index)),data=df_NUAGE_genus_uniqueness_1,na.action=na.omit))
	DiseaseAssociationMatrix[i,28] <- summary_fit$coefficients[2,1]	
	DiseaseAssociationMatrix[i,29] <- summary_fit$coefficients[2,3]
	DiseaseAssociationMatrix[i,30] <- summary_fit$coefficients[2,4]
}

DiseaseAssociationSummary <- matrix(NA,nrow(DiseaseAssociationMatrix),4)
rownames(DiseaseAssociationSummary) <- rownames(DiseaseAssociationMatrix)
colnames(DiseaseAssociationSummary) <-c("Uniqueness","Diversity","Pathobiont","Commensal")

DiseaseAssociationSummary[,1] <- ifelse(DiseaseAssociationMatrix[,3] <= 0.001,5*sign(DiseaseAssociationMatrix[,1]),ifelse(DiseaseAssociationMatrix[,3] <= 0.01,4*sign(DiseaseAssociationMatrix[,1]),ifelse(DiseaseAssociationMatrix[,3] <= 0.05,3*sign(DiseaseAssociationMatrix[,1]),ifelse(DiseaseAssociationMatrix[,3] <= 0.10,2*sign(DiseaseAssociationMatrix[,1]),sign(DiseaseAssociationMatrix[,1])))))
DiseaseAssociationSummary[,2] <- ifelse(DiseaseAssociationMatrix[,6] <= 0.001,5*sign(DiseaseAssociationMatrix[,4]),ifelse(DiseaseAssociationMatrix[,6] <= 0.01,4*sign(DiseaseAssociationMatrix[,4]),ifelse(DiseaseAssociationMatrix[,6] <= 0.05,3*sign(DiseaseAssociationMatrix[,4]),ifelse(DiseaseAssociationMatrix[,6] <= 0.10,2*sign(DiseaseAssociationMatrix[,4]),sign(DiseaseAssociationMatrix[,4])))))
DiseaseAssociationSummary[,3] <- ifelse(DiseaseAssociationMatrix[,9] <= 0.001,5*sign(DiseaseAssociationMatrix[,7]),ifelse(DiseaseAssociationMatrix[,9] <= 0.01,4*sign(DiseaseAssociationMatrix[,7]),ifelse(DiseaseAssociationMatrix[,9] <= 0.05,3*sign(DiseaseAssociationMatrix[,7]),ifelse(DiseaseAssociationMatrix[,9] <= 0.10,2*sign(DiseaseAssociationMatrix[,7]),sign(DiseaseAssociationMatrix[,7])))))
DiseaseAssociationSummary[,4] <- ifelse(DiseaseAssociationMatrix[,30] <= 0.001,5*sign(DiseaseAssociationMatrix[,28]),ifelse(DiseaseAssociationMatrix[,30] <= 0.01,4*sign(DiseaseAssociationMatrix[,28]),ifelse(DiseaseAssociationMatrix[,30] <= 0.05,3*sign(DiseaseAssociationMatrix[,28]),ifelse(DiseaseAssociationMatrix[,30] <= 0.10,2*sign(DiseaseAssociationMatrix[,28]),sign(DiseaseAssociationMatrix[,28])))))