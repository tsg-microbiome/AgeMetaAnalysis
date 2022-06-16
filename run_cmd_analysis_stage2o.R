#The sub-pipeline creates the processed dataframes for the disease-microbiome studies for the studies in the curatedMetagenomicData3 (CMD) repository.
#The sub-pipeline generates three data frames pertaining to the studies in the CMD studies in Table 1A and 1B.
#The three data frames are: cmd3_select_disease_final_species (containing the species profile), cmd3_select_disease_final_genus (containing the genus profile) and cmd3_select_disease_metadata (containing the metadata variables)

load("C:\\Projects\\curatedMetagenomicData\\curatedMetagenomicData3\\cmd3.RData")
load("C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\QinJ_Table.RData")
cmd3_select_metadata[rownames(QinJ_Table),"age"] <- QinJ_Table[rownames(QinJ_Table),"age"]
cmd3_select_disease_metadata <- cmd3_select_metadata[(!is.na(cmd3_select_metadata$age))&(cmd3_select_metadata$age>=18),]
cmd3_select_disease_metadata <- cmd3_select_disease_metadata[(cmd3_select_disease_metadata$study_name != "RaymondF_2016")&(cmd3_select_disease_metadata$study_condition != "carcinoma_surgery_history"),]
disease_control_matched_studies <- intersect(unique(cmd3_select_disease_metadata[(!is.na(cmd3_select_disease_metadata$age))&(cmd3_select_disease_metadata$study_condition!="control"),"study_name"]),unique(cmd3_select_disease_metadata[(!is.na(cmd3_select_disease_metadata$age))&(cmd3_select_disease_metadata$study_condition=="control"),"study_name"]))
disease_control_matched_studies <- disease_control_matched_studies[!is.na(disease_control_matched_studies)]

cmd3_select_disease_studies_details <- as.data.frame(matrix(NA,length(disease_control_matched_studies),4))
rownames(cmd3_select_disease_studies_details) <- disease_control_matched_studies
colnames(cmd3_select_disease_studies_details) <- c("minimum_age","maximum_age","study_conditions","nationalities")

for(i in 1:length(disease_control_matched_studies))
{
	study_name = disease_control_matched_studies[i]
	cmd3_select_disease_studies_details[i,1] <- min(cmd3_select_disease_metadata[cmd3_select_disease_metadata$study_name == study_name,"age"])
	cmd3_select_disease_studies_details[i,2] <- max(cmd3_select_disease_metadata[cmd3_select_disease_metadata$study_name == study_name,"age"])
	cmd3_select_disease_studies_details[i,3] <- paste0(unique(cmd3_select_disease_metadata[cmd3_select_disease_metadata$study_name == study_name,"study_condition"]),collapse=",")
	cmd3_select_disease_studies_details[i,4] <- paste0(unique(cmd3_select_disease_metadata[cmd3_select_disease_metadata$study_name == study_name,"country"]),collapse=",")	
}

cmd3_select_disease_final_species <- cmd3_stool_species_profile[rownames(cmd3_select_disease_metadata),]
cmd3_select_disease_final_genus <- cmd3_stool_genus_profile[rownames(cmd3_select_disease_metadata),]
cmd3_select_disease_final_metadata <- cmd3_select_disease_metadata

save(cmd3_select_disease_studies_details,cmd3_select_disease_final_species,cmd3_select_disease_final_genus,cmd3_select_disease_final_metadata,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\CMD\\cmd3_disease_analysis.RData")
