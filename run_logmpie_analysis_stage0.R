#This is the first sub-pipeline for the logmpie Data repository.
#The sub-pipeline processes the raw species, genus and metadata profiles of the logmpie data repository and generates and saves the processed species profile (logmpie_select_age_final_species), processed genus profile (logmpie_select_age_final_genus) and the processed metadata (logmpie_select_age_final_metadata) into the workspace logmpie_age_analysis.RData
#The logmpie_age_analysis.RData has already been uploaded to this github.
logmpie_species <- read.table("C:\\Projects\\ELDERMET\\NatureAgingRevision\\LogMPie\\logmpie_species.txt",sep="\t",row.names=1,header=TRUE)
logmpie_genus <- read.table("C:\\Projects\\ELDERMET\\NatureAgingRevision\\LogMPie\\logmpie_genus.txt",sep="\t",row.names=1,header=TRUE)
logmpie_metadata <- read.table("C:\\Projects\\ELDERMET\\NatureAgingRevision\\LogMPie\\logmpie_metadata.txt",sep="\t",row.names=1,header=TRUE)

logmpie_select_age_rows <- intersect(rownames(logmpie_species),rownames(logmpie_metadata[!is.na(logmpie_metadata$Age),]))

logmpie_select_age_final_species <- logmpie_species[logmpie_select_age_rows,]
logmpie_select_age_final_species <- t(apply(logmpie_select_age_final_species,1,function(x)(x/sum(x))))
logmpie_select_age_final_genus <- logmpie_genus[logmpie_select_age_rows,]
logmpie_select_age_final_genus <- t(apply(logmpie_select_age_final_genus,1,function(x)(x/sum(x))))
logmpie_select_age_final_metadata <- logmpie_metadata[logmpie_select_age_rows,]

save(logmpie_select_age_final_species,logmpie_select_age_final_genus,logmpie_select_age_final_metadata,file="C:\\Projects\\ELDERMET\\NatureAgingRevision\\LogMPie\\logmpie_age_analysis.RData")
rm(list=ls())





