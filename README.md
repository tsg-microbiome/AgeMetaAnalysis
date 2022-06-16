# AgeMetaAnalysis
The repository contains the codes and the initial starting data corresponding to the meta-analysis investigating interactions between age and different gut microbiome properties across 21000 gut microbiomes from seven major microbiome data repositories. The repository contains codes specific for each repository as well as codes pertaining to the combined analysis of data from multiple repositories. The seven repositories are: 
1. CMD3 (curatedMetagenomicData3) (codes and data pertaining to this repository are tagged with the key word cmd or cmd3).
2. AG (or AGP) (AmericanGut) (codes and data pertaining to this repository are tagged with the key word ag).
3. NU-AGE (codes and data pertaining to this repository are tagged with the key word nuage).
4. Irish Shotgun Cohorts (ISC) (codes and data pertaining to this repository are tagged with the key word isc). The sub-part of this cohort that contains only analysis or data originating from the ELDERMET project are tagged with the keyword em
5. Odamaki (codes and data pertaining to this repository are tagged with the key word odamaki)
6. He et al cohort (codes and data pertaining to this repository are tagged with the key word he)
7. LogMPie cohort (codes and data pertaining to this repository are tagged with the key word logmpie)

The different workspaces added here contain the starting data for the seven data repositories.
The details of the different workspaces are:

a. EM_Pathabundance.RData: The workspace contains the relative abundances of the different MetaCyc annotation based pathways for the different samples belonging to the ELDERMET sub-cohort of the Irish Shotgun Dataset data repository.

b. ExM_Pathabundance.RData: The workspace contains the relative abundances of the different MetaCyc annotation based pathways for the different samples belonging to the ExerciseMET sub-cohort of the Irish Shotgun Dataset data repository.

c. IBS4DPathabundance.RData: The workspace contains the relative abundances of the different MetaCyc annotation based pathways for the different samples belonging to the Jeffery et al IBS plus Control sub-cohort of the Irish Shotgun Dataset data repository.

d. IBS4DPathAbundanceAnalysis.RData: This workspace contains the unprocessed data corresponding to the previous workspace described here in.

e. ISC_Metadata.RData: The workspace contains the select descriptors (age, study condition of the subject, study name) of all the 464 samples collated as part of the Irish Shotgun Cohort data repository.

f. isc_age_analysis.RData: The workspace contains three data frames that act as the starting data for all taxonomy linked analyses corresponding to the Irish Shotgun Dataset (in the current study). The isc_select_age_final_genus containing the genus-level profiles; isc_select_age_final_species containing the species-level profiles and; isc_select_age_final_metadata containing all the relevant metadata.

g. em_aging_data.RData: The workspace contains the data frame df_em_aging_data that contains five different measures used in the current study used for measuring the extent of unhealthy aging in the current study. The features are Charlsson Comorbidity, Geriatric Depression, Inverse of Functional Independence Measure (FIM) (that is computed as -1 * FIM); inverse of Barthel Score (that is computed as -1 * Barthel Score); inverse of Mini Mental State Examination (MMSE) score (that is computed as -1 * MMSE).

h. ag_age_analysis.RData: The workspace contains three data frames that act as the starting data for all taxonomy linked analyses corresponding to the American Gut Dataset (in the current study). The ag_select_age_final_genus containing the genus-level profiles; ag_select_age_final_species containing the species-level profiles and; ag_select_age_final_metadata containing all the relevant metadata. It also additionally contains some tables corresponding to the subject IDs of the individual samples.

i. ag_select_rows.RData: The workspace contains the list of 3812 subjects from US and UK included from the American Gut Data Repository.

j. cmd3_age_analysis.RData: The workspace contains the data frames that act as the starting data for all taxonomy linked age-specific analyses corresponding to the curatedMetagenomicData3 (in the current study). cmd3_final_set_sample_details contains the initial list of studies selected for this purpose. The cmd3_select_age_final_genus contains the genus-level profiles; cmd3_select_age_final_species contains the species-level profiles and; cmd3_select_age_final_metadata contains all the relevant metadata. 

k. cmd3_disease_analysis.RData: The workspace contains the data frames that act as the starting data for all taxonomy linked disease versus control-specific analyses corresponding to the curatedMetagenomicData3 (in the current study). It includes the studies in Table 1A as well as those Table 1B. cmd3_select_disease_studies_details contains the initial list of studies selected for this purpose. The cmd3_select_disease_final_genus contains the genus-level profiles; cmd3_select_disease_final_species contains the species-level profiles and; cmd3_select_disease_final_metadata contains all the relevant metadata. 

l. QinJ_Table.RData: The workspace contains the list of samples in the Qin J et al 2012 dataset in curatedMetagenomicData3.

m. he_age_analysis.RData: The workspace contains three data frames that act as the starting data for all taxonomy linked analyses corresponding to the He et al Dataset (in the current study). The he_select_age_final_genus contains the genus-level profiles; he_select_age_final_species contains the species-level profiles and; he_select_age_final_metadata contains all the relevant metadata. 

n. odamaki_age_analysis.RData: The workspace contains three data frames that act as the starting data for all taxonomy linked analyses corresponding to the Odamaki et al Dataset (in the current study). The odamaki_select_age_final_genus contains the genus-level profiles; odamaki_select_age_final_species contains the species-level profiles and; odamaki_select_age_final_metadata contains all the relevant metadata. 

o. logmpie_age_analysis.RData: The workspace contains three data frames that act as the starting data for all taxonomy linked analyses corresponding to the LogMPie Dataset (in the current study). The logmpie_select_age_final_genus contains the genus-level profiles; logmpie_select_age_final_species contains the species-level profiles and; logmpie_select_age_final_metadata contains all the relevant metadata.

Each code is annotated at the beginning. For each data repository, the running workspace (workspace that is loaded at the beginning of the execution of each sub-pipeline and updated and saved once the sub-pipeline completes execution) and the order of execution of the individual codes along with the outputs of each sub-pipeline is listed below:

1. American Gut: The running workspace is ag_Analysis_2021_Revision.RData. The order of execution of the sub-pipelines along with intermediate outputs are: run_ag_analysis_stage1a.R > run_ag_analysis_stage1b.R (output produced: ag_stage1_results.RData) > run_ag_analysis_stage2a_1.R (output produced: ag_stage2a_results.RData) > run_ag_analysis_stage2b_1.R (output produced: ag_stage2b_results.RData) > run_ag_analysis_stage2c_1_new.R (output produced: ag_stage_2c_results.RData)

2. He et al: The running workspace is he_Analysis_2021_Revision.RData. The order of execution of the sub-pipelines along with intermediate outputs are: run_he_analysis_stage1a.R > run_he_analysis_stage1b.R (output produced: he_stage1_results.RData) > run_he_analysis_stage2a_1.R (output produced: he_stage2a_results.RData) > run_he_analysis_stage2b_1.R (output produced: he_stage2b_results.RData) > run_he_analysis_stage2c_1_new.R (output produced: he_stage_2c_results.RData)

3. Irish Shotgun Cohorts (ISC): The running workspace is isc_analysis_2021_Revision.RData. The order of execution of the sub-pipelines along with intermediate outputs are: run_isc_analysis_stage1a.R > run_isc_analysis_pathway_stage0.R (output produced: isc_stage1_pathway_results.RData) > run_isc_analysis_stage1b.R (output produced: isc_stage1_results.RData) > run_isc_analysis_stage2a_1.R (output produced: isc_stage2a_results.RData) > run_isc_analysis_stage2b_1.R (output produced: isc_stage2b_results.RData) > run_isc_analysis_stage2c_1_new.R (output produced: isc_stage_2c_results.RData)

4. LogMPie: The running workspace is logmpie_analysis_2021_Revision.RData. The order of execution of the sub-pipelines along with intermediate outputs are: run_logmpie_analysis_stage0.R > run_logmpie_analysis_stage1a.R > run_logmpie_analysis_stage1b.R (output produced: isc_stage1_results.RData) > run_logmpie_analysis_stage2a_1.R (output produced: logmpie_stage2a_results.RData) > run_logmpie_analysis_stage2b_1.R (output produced: logmpie_stage2b_results.RData)

5. 

Please contact tarini.ghosh@ucc.ie for details.
