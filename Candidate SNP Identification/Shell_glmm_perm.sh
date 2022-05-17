

#glmm perms were run as array jobs. How this will done may vary from system to system.

#high cont jobs
SGE_Batch -c 'Rscript High_Control_GLMM_perms.R $SGE_TASK_ID' -t 1-200 -b 20 -r HC_GLMM -q beagle -m 1.5G - f 1.5G

SGE_Batch -c 'Rscript High_Moderate_GLMM_perms.R $SGE_TASK_ID' -t 1-200 -b 20 -r HM_GLMM -q beagle -m 1.5G - f 1.5G

SGE_Batch -c 'Rscript Moderate_Control_GLMM_perms.R $SGE_TASK_ID' -t 1-200 -b 20 -r MC_GLMM -q beagle -m 1.5G - f 1.5G
