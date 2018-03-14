this directory contains all the simulations and behavioral data for the Ultimatum Reinforcement learning study conducted by Giffin, Lebreton, Gross, and De Dreu.

Directory tour:
arxiv_simus: archive of old simulations
data: data from the behavioral study (henceforth study 1) in both raw and processed forms
data_analysis_ML: Mael's scripts for analyzing behavioral data for both study 1 and the fMRI follow-up (henceforth study 2).
fMRI_stimulus: all the scripts used to present the 

data_analysis_ML: where Mael has done the modeling analysis, contains the following important directories:
	BEHAVcohortwithout_priors: this directory has the first tier of the modeling analysis (comparing our 4 main models), BEHAVcohordwithout_piors/run_model_learning_all.m creates a .mat file in BEHAVcohordwithout_piors_integr called Learning_All_2017_12_08.mat, which is what's read in with analysis_learning_models_all_cond_ModelSpace1.m. 
	BEHAVcohordwithout_piors_integr: contains the second tier of the model based analysis in which parameters are allowed to vary randomly across conditions, contains the aforementioned .mat file along with Learning_All_2018_03_05.mat, which contains the results of the model-fitting of the second tier of the model, and is created by BEHAVcohordwithout_piors_integr/run_model_learning_all.m
	

The above two distinctions are a bit subtle but very important to remember.

simus_recovery_learning_without_piors is the directory Mael made after we realized that the prior estimation task didn't work and we had to estimate the priors from the data, (I still need to run some simulations to see if the PPG works with this framework)