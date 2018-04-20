This folder has everything related to the predator prey game and reinforcement learning

the empirical_priors are scripts actually get the priors from the data, but I still need to use the Laplace estimation technique in order to get these. They use the Prior_Estimation funcion as well to actually do the estimating, the empirical_priors scripts do all the looping over subjects and that heavy lifting

fit_data and get_data change the R formatted data files into something more matlabby

model_vs_subject both fit the data to subjects and simulate game play using the same parameters estimated from the subject and playing agains the same sequence of feedback the subject received. 

The ModelEstimation_2params scripts actually estimate the choice prediction error models when we already set the priors (so the priors are known and don't have to be estimated, hence the 2 parameters)

recovery_mg_predprey tests to see if the models are recoverable, they contain a lot of the same information that is in the model_vs_subject script, except these actually simulate data with randomly chosen parameters and then attempts to estimate those paramters and creates plots showing the correlation between the true parameters and the estimated parameters.

test_recovery_models_MG_2017_09_21 is very important, this is the script that simulates PPG behavior of 4 different models: choice prediction error, choice prediction error with risk parameter, reward prediction error, reward prediciton error with risk parameter, and sees if the models are recoverable and identifiable from one another. This script uses the function laplace_priors_learning2_MG, learning_model_estim, and learning_models_timeseris


the directory scripts_big_hor is where I started working after the UG so the scripts are much more refined. In that directory, I fit the new models (with 4 and 5 free parameters in which we estimated the priors from the subject data), and in these models the recovered parameters made no sense, and then after running simulations I found that the parameter recovery is inconsistent, and model identifiability is also pretty bad. Therefore I'm going back to using the priors as fixed effects estimated from the first decision made by each participant, am running more similations to make sure that these look good, and will then fit that to my dataset. I think that the simulations look good, but I'm wondering if, when fit to the real data, it fits well. I'm going to start with the big data set, the big-hor dataset, and see if there's a clear winner of the models.
