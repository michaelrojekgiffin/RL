# for some model free exploring
# this is the second and more correct version of this script
# with better models. 
# It's importat to not here that the first part of this script is data shaping,
# and that I'm planning to add something akin to this in the matlab script 
# that originally shapes this stuff, which will make the first portion of this script
# obsolete. However, the lme4 models are the meat of this script, and that's
# what I use for the results shown in the summary powerpoint.


UG = read.table('/Users/michaelgiffin/Dropbox/RL/data/processed_R/UG_all.txt', header=TRUE, sep = '\t')
UG = UG[which (UG$age < 35),]

library(Hmisc)
library(ez)
library(car)
library("lme4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library(lmerTest)
library(ggplot2)

UG$social = as.factor(UG$social)
UG$opponent = as.factor(UG$opponent)
UG$explicit = as.factor(UG$explicit)

ncol(UG)

#### adding first and second half ####
UG[,10] <- rep(-999, nrow(UG))
colnames(UG)[10] <- 'half'
trial_half = length(unique(UG$trial))/2
UG[which (UG$trial <= trial_half),]$half <- 1
UG[which (UG$trial > trial_half),]$half <- 2

#### adding timed quartiles ####
UG[,11] <- rep(-999, nrow(UG))
# time quartile, is what this stands for
colnames(UG)[11] <- 'tq'
trial_quarter = length(unique(UG$trial))/4
UG[which (UG$trial <= trial_quarter),]$tq <- 1
UG[which (UG$trial > trial_quarter & UG$trial <= trial_quarter*2),]$tq <- 2
UG[which (UG$trial > trial_quarter*2 & UG$trial <= trial_quarter*3),]$tq <- 3
UG[which (UG$trial > trial_quarter*3),]$tq <- 4

UG$tq <- as.factor(UG$tq)
UG$trial = as.numeric(UG$trial)
UG$half = as.factor(UG$half)

###### Adding trial that corresponds to subject's history with specific opponent #####
ncol(UG)
UG[,12] <- rep(-999, nrow(UG))
# this is an indicator of number of trials for which the participants have been
# exposed to the same opponent
colnames(UG)[12] <- 'opp_trial'
all_subs <- unique(UG$sub_name)
uopps <- unique(UG)

for (ii in 1:length(all_subs)) {
  cur_sub = all_subs[ii]
  cur_rows = as.numeric(row.names(UG[which (UG$sub_name == cur_sub),]))
  # temp_data = UG[which (UG$sub_name == cur_sub),]
  for (ss in 1:length(cur_rows)) {
    if (UG[cur_rows[ss],]$trial == 1) {
      opp1count = 0
      opp2count = 0
      opp3count = 0 
    }
    if (UG[cur_rows[ss],]$opponent == 1) {
      opp1count = opp1count + 1
      UG[cur_rows[ss],]$opp_trial = opp1count
    }
    if (UG[cur_rows[ss],]$opponent == 2) {
      opp2count = opp2count + 1
      UG[cur_rows[ss],]$opp_trial = opp2count
    }
    if (UG[cur_rows[ss],]$opponent == 3) {
      opp3count = opp3count + 1
      UG[cur_rows[ss],]$opp_trial = opp3count
    }
  }
}


### propper models##
uglm = lmer(offer ~ ( social * opponent  * opp_trial * explicit) + (1  | sub_name), data=UG)
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)
 # significcant 4-way interaction that justifies breaking up explicit and implict conditions


uglm = lmer(offer ~ ( social * opponent  * opp_trial) + (1  | sub_name), data=UG[which (UG$explicit == 1),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)
output_table <- data.frame(summary(uglm)$coefficients)
write.table(output_table, file = "/Users/michaelgiffin/Dropbox/RL/UG_data_fits/figures/exp_regression_output.txt", sep='\t')

# in the explicit condition, an individual in the social condition playing against oppponet 3 (accepts anything)
# ofers significantly less over trials than when that person is in the non-social condition


uglm = lmer(offer ~ ( social * opponent * opp_trial ) + (1  | sub_name), data=UG[which (UG$explicit == 0),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)
# this is not the case in the implicit condition


uglm = lmer(offer ~ ( explicit * opponent * opp_trial) + (1  | sub_name), data=UG[which (UG$social == 1),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)
# in the social condition, an individual in the explici condition offers opponent 3 (accepts anything)
# significantly less over trials than if they were in the implicit condition

uglm = lmer(offer ~ ( explicit * opponent * opp_trial) + (1  | sub_name), data=UG[which (UG$social == 0),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)
# this effect is not present in the non-social condition
# 


## What this seems like is that the social effect is in the opposite direction that
## i predicted, but probably for different reasons, I don't think that the fairness
## norm is impeding learning, I think that people are actually learning better in
## the social explicit condition becaus they think that they are playing aginst 
## indivioduals who are using the fairness norm and are more predictable. 
## Perhaps if we hadn't used terminlogy like "lottery" in the instructions we
## would not have this effect, becuase I suspect that people suspect that the 
## computers are simply less predictable than humans or something to that effect.
## If that's the case, then we need to rethink some things, but I need to plot this 
## stuff out first, I think that the tail ends of these opp_trials are quite
## noisy and may have to be taken out. 

# but the effect totally disappears if I do that, here is the nice version with all the data
uglm = lmer(offer ~ ( explicit * opponent * opp_trial) + (1  | sub_name), data=UG[which (UG$social == 1),])
vif.lme(uglm)
summary(uglm)
# highly significantl 3-way interaction

# here is the version with the restricted hisotry - so only condsidering 
# the history that ALL subject's shared (which I don't like as a criterion anyway)
uglm = lmer(offer ~ ( explicit * opponent * opp_trial) + (1  | sub_name), data=UG[which (UG$social == 1 & UG$opp_trial < 40),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)
# effect disappears completely! 




# but the effect totally disappears if I do that, here is the nice version with all the data
uglm = lmer(offer ~ ( social * opponent * opp_trial) + (1  | sub_name), data=UG[which (UG$explicit == 1),])
vif.lme(uglm)
summary(uglm)
# significantl 3-way interaction

# here is the version with the restricted hisotry - so only condsidering 
# the history that ALL subject's shared (which I don't like as a criterion anyway)
uglm = lmer(offer ~ ( social * opponent * opp_trial) + (1  | sub_name), data=UG[which (UG$explicit == 1 & UG$opp_trial < 25),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)
# the effects reverse and become marginal, what is going on?