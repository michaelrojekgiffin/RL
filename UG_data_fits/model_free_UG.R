# for some model free exploring

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

UG[,10] <- rep(-999, nrow(UG))
colnames(UG)[10] <- 'half'
trial_half = length(unique(UG$trial))/2
UG[which (UG$trial <= trial_half),]$half <- 1
UG[which (UG$trial > trial_half),]$half <- 2


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

ncol(UG)
UG[,12] <- rep(-999, nrow(UG))
# this is an indicator of number of trials for which the participants have been
# exposed to the same opponent
colnames(UG)[12] <- 'opp_trial'



ezANOVA(data = UG, dv = offer, wid = sub_name, within = .(social, opponent, half), between = explicit)

#### full model with trial as interaction ###
uglm = lmer(offer ~ ( social * opponent  * explicit ) + (tq | sub_name), data=UG)
vif.lme(uglm)
Anova(uglm, type = 3)
# summary(uglm)


### full model with quartiles of time###
uglm = lmer(offer ~ ( social * opponent  * explicit * tq) + (1 | sub_name), data=UG)
vif.lme(uglm)
Anova(uglm, type = 3)
# summary(uglm)
## I am now going to attempt to make this in matlab so that at least all the plots we make are in the same general format


### full model with half of time###
uglm = lmer(offer ~ ( social * opponent  * explicit * half) + (1 | sub_name), data=UG)
vif.lme(uglm)
Anova(uglm, type = 3)
# summary(uglm)



uglm = lmer(offer ~ ( social * opponent  * tq ) + (1 | sub_name), data=UG[which (UG$explicit == 1),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)


uglm = lmer(offer ~ ( social * opponent  * tq ) + (1 | sub_name), data=UG[which (UG$explicit == 0),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)




uglm = lmer(offer ~ ( social * opponent  * half ) + (1 | sub_name), data=UG[which (UG$explicit == 1),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)


uglm = lmer(offer ~ ( social * opponent  * half ) + (1 | sub_name), data=UG[which (UG$explicit == 0),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)


aggregate(UG$offer, by=list(opponent = UG$opponent, social =  UG$social, trial = UG$trial), mean)





uglm = lmer(offer ~ ( social * opponent  * trial) + (1  | sub_name), data=UG[which (UG$explicit == 1),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)


uglm = lmer(offer ~ ( social * opponent * trial ) + (1  | sub_name), data=UG[which (UG$explicit == 0),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)



uglm = lmer(offer ~ ( explicit * opponent * trial) + (1  | sub_name), data=UG[which (UG$social == 1),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)
