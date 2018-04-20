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
uglm = lmer(offer ~ ( social * opponent * opp_trial) + (1  | sub_name), data=UG[which (UG$explicit == 1 & UG$opp_trial < 50),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)
# the effects reverse and become marginal, what is going on?
# 
# 
# 
# 
# here is the version with the restricted hisotry - so only condsidering 
# the history that ALL subject's shared (which I don't like as a criterion anyway)
uglm = lmer(offer ~ ( social * opponent * opp_trial) + (1  | sub_name), data=UG[which (UG$explicit == 1 & UG$opp_trial < 25),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)




uglm = lmer(payoff ~ ( social * opponent + explicit) + (1  | sub_name), data=UG[which (UG$opp_trial < 25),])
vif.lme(uglm)
summary(uglm)




## do subjects earn more in the nonsocial?
uglm = lmer(payoff ~ ( social * opponent ) + (1  | sub_name), data=UG[which (UG$explicit == 1 & UG$opp_trial < 25),])
vif.lme(uglm)
Anova(uglm, type = 3)
summary(uglm)


t.test(UG[which (UG$explicit == 1 & UG$opp_trial < 25 & UG$social== 0),]$payoff, UG[which (UG$explicit == 1 & UG$opp_trial < 25 & UG$social== 1),]$payoff, paired = TRUE)

t.test(UG[which ( UG$opp_trial < 25 & UG$social== 0),]$payoff, UG[which (UG$opp_trial < 25 & UG$social== 1),]$payoff, paired = TRUE)

t.test(UG[which ( UG$social== 0),]$payoff, UG[which (UG$social== 1),]$payoff, paired = TRUE)


###### plotting #####
uglm = lmer(scale(payoff) ~ ( social * opponent + explicit) + (1  | sub_name), data=UG[which (UG$opp_trial < 25),])
vif.lme(uglm)
summary(uglm)

z1 <- as.integer(c(0, 1))
z2 <- as.integer(c( 1, 2, 3))
mynewdf <- expand.grid(social=z1,opponent=z2)
mynewdf[,3] <- rep(0, nrow(mynewdf))

colnames(mynewdf) = c("social", "opponent", "explicit")
mynewdf[,1] <- as.factor(mynewdf[,1])
mynewdf[,2] <- as.factor(mynewdf[,2])
mynewdf[,3] <- as.factor(mynewdf[,3])
## re.form=NA apparently means that the random effects are all set to 0, according to this: https://stats.stackexchange.com/questions/29690/getting-fixed-effect-only-predictions-from-mixed-model-on-new-data-in-r
the_predict <- (predict(uglm, mynewdf, re.form=NA))
mydata <- transform(mynewdf, res = the_predict)

#575 X 520
p <- ggplot(data = mydata, aes(y = res, x = opponent, color=factor(social))) + stat_smooth(method=lm)
p + scale_colour_discrete(name="social") + scale_x_continuous(breaks=seq(0, 1)) + theme_bw() + labs(title = "thesis data T X C & expected reward", y = "expected reward") 






tclm <- lmer(exp_rew ~ ( ztest1 * zcort1 + role + scale(bmi)) + (1 | sub_name), data=old_hor)
z1 <- z2 <- seq(-1,1)
mynewdf <- expand.grid(ztest=z1,zcort=z2)
mynewdf[,3] <- rep("predator", nrow(mynewdf))
mynewdf[,4] <- (rep(0, nrow(mynewdf)))
colnames(mynewdf) = c("ztest1", "zcort1", "role", "bmi")
## re.form=NA apparently means that the random effects are all set to 0, according to this: https://stats.stackexchange.com/questions/29690/getting-fixed-effect-only-predictions-from-mixed-model-on-new-data-in-r
the_predict <- (predict(tclm, mynewdf, re.form=NA))
mydata <- transform(mynewdf, res = the_predict)

#575 X 520
p <- ggplot(data = mydata, aes(y = res, x = zcort1, color=factor(ztest1))) + stat_smooth(method=lm)
p + scale_colour_discrete(name="ztest1") + scale_x_continuous(breaks=seq(-1,1)) + theme_bw() + labs(title = "thesis data T X C & expected reward", y = "expected reward") 
