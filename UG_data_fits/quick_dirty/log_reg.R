# quick and dirty regression
source('/Users/michaelgiffin/Carsten PhD/hormones/scripts/expectation_functions.R')

library(Hmisc)
library(ez)
library(car)
library("lme4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library(lmerTest)

ug <- read.table('/Users/michaelgiffin/Dropbox/RL/UG_data_fits/quick_dirty/easy_sub_params.txt', header = TRUE, sep = '\t')

ugm = data.frame(sub = rep(1:50, 2), social = c(rep(1, 50), rep(0, 50)), 
                 a1 = c(ug$soc_exp_a1, ug$nonsoc_exp_a1), a2 = c(ug$soc_exp_a2, ug$nonsoc_exp_a2),
                 a0 = c(ug$soc_exp_a0, ug$nonsoc_exp_a0), beta = c(ug$soc_exp_beta, ug$nonsoc_exp_beta))

uglm <- glmer(social ~ a1 + a2 + a0 + beta + (1 | sub), data = ugm, family = binomial)
summary(uglm)

t.test(ugm[which (ugm$social == 1),]$a0, ugm[which (ugm$social == 0),]$a0, paired=TRUE)

# if this can be trusted, it means that in the nonsocial, subjects start lower (although not significantly, which worries me),
# have a lower beta (so are more exploratory), and have a lower learning rate - so take longer to learn. This is the opposite of what we 
# hypothesized, I don't really trust these results....