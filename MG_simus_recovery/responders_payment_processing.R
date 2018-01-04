# routine to pay the responders
ug <- read.table('~/Dropbox/RL/ultimatum_responders.csv', header=TRUE, sep='')

# 20 MU = €2.50
mu <- 2.5/20
all_email <- unique(ug$email)

ug_pay <- data.frame(email =  unique(ug$email), payment = rep(-999, length(all_email)))

for (ii in 1:length(all_email)) {
  samp_trial <- sample(c(1:nrow(ug[which (ug$email == all_email[ii]),])), 1)
  ug_pay[ii,]$email   <- as.character(all_email[ii])
  ug_pay[ii,]$payment <- ug[which (ug$email == all_email[ii]),][samp_trial,]$payoff * mu
}

write.table(ug_pay, '~/Dropbox/RL/MG_simus_recovery/responders_payment.txt', sep='\t', row.names=FALSE)
