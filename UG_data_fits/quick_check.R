# for some quick tests
UG <- read.table('/Users/michaelgiffin/Dropbox/RL/data/pilot/3cond_blocked/interim/my_ultimatum_3cond_block.txt', header= TRUE, sep='\t')
nrow(UG[which (UG$participant.code == 'kwi91ww8' & UG$group.current_opponent == 'circle'),])
nrow(UG[which (UG$participant.code == 'kwi91ww8' & UG$group.current_opponent == 'square'),])
nrow(UG[which (UG$participant.code == 'kwi91ww8' & UG$group.current_opponent == 'triangle'),])


t.test(UG[which (UG$group.soc_or_no == 'Non-social'),]$group.amount_offered, UG[which (UG$group.soc_or_no == 'Social'),]$group.amount_offered)



lr <- read.table('/Users/michaelgiffin/Desktop/learningrates.txt', header=FALSE, sep='\t')
lr <- stack(lr)

lr[,3] <- rep(-999, nrow(lr))
colnames(lr) <- c("lr", "cond", "sub_name")

counter = 0
condition = "social"
for (ii in 1:nrow(lr)) {
  counter = counter + 1
  # lr[ii, 2] <- sprintf('%s', condition)
  lr[ii, 3] <- counter
  if (counter == nrow(lr)/2) {
    counter = 0
    condition = "nonsocial"
  }
}


ezANOVA(data = lr, dv = lr, wid = sub_name, within = cond)
