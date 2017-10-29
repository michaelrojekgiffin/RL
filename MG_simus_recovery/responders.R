###### Looks at the responders data #####
ug <- read.table('~/Dropbox/RL/ultimatum_responders.csv', header=TRUE, sep='')
write.table(ug, '~/Dropbox/RL/MG_simus_recovery/ultimatum_responders.txt', sep='\t', row.names=FALSE)

length(unique(ug$subject))
length(unique(ug$group))
length(unique(ug$email))
length(unique(ug$endowment))

# want to see if there are different resistance points depending on the starting endowment.
# In order to do this, I need to see if there's a difference in the acceptance rate of unfair offers
# depending on the starting endowment. So I need to specify which offers are unfair

######histograms of the acceptace rates#######

end_0 = rep(-999, length(unique(ug$offer)))
end_5 = rep(-999, length(unique(ug$offer)))
end_10 = rep(-999, length(unique(ug$offer)))
end_15 = rep(-999, length(unique(ug$offer)))
end_20 = rep(-999, length(unique(ug$offer)))
for (ii in 0:length(off_hist)-1) {
  end_0[ii+1] = eval( parse( text = sprintf( 'nrow(ug[ which (ug$endowment == 0 & ug$offer == %d & ug$response == "accept"),]) / nrow(ug[ which (ug$endowment == 0 & ug$offer == %d),])', ii, ii)))
  end_5[ii+1] = eval( parse( text = sprintf( 'nrow(ug[ which (ug$endowment == 5 & ug$offer == %d & ug$response == "accept"),]) / nrow(ug[ which (ug$endowment == 5 & ug$offer == %d),])', ii, ii)))
  end_10[ii+1] = eval( parse( text = sprintf( 'nrow(ug[ which (ug$endowment == 10 & ug$offer == %d & ug$response == "accept"),]) / nrow(ug[ which (ug$endowment == 10 & ug$offer == %d),])', ii, ii)))
  end_15[ii+1] = eval( parse( text = sprintf( 'nrow(ug[ which (ug$endowment == 15 & ug$offer == %d & ug$response == "accept"),]) / nrow(ug[ which (ug$endowment == 15 & ug$offer == %d),])', ii, ii)))
  end_20[ii+1] = eval( parse( text = sprintf( 'nrow(ug[ which (ug$endowment == 20 & ug$offer == %d & ug$response == "accept"),]) / nrow(ug[ which (ug$endowment == 20 & ug$offer == %d),])', ii, ii)))
}
barplot(end_0, names=unique(ug$offer), ylim = c(0, 1), xlab = "offers", ylab = "acceptance frequency", main = "endowment = 0")
barplot(end_5, names=unique(ug$offer), ylim = c(0, 1), xlab = "offers", ylab = "acceptance frequency", main = "endowment = 5")
barplot(end_10, names=unique(ug$offer), ylim = c(0, 1), xlab = "offers", ylab = "acceptance frequency", main = "endowment = 10")
barplot(end_15, names=unique(ug$offer), ylim = c(0, 1), xlab = "offers", ylab = "acceptance frequency", main = "endowment = 15")
barplot(end_20, names=unique(ug$offer), ylim = c(0, 1), xlab = "offers", ylab = "acceptance frequency", main = "endowment = 20")

