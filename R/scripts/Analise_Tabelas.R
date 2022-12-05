setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

load("antprob_var.RData")
low = p_data[which(p_data$antprob == 0.01), ]
med = p_data[which(p_data$antprob == 0.5), ]
high = p_data[which(p_data$antprob == 0.99), ]

mean(low$opt_clusters)
sd(low$opt_clusters)

mean(med$opt_clusters)
sd(med$opt_clusters)

mean(high$opt_clusters)
sd(high$opt_clusters)

###
pol = p_data[which(p_data$type == "Pollination"),]
seed = p_data[which(p_data$type == "Seed dispersal"),]
ant = p_data[which(p_data$type == "Ant-Plant"),]

mean(ant$mpd)
sd(ant$mpd)

summary(lm(pol$mpd ~ pol$antprob))

summary(lm(pol$mpd ~ pol$antprob))

summary(lm(ant$mpd ~ ant$antprob))


####
load("minor_data_structure.RData")
pol = final_fl[which(final_fl$type == "Pollination"),]
seed = final_fl[which(final_fl$type == "Seed dispersal"),]
ant = final_fl[which(final_fl$type == "Ant-Plant"),]

mean(ant$dnest_control)
sd(ant$dnest_control)
mean(ant$dnest_coevo)
sd(ant$dnest_coevo)

mean(ant$dmod_control)
sd(ant$dmod_control)
mean(ant$dmod_coevo)
sd(ant$dmod_coevo)

####

load("central_results.RData")
pol = central_results[which(central_results$type == "Pollination"),]
seed = central_results[which(central_results$type == "Seed dispersal"),]
ant = central_results[which(central_results$type == "Ant-Plant"),]

mean(ant$mpd)
sd(ant$mpd)

central = ant[which(ant$c_ch == "Central"), ]
distributed = ant[which(ant$c_ch == "Distributed"), ]

t.test(central$mpd, distributed$mpd)

