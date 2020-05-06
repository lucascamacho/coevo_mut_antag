# simple code to generate coloured matrices like the one
# we used in our Figure 4
#
# You will need a binary matrix, here called last_m.
#
# load packages and functions

library(ggplot2)
library(dplyr)

dat_long = melt(last_m) # melt matrix
dat_long$value = factor(dat_long$value) # discrete vs continuous
gg = ggplot(dat_long)
gg = gg + geom_tile(aes(x=Var1, y=Var2, fill=value), 
                     show.legend = FALSE) # fill + legend, gray border
gg = gg + scale_fill_manual(values=c("white", "#7570B3")) # custom fill colors
gg = gg + coord_equal() # squares
gg = gg + labs(x = NULL, y = NULL) # no labels
gg = gg + theme_bw() # remove some chart junk
gg = gg + theme(panel.grid = element_blank())
gg = gg + theme(panel.border = element_blank())
gg

#save the matrix
ggsave(gg, filename = "Seed_High_Antprob.pdf", dpi = 600,
       width = 10, height = 10, units = "cm",  bg = "transparent")