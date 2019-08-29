# Figure 2
# Linear regression between p and trait variance

# Read data
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")
library(ggplot2)
library(cowplot)

dados = read.table("planilha_final.txt", header = TRUE, dec = ",")

plot = ggplot(data = dados) +
  geom_point(aes(x = antprob, y = variance, color = n)) +
#  scale_color_gradient(low="gray40", high="black") +
  xlab("p") +
  ylab("Species trait variance") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 18), 
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11))

plot
