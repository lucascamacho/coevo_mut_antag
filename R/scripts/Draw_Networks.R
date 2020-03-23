# This script generates the network figures from the 2 to 4 figures
# presented in the paper.
#
# All the networks are theoretical, only representing the shape of
# real empirical networks.

# define work directory
setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

#load packages and functions
library(igraph)

# centrality network
#e <-  c(5,1, 5,4, 5,3, 5,2, 5,7, 6,3, 2,6, 3,7, 4,7, 1,7, 
#        1,5, 4,5, 3,5, 2,5, 7,5, 6,2, 7,3, 7,4, 7,1, 3,6)

# random network
e <-  c(1,5, 6,3, 2,6, 3,7, 3,5, 4,7, 1,7, 2,5, 7,5, 5,4, 
        5,1, 6,2, 7,3, 5,3, 7,4, 7,1, 5,2, 3,6, 5,7, 4,5)

# make graph using e
g <- graph(e, n = 7, directed = TRUE)

# define colors and types (bipartite)
col = c(rep("red", 5), rep("blue", 15))
type = c(rep(FALSE, 4), rep(TRUE, 3))

# use igraph to plot networks
curve.reciprocal.edges <- function(g, curve=.3){
  # Return a graph where the edge-attribute $curved is reset to highlight reciprocal edges
  el <- t(apply(get.edgelist(g), 1, sort))
  #V(g)$color = "lightblue"
  V(g)$color = "gray18"
  #V(g)$color[5] = "red" # color of single exploiter
  E(g)$color = col # node colors
  V(g)$shape <- ifelse(type, "circle", "square") # bipartite
  E(g)$curved <- 0
  E(g)[duplicated(el) | duplicated(el,fromLast =TRUE)]$curved <- curve
  (g)
}

# save the current plotted network
pdf("1_not_central_net.pdf", width = 10, height = 10, bg = "transparent")
par(bg = NA)
plot(curve.reciprocal.edges(g), layout=layout.circle, vertex.size=15,
     edge.arrow.size=0.4, vertex.label = NA)
dev.off()