# draw networks to figures in paper

setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")


library(igraph)
#e <-  c(5,1, 5,4, 5,3, 5,2, 5,7, 6,3, 2,6, 3,7, 4,7, 1,7, 
#        1,5, 4,5, 3,5, 2,5, 7,5, 6,2, 7,3, 7,4, 7,1, 3,6)

e <-  c(1,5, 6,3, 2,6, 3,7, 3,5, 4,7, 1,7, 2,5, 7,5, 5,4, 
        5,1, 6,2, 7,3, 5,3, 7,4, 7,1, 5,2, 3,6, 5,7, 4,5)

g <- graph(e, n = 7, directed = TRUE)

col = c(rep("red", 5), rep("blue", 15))
type = c(rep(FALSE, 4), rep(TRUE, 3))

curve.reciprocal.edges <- function(g, curve=.3){
  # Return a graph where the edge-attribute $curved is reset to highlight reciprocal edges
  el <- t(apply(get.edgelist(g), 1, sort))
  #V(g)$color = "lightblue"
  V(g)$color = "gray18"
#  V(g)$color[5] = "red"
  E(g)$color = col
  V(g)$shape <- ifelse(type, "circle", "square")
  E(g)$curved <- 0
  E(g)[duplicated(el) | duplicated(el,fromLast =TRUE)]$curved <- curve
  (g)
}

tiff("distributed_net.tiff", width = 10, height = 10, units = "cm", res = 300)
par(bg=NA)
plot(curve.reciprocal.edges(g), layout=layout.circle, vertex.size=15,
     edge.arrow.size=0.4, vertex.label = NA)
dev.off()


#####

library(bipartite)
library(igraph)
source("~/Dropbox/Master/Code/coevo_mut_antag/R/functions/SquareMatrix.R")


setwd("~/Dropbox/Master/Code/coevo_mut_antag/data/")

temp = list.files(pattern = "*.txt")
redes = lapply(temp, read.table)
names(redes)  = gsub(".txt", replacement = "", temp)

pp <- as.matrix(redes[[22]])
pp = SquareMatrix(pp)
ppIg<-graph_from_adjacency_matrix(pp, mode = "undirected")

vertex_attr(ppIg)
vertex_attr(ppIg)$type = c(rep(TRUE, 10), rep(FALSE, 8))

vertex_attr(ppIg)$color<-rep("#FF9C0B", length(V(ppIg))) #coloring all the nodes orange
#Selecting everyone that is "TRUE" in $type to change it to green (take a look at the help of the grep function in case you don't know it).
vertex_attr(ppIg)$color[grep(pattern = "FALSE", vertex_attr(ppIg)$type)]<-"#1C912C"
vertex_attr(ppIg)$color # checking the color vector

wc1 <- which(V(ppIg)$type=="TRUE")
wc2 <- which(V(ppIg)$type=="FALSE")
V(ppIg)$shape[wc1] <- "square"
V(ppIg)$shape[wc2] <- "circle"

tiff("modular_net.tiff", width = 10, height = 10, units = "cm", res = 300)
par(bg=NA)
plot(ppIg, vertex.color="black",vertex.label=NA,
     vertex.size=8, edge.width=1.5, vertex.shapes = V(ppIg)$shape,
     layout=layout_as_bipartite)
dev.off()
