# APPLY FUNCTIONS NEW
traits = t(apply(z_mat, 1, TraitDegreeBalanced))
colnames(traits) = c("AM", "MM")
var_traits = t(apply(z_mat, 1, VarTraitDegreeBalanced))
colnames(var_traits) = c("AM", "MM")

#OTIMIZAR ------
traitsAM = traits[,1]
type = rep("AM", times = length(traitsAM))
traitsAM = data.frame(traitsAM, type)

traitsMM = traits[,2]
type = rep("MM", times = length(traitsMM))
traitsMM = data.frame(traitsMM, type)

vartraitsAM = var_traits[,1]
type = rep("AM", times = length(vartraitsAM))
vartraitsAM = data.frame(vartraitsAM, type)

vartraitsMM = var_traits[,2]
type = rep("MM", times = length(vartraitsMM))
vartraitsMM = data.frame(vartraitsMM, type)

traitsAM = as.matrix(traitsAM)
traitsMM = as.matrix(traitsMM)

colnames(traitsAM) = NULL
colnames(traitsMM) = NULL

new_traits = rbind(traitsAM, traitsMM)

vartraitsAM = as.matrix(vartraitsAM)
vartraitsMM = as.matrix(vartraitsMM)

colnames(vartraitsAM) = NULL
colnames(vartraitsMM) = NULL

new_vartraits = rbind(vartraitsAM, vartraitsMM)

data = cbind(new_traits, new_vartraits)
data = data[,-2]

data = data.frame(data)

data[,1] = as.numeric(as.character(data[,1]))
data[,2] = as.numeric(as.character(data[,2]))

# TRAIT K DIFF ANTPROB 0-10 PROBLEM SOLVER




#


# second data frame for last line of z_mat for each simulation
last_traits = matrix()
colnames(last_traits) = c("mean", "var", "type")

# loop to coevolution simulation and get the last line of z_mat for each simulation
for(a in 1:length(antprob_vec)){ #
  antprob = antprob_vec[a] # define a value fo antprob
  print(antprob) # print this value to follow the simulation process
  
  # OTIMIZAR -----------  
  col_variables = matrix()
  colnames(col_variables) = NULL
  
  # for each value of antprob, simulate 10 times
  for(y in 1:100){ 
    source("Trait_Diff_K_AM_MM.R")
    g.data = data[!duplicated(data$X3, fromLast=TRUE),]
    g.data[1,2] = g.data[1,2] / c[[2]]
    g.data[2,2] = g.data[2,2] / c[[3]]
    
    col_variables = bind_rows(as.data.frame(col_variables), g.data)
    
  }
  
  col_variables = col_variables[-1,-1]
  colnames(col_variables) = c("mean", "var", "type")
  
  final_m_AM = mean(col_variables$mean[which(col_variables$type == "AM")])
  final_var_AM = mean(col_variables$var[which(col_variables$type == "AM")])
  final_m_MM = mean(col_variables$mean[which(col_variables$type == "MM")])
  final_var_MM = mean(col_variables$var[which(col_variables$type == "MM")])
  
  quase = matrix(c(final_m_AM, final_m_MM, final_var_AM, final_var_MM), nrow = 2, ncol = 2)
  
  last_traits = bind_rows(as.data.frame(last_traits), as.data.frame(quase))
  
}

dados = read.table("dados.txt")
last_traits = as.data.frame(dados)
#last_traits = last_traits[1,] = NA
last_traits[,3] = rep(c("MM", "AM"), times = 100)
last_traits[,4] = rep(1:100, each=2) / 100


## ggplot
plot = ggplot(data = last_traits) + geom_point(aes(x = last_traits$V4, y = last_traits$V1, 
                                                   color = last_traits$V3), size = 5, alpha = 0.6) +
  theme_classic(base_size = 14) +
  labs(x = "Probabilidade P", 
       y = "Z-Score dos traços médios das espécies balanceados pelo grau",
       title = "Z-Score dos traços médios de Antagonismos e Mutualismos ao mudarmos P",
       color = "Interações")

plot
ggsave(plot, file = "medias_am_mm.png", dpi = 600, width = 12, height = 8, units = "in")

plot2 = ggplot(data = last_traits) + geom_point(aes(x = last_traits$V4, y = last_traits$V2, 
                                                    color = last_traits$V3), size = 5, alpha = 0.6) +
  theme_classic(base_size = 14) +
  labs(x = "Probabilidade P",
       y = "Módulo da soma das diferenças dos traços entre as espécies",
       title = "Módulo da soma das diferenças dos traços entre Antagonismos e Mutualismos ao mudarmos P",
       color = "Interações")

plot2
ggsave(plot2, file = "vars_am_mm.png", dpi = 600, width = 12, height = 8, units = "in")