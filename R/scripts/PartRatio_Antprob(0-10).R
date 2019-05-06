setwd("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/")

res_01 = vector()
res_05 = vector()
res_09 = vector()


for(i in 1:101){
  antprob = 0.1
  print(i)
  source(("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/Basic_Traits.R"))
  res_01 = append(res_01, partratio)
}

for(i in 1:101){
  antprob = 0.5
  print(i)
  source(("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/Basic_Traits.R"))
  res_05 = append(res_05, partratio)
}

for(i in 1:101){
  antprob = 0.9
  print(i)
  source(("~/Dropbox/Master/Code/coevo_mut_antag/R/scripts/Basic_Traits.R"))
  res_09 = append(res_09, partratio)
}

mean(res_09)
var(res_09)

x = seq(0, 10, 0.1)
plot(x, res_01, col = "blue")
plot(x, res_05, col = "yellow")
plot(x, res_09, col = "red")
