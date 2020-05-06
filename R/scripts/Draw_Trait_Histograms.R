# Simple code to draw the histograms
# from the Figure 2 of the paper.
# You need to have an z_mat matriz of species traits

par(bg = NA)
hist(z_mat, # histogram
      col = "peachpuff", # column color
      border = "black",
      prob = TRUE, # show densities instead of frequencies
      xlab = "Trait values",
      main = "", ylim = c(0,0.15))
 lines(density(z_mat), # density plot
       lwd = 2, # thickness of line
       col = "black")
