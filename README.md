# Insert cheaters interaction outcomes in mutualistic networks changes the coevolutionary process in unexpected ways

Links of the code that generate my figures. IN PROGRESS:

- [Figure 2](https://github.com/lucascamacho/coevo_mut_antag/blob/master/R/scripts/Empirical_Antprob(0-10)_Figure2.R)
- [Figure 3](https://github.com/lucascamacho/coevo_mut_antag/blob/master/R/scripts/Empirical_CentralAntprob_Figure3.R)
- [Figure 4](https://github.com/lucascamacho/coevo_mut_antag/blob/master/R/scripts/Empirical_ForbLinks_Antprob_Figure4.R)

## Technical specifications

I run my simulations in two different computers:
- Desktop MacPro (Mid 2010) 2x2.4 GHz Quad-Core Intel Xeon, 6GB DDR3, macOS High Sierra 10.13.6
- Notebook MacBook Pro (13-inch, Mid 2012) 2.9 GHz Intel Core i7, 8GB DDR3, macOS High Sierra 10.13.6

The softwares used was:
- Dropbox 70.4.93
- GitHub Desktop 1.6.5
- MODULAR 0.1(1)
- R 3.5.3
- RStudio 1.1.453
- TexMaker 5.0.2
- Zotero 5.0.61


### Repository organization

This repository is organized in folders which contain different aspects of my master's project:

- R: separate in functions and scripts. You will find all the .R files to run my simulations
- Data: the networks of interactions that i use and the .RData files with the results of all my simulations
- Output: figures, pdf's, graphs, etc. All the visual results of my project are here
- Manuscript: the LaTeX files that compose my master's dissertation.
- Tutorials: Rmd and PDF files to help the understanding of the model and the simulation process

### Tutorial and guides

To improve the understanding of the coevolutionary model and the simulations, i prepare some small guides and tutorials to help people know how my model and code works. You 
will find a guide to my model and R code in the main page of this repository and some tutorials together with some scripts in R/scripts folder. I hope this guide helps people 
that are interested in my work and my results. Also, and equally important, i hope this guides and tutorials will increase the understanding about what i've been doing in my master's.
Please, let me know if these tutorials and guides are usefull or has some erros or inconsistencies sending me an email (lucas.camacho@usp.br). Check below the lists of 
guides and tutorials available:

- [Guide_CoevoMutAntag](https://github.com/lucascamacho/coevo_mut_antag/blob/master/tutorials/Guide_CoevoMutAntag.pdf): Guide showing my project basic ideas, my model and how do 
i implement it on R. These guides ends in a single plot showing the average species traits changing in time.

- [Running_CoevoMutAntag](https://github.com/lucascamacho/coevo_mut_antag/blob/master/tutorials/Running_CoevoMutAntag.pdf): Tutorial to explain what i've being doing in the scripts where i explore how the presence of cheaters influences the coevolutionary process 
in mutualistic networks. You will that i try to explain step-by step the basic ideas of my scripts and how do i visualize the results of my simulations.


### Running my simulations

My scripts are directed to my own working directory. Please, remember to use your own working directory when run the simulations. The scripts already have an automate way of installing the 
required packages listed below:

List of R packages used:
- [ggplot2](https://ggplot2.tidyverse.org)
- [cowplot](https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html)
- [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
- [bipartite](https://cran.r-project.org/web/packages/bipartite/index.html)
- [plyr](https://www.rdocumentation.org/packages/plyr/versions/1.8.4)
- [viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)



