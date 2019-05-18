# From mutualism to antagonism: coevolutionary dynamics of cheaters and context dependent interactions

Coevolution is a fundamental process on the characterization of ecological communities, creating and breaking interaction between species. Nowadays, we know that different interaction outcomes, like mutualism or antagonism, can influence the species coevolution in distinct ways. However, the ecological interaction outcomes are variable in space and time. The interactions which the outcomes vary in space and time are called context-dependent and can influence the arrange of antagonism/mutualism interaction in time and consequently, change the coevolutionary dynamics os species. Thus, take into account the context-dependent interactions in coevolutionary models and explore how the interactions outcomes shifts changes the coevolutionary process in a next step on the interface between ecology and evolution. Here, we combine a single trait coevolutive model, empirical networks of interactions and numerical simulations to contribute to this next step.

## Technical specifications

I run my simulations in two different computers:
- Desktop MacPro (Mid 2010) 2x2.4 GHz Quad-Core Intel Xeon, 6GB DDR3, macOS High Sierra 10.13.6
- Notebook MacBook Pro (13-inch, Mid 2012) 2.9 GHz Intel Core i7, 8GB DDR3, macOS High Sierra 10.13.6

The softwares used was:
- R 3.5.3
- RStudio 1.1.453
- Zotero 5.0.61
- TexMaker 5.0.2
- GitHub Desktop 1.6.5
- Dropbox 70.4.93

### Repository organization

This repository is organized in folders which contain different aspects of my master's project:

- R: separate in functions and scripts. You will find all the .R files to run my simulations
- Data: the networks of interactions that i use and the .RData files with the results of all my simulations
- Output: figures, pdf's, graphs, etc. All the visual results of my project are here
- Manuscript: the LaTeX files that compose my master's dissertation.

### Tutorial and guides

To improve the understanding of the coevolutionary model and the simulations, i prepare some small guides and tutorials to help people know how my model and code works. You 
will find a guide to my model and R code in the main page of this repository and some tutorials together with some scripts in R/scripts folder. I hope this guide helps people 
that are interested in my work and my results. Also, and equally important, i hope this guides and tutorials will increase the understanding about what i've been doing in my master's.
Please, let me know if these tutorials and guides are usefull or has some erros or inconsistencies sending me an email (lucas.camacho@usp.br). Check below the lists of 
guides and tutorials available:

- [Guide_coevo_mut_antag](https://github.com/lucascamacho/coevo_mut_antag/blob/master/Guide_coevo_mut_antag.pdf): Guide showing my project basic ideas, my model and how do 
i implement it on R. These guides ends in a single plot showing the average species traits changing in time.

- [Tutorial_MutAntNet](): IN PROGRESS. Tutorial to explain what i've being doing in the scripts where i explore how the presence of antagonism (like cheaters) influences the coevolutionary process 
in mutualistic networks. You will that i try to explain step-by step the basic ideas of my scripts and how do i visualize the results of my simulations.

- [Tutorial_ContDepCoevo](): IN PROGRESS. This tutorial shows the base ideas of the scripts which explores the coevolutionary process in context-dependent interactions. 
You will see that the basic ideias of this scripts and those from the MutAntNet tutorial are pretty similar.


### Running my simulations

My scripts are directed to my own working directory. Please, remember to use your own working directory when run the simulations. The scripts already have an automate way of installing the 
required packages listed below:

List of R packages used:
- [ggplot2](https://ggplot2.tidyverse.org)
- [cowplot](https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html)
- [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)

## Authors

* **Lucas A. Camacho** - *Grad student at the Departamento de Ecologia, IB - USP*
* **Paulo Roberto Guimarães Junior** - *Associate Professor at the Departamento de Ecologia, IB - USP* - [Paulo's website](http://guimaraeslab.weebly.com)

## Acknowledgments

* My supervisor Paulo Guimarães.
* To all the actual and pass students of MiudoLab and Lage's at USP.
* Charles Darwin for dive in that crazy trip in a crazy ship and write about it.
* Alan Turing for develop the base of computers that we have nowadays.
* Electric Octopus, Jimmy Hendrix, Naxatras and several other psychodelic bands that help me code better.
