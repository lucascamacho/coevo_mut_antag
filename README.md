# From mutualism to antagonism: coevolutionary dynamics of cheaters and context dependent interactions

My Master's dissertation that explore how the coevolution happens when cheaters shows up in mutualistic networks. Also, i wanna know how the fluctuations in
interactions signals between positive and negative influence the coevolutionary process of species. For that, i use a mathematical model that describes how
a certain average trait of species changes due to mutualism and antagonism. Latter i insert fluctuations in interaction signal's in time.

## Technical specifications

I run my simulations in two different computers:
- Desktop MacPro (Mid 2010) 2x2.4 GHz Quad-Core Intel Xeon, 6GB DDR3, macOS High Sierra 10.13.6
- Notebook MacBook Pro 

The softwares used was:
- R 3.5.0
- Zotero 5.0.61
- TexMaker 5.0.2
- GitHub Desktop 1.6.5

### Repository organization

This repository is organized in folders which contain different aspects of my master's project:

- Code: separate in functions and scripts. You will find all the R files to run my simulations.
- Data: the networks of interactions that i use and the RData files which the results of all my simulations
- Output: figures, pdf's, graphs, etc. All the visual results of my projects are here
- Manuscript: the RMarkdown files that compose my master's dissertation.


### Running my simulations

I deeply recommend to run my scripts looking at the working directory that i'am using. Some parts of my script
depends on R packages that must be previous installed (see the used packages list below) and my scripts are directed
to my own working directory. Please, remember to use your own working directory when run the simulations.

List of R packages used:
- ggplot2
- reshape2
- cowplot
- MASS

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
