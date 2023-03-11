# A Non-stationary Dependence Model for Extreme European Windstorms

This repository contains R code used in my master thesis. 
This code implements non-stationary dependence in the generalized r-Pareto framework using locally anisotropic SPDEs.
The optimization process uses the recent parallel implementation of Marquardt–Levenberg Algorithm in R.

*[ **The project is not under development anymore, please fork for functionalities improvement**. If you see any bug please submit a pull/merge request.]*

### Documentation

To understand the approach and the implementation you can read the following documents:
* [Final presentation](Documentation/Presentation_02_2021.pdf) - A presentation giving a quick understanding of the project
* [Thesis report](Documentation/Master_Thesis.pdf) - The detailed master thesis report

### Background

The code implementation is mainly based on the following articles:
* de Fondeville, R., and Davison, A. C. (2020), “Functional Peaks-over-threshold Analysis,” arXiv:2002.02711 [stat].
* de Fondeville, R., and Davison, A. C. (2018), “High-dimensional peaks-over-threshold inference,” Biometrika, 105, 575–592. https://doi.org/10.1093/biomet/asy026.
* Fuglstad, G.-A., Lindgren, F., Simpson, D., and Rue, H. (2015), “Exploring a New Class of Non-stationary Spatial Gaussian Random Fields with Varying Local Anisotropy,” Statistica Sinica, 25, 115–133. https://doi.org/10.5705/ss.2013.106w.
* Philipps, V., Hejblum, B. P., Prague, M., Commenges, D., and Proust-Lima, C. (2020), “Robust and Efficient Optimization Using a Marquardt-Levenberg Algorithm with R Package marqLevAlg,” arXiv:2009.03840 [stat].


## Repository structure

The repository contains the following folders:
* Documentation - Contains the thesis report and a presentation on the project.
* Data - Contains windgust data from ERA-Interim model.
* Code - Contains scripts, functions and a package to perform the inference on the data.
* Tmp - Temporary data folder, contains intermidiate results.
* Debug - Default debug directory.

## Getting Started


### Prerequisites

You will need a few libraries to run the scripts.
You can install them with
```
install.packages(c('sp','rgdal','patchwork','parallel','ggplot2','ggmap','EnvStats','evd','Matrix','marqLevAlg','devtools','grid'))
```

To perform parallel optimization using the 'marqLevAlg' package we need to package some functions to send them to parallel workers. This is done with a custom package called 'nonStatInf' located in 'Code/nonStatInf' that is rebuild everytime the inference script is run using 'devtools' functionalities.

### Scripts and functions

The two main scripts are located in 'Code' and will take you through the inference process.

Most of the functions used are located in 'Code/Functions' and a few of them (that need to be packaged) are located in 'Code/nonStatInf/R'

### Pre-computed Results

Some temporary files are already located in the 'Tmp' directory and can be used to save-up on computationally intensive part of the scripts.

## Author and acknowledgements

**Paul Castelain** - <paul.castelain@alumni.epfl.ch>

The implementation of Fuglstad et al. (2015) was kindly provided by Raphaël de Fondeville.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
