# ERC (Evolutionary Rate Covariation)
This software is for determining Evolutionary Rate Covariation (ERC) between sets of genes, using R.

## Input:
To run it, you will need a tree of genes you want to examine

## Output:
The code will output a correlation residual matrix and a matrix with branch values. The residual matrix has numbers representing ERC correlations. A positive value indicates that the genes are above the expected variation, while a negative value represents that they are below the expected variation. A zero means that the genes have the expected evolutionary relationship.


## Installation:
This software is built in R.
Once you download this software, you will need to install its dependencies:
```
devtools
RColorBrewer
gplots
phytools
├──ape
├──maps
├──Rcpp
geiger
knitr
RcppArmadillo
weights
phangorn
```
Once the package is downloaded and its dependencies installed, the main code you will need to run is RunERC.R. For the first runthrough, uncomment line 3 (which will grab code from github), but then comment it back out for the remainder of your session.



## Use Instructions:
First, change the file extensions to match your folders. This may include the other R files in this package, such as Updates2022.R. Make sure you set a valid output file name as well. Finally, you're ready to run the code!


At the end of execution, you will have a data object called "corres" which will also be saved as your output .RDS file.
Use ``corres[["cor"]]`` to see and operate on the correlation matrix. ``corres[["count"]]`` represents the number of observations/branches that went into each correlation.
