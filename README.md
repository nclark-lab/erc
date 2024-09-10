# ERC (Evolutionary Rate Covariation)
This software is for determining Evolutionary Rate Covariation (ERC) between sets of genes, using R.

## Input:
You will need a .tre file containing Newick format trees of the genes you want to examine. Optionally, you can provide a master tree, but if you do not have one the program will make one.

## Output:
The code will output a correlation residual matrix and a matrix with branch values. The residual matrix has numbers representing ERC correlations. A positive value indicates that the genes are above the expected variation, while a negative value represents that they are below the expected variation. A zero means that the genes have the expected evolutionary relationship.


## Installation:
Check our [installation](https://github.com/nclark-lab/erc/blob/main/install.md) page for more details.


## Use Instructions:
First, change the file paths to match your folders. This may include the other R files in this package, such as Updates2022.R. Make sure you set a valid output file name as well. Finally, you're ready to run the code!


At the end of execution, you will have a data object called "corres" which will also be saved as your output .RDS file.
Use ``corres[["cor"]]`` to see and operate on the correlation matrix. ``corres[["count"]]`` represents the number of observations/branches that went into each correlation.


## Quickstart Guide
0. Install dependencies ([visit install page](https://github.com/nclark-lab/erc/blob/main/install.md))
1. Download package (either download from github or use git clone)
2. Download [yeast tree example](https://github.com/nclark-lab/erc/blob/main/physical_interaction_paper/domains_trees.tre)
3. Set output file name on line 20 of RunERC.R
4. Uncomment line 3, then run the program
   -  You may need to be explicit about source file paths if the code gives you an error.
5. Comment line 3 until you next boot up R
6. Manipulate corres matrices (exapmles coming soon.)

For example:

`corres[["cor"]][0:10,0:10]` shows the first 10x10 of the correlation matrix.

`corres[["cor"]][300:310,300:310]` shows elements 300-310 of the axes of the correlation matrix.
