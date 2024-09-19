# Function Overview (mk. I)
On this page we will summarize the parameters, inputs, and outputs of the main functions involved in calculating evolutionary rate covariation.
- Input parameters will be in **bold** (`with defaults in parentheses`) and required parameters will be <ins>**underlined.**</ins>


## readTrees
`readTrees` takes a `.tre` input file in newick format, and creates a master tree with branch lengths according to the supplied input trees.
### Input
- <ins>**file**</ins>: the file you want to open
- max.read: (`NULL`) maximum number of trees to read. It is useful for limiting the running time of the program for large tree files
- masterTree: (`NULL`) optional input of all species present in the trees file
### Output
- treesObj object that contains all trees inputted, their branch lengths, and a master tree


## transformPaths
`transformPaths` applies a transformation to the paths of the trees, normalizing them. By default it is set to apply a logarithmic transformation but it can also be set to apply a square root transformation, or no transformation.

### Input
- <ins>**treesObj**</ins>: treesObj object to transform
- **transform**: (`"log"`) transformation to apply. `"log"` by default but can be set to `"sqrt"` or `"none"`
- **impute**: (`T`) whether or not to impute for missing data.
### Output
- A treesObj object, modified as specified


## getAllResiduals
`getAllResiduals` takes a tree object and several optional parameters and returns the residuals of the tree branches.

### Input
- <ins>**treesObj**</ins>: treesObj to get residuals from
- **nvMod** (`NULL`): model or design matrix used for linear regression. If one is not supplied, the function creates one from imputed path lengths. Generally, the user should not supply this, unless they are trying to model with something other than a linear regression. The user should also not be doing this. However, the functionality remains in case they know exactly what they're doing.
- **n.pcs** (`0`): number of principal components for PCA dimensionality reduction
- **cutoff** (`NULL`): Cutoff quantile for determining bad paths (to be thrown out before regression). If left null, it will default to 0.05
- **useSpecies** (`NULL`): list of species to focus on
- **min.sp** (`10`): minimum number of species that must be present in a tree
- **min.valid** (`20`): decides what the minimum score is for a branch to be valid
    - what is the score? not understanding the iigood/iibad/allbranch setup
- **scale** (`F`): unused (?)
- **doOnly** (`NULL`): do only n trees
- **maxT** (`NULL`): Maximum number of trees to use (will default to however many there are in treesObj)
- **block<span>.d</span>o** (`F`): ???
- **weights** (`NULL`): Weights for the paths. If not provided it grabs them from treesObj
- **do.loess** (`F`): toggle for using loess regression. At present unused, but if you know what you're doing you can remove the code preventing its use.
- **family** (`"gaussian"`): Used for the loess regression
- **span** (`0.7`): Also for loess regression
- **interaction**: (`F`): toggle for user interactivity: creates a data frame object

### Output
`getAllResiduals` returns a list with four elements:
- **allresiduals**: residual values
- **index**: ?????
- **preds**: predicted branch values via regression
- **tpaths**: transformed paths

## getRMat
`getRMat` gets creates a residual matrix from an input of residuals.

### Input
- <ins>**resOut**</ins>: input of residual values: the output from `getAllResiduals`.
- **all** (`F`): toggle for using all rows of residual input. (Seems redundant with use.rows?)
- **weights** (`NULL`): weights for residuals
- **scale** (`F`): toggle to T to normalize the residuals (divide by the mean and divide by SD)
- **use.rows** (`NULL`): list of which rows to use. Will default to using all of them.

### Output
- residual matrix

## getClusterList
`getClusterList` takes a treesObj object and returns a list of the trees, clustered (by **???**)

### Input
- <ins>**treesObj**</ins>: tree object to read in data from
### Output
`getClusterList` outputs a list with two elements:
- **idList**: a list of clusters, each cluster containing the IDs of the trees within it
- **count1**: a count of the number of genes seen in each cluster/number of unique tips in the cluster

## runERC

### Input
- <ins>**rr**</ins>:
- <ins>**treesObj**</ins>: Input the trees you want to find ERC for
- **minSp** (`NULL`): minimum species in a tree/cluster
- **doOnly** (`NULL`): list of genes to restrict output to
- **againstAll** (`F`): Whether to compare genes against all genes present (instead of just themselves)
- **parallel** (`F`): Toggle for execution in parallel
- **clusterListOutput** (`NULL`): input for a list of clustered trees. If this is not provided it basically runs getClusterList on its own. (*should we merge them, then? Or at least have one function call the other?*)
- **saveFile** (`NULL`): Filename to save to. If left null(the default) it will not save a file.
- **plot** (`F`): Toggle for displaying a plot.
- **win_value** (`3`):
### Output
`runERC` outputs a list of two matrices:
- A matrix containing ERC values for gene pairs
- A matrix containing number of branches that went into each ERC value
