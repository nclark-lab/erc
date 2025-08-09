require(devtools)
source("ERC_functions.R")
source("ERC.R")
Rcpp::sourceCpp("cppFuncs.cpp")


##put in the file with all gene trees here
treefile = "physical_interaction_paper/domains_trees.tre"
outputfile = "out.RDS"


## readTrees reads in your file of trees and constructs a master tree based on them
trees=readTrees(treefile)


compTrees = transformPaths(trees, transform = "sqrt",impute = F)
compResid = coreGetResiduals(compTrees, n.pcs=0)
residuals = getRMat(compResid, all = T, rmatweights = compTrees$weights)

clusterList=getClusterList(compTrees)



##Here is where you will edit the threshold you want, minSp is the number of species two genes have to share
##If you only want to run a few genes you can set doOnly = c("genea","geneb")
##If you want the plot of the RERs set plot = T (I would only recommend doing this for a few genes because it uses up a lot of space)

corres=computeERC(residuals, compTrees, parallel = F, clusterListOutput = clusterList,  minSp = 15, saveFile = outputfile)


# We strongly recommend you Fisher transform the data; this is what we usually operate on.
ft_data = fisherTransform(corres)
