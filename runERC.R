require(devtools)
#install TreeTools
#remotes::install_github("ms609/TreeTools")
#this takes a few seconds and only needs to be done once per session
source("updates2022.R")
source("ERC.R")
Rcpp::sourceCpp("cppFuncs.cpp")



##put in the file with all gene trees here
comptrees=readTrees("your/tree/here.tre")
comptrees=transformPaths(comptrees, transform = "sqrt",impute = F)
compResid=getAllResiduals(comptrees, n.pcs = 0)
rMat=getRMat(compResid, all = T, weights = comptrees$weights)
clusterList=getClusterList(comptrees)
##Here is where you will edit the threshold you want, minSp is the number of species two genes have to share
##If you only want to run a few genes you can set doOnly = c("genea","geneb")
##If you want the plot of the RERs set plot = T (I would only recommend doing this for a few genes because it uses up a lot of space)
corres=computeERC(rMat, comptrees, parallel = F, clusterListOutput = clusterList,  minSp = 15, saveFile = "your/outfile/.RDS")

