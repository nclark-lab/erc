library(TreeTools)
library(phangorn)
library(progress)
library(RcppArmadillo)
library(rsvd)
library(impute)
library(pheatmap)
library(Rcpp)
#sourceCpp("cppFuncs.cpp")
computeWeights=function(treesObj, plot=T){
  if(min(treesObj$paths, na.rm = T)<1e-6){
    offset=1e-6
  }
  else{
    offset=0
  }
  mml=apply(log(treesObj$paths+offset),2,mean, na.rm=T)
  varl=apply(log(treesObj$paths+offset),2,var, na.rm=T)
  #use spline fitting
  l=smooth.spline(x = mml, y = varl, w = (exp(mml)), nknots = 6, all.knots = F, spar=0.25)
  f=approxfun(l, rule=2)
  if(plot){
    plot(mml, varl, xlab="mean log", ylab="var log")
    lines(l, lwd=2, col=2)
  }
  #weights are 1/var where the var is the predicted variance
  matrix(1/(f(log(treesObj$paths))), nrow=nrow(treesObj$paths))
}






evalBias=function(x,y,n=20, plot=T){
  
  x=as.vector(x)
  y=as.vector(y)
  
  
  ii=which(!is.na(x)&!is.na(y))
  x=x[ii]
  y=y[ii]
  x=cut(rank(x), breaks = n)
  
  res=boxplot(y~x, outline=F, plot=plot)
  thisMed=res$stats[3,]
  thisSpread=res$stats[4,]-res$stats[2,]
  intStats=mean(thisMed^2)*100
  medStats=summary(lm(thisMed~1+as.vector(1:20)))$coef[2,4]
  spreadStats=summary(lm(thisSpread~1+as.vector(1:20)))$coef[2,4]
  allstats=round(c(intStats,-log10(c(medStats, spreadStats))),2)
  names(allstats)=c("medianInt", "medianSlope", "spreadSlope")
  if(plot){
    abline(h=0, col=2, lwd=2)
    title(paste(paste(names(allstats), allstats), collapse=", "))
  }
  return(allstats)
}


evalResid=function( residOut, treesObj, maxUse=5000, weighted=F, plot=F){
  if(weighted){
    rr=getRMat(residOut, weights = treesObj$weights, use.rows = 1:maxUse)
  }
  else{
    rr=getRMat(residOut, use.rows = 1:maxUse)
  }
  evalBias(treesObj$paths[1:maxUse,], rr[1:maxUse,], plot=plot)
}
copyMat=function(mat, names=T){
  newmat=matrix(nrow=nrow(mat), ncol=ncol(mat))
  if(names){
    rownames(newmat)=rownames(mat)
    colnames(newmat)=colnames(mat)
  }
  newmat
}

allPaths=function(tree, needIndex=T){
  pres=PathLengths(tree)
  
  allD=pres[,3]
  
  nA=length(tree$tip.label)+tree$Nnode
  if(needIndex){
    matIndex=matrix(nrow=nA, ncol=nA)
    for( j in 1:nrow(pres)){
      matIndex[pres[j,2], pres[j,1]]=j
    }
    
    
    return(list(dist=allD, nodeId=pres[,c(2,1)], matIndex=matIndex))
  }
  else{
    return(list(dist=allD, nodeId=pres[,c(2,1)]))
  }
}




matchAllnodes=function(tree, masterTree){
  index = KeptVerts(masterTree, TipLabels(masterTree) %in% tree$tip.label)
  key = which(index)
  map = cbind(seq_along(key), key)
  map
}


edgeIndexRelativeMaster=function(tree, masterTree){
  map=matchAllnodes(tree,masterTree)
  newedge=tree$edge
  newedge[,1]=map[newedge[,1],2]
  newedge[,2]=map[newedge[,2],2]
  newedge
}

allPathsMasterRelative=function(tree, masterTree, masterTreePaths=NULL,i=NULL){
  if(! is.list(masterTreePaths)){
    masterTreePaths=allPaths(masterTree)
  }
  
  treePaths=allPaths(tree, needIndex = F)
  map=matchAllnodes(tree,masterTree)
  
  #remap the nodes
  treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
  treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]
  
  
  #ii=masterTreePaths$matIndex[(treePaths$nodeId[,2]-1)*nrow(masterTreePaths$matIndex)+treePaths$nodeId[,1]]
  ii=masterTreePaths$matIndex[cbind(treePaths$nodeId[,1],treePaths$nodeId[,2])]
  vals=double(length(masterTreePaths$dist))
  vals[]=NA
  if(sum(is.na(ii))>0 & !is.null(i)) {
    message("warning: discordant tree topology in tree ", i,", returning NA row",sep="")
    return(vals)
  }
  vals[ii]=treePaths$dist
  vals
}

readTrees=function(file, max.read=NA, masterTree=NULL, minTreesAll=20, reestimateBranches=F, minSpecs=NULL, useSpecies=NULL){
  message("Using readTrees 2")
  tmp=scan(file, sep="\t", what="character", quiet = T)
  message(paste("Read ",length(tmp)/2, " items", collapse=""))
  trees=vector(mode = "list", length = min(length(tmp)/2,max.read, na.rm = T))
  keeptrees=rep(TRUE,length(trees))
  treenames=character()
  maxsp=0; # maximum number of species
  
  for ( i in 1:min(length(tmp),max.read*2, na.rm = T)){
    if (i %% 2==1){
      treenames=c(treenames, tmp[i])
    }
    else{
      trees[[i/2]]=unroot(read.tree(text=tmp[i]))
      if (i == 2){  #initial setup
        if(!is.null(useSpecies)){allnames=intersect(trees[[i/2]]$tip.label,useSpecies)
        } else {allnames=trees[[i/2]]$tip.label}
      }
      
      #check useSpecies parameter
      if(!is.null(useSpecies)){
        if(length(intersect(trees[[i/2]]$tip.label,useSpecies)) < 3){
          keeptrees[[i/2]] = FALSE; next;
        }
        trees[[i/2]] = unroot(keep.tip(trees[[i/2]],intersect(trees[[i/2]]$tip.label,useSpecies)))
      }
      
      #check if it has more species
      if(length(trees[[i/2]]$tip.label)>maxsp){
        allnames = unique(c(allnames,trees[[i/2]]$tip.label))
        maxsp = length(allnames)
      }
    }
    
  }
  trees = trees[keeptrees]
  treenames = treenames[keeptrees]
  names(trees)=treenames
  treesObj=vector(mode = "list")
  treesObj$trees=trees
  treesObj$numTrees=length(trees)
  treesObj$maxSp=maxsp
  message("Done")
  message(paste("Max number of species is ", maxsp))
  
  report=matrix(nrow=treesObj$numTrees, ncol=maxsp)
  colnames(report)=allnames
  
  rownames(report)=treenames
  for ( i in 1:nrow(report)){
    ii=match(allnames, trees[[i]]$tip.label)
    report[i,]=1-is.na(ii)
    
  }
  treesObj$report=report
  
  
  if(is.null(masterTree)){
    
    
    ii=which(rowSums(report)==maxsp)
    
    #Create a master tree with no edge lengths
    master=trees[[ii[1]]]
    master$edge.length[]=1
  }
  else{
    master=masterTree
  }
  
  
  
  master=Preorder(master)
  treesObj$masterTree=master
  
  for ( i in 1:treesObj$numTrees){
    treesObj$trees[[i]]=RenumberTips(treesObj$trees[[i]], master$tip.label)
    treesObj$trees[[i]]=Preorder(treesObj$trees[[i]])
    
  }
  
  ap=allPathsTT(master)
  treesObj$ap=ap
  matAnc=(ap$matIndex>0)+1-1
  matAnc[is.na(matAnc)]=0
  
  paths=matrix(nrow=treesObj$numTrees, ncol=length(ap$dist))
  
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = treesObj$numTrees,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  message("Extracting paths")
  for( i in 1:treesObj$numTrees){
    pb$tick()
    paths[i,]=allPathsMasterRelativeTT(treesObj$trees[[i]], master, ap,i)
  }
  
  
  #  paths=paths+min(paths[paths>0], na.rm=T)
  treesObj$paths=paths
  treesObj$matAnc=matAnc
  treesObj$matIndex=ap$matIndex
  treesObj$lengths=unlist(lapply(treesObj$trees, function(x){sqrt(sum(x$edge.length^2))}))
  
  #require all species and tree compatibility
  #ii=which(rowSums(report)==maxsp)
  ii=intersect(which(rowSums(report)==maxsp),which(is.na(paths[,1])==FALSE))
  
  
  
  #if masterTree is provided by user, must use minSpecs<maxsp
  #if no user supplied tree and not minSpec, calculate branch lengths from trees with all species
  #if minSpecs<maxsp, calculate branch lengths from trees with minSpecs species
  if(is.null(minSpecs)){
    #if minimum species not specified,
    #minimum is all species
    minSpecs=maxsp
  }
  
  if(!is.null(masterTree) && !reestimateBranches){
    message("Using user-specified master tree")
  }
  
  
  if(minSpecs==maxsp){ #if we're using all species
    if (is.null(masterTree)) { #and if the user did not specify a master tree
      if(length(ii)>=minTreesAll){
        message (paste0("Estimating master tree branch lengths from ", length(ii), " genes"))
        tmp=lapply( treesObj$trees[ii], function(x){x$edge.length})
        
        allEdge=matrix(unlist(tmp), ncol=2*maxsp-3, byrow = T)
        allEdge=scaleMat(allEdge)
        allEdgeM=apply(allEdge,2,mean)
        treesObj$masterTree$edge.length=allEdgeM
      }else {
        message("Not enough genes with all species present: master tree has no edge.lengths")
      }
    }else{
      message("Must specify minSpecs when supplying a master tree: master tree has no edge.lengths")
    }
  }else{ #if we are not using all species
    #estimating from trees with minimum number of species
    treeinds=which(rowSums(report)>=minSpecs) #which trees have the minimum species
    message (paste0("estimating master tree branch lengths from ", length(treeinds), " genes"))
    
    
    if(length(treeinds)>=minTreesAll){
      pathstouse=treesObj$paths[treeinds,] #get paths for those trees
      
      colnames(pathstouse) = ap$destinNode
      colBranch = vector("integer",0)
      unq.colnames = unique(colnames(pathstouse))
      
      for (i in 1:length(unq.colnames)){
        ind.cols = which(colnames(pathstouse) == unq.colnames[i])
        colBranch = c(colBranch,ind.cols[1])
      }
      
      allEdge = pathstouse[,colBranch]
      allEdgeScaled = allEdge
      for (i in 1:nrow(allEdgeScaled)){
        allEdgeScaled[i,] = scaleDistNa(allEdgeScaled[i,])
      }
      colnames(allEdgeScaled) = unq.colnames
      
      edgelengths = vector("double", ncol(allEdgeScaled))
      
      edge.master = treesObj$masterTree$edge
      
      for (i in 1:nrow(edge.master)){
        destinNode.i = edge.master[i,2]
        col.Node.i = allEdgeScaled[,as.character(destinNode.i)]
        edgelengths[i] = mean(na.omit(col.Node.i))
      }
      
      treesObj$masterTree$edge.length = edgelengths
    }else{
      message("Not enough genes with minSpecs species present: master tree has no edge.lengths")
    }
  }
  
  
  message("Naming columns of paths matrix")
  colnames(treesObj$paths)=namePathsWSpeciesTT(treesObj)
  
  
  class(treesObj)=append(class(treesObj), "treesObj")
  treesObj
}

allPathsMasterRelativeTT =function (tree, masterTree, masterTreePaths=NULL,i=NULL){
  
  if(! is.list(masterTreePaths)){
    masterTreePaths=allPathsTT(masterTree)
  }
  
  treePaths=allPathsTT(tree, needIndex = F)
  map=matchAllnodesTT(tree,masterTree)
  
  #remap the nodes
  treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
  treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]
  
  
  #ii=masterTreePaths$matIndex[(treePaths$nodeId[,2]-1)*nrow(masterTreePaths$matIndex)+treePaths$nodeId[,1]]
  ii=masterTreePaths$matIndex[cbind(treePaths$nodeId[,1],treePaths$nodeId[,2])]
  vals=double(length(masterTreePaths$dist))
  
  if(sum(is.na(ii))>0 & !is.null(i)) {
    message("warning: discordant tree topology in tree ", i,", returning NA row",sep="")
    return(vals)
  }
  
  vals[]=NA
  vals[ii]=treePaths$dist
  vals
}
namePathsWSpeciesTT =function (treesObj){
  cnames=vector("character", ncol(treesObj$paths))
  
  for(i in 1:ncol(treesObj$paths)){
    tip=which(treesObj$matIndex==i, arr.ind = T)[,1]
    #  show(tip)
    if(tip<=treesObj$maxSp){
      cnames[i]=treesObj$masterTree$tip.label[tip]
    }
  }
  cnames
}
allPathsTT =function (tree, needIndex=T){
  pres=PathLengths(tree)
  
  allD=pres[,3]
  
  nA=length(tree$tip.label)+tree$Nnode
  if(needIndex){
    matIndex=matrix(nrow=nA, ncol=nA)
    for( j in 1:nrow(pres)){
      matIndex[pres[j,2], pres[j,1]]=j
    }
    
    
    return(list(dist=allD, nodeId=pres[,c(2,1)], matIndex=matIndex))
  }
  else{
    return(list(dist=allD, nodeId=pres[,c(2,1)]))
  }
}
matchAllnodesTT =function (tree, masterTree){
  index = KeptVerts(masterTree, TipLabels(masterTree) %in% tree$tip.label)
  key = which(index)
  map = cbind(seq_along(key), key)
  map
}
edgeIndexRelativeMasterTT =function (tree, masterTree){
  map=matchAllnodesTT(tree,masterTree)
  newedge=tree$edge
  newedge[,1]=map[newedge[,1],2]
  newedge[,2]=map[newedge[,2],2]
  newedge
}


scaleMat=function(mat){t(apply(mat,1,scaleDist))}
scaleDist=function(x, na.rm=T){
  x/sqrt(sum(x^2))
}





transformMat=function(tree){
  
  nA=length(tree$tip.label)+tree$Nnode
  matIndex=matrix(nrow=nA, ncol=nA)
  mat=matrix(nrow=nrow(tree$edge), ncol=0)
  index=1
  for ( i in 1:nA){
    ia=getAncestors(tree,i)
    if(length(ia)>0){
      thisindex=double()
      ia=c(i,ia)
      for (k in 2:length(ia)){
        j=ia[k]
        
        thisindex=c(thisindex, which(tree$edge[,2]==ia[k-1]&tree$edge[,1]==ia[k]))
        
        vals=rep(0, nrow(mat))
        vals[thisindex]=1
        mat=cbind(mat, vals)
      }
    }
  }
  mat
}




# internal
getAncestors=function(tree, nodeN){
  if(is.character(nodeN)){
    nodeN=which(tree$tip.label==nodeN)
  }
  im=which(tree$edge[,2]==nodeN)
  ianc=tree$edge[im,1]
  ancVec=double()
  while(length(ianc)>0){
    ancVec=c(ancVec,ianc)
    im=which(tree$edge[,2]==ianc)
    ianc=tree$edge[im,1]
  }
  ancVec
}


getLogPaths=function(treesObj){
  set.seed(1);  iis=sample(nrow(treesObj$paths),500)
  tmpdat=treesObj$paths[iis,]
  mval=min(tmpdat[tmpdat>0], na.rm = T)
  log(treesObj$paths+mval)
}

#various transformation functions, using log or arc hyperbolic sine or square root
transformData=function(x, transform=c("log", "asinh")){
  if(transform=="log"){
    if(length(x)>10000){
      set.seed(1)
      iis=sample(length(x), 10000)
      tmpdat=as.vector(x[iis])
      mval=min(tmpdat[tmpdat>0], na.rm = T)
    }
    else{
      mval=min(x[x>0], na.rm=T)
    }
    log(x+mval)
  }
  else{
    asinh(x)
  }
}

transformPaths=function(treesObj, transform="sqrt", impute=T){
  transform=match.arg(transform, c("sqrt", "log", "none"))
  if(transform=="log"){
    set.seed(1);  iis=sample(nrow(treesObj$paths),500)
    tmpdat=treesObj$paths[iis,]
    mval=min(tmpdat[tmpdat>0], na.rm = T)
    treesObj$paths= log(treesObj$paths+mval)
    treesObj$logOffset=mval
    treesObj$transform="log"
  }
  else if (transform=="sqrt"){
    treesObj$paths= sqrt(treesObj$paths)
    treesObj$transform="sqrt"
  }    
  
  else if (transform=="none"){
    
    treesObj$transform="none"
  }    
  
  
  
  if(impute){
    set.seed(1);kres=impute.knn(treesObj$paths, colmax = 95)
    treesObj$pathsImputed=kres$data
    
  }
  treesObj$weights=computeWeightsAllVar(treesObj$paths, transform = "none")
  treesObj
  
}

unlogPaths=function(treesObj){
  if(is.null(treesObj$isLogged)){
    message("paths are not log transformed")
  }
  else{
    exp(treesOb$paths)-treesObj$logOffset
  }
}

getTrimmedAverage=function(x, trim=0.05){
  apply(x,2, mean, trim=trim, na.rm=T)
}

coreGetResiduals=function(treesObj, nvMod=NULL, n.pcs=0,cutoff=NULL, useSpecies=NULL,  min.sp=10, min.valid=20,  doOnly=NULL, maxT=NULL, block.do=F, residweights=NULL, do.loess=F, family="gaussian", span=0.7, interaction=F){
  
  
  
  if (is.null(cutoff)) {
    cutoff = quantile(treesObj$paths, 0.05, na.rm = T)
    message(paste("cutoff is set to", cutoff))
  }
  
  if(is.null(treesObj$pathsImputed) && n.pcs>0){
    message("PC normalization not supported without imputation\n
            Run transformPaths with imputation or set n.pcs=0")
    return()
  }
  if(block.do){
    pStr=apply(treesObj$report, 1, paste, collapse="")
    uClust=unique(pStr)
    blockId=match(pStr, uClust)
  }

  tPaths= treesObj$paths
  if(is.null(nvMod)){
    if(n.pcs>0){
      # dimensionality reduction
      svdres=rsvd(treesObj$pathsImputed, k=n.pcs+1)

      if(!interaction){
        nvMod=model.matrix(~1+svdres$v[, 1:n.pcs])
      }
      else{
        nvMod=model.matrix(~1+.^2, data = as.data.frame(svdres$v[, 1:n.pcs]))
      }

    }
    else {
      nvAve=apply(tPaths,2, mean, trim=0.05, na.rm=T)
      nvMod=model.matrix(~1+nvAve)
    }
  }
  if(any(is.na(nvMod))){
    message("NA values in model: something is wrong")
    stop()
  }
  if(ncol(nvMod)>5 && do.loess){
    message("Model is too big for loess. Set n.pcs<=4 or n.pcs<=2 if interaction=T")
    stop()
  }
  
  if (is.null(useSpecies)){
    useSpecies=treesObj$masterTree$tip.label
    mappedEdges=treesObj$mappedEdges
  }
  if(is.null(maxT)){
    maxT=treesObj$numTrees
  }
  
  
  if(is.null(residweights)){
    residweights=treesObj$weights
  }
  
  
  #useSpecies has to exist in masterTree
  useSpecies=intersect(treesObj$masterTree$tip.label, useSpecies)
  
  #maximum number of present species
  maxSpecies=rowSums(treesObj$report[,useSpecies])
  
  #this will hold the predictions  
  preds=copyMat(tPaths)
  preds[]=NA
  
  
  if(is.null(doOnly)){
    doOnly=1
  }
  else{
    maxT=1
  }
  isDone=vector("logical", nrow(preds))
  isDone[]=F
  
  #this will hold the indecies used to construct the regression
  #which is a subset of the indecies for predictions
  modelIndexList=vector("list", nrow(preds))
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = treesObj$numTrees,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  for (i in doOnly:(doOnly+maxT-1)){
    
    pb$tick()
    
    if(!isDone[i]){
      
      
      #get the ith tree
      tree1=treesObj$trees[[i]]
      
      #get the common species, prune and unroot
      thisUseSpecies=intersect(tree1$tip.label, useSpecies)
      
      if(length(thisUseSpecies)<min.sp){
        
        next
      }

      
      #find all the genes that that whose maximal species set is the same as tree1
      
      if(block.do){
        #      thisMaxSpecies=rowSums(treesObj$report[,thisUseSpecies])
        #     iido=which(maxSpecies==thisMaxSpecies&thisMaxSpecies==length(thisUseSpecies))
        iido=which(blockId==blockId[i])
      }
      else{
        iido=i
      }

      ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)
      
      #these are the indecies for the paths that span the tree 
      #these are used to build the regression
      
      iiPaths= treesObj$matIndex[ee[, c(2,1)]]


      #extract the tPaths and corresponding weights
      
      allbranch=tPaths[iido,iiPaths,drop=F]
      allbranchw=residweights[iido,iiPaths, drop=F]
      iibad=which(tPaths[iido, iiPaths, drop=F]<cutoff)
      
      allbranch[iibad]=NA
      
      iidoAll=iido
      sumCanUse=rowSums(!is.na(allbranch))
      iigood=which(sumCanUse>=min.valid)
      
      if(length(iigood)<1){
        next
      }
      allbranch=allbranch[iigood, ,drop=F]
      allbranchw=allbranchw[iigood, ,drop=F]
      iido=iido[iigood]
      
      #these are all the indecies that are not NA and thus can get predictions
      iiallnotna=which(!is.na(treesObj$paths[iido[1],]))
      
      
      
      for(j in 1:length(iido)){
        
        modelIndexList[[iido[j]]]=iiPaths[which(!is.na(allbranch[j,]))]
      }

      
      if(!do.loess){
        
        preds[iido,iiallnotna]= fastLmResidMatWeightedPredict(allbranch, nvMod[iiPaths,], allbranchw, nvMod[iiallnotna,])
        
      }
      else{ #this is loess
        message("loess should not be used since it produces weird results and also cannot predict outside the input range")
        stop()
        for(id in 1:length(iido)){
          lres=loess(allbranch[id,]~nvMod[iiPaths,-1], weights =allbranchw[id,], span=span, family=family)
          preds[iido[id],iiallnotna]=predict(lres, nvMod[iiallnotna,-1])
        }
      }
      #plotting function. change to T for plotting functionality
      if(F){
        plot(preds[iido[1],], tPaths[iido[1],])
        
        abline(a=0,b=1, lwd=3, col=2)
      }
      isDone[iidoAll]=T
    }
  }
  allresiduals=(tPaths-preds)
  
  rownames(allresiduals)=names(treesObj$trees)
  #  allresiduals[treesObj$paths<cutoff]=0
  return(list(allresiduals=allresiduals, index=modelIndexList, preds=preds, tpaths=tPaths))
  
  
}


#' Provides names for paths/RERs representing terminal branches for plotting
#' Originally an internal function but necessary for the vignette/walk-through
#' @param  masterTree The master tree used for analysis
#' @return  Names corresponding to the paths/RERs for terminal branches
#' @export
namePathsWSpecies=function(masterTree){
  mat=transformMat(masterTree)
  n=length(masterTree$tip.label)
  #these are the tip edges in the master tree
  iim=match(1:n, masterTree$edge[,2])
  #each column in the mat is composed of at most one tip edge
  tip.edge=apply(mat[iim,],2,function(x){if(max(x)>0){which(x==1)} else{NA}})
  return(masterTree$tip[tip.edge])

}


flatIndex=function(twocol,n){
  
  (twocol[,2]-1)*n+twocol[,1]
}


myscale=function(x ){
  s = apply(x, 2, sd, na.rm=T)
  m = apply(x, 2, mean, na.rm=T)
  x = sweep(x, 2, m)
  x = sweep(x, 2, s, "/")
  
  x
}

getRMat=function(resOut, all=F, rmatweights=NULL, scale=F, use.rows=NULL){
  
  allres=copyMat(resOut$allresiduals)
  if(is.null(use.rows)){
    use.rows=1:nrow(allres)
  }
  resIn=resOut$allresiduals
  if(!is.null(rmatweights)){
    resIn=resIn*sqrt(rmatweights)
  }
  if(scale){
    resIn=myscale(resIn)
  }
  
  if(!all){
    for(i in use.rows){
      
      ii=resOut$index[[i]]
      allres[i,ii]=resIn[i,ii]
    }
    allres
  }
  else{
    allres=resIn
  }
  
  
  allres
}


getAllResiduals=function(treesObj, transform="sqrt", impute=T,  # transformPaths arguments
                         #--------------------------------------------------------#
                         nvMod=NULL, n.pcs=0, cutoff=NULL,      #
                         useSpecies=NULL,  min.sp=10,           #
                         min.valid=20, doOnly=NULL, maxT=NULL,  # getAllResiduals core function arguments
                         block.do=F, residweights=NULL,              #
                         do.loess=F, family="gaussian",         #
                         span=0.7, interaction=F,               #
                         #--------------------------------------------------------#
                         all=F, rmatweights=NULL, scale=F, use.rows=NULL)   # getRMat arguments
{  
  
  
  tree2 = transformPaths(treesObj, transform = transform, impute = impute)
  
  resids = coreGetResiduals(tree2, nvMod=nvMod, n.pcs=n.pcs, cutoff=cutoff,
                    useSpecies=useSpecies, min.sp=min.sp, min.valid=min.valid,
                    doOnly=doOnly,maxT=maxT, block.do=block.do, residweights=residweights,
                    do.loess=do.loess, family=family, span=0.7, interaction=interaction)
  if(is.null(rmatweights)){rmatweights=tree2$weights}
  r = getRMat(resids, all=all,rmatweights=rmatweights, scale=scale,use.rows=use.rows)
  
  return(r)
}

win=function(x,w){
  xs=sort(x[!is.na(x)], decreasing = T)
  xmax=xs[w]
  xmin=xs[length(xs)-w+1]
  
  x[x>xmax]=xmax
  x[x<xmin]=xmin
  x
}

OldgetAllCor=function(RERmat, charP, method="auto",min.sp=10, min.pos=2, winsorize=NULL,weights=NULL){
  if (method=="auto"){
    lu=length(unique(charP))
    if(lu==2){
      method="k"
      message("Setting method to Kendall")
    }
    else if (lu<=5){
      method="s"
      message("Setting method to Spearman")
    }
    else{
      method="p"
      message("Setting method to Pearson")
      if(is.null(winsorize)){
        message("Setting winsorise=3")
      }
    }
  }
  
  corout=matrix(nrow=nrow(RERmat), ncol=3)
  rownames(corout)=rownames(RERmat)
  
  colnames(corout)=c("Rho", "N", "P")
  if(!is.null(winsorize)){
    charP=win(charP, winsorize)
  }
  
  for( i in 1:nrow(corout)){
    
    if(((nb<-sum(ii<-(!is.na(charP)&!is.na(RERmat[i,]))))>=min.sp)&sum(charP[ii]!=0)>=min.pos){
      
      if(is.null(weights)){
        
        
        if (!is.null(winsorize)){
          x=win(RERmat[i,], winsorize)
        }
        else{
          x=RERmat[i,]
        }
        cres=cor.test(x, charP, method=method)
        corout[i,1:3]=c(cres$estimate, nb, cres$p.value)
      }
      else{
        
        cres=wtd.cor(rank(RERmat[i,ii]), rank(charP[ii]), weight = weights[rownames(RERmat)[i],ii], mean1 = T)
        
        corout[i, 1:3]=c(cres[1], nb, cres[4])
      }
    }
    else{
      #show(i)
      #show(c(nb, charP[ii]))
    }
    
  }
  as.data.frame(corout)
}

mypr=function(vals, genes, foldChange=T){
  ii=match(names(vals), genes)
  fg=(!is.na(ii))+1-1
  oo=order(-vals)
  fgo=fg[oo]
  prc=matrix(nrow=length(fgo), ncol=2)
  rec=0
  for (i in 1:length(fgo)){
    if(fgo[i]==1){
      rec=rec+1
    }
    prc[i,1]=rec
    prc[i,2]=rec/i
  }
  
  
  prc[,2]=rev(cummax(rev(prc[,2])))
  if(foldChange){
    message("FoldChange")
    bg=sum(fg)/length(vals)
    prc[,2]=prc[,2]/bg
    colnames(prc)=c("x","y")
  }
  else{
    colnames(prc)=c("x","y")
  }
  
  as.data.frame(prc)
}

plotAllPRggplot=function(resAll, genes, subset=1:length(resAll), grp1, grp2=NULL, grp3=NULL, invert=F, sz=1.7){
  grp1=as.character(grp1)
  if(!is.null(grp2)){
    grp2=as.character(grp2)
  }
  else{
    grp2=as.factor(rep(1, length(grp1)))
  }
  
  if(!is.null(grp3)){
    grp3=as.character(grp3)
  }
  
  
  df=list()
  
  for(i in 1:length(resAll)){
    if(i  %in% subset){
      show(i)
      
      stats=getStat(resAll[[i]])
      if(invert){
        pr=mypr(-stats, genes) 
      }
      else{
        pr=mypr(stats, genes)
      }
      ii=which(pr$x>=10)
      pr=pr[ii,]
      df$x=c(df$x, pr$x)
      df$y=c(df$y,pr$y)
      df$grp1=c(df$grp1, rep(grp1[i], length(pr$x)))
      df$grp2=c(df$grp2, rep(grp2[i], length(pr$x)))
      df$grp3=c(df$grp3, rep(grp3[i], length(pr$x)))
      df$id=c(df$id, rep(i, length(pr$x)))
    }
  }
  
  df=as.data.frame(df)
  # show(df[1:10,])
  #  show(table(df$grp1))
  df$id=as.factor(df$id)
  df$grp1=as.factor(df$grp1)
  df$grp2=as.factor(df$grp2)
  if(length(unique(df$grp2))>1){
    p=ggplot(df, aes(x=x, y=y,color=grp1, linetype=grp2, group=id))+geom_line(size=sz)+scale_x_log10()+theme_bw()
  }
  else{
    p=ggplot(df, aes(x=x, y=y,color=grp1,  group=id))+geom_line(size=sz)+scale_x_log10()+theme_bw()
  }
  
  if(!is.null(grp3)){
    p=p+facet_wrap(~grp3)
  }
  p+theme_bw(base_size = 20)+xlab("recall")+ylab("precision fold change")
}

naresidCPP=function(data, mod, weights=NULL){
  # message("CPP")
  if(is.null(weights)){
    out=fastLmResidMat(data, mod)
  }
  else{
    out=fastLmResidMatWeighted(data,mod, weights)
    out=out*sqrt(weights)
  }
  rownames(out)=rownames(data)
  colnames(out)=colnames(data)
  out
}

# myRfastNorm=function(x){
#   mat <- (x) - Rfast::rowmeans(x)
#   mat <- mat/sqrt(Rfast::rowsums(mat^2))
# mat
#   }
# mycor=function (x, y){
#   if(nrow(x)>100&&nrow(y)>100){
#   matx <- myRfastNorm(x)
# maty=myRfastNorm(y)
# 
#     Rfast::Tcrossprod(matx, maty)
# }
#   else cor(t(x), t(y))
# }







computeWeightsAllVarOld=function (mat, nv = NULL, transform = "none", plot = T, predicted = T) 
{
  if (is.null(nv)) {
    nv = apply(mat, 2, mean, na.rm = T, trim = 0.05)
  }
  transform = match.arg(transform, choices = c("none", "sqrt", 
                                               "log"))
  if (transform == "log") {
    offset = 0
    if (min(mat, na.rm = T) < 1e-08) {
      offset = min(mat[mat > 1e-08])
    }
    mat = log(mat + offset)
    nv = log(nv + offset)
  }
  if (transform == "sqrt") {
    mat = sqrt(mat)
    nv = sqrt(nv)
  }
  matsub = mat
  matr = naresidCPP(matsub, model.matrix(~1 + nv))
  matpred = fastLmPredictedMat(matsub, model.matrix(~1 + nv))
  mml = as.vector(matsub)
  varlall=varl = as.vector(log(matr^2))
  ii = which(!is.na(mml))
  mml = mml[ii]
  varl = varl[ii]
  set.seed(123)
  iis = sample(length(mml), min(5e+05, length(mml)))
  mml = mml[iis]
  varl = varl[iis]
  l = lowess(mml, varl, f = 0.7, iter = 2)
  f = approxfun(l, rule = 2)
  if (plot) {
    par(mfrow = c(1, 2), omi = c(1, 0, 0, 0))
    nbreaks = 20
    qq = quantile(mml, seq(0, nbreaks, 1)/nbreaks)
    qqdiff = diff(qq)
    breaks = qq[1:nbreaks] + qqdiff/2
    rr = quantile(mml, c(1e-04, 0.99))
    breaks = unique(round(breaks, 3))
    nbreaks = length(breaks)
    cutres <- cut(mml, breaks = breaks)
    cutres_tt = table(cutres)
    boxplot((varl) ~ cutres, xlab = "", ylab = "log var", 
            outline = F, log = "", las = 2)
    title("Before")
    xx = (qq[1:nbreaks] + breaks)/2
    lines(1:length(xx), (f(qq[1:nbreaks])), lwd = 2, col = 2)
  }
  wr = 1/exp(f(mml))
  if (!predicted) {
    weights = (matrix(1/exp(f(mat)), nrow = nrow(mat)))
  }
  else {
    weights = (matrix(1/exp(f(matpred)), nrow = nrow(mat)))
  }
  
  if (plot) {
    matr = naresidCPP(matsub[1:2000,], model.matrix(~1 + nv), weights[1:2000,])
    iisub=ii[iis]
    iisub=iisub[iisub<=length(matr)]
    varl = (as.vector(log(matr^2))[ii])[iis]
    boxplot((varl) ~ cutres, ylab = "log var", outline = F, 
            log = "", main = "After", las = 2)
    abline(h = 0, col = "blue3", lwd = 2)
    mtext(side = 1, text = "bins", outer = T, line = 2)
  }
  save(matr, weights, matsub, ii, iis, nv, file="matrOld.Rdata")
  weights
}

naresid=function(data, X,  weights=NULL, covar=NULL, numiter=0, useZero=F){
  if(is.vector(X)){
    mod=model.matrix(~1+as.vector(X));
  }
  else{
    mod=X
  }
  #show(mod)
  resid=matrix(nrow=nrow(data), ncol=ncol(data))
  if(useZero){
    resid[]=0;
  }
  else{
    resid[]=NA;
  }
  for ( i in 1:nrow(data)){
    #  show(i)
    ii=which(!is.na(data[i,]))
    if(length(ii)>2){
      iiinv=which(is.na(data[i,]))
      
      dat=data[i,ii,drop=F]
      
      modtmp=mod[ii,]
      if(!is.null(weights)){
        W=diag(weights[i,ii])
      }
      else{
        W=NULL
      }
      if (!is.null(covar)){
        if(!is.null(W)){
          W=W%*%covar
        }
        else{
          W=covar
        }
      }
      n=dim(dat)[2]
      Id=diag(n)
      if(!is.null(W)){
        #  message("here")
        coeff=dat%*%W%*%modtmp %*% solve(t(modtmp) %*% W %*% modtmp)
        #check there is no error with lm
        # lmres=lm(t(dat)~0+modtmp, weights = diag(W))
        #  show(coeff)
        #  show(coefficients(lmres))
        resid[i, ii] = dat -(coeff %*% t(modtmp))
        resid[i, ii]=resid[i,ii]*sqrt(diag(W))
        if(numiter>0){
          for(iter in 1:numiter){
            ii.use=isNotOutlier(resid[i,], quant)
            Wq=diag(weights[i, ii.use])
            modq=mod[ii.use,]
            coeff=dat[i, ii.use]%*%Wq%*%modq %*% solve(t(modq) %*% Wq %*% modq)
            
            resid[i, ] = dat[i,] -(coeff %*% t(mod))
            resid[i, ]=resid[i,]*sqrt(diag(W))
          }
        }
        
      }
      else{
        coeff=dat%*%modtmp %*% solve(t(modtmp)  %*% modtmp)
        #check there is no error with lm
        # lmres=lm(t(dat)~0+modtmp)
        #  show(coeff)
        #  show(coefficients(lmres))
        resid[i, ii] = dat -(coeff %*% t(modtmp))
        
      }
      
      
      
    }
    else{
      # message("Cannot compute residuals")
      
    }
    
  }
  rownames(resid)=rownames(data)
  colnames(resid)=colnames(data)
  resid
}



plotAsBoxMeanVar=function(mat,matr){
  mml = as.vector(mats)
  varl = as.vector(log(matr^2))
  ii = which(!is.na(mml))
  mml = mml[ii]
  varl = varl[ii]
  set.seed(123)
  iis = sample(length(mml), min(5e+05, length(mml)))
  mml = mml[iis]
  varl = varl[iis]
  plotAsBox(mml, varl, 20)
}


plotAsBox=function(x,y,n=20,intercept=0, ...){
  
  x=as.vector(x)
  y=as.vector(y)
  #ii=which(!is.na(x)&!is.na(y))
  #x=x[ii]
  #y=y[ii]
  
  if(length(x)>1e5){
    set.seed(123);
    iis=sample(length(x), 1e5)
    qq = quantile(x[iis], seq(0, n, 1)/n, na.rm = T)
  }
  else{
    qq = quantile(x, seq(0, n, 1)/n, na.rm = T)
  }
  qqdiff = diff(qq)
  breaks = qq[1:n] + qqdiff/2
  
  breaks = unique(round(breaks, 3))
  
  x=cut(x, breaks=breaks)
  boxplot(y~x, outline=F, las=2, ...)
  abline(h = intercept, col = "blue3", lwd = 2)
}

getNV=function(mat){
  nv = apply(mat, 2, mean, na.rm = T, trim = 0.05)
}


computeWeightsAllVarOriginal=function (mat, nv = NULL, transform = "none", plot = T, predicted = T) 
{
  if (is.null(nv)) {
    nv = apply(mat, 2, mean, na.rm = T, trim = 0.05)
  }
  transform = match.arg(transform, choices = c("none", "sqrt", 
                                               "log"))
  if (transform == "log") {
    offset = 0
    if (min(mat, na.rm = T) < 1e-08) {
      offset = min(mat[mat > 1e-08])
    }
    mat = log(mat + offset)
    nv = log(nv + offset)
  }
  if (transform == "sqrt") {
    mat = sqrt(mat)
    nv = sqrt(nv)
  }
  matsub = mat
  #matr = naresidCPP(matsub, model.matrix(~1 + nv))
  message("computing unweighted predictions")
  matpred = fastLmPredictedMat(matsub, model.matrix(~1 + nv))
  message("Done")
  matr=matsub-matpred
  mml = as.vector(matsub)
  varl = as.vector(log(matr^2))
  ii = which(!is.na(mml))
  mml = mml[ii]
  varl = varl[ii]
  set.seed(123)
  iis = sample(length(mml), min(1e+07, length(mml)))
  mml = mml[iis]
  varl = varl[iis]
  show(length(ii))
  show(ii[iis[1:10]])
  show(varl[1:10])
  show(mml[1:10])
  
  l = lowess(mml[!is.na(varl)], varl[!is.na(varl)], f = 0.7, iter = 2)
  f = approxfun(l, rule = 2)
  if (plot) {
    par(mfrow = c(1, 2), omi = c(1, 0, 0, 0))
    nbreaks = 20
    qq = quantile(mml, seq(0, nbreaks, 1)/nbreaks)
    qqdiff = diff(qq)
    breaks = qq[1:nbreaks] + qqdiff/2
    rr = quantile(mml, c(1e-04, 0.99))
    breaks = unique(round(breaks, 3))
    nbreaks = length(breaks)
    cutres <- cut(mml, breaks = breaks)
    cutres_tt = table(cutres)
    boxplot((varl) ~ cutres, xlab = "", ylab = "log var", 
            outline = F, log = "", las = 2)
    title("Before")
    xx = (qq[1:nbreaks] + breaks)/2
    lines(1:length(xx), (f(qq[1:nbreaks])), lwd = 2, col = 2)
  }
  wr = 1/exp(f(mml))
  if (!predicted) {
    weights = (matrix(1/exp(f(mat)), nrow = nrow(mat)))
  }
  else {
    weights = (matrix(1/exp(f(matpred)), nrow = nrow(mat)))
  }
  if (plot) {
    mmltmp=(as.vector(matsub)[ii])[iis]
    show(all(mml==mmltmp))
    
    matr = naresidCPP(matsub, model.matrix(~1 + nv), weights)
    varl = (as.vector(log(matr^2))[ii])[iis]
    
    
    cutres <- cut(mml, breaks = breaks)
    boxplot((varl) ~ cutres, ylab = "log var", outline = F, 
            log = "", main = "After", las = 2)
    abline(h = 0, col = "blue3", lwd = 2)
    mtext(side = 1, text = "bins", outer = T, line = 2)
  }
  weights
}


computeWeightsAllVar=function (mat, nv = NULL, transform = "none", plot = T, predicted = T) 
{
  message("new function")
  if (is.null(nv)) {
    nv = apply(mat, 2, mean, na.rm = T, trim = 0.05)
  }
  transform = match.arg(transform, choices = c("none", "sqrt", 
                                               "log"))
  if (transform == "log") {
    offset = 0
    if (min(mat, na.rm = T) < 1e-08) {
      offset = min(mat[mat > 1e-08])
    }
    mat = log(mat + offset)
    nv = log(nv + offset)
  }
  if (transform == "sqrt") {
    mat = sqrt(mat)
    nv = sqrt(nv)
  }
  matsub = mat
  #matr = naresidCPP(matsub, model.matrix(~1 + nv))
  message("computing unweighted predictions")
  matpred = fastLmPredictedMat(matsub, model.matrix(~1 + nv))
  message("Done")
  matr=matsub-matpred
  mml = as.vector(matsub)
  varl = as.vector(log(matr^2))
  ii = which(!is.na(mml))
  mml = mml[ii]
  varl = varl[ii]
  set.seed(123)
  iis = sample(length(mml), min(1e+07, length(mml)))
  mml = mml[iis]
  varl = varl[iis]
  # show(length(ii))
  #  show(ii[iis[1:10]])
  #  show(varl[1:10])
  #  show(mml[1:10])
  
  l = lowess(mml[!is.na(varl)], varl[!is.na(varl)], f = 0.7, iter = 2)
  f = approxfun(l, rule = 2)
  if (plot) {
    par(mfrow = c(1, 2), omi = c(1, 0, 0, 0))
    nbreaks = 100
    qq = quantile(mml, seq(0, nbreaks, 1)/nbreaks)
    qqdiff = diff(qq)
    breaks = qq[1:nbreaks] + qqdiff/2
    rr = quantile(mml, c(1e-04, 0.99))
    breaks = unique(round(breaks, 3))
    nbreaks = length(breaks)
    cutres <- cut(mml, breaks = breaks)
    cutres_tt = table(cutres)
    boxplot((varl) ~ cutres, xlab = "", ylab = "log var", 
            outline = F, log = "", las = 2)
    title("Before")
    xx = (qq[1:nbreaks] + breaks)/2
    lines(1:length(xx), (f(qq[1:nbreaks])), lwd = 2, col = 2)
  }
  wr = 1/exp(f(mml))
  if (!predicted) {
    weights = (matrix(1/exp(f(mat)), nrow = nrow(mat)))
  }
  else {
    weights = (matrix(1/exp(f(matpred)), nrow = nrow(mat)))
  }
  if (plot) {
    message("computing weighted residuals on a subset")
    if(nrow(matsub)<2000){
      set.seed(1);iis=1:nrow(matsub)
    }
    else{
      iis=1:2000
    }
    
    matin=matsub[iis,]
    win=weights[iis,]
    matr = naresidCPP(matin, model.matrix(~1 + nv), weights[iis,])
    save(matr, matin, nv, win, file="matr.RData")
    message("done")
    #    varl = (as.vector(log(matr^2))[ii])[iis]
    #   boxplot((varl) ~ cutres, ylab = "log var", outline = F, 
    #         log = "", main = "After", las = 2)
    plotAsBox( matsub[iis,], log(matr^2), main="After", n = 100)
    abline(h = 0, col = "blue3", lwd = 2)
    mtext(side = 1, text = "bins", outer = T, line = 2)
  }
  weights
}
correlateWithBinaryPhenotype=function (RERmat, charP, min.sp = 10, min.pos = 2, weighted = "auto") 
{
  if (weighted == "auto") {
    if (any(charP > 0 & charP < 1, na.rm = TRUE)) {
      message("Fractional values detected, will use weighted correlation mode")
      weighted = T
    }
    else {
      weighted = F
    }
  }
  getAllCor(RERmat, charP, min.sp, min.pos, method = "k", weighted = weighted)
}
getAllCor=function (RERmat, charP, method = "auto", min.sp = 10, min.pos = 2, 
                    winsorizeRER = NULL, winsorizetrait = NULL, weighted = F) 
{
  RERna = (apply(is.na(RERmat), 2, all))
  iicharPna = which(is.na(charP))
  if (!all(RERna[iicharPna])) {
    warning("Species in phenotype vector are a subset of the those used for RER computation. For best results run getAllResiduals with the useSpecies")
  }
  if (method == "auto") {
    lu = length(unique(charP))
    if (lu == 2) {
      method = "k"
      message("Setting method to Kendall")
    }
    else if (lu <= 5) {
      method = "s"
      message("Setting method to Spearman")
    }
    else {
      method = "p"
      message("Setting method to Pearson")
      if (is.null(winsorizeRER)) {
        message("Setting winsorizeRER=3")
        winsorizeRER = 3
      }
      if (is.null(winsorizetrait)) {
        message("Setting winsorizetrait=3")
        winsorizetrait = 3
      }
    }
  }
  win = function(x, w) {
    xs = sort(x[!is.na(x)], decreasing = T)
    xmax = xs[w]
    xmin = xs[length(xs) - w + 1]
    x[x > xmax] = xmax
    x[x < xmin] = xmin
    x
  }
  corout = matrix(nrow = nrow(RERmat), ncol = 3)
  rownames(corout) = rownames(RERmat)
  colnames(corout) = c("Rho", "N", "P")
  for (i in 1:nrow(corout)) {
    if (((nb <- sum(ii <- (!is.na(charP) & !is.na(RERmat[i, 
    ])))) >= min.sp)) {
      if (method != "p" && sum(charP[ii] != 0) < min.pos) {
        next
      }
      if (!weighted) {
        x = RERmat[i, ]
        indstouse = which(!is.na(x) & !is.na(charP))
        if (!is.null(winsorizeRER)) {
          x = win(x[indstouse], winsorizeRER)
        }
        else {
          x = x[indstouse]
        }
        if (!is.null(winsorizetrait)) {
          y = win(charP[indstouse], winsorizetrait)
        }
        else {
          y = charP[indstouse]
        }
        cres = cor.test(x, y, method = method, exact = F)
        corout[i, 1:3] = c(cres$estimate, nb, cres$p.value)
      }
      else {
        charPb = (charP[ii] > 0) + 1 - 1
        weights = charP[ii]
        weights[weights == 0] = 1
        cres = wtd.cor(RERmat[i, ii], charPb, weight = weights, 
                       mean1 = F)
        corout[i, 1:3] = c(cres[1], nb, cres[4])
      }
    }
    else {
    }
  }
  corout = as.data.frame(corout)
  corout$p.adj = p.adjust(corout$P, method = "BH")
  corout
}

##supplementary ERC functions
##updated GH 10/28/24


##make sure the gene list is contained in the erc matrix
clean_list = function(list, names){
  list <- unique(list)
  list <- list[which(match(list,names) != "NA")]
  list
}

pair_list <- function(list, erc_matrix , na.val = -2) {
  #provide vector list of genes (row and col names) and erc_matrix
  #returns vector list of gene x gene ERC values with no repeats or 'self x self' entries
  # good for computing averages and making histograms, etc...
  list <- clean_list(list , colnames(erc_matrix))
  mat <- erc_matrix[list,list]
  mat[which(mat == na.val)] <- NA
  mat[lower.tri(mat)]
}
##Make the ERC matrix symmetrical
make_symmetric = function(erc_matrix){
  size = dim(erc_matrix)[1]
  hold = matrix(NA,size,size)
  nindex = is.na(erc_matrix)   # find NA indices, make them zero
  hold[nindex] <- 0
  hold <- hold + t(erc_matrix) # add transpose to zeroes
  erc_matrix[nindex] <- hold[nindex] # combine matrices
  erc_matrix
}

##get permutations for a list of genes against the genome
permTestMat <- function(list, erc_matrix, perms=10000) {
  cleannames <- clean_list(colnames(erc_matrix),colnames(erc_matrix))
  list <- clean_list(list, cleannames)
  if(length(list) == 1){return(NaN)}
  obs <- mean(erc_matrix[list,list] , na.rm=TRUE)
  print(obs)
  if (is.nan(obs)) { return(NaN) }
  
  null <- c()
  while (length(null) < perms) {
    s <- sample(cleannames , length(list))
    m <- mean(pair_list(s , erc_matrix) , na.rm=TRUE)
    if (is.nan(m)) next
    null <- c(null , m)
  }
  pval <- sum(as.numeric(obs <= null)) / perms
  list( obs=obs , p=pval , null=null)
}
##Supply the output from computeERC to get the Fisher transformed correlation values
fisherTransform <- function(erc_matrix){
  out <- erc_matrix$cor
  out <- atanh(erc_matrix$cor)*sqrt((erc_matrix$count)-3)
  ##out[upper.tri(out)] <- erc_matrix[upper.tri(erc_matrix)]
  out[is.nan(out)] <- NA
  diag(out)<-0
  out
}

##Get a subset matrix based on two gene lists
betweencomplex <- function(list1 , list2 , erc_matrix) {
  list1 <- clean_list(list1 , colnames(erc_matrix))
  list2 <- clean_list(list2 , colnames(erc_matrix))
  mat <- erc_matrix[list1,list2]
  mat
}

##get mean ERC b/w 2 gene lists and perform permutation tests
#	pval1 is pval after permuting genes in list 1.  Similar for pval2.
#	In practice, take the mean or max of pval1,pval2. It's up to the user to decide which.
permTestPair <- function(list1, list2, erc_matrix, perms=10000) {
  cleannames <- clean_list(colnames(erc_matrix),colnames(erc_matrix))
  list1 <- clean_list(list1, cleannames)
  list2 <- clean_list(list2, cleannames)
  if(length(list1) == 1 | length(list2) == 1){return(NaN)}
  mat <- erc_matrix[list1,list2]
  obs <- mean(mat[upper.tri(mat)] , na.rm=TRUE)
  print(obs)
  if (is.nan(obs)) { return(NaN) }
  
  null1 <- c()
  null2 <- c()
  while (length(null1) < perms) {
    s1 <- sample(cleannames , length(list1))
    s2 <- sample(cleannames , length(list2))
    
    mat1 <- erc_matrix[s1,list2]
    mat2 <- erc_matrix[list1,s2]
    m1 <- mean(mat1[upper.tri(mat1)] , na.rm=TRUE)
    m2 <- mean(mat2[upper.tri(mat2)] , na.rm=TRUE)
    
    if (is.nan(m1) | is.nan(m2)) next
    null1 <- c(null1 , m1)
    null2 <- c(null2 , m2)
  }
  pval1 <- sum(as.numeric(obs <= null1)) / perms
  pval2 <- sum(as.numeric(obs <= null2)) / perms
  list( obs=obs , p=c(pval1,pval2) , null=c(null1,null2))
}

