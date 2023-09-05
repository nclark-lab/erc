library(doParallel)
computeERC=function(rr, treesObj, minSp=NULL, doOnly=NULL, againstAll=F, parallel=F, clusterListOutput=NULL, saveFile=NULL, plot = F, win_value = 3){
  
  reportBin=treesObj$report[, TipLabels(treesObj$masterTree)]
  if(is.null(clusterListOutput)){
    pStr=apply(treesObj$report, 1, paste, collapse="")
    uClust=unique(pStr)
    char1=iconv("1", toRaw=T)[[1]]
    charL=(lapply(uClust, function(x){iconv(x, toRaw=T)[[1]]}))
    count1=unlist(lapply(charL, function(x){sum(x==char1)}))
    oo=order(-count1)
    uClust=uClust[oo]
    count1=count1[oo]
    idList=list()
    for(i in 1:nclust){
      idList[[i]]=which(pStr==uClust[i])
    }
  }
  else{
    count1=clusterListOutput$count1
    idList=clusterListOutput$idList
  }
  
  
  #set min species
  maxSp=ncol(treesObj$report)
  if(is.null(minSp)){
    minSp=round(maxSp*0.25)
    message(paste("minimum species is set to", minSp))
  }
  #initialize output
  n=treesObj$numTrees
  corout=matrix(nrow=n, ncol=n)
  colnames(corout)=rownames(corout)=names(treesObj$trees)
  countout=matrix(nrow=n, ncol=n)
  colnames(countout)=rownames(countout)=names(treesObj$trees)
  
  #fine the total number of operations
  nclust=sum(count1>=minSp)
  
  if(is.null(doOnly)){
    update=100
    total=(nclust*(nclust-1)/2+nclust)
    if(parallel){
      total=nclust
    }
    else{
      total=(nclust*(nclust-1)/2+nclust)
    }
    idListI=idList
    idListJ=idList
  }
  else{
    if(class(doOnly)=="character") {#there are genes and not indecied
      tmpdoOnly=double()
      #map indecies to clusterIDs
      map=double(nrow(treesObj$paths))
      for(i in 1:length(idList)){
        map[idList[[i]]]=i
      }

      for(gname in doOnly){
        gIndex=match(gname, names(treesObj$trees))
        if(is.na(gIndex)){
          message(paste("couldn't find", gname))
        }
       else{
         tmpdoOnly=c(tmpdoOnly, map[gIndex])
       }
      }
      doOnly=tmpdoOnly
    }
  
    update=500
    idListI=idList[doOnly]
    if(againstAll){
      idListJ=idList
    }
    else{
      idListJ=idListI
    }
    tmpn=length(doOnly)
    if(parallel){
      total=tmpn
    }
    else if (!againstAll){
      total=(tmpn*(tmpn-1)/2+tmpn)
    }
    else{
     total=length(idListI)*length(idListJ)
    }
  }
  if (parallel){
    update=1
  }
  pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                        total = total/update,
                        complete = "=",   # Completion bar character
                        incomplete = "-", # Incomplete bar character
                        current = ">",    # Current bar character
                        clear = F,    # If TRUE, clears the bar when finish
                        width = 100)      # Width of the progress bar
  
  
  
  
  mpaths = PathLengths(treesObj$masterTree)
  mpathsMat=as.matrix(mpaths[, c(2,1)])
  
  mympaths=mpaths
  mympaths$id=1:nrow(mpaths)
  mympaths=mympaths[nrow(mympaths):1,] #reverse
  done=0;
  for (i in 1:length(idListI)){
    if(parallel && ! pb$finished){
      pb$tick()
    }
    if(count1[i]<minSp){
      message("Remaining clusters too small")
      break
    }
    ii=idListI[[i]]
    
    # message(paste("cluster", i, "has length",length(ii)))
    
    rriall=rr[ii,,drop=F]
    tmpout=matrix(nrow=nrow(rriall), ncol=ncol(corout))
    if(!is.null(doOnly) && againstAll){
      istart=1
    }
    else{
      istart=i
    }
    if(!parallel){
      for (j in istart:length(idListJ)){
done=done+1
        if( j%% update ==0 && !pb$finished){
          pb$tick()
    
            }
        
        jj=idListJ[[j]]
        if(i==j && length(ii)==1){ #no point in doing self correlation
          next
        }
      #  bothIndex=which(colSums(treesObj$report[c(ii[1], jj[1]),])==2)
        bothIndex=reportBin[ii[1], ]&reportBin[jj[1],]
       # toKeep=colnames(treesObj$report)[bothIndex]
        if(sum(bothIndex)>=minSp){
       #   show(j)
         # mindex = KeptVerts(treesObj$masterTree, TipLabels(treesObj$masterTree) %in% toKeep)
        mindex=TreeTools:::kept_vertices(treesObj$masterTree$edge, bothIndex)[-1]>1L
     
        mindexI=which(mindex)
        
        tmp=mympaths[mympaths$start %in% mindexI,]
        tmpi=match(mindexI, tmp$end)
        pathsIndex=tmp$id[tmpi]
        pathsIndex=pathsIndex[!is.na(pathsIndex)]
        
        #previous indexing code can be used for debuggin
        #  pathsIndexThis = KeptPaths(mpaths, mindex, all = FALSE)
        #  pathsThis= mpathsMat[pathsIndexThis, ]
         # pathsIndex=  treesObj$matIndex[as.matrix(pathsThis[, ])]
        #  stopifnot(length(setdiff(pathsIndex, pathsIndex2))==0)
          
          #show(mean(!is.na(rr[c(ii,jj), pathsIndex])))
          rri=rriall[,pathsIndex,drop=F]
          rrj=rr[jj,pathsIndex,drop=F]
          tmp=cor(win(t(rri),win_value), win(t(rrj),win_value))
          if(plot == T){
            plot(win(t(rri),win_value), win(t(rrj),win_value), xlab = rownames(rri), ylab = rownames(rrj))
          }
          ##tmp=cor(t(rri), t(rrj))
          corout[ii,jj]=tmp
          countout[ii,jj]=length(pathsIndex)
        }
      }
    } #end not parallel
    else{
      outlist=foreach (j=istart:length(idListJ)) %dopar%{
        jj=idListJ[[j]]
#        bothIndex=which(colSums(treesObj$report[c(ii[1], jj[1]),])==2)
 #       toKeep=colnames(treesObj$report)[bothIndex]
        bothIndex=reportBin[ii[1], ]&reportBin[jj[1],]
        if(sum(bothIndex)>=minSp){
        
          
          mindex=TreeTools:::kept_vertices(treesObj$masterTree$edge, bothIndex)[-1]>1L
          
          mindexI=which(mindex)
          
          tmp=mympaths[mympaths$start %in% mindexI,]
          tmpi=match(mindexI, tmp$end)
          pathsIndex=tmp$id[tmpi]
          pathsIndex=pathsIndex[!is.na(pathsIndex)]
      
          #show(mean(!is.na(rr[c(ii,jj), pathsIndex])))
          rri=rriall[,pathsIndex,drop=F]
          rrj=rr[jj,pathsIndex,drop=F]
          ##JL added winsorization 7/21/22
          tmp=cor(win(t(rri),win_value), win(t(rrj),win_value))
          if(plot == T){
            plot(win(t(rri),win_value), win(t(rrj),win_value))
          }
          ##tmp=cor(t(rri), t(rrj))
          return(list(cors=tmp, lengths=length(pathsIndex), jj=jj))
        }
      }
      for (k in 1:length(outlist)){
        jj=outlist[[k]]$jj
        corout[ii,jj]=outlist[[k]]$cors
        countout[ii,jj]=outlist[[k]]$lengths
      }
    }
    
  }
  message("Done!")
  result=list(cor=corout, count=countout)
if(!is.null(saveFile)){
  saveRDS(result, saveFile)
}
return(result)
}

getClusterList=function(treesObj){
  pStr=apply(treesObj$report, 1, paste, collapse="")
  uClust=unique(pStr)
  char1=iconv("1", toRaw=T)[[1]]
  charL=(lapply(uClust, function(x){iconv(x, toRaw=T)[[1]]}))
  count1=unlist(lapply(charL, function(x){sum(x==char1)}))
  oo=order(-count1)
  uClust=uClust[oo]
  count1=count1[oo]
  nclust=length(uClust)
  idList=list()
  for(i in 1:nclust){
    idList[[i]]=which(pStr==uClust[i])
  }
  
  return(list(idList=idList, count1=count1))
}

findYeastGenes=function(inGenes){
  names(which(unlist(lapply(og2geneList, function(x){length(intersect(x, inGenes))>0}))))
}
