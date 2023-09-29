##Pipeline to compare the covariation of physically interacting domains versus non-interacting domains
##Updated by Jordan Little 9/28/2023

##package dependencies
library(PRROC)
library(phangorn)
library(plyr)
library(reshape2)
library(RERconverge)


prep_complexes = function(name, full_mat,complex_association, tp_list){
  ##make complex matrix symmetrical and only take the bottom half so we don't get
  ##duplicate edges
  complex_list = as.character(complex_association$domain[complex_association$complex == name])
  tp_list = tp_list[tp_list$GENEA %in% complex_list,]
  
  complex_mat = full_mat[complex_list, complex_list]
  complex = make_symmetric(complex_mat)
  complex[upper.tri(complex)] = NA
  
  raw_genes = unique(sapply(strsplit(complex,'_'),'[',1))
  ##set self x self cells to NA
  for(i in 1:length(raw_genes)){
    complex[grepl(paste(raw_genes[i],"_"), rownames(complex)), grepl(paste(raw_genes[i],"_"), colnames(complex))] = NA
  }
  
  
  domain_edges = reshape2::melt(complex, na.rm = TRUE)
  colnames(domain_edges) = c("GENEA","GENEB","ftERC")
  
  
  ##Turn the true positive list into a matrix so that we can smoosh it with the 
  ##erc edges 
  tp_mat = matrix(data = NA, nrow = length(complex_list), ncol = length(complex_list))
  rownames(tp_mat) = complex_list
  colnames(tp_mat) = complex_list
  
  for(i in 1:nrow(tp_list)){
    genea = tp_list[,1][i]
    geneb = tp_list[,2][i]
    tp_mat[genea, geneb] = 1
  }
  
  ##make the true positive matrix symmetrical and only take the bottom half so we
  ##don't have duplicate edges
  tp_mat[is.na(tp_mat)] = 0
  tp_mat[is.na(domain_mat)] = NA
  tp_mat = make_symmetric(tp_mat)
  tp_mat[upper.tri(tp_mat)] = NA
  tp_edges = reshape2::melt(tp_mat, na.rm = TRUE)
  colnames(tp_edges) = c("GENEA","GENEB","true_positive")
  
  ##merge the two edgelists
  domain_tp = merge(domain_edges, tp_edges)
  
  domain_tp
}

plot_roc_domains_only = function(edgelist, title){
  res=edgelist[order(edgelist$ftERC,decreasing = T),]
  
  roc_d=roc.curve(scores.class0 =res[res[,4]=="1",]$ERC,
                  scores.class1 =res[res[,4]=="0",]$ERC,
                  curve=T,)
  plot(roc_d, auc.main = TRUE,main = paste(title, "domain interactions"), color = 'black')
  roc_d
}

plot_all_complexes = function(ROC){
  plot.new()
  par(bg = "#a9a9a9")
  
  par(mar = c(5,4,4,11), xpd = TRUE) ##these dimensions have to be changed based on the plot window size in Rstudio
  all_roc = plot(ROC[[1]], main = "Complex ROC curves",xlab = "False positive rate", max.plot =TRUE, min.plot = TRUE, rand.plot = TRUE, fill.area = T, color =  colors[9], auc.main = FALSE)
  aucs = as.matrix(c(sapply(ROC,'[[',2)))
  aucs = format(round(aucs, 3))
  names = rownames(aucs)
  list_for_legend = mapply(function(x,y) paste(x,"(", y,")", sep = ""), names, aucs) 
  legend("topright", inset = c(-0.85,0.1), title = "Complex (ROC-AUC)",legend = list_for_legend, col = colors, cex=0.8, lty = 1, lwd = 3)
  for(i in 2:length(ROC)){
    plot(ROC[[i]], color = colors[i], add = TRUE)
    op = par(cex = 1.7)
  }

  
}

##get the permutation p-value and False positive rate for each complex
pairwise_rank_permutation = function(erc_mat, complex_edgelist, perms = 1000){
  complex_tp = complex_edgelist[,c("GENEA","GENEB")][complex_edgelist$true_positive ==1]
  complex_domains = unique(as.character(c(complex_edgelist$GENEA, complex_edgelist$GENEB)))
  complex_mat = make_symmetric(erc_mat[complex_domains, complex_domains])
  complex_pval = list()
  obs = NA
  for(true_pos in 1:nrow(complex_tp)){
    genea = sapply(strsplit(as.character(complex_tp$GENEA[true_pos]),'_'),'[',1)
    geneb = sapply(strsplit(as.character(complex_tp$GENEB[true_pos]),'_'),'[',1)
    if(length(complex_pval[[paste(genea,geneb)]]) == 3){
      obs = obs
      cnt = cnt + 1
    }else{
      obs = NA
      cnt = 1
    }
    rows = rownames(complex_mat)[grep(paste(genea,"_",sep=''), rownames(complex_mat))]
    cols = rownames(complex_mat)[grep(paste(geneb,"_",sep=''), rownames(complex_mat))]
    pair_mat = as.matrix(complex_mat[rows,cols])
    
    if(sum(rownames(pair_mat) == rows) != length(rows)){
      pair_mat = t(pair_mat)
      rownames(pair_mat) = rows
      colnames(pair_mat) = cols
    }

    pair_edgelist = reshape2::melt(pair_mat, na.rm = TRUE)
    if(nrow(pair_edgelist) <1){
      next
    }
    pair_edgelist$true_pos = 0
    pair_edgelist = pair_edgelist[order(pair_edgelist$value, decreasing = TRUE),]
    row.names(pair_edgelist) = NULL
    tp = which(pair_edgelist$Var1 == as.character(complex_tp$GENEA[true_pos]) & pair_edgelist$Var2 == as.character(complex_tp$GENEB[true_pos]))
    if(length(tp) ==0){
      tp = which(pair_edgelist$Var1 == as.character(complex_tp$GENEB[true_pos]) & pair_edgelist$Var2 == as.character(complex_tp$GENEA[true_pos]))
    }
    
    pair_edgelist$true_pos[tp] = 1
    temp = (1-tp)/(1-nrow(pair_edgelist))
    if(length(temp) == 0){
      next
    }
    null = c()
    m = c()
    while (length(null) < perms) {
      s <- sample(x=1:nrow(pair_edgelist) , cnt)
      for(k in 1:cnt){
        m[k] = (1-s[k])/(1-nrow(pair_edgelist))
      }
      m = mean(m)
      if (is.nan(m)) next
      null <- c(null , m)
    }
    
    obs = mean(c(temp, obs), na.rm = TRUE)
    pval <- sum(as.numeric(obs >= null)) / perms
    print(pval)
    complex_pval[[paste(genea,geneb)]] = list( obs=obs , p=pval , null=null)
  }
  complex_pval
}

complex_fpr_permutations = function(complex_fpr, perms = 1000){
  obs = mean(sapply(complex_fpr, "[[",1))
  null = rowMeans(as.matrix(sapply(complex_fpr, "[[",3)))
  pval = sum(as.numeric(obs>=null)) / perms
  complex_pval = list(obs=obs, p = pval, null=null)
  complex_pval
}

##generate domain gene trees with Phangorn
estimatePhangornTreeAll(alndir = "domains", treefile = "343yeast_master.tre", output.file = "domain_trees.tre")
##merge domain gene trees with full gene trees
domains = read.csv("domain_trees.tre", header = FALSE, sep = "\t")
domain_list = as.character(domains$V1)
full = read.csv("full_trees.tre", header = FALSE, sep = "\t")
all = rbind(domains, full)
write.table(all, "domain_trees.tre", col.names = FALSE, sep = "\t", quote = FALSE, row.names = FALSE)

##Run ERC on concatenated domain/full gene trees - make sure the correct file names are put in to run_ERC
##This part takes approximately 3.5 hrs

detach("package:RERconverge", unload = TRUE)
source("run_ERC.R")
ERC = readRDS("domain_by_domain_ERC.RDS")
##Fisher transform ERC matrix and get the bottom half of the matrix so that we don't get duplicate edges further down
ftERC = fishertransformed_updated(ERC)
ftERC = make_symmetric(ftERC)
ftERC[upper.tri(ftERC)] = NA

##Subset matrix on domains only
domain_fterc = ftERC[domain_list, domain_list]
##Move the domains into their complexes
complex_list = c("EIF3","MCM","NUP84","ORC","PAN1","SMC5_6","TREX","EXOCYST","COMA","SWI_SNF","CUL8_MMS1_MMS22_CTF4","GET4_GET5","ATG17_ATG31_ATG29","ESCRT_I","MITO_ATP","SEC23_24","EXPORTIN")
tp_list = read.csv("true_positives.csv", header = FALSE)
complex_association = read.csv("complex_association.csv")

complex_complete = vector(mode = "list", length = length(complex_list))
complex_complete = lapply(seq_len(length(complex_list)),function(X) domain_fterc )
names(complex_complete) = complex_list


##Generate complex edgelists with [GENEA, GENEB, ERC, true_positive]
complex_complete = lapply(complex_complete, prep_complexes, full_mat=complex_complete,name=name(complex_complete),complex_association=complex_association, tp_list=tp_list)



##Get the ROC curve and AUC for each complex, plot all curves on one plot
ROC = vector(mode = "list", length = complex_list)
names(ROC) = complex_list
ROC = lapply(complex_complete,plot_roc_domains_only, title = name(complex_complete))


##Get each protein pair's FPR 
individual_fpr = vector(mode = "list", length = length(complex_list))
names(individual_fpr) = complex_list
individual_fpr = lapply(complex_complete, pairwise_rank_permutation, erc_mat = domain_fterc, complex_edgelist = complex_complete)

##Get the complex average FPR and p-value
complex_fpr = vector(mode = "list", length = length(complex_list))
complex_fpr = lapply(individual_fpr, complex_fpr_permutations)
complex_fpr = unlist(complex_fpr)

