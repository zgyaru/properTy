
#' @import SingleCellExperiment
#' @import dplyr
#' @import cowplot
#' @import ggplot2
#' @export
properTy = function(scMat,
                    clusterNames){
  clusterNames = data.frame(clusterNames[colnames(scMat),1],
                            row.names = colnames(scMat))
  if(nrow(clusterNames) != ncol(scMat)){
    stop("Please ensure that clusterNames dataset includes all cells in scMat.")
  }
  immune_path = system.file('immunePropertyGenes.txt',package = "properTy")
  geneSets = read.csv(immune_path,sep = '\t')
  immuneGenes = split(geneSets$Genes, geneSets$Property)
  rm(geneSets)
  ## aucell score
  ranks = AUCell_buildRankings(scMat,plotStats=FALSE)
  aucellRes = AUCell_calcAUC(immuneGenes, ranks)
  aucellScore = t(aucellRes)
  aucellScore = data.frame(aucellScore)
  aucellScore$cellCluster = as.character(clusterNames[,1])
  return(aucellScore)
}




#' @export
DensityProp = function(score,
                       celltypes){
  testPlot = score[which(score$cellCluster %in% celltypes),]
  g1 = ggplot(testPlot, aes(x=activation,fill=cellCluster))+
    geom_density(alpha=0.5)+theme_bw()+
    theme(panel.grid = element_blank())+
    ggtitle('Activation')+
    xlab("properTy Score")
  g2 = ggplot(testPlot, aes(x=cytotoxic,fill=cellCluster))+
    geom_density(alpha=0.5)+theme_bw()+
    theme(panel.grid = element_blank())+
    ggtitle('Cytotoxic')+
    xlab("properTy Score")
  g3 = ggplot(testPlot, aes(x=exhausted,fill=cellCluster))+
    geom_density(alpha=0.5)+theme_bw()+
    theme(panel.grid = element_blank())+
    ggtitle('Exhausted')+
    xlab("properTy Score")
  g4 = ggplot(testPlot, aes(x=inhibitory,fill=cellCluster))+
    geom_density(alpha=0.5)+theme_bw()+
    theme(panel.grid = element_blank())+
    ggtitle('Inhibitory')+
    xlab("properTy Score")
  g5 = ggplot(testPlot, aes(x=`co.simulatory`,fill=cellCluster))+
    geom_density(alpha=0.5)+theme_bw()+
    theme(panel.grid = element_blank())+
    ggtitle('Co-simulatory')+
    xlab("properTy Score")
  g6 = ggplot(testPlot, aes(x=`Na.ve.memory`,fill=cellCluster))+
    geom_density(alpha=0.5)+theme_bw()+
    theme(panel.grid = element_blank())+
    ggtitle('Naive/memory')+
    xlab("properTy Score")
  return(plot_grid(g6,g1,g5,g2,g3,g4, nrow=3))
}





#' @export
PieProp = function(score,
                   celltypes){
  cluster_property = group_by(score, cellCluster)
  cluster_property = summarise(cluster_property,
                               `Naive/memory` = mean(`Na.ve.memory`),
                               Activation = mean(activation),
                               `Co-simulatory` = mean(co.simulatory),
                               Cytotoxic = mean(cytotoxic),
                               Exhausted = mean(exhausted),
                               Inhibitory = mean(inhibitory)
  )
  cluster_property = as.data.frame(cluster_property)
  rownames(cluster_property) = cluster_property$cellCluster
  cluster_property = cluster_property[,which(colnames(cluster_property)!='cellCluster')]
  cluster_property = cluster_property/rowSums(cluster_property)

  par(mfrow=c(ceiling(length(celltypes)/2),2))
  for(c in celltypes){
    pie(as.numeric(cluster_property[`c`,]),
        colnames(cluster_property),main = `c`)
  }
}



#' @export
loadExample = function(){
  f = system.file("testData.rds",package = "properTy")
  if(f == ""){
    stop("Could not find example data directory, try-re-installing
         'testSctpa properTy'")
  }
  readRDS(f)
}

















