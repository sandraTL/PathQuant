
allMetabolitesAtDistanceX <- function(geneList){

    graph <- createGraphFromPathway("hsa01100")
    neighborsKeggIds <- c();
    absentyGene <- c();
    for(l in 1:length(geneList)){

     m1 <- getHeadTailKgmlIdOfEdge(graph@graph , geneList[l], graph@edgeDF);

    if(length(which(rowSums(is.na(m1))==2)) == 0){

       geneGraphIds <- c(as.vector(m1[,2]), as.vector(m1[,3]))
       neighbors <- c();

     for(i in 1:length(geneGraphIds)){
       tempV <-igraph::ego(graph@graph,i,
               nodes = geneGraphIds[i], "all", 0)
       neighbors <- c(neighbors, tempV)
     }

    neighbors <- unique(Reduce(c,neighbors))
    print(geneList[l])
    print(neighbors)

   # neighborsKeggIds <- c(neighborsKeggIds, geneList[l], neighbors,
   #                       as.vector(getNodeKeggIdbyId(neighbors,graph@graph)))

   # neighborsKeggIds <- c(neighborsKeggIds,
   # as.vector(getNodeKeggIdbyId(neighbors,graph@graph)))

   # print(as.vector(getNodeKeggIdbyId(neighbors,graph@graph)))



    }else{
        absentyGene <- c(absentyGene, geneList[l]);
    }

    }
 print(neighborsKeggIds)
 print(unique(Reduce(c,neighborsKeggIds)))
 #print(absentyGene)
 return <- neighborsKeggIds;
}

getNodeKeggIdbyId <- function(nodeList, graph){

    keggIdList <- c();

    for(j in 1:length(nodeList)){

        keggId <- igraph::V(graph)[nodeList[j]]$keggId

        keggIdList <- c(keggIdList, keggId)
    }

    return <- keggIdList;
}
