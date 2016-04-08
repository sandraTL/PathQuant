
# function returning the nodes (sub, prod ) of the egdes representing the gene
#  of interest

getHeadTailKgmlIdOfEdge <- function(g, hsaGene,  reactionDF){

    # return lines of DF that contains hsa gene

    listId <-grep(hsaGene, reactionDF$ko)
    nodesVector1 <- data.frame();
    if(length(listId) > 1){
      f<- lapply(listId, function(x) {

        nodesVector <- as.vector(igraph::get.edges(g, igraph::E(g)[x[1]]));
        gene <- c(as.character(hsaGene))
        sub <- c(as.character(igraph::V(g)[nodesVector[1]]$name));
        prod <- c(as.character(igraph::V(g)[nodesVector[2]]$name));

        nodesVector1 <- data.frame(gene, sub, prod);

            return <- nodesVector1;
        })

         nodesVector1 <- f;
         nodesVector1 <- do.call(rbind.data.frame, nodesVector1)

    }else if(length(listId) == 1){

      nodesVector <- as.vector(igraph::get.edges(g, igraph::E(g)[listId[1]]));
      # from data frame id toi kgmlId(used in graph)
      gene <- c(as.character(hsaGene))
      sub <- c(as.character(igraph::V(g)[nodesVector[1]]$name));
      prod <- c(as.character(igraph::V(g)[nodesVector[2]]$name));


      nodesVector1 <- data.frame(gene, sub, prod);

    }else{

        gene <- c(as.character(hsaGene))
        sub <- c(NA)
        prod <- c(NA)
        nodesVector1 <-data.frame(gene, sub,prod)

    }
    return <- nodesVector1;

}

# function returning the kgmlId of the compound of interest

getCompoundNodeKgmlId <- function(g, compoundKeggId, nodeDF){


    listId <-grep(compoundKeggId, nodeDF$keggId)


    if(length(listId) > 1 ){

         f<- lapply(listId, function(x) {

             temp <- igraph::V(g)[x[1]];
             temp <- temp$name
             return <- temp;
         })

         temp<- f;
    }else if(length(listId) == 1){
        temp <- igraph::V(g)[listId[1]];
        temp <- temp$name;

    } else{
        temp <- NA;
    }

    return <- temp;

}


