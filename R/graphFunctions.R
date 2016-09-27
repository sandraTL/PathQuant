
# function returning the nodes (sub, prod ) of the egdes representing the gene
#  of interest

getHeadTailKgmlIdOfEdge <- function(g, hsaGene,  reactionDF){

    # return lines of DF that contains hsa gene
    hsaGene_Grep<- paste("\\",hsaGene,"\\b",sep="")
    listId <-grep(hsaGene_Grep, reactionDF$ko)


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

    compoundKeggId1 <- paste("\\",compoundKeggId,"\\b",sep="")

    listId <-grep(compoundKeggId1, nodeDF$keggId)


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

grep2DF <- function(df1, df2){
   j <- 0 ;
    for(i in 1:length(df1[,1])){

    item <- paste("\\",df1[i,1],"\\b",sep="")

    listId <-grep(item, df2)
    if(length(listId) > 0 ){
        j <- j+1;
    match <- paste(j,df1[i,1], listId, sep = " ");

    }

    }

}

