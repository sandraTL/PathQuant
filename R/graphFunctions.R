#' A igraph function
#'
#' function returning the nodes (sub, prod ) of the egdes representing the gene
#'  of interest
#'
#' @param graphe, gene of interest, dataframe of reaction (edges) in the graph
#' @keywords  graph, edge
#' @examples getHeadTailKgmlIdOfEdge(g, hsa, reactionDF)

getHeadTailKgmlIdOfEdge <- function(g, hsaGene,  reactionDF){

    # return lines of DF that contains hsa gene

    listId <-grep(hsaGene, reactionDF$ko)


    nodesVector1 <- data.frame();
    if(length(listId) > 1){
       # print("getHeadTailKgmlIdOfEdge---------1")
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

       # print("getHeadTailKgmlIdOfEdge---------2")
      nodesVector <- as.vector(igraph::get.edges(g, igraph::E(g)[listId[1]]));
      # from data frame id toi kgmlId(used in graph)
      gene <- c(as.character(hsaGene))
      sub <- c(as.character(igraph::V(g)[nodesVector[1]]$name));
      prod <- c(as.character(igraph::V(g)[nodesVector[2]]$name));


      nodesVector1 <- data.frame(gene, sub, prod);

    }else{

      #  print("getHeadTailKgmlIdOfEdge---------3")
        gene <- c(as.character(hsaGene))
        sub <- c(NA)
        prod <- c(NA)
        nodesVector1 <-data.frame(gene, sub,prod)

    }
    #print(nodesVector1)
    return <- nodesVector1;

}


#' A igraph function
#'
#' function returning the kgmlId of the compound of interest
#' @param graphe, kegg id of compound of interest, dataframe of compounds
#' (nodes)
#' @keywords  graph, node, kgmlId
#' @examples getCompoundNodeKgmlId(g, compoundKeggId, nodeDF)


getCompoundNodeKgmlId <- function(g, compoundKeggId, nodeDF){


    listId <-grep(compoundKeggId, nodeDF$keggId)

    #########################################################################
    #########################################################################
    ###                   PUT THIS IN Graph class                        ####
    ###           ADD CONDITION FOR WHEN IS IN MULTIPLE NODES            ####
    #########################################################################
    #########################################################################

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
#      AM I USING THIS
# #' A igraph function
# #'
# #' fucntion return distance between 2 compound of the graph, in our case
# #' the first (c1) compound is one of the edge of our gene of interest and the second
# #' (c2) is the target metabolite
# #' @param graphe, kegg id of compound of interest, dataframe of compounds
# #' (nodes)
# #' @keywords  graph, node, kgmlId
# #' @examples getCompoundNodeKgmlId(g, compoundKeggId, nodeDF)
#
# getShortestDistanceBetween2Compounds <- function(g, c1, c2){
#
#
#     ## fonctionne seulement lorsqu'il y a une seul distance de calculé
#     # retourne la première distance. Si on veut calculer plusieurs distance
#     # il faudra ajusté...
#     ###                   PUT THIS IN Graph class                        ####
#     distance <- igraph::shortest.paths(g,
#         igraph::V(g)[as.character(c1)], igraph::V(g)[as.character(c2)]);
#
#     return <- distance[1];
#
#
#
# }

