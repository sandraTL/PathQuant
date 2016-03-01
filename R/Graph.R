
##############################   IGRAPH S4   ###################################
#' object igraph instanciated as S4 classe to be able to use igraph objects in
#' other S4 classes (Graphe)
#'
setClass(Class = "igraph")

##########################   GRAPH ELEMENTS S4   ###############################

setClass("GraphElements" ,
         representation = representation(
             pathwayId = "character",
             nodeDF = "list",
             edgeDF = "list"
         ),  contains="list")

setGeneric("createIGraph", function(object)
{
    standardGeneric("createIGraph")
}
)

#' Creation of object igraph with data from XML file of pathway
#' of interest


setMethod("createIGraph", "GraphElements", function(object) {

    # print("createIGraph");
    g <- igraph::graph.data.frame(object@edgeDF, directed=FALSE,
                                  vertices=object@nodeDF);
    return <- g;

})

##########################   GRAPHELEMENTS S4   ################################

#'S4 classe Graph, contains all information important to create graphe argument
#'and to calcul distances for it
setClass(
    Class= "Graph",
    representation = representation(
        graph = "igraph"
    ),
    contains=c("igraph", "GraphElements")
)

setGeneric("allShortestPaths", function(object, associatedGeneMetaDF,
                                        completeMetaboliteDF)
{
    standardGeneric("allShortestPaths")
}
)

#' function that uses the object Graph and calculs distances for every pair of
#' gene - metabolite from data.
#'
#' Note that for 1 gene is attached 1 enzyme which is an edge in our graph.
#' In order to calcul a distance from vertice to vertice this function is
#' called twice in function 'getFinalDFSHortestDistance'.
#'
#' First call is to calcul all the distance for the one vertice of the gene
#' to the all metabolites, the second call is to calcul the distance
#' from the other vertice to all metabolites and then
#' choose the smallest distance between the 2.
#'
#' The selection of the first vertice and the second vertice is done in a other
#' function 'getIdGeneInGraph'
#'
#' param data where
#'  mg = id (from igraph in object Graph) of 1 metabolite related to the gene
#'  m = id (igraph in object Graph)
#'  genes = hsa:... id from KEGG
#'  metabolites = C.... id from KEGG
#'
#' @importFrom igraph shortest_paths
#' @param object Graph, data(mg, m, genes, metabolites)
#' @keywords  kegg
#' @examples allShortestPaths(Graph, data)


setMethod("allShortestPaths","Graph", function(object, associatedGeneMetaDF,
                                               completeMetaboliteDF){

    # removing duplicate in metaboliteKEGGId's
    completeMetaboliteDF <-  unique(completeMetaboliteDF)

    #'Name column names by metabolite ids
    repeatedGeneVector <- paste(associatedGeneMetaDF[,1], sep="")
    metaboliteVector <- paste(completeMetaboliteDF[,2], sep="")

    #'calcul all distances
    pl <-  apply( associatedGeneMetaDF, 1, function(x){

        op <- options(warn=2)
        tt <- tryCatch(igraph::shortest.paths(object@graph , x[2],
                                              as.vector(unlist(completeMetaboliteDF[,1])))
                       ,error=function(e) e,
                       warning=function(w) w)
        #catch warnings when ther is not pat between 2 nodes
        if(is(tt,"warning")) {}
        else if(is(tt,"error")) {}
        else
            return <- tt;
    })
    pl <- data.frame(pl);
    # combine all vectors of distances
    output <- do.call(cbind.data.frame, pl);

    ##### choose smallest values for each metabolites ######
    # adding column with metaboliteKEGGId
    output <- cbind(output, KEGGId = c(metaboliteVector));
    output <- mergeRowsWithSmallestValueByKEGGId(output)

    finalColNames <- output$KEGGId;
    output <- output[-ncol(output)]

    # transposing output to do the same with genes
    output <- data.frame(t(output))
    rownames(output) <- c(1:nrow(output))

    ##### choose smallest values for each genes ######
    # adding column with geneKEGGId
    output <- cbind(output, KEGGId = c(repeatedGeneVector));
    output <-mergeRowsWithSmallestValueByKEGGId(output)


    # adding final rows and columns names (genes, metbaolites)
    rownames(output) <- output$KEGGId;
    output <- output[-ncol(output)]
    colnames(output) <- finalColNames

    return <- output;
})

setGeneric("associatedShortestPaths", function(object, data)
{
    standardGeneric("associatedShortestPaths")
}
)

#' function that uses the object Graph and calculs distances for every pair of
#' gene - metabolite from data.
#'
#' Note that for 1 gene is attached 1 enzyme which is an edge in our graph.
#' In order to calcul a distance from vertice to vertice this function is
#' called twice in function 'getFinalDFSHortestDistance'.
#'
#' First call is to calcul all the distance for the one vertice of the gene
#' to the all metabolites, the second call is to calcul the distance
#' from the other vertice to all metabolites and then
#' choose the smallest distance between the 2.
#'
#' The selection of the first vertice and the second vertice is done in a other
#' function 'getIdGeneInGraph'
#'
#' param data where
#'  mg = id (from igraph in object Graph) of 1 metabolite related to the gene
#'  m = id (igraph in object Graph)
#'  genes = hsa:... id from KEGG
#'  metabolites = C.... id from KEGG
#'
#' @importFrom igraph shortest_paths
#' @param object Graph, data(mg, m, genes, metabolites)
#' @keywords  kegg
#' @examples associatedShortestPaths(Graph, data)

setMethod("associatedShortestPaths","Graph", function(object, data){

    #'calcul al distances

    pl <-  apply(data,1, function(x){

        dfTemp <- data.frame();

        if(is.na(x['metaboliteGraphId']) || is.na(x['geneGraphId'])){

            dfTemp <- NA;
        }else{


            op <- options(warn=2)
            tt <- tryCatch((igraph::shortest.paths(object@graph
                                                   ,x['geneGraphId'], x['metaboliteGraphId'])),
                           error=function(e) e,
                           warning=function(w) w)
            #catch warnings when ther is not pat between 2 nodes
            if(is(tt,"warning")) {}
            else if(is(tt,"error")) {}
            else{
                dfTemp <- tt;
            }

        }

        return <- dfTemp;
    })

    #' choosing smallest distance between the two metabolites of gene and
    #' all metabolites
    outputFinal <- data.frame();

#     pl1 <- lapply(pl, function(x){
#
#         if(length(x)>0){
#             lengthPath <- length(x$vpath[[1]])
#
#             output <- lengthPath;
#             #' if length of path is 0 -> couldnt reach a path
#             #' if length of path is >0 i have to do length -1
#             #' for example path i am lookinf for path from 1523 to 1523
#             #' path is 1523 and the length is 1. But the real length is 0.
#             if(lengthPath  != 0){
#                 lengthPath <-  (lengthPath -1);
#             }else if(lengthPath  == 0)
#                 lengthPath  <-  NA;
#
#             output <- lengthPath;
#         }else output <- NA
#         return <- output;
#
#     })

  #  print(pl)

    outputFinal <- rbind(outputFinal, pl)
    outputFinal <- t(outputFinal)
    colnames(outputFinal) <- c("lengthShortestPath")

    return <- outputFinal;
})



#' Fonction that calculates distance between each gene-metabolite pairs.
#'
#' The igrpah created to simulate the KEGG pathways as metabolties as nodes
#' and genes (related gene enzymes and reaction) as egdes.
#'
#' The shortest distance is taking from calculation from both vertices related
#' to a gene to the metabolites of interest.
#'
#' for param data:
#'      gene = KEGGid of gene hsa:...
#'      metabolites : KEGGid of metabolites C....
#' for param pathwayId : KEGG id of pathways without ':' ex: hsa01100
#'
#' @param data(gene, metabolites )
#' @keywords KEGG
#' @export
#' @examples getDistanceAsso(pathwayId, data, ouput)
getDistanceAsso <- function(pathwayId, data, ordered = FALSE,
                            output = c("xslx","data.frame")){

    #if the xml file was already dowmloaded
    if(isFileInDirectory(pathwayId) == FALSE){

        file <-  getPathwayKGML(pathwayId)
        # op <- options(warn=2)
        #        file <- tryCatch(getPathwayKGML(pathwayId),error=function(e) e,
        #                       warning=function(w) w)

        #        if(is(file,"warning")){
        #            if(file[1]$message == "download had nonzero exit status"){
        #            stop("pathway doesn't exist in KEGG database",call. = FALSE )
        #            }
        #        }

    }



    finalDF <- data.frame();

    #graph creation

    if(!exists("graphe")){

    graphe <-  createGraphFromPathway(pathwayId);
    }


    #modify function calculate distance directly for association
    finalDF <- getFinalAssoDfSd(graphe, data);


    #print("changeDFassosToRigthDistances")
     #Change Na in finalDF to Inf value
    finalDF <- changeDFassoToRigthDistances(finalDF);

    # order result by increasing distances
    finalDF$distance[is.na(finalDF$distance)] <- NaN;



    if(ordered == TRUE){
        finalDF <- finalDF[ order(finalDF[,7]), ]
    }

    ######################################################################
    ######################################################################
    ######## Could had information on the nodes or the multiple   ########
    ########               path of shortest paths                 ########
    ######################################################################
    ######################################################################

    # Remove rows with distance between same gene and metbolites choosing
    # the smallest distance.
    finalDF <- removeRowsDistanceAsso(finalDF)

    rowNumbers <- 1:length(finalDF[,1])
    row.names(finalDF) <- row.names(1:length(finalDF[,1]))
    finalDF <- subset(finalDF, , c(2,3,5,6,7))


    ##Normal use#########


#Adding common names for genes and emtabolites
#     geneCommonName <- getCommonNames(as.vector(unlist(finalDF[,1])), "gene")
#     geneCommonName <- as.vector(unlist(geneCommonName))
#
#     metaboliteCommonName <- getCommonNames(as.vector(unlist(finalDF[,3])),
#                                            "metabolite")
#     metaboliteCommonName <- as.vector(unlist(metaboliteCommonName))

#     finalDF1 <- data.frame("geneCommonName" = geneCommonName,
#                            "geneKEGGId" = finalDF[,1],
#                            "isGeneInMap" = finalDF[,2],
#                            "metaboliteCommonName" = metaboliteCommonName,
#                            "metaboliteKEGGId" = finalDF[,3],
#                            "isMetaboliteInMap" = finalDF[,4],
#                            "distance" = finalDF[,5]);

    ##permutation use
    finalDF1 <- data.frame("geneKEGGId" = finalDF[,1],
                           "isGeneInMap" = finalDF[,2],
                           "metaboliteKEGGId" = finalDF[,3],
                           "isMetaboliteInMap" = finalDF[,4],
                           "distance" = finalDF[,5]);  # return <- finalDF1;
    #     write.table(finalDF1, file = "AssociatedDataDistance.txt",sep="\t"
    #                 ,row.names=FALSE);
    if(output == "xslx"){assoDataXlsx(finalDF1)}
    else if(output == "data.frame"){return <- finalDF1;}

}

changeDFassoToRigthDistances <- function(associatedShortestPathsDF){

    for(row in 1:nrow(associatedShortestPathsDF)){

        # test if both gene and metabolite are in map
        # but no distance is found
        if(associatedShortestPathsDF[row,3] == TRUE &&
           associatedShortestPathsDF[row,6] == TRUE &&
           is.na(associatedShortestPathsDF[row,7])){

            associatedShortestPathsDF[row,7] <- Inf;
        }
    }
    return <- associatedShortestPathsDF;

}




createGraphFromPathway <- function(pathwayId){
    #print("CreateGraphFromPathway")
    #' create df for vertices
    nodeDF <- getListNodeFromKGML(pathwayId);

    #' create df edges
    edgeDF <- finalReactionEdgeDF(pathwayId);

    #' create graphEl objects
    graphEl <- new("GraphElements", nodeDF= nodeDF,
                   edgeDF= edgeDF, pathwayId = pathwayId);


    #' create igraph with graphEl object elements
    igraphe <- createIGraph(graphEl);

    #' create Graph object
    graphe <- new("Graph", graph = igraphe, graphEl);

    return <- graphe;
}



setGeneric("getIdGeneInGraph", function(object, associatedGeneMetaDF,
                                        indexMetabolite) {
    standardGeneric("getIdGeneInGraph");
}
)

#' A igraph function
#'
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to a dataFrame contaning one metabolite associated
#' with the gene edge in the Graph by adding igraph id from Graph object
#'
#' indcieMetabolite selects either the first or second metabolite related to
#' that gene
#'
#' if indexMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' if indexMetabolite = 2
#' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
#'
#' where Graph param is the Graph object
#' where data is (gene, metabolite)
#' where indice metabolites is only 1 or 2
#'
#' @param Graph, data, indexMetabolites
#' @keywords  igraph, node, kgmlId
#' @examples getIdGeneInGraph(g, data, indexMetabolite)

setMethod("getIdGeneInGraph", "Graph", function(object,
                                                associatedGeneMetaDF, indexMetabolite){
    #########################################################################
    ##### Add condition to insure indexMetabolite can only be 1 or 2    #####
    #########################################################################
    # print("getIdGeneInGraph")
    f <- apply(associatedGeneMetaDF,1, function(x){

        # ' get both metabolites id from Graph related to the gene of data
        m1 <- getHeadTailKgmlIdOfEdge(object@graph , x[1], object@edgeDF);

    })

    f <- do.call(rbind, f)


    f <- f[, c(indexMetabolite,3)]

    f <- f[rowSums(is.na(f)) != 1,] # delete rows with 2 NA or more
    colnames(f) <- c("mg","genes") #set colnames

    return <- f;

})

setGeneric("fromAssosDFEntryToIGraphIdDF", function(object, data,
                                                    indexMetabolite) {
    standardGeneric("fromAssosDFEntryToIGraphIdDF");
}
)

#' A igraph function
#'
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to a dataFrame contaning one metabolite associated
#' with the gene edge in the Graph by adding igraph id from Graph object
#'
#' indcieMetabolite selects either the first or second metabolite related to
#' that gene
#'
#' if indexMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' if indexMetabolite = 2
#' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
#'
#' where Graph param is the Graph object
#' where data is (gene, metabolite)
#' where indice metabolites is only 1 or 2
#'
#' @param Graph, data, indexMetabolites
#' @keywords  igraph, node, kgmlId
#' @examples fromAssosDFEntryToIGraphIdDF(g, data, indexMetabolite)

setMethod("fromAssosDFEntryToIGraphIdDF", "Graph", function(object, data,
                                                            indexMetabolite){
    #########################################################################
    #####  Add condition to insure indexMetabolite can only be 1 or 2   #####
    #########################################################################

    completeAssoDF <- data.frame();
    f <- apply(data,1, function(x){
        #  print("fromAssosDFEntryToIGraphIdDF---------1")
        # ' get both metabolites id from Graph related to the gene of data
        m1 <- getHeadTailKgmlIdOfEdge(object@graph , x[1], object@edgeDF);

        # print("fromAssosDFEntryToIGraphIdDF---------2")
        # ' get metabolites id from Graph of data
        m2 <- getCompoundNodeKgmlId(object@graph, x[2], object@nodeDF);

        temp <- x;
        colnames(m1) <- c("gene","geneGraphId", "geneGraphId")

        if(length(m2) > 1){
            #  print("fromAssosDFEntryToIGraphIdDF---------3")
            f1 <- lapply(m2, function(x){
                tempMetabolite <- x;
                if(length(m1[,1]) > 1){
                    #  print("fromAssosDFEntryToIGraphIdDF---------4")
                    f2 <- lapply(m1[,indexMetabolite+1], function(x){

                        geneGraphId <- c(as.character(x))
                        geneKEGGId <- c(temp[1])
                        metaboliteGraphId <- c(tempMetabolite)
                        metaboliteKEGGId <- c(temp[2])

                        completeAssoDF <- data.frame(geneGraphId, geneKEGGId,
                                         metaboliteGraphId, metaboliteKEGGId)

                        return<- completeAssoDF;
                    })

                    f2 <- do.call(rbind, f2)

                }else if(nrow(m1) == 1){
                    # print("fromAssosDFEntryToIGraphIdDF---------5")

                    geneGraphId <- c(as.character(m1[,indexMetabolite+1]))
                    geneKEGGId <- c(temp[1])
                    metaboliteGraphId <- c(tempMetabolite)
                    metaboliteKEGGId <- c(temp[2])

                    completeAssoDF <- data.frame(geneGraphId, geneKEGGId,
                                         metaboliteGraphId, metaboliteKEGGId)

                    return <-  completeAssoDF;
                }

            })
            f1 <- do.call(rbind, f1)
        }else {
            #  print("fromAssosDFEntryToIGraphIdDF---------7")

            if(nrow(m1) > 1){
                #print("fromAssosDFEntryToIGraphIdDF---------8")
                f2 <- lapply(m1[,indexMetabolite+1], function(x){

                    geneGraphId <- c(as.character(x))
                    geneKEGGId <- c(temp[1])
                    metaboliteGraphId <- c(m2)
                    metaboliteKEGGId <- c(temp[2])

                    completeAssoDF <- data.frame(geneGraphId, geneKEGGId,
                                         metaboliteGraphId, metaboliteKEGGId)

                    return<- completeAssoDF;
                })
                f2 <- do.call(rbind, f2)

            }else{

                #  print("fromAssosDFEntryToIGraphIdDF---------9")

                geneGraphId <- c(m1[indexMetabolite+1])
                geneKEGGId <- c(temp[1])
                metaboliteGraphId <- c(m2)
                metaboliteKEGGId <- c(temp[2])

                completeAssoDF <- data.frame(geneGraphId, geneKEGGId,
                                         metaboliteGraphId, metaboliteKEGGId)

                f2 <- completeAssoDF;
                return <- f2
            }

            f1 <- f2;


        }

        return<- f1


    })
    f <- do.call(rbind, f)
    f <- data.frame(f)

    return <- f;


})


setGeneric("getFinalAssoDfSd", function(object, data)
{
    standardGeneric("getFinalAssoDfSd");
}
)

#' A igraph function
#'
#'
#' Calls getIdGeneInGraph() twice to get :
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to 2 dataFrame
#' where indexMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' #' where indexMetabolite = 2
#' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
#'
#' This is because genes are edges and we calculate vertex to vertex
#' the distance between genes and metabolites.
#' From each gene we want to take the shortest distance between the 2 vertex
#' attached, so we calcul both and in an other function we take the smallest
#' and that will be the output.
#'
#' This function keeps, between the 2 metabolites of each gene vertexes,
#' the distance that is the smallest between the 2 to the same metabolite.
#'
#' @param graphe, kegg id of compound of interest, dataframe of compounds
#' (nodes)
#' @keywords  igraph, node, kgmlId
#' @examples getCompoundNodeKgmlId(g, data, indexMetabolite)

setMethod("getFinalAssoDfSd", "Graph", function(object, data){

    finalDF <- data.frame();


    # indexMetabolite = 1 to get the first metabolite attached to gene
    m1g <- fromAssosDFEntryToIGraphIdDF(object, data, 1);


    # indexMetabolite = 2 to get the second metabolite attached to gene
    m2g <- fromAssosDFEntryToIGraphIdDF(object, data, 2);

    m3g <- rbind(m1g,m2g)
    m3g <- m3g[!duplicated(m3g),]

    r3 <- associatedShortestPaths(object, m3g);

    final <- cbind(m3g,"distance" = r3)

    final <- removeRowsDistanceAsso(final)

    final[,6] <- !(is.na(final[,1]));
    final[,7] <- !(is.na(final[,3]));
    final <- final[,as.vector(c(1,2,6,3,4,7,5))]

    colnames(final) <- c("geneGraphId", "geneKEGGId", "geneInGraph",
                       "metaboliteGraphId", "metaboliteKEGGId",
                       "metaboliteInGraph" ,"distance")

    return <- final;
})

#'function merge 2 numeric vectors and return a vector with the smalest
#'values for each position of the vectors
#'
#' @param 2 numeric vectors
#' @keywords  kegg
#' @examples d1 <- c(1,2,1,2,3,4), d2 <- c(2,3,1,1,5,2)
#'  mergeVectorsLowerValues(d1,d2) returns -> 1,2,1,1,3,2

mergeVectorsLowerValues <- function(A,B) {

    dataCombineDF <- rbind(A,B);
    dataCombineDF <- t(dataCombineDF);

    apply(dataCombineDF,1,function(x) {
        result <- data.frame();
        if(is.na(x[1]) && is.na(x[2])){
            result <- cbind(NA);
        }
        else if (x[1] <= x[2]) {
            result <- cbind(x[1]);
        } else
            result <- cbind(x[2]);

        return <- result;
    })
}

#' Fonction that calculates every shortest distances between each gene and
#' all metabolites in KEGG pathway of choice
#'
#' The igrpah created to simulate the KEGG pathways as metabolties as nodes
#' and genes (related gene enzymes and reaction) as egdes.
#'
#' The shortest distance is taking from calculation from both vertices related
#' to a gene to the metabolites of interest.
#'
#' for param data:
#'      gene = KEGGid of gene hsa:...
#'      metabolites : KEGGid of metabolites C....
#' for param pathwayId : KEGG id of pathways without ':' ex: hsa01100
#'
#' @param data(gene, metabolites )
#' @keywords KEGG
#' @export
#' @examples getDistanceAll(pathwayId, data)
getDistanceAll <- function(pathwayId, associatedGeneMetaDF,
                           completeMetaboliteDF){


    finalDF <- data.frame();
    #if the xml file was already dowmloaded
    # look when it was downloaded if it has been to long redownload
    if(isFileInDirectory(pathwayId) == FALSE){
        getPathwayKGML(pathwayId);
    }

    if(!is.data.frame(associatedGeneMetaDF) || length(associatedGeneMetaDF[1,])< 2 ||
       length(associatedGeneMetaDF[1,])> 3){
        e <- simpleError("dataframe dimension is wrong, please enter you data
                         frame with KEGG id of genes (ex : hsa:00001) in first
                         column and associated KEGG id metabolites (ex: C00001)
                         in second column")
        tryCatch(stop(e), finally = print("please try again"))

    }else{
        if(is.null(getKGMLRootNode(pathwayId))){
            print("path you entered do not exist, enter valid hsa number without :");
        }else{
            #graph creation
            graphe <-  createGraphFromPathway(pathwayId);
            finalDF <- getFinalDFSHortestDistance(graphe, associatedGeneMetaDF,
                                                  completeMetaboliteDF );

            return <- finalDF;
        }
    }

}

setGeneric("getFinalDFSHortestDistance", function(object,  associatedGeneMetaDF,
                                                  completeMetaboliteDF)
{
    standardGeneric("getFinalDFSHortestDistance");
}
)

#' A igraph function
#'
#'
#' Calls getIdGeneInGraph() twice to get :
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to 2 dataFrame
#' where indexMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' #' where indexMetabolite = 2
#' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
#'
#' This is because genes are edges and we calculate vertex to vertex
#' the distance between genes and metabolites.
#' From each gene we want to take the shortest distance between the 2 vertex
#' attached, so we calcul both and in an other function we take the smallest
#' and that will be the output.
#'
#' This function keeps, between the 2 metabolites of each gene vertexes,
#' the distance that is the smallest between the 2 to the same metabolite.
#'
#' @param graphe, kegg id of compound of interest, dataframe of compounds
#' (nodes)
#' @keywords  igraph, node, kgmlId
#' @examples getCompoundNodeKgmlId(g, data, indexMetabolite)

setMethod("getFinalDFSHortestDistance", "Graph", function(object,
                                                          associatedGeneMetaDF, completeMetaboliteDF){
    # print("getFinalDFShortestDistance")
    finalDF <- data.frame();

    # indexMetabolite = 1 to get the first metabolite attached to gene
    idM1g <- getIdGeneInGraph(object,associatedGeneMetaDF, 1);
    # indexMetabolite = 2 to get the second metabolite attached to gene
    idM2g <- getIdGeneInGraph(object,associatedGeneMetaDF, 2);

    #bind all node related to a gene to
    idMg <- data.frame(rbind(idM1g,idM2g))

    idMg <- unique(idMg)
    # print(idMg)

    idM <- getIdMetabolitesInGraph(object, completeMetaboliteDF)
    idM <- na.omit(idM)

    #' get all shortest paths for both ends of gene to all metabolites
    r <- allShortestPaths(object, idMg  , idM);

    return <- r;
})


setGeneric("getIdGeneInGraph", function(object, associatedGeneMetaDF,
                                        indexMetabolite) {
    standardGeneric("getIdGeneInGraph");
}
)

#' A igraph function
#'
#' Modification of dataframe giving in entry by user, from
#' c("gene", "metabolite") to a dataFrame contaning one metabolite associated
#' with the gene edge in the Graph by adding igraph id from Graph object
#'
#' indcieMetabolite selects either the first or second metabolite related to
#' that gene
#'
#' if indexMetabolite = 1
#' d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
#' if indexMetabolite = 2
#' d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")
#'
#' where Graph param is the Graph object
#' where data is (gene, metabolite)
#' where indice metabolites is only 1 or 2
#'
#' @param Graph, data, indexMetabolites
#' @keywords  igraph, node, kgmlId
#' @examples getIdGeneInGraph(g, data, indexMetabolite)

setMethod("getIdGeneInGraph", "Graph", function(object,
                                                associatedGeneMetaDF, indexMetabolite){
    #########################################################################
    ##### Add condition to insure indexMetabolite can only be 1 or 2    #####
    #########################################################################
    # print("getIdGeneInGraph")

    # print(associatedGeneMetaDF)
    f <- apply(associatedGeneMetaDF,1, function(x){

        # ' get both metabolites id from Graph related to the gene of data
        m1 <- getHeadTailKgmlIdOfEdge(object@graph , x[1], object@edgeDF);

    })

    f <- do.call(rbind, f)


    f <- f[, c(1,indexMetabolite+1)]

    f <- f[rowSums(is.na(f)) != 1,] # delete rows with 2 NA or more
    colnames(f) <- c("geneKEGGId","geneGraphId ") #set colnames
    return <- f;

})


setGeneric("getIdMetabolitesInGraph", function(object,completeMetaboDF) {
    standardGeneric("getIdMetabolitesInGraph");
}
)

#' A igraph function
#'
#' where Graph param is the Graph object
#' where data is (gene, metabolite)
#' where indice metabolites is only 1 or 2
#'
#' @param Graph, data, indexMetabolites
#' @keywords  igraph, node, kgmlId
#' @examples getIdGeneInGraph(g, data, indexMetabolite)

setMethod("getIdMetabolitesInGraph", "Graph", function(object,
                                                       completeMetaboDF){

    #print("getIdMetabolitesInGraph")
    f <- apply(completeMetaboDF,1, function(x){

        # ' get both metabolites id from Graph related to the gene of data
        m1 <- getCompoundNodeKgmlId(object@graph, x[1], object@nodeDF);

        temp <- x;

        if(length(m1) > 1){

            f1 <- lapply(m1, function(x){

                mgmDF<- c(m=x, temp[1]);
                return<- mgmDF;

            })
        }else {

            f1<- c(m=  m1, temp[1]);

        }

        if(typeof(f1) == 'list'){
            f1 <- do.call(rbind, f1)

        }else{

            f1 <- data.frame(f1);
            f1 <- t(f1)
            f1 <- data.frame(f1);

        }

        return <- f1;
    })


    f <- do.call(rbind, f)

    f <- data.frame(f)

    return <- f;
})



