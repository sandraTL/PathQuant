
##############################   IGRAPH S4   ###################################
# object igraph instanciated as S4 classe to be able to use igraph objects in
# other S4 classes (Graphe)

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

# Creation of object igraph with data from XML file of pathway
# of interest
setMethod("createIGraph", "GraphElements", function(object) {


    g <- igraph::graph.data.frame(object@edgeDF, directed=FALSE,
                                  vertices=object@nodeDF);
    return <- g;

})

##########################   GRAPHELEMENTS S4   ################################

# S4 classe Graph, contains all information important to create graphe argument
# and to calcul distances for it
setClass(
    Class= "Graph",
    representation = representation(
        graph = "igraph"
    ),
    contains=c("igraph", "GraphElements")
)

setGeneric("allShortestPaths", function(object, data,
                                        metabolite)
{
    standardGeneric("allShortestPaths")
}
)

# function that uses the object Graph and calculs distances for every pair of
# gene - metabolite from data.

setMethod("allShortestPaths","Graph", function(object, data,
                                               metabolite){

    # removing duplicate in metaboliteKEGGId's
    metabolite <-  unique(metabolite)

    # name column names by metabolite ids
    repeatedGeneVector <- paste(data[,1], sep="")
    metaboliteVector <- paste(metabolite[,2], sep="")

    # calcul all distances
    pl <-  apply( data, 1, function(x){

        op <- options(warn=2)
        tt <- tryCatch(igraph::shortest.paths(object@graph , x[2],
                                              as.vector(unlist(metabolite[,1])))
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

# function that uses the object Graph and calculs distances for every pair of
# gene - metabolite from data.
#
# Note that for 1 gene is attached 1 enzyme which is an edge in our graph.
# In order to calcul a distance from vertice to vertice this function is
# called twice in function 'getFinalDFSHortestDistance'.
#
# First call is to calcul all the distance for the one vertice of the gene
# to the all metabolites, the second call is to calcul the distance
# from the other vertice to all metabolites and then
# choose the smallest distance between the 2.
#
# The selection of the first vertice and the second vertice is done in a other
# function 'getIdGeneInGraph'
#

setMethod("associatedShortestPaths","Graph", function(object, data){

    #calcul al distances

    pl <-  apply(data,1, function(x){

        dfTemp <- data.frame();

        if(is.na(x['metaboliteGraphId']) || is.na(x['geneGraphId'])){

            dfTemp <- NA;
        }else{


            op <- options(warn=2)
            tt <- tryCatch((igraph::shortest.paths(object@graph,
                           x['geneGraphId'], x['metaboliteGraphId'])),
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

    # choosing smallest distance between the two metabolites of gene and
    # all metabolites
    outputFinal <- data.frame();



    outputFinal <- rbind(outputFinal, pl)
    outputFinal <- t(outputFinal)
    colnames(outputFinal) <- c("lengthShortestPath")

    return <- outputFinal;
})



#' Function calculating shortest distance between each gene-metabolite pairs.
#'
#' Function calculating shortest distance between each gene-metabolite
#' associations on your selected KEGG pathway.
#'
#'
#' If a gene or a metabolite is present on multiple edges or nodes, then the
#' shortest distance is selected.
#' Output : dataframe with the following columns : geneCommonName, geneKEGGId,
#' isGeneInMap, metaboliteCommonName, metaboliteKEGGId,
#' isMetaboliteInMap, distance
#'
#' @param pathwayId KEGG Id of selected pathway.
#' @param association Dataframe with 2 columns, where each line reprensents an
#'        associations. First column are the genes and the sencond column as the
#'        metabolites. Only use KEGG Ids.
#' @param ordered [option] ascendent ordering of distance
#' @keywords graph, shortestDistance, KEGG
#' @export
#' @examples getDistanceAsso("hsa01100",shinAndAlDF)

getDistanceAsso <- function(pathwayId, association, ordered = FALSE){

    pathwayId <- gsub("hsa:", "hsa", pathwayId)
    finalDF <- data.frame();

    # test input parameters
    test_getDistanceAsso(pathwayId,association)

    # graph creation
    if(!exists("graphe")){
    graphe <-  createGraphFromPathway(pathwayId);
    }

    # modify function calculate distance directly for association
    finalDF <- getFinalAssoDfSd(graphe, association);

    # Change Na in finalDF to Inf value
    finalDF <- changeDFassoToRigthDistances(finalDF);

    # order result by increasing distances
    finalDF$distance[is.na(finalDF$distance)] <- NaN;

    if(ordered == TRUE){
        finalDF <- finalDF[ order(finalDF[,7]), ]
    }

    #****** Could had information on the nodes or the multiple
    #****** path of shortest paths finalDF col 1 and 4

    # Remove rows with distance between same gene and metbolites choosing
    # the smallest distance.
    finalDF <- removeRowsDistanceAsso(finalDF)


    # Adding common names for genes and emtabolites
    geneCommonName <- as.vector(unlist(getCommonNames(as.vector
                                        (unlist(finalDF[,2])), "gene")))
    metaboliteCommonName <- as.vector(unlist(getCommonNames(as.vector
                                        (unlist(finalDF[,5])),"metabolite")))

     # create final output
    finalDF1 <- data.frame("geneCommonName" = geneCommonName,
                            "geneKEGGId" = finalDF[,2],
                            "isGeneInMap" = finalDF[,3],
                            "metaboliteCommonName" = metaboliteCommonName,
                            "metaboliteKEGGId" = finalDF[,5],
                            "isMetaboliteInMap" = finalDF[,6],
                            "distance" = finalDF[,7]);
}



# Function calculating shortest distance between each gene-metabolite pairs.
#
# Function calculating shortest distance between each gene-metabolite pairs.
# on a graph model of KEGG map selected, where nodes are metabolites and
# reactions are edges.
#
# If a gene or a metabolite is present on multiple edges or nodes, then
# shortest distance are calculated for every combinaison possible and the
# shortest distance is selected.

getDistanceAssoPerm <- function(pathwayId, data, ordered = FALSE){


    # test input parameters
    test_getDistanceAssoPerm(pathwayId, data);

    finalDF <- data.frame();

    #graph creation
    if(!exists("graphe")){
        graphe <-  createGraphFromPathway(pathwayId);
    }

    #modify function calculate distance directly for association
    finalDF <- getFinalAssoDfSd(graphe, data);

    #Change Na in finalDF to Inf value
    finalDF <- changeDFassoToRigthDistances(finalDF);

    # order result by increasing distances
    finalDF$distance[is.na(finalDF$distance)] <- NaN;

    if(ordered == TRUE){
        finalDF <- finalDF[ order(finalDF[,7]), ]
    }

    # Remove rows with distance between same gene and metbolites choosing
    # the smallest distance.
    finalDF <- removeRowsDistanceAsso(finalDF)

    # output for permutation use
    finalDF1 <- data.frame("geneKEGGId" = finalDF[,2],
                               "isGeneInMap" = finalDF[,3],
                               "metaboliteKEGGId" = finalDF[,5],
                               "isMetaboliteInMap" = finalDF[,6],
                               "distance" = finalDF[,7]);
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

    #test if KGML file was downloaded already
    if(isFileInDirectory(pathwayId) == FALSE){

        getPathwayKGML(pathwayId);
    }
    # create df for vertices
    nodeDF <- getListNodeFromKGML(pathwayId);

    # create df edges
    edgeDF <- finalReactionEdgeDF(pathwayId);

    # create graphEl objects
    graphEl <- new("GraphElements", nodeDF= nodeDF,
                   edgeDF= edgeDF, pathwayId = pathwayId);


    # create igraph with graphEl object elements
    igraphe <- createIGraph(graphEl);

    # create Graph object
    graphe <- new("Graph", graph = igraphe, graphEl);

    return <- graphe;
}



setGeneric("getIdGeneInGraph", function(object, data,
                                        indexMetabolite) {
    standardGeneric("getIdGeneInGraph");
}
)

# A igraph function
#
# Modification of dataframe giving in entry by user, from
# c("gene", "metabolite") to a dataFrame contaning one metabolite associated
# with the gene edge in the Graph by adding igraph id from Graph object
#
# indcieMetabolite selects either the first or second metabolite related to
# that gene
#
# if indexMetabolite = 1
# d1 = c("metabolite1giIdIniGraph", "metabolitesIdIniGraph")
# if indexMetabolite = 2
# d2 = c("metabolite2giIdIniGraph", "metabolitesIdIniGraph")

setMethod("getIdGeneInGraph", "Graph", function(object,
                                                data, indexMetabolite){
    #########################################################################
    ##### Add condition to insure indexMetabolite can only be 1 or 2    #####
    #########################################################################

    f <- apply(data,1, function(x){

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


setMethod("fromAssosDFEntryToIGraphIdDF", "Graph", function(object, data,
                                                            indexMetabolite){
    #########################################################################
    #####  Add condition to insure indexMetabolite can only be 1 or 2   #####
    #########################################################################
    testDataDF <- apply(data,1,function(x){
        if(substr(x[1],0,4) != "hsa:"){

            stop("genes are not all valide KEGGId's",call. = FALSE )
        }


        if(substr(x[2],0,1) != "C" && length(x[2]) != 6){

            stop("metabolites are not all valide KEGGId's",call. = FALSE )
        }

    })


    completeAssoDF <- data.frame();
    f <- apply(data,1, function(x){

        #  get both metabolites id from Graph related to the gene of data
        m1 <- getHeadTailKgmlIdOfEdge(object@graph , x[1], object@edgeDF);


        #  get metabolites id from Graph of dat
        m2 <- getCompoundNodeKgmlId(object@graph, x[2], object@nodeDF);

        temp <- x;
        colnames(m1) <- c("gene","geneGraphId", "geneGraphId")

        if(length(m2) > 1){
            f1 <- lapply(m2, function(x){
                tempMetabolite <- x;
                if(length(m1[,1]) > 1){
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

            if(nrow(m1) > 1){
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

# function merge 2 numeric vectors and return a vector with the smalest
# values for each position of the vectors
#

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



#' Function calculating shortest distance between every gene in a
#' gene-metabolite pairs and all metabolite.
#'
#' Function calculating shortest distance between every gene in a
#' gene-metabolite pairs of your association parameter and every metabolite
#' (in a pair or not) on a graph model of KEGG map selected, where nodes are
#' metabolites and reactions are edges.
#'
#' If a gene or a metabolite is present on multiple edges or nodes, then
#' shortest distance are calculated for every combinaison possible and the
#' shortest distance is selected.
#'
#' Output : dataframe with metabolite in columns and gene in rows with
#'          shortest distance values.
#'
#' @param pathwayId KEGG Id of selected pathway.
#' @param gene Dataframe of 1 column, representing all genes reported. Only use
#'        KEGG Ids.
#' @param metabolite Dataframe of 1 column, representing all the measured
#'        metabolites. Only use KEGG Ids.
#' @keywords graph, shortestDistance, KEGG
#' @export
#' @examples getDistanceAll(metabolismOverviewMapKEGGId,completeGeneDF,
#'           completeMetaboDF)

getDistanceAll <- function(pathwayId, gene,
                           metabolite){

    pathwayId <- gsub("hsa:", "hsa", pathwayId)

    # test inputs
    test_getDistanceAll(pathwayId, gene, metabolite)

    finalDF <- data.frame();

    # look when it was downloaded if it has been to long redownload
    if(isFileInDirectory(pathwayId) == FALSE){
        getPathwayKGML(pathwayId);
    }

    #graph creation
    graphe <-  createGraphFromPathway(pathwayId);

    finalDF <- getFinalDFSHortestDistance(graphe, gene,
                                                  metabolite);
    return <- finalDF;


}

setGeneric("getFinalDFSHortestDistance", function(object,  data,
                                                  metabolite)
{
    standardGeneric("getFinalDFSHortestDistance");
}
)

setMethod("getFinalDFSHortestDistance", "Graph", function(object,
                                data, metabolite){
    # print("getFinalDFShortestDistance")
    finalDF <- data.frame();

    # indexMetabolite = 1 to get the first metabolite attached to gene
    idM1g <- getIdGeneInGraph(object,data, 1);
    # indexMetabolite = 2 to get the second metabolite attached to gene
    idM2g <- getIdGeneInGraph(object,data, 2);

    # bind all node related to a gene to
    idMg <- data.frame(rbind(idM1g,idM2g))

    idMg <- unique(idMg)
    # print(idMg)

    idM <- getIdMetabolitesInGraph(object, metabolite)
    idM <- na.omit(idM)

    if(length(idMg) == 0){

        stop("Sorry no gene of your entry data are on the selected map, thus
              no distance was calculated", call. = FALSE)

    }
    if(length(idM) == 0){

        stop("Sorry no metabolites of your entry data are on the selected map, thus
              no distance was calculated", call. = FALSE)
    }

    # get all shortest paths for both ends of gene to all metabolites
    r <- allShortestPaths(object, idMg  , idM);

    return <- r;
})


setGeneric("getIdGeneInGraph", function(object, data,
                                        indexMetabolite) {
    standardGeneric("getIdGeneInGraph");
}
)


setMethod("getIdGeneInGraph", "Graph", function(object,
                                                data, indexMetabolite){
    #########################################################################
    ##### Add condition to insure indexMetabolite can only be 1 or 2    #####
    #########################################################################
    # print("getIdGeneInGraph")


    f <- apply(data,1, function(x){

        # ' get both metabolites id from Graph related to the gene of data
        m1 <- getHeadTailKgmlIdOfEdge(object@graph , x, object@edgeDF);

    })

    f <- do.call(rbind, f)

    if(!is.null(f)){

    f <- f[, c(1,indexMetabolite+1)]

    f <- f[rowSums(is.na(f)) != 1,] # delete rows with 2 NA or more
    colnames(f) <- c("geneKEGGId","geneGraphId ") #set colnames
    }
    return <- f;

})


setGeneric("getIdMetabolitesInGraph", function(object,completeMetaboDF) {
    standardGeneric("getIdMetabolitesInGraph");
}
)



# where Graph param is the Graph object
# where data is (gene, metabolite)
# where indice metabolites is only 1 or 2

setMethod("getIdMetabolitesInGraph", "Graph", function(object,
                                                       completeMetaboDF){

    f <- apply(completeMetaboDF,1, function(x){

        #  get both metabolites id from Graph related to the gene of data
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



