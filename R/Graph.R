
##############################   IGRAPH S4   ###################################
# object igraph instanciated as S4 classe to be able to use igraph objects in
# other S4 classes (Graphe)

setClass(Class = "igraph")

################################################################################
##########################   GRAPH ELEMENTS S4   ###############################
################################################################################

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

    # print("createIGraph")
    g <- igraph::graph.data.frame(object@edgeDF, directed=FALSE,
                                  vertices=object@nodeDF);
    return <- g;

})

################################################################################
##########################   GRAPHELEMENTS S4   ################################
################################################################################

# S4 classe Graph, contains all information important to create graphe argument
# and to calcul distances for it

setClass(
    Class= "Graph",
    representation = representation(
        graph = "igraph"
    ),
    contains=c("igraph", "GraphElements")
)


setGeneric("getPath", function(object, path.id.vec){
    standardGeneric("getPath")
}
)

setMethod("getPath", "Graph", function(object, path.id.vec) {

     print("getPath")

    path <- NULL

    for (i in 1:length(path.id.vec)) {

        metabolite.g.id <- igraph::V(object@graph)[path.id.vec[i]]$name
        metabolite.kegg.id <- which(object@nodeDF$kgmlId == metabolite.g.id)

        metabolite <- object@nodeDF$keggId[metabolite.kegg.id]


        if (length(path.id.vec) == 1) {

            path <- NA

        } else if (i+1 <= length(path.id.vec)) {

            # Create a separated function

            # Get the kgml id of the graph node id
            edge.g.id.1 <- igraph::V(object@graph)[path.id.vec[i]]$name
            edge.g.id.2 <- igraph::V(object@graph)[path.id.vec[i+1]]$name

            # Find the edges with the metabolites as products and substrates
            # Each metabolite node might be related to multiple edges
            # Find the edge that has both metabolite as substrate and porducts
            edge.kegg.id.1 <- which(object@edgeDF$substrateId == edge.g.id.1)
            edge.kegg.id.2 <- which(object@edgeDF$productId == edge.g.id.1)
            edge.kegg.id.3 <- which(object@edgeDF$substrateId == edge.g.id.2)
            edge.kegg.id.4 <- which(object@edgeDF$productId == edge.g.id.2)

            # Intersect subtrates and products
            edge.id.1 <- intersect(edge.kegg.id.1, edge.kegg.id.4)
            edge.id.2 <- intersect(edge.kegg.id.2, edge.kegg.id.3)

            # Choose the non-empty intersection
            if (length(edge.id.1) > 0) {

                edge.id <- edge.id.1

            } else {

                edge.id <- edge.id.2

            }

            # format the KEGG ids of the genes on the edge for print
            edge <- object@edgeDF$ko[edge.id]
            edge <- paste("[" , edge, "]", collapse = "")

            # Handle the begining of print for metabolite
            if (i == 1) {

                path <- paste(path, metabolite, sep = "")

            } else {

                path <- paste(path, metabolite, sep = " -> ")

            }
            # concat gene
            path <- paste(path, edge, sep = " -> ")

        } else {
            # concat last metabolite
            path <- paste(path, metabolite, sep = " -> ")

        }

    }

    if (is.data.frame(path) == TRUE) {
        if (nrow(path) == 0) {
            path <- NA;
        }
    }

    return <- path

})

################################################################################
######################## INPUT  DATA AND RESULT S4 #############################
################################################################################

setClass("data",
         representation = representation(
             id = "list",
             gene = "list",
             metabolite = "list"
         ), contains = "list")

setMethod("show","data",
          function(object){
              cat("data: ")
              print(object@id)
              print(object@gene)
              print(object@metabolite)
          })

setClass("data_annotated",
         representation = representation(
             data = "data",
             geneType = "list",
             geneBrite = "list",
             geneEnzyme = "list",
             metabolite_classification = "list"
         ), contains = c("list", "data"))


setMethod("show","data_annotated",
          function(object){
              cat("data_annotated: ")
              show(object@data)
              print(object@geneType)
              print(object@geneBrite)
              print(object@geneEnzyme)
              print(object@metabolite_classification)
          })

setGeneric("annotate_data", function(object){
    standardGeneric("annotate_data")
})


setMethod("annotate_data", "data" ,function(object){

    # print("annotate_data")

    data.df <- data.frame(gene = object@gene,
                          metabolite = object@metabolite)

    gene.annotation <- annotateAssociationData(data.df)

    metabolite.annotation <-
        getSuperClassByKEggId(as.vector(object@metabolite))
    # print(metabolite.annotation)
    data.annotated <- new("data_annotated",
                          data = object,
                          geneType = list(gene.annotation[,4]),
                          geneBrite = list(gene.annotation[,3]),
                          geneEnzyme = list(gene.annotation[,2]),
                          metabolite_classification = list(metabolite.annotation)
    )

    return <- data.annotated
})


setClass("data_result",
         representation = representation(
             data_annotated = "data_annotated",
             distance = "list",
             pathway = "list",
             path = "list"
         ), contains = "data_annotated")


setGeneric("create.data.result", function(object, distance){
    standardGeneric("create.data.result")
})

setMethod("create.data.result" ,"data_annotated", function(object, distance){

    # print("create.data.result")

    distance.vec <- NULL;
    pathway.vec <- NULL;
    path.vec <- NULL;

    for(i in 1:length(object@data@gene[[1]])) {

        match <- distance[distance$geneKEGGId == object@data@gene[[1]][i]  &
                              distance$metaboliteKEGGId == object@data@metabolite[[1]][i],]

        if (nrow(match) > 0) {

            # get weird results of all NA rows...
            distance.vec <- c(distance.vec, match[1,]$distance)
            pathway.vec <- c(pathway.vec, as.character(match[1,]$pathway))
            path.vec <- c(path.vec, as.character(match[1,]$path))

        } else {

            distance.vec <- c(distance.vec, NA)
            pathway.vec <- c(pathway.vec, NA)
            path.vec <- c(path.vec, NA)

        }

    }

    data.result.ob <- new("data_result",
                          data_annotated = object,
                          distance = list(distance.vec),
                          pathway = list(pathway.vec),
                          path = list(path.vec))


    return <- data.result.ob
})


setMethod("show","data_result",
          function(object){
              cat("data_result: ")
              show(object@data_annotated)
              print(object@distance)
              print(object@pathway)
              print(object@path)
          })


setGeneric("df.sub.data.result.numeric", function(object, numeric){

    standardGeneric("df.sub.data.result.numeric")

})

setMethod("df.sub.data.result.numeric" , "data_result",
          function(object, numeric = FALSE){

              # print("df.sub.data.result.numeric")

              df.out <- data.frame()

              id <- as.vector(object@data_annotated@data@id[[1]])
              gene <- as.vector(object@data_annotated@data@gene[[1]])
              metabolite <- as.vector(object@data_annotated@data@metabolite[[1]])
              geneType <- as.vector(object@data_annotated@geneType[[1]])
              geneEnzyme <- as.vector(object@data_annotated@geneEnzyme[[1]])
              geneBrite <- as.vector(object@data_annotated@geneBrite[[1]])
              metaboliteClassification <-
                  as.vector(object@data_annotated@metabolite_classification[[1]])
              distance = as.vector(object@distance[[1]])
              pathway  = as.vector(object@pathway[[1]])
              path = as.vector(object@path[[1]])

              df.out <- data.frame(
                  "id" = id,
                  "gene" = gene,
                  "metabolite" = metabolite,
                  "geneType" = geneType,
                  "geneEnzyme" = geneEnzyme,
                  "geneBrite" = geneBrite,
                  "metaboliteClassification" = metaboliteClassification,
                  "distance" = distance,
                  "pathway"  = pathway,
                  "path" = path
              )

              if (numeric == TRUE) {
                  df.out <- df.out[complete.cases(df.out[,8]),]
              }

              return <- df.out
          })

createGraphFromPathway <- function(pathwayId){

    # print("createGraphFromPathway")
    #test if KGML file was downloaded already

    if(isFileInDirectory(pathwayId) == FALSE) {

        getPathwayKGML(pathwayId);

    }
    # create df for vertices
    nodeDF <- getListNodeFromKGML(pathwayId);

    # create df edges
    edgeDF <- finalReactionEdgeDF(pathwayId);

    # test if nodeDF and edgeDF aren't empty
    if(!(is.null(nodeDF)) && !(is.null(edgeDF))) {

        # create graphEl objects
        graphEl <- new("GraphElements", nodeDF= nodeDF,
                       edgeDF= edgeDF, pathwayId = pathwayId);

        # create igraph with graphEl object elements
        igraphe <- createIGraph(graphEl);

        # create Graph object
        graphe <- new("Graph", graph = igraphe, graphEl);

    } else {

        graphe <- NULL;

    }

    return <- graphe;
}


setGeneric("associatedShortestPaths", function(object, data, path)
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

setMethod("associatedShortestPaths","Graph",
          function(object, data, path = FALSE){

              # print("associatedShortestPaths")

              # Intantiate empty data frames for results
              final.df <- data.frame()
              distances.df <- data.frame()
              paths.df <- data.frame()

              # Test if it is faster with cluster (not very long already)
              # cl <- setCluster()
              # doParallel::registerDoParallel(cl)

              # First apply get distance between each associations

              if (path == T) {
                  distances.df  <- get.distances(object, data)
                  tempDf <- as.data.frame(distances.df)

                  if(all(is.na(distances.df[,1])) == F){
                     paths.df <- get.paths(object, data)
                     final.df <- cbind("distance" = distances.df , "path" = paths.df)
                  } else {
                      distances.df  <- get.distances(object, data)
                      paths.df <- c(rep("not computed", length(distances.df)))
                      final.df <- cbind("distance" = distances.df, "path" = paths.df)
                  }
              } else {
                  distances.df  <- get.distances(object, data)
                  paths.df <- c(rep("not computed", length(distances.df)))
                  final.df <- cbind("distance" = distances.df, "path" = paths.df)
              }

              distances.df  <- get.distances(object, data)


              return <- final.df;

          })

setGeneric("get.distances", function(object, data)
{
    standardGeneric("get.distances")
}
)

setMethod("get.distances","Graph", function(object, data){

    distances.df <- data.frame()

    distances <-  apply(data,1, function(x) {

        dfTemp <- data.frame();
        # check if gene and metabolite have Graph Ids for compute

        if (is.na(x['metaboliteGraphId']) || is.na(x['geneGraphId'])) {
            df.temp <- NA;
        } else {

            op <- options(warn=2)
            tt <- tryCatch((igraph::shortest.paths(object@graph,
                                                   as.character(x[2]), as.character(x[4]))),
                           error=function(e) e,
                           warning=function(w) w)

            #catch warnings when ther is not pat between 2 nodes
            if (is(tt,"warning")) {}
            else if (is(tt,"error")) {}
            else {
                df.temp <- tt
            }
        }
        return <- as.numeric(df.temp);

    })

    # Form result in one data.frame row
    distances.df <- rbind(distances.df, distances)

    # Change to column
    distances.df <- t(distances.df)

    return <- distances.df

})

setGeneric("get.paths", function(object, data)
{
    standardGeneric("get.paths")
}
)

setMethod("get.paths","Graph", function(object, data){

    print("get.paths")
    paths.df <- data.frame()

    paths <- apply(data, 1, function(x) {

        df.temp <- data.frame();
        # check if gene and metabolite have Graph Ids for compute

        if (is.na(x['metaboliteGraphId']) || is.na(x['geneGraphId'])) {

            df.temp <- NA;

        } else {
            op <- options(warn=2)
            tt <- tryCatch((igraph::shortest_paths(object@graph,
                                                   as.character(x[2]), as.character(x[4]))),
                           error=function(e) e,
                           warning=function(w) w)

            #catch warnings when ther is not pat between 2 nodes

            if (is(tt,"warning")) { df.temp <- "no path exist"
            } else if (is(tt,"error")) { df.temp <- NA
            } else {
                df.temp <- tt$vpath[[1]];
                print(df.temp)
                # getPath <- keggIds of each node reported un vpath
                df.temp <- getPath(object, tt$vpath[[1]])
                print(df.temp)
            }
        }

        return <- df.temp;
    })

    # Form result in one data.frame row
    paths.df <- rbind(paths.df, paths)

    # Change to column
    paths.df <- t(paths.df)

    return <- paths.df

})

setGeneric("getDistanceAsso", function(object, pathwayId, path)
{
    standardGeneric("getDistanceAsso")
}
)

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
#'        association. First column are the genes and the sencond column are the
#'        metabolites. Only use KEGG Ids.
#' @param ordered [option] ascendent ordering of distance
#' @param commonNames get KEGG's Common Names of the KEGG Id in the results.
#' @keywords graph, shortestDistance, KEGG
#' @export
#' @examples getDistanceAsso("hsa01100",shinAndAlDF)

setMethod("getDistanceAsso", "data_annotated", function(object,
                                                        pathwayId,
                                                        path = FALSE){

    # print("GetDistanceAsso")

    # Arrange pathway name to get the KGML file
    pathwayId <- gsub(":", "", pathwayId)
    finalResult <- data.frame()
    # test input parameters
    # test_getDistanceAsso(pathwayId, association)

    # graph creation
    if (!exists("graphe")) {
        graphe <- createGraphFromPathway(pathwayId);
    }

    # test if the graph is null
    if (!is.null(graphe)) {

        # modify function calculate distance directly for association
        # result include, gene and metabolite KEGG Id
        # also the node or edge id, only the distance was used here
        # the node and or edges ids could be retrieved for other purposes
        # put condition to only compute on associations with enzyme

        if(path == FALSE){
            finalDF <- getFinalAssoDfSd(graphe, data.frame(
                gene = object@data@gene,
                metabolite = object@data@metabolite),
                path = FALSE)
        } else {
            finalDF <- getFinalAssoDfSd(graphe, data.frame(
                gene = object@data@gene,
                metabolite = object@data@metabolite),
                path = TRUE)
        }

        pathway.vec <- rep(pathwayId, nrow(finalDF))

        finalResult <- data.frame("geneKEGGId" = finalDF$geneKEGGId,
                                  "metaboliteKEGGId" = finalDF$metaboliteKEGGId,
                                  "distance" = finalDF$distance,
                                  "pathway" = pathway.vec,
                                  "path" = finalDF$path)
    } else {
        finalResult <- NULL
    }

    return <- finalResult;

})


#' isGeneInMap, metaboliteCommonName, metaboliteKEGGId,
#' isMetaboliteInMap, distance
#'
#' @param pairs dataframe with 2 columns, where each line reprensents a pair.
#'        The first column is genes and the sencond column is metabolites.
#'         Only use KEGG Ids.
#' @param pathway list of selected metabolic pathway in
#'        KEGG. Default = "All". Only use KEGG Ids.
#' @param ordered [option] ascendent ordering of distance
#' @param commonNames [option] get KEGG's Common Names of the KEGG Id in the
#'        results.
#' @param path [option] get reactional path used for srd computation.
#' @keywords graph, shortest path, shortest reactional path, KEGG
#' @export
#' @examples get.srd(pairs.df,"All")


get.srd <- function(pairs,
                    organism_code,
                    pathway = "All",
                    ordered = FALSE,
                    commonNames = FALSE,
                    path = FALSE){

    data.result.ob <- create.data.restul.ob(pairs, organism_code, pathway,
                                            ordered, commonNames, path)
    # object to df
    data.result.df <- df.sub.data.result.numeric(data.result.ob, F)

    if(commonNames == TRUE){
        # Get Common names could be tested with cluster
        geneCommonName <- as.vector(unlist(getCommonNames(as.vector
                                (unlist(data.result.df[,2])),"gene")))
        metaboliteCommonName <- as.vector(unlist(getCommonNames(as.vector
                                (unlist(data.result.df[,3])), "metabolite")))

        data.result.df <- cbind(data.result.df,
                          "geneName" = as.vector(geneCommonName),
                          "metaboliteName" = as.vector(metaboliteCommonName))
    }

    return <- data.result.df
}

create.data.restul.ob <- function(pairs,
                                  organism_code,
                                  pathway = "All",
                                  ordered = FALSE,
                                  commonNames = FALSE,
                                  path = FALSE){


    finalDF <- data.frame();
    pathway.list <- list();
    geneCommonName <- list();
    metaboliteCommonName <- list();

    # Create data object
    data.ob <- new("data",  id = list(c(1:nrow(pairs))),
                   gene = list(pairs[,1]),
                   metabolite = list(pairs[,2]))

    # Create data.annotated object for info on genes and metabolites
    data.annotated.ob <- annotate_data(data.ob)

    # Get List of all pathway maps in KEGG
    if (length(pathway) == 1 && pathway == "All") {
        pathway.list <- get.complete.list.pathways(organism_code)
        pathway.list <- pathway.list[-22] # problem with the encoding of this
    } else {
        pathway.list <- pathway;
    }

    # get distance for each pathway
    for(i in 1:length(pathway.list)){

        if (path == FALSE) {
            r <-getDistanceAsso(data.annotated.ob, pathway.list[i], FALSE);
        } else {
            r <-getDistanceAsso(data.annotated.ob, pathway.list[i], TRUE);
        }
        if (is.data.frame(r)) {

            # Combine the computed distance of each pathway
            if (nrow(r) >0) {
                finalDF <- rbind(finalDF, r)
            }

        }
    }

    # Keeps associations only when the gene and the metabolite are present in
    # the same map
    finalDF <- finalDF[!duplicated(finalDF),]

    # Choose the smallest distance for association between all pathways
    finalDF <- removeRowsDistanceAsso(finalDF)

    # Create final object data.results with annotation + distance and path
    data.result.ob <- create.data.result(data.annotated.ob, finalDF)


    return <- data.result.ob;


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

    # print("getDistanceAssoPerm")
    # test input parameters
    # test_getDistanceAssoPerm(pathwayId, data);

    finalDF <- data.frame();

    #graph creation
    if(!exists("graphe")){
        graphe <-  createGraphFromPathway(pathwayId);
    }

    #modify function calculate distance directly for association
    finalDF <- getFinalAssoDfSd(graphe, data, F);

    #Change Na in finalDF to Inf value
    finalDF <- changeDFassoToRigthDistances(finalDF);

    # order result by increasing distances
    finalDF$distance[is.na(finalDF$distance)] <- NaN;



    # Remove rows with distance between same gene and metbolites choosing
    # the smallest distance.
    finalDF <- removeRowsDistanceAsso(finalDF)


    # order found distances from small to big
    if(ordered == TRUE){
        finalDF <- finalDF[ order(finalDF[,7]), ]
    }

    # output for permutation use
    finalDF1 <- data.frame("geneKEGGId" = finalDF[,1],

                           "metaboliteKEGGId" = finalDF[,2],

                           "distance" = finalDF[,3]);
    # print(Sys.time())

    return <- finalDF1
}



changeDFassoToRigthDistances <- function(associatedShortestPathsDF){

    # print("changeDFassoToRigthDistances")

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
    #  print("GETIDGENEINGRAPH")

    f <- apply(data,1, function(x){


        hsaGene_Grep<- paste("\\",x[1],"\\b",sep="")
        listId <-grep(hsaGene_Grep, object@reactionDF$ko)


        # ' get both metabolites id from Graph related to the gene of data
        #m1 <- getHeadTailKgmlIdOfEdge(object@graph , x[1], object@edgeDF);

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

# Get graph ids for genes and metabolites
setMethod("fromAssosDFEntryToIGraphIdDF", "Graph", function(object, data,
                                                            indexMetabolite){
    #########################################################################
    #####  Add condition to insure indexMetabolite can only be 1 or 2   #####
    #########################################################################
    #
    #      testDataDF <- apply(data,1,function(x){
    #         if(substr(x[1],0,4) != "hsa:"){
    #
    #             stop("genes are not all valide KEGGId's",call. = FALSE )
    #         }
    #
    #
    #         if(substr(x[2],0,1) != "C" && length(x[2]) != 6){
    #
    #             stop("metabolites are not all valide KEGGId's",call. = FALSE )
    #         }
    #
    #     })
    colnames(data) <- c("gene", "metabolite")
    # print(data)


    # print("fromAssosDFEntryToIGraphIdDF")

    testDF<- list();
    completeAssoDF <- data.frame();
    f <- apply(data,1, function(x){

        #  print(x)
        # Get ids for the gene
        hsaGene_Grep<- paste("\\",x[1],"\\b",sep="")
        #  listId <-grep(x[1], object@edgeDF$ko)
        listId <-grep(hsaGene_Grep, object@edgeDF$ko)
        #print(object@edgeDF$ko)
        #print(listId)
        # Get ids for the metabolite
        hsaM_me<- paste("\\",x[2],"\\b",sep="")
        listId_M <-grep(hsaM_me, object@nodeDF$keggId)
        # print(listId_M)
        # Create a row for each substrate and product related to the gene
        # and for every node containing the metabolite
        # To calculate the distance between every option avaliable for
        # the gene and metabolite of an association
        if((length(listId_M)==0) || (length(listId)==0)){
            testDF2 <- list(as.character(x[1]),
                            NA,
                            as.character(x[2]),
                            NA)
            testDF <- append(testDF, list(testDF2))
        }else{
            for(i in 1:length(listId)){
                for(j in 1:length(listId_M)){

                    testDF1 <- list(x[1],
                                    as.numeric(as.character(object@edgeDF$substrateId[listId[i]])),
                                    as.character(x[2]),
                                    as.numeric(as.character(object@nodeDF$kgmlId[listId_M[j]])))

                    testDF2 <- list(gsub(" ","",as.character(x[1])),
                                    as.numeric(as.character(object@edgeDF$productId[listId[i]])),
                                    as.character(x[2]),
                                    as.numeric(as.character(object@nodeDF$kgmlId[listId_M[j]])))

                    testDF <- append(testDF, list(testDF1))
                    testDF <- append(testDF, list(testDF2))

                }

            }
        }

        # as.data.frame important to remove quotes from
        # entering in the string of the data.frame
        testDF <- as.data.frame(sapply(testDF, rbind))
        testDF <- as.data.frame(t(testDF))


        colnames(testDF) <- c("geneKEGGId", "geneGraphId",
                              "metaboliteKEGGId", "metaboliteGraphId")

        return<- testDF
    })

    f <- do.call(rbind, f)
    rownames(f) <- c(1:nrow(f))

    return <- f;


})


setGeneric("getFinalAssoDfSd", function(object, data, path)
{
    standardGeneric("getFinalAssoDfSd");
}
)

# get ids for every gene and metabolites and compute distances
setMethod("getFinalAssoDfSd", "Graph", function(object, data, path=FALSE){

    # print("getFinalAssoDfSd")

    finalDF <- data.frame();

    # get graph id for genes and metabolites of associations
    # get subsutrate and product of genes to compute distance
    m1g <- fromAssosDFEntryToIGraphIdDF(object, data, 1);

    # compute associations for every pair
    r3 <- associatedShortestPaths(object, m1g, path);

    # combine results with initial data.frame
    final <- cbind(m1g,r3)

    colnames(final) <- c("geneKEGGId", "geneGraphId",
                         "metaboliteKEGGId", "metaboliteGraphId","distance", "path")

    # Remove rows with distance between same gene and metbolites choosing
    # the smallest distance.

    final <- removeRowsDistanceAsso(final)

    return <- final;
})

# function merge 2 numeric vectors and return a vector with the smalest
# values for each position of the vectors
#

mergeVectorsLowerValues <- function(A,B) {

    # print("mergeVectorsLowerValues")

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



#' Function calculating shortest distance between every gene and every
#' metabolite.
#'
#' Function calculating shortest distance between every gene and every
#' metabolites on a graph model of KEGG map selected, where nodes are
#' metabolites and reactions (genes) are edges.
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

#' @examples getDistanceAll(metabolismOverviewMapKEGGId,completeGeneDF,
#'           completeMetaboDF)

getDistanceAll <- function(pathwayId, gene,
                           metabolite){

    # print("getDistanceAll")

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

    # print("getFinalDFSHortestDistance")
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

        stop("Sorry no metabolites of your entry data are on the selected map,
             thus no distance was calculated", call. = FALSE)
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

    # print("getIdMetabolitesInGraph")

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


#
# getAllShortestPath <- function(){
#
#
#
# print("alllo")
#     ##get list of all nodes implacated in an 1 edge at least
#     graphe <- createGraphFromPathway("hsa01100");
#     print("allo2")
#     p<- list();
#     for(i in 1:1002){
#         nodesVector <- as.vector(igraph::get.edges(graphe@graph, igraph::E(graphe@graph)[i]))
#
#         p<-c(p,nodesVector)
#     }
#     print("allo3")
#      p<- as.vector(unique(unlist(p)))
#
#
#    print("allo4")
# #     c <- list()
# #     sandp <- mapply(c, substrate, product, SIMPLIFY=FALSE)
# #     sandp <- as.data.frame(sandp)
# #     sandp <- sandp[!duplicated(sandp),]
# # print("allo3")
#     #exempel d'accès à un sommet
#
#  result <- list()
# # print(length(sandp))
# # for(i in 1:length(sandp)){
#     for(i in 200:899){
#       #  for(i in 900:950){
#         print(i)
#         print(igraph::V(graphe@graph)[i]$keggId);
#       for(j in 1:length(p)){
#
#         op <- options(warn=2)
#         tt <- tryCatch(igraph::shortest.paths(graphe@graph , as.character(igraph::V(graphe@graph)[i]$name),as.character(igraph::V(graphe@graph)[j]$name))
#
#                        ,error=function(e) e,
#                        warning=function(w) w)
#       # print(paste(igraph::V(graphe@graph)[i]$keggId, " ", igraph::V(graphe@graph)[j]$keggId, " ", tt))
#         #catch warnings when ther is not pat between 2 nodes
#        # print(tt)
#         if(is(tt,"warning")) {}
#         else if(is(tt,"error")) {}
#         else
#            result <- c(result, paste(as.character(igraph::V(graphe@graph)[i]$keggId), as.character(igraph::V(graphe@graph)[j]$keggId), tt))
#        # print(igraph::V(graphe@graph)[i]$keggId);
#   }
#     # print(result)
#  }
#
#   return <- result
# }


#
# setGeneric("allShortestPaths", function(object, data,
#                                         metabolite)
# {
#     standardGeneric("allShortestPaths")
# }
# )
#
# # function that uses the object Graph and calculs distances for every pair of
# # gene - metabolite from data.
#
# setMethod("allShortestPaths","Graph", function(object, data,
#                                                metabolite){
#
#     # removing duplicate in metaboliteKEGGId's
#     metabolite <-  unique(metabolite)
#
#     # name column names by metabolite ids
#     repeatedGeneVector <- paste(data[,1], sep="")
#     metaboliteVector <- paste(metabolite[,2], sep="")
#
#     # calcul all distances
#     pl <-  apply( data, 1, function(x){
#
#         op <- options(warn=2)
#         tt <- tryCatch(igraph::shortest.paths(object@graph , x[2],
#                                               as.vector(unlist(metabolite[,1])))
#                        ,error=function(e) e,
#                        warning=function(w) w)
#         #catch warnings when ther is not pat between 2 nodes
#         if(is(tt,"warning")) {}
#         else if(is(tt,"error")) {}
#         else
#             return <- tt;
#     })
#     pl <- data.frame(pl);
#
#     # combine all vectors of distances
#     output <- do.call(cbind.data.frame, pl);
#
#     ##### choose smallest values for each metabolites ######
#     # adding column with metaboliteKEGGId
#     output <- cbind(output, KEGGId = c(metaboliteVector));
#
#     output <- mergeRowsWithSmallestValueByKEGGId(output)
#
#     finalColNames <- output$KEGGId;
#
#     output <- output[-ncol(output)]
#
#     # transposing output to do the same with genes
#     output <- data.frame(t(output))
#
#     rownames(output) <- c(1:nrow(output))
#
#     ##### choose smallest values for each genes ######
#     # adding column with geneKEGGId
#
#     output <- cbind(output, KEGGId = c(repeatedGeneVector[1:nrow(output)]));
#
#     output <-mergeRowsWithSmallestValueByKEGGId(output)
#
#     # adding final rows and columns names (genes, metbaolites)
#     rownames(output) <- output$KEGGId;
#     output <- output[-ncol(output)]
#     colnames(output) <- finalColNames
#
#     return <- output;
# })




