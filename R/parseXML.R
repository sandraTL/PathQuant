
# This function extract from KGML a list of metabolites (id and name) that
# makes the nodes of the graph ** some metablite can from more than 1 node,
# which is why the id of the nodes is not the name of the metabolite.
# data.frame (id, keggId)

getListNodeFromKGML <- function(pathwayId) {

    xmltop <- getKGMLRootNode(pathwayId);


    nodeListId <- XML::xpathSApply(xmltop, "//entry[@type = 'compound']",
                                   function(x) (XML::xmlAttrs(x))['id']);
    nodeListCpd <- XML::xpathSApply(xmltop, "//entry[@type = 'compound']",
                                   function(x) (XML::xmlAttrs(x))['name']);
    nodeListColor <- XML::xpathSApply(xmltop, "//entry[@type = 'compound']//graphics",
                                    function(x) (XML::xmlAttrs(x))['fgcolor'])
    nodeDF <- data.frame(kgmlId = as.vector(unlist(nodeListId)),
                         keggId = as.vector(unlist(nodeListCpd)),
                         color = as.vector(unlist(nodeListColor)));

    nodeDF <- correctKeggIdString(nodeDF);

    return <- nodeDF;

}


# This function creates a data.frame containing the first part of our edge list
# the "gene" reactions
# data.frame(kgmlIdSubstratre, kgmlIdProduct, substrateName, productName,
# kgmlIdReaction, reactionName, reactionType)

getListReactionFromKGML <- function(pathwayId) {

    #gives content of root
    xmltop <- getKGMLRootNode(pathwayId);


    reactionIdNodes <- XML::getNodeSet(xmltop, "//reaction");

    reactionId <- lapply(reactionIdNodes,
                         function(x) XML::xmlAttrs(x)['id']);
    reactionName <- lapply(reactionIdNodes,
                           function(x) XML::xmlAttrs(x)['name']);
    reactionType <- lapply(reactionIdNodes,
                           function(x) XML::xmlAttrs(x)['type']);
    substrateId <- lapply(reactionIdNodes, XML::xpathApply,path = './substrate',
                          function(x) XML::xmlAttrs(x)['id']);
    substrateName <- lapply(reactionIdNodes, XML::xpathApply,
                            path = './substrate',
                            function(x) XML::xmlAttrs(x)['name']);
    productId <- lapply(reactionIdNodes, XML::xpathApply,
                        path = './product',
                        function(x) XML::xmlAttrs(x)['id']);
    productName <- lapply(reactionIdNodes, XML::xpathApply,
                          path = './product',
                          function(x) XML::xmlAttrs(x)['name']);
    reactionList <- do.call(rbind.data.frame,
                            mapply(cbind,
                                   "substrateId" = substrateId,
                                   "productId" = productId,
                                   "substrateName" = substrateName,
                                   "productName" = productName,
                                   "reactionId" = reactionId,
                                   "reactionType" = reactionType,
                                   "reactionName" = reactionName));

    return <- reactionList;

}


# This function extract from KGML a list of reaction of entry-type "gene"
# Returns a dataFrame object : id, entryId, reaction, ko

getListEdgeFromGeneKGML <- function(pathwayId) {

    # get the root of the KGML document
    xmltop <- getKGMLRootNode(pathwayId);

    # Get value of atributes (id, name(ko) and reaction) of entry of type
    # ortholog, some don't have reaction argument thus won't have an edge in
    # our final graph they will have NA in data.frame.
    edgeListId <- XML::xpathSApply(xmltop, "//entry[@type = 'gene']",
                                   function(x) (XML::xmlAttrs(x))['id']);
    edgeListKo <- XML::xpathSApply(xmltop, "//entry[@type = 'gene']",
                                   function(x) (XML::xmlAttrs(x))['name']);
    edgeListReaction <- XML::xpathSApply(xmltop, "//entry[@type = 'gene']",
                                         function(x) (XML::xmlAttrs(x))['reaction']);


    edgeDF <- data.frame("reactionId" = as.vector(as.character(edgeListId)),
                         "reactions" = as.vector(as.character(edgeListReaction)),
                         "type" = as.vector("gene"),
                         "ko" = as.vector(as.character(edgeListKo)));
    return <- edgeDF;

}

# This function extract from KGML a list of reaction of entry-type "gene"
# Returns a dataFrame object : id, entryId, reaction, ko

getListOrthologGeneFromKGML <- function(pathwayId) {

    # get the root of the KGML document
    xmltop <- getKGMLRootNode(pathwayId);

    # Get value of atributes (id, name(ko) and reaction) of entry of type
    # ortholog, some don't have reaction argument thus won't have an edge in
    # our final graph they will have NA in data.frame.
    orthologListId <- XML::xpathSApply(xmltop, "//entry[@type = 'ortholog']",
                                   function(x) (XML::xmlAttrs(x))['id']);
    print(length(orthologListId))
    orthologListKo <- XML::xpathSApply(xmltop, "//entry[@type = 'ortholog']",
                                   function(x) (XML::xmlAttrs(x))['name']);
    print(length(orthologListKo))
    orthologListReaction <- XML::xpathSApply(xmltop, "//entry[@type = 'ortholog']",
                                         function(x) (XML::xmlAttrs(x))['reaction']);
    print(length(orthologListReaction))

    print(length(orthologListCoords))
    orthologDF <- data.frame("reactionId" = as.vector(as.character(orthologListId)),
                             "reactions" = as.vector(as.character(orthologListReaction)),
                             "ko" = as.vector(as.character(orthologListKo)),
                             "x" = as.vector(rep(-1, length(orthologListReaction))),
                             "y" = as.vector(rep(-1, length(orthologListReaction))));

    for(row in 1:length(orthologDF[,1])){
      id <- orthologDF[row,1]
      coords <- XML::xpathSApply(xmltop, "//entry[@id = id]//graphics",
                         function(x) (XML::xmlAttrs(x))['coords']);
      print(length(coords))
     if(length(coords) == 1){
          orthologDF[row,4] <- coords;
     }
    }

    orthologListNameCoords <- XML::xpathSApply(xmltop, "//entry[@type = 'ortholog']//graphics",
                                           function(x) (XML::xmlAttrs(x))['name']);
    print(length(orthologListNameCoords))
    orthologListCoords <- XML::xpathSApply(xmltop, "//entry[@type = 'ortholog']//graphics",
                                         function(x) (XML::xmlAttrs(x))['coords']);


    print(length(orthologDF[,1]))
    orthologCoords <- data.frame("nameCoords" =as.vector(as.character(orthologListNameCoords)),
                                 "coords" = as.vector(as.character(orthologListCoords)));
    print(length(orthologCoords[,1]))
    #return <- orthologDF;

}



getKGMLRootNode <- function(pathwayId){
    # get the root of the KGML document
    pathFile <- toStringPathFile(pathwayId);

    if(is.na(file.info(pathFile)$size)== FALSE){

        xmlfile <- XML::xmlParse(pathFile);
        xmltop <- XML::xmlRoot(xmlfile); # gives content of root
    }else{
        getPathwayKGML(pathwayId)
        xmlfile <- XML::xmlParse(pathFile);
        xmltop <- XML::xmlRoot(xmlfile); # gives content of root


    }

    return <- xmltop;

}



toStringPathFile <- function(pathwayId){

    # concatenation of pathwayId to set swdir for the xml

    s2 <-  toString(pathwayId);
    s3 <- ".txt"

    pathFile <- paste(s2, s3, sep="");

    return <- pathFile;
}

getCommonNames <- function(vectorOfKEGGIds, type = c("gene","metabolite")){


        count <- 1;

    ### Vérifiez la connection internet

    if(length(vectorOfKEGGIds) > 10 ){
        count <- 0;
    # analysis of submaps - beggining
        names <- character();

        while(count <= length(vectorOfKEGGIds)){
            names1<- getNames(vectorOfKEGGIds[count])
            names <- append(names,names1)
            count <- count + 1;
        }
}
       return <- names;

}

getNames <- function(geneId){

    ### Vérifiez la connection internet
     url <- getGeneInfoUrl(geneId)
     foundName <- FALSE;
     allLines <- readLines(url);
     i <- 1;
     name<- NULL;

     while(foundName == FALSE){

         allLines[i] <- stringr::str_trim(allLines[i], "both")
         tmp <- strsplit(allLines[i], "\\s+|,|;")

         if(!is.null(tmp[[1]][1])){
            if(tmp[[1]][1] =="NAME"){

             name <- tmp[[1]][2]
             foundName <- TRUE;

         }}

         i <- i+1;
     }
     return <- name
}

getGeneInfoUrl <- function(geneId){

 url <- "http://rest.kegg.jp/get/"
 url <- paste(url, geneId, sep = "")
 return <- url;
}


getFirstLevelSubMaps <- function(pathwayId) {

    xmltop <- getKGMLRootNode(pathwayId);


    mapId <- XML::xpathSApply(xmltop, "//entry[@type = 'map']",
                                   function(x) (XML::xmlAttrs(x))['name']);

return <- mapId
}






isAllGeneOfSubGraphsInKEGGOverview <- function(subGraphId){

    overviewGeneList <- KEGGREST::keggLink("genes", "hsa01100")


    subGraphGeneList <- KEGGREST::keggLink("genes", subGraphId)

    notInOverview <- c();
    for(i in 1:length(subGraphGeneList)){

        subGraphGene_Grep<- paste("\\",subGraphGeneList[[i]],"\\b",sep="")

        r <- grep(subGraphGene_Grep, overviewGeneList);

        if(length(r) == 0){
            notInOverview <- c(notInOverview,subGraphGeneList[[i]])
        }

    }

    return(notInOverview)

}

allSubMapsMetaboliteAnalysis <- function(){

    allSubMaps <- getFirstLevelSubMaps("hsa01100")
    analysisResults<- NULL;

    for(i in 1:length(allSubMaps)){
        allSubMaps[[i]] <- gsub("path:", "", allSubMaps[[i]])

        if(substr(allSubMaps[[i]],1,3) == "hsa"){
            mapSubMapId <- gsub("hsa", "map", allSubMaps[[i]])
       numberGeneOfSubMap <- length(KEGGREST::keggLink("cpd",mapSubMapId))

     notInOverview  <- isAllMetaboliteOfSubGraphsInKEGGOverview(mapSubMapId)
     subMapName <- KEGGREST::keggGet(allSubMaps[[i]])[[1]]$NAME
     subMapClass <- KEGGREST::keggGet(allSubMaps[[i]])[[1]]$CLASS

        if(length(notInOverview) == 0){
            notInOverview <- NA;
            numberGeneNotInOverview <- 0;
        }else{
        numberGeneNotInOverview <- length(notInOverview)
        notInOverview <- paste(notInOverview, collapse = ',')
        }
        analysisResults <-rbind(analysisResults,
              c(subMapKEGGId = as.vector(allSubMaps[[i]]),
                subMapName = as.vector(subMapName),
                subMapClass = as.vector(subMapClass),
                numberMetaboliteOfSubMap = as.vector(numberGeneOfSubMap),
                numberMetaboliteNotInOverview = as.vector(numberGeneNotInOverview),
                metaboliteNotInOverviewKEGGId = as.vector(notInOverview)))


        }
    }

   return <- analysisResults;

}

isAllMetaboliteOfSubGraphsInKEGGOverview <- function(subGraphId){

    overviewMetaboliteList <- KEGGREST::keggLink("cpd", "map01100")


    subGraphMetaboliteList <- KEGGREST::keggLink("cpd", subGraphId)

    notInOverview <- c();

    if(length(subGraphMetaboliteList)>0){
    for(i in 1:length(subGraphMetaboliteList)){

        subGraphMetabolite_Grep<-
            paste("\\",subGraphMetaboliteList[[i]],"\\b",sep="")

        r <- grep(subGraphMetabolite_Grep, overviewMetaboliteList);

        if(length(r) == 0){
            notInOverview <- c(notInOverview,subGraphMetaboliteList[[i]])
        }

    }
   }
    return(notInOverview)

}

allSubMapsGeneAnalysis <- function(){

    allSubMaps <- getFirstLevelSubMaps("hsa01100")
    analysisResults<- NULL;
    for(i in 1:length(allSubMaps)){
        allSubMaps[[i]] <- gsub("path:", "", allSubMaps[[i]])
        if(substr(allSubMaps[[i]],1,3) == "hsa"){
            numberGeneOfSubMap <- length(KEGGREST::keggLink("genes",allSubMaps[[i]]))

            notInOverview  <- isAllGeneOfSubGraphsInKEGGOverview(allSubMaps[[i]])
            subMapName <- KEGGREST::keggGet(allSubMaps[[i]])[[1]]$NAME
            subMapClass <- KEGGREST::keggGet(allSubMaps[[i]])[[1]]$CLASS

            if(length(notInOverview) == 0){
                notInOverview <- NA;
                numberGeneNotInOverview <- 0;
            }else{
                numberGeneNotInOverview <- length(notInOverview)
                notInOverview <- paste(notInOverview, collapse = ',')
            }
            analysisResults <-rbind(analysisResults,
                                    c(subMapKEGGId = as.vector(allSubMaps[[i]]),
                                      subMapName = as.vector(subMapName),
                                      subMapClass = as.vector(subMapClass),
                                      numberGeneOfSubMap = as.vector(numberGeneOfSubMap),
                                      numberGeneNotInOverview = as.vector(numberGeneNotInOverview),
                                      geneNotInOverviewKEGGId = as.vector(notInOverview)))



        }
    }

    return <- analysisResults;

}


