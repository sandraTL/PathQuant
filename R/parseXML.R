
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
    nodeDF <- data.frame("kgmlId" = as.vector(nodeListId),
                         "keggId" = as.vector(nodeListCpd));

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




getKGMLRootNode <- function(pathwayId){
    # get the root of the KGML document
    pathFile <- toStringPathFile(pathwayId);

    if(is.na(file.info(pathFile)$size)== FALSE){
        xmlfile <- XML::xmlParse(pathFile);
        xmltop <- XML::xmlRoot(xmlfile); # gives content of root
    }else
        xmltop = NULL;

    #
    return <- xmltop;

}



toStringPathFile <- function(pathwayId){

    # concatenation of pathwayId to set swdir for the xml
    s1 <- "~/";
    s2 <-  toString(pathwayId);
    s3 <- ".xml"
    s4 <- paste(s1,s2, sep= "");
    pathFile <- paste(s4, s3, sep="");

    return <- pathFile;
}

getCommonNames <- function(vectorOfKEGGIds, type = c("gene","metabolite")){

        count <- 1;
        names <- character();

        while(count <= length(vectorOfKEGGIds)){
            names1<- getNames(vectorOfKEGGIds[count])
            names <- append(names,names1)
            count <- count + 1;
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

         allLines[i] <- str_trim(allLines[i], "both")
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
