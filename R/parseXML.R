
# This function extract from KGML a list of metabolites (id and name) that
# makes the nodes of the graph ** some metablite can from more than 1 node,
# which is why the id of the nodes is not the name of the metabolite.
# data.frame (id, keggId)

getListNodeFromKGML <- function(pathwayId) {

    # print("getListNodeFromKGML")

    xmltop <- getKGMLRootNode(pathwayId);

    nodeListId <- XML::xpathSApply(xmltop,
                                   "//entry[@type = 'compound']",
                                   function(x) (XML::xmlAttrs(x))['id']);
    nodeListCpd <- XML::xpathSApply(xmltop,
                                    "//entry[@type = 'compound']",
                                    function(x) (XML::xmlAttrs(x))['name']);
    nodeListColor <- XML::xpathSApply(xmltop,
                                      "//entry[@type = 'compound']//graphics",
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

    #print("getListReactionFromKGML")
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

    if(length(reactionIdNodes) > 0) {
        max1 <- max(lengths(substrateId));
        max2 <- max(lengths(productId));
        max3 <- max(lengths(substrateName));
        max4 <- max(lengths(productName));
        max5 <- max(lengths(reactionId));
        max6 <- max(lengths(reactionType));
        max7 <- max(lengths(reactionName));



        sum <- (max1 + max2 + max3 +max4 + max5 + max6 +max7);

        # print(substrateId)
        # print(productId)
        # print(substrateName)
        # print(productName)
        # print(reactionId)
        # print(reactionType)
        # print(reactionName)

        # print(length(substrateId));
        # print(length(productId));
        # print(length(substrateName));
        # print(length(productName));
        # print(length(reactionId));
        # print(length(reactionType));
        # print(length(reactionName));

        # print(pathwayId)

        if(sum > 7) {
            reactionList <- do.call(rbind.data.frame,
                                    mapply(cbind,
                                           "substrateId" = substrateId,
                                           "productId" =  productId,
                                           "substrateName" = substrateName,
                                           "productName" = productName,
                                           "reactionId" =  reactionId,
                                           "reactionType" = reactionType,
                                           "reactionName" = reactionName
                                    ));
        } else if(sum == 7) {
            reactionList <- data.frame(cbind(
                "substrateId" = as.numeric(unlist(substrateId)),
                "productId" = as.numeric(unlist(productId)),
                "substrateName" = as.character(unlist(substrateName)),
                "productName" = as.character(unlist(productName)),
                "reactionId" = as.numeric(reactionId),
                "reactionType" = reactionType,
                "reactionName" = reactionName
            ));
        }
    } else {
        reactionList <- data.frame();
    }

    return <- reactionList;

}


getListOfUniqueHSAGeneId <- function(pathwayId){

    print("getListOfUniqueHSAGeneId")

    gene.list <- getListEdgeFromGeneKGML(pathwayId)$ko
    gene.list <- paste(gene.list, collapse = " ")
    gene.list <- strsplit(gene.list, " ")
    gene.list <- gene.list[[1]][!duplicated(gene.list[[1]])]

    return <- gene.list
}

getListOfUniqueMetabolite <- function(pathwayId){

    print("getListOfUniqueMetabolite")

    metabo.list <- getListNodeFromKGML(pathwayId)$keggId
    metabo.list <- paste(metabo.list, collapse = " ")
    metabo.list <- strsplit(metabo.list, " ")
    metabo.list <- metabo.list[[1]][!duplicated(metabo.list[[1]])]

    return <- metabo.list
}

getListUniqueMetaboliteInReactions <- function(pathwayId){

    print("getListUniqueMetaboliteInReactions")

    reac.sub.list <- getListReactionFromKGML(pathwayId)$substrateName
    reac.prod.list <- getListReactionFromKGML(pathwayId)$productName
    reac.sub.list <- paste(reac.sub.list, collapse = " ")
    reac.prod.list <- paste(reac.prod.list, collapse = " ")

    reac.me.list <- paste(reac.sub.list, reac.prod.list, collapse = " ")

    reac.me.list <- strsplit(reac.me.list, " ")
    reac.me.list <- reac.me.list[[1]][!duplicated(reac.me.list[[1]])]

    return <- reac.me.list

}


# This function extract from KGML a list of reaction of entry-type "gene"
# Returns a dataFrame object : id, entryId, reaction, ko

getListEdgeFromGeneKGML <- function(pathwayId) {


    # print("getListEdgeFromGeneKGML")
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

     print("getListOrthologGeneFromKGML")
    # get the root of the KGML document
    xmltop <- getKGMLRootNode(pathwayId);

    # Get value of atributes (id, name(ko) and reaction) of entry of type
    # ortholog, some don't have reaction argument thus won't have an edge in
    # our final graph they will have NA in data.frame.
    orthologListId <- XML::xpathSApply(xmltop,
                                       "//entry[@type = 'ortholog']",
                                       function(x) (XML::xmlAttrs(x))['id']);

    orthologListKo <- XML::xpathSApply(xmltop,
                                       "//entry[@type = 'ortholog']",
                                       function(x) (XML::xmlAttrs(x))['name']);

    orthologListReaction <- XML::xpathSApply(xmltop,
                                             "//entry[@type = 'ortholog']",
                                             function(x) (XML::xmlAttrs(x))['reaction']);

    orthologDF <- data.frame(
        "reactionId" = as.vector(as.character(orthologListId)),
        "reactions" = as.vector(as.character(orthologListReaction)),
        "ko" = as.vector(as.character(orthologListKo)),
        "x" = as.vector(rep(-1, length(orthologListReaction))),
        "y" = as.vector(rep(-1, length(orthologListReaction)))
    );

    for(row in 1:length(orthologDF[,1])){
        id <- orthologDF[row,1]

        coords <- XML::xpathSApply(xmltop,
                                   "//entry[@id = id]//graphics",
                                   function(x) (XML::xmlAttrs(x))['coords']);

        if(length(coords) == 1){
            orthologDF[row,4] <- coords;
        }
    }

    orthologListNameCoords <- XML::xpathSApply(xmltop, "//entry[@type = 'ortholog']//graphics",
                                               function(x) (XML::xmlAttrs(x))['name']);

    orthologListCoords <- XML::xpathSApply(xmltop, "//entry[@type = 'ortholog']//graphics",
                                           function(x) (XML::xmlAttrs(x))['coords']);

    orthologCoords <- data.frame("nameCoords" =as.vector(as.character(orthologListNameCoords)),
                                 "coords" = as.vector(as.character(orthologListCoords)));

}



getKGMLRootNode <- function(pathwayId){
    # get the root of the KGML document
    # print("getKGMLRootNode")

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

    # print("toStringPathFile")
    # concatenation of pathwayId to set swdir for the xml

    s2 <-  toString(pathwayId);
    s3 <- ".txt"

    pathFile <- paste(s2, s3, sep="");

    return <- pathFile;
}

getPathCommonNames <- function(path){

     print("getPathCommonNames")

    pathSplit <- strsplit(as.character(path[1,]), " ")

    for(i in 1:length(pathSplit[[1]])){

        if(substring(pathSplit[[1]][i],1,1) == "C" && nchar(pathSplit[[1]][i]) == 6){
            pathSplit[[1]][i] <- getCommonNames(pathSplit[[1]][i], "metabolite")
        }else if(substring(pathSplit[[1]][i],1,4) == "hsa:"){
            pathSplit[[1]][i] <- getCommonNames(pathSplit[[1]][i], "gene")
        }
    }

    pathSplit <- do.call(paste, c(as.list(pathSplit[[1]]), sep=" "));

    return <- pathSplit;
}



getCommonNames <- function(vectorOfKEGGIds, type = c("gene","metabolite")){

    # print("getCommonNames")
    count <- 1;

    ### Vérifiez la connection internet

    if(length(vectorOfKEGGIds) > 0 ){

        # analysis of submaps - beggining
        names <- character();

        while(count <= length(vectorOfKEGGIds)){

            names1<- getNames(vectorOfKEGGIds[count], type)

            names <- append(names,names1)

            count <- count + 1;
        }
    }


    return <- names;

}

getNames <- function(keggId, type){

    # print("getNames")
    # print(keggId)
    ### Vérifiez la connection internet

    name <- NULL;

    if (is.na(keggId)) {

        name <- NA;



    } else {

        url <- getKEGGInfoUrl(keggId)

        foundName <- FALSE;
        allLines <- readLines(url);
        i <- 1;

        while(foundName == FALSE){

            allLines[i] <- stringr::str_trim(allLines[i], "both")

            if(type == "gene"){

                tmp <- strsplit(allLines[i], "\\s+|,|;")



                if(!is.null(tmp[[1]][1])){
                    if(tmp[[1]][1] =="NAME"){

                        name <- tmp[[1]][2]
                        foundName <- TRUE;

                    }}

                i <- i+1;
            }
            else if(type == "metabolite"){

                tmp <- strsplit(allLines[i], "\\s+|;")

                if(!is.null(tmp[[1]][1])){
                    if(tmp[[1]][1] =="NAME"){

                        name <- tmp[[1]][2]
                        foundName <- TRUE;

                    }}

                i <- i+1;

            }

        }
    }
    return <- name
}

# trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# getBrites <- function(Associations){
#      print("getBrites")
#
#     brite<- character();
#     # Associations <- Associations[,c(2,5)]
#     # Associations <- Associations[!duplicated(Associations),]
#     for(j in 1:nrow(Associations)){
#
#         url <- getKEGGInfoUrl(Associations[j,3])
#
#         foundName <- FALSE;
#         allLines <- readLines(url);
#         i <- 1;
#
#         for(i in 1:length(allLines)){
#             temp <- character();
#             tmp <- strsplit(allLines[i], "\\s+|;")
#             if(!is.null(tmp[[1]][1])){
#
#                 if(tmp[[1]][1] =="BRITE"){
#                     foundName <- TRUE;
#                     if(tmp[[1]][2] == "Compounds"){
#                         temp <- paste(trim(allLines[i+1]),
#                                       "-", trim(allLines[i+2]),sep = " ");
#
#                     }else if(tmp[[1]][2] == "Anatomical"){
#                         temp <-paste(tmp[[1]][2],"-",trim(allLines[i+4]),sep = " ")
#                     }else if(tmp[[1]][2] == "Pharmaceutical"){
#                         temp <-paste(tmp[[1]][2],tmp[[1]][3],
#                                      sep = " ")
#                     }else if(tmp[[1]][2] == "Phytochemical"){
#                         temp <-paste(tmp[[1]][2],"-",trim(allLines[i+1]),
#                                      sep = " ")
#                     }else if(tmp[[1]][2] == "Pesticides"){
#                         temp <-paste(tmp[[1]][2],"-",trim(allLines[i+1]),
#                                      sep = " ")
#                     }else{
#                         temp <-paste(tmp[[1]][2],"-",trim(allLines[i+1]), sep = " ")
#                     }
#
#                 }
#                 brite <- append(brite,temp)
#             }
#
#         }
#
#         if(foundName==FALSE){
#             temp <- "not classified by KEGG yet"
#
#             brite <- append(brite,temp)
#         }
#
#     }
#
#     Associations <- cbind(Associations, "Chemical Class" = as.vector(brite))
#     return <- Associations;
# }


getKEGGInfoUrl <- function(keggId){

    # print("getKEGGInfoUrl")

    url <- "http://rest.kegg.jp/get/"
    url <- paste(url, keggId, sep = "")
    return <- url;
}

#
# getFirstLevelSubMaps <- function(pathwayId) {
#
#       print("getFirstLevelSubMaps")
#
#     xmltop <- getKGMLRootNode(pathwayId);
#
#
#     mapId <- XML::xpathSApply(xmltop, "//entry[@type = 'map']",
#                               function(x) (XML::xmlAttrs(x))['name']);
#
#     return <- mapId
# }
#
#
#
#
#
#
# isAllGeneOfSubGraphsInKEGGOverview <- function(subGraphId){
#
#     print("isAllGeneOfSubGraphsInKEGGOverview")
#
#     overviewGeneList <- KEGGREST::keggLink("genes", "hsa01100")
#
#
#     subGraphGeneList <- KEGGREST::keggLink("genes", subGraphId)
#
#     notInOverview <- c();
#     for(i in 1:length(subGraphGeneList)){
#
#         subGraphGene_Grep<- paste("\\",subGraphGeneList[[i]],"\\b",sep="")
#
#         r <- grep(subGraphGene_Grep, overviewGeneList);
#
#         if(length(r) == 0){
#             notInOverview <- c(notInOverview,subGraphGeneList[[i]])
#         }
#
#     }
#
#     return(notInOverview)
#
# }
#
# allSubMapsMetaboliteAnalysis <- function(){
#
#     print("allSubMapsMetaboliteAnalysis")
#
#     allSubMaps <- getFirstLevelSubMaps("hsa01100")
#     analysisResults<- NULL;
#
#     for(i in 1:length(allSubMaps)){
#         allSubMaps[[i]] <- gsub("path:", "", allSubMaps[[i]])
#
#         if(substr(allSubMaps[[i]],1,3) == "hsa"){
#             mapSubMapId <- gsub("hsa", "map", allSubMaps[[i]])
#             numberGeneOfSubMap <- length(KEGGREST::keggLink("cpd",mapSubMapId))
#
#             notInOverview  <- isAllMetaboliteOfSubGraphsInKEGGOverview(mapSubMapId)
#             subMapName <- KEGGREST::keggGet(allSubMaps[[i]])[[1]]$NAME
#             subMapClass <- KEGGREST::keggGet(allSubMaps[[i]])[[1]]$CLASS
#
#             if(length(notInOverview) == 0){
#                 notInOverview <- NA;
#                 numberGeneNotInOverview <- 0;
#             }else{
#                 numberGeneNotInOverview <- length(notInOverview)
#                 notInOverview <- paste(notInOverview, collapse = ',')
#             }
#             analysisResults <-rbind(analysisResults,
#                                     c(subMapKEGGId = as.vector(allSubMaps[[i]]),
#                                       subMapName = as.vector(subMapName),
#                                       subMapClass = as.vector(subMapClass),
#                                       numberMetaboliteOfSubMap = as.vector(numberGeneOfSubMap),
#                                       numberMetaboliteNotInOverview = as.vector(numberGeneNotInOverview),
#                                       metaboliteNotInOverviewKEGGId = as.vector(notInOverview)))
#
#
#         }
#     }
#
#     return <- analysisResults;
#
# }
#
# isAllMetaboliteOfSubGraphsInKEGGOverview <- function(subGraphId){
#
#      print("isAllMetaboliteOfSubGraphsInKEGGOverview")
#
#     overviewMetaboliteList <- KEGGREST::keggLink("cpd", "map01100")
#
#
#     subGraphMetaboliteList <- KEGGREST::keggLink("cpd", subGraphId)
#
#     notInOverview <- c();
#
#     if(length(subGraphMetaboliteList)>0){
#         for(i in 1:length(subGraphMetaboliteList)){
#
#             subGraphMetabolite_Grep<-
#                 paste("\\",subGraphMetaboliteList[[i]],"\\b",sep="")
#
#             r <- grep(subGraphMetabolite_Grep, overviewMetaboliteList);
#
#             if(length(r) == 0){
#                 notInOverview <- c(notInOverview,subGraphMetaboliteList[[i]])
#             }
#
#         }
#     }
#     return(notInOverview)
#
# }
#
# allSubMapsGeneAnalysis <- function(){
#
#      print("allSubMapsGeneAnalysis")
#
#     allSubMaps <- getFirstLevelSubMaps("hsa01100")
#     analysisResults<- NULL;
#     for(i in 1:length(allSubMaps)){
#         allSubMaps[[i]] <- gsub("path:", "", allSubMaps[[i]])
#         if(substr(allSubMaps[[i]],1,3) == "hsa"){
#             numberGeneOfSubMap <- length(KEGGREST::keggLink("genes",allSubMaps[[i]]))
#
#             notInOverview  <- isAllGeneOfSubGraphsInKEGGOverview(allSubMaps[[i]])
#             subMapName <- KEGGREST::keggGet(allSubMaps[[i]])[[1]]$NAME
#             subMapClass <- KEGGREST::keggGet(allSubMaps[[i]])[[1]]$CLASS
#
#             if(length(notInOverview) == 0){
#                 notInOverview <- NA;
#                 numberGeneNotInOverview <- 0;
#             }else{
#                 numberGeneNotInOverview <- length(notInOverview)
#                 notInOverview <- paste(notInOverview, collapse = ',')
#             }
#             analysisResults <-rbind(analysisResults,
#                                     c(subMapKEGGId = as.vector(allSubMaps[[i]]),
#                                       subMapName = as.vector(subMapName),
#                                       subMapClass = as.vector(subMapClass),
#                                       numberGeneOfSubMap = as.vector(numberGeneOfSubMap),
#                                       numberGeneNotInOverview = as.vector(numberGeneNotInOverview),
#                                       geneNotInOverviewKEGGId = as.vector(notInOverview)))
#
#
#
#         }
#     }
#
#     return <- analysisResults;
#
# }
#
# getHMDBByKEGGId <-function(keggId) {
#
#      print("getHMDBByKEGGId")
#
#     files <- list.files("./hmdb_metabolites")
#     files <- files[-1]
#
#     res <-  lapply(keggId, function(i)
#
#         lapply(files, function(j)
#             compareWihtNodeKeggId(j,i)))
#
#
#     return <- res;
#
# }
#
# getHMDBInfo <- function(){
#
#      print("getHMDBInfo")
#
#     files <- list.files("/Users/sandra/Documents/workspaceMetabolomics/hmbdParser/hmdb_metabolites")
#     files <- files[-1]
#
#     hmdb <- getHMDBIds(files)
#
#     kegg_ids <- lapply(files, function(j) getKeggIdByHMDB(j))
#     super_class <- lapply(files, function(j) getSuperClassByHMDB(j))
#     class <- lapply(files, function(j) getClassByHMDB(j))
#     hmdb_infoDF <-data.frame(cbind("hmdb_Id" = as.vector(hmdb),
#                                    "kegg_Id" = as.vector(kegg_ids),
#                                    "super_class" = as.vector(super_class),
#                                    "class" = as.vector(class)))
#     exportHMDBInfo(hmdb_infoDF)
#     return <- hmdb_infoDF
# }
#
#
# getKeggIdByHMDB <- function(fileName){
#
#      print("getKeggIdByHMDB")
#
#     fileName <- paste("/Users/sandra/Documents/workspaceMetabolomics/hmbdParser/hmdb_metabolites/", fileName, sep="")
#     doc1 <- XML::xmlParse(fileName);
#     node_kegg_id <- XML::xpathApply(doc1, "//kegg_id")
#
#     kegg_id <- XML::xmlSApply(node_kegg_id, XML::xmlValue)
#
#     if(nchar(kegg_id) == 0) kegg_id <- NA;
#
#     return <- kegg_id;
# }
#
# getSuperClassByHMDB <- function(fileName){
#
#      print("getSuperClassByHMDB")
#
#     fileName <- paste("/Users/sandra/Documents/workspaceMetabolomics/hmbdParser/hmdb_metabolites/", fileName, sep="")
#     doc1 <- XML::xmlParse(fileName);
#
#     node_super_class <- XML::xpathApply(doc1, "//super_class")
#
#     super_class <- XML::xmlSApply(node_super_class, XML::xmlValue)
#
#     if(nchar(super_class) == 0) super_class <- NA
#
#     return <- super_class
# }
#
# getClassByHMDB <- function(fileName){
#
#      print("getClassByHMDB")
#
#     fileName <- paste("/Users/sandra/Documents/workspaceMetabolomics/hmbdParser/hmdb_metabolites/", fileName, sep="")
#     doc1 <- XML::xmlParse(fileName);
#
#     node_class <- XML::xpathApply(doc1, "//class")
#
#     class <- XML::xmlSApply(node_class, XML::xmlValue)
#
#     if(nchar(class) == 0) class <- NA;
#
#     return <- class;
# }



#Data <- ldply(files,parse_xml)
#
# getHMDBIds <- function(fileList){
#
#      print("getHMDBIds")
#     finalList <- lapply(fileList, function(x) gsub( ".xml", "",x))
#     return <- finalList;
#
# }
#
# exportHMDBInfo <- function(hmdb_infoDF){
#
#      print("exportHMDBInfo")
#     hmdb_infoDF <- data.frame(lapply(hmdb_infoDF, as.character), stringsAsFactors=FALSE)
#     exportDFtoTxt(hmdb_infoDF, "hmdb_info")
# }
#
#
# getHMDBIdByKEggId <- function(keggIdList){
#
#      print("getHMDBIdByKEggId")
#
#     res <- importTXTtoDF("hmdb_info.txt")
#
#     r <-lapply(keggIdList, function(x) subset(res, res[2] == x));
#
#     r <- do.call(rbind, r)
#
#     r <- prepOutputHMDBIdByKEggId(r);
#
#     return <- r[,c(2,1)];
# }


## could also annotate from HMDB database
getSuperClassByKEggId <- function(kegg.id.list){

    # print("getSuperClassByKEggId")

    # retrieve data from HMDB xml files
    hmdb.data <- data.frame(importTXTtoDF("hmdb_info.txt"))
   # hmdb.data <- load("")
   # print(typeof(hmdb.data))
   # print(hmdb.data)
    # unlist the data to for lapply
    kegg.id.list <- unlist(kegg.id.list)
    # print(kegg.id.list)
    r <- lapply(kegg.id.list, function(x) {

        x.1 <- paste("\\b", x,"\\b",sep="")

        if (nrow(hmdb.data[grep(x.1, hmdb.data$kegg_Id),]) == 0) {

            data.frame("hmdb_Id" = NA,
                       "kegg_Id" = x,
                       "super_class" = NA,
                       "class" = NA)

            # carefull I only took the first line
            # but multiple line mean multiple classes that should be
            # bind together
        } else {
            hmdb.data[grep(x.1, hmdb.data$kegg_Id),][1,]
        }
    });

    #bind the results for each compound
    r <- do.call(rbind, r)

    return <- r[,3];
}
#
# getClassByKEggId <- function(keggIdList){
#
#      print("getClassByKEggId")
#
#     res <- importTXTtoDF("hmdb_info.txt")
#
#     r <-lapply(keggIdList, function(x) subset(res, res[2] == x));
#
#     r <- do.call(rbind, r)
#
#     #  r <- prepOutputHMDBIdByKEggId(r);
#
#     return <- r[,c(2,4)];
# }
#
# getAllHMDBbyKeggId <- function(keggIdList){
#
#      print("getAllHMDBbyKeggId")
#
#     res <- importTXTtoDF("hmdb_info.txt")
#
#     r <-lapply(keggIdList, function(x) subset(res, res[2] == x));
#
#     r <- do.call(rbind, r)
#
#     #  r <- prepOutputHMDBIdByKEggId(r);
#
#     return <- r;
# }
#
# ## put this in a general data frame help function
# prepOutputHMDBIdByKEggId <- function(metabolitesNamesDF){
#
#      print("prepOutputHMDBIdByKEggId")
#
#     metabo <- ""
#     hmdb_ids <- ""
#     df <- data.frame();
#
#     for(x in 1:nrow(metabolitesNamesDF)){
#         # conditino to get first row
#         if(metabo == ""){
#             metabo <- metabolitesNamesDF[x,2];
#             hmdb_ids <- paste(hmdb_ids, metabolitesNamesDF[x,1])
#
#             # conditino get inside rows and combine first colinfo info based on
#             # duplicated 2 col
#         }else if(!(metabo == as.character(metabolitesNamesDF[x,2]))){
#
#             df <- rbind(df, data.frame("keggId" = as.vector(as.character(metabo)),
#                                        "hmdbIds" = as.vector(hmdb_ids)))
#             hmdb_ids <- "";
#             metabo <- metabolitesNamesDF[x,2]
#             hmdb_ids <- paste(hmdb_ids, metabolitesNamesDF[x,1])
#
#         }else{
#             hmdb_ids <- paste(hmdb_ids, metabolitesNamesDF[x,1])
#         }
#
#         # ctach last row
#         if(x == nrow(metabolitesNamesDF)){
#             df <- rbind(df, data.frame("keggId" = as.vector(as.character(metabo)),
#                                        "hmdbIds" = as.vector(hmdb_ids)))
#         }
#     }
#
#
#     return <- df;
#
# }


