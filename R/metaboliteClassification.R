
# metabolite association using HMDB

metaboliteClasseCount <- function(associations){

   # print("metaboliteClasseCount")

    metaboTable <-table(associations$metaboliteSuperClasse,useNA = "always")
    metaboNames <- as.vector(names(metaboTable))
    metaboCount <- as.vector(metaboTable)
    metaboClasse <- data.frame("classNames" = metaboNames,
                               "count" = metaboCount)
    metaboClasse <- metaboClass[apply(metaboClass[2], 1, function(x) !any(x==0)),]
    return <- metaboClasse;

}

annotateMetabolite <- function(associations){

  #  print("annotateMetabolite")

    briteDF <- data.frame("brite" =
                             as.vector(getMetaboBriteListWithNa(associations)))

    df <- data.frame("gene" = as.vector(associations[,1]),
                     "metabolite" = as.vector(associations[,2]),
                     "brite" = briteDF)

    classDF <- data.frame("brite" = as.vector(getMetaboClassification(df)))
   # defDF <- data.frame("def" = as.vector(getBriteDefDF(df, associations)))
   # print(defDF)

    df <- data.frame("gene" = as.vector(associations[,1]),
                     "metabolite" = as.vector(associations[,2]),
                     "brite" = briteDF,
                     "classification"= classDF)

    return <- df
}


# classification for pathways enrichment using compound brite
getMetaboliteBriteList <- function(){

  #  print("getMetaboliteBriteList")

    url = "http://rest.kegg.jp/link/br/cpd"
    urlFile <- downloadFileByUrl(url)

    metaboBriteDF <- strsplit(urlFile, "[\n]");

    metaboBriteDF <- lapply(metaboBriteDF, function(x)
        gsub("cpd:", "", x))
    metaboBriteDF <- lapply(metaboBriteDF, function(x) strsplit(x, "[\t]"))

    metaboBriteDF <- do.call(rbind, metaboBriteDF[[1]])



    return <- metaboBriteDF

}

# get the brites for every gene in the associations data set
# adding line if gene as more than one brite
getMBriteListByKEGGId <- function(associations, metaboBriteDF){

   # print("getMBriteListByKEGGId")

    if(nrow(metaboBriteDF) == 0){ metaboBriteDF <- getMetaboliteBriteList() }

    l <- lapply(associations[,2], function(x){
        subset(metaboBriteDF, metaboBriteDF[,1] == x)})

    l<- do.call(rbind,l)


    return <- l

}

# get the list of brites without duplicates for the input list of metabolite in
# dataset to get url data info.
getBritesListForUrl <- function(associations){

  #  print("getBritesListForUrl")

    metaboBriteDF <- getMBriteListByKEGGId(associations, data.frame())

    briteList <- data.frame("brite" = metaboBriteDF[,2])

    briteList <- data.frame("brite" = briteList[!duplicated(briteList),])

    return <- briteList
}


getMetaboBriteListWithNa <- function(associations){

   # print("getMetaboBriteListWithNa")

    df <- list()

    metaboBriteDF <- getMetaboliteBriteList()



    for(i in 1:nrow(associations)){

        r <- getMBriteListByKEGGId(associations[i,], metaboBriteDF)

        if(length(r) == 0){
            row <- data.frame("brite"   = NA);
        }else if(length(r[,1]) > 1){

            r<- as.data.frame(r)
            r <- concatDfColInfoFromDuplicate(r, 1, 2)
            row <- data.frame("brite"   = r[1,2]);
        }else if(length(r[,1] == 1)){
            row <- data.frame("brite"   = r[1,2]);
        }

        df <- rbind(df, row)
    }

    return <- df

}

getMetaboClassification <- function(associations_annotated){

   # print("getMetaboClassification")

    briteDF <- NULL;

    for(i in 1:nrow(associations_annotated)){

        #r <- strsplit(as.character(associations_annotated[i,3]), "\\s+")
        # print(grep("br:br08001", associations_annotated[i,3]))
        # print(grep("br:br08002", associations_annotated[i,3]))
        briteCN <- character();

        re <- NULL;
       # if(is.na(r[[1]][1])){ re <- c(re, "no brite classification")

       # }

         if(length(grep("br:br08002", associations_annotated[i,3])) == 1){
             re <-c(re, "lipid")
        }
        else if(length(grep("br:br08001", associations_annotated[i,3])) == 1){
                      print(associations_annotated[i,2])

             meta <- getBriteDefByGeneKeggId(associations_annotated[i,2],
                                             "br:br08001")

             re <-c(re, as.vector(meta))

         }
         else if(length(grep("br:br08003", associations_annotated[i,3])) == 1){
            re <-c(re, "Phytochemical Compounds")
         }
         else if(length(grep("br:br08007", associations_annotated[i,3])) == 1){
            re <-c(re, "Pesticides")
         }
         else if(length(grep("br:br08303", associations_annotated[i,3])) == 1){
            re <-c(re, "Anatomical Therapeutic Chemical")
         }
         else if(length(grep("br:br08323", associations_annotated[i,3])) == 1){
            re <-c(re, "Major components of natural products")
         }
         else re <- c(re, "no brite classification")

        # for(j in 2:length(r[[1]])){
        #     if(!is.na(r[[1]][j])){
        #         if(r[[1]][j] == "br:br08002"){re <-c(re, "lipid")}
        #         else if(r[[1]][j] == "br:br08001"){
        #             print(associations_annotated[i,2])
        #
        #             meta <- getBriteDefByGeneKeggId(associations_annotated[i,2], r[[1]][j])
        #             print(meta)
        #             re <-c(re, as.vector(meta))
        #         }
        #
        #
        #
        #     }
        #
        # }
        if(is.null(re)) re <- NA

        re <-  paste(as.vector(unlist(re)), collapse=" ")
        briteDF <- rbind(briteDF, "Classification" = data.frame(as.vector(re)))
    }

    briteDF <- data.frame("Classification" = as.vector(briteDF))
    return <- briteDF
}
