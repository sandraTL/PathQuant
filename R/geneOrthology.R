
getGeneBriteList <- function(){

    # print("getGeneBriteList")

     url = "http://rest.kegg.jp/link/br/hsa"
     urlFile <- downloadFileByUrl(url)

     geneBriteDF <- strsplit(urlFile, "[\n]");

     geneBriteDF <- lapply(geneBriteDF, function(x) strsplit(x, "[\t]"))

     geneBriteDF <- do.call(rbind, geneBriteDF[[1]])

     return <- geneBriteDF

 }

# get the brites for every gene in the associations data set
# adding line if gene as more than one brite
getBriteListByKEGGId <- function(associations, geneBriteDF){

    #print("getBriteListByKEGGId")

    if(nrow(geneBriteDF) == 0){ geneBriteDF <- getGeneBriteList() }

    l <- lapply(associations[,1], function(x){
            subset(geneBriteDF, geneBriteDF[,1] == x)})

    l<- do.call(rbind,l)

    return <- l

}

getBriteListWithNa <- function(associations){

    #print("getBriteListWithNa")

    df <- list()

    geneBriteDF <- getGeneBriteList()

    for(i in 1:nrow(associations)){

        r <- getBriteListByKEGGId(associations[i,], geneBriteDF)

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
