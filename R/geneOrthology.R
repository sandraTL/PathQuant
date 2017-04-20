

getGeneBriteList <- function(){

    #print("getGeneBriteList")

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

   # print("getBriteListByKEGGId")

    if(nrow(geneBriteDF) == 0){ geneBriteDF <- getGeneBriteList() }

    l <- lapply(associations[,1], function(x){
            subset(geneBriteDF, geneBriteDF[,1] == x)})

    l<- do.call(rbind,l)

    return <- l

}

# get the list of brites without duplicates for the input list of gene in
# dataset to get url data info.
getBritesListForUrl <- function(associations){

   # print("getBritesListForUrl")

    geneBriteDF <- getBriteListByKEGGId(associations, data.frame())

    briteList <- data.frame("brite" = geneBriteDF[,2])


    briteList <- data.frame("brite" = briteList[!duplicated(briteList),])

    return <- briteList
}

getBriteListWithNa <- function(associations){

   # print("getBriteListWithNa")

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

getBriteDefByGeneKeggId <- function(gene, brite){

   # print("getBriteDefByGeneKeggId")

    file_name <- paste(brite, ".txt", sep="")

    d <- read.delim(file_name, sep="\n")

    d <- data.frame("line" = d[,1])

    gene <- gsub("hsa:", "", gene)

    gene <- paste(" ", gene, " ", sep="")

    gene_line_num <- grep(gene,d$line)

    defineBriteByLineNum(gene_line_num, d)

}



## I got change it here!!!!
getBriteDefDF <- function(associations_annotated, association){

    # print("getBriteDefDF")

   # annoted Association with brite and EC
   # associations_annotated <- annotateAssociationData(association)
   # print(nrow(associations_annotated))
     # get the list of brite (no duplicate) used the list
    list <- getBritesListForUrl(association)

    briteDF <- data.frame()
    # download and save text of the brite rest url
    apply(list,1, function(x){
        url <- paste("http://rest.kegg.jp/get/", x, sep = "")
        file <- downloadFileByUrl(url)
        file_name <- paste(x,".txt", sep="")
        write(file, file_name)
    })

    for(i in 1:nrow(associations_annotated)){

        r <- strsplit(as.character(associations_annotated[i,3]), "\\s+")

        briteCN <- character()
        for(j in 2:length(r[[1]])){
          re <- "";
         if(!is.na(r[[1]][j])){
          re <- getBriteDefByGeneKeggId(associations_annotated[i,2], r[[1]][j])
          re <-  paste(as.vector(unlist(re)), collapse=" ")
          briteCN <- paste(briteCN, re, sep = " | ")

          }
        }
        if(length(briteCN) == 0) briteCN <- NA
        briteDF <- rbind(briteDF, "briteName" = data.frame(as.vector(briteCN)))

    }

   briteDF <- data.frame("briteName" = as.vector(briteDF))
 #  return <- briteDF
}

defineBriteByLineNum <- function(geneLN, fileDF){

    # print("defineBriteByLineNum")

    lines_Beg_A <- grep("^A.*",fileDF$line)
    lines_Beg_B <- grep("^B.*",fileDF$line)
   # lines_Beg_C <- grep("^C.*",fileDF$line)

    geneLN_split <- strsplit(as.character(geneLN), "[\\s+]")
    geneLN_split <- do.call(rbind,geneLN_split)

    list <- data.frame()

     for(i in 1:nrow(geneLN_split)){

         line_A <- findPlaceOfNumInList(lines_Beg_A, geneLN_split[i,])

        # line_B <- findPlaceOfNumInList(lines_Beg_B, geneLN_split[i,])

         # def <- paste(as.character(fileDF[line_A,]),
         #              as.character(fileDF[line_B,], sep= " "))
         def <- paste(as.character(fileDF[line_A,]), sep= " ")

         list <- rbind(list,  data.frame("brite" = as.vector(def)))

     }

   list <- list[!duplicated(list),]

   return <- list
}


findPlaceOfNumInList <- function(list, brite_num){

   # print("findPlaceOfNumInList")
    place <- length(list);

    for(i in 1:(length(list))){

        if(as.numeric(as.character(brite_num)) < list[i]){
            place <- i;
            break;
        }

    }

   return <- list[place-1];

}




