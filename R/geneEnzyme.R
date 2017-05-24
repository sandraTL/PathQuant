getEcListWithNa <- function(associations){

    # print("getEcListWithNa")

    df <- list()

    url <- "http://rest.kegg.jp/get/br:hsa01000"
    urlFile <- downloadFileByUrl(url)

    briteEcDF <- strsplit(urlFile, "[\n]");

    briteEcDF <- do.call(rbind, briteEcDF)

    briteEcDF <- data.frame(t(briteEcDF))

    enzymeVec <- NULL

    for(i in 1:nrow(associations)){

        # See if the gene is annotated to the br:hsa01000
        # The br:hsa01000 is the brite for enzymes
        isEc <- grep("br:hsa01000", as.character(associations[i,2]))

        # if the gene encodes an enzymes, than find the EC number
        if(length(isEc) == 1){

            #get the numerical part of hsa:XXXXX
            geneId <- gsub("hsa:", "", associations[i,1])
            geneId<- paste("E        ", geneId,"\\b",sep="")

            ec.line.df <- data.frame(briteEcDF[grep(geneId, briteEcDF[,1]),])

            ecList <- NULL

            # extract br:hsa01000 file to find ec number
            # for each gene
            if(nrow(ec.line.df) > 0){

             for(j in 1:nrow(ec.line.df)){

                    # get string of this line in one String
                    line <- paste(ec.line.df[j,], sep="", collapse="")

                    # get the start and end position of EC number [ ]
                    startCharEc <-  gregexpr('\\[', line)
                    endCharEc <-  gregexpr('\\]', line)

                    # get the number
                    ecNumber <- substr(line,
                                       startCharEc[[1]][1],
                                       endCharEc[[1]][1])

                    # put the ec number in a list
                    ecList <- c(ecList, ecNumber)

                }
            }

            # remove replicat enzymes in the list
            ecList <- unique(ecList)
            ecList <- paste(ecList, sep="", collapse="")

            #
            enzymeVec <- c(enzymeVec,ecList)

       } else {
           enzymeVec <- c(enzymeVec, "NA")
       }

    }

    return <- enzymeVec ;

}
