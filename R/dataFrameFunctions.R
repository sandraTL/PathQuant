#'merge rows with smallest values in each col for repeated metabolites or genes
#'Genes and metabolites are likely to be repeated since we can find them in
#'multiple places in the graph reprensentation of a map

mergeRowsWithSmallestValueByKEGGId <- function(df){

    metaboName <- as.vector(df$KEGGId)

    replicat <- unique(metaboName[duplicated(metaboName)])

    #colIds <- colnames(df)
    for(x in replicat){

        #find rows ID of this replicated element
        rowReplicat <- which(df$KEGGId == x)
      #  print("dataframe 15")

        keggId <- data.frame(KEGGId = df$KEGGId);

       # print(df)
        if(ncol(df) > 2){
           # print("1")
            df <- df[,-ncol(df)]

        rowReplicatToReplace <- rowReplicat[1]

        # row to remove
        rowReplicatToRemove <- rowReplicat[-1]

        # subset of rows with the replicated value
        rowReplicatSubDF <- df[rowReplicat,]
        rowReplicatSortedSubDF <- apply(rowReplicatSubDF,2,sort)
        rowMinValues <- rowReplicatSortedSubDF[1,]
        }else{
       # print("2")
        df <- data.frame("pl" = df[,-ncol(df)])
        rowReplicatToReplace <- rowReplicat[1]

            # row to remove
        rowReplicatToRemove <- rowReplicat[-1]
            # subset of rows with the replicated value
        rowReplicatSubDF <- data.frame("pl" = df[rowReplicat,])
       # print("2_1")
       # print(dim(rowReplicatSubDF))
        rowReplicatSortedSubDF <- data.frame("pl"=
                                rowReplicatSubDF[order(rowReplicatSubDF[,1]),])
       # print("2_2")
       # print(rowReplicatSortedSubDF)
       # print("2_3")
      #  print(rowReplicatSubDF)
        rowMinValues <- rowReplicatSortedSubDF[1,]
        }

        # row ID to be replaced with smallest value

        # sort all columns by ascending order (more than 1 column)

        # keep only the first row -> the one with smallest vlaues


        df[rowReplicatToReplace,] <- rowMinValues;

        df <- cbind(df, keggId)
        df <- df[-rowReplicatToRemove,]
       # print(df)

     }

    return <- df;

}

merge2DFWithSmallestValue <- function(df1, df2){

for (row in 1:nrow(df)) {
    r <- mergeVectorsLowerValues(df1[row,], df2[row,]);
    rf <- t(data.frame(r));
    finalDF <- rbind(finalDF,rf);
    return <- finalDF;
}
}
removeDuplicatedColumnDF <- function(df){

df[is.na(df)] <- Inf;

    df<- t(df)


    df<- aggregate(df,list(row.names(df)),function(x) x[which.min(abs(x))]);

    row.names(df) <- df[,1]
    df <- data.frame(t(df))

    df <- df[-c(1), ]

return <- df;
}

removeRowsDistanceAsso <- function(df){

    # df$geneKEGGId <- factor(df$metabolites, levels=unique(df$metabolites))
    df$geneKEGGId <- I(df$geneKEGGId)
    df$metaboliteKEGGId <- I(df$metaboliteKEGGId)

    aa <- split(df, list(df$metaboliteKEGGId,df$geneKEGGId),drop = TRUE)

    f <- lapply(aa ,function(x){

        tempDF <- data.frame(x)

        tempDF <- tempDF[order(tempDF[,5]),]
        return <-tempDF[1,]
    })

    finalDF <- do.call(rbind.data.frame, f)
    row.names(finalDF) <- row.names(1:length(finalDF[,1]))
    return <- finalDF;

}

#' Thsi function also sorts the data by metabolite or gene ids... we don't want
#' that
removeRowsDistanceAll <- function(df){

    df$geneKEGGId <- I(df$geneKEGGId)
    aa <- split(df, list(df$geneKEGGId),drop = TRUE)
    nbreOfCol <- ncol(df) -1

    f <- lapply(aa ,function(x){

        tempDF <- data.frame(x)
        tempDF <- apply(tempDF,2,sort,decreasing=F)

        if((length(tempDF)/ncol(df)) == 1){
            tempDF <- as.data.frame(tempDF)
            tempDF <- as.data.frame(t(tempDF))
        }else tempDF <- as.data.frame(tempDF)
        return <-tempDF[1,]
    })

    finalDF <- do.call(rbind.data.frame, f)
    return <- finalDF;

}

getMaxValIgnoreInfVal <- function(df){

    infVal <- which(df[,1] == Inf)
    temp <- df
    temp[infVal,] <- -1
    maxVal <- max(temp[,1])

    return <- maxVal;

}
