# merge rows with smallest values in each col for repeated metabolites or genes
# Genes and metabolites are likely to be repeated since we can find them in
# multiple places in the graph reprensentation of a map

mergeRowsWithSmallestValueByKEGGId <- function(df){


    # print("mergeRowsWithSmallestValueByKEGGId")

    metaboName <- as.vector(df$KEGGId)

    replicat <- unique(metaboName[duplicated(metaboName)])

    #colIds <- colnames(df)
    for(x in replicat){

        #find rows ID of this replicated element
        rowReplicat <- which(df$KEGGId == x)


        keggId <- data.frame(KEGGId = df$KEGGId);


        if(ncol(df) > 2){

            df <- df[,-ncol(df)]

        rowReplicatToReplace <- rowReplicat[1]

        # row to remove
        rowReplicatToRemove <- rowReplicat[-1]

        # subset of rows with the replicated value
        rowReplicatSubDF <- df[rowReplicat,]
        rowReplicatSortedSubDF <- apply(rowReplicatSubDF,2,sort)
        rowMinValues <- rowReplicatSortedSubDF[1,]
        }else{

        df <- data.frame("pl" = df[,-ncol(df)])
        rowReplicatToReplace <- rowReplicat[1]

        # row to remove
        rowReplicatToRemove <- rowReplicat[-1]
        # subset of rows with the replicated value
        rowReplicatSubDF <- data.frame("pl" = df[rowReplicat,])


        rowReplicatSortedSubDF <- data.frame("pl"=
                                rowReplicatSubDF[order(rowReplicatSubDF[,1]),])

        rowMinValues <- rowReplicatSortedSubDF[1,]
        }

        # row ID to be replaced with smallest value

        # sort all columns by ascending order (more than 1 column)

        # keep only the first row -> the one with smallest vlaues


        df[rowReplicatToReplace,] <- rowMinValues;

        df <- cbind(df, keggId)
        df <- df[-rowReplicatToRemove,]

     }
    return <- df;

}

merge2DFWithSmallestValue <- function(df1, df2){

   #  print("merge2DFWithSmallestValue")

     for (row in 1:nrow(df)) {
        r <- mergeVectorsLowerValues(df1[row,], df2[row,]);
        rf <- t(data.frame(r));
        finalDF <- rbind(finalDF,rf);
        return <- finalDF;
     }
}

removeDuplicatedColumnDF <- function(df){

    # print("removeDuplicatedColumnDF")

    df[is.na(df)] <- Inf;

    df<- t(df)

    df<- aggregate(df,list(row.names(df)),function(x) x[which.min(abs(x))]);

    row.names(df) <- df[,1]

    df <- data.frame(t(df))

    df <- df[-c(1), ]

  return <- df;

}

removeRowsDistanceAsso <- function(df){

   #  print("removeRowsDistanceAsso")

   # print(df)
    row.names(df) <- c(1:nrow(df))

    if(is.data.frame(df) && nrow(df) > 0){
    df$distance <- as.numeric(as.character(df$distance))
    df$geneKEGGId <- as.character(df$geneKEGGId)
    df$metaboliteKEGGId <- as.character(df$metaboliteKEGGId)

    df <- df[order(df$geneKEGGId, df$metaboliteKEGGId),]

    # Get the smallest value for each pair of gene-metabolite. from
    # the multiplepossibilities of graph nodes

    # Aggregates as problem with NA values, so i will give all NA distance value
    # a value of 100 000 which can't be reach within the biological graph and
    # give back a value of NA Afterwards. It is not very nice fix but works
    # for now

    df$distance[is.na(df$distance)] <- 10000

    dft <- aggregate(df$distance,
             by = list(df$geneKEGGId, df$metaboliteKEGGId),FUN = min)
        #     na.action=na.pass, na.rm=TRUE)

    colnames(dft) <- c("geneKEGGId", "metaboliteKEGGId", "distance")

    dft$distance[dft$distance == 10000] <- NA
     dft1 <- merge(dft, df)

     dft1 <- dft1[!duplicated(dft1[c(1,2,3)]),]
     end.time <- Sys.time()

    } else {
        dft1 <- df;
    }

     return <- dft1;

}

# This function also sorts the data by metabolite or gene ids... we don't want
# that
removeRowsDistanceAll <- function(df){

    # print("removeRowsDistanceAll")

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

    # print("getMaxValIgnoreInfVal")

    infVal <- which(df$distance == Inf)

    temp <- df

    temp[infVal,3] <- -1

    maxVal <- max(temp$distance)

    return <- maxVal;

}


## function to remove duplicated rows of one col and concat data from another
##  col
## exemple
#    df : a   b
#         1   1
#         1   2
#         1   3
#         2   4
#         2   1
#         3   1
#
# df <- concatDfColInfoFromDuplicate(df, a, b)
# output
#    df : a   b
#         1   1 2 3
#         2   4 1
#         3   1

concatDfColInfoFromDuplicate <- function(df, duplicatedCol, concatCol){

    # print("concatDfColInfoFromDuplicate")
    # for this algorithm to work well, sort the dataframe by duplicate col
    df <- df[order(df[,duplicatedCol]),]

    # first remove duplicated rows from df
    df <- df[!duplicated(df),]

    concat <- ""
    duplicat <- ""
    concatColName <- colnames(df)[concatCol]
    duplicatColName <- colnames(df)[duplicatedCol]


    df_new <- data.frame();

    for(x in 1:nrow(df)){

       # condition to get first row
        if(duplicat == ""){
            duplicat <- df[x,duplicatedCol];
            concat <- paste(concat, df[x,concatCol])

       # conditino get inside rows and combine first colinfo info based on
       # duplicated 2 col
        }else if(!(duplicat == as.character(df[x,duplicatedCol]))){

            df_new <- rbind(df_new,
                data.frame(duplicatColName = as.vector(as.character(duplicat)),
                           concatColName = as.vector(concat)))
            concat <- "";
            duplicat <- df[x,duplicatedCol]
            concat <- paste(concat, df[x,concatCol])

        }else{
         concat <- paste(concat, df[x,concatCol])
        }

        # catch last row
        if(x == nrow(df)){
            df_new <- rbind(df_new,
                data.frame(duplicatColName = as.vector(as.character(duplicat)),
                          concatColName = as.vector(concat)))
        }
    }

    return <- df_new;

}

