# isGeneInMap, metaboliteCommonName, metaboliteKEGGId,
# isMetaboliteInMap, distance
#
#  association Dataframe with 2 columns, where each line reprensents an
#        association. First column are the genes and the sencond column are the
#        metabolites. Only use KEGG Ids.
# ordered [option] ascendent ordering of distance
# commonNames get KEGG's Common Names of the KEGG Id in the
#        results.
# graph, shortestDistance, KEGG

# getDistanceAssoAllMaps(shin,ordered = T, commonNames = F)


getDistanceAssoAllMaps <- function(association, ordered =FALSE,
                                   commonNames = TRUE){


     l <- getListAllHumanMaps();

    #bug with map hsa00250....
    l[22] <- "hsa00250";
    finalDF <- data.frame();
    for(i in 1:length(l)){

      r <-getDistanceAsso(l[i], association, F, commonNames=FALSE);


      if(is.data.frame(r)){

         if(nrow(r) >0){
             finalDF <- rbind(finalDF,r )
         }
      }

    }

    finalDF <- removeRowsDistanceAsso(finalDF)

    # Keeps associations only when the gene and the metabolite are present in
    # the same map.
    #     if(mapped == 1){
    #        finalDF <- subset(finalDF , finalDF[3] == TRUE & finalDF[6] == TRUE)
    #        finalDF <- finalDF[!duplicated(finalDF),]
    #     }else if(mapped ==2){
    #
    #         finalDF <- subset(finalDF , finalDF[3] == FALSE | finalDF[6] == FALSE)
    #         finalDF <- finalDF[!duplicated(finalDF),]
    # }



    # finalDF <- subset(finalDF , finalDF[2] == TRUE & finalDF[4] == TRUE)
    # finalDF <- finalDF[!duplicated(finalDF),]

    for(i in 1:nrow(finalDF)){
        if(finalDF[i,5] =="NaN"){
            finalDF[i,6] = NaN;
        }
    }

    if(ordered == TRUE){
        finalDF <- finalDF[order(as.numeric(as.character(finalDF[,5])), decreasing = FALSE, na.last = TRUE), ]
    }

    if(commonNames == TRUE){

         geneCommonName <- as.vector(unlist(getCommonNames(as.vector
                                    (unlist(finalDF[,1])),"gene")))
         metaboliteCommonName <- as.vector(unlist(getCommonNames(as.vector
                                    (unlist(finalDF[,3])), "metabolite")))


         finalDF1 <- data.frame(
             "geneCommonName" = geneCommonName,
             "geneKEGGId" = finalDF[,1],
             "isGeneInMap" = finalDF[,2],
             "metaboliteCommonName" = metaboliteCommonName,
             "metaboliteKEGGId" = finalDF[,3],
             "isMetaboliteInMap" = finalDF[,4],
             "distance" = finalDF[,5],
             "pathwayId" = finalDF[,6],
             "path" = finalDF[,7]);
         rownames(finalDF1) <- c(1:nrow(finalDF1))
    }else{
         finalDF1 <- finalDF
     }

    return <- finalDF1;
}



getUnmappedAsso <- function(aM){


    # aM <- getDistanceAssoAllMaps(shinAndAlDF,2)

    aMAsso <- aM[,c(1,2,3,4)]
    aMAsso <- aMAsso[!duplicated(aMAsso),]
    print(aMAsso)

    # v <- getAssociationForHeatmap(aMAsso,shinAndAlDF)

    newAsso <- aMAsso
    geneCommonName <- getCommonNames(as.vector(unlist(newAsso[,1])), "gene")

    geneCommonName <- as.vector(unlist(geneCommonName))

    metaboliteCommonName <- getCommonNames(as.vector(unlist(newAsso[,3])),
                                           "metabolite")
    metaboliteCommonName <- as.vector(unlist(metaboliteCommonName))
    newAsso <- cbind(newAsso, "geneCommonName" = geneCommonName,
                     "metaboliteCommonNames" = metaboliteCommonName)

    newAsso <- getBrites(newAsso)

    ## Ouput in xlsx file

    WriteXLS::WriteXLS(newAsso, "AssoNotInKEGG.xlsx")

    return <- newAsso;
}




