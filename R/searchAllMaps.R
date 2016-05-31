getDistanceAssoAllMaps <- function(association, mapped=c(1,2)){
    l <- getListAllHumanMaps();

    #bug with map hsa00250....
    l[22] <- "hsa00250";
    finalDF <- data.frame();
    for(i in 1:length(l)){

      r <-getDistanceAsso(l[i], association);

      if(is.data.frame(r)){

         if(nrow(r) >0){
             finalDF <- rbind(finalDF,r )
         }
      }

    }
   # Keeps associations only when the gene and the metabolite are present in
   # the same map.
    if(mapped == 1){
 finalDF <- subset(finalDF , finalDF[2] == TRUE & finalDF[4] == TRUE)
     finalDF <- finalDF[!duplicated(finalDF),]
    }else if(mapped ==2){

        finalDF <- subset(finalDF , finalDF[2] == FALSE | finalDF[4] == FALSE)
        finalDF <- finalDF[!duplicated(finalDF),]
}

    return <- finalDF;
}



getUnmappedAsso <- function(aM){

  #  aM <- getDistanceAssoAllMaps(shinAndAlDF,2)

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
    print(newAsso)
    WriteXLS::WriteXLS(newAsso, "AssoNotInKEGG.xlsx")

    return <- newAsso;
}




