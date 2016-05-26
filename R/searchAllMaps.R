getDistanceAssoAllMaps <- function(association){
    l <- getListAllHumanMaps();
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

    finalDF <- subset(finalDF , finalDF[2] == TRUE & finalDF[4] == TRUE)
    finalDF <- finalDF[!duplicated(finalDF),]

    return <- finalDF;
}
