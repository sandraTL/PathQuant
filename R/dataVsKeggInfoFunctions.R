
metabolitesInKEGGandInMap <- function(completeMetaboliteDF, mapId){


    allMetaboliteInKEGG <- getAllMetaboliteInKEGG();
    allMetaboliteInMap <- getAllMetaboliteInMap(mapId);

    for(x in 1:length(completeMetaboliteDF[,1])){

        p1 <- paste("^",completeMetaboDF[x,1],sep="")
        p1 <- paste(p1,"$",sep="")

        print(any(grepl(p1, allMetaboliteInMap)));

    print(p1)



    }


}


genesInKEGGandInMap <- function(completeGeneDF, hsaId){

    allgeneInKEGG <- getAllHumanGeneInKEGG();
    allgeneInMap <- getAllGeneInMap(hsaId)
    for(x in 1:length(completeGeneDF[,1])){

        p1 <- paste("^",completeGeneDF[x,1],sep="")
        p1 <- paste(p1,"$",sep="")

        print(any(grepl(p1, allgeneInMap)));

        print(p1)
    }

}
