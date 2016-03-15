
# for supplement files...
metabolitesInKEGGandInMap <- function(completeMetaboliteDF, mapId){


    allMetaboliteInKEGG <- getAllMetaboliteInKEGG();
    allMetaboliteInMap <- getAllMetaboliteInMap(mapId);

     numberOfMetabolites(allMetaboliteInMap, shinAndAlDF[,2])

}

# for supplement files..
genesInKEGGandInMap <- function(completeGeneDF, hsaId){

    allgeneInKEGG <- getAllHumanGeneInKEGG();
    allgeneInMap <- getAllGeneInMap(hsaId)
    for(x in 1:length(completeGeneDF[,1])){

        p1 <- paste("^",completeGeneDF[x,1],sep="")
        p1 <- paste(p1,"$",sep="")

    }

}

associationsInKegg <- function(associationDF){




}
