
getListAllHumanMaps <- function(){

    url <- AllHumanMapsFile();
    file <- downloadFileByUrl(url);
    file <- gsub("path:","",file);
    file <- gsub(":","",file);

    file <- strsplit(file,"[\\\\]|[^[:print:]]",fixed=FALSE);

    file <- as.vector(unlist(file));
    allMapsIds <- file[grep("hsa", file)];

    return <- allMapsIds;

}


getHumanGeneByMap <- function(pathwayId){

    url <- getHumanGeneByMapUrl(pathwayId);
    file <- downloadFileByUrl(url);

}

getMetaboliteByMap <- function(pathwayId){

    url <- getMetaboliteByMapUrl(pathwayId);
    file <- downloadFileByUrl(url);


}

getHumanGeneByMapUrl <- function(pathwayId){

    url <-  paste("http://rest.kegg.jp/link/genes/", pathwayId, sep = "")
    return <- url;

}

getMetaboliteByMapUrl <- function(pathwayId){

    pathwayId <- gsub("hsa","map",pathwayId);
    url <-  paste("http://rest.kegg.jp/link/cpd/",pathwayId,sep = "")
    return <- url;

}



