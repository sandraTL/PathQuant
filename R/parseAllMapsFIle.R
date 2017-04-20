
getListAllHumanMaps <- function(){

   # print("getListAllHumanMaps")

    url <- AllHumanMapsFile();
    file <- downloadFileByUrl(url);
    file <- gsub("path:","",file);
    file <- gsub(":","",file);

    file <- strsplit(file,"[\\\\]|[^[:print:]]",fixed=FALSE);

    file <- as.vector(unlist(file));
    allMapsIds <- file[grep("hsa", file)];
    #print(allMapsIds)
    return <- allMapsIds;

}


getHumanGeneByMap <- function(pathwayId){

   # print("getHumanGeneByMap")

    url <- getHumanGeneByMapUrl(pathwayId);
    file <- downloadFileByUrl(url);

}

getMetaboliteByMap <- function(pathwayId){

   # print("getMetaboliteByMap")

    url <- getMetaboliteByMapUrl(pathwayId);
    file <- downloadFileByUrl(url);


}

getHumanGeneByMapUrl <- function(pathwayId){

   # print("getHumanGeneByMapUrl")

   url <-  paste("http://rest.kegg.jp/link/genes/", pathwayId, sep = "")
   return <- url;

}

getMetaboliteByMapUrl <- function(pathwayId){

   # print("getMetaboliteByMapUrl")

    pathwayId <- gsub("hsa","map",pathwayId);
    url <-  paste("http://rest.kegg.jp/link/cpd/",pathwayId,sep = "")
    return <- url;

}



