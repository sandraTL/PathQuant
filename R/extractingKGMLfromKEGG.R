
# A KEGG Function
#
# This function allows you to get a KGML format of a pathway and save
# xml format in data directory
# getPathwayKGML("hsa01100")

getPathwayKGML <- function(pathwayId) {

    adressfile <- toStringAdressfile(pathwayId)
    destfile <- toStringDestfile(pathwayId)

    URL_S <- "http://rest.kegg.jp/get/"
    URL_E <- "/kgml"

    URL_FINAL<- paste(URL_S,pathwayId, URL_E, sep="");

    xmlFile <- RCurl::getURL(URL_FINAL);

    if (xmlFile == "") {
        stop("pathway do not exist in KEGG database", call. = FALSE)
    } else{
        xmlFile <- XML::xmlParse(xmlFile)

        XML::saveXML(xmlFile, file = destfile)
    }
    return <- xmlFile

}

#set path to store downloaded file
toStringDestfile <- function(pathwayId){
    #concatenation of pathwayId to set swdir for the xml

    s2 <-  toString(pathwayId);
    s3 <- ".txt"

    destfile <- paste(s2, s3, sep="");

    return <- destfile;
}

#set path for download
toStringAdressfile <- function(pathwayId){

    s1 <- "rest.kegg.jp/get/";
    s2 <-  toString(pathwayId);
    s3 <- "/kgml"
    s4 <- paste(s1,s2, sep= "");
    adressfile <- paste(s4, s3, sep="");

    return <- adressfile;
}

# see if file was already dowloaded
isFileInDirectory <- function(pathwayId){
    #concatenation of pathwayId to set swdir for the xml

    pathFile <- toStringDestfile(pathwayId)
    bool <- FALSE;
      res <-  tryCatch({
             XML::xmlParse(pathFile)
             bool <- TRUE
         }, error = function(e) {
             bool <- FALSE;

         }, finally = {
             return <- bool;
         })

 return <- res;
}

# set adress to download compound kgml file
toCompoundAdressfile <- function(compoundKeggId){

    s1 <- "rest.kegg.jp/list/";
    s2 <-  toString(compoundKeggId);

    adressfile <- paste(s1,s2, sep= "");

    return <- adressfile;
}

AllHumanMapsFile <- function(){

    url_Address <-  "http://rest.kegg.jp/list/pathway/hsa";
    return <- url_Address;

}

getPathwayKGML <- function(pathwayId) {

    adressfile <- toStringAdressfile(pathwayId)
    destfile <- toStringDestfile(pathwayId)

    URL_S <- "http://rest.kegg.jp/get/"
    URL_E <- "/kgml"

    URL_FINAL<- paste(URL_S,pathwayId, URL_E, sep="");

    xmlFile <- RCurl::getURL(URL_FINAL);

    if (xmlFile == "") {
        stop("pathway do not exist in KEGG database", call. = FALSE)
    } else{
        xmlFile <- XML::xmlParse(xmlFile)

        XML::saveXML(xmlFile, file = destfile)
    }
    return <- xmlFile

}

downloadFileByUrl <- function(url){

    file <- RCurl::getURL(url);


    return <- file;

}
