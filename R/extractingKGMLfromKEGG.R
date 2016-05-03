
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

# # get list of all kinds metabolite in KEGG.
# getAllMetaboliteInKEGG <- function(){
#
#
#     metaboliteList <- as.vector(names(KEGGREST::keggList("cpd")))
#
#     metaboliteList <- gsub("cpd:", "", metaboliteList)
#    # adressfile <- "rest.kegg.jp/list/cpd";
#     return <- metaboliteList;
# }
#
# # get list of all human genes in KEGG
# getAllHumanGeneInKEGG<- function(){
#
#     geneList <- as.vector(names(KEGGREST::keggList("hsa")))
#     #adressfile <- "rest.kegg.jp/list/hsa";
#     return <- geneList;
# }
#
# # get a list of all metbolites on a specific map
# getAllMetaboliteInMap <- function(mapId){
#
#     metaboliteList <- as.vector(KEGGREST::keggLink("cpd",mapId))
#     metaboliteList <- gsub("cpd:", "", metaboliteList)
#
#     return <- metaboliteList;
# }
#
# # get list of all genes in a specific map
# getAllGeneInMap<- function(hsaId){
#
#     geneList <- as.vector(KEGGREST::keggLink("genes",hsaId))
#
#     return <- geneList;
# }



