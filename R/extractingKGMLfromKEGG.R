
#' A KEGG Function
#'
#' This function allows you to get a KGML format of a pathway and save
#' xml format in data directory
#' @param this function parses all data from kegg
#' @keywords  kegg
#' @examples
#' getPathwayKGML("hsa01100")

getPathwayKGML <- function(pathwayId) {

    adressfile <- toStringAdressfile(pathwayId);
    destfile <- toStringDestfile(pathwayId);

    op <- options(warn=2)
    file <- tryCatch( download.file(adressfile, destfile, "libcurl")
                      ,error=function(e) e,
                     warning=function(w) w)

    if(is(file,"warning")){

        if(file[1]$message == "download had nonzero exit status"){

            stop("pathway do not exist in KEGG database",call. = FALSE )
        }
    }

return <-file;
}

#set path to store downloaded file
toStringDestfile <- function(pathwayId){
    #concatenation of pathwayId to set swdir for the xml
    s1 <- "~/";
    setwd(s1);
    s2 <-  toString(pathwayId);
    s3 <- ".xml"
    s4 <- paste(s1,s2, sep= "");
    destfile <- paste(s4, s3, sep="");

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

    bool = FALSE;
    files <- list.files("~/")
    s2 <-  toString(pathwayId);
    s3 <- ".xml"

    file <- paste(s2, s3, sep="");
    m <- match(file, files, nomatch = NA, incomparables = NULL)

    if(is.na(m) == FALSE){
        bool = TRUE;
    }

    return <- bool;
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



