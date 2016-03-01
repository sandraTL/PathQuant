
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
    file <- tryCatch( download.file(adressfile, destfile,
                                    quiet = TRUE,method = "curl"),error=function(e) e,
                     warning=function(w) w)

    if(is(file,"warning")){

        if(file[1]$message == "download had nonzero exit status"){

            stop("pathway doesn't exist in KEGG database",call. = FALSE )
        }
    }

return <-file;
}

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


toStringAdressfile <- function(pathwayId){

    s1 <- "rest.kegg.jp/get/";
    s2 <-  toString(pathwayId);
    s3 <- "/kgml"
    s4 <- paste(s1,s2, sep= "");
    adressfile <- paste(s4, s3, sep="");

    return <- adressfile;
}

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

toCompoundAdressfile <- function(compoundKeggId){

    s1 <- "rest.kegg.jp/list/";
    s2 <-  toString(compoundKeggId);

    adressfile <- paste(s1,s2, sep= "");

    return <- adressfile;
}

getAllMetaboliteInKEGG <- function(){


    metaboliteList <- as.vector(names(KEGGREST::keggList("cpd")))

    metaboliteList <- gsub("cpd:", "", metaboliteList)
   # adressfile <- "rest.kegg.jp/list/cpd";
    return <- metaboliteList;
}

getAllHumanGeneInKEGG<- function(){

    geneList <- as.vector(names(KEGGREST::keggList("hsa")))
    #adressfile <- "rest.kegg.jp/list/hsa";
    return <- geneList;
}

getAllMetaboliteInMap <- function(mapId){

    metaboliteList <- as.vector(KEGGREST::keggLink("cpd",mapId))
    metaboliteList <- gsub("cpd:", "", metaboliteList)

    return <- metaboliteList;
}

getAllGeneInMap<- function(hsaId){

    geneList <- as.vector(KEGGREST::keggLink("genes",hsaId))

    return <- geneList;
}


