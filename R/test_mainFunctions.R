test_getDistanceAsso <- function(pathwayId, data){
    mError <- "error in data,
    where colnames(df) <- c(gene,metabolite) frame with
    KEGG ids of genes (ex : hsa:00001) in first
    column and associated KEGG ids metabolites (ex: C00001)
    in second column"

    if(isFileInDirectory(pathwayId) == FALSE){
        file <-  getPathwayKGML(pathwayId)
    }
    if(is.data.frame(data) && nrow(data)==0){
        stop(mError,call. = FALSE )
    }
    if(is.data.frame(data) &&
       !ncol(data) == 2){
        stop(mError,call. = FALSE )
    }
    for(row in 1:nrow(data)){

        if(substr(data[row,1],1,4)!="hsa:")
            stop(mError, call. = FALSE);
        if(substr(data[row,2],0,1) != "C"
           && length(data[row,2]) != 5)
            stop(mError, call. = FALSE);
    }
}

test_getDistanceAssoPerm <- function(pathwayId, data){

    mError <- "error in data,
    where colnames(df) <- c(gene,metabolite) frame with
    KEGG ids of genes (ex : hsa:00001) in first
    column and associated KEGG ids metabolites (ex: C00001)
    in second column"

    if(isFileInDirectory(pathwayId) == FALSE){
        file <-  getPathwayKGML(pathwayId)
    }
    if(is.data.frame(data) && nrow(data)==0){
        stop(mError,call. = FALSE )
    }
    if(is.data.frame(data) &&
       !ncol(data) == 2){
        stop(mError,call. = FALSE )
    }
    for(row in 1:nrow(data)){

        if(substr(data[row,1],1,4)!="hsa:")
            stop(mError, call. = FALSE);
        if(substr(data[row,2],0,1) != "C"
           && length(data[row,2]) != 5)
            stop(mError, call. = FALSE);

    }

}

test_getDistanceAll <- function(pathwayId, gene, metabolite){

    mError1 <- "error in metabolite, please input a dataframe of 1 column
    with a list of KEGG ids metabolites (ex: C00001)"

    mError2 <- "error in data,
    where colnames(df) <- c(gene,metabolite) frame with
    KEGG ids of genes (ex : hsa:00001) in first
    column and associated KEGG ids metabolites (ex: C00001)
    in second column"

    # test metabolite
    if(is.data.frame(metabolite) && nrow(metabolite) == 0){
        stop(mError1, call. = FALSE);
    }
    # test gene
    if(is.data.frame(gene) && nrow(gene) == 0){
        stop(mError1, call. = FALSE);
    }

    for(row in 1:nrow(metabolite)){

        if(substr(metabolite[row,1],0,1) != "C"
           && length(metabolite[row,1]) != 5)
            stop(mError1, call. = FALSE);

    }

    for(row in 1:nrow(gene)){

        if(substr(gene[row,1],1,4)!="hsa:")
            stop(mError2, call. = FALSE);

    }

}
