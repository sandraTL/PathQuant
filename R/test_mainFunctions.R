


mError <<- "error in association, KEGG ids of genes (ex : hsa:00001) in first column and associated KEGG ids metabolites (ex: C00001) in second column";
mError1 <<- "error in metabolite, please input a dataframe of 1 column with a list of KEGG ids metabolites (ex: C00001)";
mError2 <<- "error in gene, please input a dataframe of 1 column with a list of KEGG ids gene (ex: hsa:00001)";
mError3 <<- "error in argument gene, the gene entered doesn't match any gene in association";

test_getDistanceAsso <- function(pathwayId, association){


    if(isFileInDirectory(pathwayId) == FALSE){
        file <-  getPathwayKGML(pathwayId)
    }
    if(is.data.frame(association) && nrow(association)==0){
        stop(mError,call. = FALSE )
    }
    if(is.data.frame(association) &&
       !ncol(association) == 2){
        stop(mError,call. = FALSE )
    }
    for(row in 1:nrow(association)){

        if(substr(association[row,1],1,4)!="hsa:")
            stop(mError, call. = FALSE);
        if(substr(association[row,2],0,1) != "C"
           && length(association[row,2]) != 5)
            stop(mError, call. = FALSE);
    }
}

test_getDistanceAssoPerm <- function(pathwayId, data){



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

    # test metabolite
    if(is.data.frame(metabolite) && nrow(metabolite) == 0){
        stop(mError1, call. = FALSE);
    }
    # test gene
    if(is.data.frame(gene) && nrow(gene) == 0){
        stop(mError2, call. = FALSE);
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

test_distributionGene <- function(pathwayId,association, metabolite, gene){

    #test association
    if(is.data.frame(metabolite) && nrow(metabolite)==0){
        stop(mError1, call. = FALSE);
    }
    if(!is.data.frame(association) ||
       length(association[1,])< 2 ||
       length(association[1,])> 3){

        stop(mError, call. = FALSE)
    }
    for(row in 1:nrow(metabolite)){

        if(substr(metabolite[row,1],0,1) != "C"
           && length(association[row,1]) != 5)
            stop(mError1, call. = FALSE);

    }
    for(row in 1:nrow(association)){

        if(substr(association[row,1],1,4)!="hsa:")
            stop(mError, call. = FALSE);
        if(substr(association[row,2],0,1) != "C"
           && length(association[row,2]) != 5)
            stop(mError, call. = FALSE);
    }

    if(length(data.frame(
        association[association$gene == gene,])[,1]) == 0){

        stop(mError3, call. = FALSE);
    }


}

test_heatmap <- function(pathwayId, association){

    if(length(association) == 0){
        stop(mError, call. = FALSE);
    }
    if(!is.data.frame(association) ||
       length(association[1,])< 2 ||
       length(association[1,])> 3){
        stop(mError, call. = FALSE);
    }
    for(row in 1:nrow(association)){

        if(substr(association[row,1],1,4)!="hsa:")
            stop(mError, call. = FALSE);
        if(substr(association[row,2],0,1) != "C"
           && length(association[row,2]) != 5)
            print(substr(association[row,2],0,1))
            stop(mError, call. = FALSE);
    }

}

test_permutationTest <- function(pathwayId, association, gene, metabolite,
                                 permutation, output){
    mError4 <- "the number of permutation you entered is wrong,
    Enter a positive integer";
    mError5 <- "output parameter should be 'medians',
    'histogram' or 'pvalue'"


    if(length(association) == 0){
        stop(mError, call. = FALSE);
    }
    if(!is.data.frame(association) ||
       length(association[1,])< 2 ||
       length(association[1,])> 3){
        stop(mError, call. = FALSE);
    }
    for(row in 1:nrow(association)){

        if(substr(association[row,1],1,4)!="hsa:")
            stop(mError, call. = FALSE);
        if(substr(association[row,2],0,1) != "C"
           && length(association[row,2]) != 5)
            stop(mError, call. = FALSE);
    }
    # test metabolite
    if(is.data.frame(metabolite) && nrow(metabolite) == 0){

        stop(mError1, call. = FALSE);
    }
    # test gene
    if(is.data.frame(gene) && nrow(gene) == 0){
        stop(mError2, call. = FALSE);
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

    if(!is.integer(permutation) &&
       permutation <= 0){

        stop(mError4, call. = FALSE);
    }

    if(!(output == "medians" ||
       output == "histogram" ||
       output == "pvalue")){

        stop(mError5, call. = FALSE);
    }

}
