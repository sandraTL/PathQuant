
#' This fucntion creates a data frame which contains the duplicated compounds
#' in the initial graph. We need this because some reaction has a subtrate or
#' a product that has a metabolite that we find more than one in the graph but
#' we don't know which node we should use (reaction type ortholog). We prefer
#' not to use those reactions
#' #### ortholog
#' @param compounds data frame
#' @keywords  compounds
#' @examples getDataFrameOfDuplicateCompounds(compoundDataFrame)

getDataFrameOfDuplicateCompounds <- function(compoundDataFrame){

    # Order and remove duplicates in dataFrame
    # Gets all the duplicated compound but if duplicated more than once,
    # it only takes out 1 duplicate (not optimal).... Would apply unique()
    # on the result but doesnt work for some reason... Got to find another
    # solution but it is not really important since the duplicate list isnt big.
    compoundDataFrame <- compoundDataFrame[order(compoundDataFrame[,2]), ];
    compoundDataFrame <- compoundDataFrame[duplicated(compoundDataFrame[,2]), ];

    # dataframe with only duplicated metabolites



    duplicatedCompoundDataFrame <- data.frame(compoundDataFrame[2]);

    return <- duplicatedCompoundDataFrame;
}


correctKeggIdString <- function(nodeDF){

    nodeDF <- lapply(nodeDF,
                     function (x) gsub("cpd:","",x))

    return <- nodeDF;

}

