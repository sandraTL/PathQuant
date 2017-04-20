
# This function creates a data frame which contains the duplicated compounds
# in the initial graph. We need this because some reaction has a subtrate or
# a product that has a metabolite that we find more than one in the graph but
# we don't know which node we should use (reaction type ortholog). We prefere
# not to use those reactions

getDataFrameOfDuplicateCompounds <- function(compoundDataFrame){

  #  print("getDataFrameOfDuplicateCompounds")
    # Order and remove duplicates in dataFrame
    compoundDataFrame <- compoundDataFrame[order(compoundDataFrame[,2]), ];
    compoundDataFrame <- compoundDataFrame[duplicated(compoundDataFrame[,2]), ];

    # dataframe with only duplicated metabolites

    duplicatedCompoundDataFrame <- data.frame(compoundDataFrame[2]);

    return <- duplicatedCompoundDataFrame;
}


correctKeggIdString <- function(nodeDF){

   # print("correctKeggIdString")

    nodeDF <- lapply(nodeDF,
                     function (x) gsub("cpd:","",x))
    return <- nodeDF;

}

