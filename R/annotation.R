

# Enzyme extracted from http://rest.kegg.jp/get/br:hsa01000 parsing
# using gene KEGG id number
#
# Brite extracted from http://rest.kegg.jp/link/br/hsa parsing
# Classification made by brite definition

#' Gene annotation classification of associations
#'
#' Gene annotation involve EC number, kegg brite annotation
#' High level classification made by involved brite
#' trnasporter, enzyme or protein
#'
#' @param association Dataframe with 2 columns, where each line reprensents an
#'        association. First column are the genes and the sencond column are the
#'        metabolites. Only use KEGG Ids.
#' @keywords gene annotation, classification, KEGG
#' @export
#' @examples annotateAssociationData(associations)

annotateAssociationData <- function(associations){

   # print("annotateAssociationData")
     # get brite ids for each gene
     # can be upgraded by remebering wich gene already searched

     briteDF <- data.frame("brite" =
                               as.vector(getBriteListWithNa(associations)));

     # form data frame with brite info for next steps
     df <- data.frame("gene" = as.vector(associations[,1]),
                      "brite" = briteDF,
                      "metabolite" = as.vector(associations[,2]))

     # get EC number of gene has br:01000
     ecDF <- data.frame("ec" = as.vector(getEcListWithNa(df)));

     # get classification for each gene
     classDF <- get.gene.product(df)

     # form final data frame
     df <- data.frame( "gene" = as.vector(associations[,1]),
                       "ec" =  ecDF,
                       "brite" = briteDF,
                       "classification" = classDF,
                       "metabolite" = as.vector(associations[,2]))

     return <- df;

}


geneClassificationForMappingMetabolicGenes <- function(associations){

   print("geneClassificationForMappingMetabolicGenes")
  #first step Get Brite KEGG Ids associated with the genes in the associations
  briteDF <- data.frame("brite" =
                            as.vector(getBriteListWithNa(associations)));

  geneType <- list()

  #find brite hsa01000 : enzyme brite
  geneType <- lapply(briteDF[,1], function(x){
      grepl(paste(c("br:hsa01000","^*br:hsa01000$"), collapse = "|"),x)
  })

  geneType <- data.frame("brite" = do.call(rbind, geneType))
  geneTypeName <- data.frame()

  # identify the associations with enzymatic genes
  for(i in 1:nrow(geneType)){
     if(geneType[i,1] == TRUE){ geneTypeName[i,1] <- "metabolic gene"
     }else geneTypeName[i,1] <- "not metabolic gene"
  }

  # rebuild data frame with metabolic genes info
  associations <- cbind("gene" = as.vector(associations[,1]),
                        "briteId" = as.vector(briteDF),
                        "geneType" = as.vector(geneTypeName),
                        "metabolite" = as.vector(associations[,2]));

  return <- briteDF

}




