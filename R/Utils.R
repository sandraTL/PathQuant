
#' Inporting data from an XSLX file, in the context of this package, one of the
#' data to import is the association file.
#'
#' @param path to XSLX file with 2 columns where each line reprensents an
#'        association. First column are genes and the sencond column are
#'        metabolites. Only use KEGG Ids.
#' @keywords data, input, XSLX, data frame, associations.
#' @export
#' @examples inputXLSXtoDF("file_path/fileName.xlsx")

inputXLSXtoDF <- function(association){

    input_data = gdata::read.xls(association);

    return <- input_data;

}

#' Export results data of a data frame into a text files, for which each lines
#' separates the data from every data frame colum by tabulation. This text file
#' can be easily opened with excel.

#' @param a dataframe.
#' @keywords data, export, txt, dataframe.
#' @export
#' @examples exportDFtoTxt(dataFrame, exportFileName)

exportDFtoTxt <- function(export_data, exportNameFile){



    exportNameFile = paste(exportNameFile, ".txt", sep = "");

    return <- write.table(export_data, exportNameFile, sep="\t")

}


#' get the shortest path between on specific association between one gene and
#' one metabolite

#' @param pathwayId KEGG Id of selected pathway.
#' @param gene KEGG Id.
#' @param Metabolite KEGG Id.
#' @param commonNames get KEGG's Common Names of the KEGG Id provided in the
#'        path results.
#' @keywords association, path.
#' @export
#' @examples getPath("hsa:01100","hsa:1373", "C00300", T)

getPath <- function(pathwayId, gene, metabolite, commonNames = F){

    asso <- data.frame(gene = as.vector(gene),
                       metabolite = as.vector(metabolite));

      res <- getDistanceAsso(pathwayId, asso,F,F);
      path <- res["path"]

      if(commonNames == T){

          path <-getPathCommonNames(path);
      }

     return <- path;

}






