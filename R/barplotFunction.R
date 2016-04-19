
#' Function calculating shortest distance between a selected gene of a
#' gene-metabolite pairs and every metabolites.
#'
#' Function calculating shortest distance between a selected gene of a
#' gene-metabolite pairs and every metabolites (in a pair or not)
#' on a graph model of KEGG map selected, where nodes are metabolites and
#' reactions are edges.
#'
#' If a gene or a metabolite is present on multiple edges or nodes, then
#' shortest distance are calculated for every combinaison possible and the
#' shortest distance is selected.
#'
#' Output: barplot showing the distribution of the calculated distance, where
#' red bars represent a distance with an associated metabolite and grey bars
#' distance with no associated metabolites.
#'
#' @param pathwayId KEGG Id of selected map
#' @param data is a dataFrame with 2 columns. Where each line
#'        reprensents an associations with the first column as
#'        gene KEGG Ids and the sencond column as metabolite
#'        KEGG Ids.
#' @param metabolite is a dataframe of 1 column with the KEGG Ids of
#'        all measured metabolites.
#' @param gene is the KEGG Id of the selected gene
#' @keywords graph, shortestDistance, KEGG
#'
#' @export
#' @examples distributionGene(metabolismOverviewMapKEGGId,
#' shinAndAlDF, completeMetaboDF, "hsa:1373")

<<<<<<< HEAD
distanceGeneToAllMetabolite <- function(pathwayId, data,metabolite, gene){
=======
distributionGene <- function(pathwayId, data,
                            metabolite, gene){
>>>>>>> refs/remotes/sandraTL/master

    pathwayId <- gsub("hsa:", "hsa", pathwayId)
    mError1 <-"error in metabolite, please input a dataframe of 1
    column with a list of KEGG ids metabolites (ex: C00001)"

    mError2 <-"error in data,
    where colnames(df) <- c(gene,metabolite) frame with
    KEGG ids of genes (ex : hsa:00001) in first
    column and associated KEGG ids metabolites (ex: C00001)
    in second column"

    mError3 <- "error in argument gene, the gene entered doesn't match any
                gene in data"


    #test data
    if(is.data.frame(metabolite) && nrow(metabolite)==0){
        stop(mError1, call. = FALSE);
    }
    if(!is.data.frame(data) ||
       length(data[1,])< 2 ||
       length(data[1,])> 3){

        stop(mError2, call. = FALSE)
    }
    for(row in 1:nrow(metabolite)){

        if(substr(metabolite[row,1],0,1) != "C"
           && length(data[row,1]) != 5)
            stop(mError1, call. = FALSE);

    }
    for(row in 1:nrow(data)){

        if(substr(data[row,1],1,4)!="hsa:")
            stop(mError2, call. = FALSE);
        if(substr(data[row,2],0,1) != "C"
           && length(data[row,2]) != 5)
            stop(mError2, call. = FALSE);
    }

    if(length(data.frame(
        data[data$gene == gene,])[,1]) == 0){

        stop(mError3, call. = FALSE);
    }


    # get all shortest paths from data entry
    shortestsPathsDF <- data.frame(t(getDistanceAll(pathwayId,
                     data[data$gene == gene,],
                     metabolite)));

    associatedMetabo <- data.frame(
                        getAssociatedMetaboByGene(data,gene))

    # adjust gene parameter
    gene1 <- gsub(":", ".", gene);

    # add metabolite row
    shortestsPathsDF[ "metabolites" ] <- rownames(shortestsPathsDF);

    # get a subset of shortestsPathsDF contaning only geneOf interest, gene
    # and metaboltie column

    maxVal <- getMaxValIgnoreInfVal(shortestsPathsDF)

    # get frequency of every value until the maxVal found + Inf val
    frequenceDF <- data.frame(table(factor(shortestsPathsDF[,1],
                                            levels=c(0:maxVal,Inf))))
    # set column names of frequenceDF
    colnames(frequenceDF) <- c("Distance", "Freq")

    # get associations values for all distance pair calculated
    associationsMetaboDF <- getAssociationsDF(shortestsPathsDF,associatedMetabo)

    # bind info into 1 DF
    shortestsPathsDF <- cbind(shortestsPathsDF,
                              Associations = associationsMetaboDF);

    # Association of the distance where there is an association
    associationsfrequencyDF <- getFrequenceAssociationsDF(frequenceDF,
                                                          shortestsPathsDF,
                                                          gene1);

    # Add a column for the coloring of the bar associated with gene to subgraph
    frequenceDF<-cbind(frequenceDF,
                           Associations = associationsfrequencyDF);

    # create barplot
    barplotFunctionGeneToAllMetabo(frequenceDF,gene)

}

# Function to output barplot.

barplotFunctionGeneToAllMetabo <- function(frequenceDF,gene){


    # initiating variable for barplotGraph
    geneCommonName <- getCommonNames(c(gene), "gene")
    numInfValue <- frequenceDF[frequenceDF$Distance == Inf,][,2]
    frequenceDF <- frequenceDF[-length(frequenceDF[,1]),]
    maxDistance <- as.numeric(as.character(frequenceDF[nrow(frequenceDF),][,1]))
    maxFrequency <- max(frequenceDF$Freq, na.rm = TRUE)

    legend_text <- paste("infinite distance count: ",numInfValue, sep = "")

    plot <- ggplot2::ggplot(frequenceDF, ggplot2::aes(
        x = factor(Distance),
        y = Freq,
        fill = Associations
    ),environment = environment())

    plot <- (plot + ggplot2::geom_bar(stat="identity", colour="black",
                                      size=0.5, position="identity",width=1)
             + ggplot2::theme_bw()
             + ggplot2::theme(panel.border = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 text = ggplot2::element_text(size=12),
                 axis.line.x = ggplot2::element_line(color="black"),
                 axis.line.y = ggplot2::element_line(color="black")
                              )
             + ggplot2::xlab("Distance from Gene")
             + ggplot2::ylab("Metabolite count")
             + ggplot2::ggtitle(geneCommonName)
             + ggplot2::coord_fixed(ratio = 1)
             + ggplot2::geom_rect(data = frequenceDF,
                                  ggplot2::aes(xmin = (maxDistance+1 -6),
                                               xmax = maxDistance+1,
                                               ymin = (maxFrequency -1),
                                               ymax = maxFrequency),
                                  fill = "grey80")
             + ggplot2::annotate("text", x = (maxDistance-2),
                                 y = (maxFrequency -0.5),
                                 label = legend_text,
                                 colour = "black",
                                 size=4 )
             + ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(2,4,6,8,10) )

             + ggplot2::scale_fill_manual(values = c("FALSE" ="grey",
                                                     "TRUE" = "red3"),
                                          guide = FALSE)

    );

    filename = paste0(geneCommonName,".png")

    ggplot2::ggsave(filename,width =10, height =7, dpi=300)
    print(plot);

}

getFrequenceAssociationsDF <- function(frequenceDistDF,shortestsPathsDF,gene){

    test = FALSE;
    results <- data.frame();

    #' creation of vector to fill bar colors automatically
    for(row1 in 1:nrow(frequenceDistDF)){
        test = FALSE;
        for(row2 in 1:nrow(shortestsPathsDF)){

            if(frequenceDistDF[row1,"Distance"] == shortestsPathsDF[row2,gene]){
                if(shortestsPathsDF[row2,"Associations"] == TRUE){
                    test<- TRUE;

                    break;
                }
            }else test<- FALSE;
        }

        results <- rbind(results,test);
        colnames(results) <- c("Associations")

    }

    return <- results;
}


getAssociationsDF <- function(assoDistanceDF, associatedMetaboDF){

    #' initiation of values
    test = FALSE;
    results <- data.frame();

    #' buiding boolean associations DF
    for(row1 in 1:nrow(assoDistanceDF)){
        test = FALSE;
        for(row2 in 1:nrow(associatedMetaboDF)){

            if(assoDistanceDF[row1,"metabolites"] == associatedMetaboDF[row2,]){
                test<- TRUE;

                break;
            }else test<- FALSE;
        }

        results <- rbind(results,test);
        colnames(results) <- c("Associations")

    }
    return <- results;
}

getAssociatedMetaboByGene <- function(data, gene){

    selectedRows <- data[data$gene == gene,];
    associatedMetabo <- selectedRows[,2];

    return <- associatedMetabo;

}


# function to modify the frequency data.frame by adding more distance that are
# then the one already defined, ulterior use.

barplot_adjustMaximalDistance <- function(maximumDistance, frequenciesDF,
                                          maxDistFrequenciesDF){
        levels(frequenciesDF$Distance) <- c(levels(frequenciesDF$Distance),
                                            c(1:25))
        if(maxDistFrequenciesDF < maximumDistance){
            for(i in (maxDistFrequenciesDF+1):25){

                newRow <- (data.frame("Distance" = i,
                                       "Freq" = 0,
                                       "Associations" = FALSE))
                frequenciesDF <- rbind(frequenciesDF, newRow)

            }
        }
    return <- frequenciesDF

}
