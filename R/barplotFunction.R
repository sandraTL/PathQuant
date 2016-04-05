

#' function that output a barplot graph related to one specific gene with all
#' the shortest distances from that gene to all metabolites
#'
#' for param data:
#'      gene = KEGGid of gene hsa:...
#'      metabolites : KEGGid of metabolites C....
#' for param pathwayId : KEGG id of pathways without ':' ex: hsa01100
#' for param gene : is a gene in data ex: hsa:8801
#' @param KEGGPathwayId, data(gene, metabolites),allMeasuredMetaboliteesDF, gene
#' @keywords  Graph, barplot, shortestDistance
#' @export
#' @examples barplotFunction(hsa01100, associatedGeneMetaDF,
#'  completeMetaboliteDF, gene)

distanceGeneToAllMetabolite <- function(pathwayId, associatedGeneMetaDF,
                            completeMetaboliteDF, gene){

    mError1 <-"error in completeMetaboliteDF, please input a dataframe of 1
    column with a list of KEGG ids metabolites (ex: C00001)"

    mError2 <-"error in associatedGeneMetaDF,
    where colnames(df) <- c(gene,metabolite) frame with
    KEGG ids of genes (ex : hsa:00001) in first
    column and associated KEGG ids metabolites (ex: C00001)
    in second column"

    mError3 <- "error in argument gene, the gene entered doesn't match any
                gene in associatedGeneMetaDF"


    #test associatedGeneMetaDF
    if(is.data.frame(completeMetaboliteDF) && nrow(completeMetaboliteDF)==0){
        stop(mError1, call. = FALSE);
    }
    if(!is.data.frame(associatedGeneMetaDF) ||
       length(associatedGeneMetaDF[1,])< 2 ||
       length(associatedGeneMetaDF[1,])> 3){

        stop(mError2, call. = FALSE)
    }
    for(row in 1:nrow(completeMetaboliteDF)){

        if(substr(completeMetaboliteDF[row,1],0,1) != "C"
           && length(associatedGeneMetaDF[row,1]) != 5)
            stop(mError1, call. = FALSE);

    }
    for(row in 1:nrow(associatedGeneMetaDF)){

        if(substr(associatedGeneMetaDF[row,1],1,4)!="hsa:")
            stop(mError2, call. = FALSE);
        if(substr(associatedGeneMetaDF[row,2],0,1) != "C"
           && length(associatedGeneMetaDF[row,2]) != 5)
            stop(mError2, call. = FALSE);
    }

    if(length(data.frame(
        associatedGeneMetaDF[associatedGeneMetaDF$gene == gene,])[,1]) == 0){

        stop(mError3, call. = FALSE);
    }


    #' get all shortest paths from data entry
    shortestsPathsDF <- data.frame(t(getDistanceAll(pathwayId,
                     associatedGeneMetaDF[associatedGeneMetaDF$gene == gene,],
                     completeMetaboliteDF)));

    associatedMetabo <- data.frame(
                        getAssociatedMetaboByGene(associatedGeneMetaDF,gene))

    #'adjust gene parameter
    gene1 <- gsub(":", ".", gene);

    #' add metabolite row
    shortestsPathsDF[ "metabolites" ] <- rownames(shortestsPathsDF);

    #' get a subset of shortestsPathsDF contaning only geneOf interest, gene
    #' and metaboltie column

    maxVal <- getMaxValIgnoreInfVal(shortestsPathsDF)

    #' get frequency of every value until the maxVal found + Inf val
    frequenceDF <- data.frame(table(factor(shortestsPathsDF[,1],
                                            levels=c(0:maxVal,Inf))))
    #' set column names of frequenceDF
    colnames(frequenceDF) <- c("Distance", "Freq")

    #' get associations values for all distance pair calculated
    associationsMetaboDF <- getAssociationsDF(shortestsPathsDF,associatedMetabo)

    #' bind info into 1 DF
    shortestsPathsDF <- cbind(shortestsPathsDF,
                              Associations = associationsMetaboDF);

    #' Association of the distance where there is an association
    associationsfrequencyDF <- getFrequenceAssociationsDF(frequenceDF,
                                                          shortestsPathsDF,
                                                          gene1);

    # Add a column for the coloring of the bar associated with gene to subgraph
    frequenceDF<-cbind(frequenceDF,
                           Associations = associationsfrequencyDF);

    # create barplot
    barplotFunctionGeneToAllMetabo(frequenceDF,gene)

}

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
                              text = ggplot2::element_text(size=12,family="Arial"),
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
                                 size=4,
                                 family="Arial" )
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

#' function to modify the frequency data.frame by adding more distance that are
#' then the one already defined

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
