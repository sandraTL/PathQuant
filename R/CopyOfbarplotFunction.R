

#' function that output a barplot graph related to one specific gene with all
#' the shortest distances from that gene to all metabolites
#'
#' for param data:
#'      gene = KEGGid of gene hsa:...
#'      metabolites : KEGGid of metabolites C....
#' for param pathwayId : KEGG id of pathways without ':' ex: hsa01100
#' for param gene : is a gene in data ex: hsa:8801
#' @param pathwayId, data(gene, metabolites), gene
#' @keywords  Graph, barplot, shortestDistance
#' @export
#' @examples barplotFunction(hsa01100, associatedGeneMetaDF,
#'  completeMetaboliteDF, gene)

distanceGeneToAllMetabolite <- function(pathwayId, associatedGeneMetaDF,
                            completeMetaboliteDF, gene){
    ############################################################
    #' Serious need to refactor this function in multiple ones
    #' It is way to long
    ############################################################

    geneCommonName <- getCommonNames(c(gene), "gene")

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
    numberUniqueMetabolites <- length( shortestsPathsDF[,1])

    infVal <- which(shortestsPathsDF[,1] == Inf)
    temp <- shortestsPathsDF
    temp[infVal,] <- -1
    maxVal <- max(temp[,1])

    #' get frequency of every value until the maxVal found + Inf val
    frequenceDistDF <- data.frame(table(factor(shortestsPathsDF[,1],
                                            levels=c(0:maxVal,Inf))))
    colnames(frequenceDistDF) <- c("Var1", "Freq")

    #' initiation of values
    test = FALSE;
    results <- data.frame();

    #' creation of vector to fill bar colors automatically
    for(row1 in 1:nrow(shortestsPathsDF)){
        test = FALSE;
        for(row2 in 1:nrow(associatedMetabo)){

            if(shortestsPathsDF[row1,"metabolites"] == associatedMetabo[row2,]){
                test<- TRUE;

                break;
            }else test<- FALSE;
        }

        results <- rbind(results,test);
        colnames(results) <- c("Associations")
        return <- results;
    }

    shortestsPathsDF.plot<-shortestsPathsDF;
    shortestsPathsDF.plot<-cbind(shortestsPathsDF.plot, Associations = results);

    #' initiation of values
    test = FALSE;
    results1 <- data.frame();

    #' creation of vector to fill bar colors automatically
    for(row1 in 1:nrow(frequenceDistDF)){
        test = FALSE;
        for(row2 in 1:nrow(shortestsPathsDF.plot)){

            if(frequenceDistDF[row1,"Var1"] == shortestsPathsDF.plot[row2,gene1]){
                if(shortestsPathsDF.plot[row2,"Associations"] == TRUE){
                    test<- TRUE;

                    break;
                }
            }else test<- FALSE;
        }

        results1 <- rbind(results1,test);
        colnames(results1) <- c("Associations")
        return <- results1;
    }

    # Add a column for the coloring of the bar associated with gene to subgraph
    frequencies<-frequenceDistDF;
    frequencies<-cbind(frequencies, Associations = results1);

    # create barplot
    barplotFunctionGeneToAllMetabo(frequencies,geneCommonName)

}

barplotFunctionGeneToAllMetabo <- function(frequencies,geneCommonName){

    # create barplot
    numInfValue <- frequencies[frequencies$Var1 == Inf,][,2]
    frequencies <- frequencies[-length(frequencies[,1]),]
    maxDistance <- as.numeric(as.character(frequencies[nrow(frequencies),][,1]))
    maxFrequency <- max(frequencies$Freq, na.rm = TRUE)

    legend_text <- paste("infinite distance count: ",numInfValue, sep = "")

    plot <- ggplot2::ggplot(frequencies, ggplot2::aes(
        x = factor(Var1),
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
                              axis.line = ggplot2::element_line(colour = "black"))
             + ggplot2::xlab("Distance from Gene")
             + ggplot2::ylab("Metabolite count")
             + ggplot2::ggtitle(geneCommonName)
             + ggplot2::coord_fixed(ratio = 1)
             + ggplot2::geom_rect(data = frequencies,
                                  ggplot2::aes(xmin = (maxDistance+1 -8),
                                               xmax = maxDistance+1,
                                               ymin = (maxFrequency -1),
                                               ymax = maxFrequency),
                                  fill = "grey80")
             + ggplot2::annotate("text", x = (maxDistance-3), y = (maxFrequency -0.5),
                                 label = legend_text,colour = "black",size=5,
                                 family="Arial" )
             + ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(2,4,6,8,10) )

             + ggplot2::scale_fill_manual(values = c("FALSE" ="grey",
                                                     "TRUE" = "red3"))
    );
    filename = paste0(geneCommonName,".png")
    ggplot2::ggsave(filename,width =10, height =7, dpi=300)
    # print plot it could change to save the graph image somewhere
    print(plot);

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
        levels(frequenciesDF$Var1) <- c(levels(frequenciesDF$Var1), c(1:25))
        if(maxDistFrequenciesDF < maximumDistance){
            for(i in (maxDistFrequenciesDF+1):25){

                newRow <- (data.frame("Var1" = i,
                                       "Freq" = 0,
                                       "Associations" = FALSE))
                frequenciesDF <- rbind(frequenciesDF, newRow)

            }
        }
    return <- frequenciesDF

}
