#' Distance distribution plots for single gene
#'
#' Function ploting the distribution of distances between a gene and all
#' measured metabolite, highlithing the distance of its associated metabolites.
#'
#' The plot is depicted as frequency bars, which represent
#' the number of metabolites at a given distance for the selected gene.
#' Frequency bars are shown in grey for metabolites that are not associated with
#' the selected gene and in red if there is at least one metabolite associated
#' with this gene.
#'
#' If a gene or a metabolite is present on multiple edges or nodes, then the
#' shortest distance is selected.
#'
#' @param pathway Selected pathway. Only use KEGG Ids.
#' @param gene Gene selected. Only use KEGG Ids.
#' @param pairs Dataframe with 2 columns, where each line reprensents an
#'        association. First column are the genes and the sencond column are the
#'        metabolites. Only use KEGG Ids.
#' @param metabolites Dataframe of 1 column containing a lsit of selected
#'        metabolites. Default = "All". All metabolite of the selected pathway
#'        map. Only use KEGG Ids.
#' @keywords graph, shortest path, KEGG.
#' @export
#' @examples gene.distribution(hsa01100, "hsa:1373" , pairs.df, metabolite.df)


gene.distribution <- function(pathway, gene, pairs, metabolites){

   pathway <- gsub(":", "", pathway);
   #  test_distributionGene(pathwayId, association, metabolite, gene);

   # get list of uniq metabolites in association
   # metabolite.uniq <-
   #    unique(association$Metabolite[!is.na(association$Metabolite)])

   if(metabolites == "All"){
   metabolites <- unique(getListNodeFromKGML(pathway)$keggId)
   }
   # Get the real association to highlight bar of their distance
   asso.gene <- pairs[pairs$Gene == gene,]

   # Remove Na values
   asso.gene <- asso.gene[complete.cases(asso.gene), ]

   geneDF <- data.frame("gene" = c(gene))

   # Prepare association to compute distance
   df.barplot <- data.frame(
        Gene = rep(gene, times= length(metabolites)),
        Metabolite = metabolites)


   data.result.ob <- get.srd(df.barplot,
                             pathway = c(pathway),
                             F,
                             F,
                             F)

      # get all shortest paths from association entry
   shortestsPathsDF <- data.result.ob[,c(2,3,8)]

   # adjust gene parameter
   gene1 <- gsub(":", ".", gene);

   shortestsPathsDF <-
        shortestsPathsDF[(complete.cases(shortestsPathsDF$distance)),]
    rownames(shortestsPathsDF) <- c(1:nrow(shortestsPathsDF))

    maxVal <- getMaxValIgnoreInfVal(shortestsPathsDF)

    # get frequency of every value until the maxVal found + Inf val
    frequenceDF <- data.frame(table(factor(shortestsPathsDF$distance,
                                            levels=c(0:maxVal,Inf))))


    # set column names of frequenceDF
    colnames(frequenceDF) <- c("Distance", "Freq")

    # get associations values for all distance pair calculated
    associationsMetaboDF <- getAssociationsDF(shortestsPathsDF,asso.gene)

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

    # print("barplotFunctionGeneToAllMetabo")

    # initiating variable for barplotGraph
    geneCommonName <- getCommonNames(c(gene), "gene")
    numInfValue <- frequenceDF[frequenceDF$Distance == Inf,][,2]
    frequenceDF <- frequenceDF[-length(frequenceDF[,1]),]
    maxDistance <- as.numeric(as.character(frequenceDF[nrow(frequenceDF),][,1]))
    maxFrequency <- max(frequenceDF$Freq, na.rm = TRUE)

    y.axis <- round(maxFrequency/2)
    breaks <- seq(1, maxFrequency, by=2)


    legend_text <- paste("Inf srd : ",numInfValue, sep = "")

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
                 axis.line.y = ggplot2::element_line(color="black"),
                 plot.title = ggplot2::element_text(hjust = 0.5)
                              )
             + ggplot2::xlab("srd from gene")
             + ggplot2::ylab("Metabolite count")
             + ggplot2::ggtitle(geneCommonName)
             + ggplot2::coord_fixed(ratio = 1)

             + ggplot2::annotate("text", x = (maxDistance-3),
                                 y = (maxFrequency -0.5),
                                 label = legend_text,
                                 colour = "black",
                                 size=4 )

              + ggplot2::scale_y_continuous(expand = c(0,0),
                                            breaks = breaks)

             + ggplot2::scale_fill_manual(values = c("FALSE" ="grey",
                                                     "TRUE" = "red3"),
                                          guide = FALSE)

    );

    filename = paste0(geneCommonName,".png")

    ggplot2::ggsave(filename,width =15, height =7, dpi=300)
    print(plot);

}




getFrequenceAssociationsDF <- function(frequenceDistDF,shortestsPathsDF,gene){

  #  print("getFrequenceAssociationsDF")

    test = FALSE;
    results <- data.frame();

    # creation of vector to fill bar colors automatically
    for(row1 in 1:nrow(frequenceDistDF)){
        test = FALSE;
        for(row2 in 1:nrow(shortestsPathsDF)){



        if(frequenceDistDF[row1,"Distance"] == shortestsPathsDF[row2,"distance"]){

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

    associatedMetaboDF[,2] <- factor(associatedMetaboDF[,2],
                                 levels = levels(assoDistanceDF[,2]))

    # initiation of values
    test = FALSE;
    results <- data.frame();

    # buiding boolean associations DF
    for(row1 in 1:nrow(assoDistanceDF)){
        test = FALSE;


        if(is.na(assoDistanceDF[row1,2])){
            test = FALSE
        }else{
          for(row2 in 1:nrow(associatedMetaboDF)){
            if(is.na(associatedMetaboDF[row2,2])){
                test = FALSE
            }else if(assoDistanceDF[row1,2] == associatedMetaboDF[row2,2]){

                test <- TRUE;

                break;
            }else test<- FALSE;
        }
      }
        results <- rbind(results,test);
        colnames(results) <- c("Associations")

    }
    return <- results;
}

getAssociatedMetaboByGene <- function(data, gene){

    selectedRows <- data[data[,1] == gene,];

    associatedMetabo <- selectedRows[,2];

    return <- associatedMetabo;

}


# function to modify the frequency data.frame by adding more distance that are
# then the one already defined, ulterior use.

barplot_adjustMaximalDistance <- function(maximumDistance, frequenciesDF,
                                          maxDistFrequenciesDF){

       # print("barplot_adjustMaximalDistance")
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

groupedBarPlot <- function(data){

   # print("groupedBarPlot")

    statsData <- getDataForBarPlot(data)
    colnames(statsData) <- c("classNames","all_uniq", "unmapped", "unmappedAsso")
   # print(statsData)
    statsData.m <- reshape2::melt(statsData, "id.vars"="classNames")
   # print(statsData.m)
    plot<- ggplot2::ggplot(statsData.m, ggplot2::aes(classNames, value)) +
    ggplot2::geom_bar(ggplot2::aes(fill = variable),
                      position = "dodge", stat="identity") +
    ggplot2::xlab("Metabolite Class")+
    ggplot2::ylab("count")+
    ggplot2::scale_x_discrete(labels =
                     function(x) stringr::str_wrap(x, width = 0.6))


   # ggplot2::ggsave("stats",width =10, height =7, dpi=300)
    print(plot);

}

groupedBarPlot <- function(data){

    #print("groupedBarPlot")

    statsData <- getDataForBarPlot(data)
    colnames(statsData) <- c("classNames","all_uniq", "unmapped", "unmappedAsso")
   # print(statsData)
    statsData.m <- reshape2::melt(statsData, "id.vars"="classNames")
   # print(statsData.m)
    plot<- ggplot2::ggplot(statsData.m, ggplot2::aes(classNames, value)) +
        ggplot2::geom_bar(ggplot2::aes(fill = variable),
                          position = "dodge", stat="identity") +
        ggplot2::theme(
            text = ggplot2::element_text(size=14),
            axis.text=ggplot2::element_text(colour="black", size = 12))+
        ggplot2::xlab("Metabolite Class")+
        ggplot2::ylab("count")+
        ggplot2::scale_x_discrete(labels =
                                      function(x) stringr::str_wrap(x, width = 0.6))


    # ggplot2::ggsave("stats",width =10, height =7, dpi=300)
    print(plot);

}

#' srd.distribution
#'
#' Function ploting the distribution of srds computed by get.srd function.
#'
#' The plot is depicted as frequency bars, which represent
#' the number of pairs at a given srd value.
#' Frequency bars are shown in heatmap colors to depict the small from larger
#' srd values.
#'
#'
#' @param col of srd values of the data.frame obtained by get.srd function
#' @keywords graph, srd, KEGG, barplot.
#' @export
#' @examples srd.distribution(res[,7])

srd.distribution <- function(distance){

    statsData.m <- getDistanceDistribution(distance)

    distance <- as.character(as.vector(statsData.m$distance))

    distanceData <- data.frame("distance" = as.character(distance),
                              "count" = statsData.m$count)

    fitColors <-  heat.colors(nrow(distanceData), 1)
    distanceData <- cbind(distanceData, "colors" = fitColors)
    p<-ggplot2::ggplot(distanceData,
                ggplot2::aes(x=distance, y=count, fill=colors)
                       ,environment = environment())
       p <- (p +  ggplot2::geom_bar(stat="identity", colour="black",
                         position="identity")+
       ggplot2::theme(
            text = ggplot2::element_text(size=14),
            axis.text=ggplot2::element_text(colour="black", size = 14))+
       ggplot2::scale_fill_manual(values=fitColors)+
       ggplot2::labs(x="srd",y="count")+
       ggplot2::ggtitle("srd distribution")+
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
       ggplot2::theme(legend.position="none")+
       ggplot2::scale_x_discrete(limits=distance))



 print(p);


}


getDistanceDistribution <- function(distance){

    #print("getDistanceDistribution")
    # test <- apply(data, 1,function(x) as.numeric(distance))

    distTable <-table(distance)


    distNames <- as.vector(names(distTable))
    distCount <- as.vector(distTable)
    distDF <- data.frame("distance" = distNames,
                         "count" = distCount)

    return <- distDF;

}


fitColors <- function(barNumber, colors){

   repeatTimes <- floor(barNumber/length(colors))

   newColors <- list()
   y <- 0;
   colNum <- 1;
   for(i in 1:barNumber){
   #print(colNum)
       if(y < repeatTimes){

           newColors <- c(newColors, colors[colNum])
           y <- y+1
       }else if(y == repeatTimes){
           y <- 0
           colNum <- colNum + 1;
           newColors <- c(newColors, colors[colNum])
           y <- y+1
       }
   }

   newColors <- unlist(newColors)

   return <- newColors
}
