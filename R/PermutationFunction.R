#' Permutation test to evauluate if gene-metabolite associations pairs
#' are significantly closer than randomly selected gene-metabolite pairs
#'
#' @param data, data.frame() colunm c(gene, metabolite) representing association
#' between them for each row geneMeasured, list of evaluated genes in study
#' (associated and not associated) metaboliteMeasured, list of measured
#' metabolites in study (associated and not associated)
#' permutation, number of permutation to execute
#' output, "medians" median value of permutated sets - "pvalue" - "histogram"
#'
#' @keywords permutation
#' @export
#' @examples permutationFunction("hsa01100", data, gene, metabolite, 1000)

permutationFunction <-
    function(pathwayId,data,geneMeasured,metaboliteMeasured,
             permutation,output = c("medians","pvalue","histogram")) {

        geneMeasured <- c(t(geneMeasured))
        metaboliteMeasured <- c(t(metaboliteMeasured))

        #create graph object
        graphe <-  createGraphFromPathway(pathwayId);

        #initalise vector for medians of all permutations
        permutatedMedians <- rep(NA,permutation)

        # change data to eliminate associations not where the gene or
        # the metabolite or both are not in the graph
        geneList <- (data[,1])

        rGeneList <- numberOfReactions(graphe@edgeDF,geneList)

        metaboliteList <- (data[,2])
        rMetaboliteList <-
            numberOfMetabolites(graphe@nodeDF, metaboliteList)

        tempDf1 <- data.frame(cbind(
            g1 = as.vector(geneList),
            g2 = as.vector(as.numeric(rGeneList)),
            m1 = as.vector(metaboliteList),
            m2 = as.vector(as.numeric(rMetaboliteList))
        ))


        tempDf1 <- removeNotInGraph(tempDf1)
        rGeneList <- as.numeric(as.vector(tempDf1[,2]))
        rMetaboliteList <- as.numeric(as.vector(tempDf1[,4]))
        data <- subset(tempDf1, select = c(1,3))

        geneList <- as.vector(unique(data[,1]))
        metaboliteList <- as.vector(unique(data[,2]))

        # making sure there is no doubles in the list of all gene measured
        geneMeasured <- unique(geneMeasured)

        # generate vectors with number of reactions for catalysed by gene
        rPossibleGeneToPermutate <-
            numberOfReactions(graphe@edgeDF,geneMeasured)

        # remove gene not in the graph from list of all genes
        tempDf <- data.frame(cbind(g1 = as.vector(geneMeasured),
                                   g2 = as.vector(
                                       as.numeric(rPossibleGeneToPermutate)
                                   )))

        tempDf <- removeNotInGraph(tempDf)

        tempDf <- tempDf[!tempDf$g1 %in% geneList,]

        possibleGeneToPermutate <- as.vector(tempDf[,1])

        rPossibleGeneToPermutate <-
            as.numeric(as.vector(tempDf[,2]))

        # making sure there is no doubles in the list of all metabolites measured
        metaboliteMeasured <- as.vector(unique(metaboliteMeasured))

        # generate vectors with number of repetition of each metabolite
        rPossibleMetaboliteToPermutate <-
            numberOfMetabolites(graphe@nodeDF,
                                metaboliteMeasured)

        # remove metabolites not in the graph from list of all metabolites
        tempDf2 <- data.frame(cbind(
            g1 = as.vector(metaboliteMeasured),
            g2 = as.vector(as.numeric(rPossibleMetaboliteToPermutate))
        ))

        tempDf2 <- removeNotInGraph(tempDf2)
        metaboliteMeasured <- as.vector(tempDf2[,1])
        rPossibleMetaboliteToPermutate <-
            as.numeric(as.vector(tempDf2[,2]))
        distPermutated <- "";

        pb <- txtProgressBar(0,permutation,style = 3)

        for (k in 1:permutation)
        {
            setTxtProgressBar(pb, k)

            metaboliteMeasured <- as.vector(tempDf2[,1])

            rPossibleMetaboliteToPermutate <-
                as.numeric(as.vector(tempDf2[,2]))

            possibleGeneToPermutate <- as.vector(tempDf[,1])

            rPossibleGeneToPermutate <-
                as.numeric(as.vector(tempDf[,2]))

            metaboliteShuffled <-
                sample(metaboliteMeasured,length(metaboliteList))

            geneShuffled <- vector()
            # for all associated, replace by permutated genes
            for (i in 1:length(geneList)) {
                # pick genes that catalyze (+/-1) the same number of reaction

                if (rGeneList[i] != 1)
                {
                    # to get position of possible genes
                    possibleGenesReaction <-
                        as.vector(which(
                            abs(rPossibleGeneToPermutate
                                - rGeneList[i]) <= 1
                        ))

                    # condition if there is no gene that has +/- 1 the number of
                    # reaction of the gene to replace then use the same gene..?
                    if (length(possibleGenesReaction) == 0) {
                        possibleGeneToPermutate <- c(possibleGeneToPermutate,
                                                     geneList[i])

                        genePositionShuffled <-
                            length(possibleGeneToPermutate)
                    }else{
                        genePositionShuffled <- sample(possibleGenesReaction,1)
                    }
                }
                if (rGeneList[i] == 1)
                {
                    # to get position of possible genes
                    possibleGenesReaction <-
                        which((rPossibleGeneToPermutate
                               - rGeneList[i]) == 0)
                    # condition if there is no gene that has +/- 1 the number of
                    # reaction of the gene to replace then use the same gene..?
                    if (length(possibleGenesReaction) == 0) {
                        possibleGeneToPermutate <- c(possibleGeneToPermutate,
                                                     geneList[i])
                        genePositionShuffled <-
                            length(possibleGeneToPermutate)
                    }else{
                        genePositionShuffled <- sample(possibleGenesReaction,1)
                    }
                }
                # add chosen gene to shuffle list
                geneShuffled <- c(geneShuffled,
                    as.character(possibleGeneToPermutate[genePositionShuffled]))

                #' take out chosen gene (and his number of reactions)
                #' from list of possible genes
                possibleGeneToPermutate <-
                    possibleGeneToPermutate[-genePositionShuffled]
                rPossibleGeneToPermutate <-
                    rPossibleGeneToPermutate[-genePositionShuffled]
            }

            permutatedData <- data.frame(data);
            rownames(permutatedData) <-
                c(1:length(permutatedData[,1]))

            #' create new 'data' where associated genes and metabolites
            #' are replace by shuffle genes and metabolites for genes

            for (j in 1:length(geneList)) {

                f <-  apply(permutatedData,1, function(x) {
                    geneListTemp <- paste("\\<", geneList[j], "\\>", sep = '')

                    permutatedData <-
                        gsub(geneListTemp,geneShuffled[j], x)
                    return <- permutatedData;

                })
                permutatedData <- data.frame(t(f));
            }
            for (m in 1:length(metaboliteList)) {
                f1 <-  apply(permutatedData,1, function(x) {
                    metaboliteListTemp <- paste("\\<", metaboliteList[m]
                                                ,"\\>", sep = '')
                    permutatedData <- gsub(metaboliteListTemp,
                                           metaboliteShuffled[m], x)
                    return <- permutatedData;

                })
                permutatedData <- data.frame(t(f1));

            }

            if (k == 1) {
                distPermutated <-
                    getDistanceAsso(pathwayId,
                                    permutatedData,F,"data.frame")[,c(1,3,5)]
                count <- c(rep(1, length(distPermutated[,1])))
                distPermutated <- cbind(distPermutated, "count" = count)

                permutatedMedians[k] <-
                    ceiling(median(distPermutated$distance))

            }else if (k > 1) {
                colnames(permutatedData) <- c("geneKEGGId","metaboliteKEGGId")


                knownDistances <-
                    data.frame(unique(merge(
                        distPermutated,permutatedData,by = c("geneKEGGId",
                                                             "metaboliteKEGGId")
                    )))

                permutatedMedians[k] <-
                    ceiling(median(distPermutated$distance))

                if (nrow(knownDistances) > 0) {

                    for (v in 1:length(knownDistances[,1])) {
                        knownDistances[,1] <-
                            factor(knownDistances[,1],
                                   levels = levels(permutatedData[,1]))
                        knownDistances[,2] <-
                            factor(knownDistances[,2],
                                   levels = levels(permutatedData[,2]))



                        #' get id of known distance rows in permutation set to
                        #' to delete from distance calculation process
                        idRowToDelete <- permutatedData[which(
                          permutatedData$geneKEGGId == knownDistances[v,1]&
                          permutatedData$metaboliteKEGGId == knownDistances[v,2]
                        ),]
                        knownDistances[,1] <-
                            factor(knownDistances[,1],
                                   levels = levels(distPermutated[,1]))
                        knownDistances[,2] <-
                            factor(knownDistances[,2],
                                   levels = levels(distPermutated[,2]))
                        #'get rows to update count
                        idRowToUpdateCount <- distPermutated[which(
                          distPermutated$geneKEGGId == knownDistances[v,1]&
                          distPermutated$metaboliteKEGGId == knownDistances[v,2]

                        ),]

                        if (nrow(idRowToUpdateCount) > 0) {
                            idn <- as.numeric(rownames(idRowToUpdateCount))


                     distPermutated[idn,]$count <- distPermutated[idn,]$count +1

                        }

                        if (nrow(idRowToDelete) > 0) {
                            idn <- as.numeric(rownames(idRowToDelete))

                            permutatedData <-
                                permutatedData[-c(idn),]
                            if (nrow(permutatedData) > 0) {
                                rownames(permutatedData) <-
                                    c(1:length(permutatedData[,1]))
                            }
                            return <- permutatedData
                        }
                    }

                }

                if (nrow(permutatedData) > 0) {
                    distkPermutated <-
                        getDistanceAsso(pathwayId,
                         permutatedData,F,"data.frame")[,c(1,3,5)]
                    count <- c(rep(1, length(distkPermutated[,1])))
                    distkPermutated <- cbind(distkPermutated,"count" = count)
                    distPermutated <-
                        rbind(distPermutated, distkPermutated)
                    distkPermutated <-
                        rbind(distkPermutated,knownDistances)
                    permutatedMedians[k] <-
                        ceiling(median(distkPermutated$distance))
                }else if (nrow(permutatedData) == 0) {
                    permutatedMedians[k] <-
                        ceiling(median(knownDistances$distance))
                }
                print(length(distPermutated$count))
            }

         # print(distPermutated)
        }
        # get median distance associated
        #print(distPermutated)
        print(distPermutated$count)
        distAssociated <-
            getDistanceAsso(pathwayId,data,F, "data.frame")

        medianAssociated <- median(distAssociated$distance)

        # output functions
        if (output == "medians") {
            return <- distPermutated
            #return <- permutatedMedians
        }
        else if (output == "pvalue") {
            pvalue <-
             (sum(permutatedMedians <= medianAssociated) + 1)/(permutation + 1);
            return <- pvalue;
        }
        else if (output == "histogram") {
            # get median distance associated

            permutatedMedians <-
                data.frame("medians" =  permutatedMedians);

            histogramFunction(permutatedMedians, medianAssociated);

        }
    }

histogramFunction <- function(permutatedMedians, medianAssociated, permutation) {

    # get the maximum median value before the Inf value
    infVal <- which(permutatedMedians$medians == Inf)
    temp <- permutatedMedians
    temp[infVal,] <- -1
    maxDistance <- max(temp)

    # get frequency of every value until the maxVal found + Inf val
    frequencies <-
        data.frame(table(factor(
            permutatedMedians$medians,
            levels = c(0:maxDistance,Inf)
        )))

    #percentage of infinite medians
    infValues <- (frequencies[nrow(frequencies),]$Freq / permutation) * 100

    #take out inf values of frequency table
    frequencies <- frequencies[-length(frequencies[,1]),]

    #change all frequency values to percentage
    for (v in 1:length(frequencies[,2])) {
        frequencies[v,2] <- (frequencies[v,2] / permutation) * 100
    }
    ## frequence of inifinite medians
    infFrequencies <- data.frame(frequencies[nrow(frequencies),])

    # maximum of frequence all medians\infinite
    maxFrequencie <- floor(max(frequencies$Freq))

    # define annotation text
    legend_text <- paste("infinite median : ",infValues,"%", sep = "")

    #define title text
    title <- paste("Number of Permutations : ", permutation, sep = "")

    fileName <-   paste("permutations", permutation,".png", sep = "")

    colnames(frequencies) <- c("medians", "Freq")

    plot1 <- ggplot2::ggplot(frequencies,
                             ggplot2::aes(x = medians,y = Freq),
                             environment = environment())

    plot1 <- (
        plot1 + ggplot2::geom_bar(
            stat = "identity",position = "stack",width = 1,
            colour = "black",fill = "grey",
            size = 0.5
        )
        + ggplot2::theme_bw()
        + ggplot2::geom_vline(xintercept = medianAssociated + 1,
                              colour ="red")
        + ggplot2::theme(
            panel.border = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            text = ggplot2::element_text(size = 12, family = "Arial"),
            axis.line = ggplot2::element_line(colour = "black")

        )
        + ggplot2::geom_rect(data = frequencies,
                             ggplot2::aes(xmin = maxDistance+1 -7,
                                          xmax = maxDistance+1,
                                          ymin = maxFrequencie -1,
                                          ymax = maxFrequencie),
                             fill = "grey80")
        + ggplot2::annotate("text", x = maxDistance +1 -3.5, y = maxFrequencie -0.5,
                 label = legend_text,colour = "black",size=5,
                family="Arial" )
        + ggplot2::xlab("Permutated Medians")
        + ggplot2::ylab("Frenquency (%)")
        + ggplot2::scale_y_continuous(expand = c(0,0))
        + ggplot2::ggtitle(title)
    );


    ggplot2::ggsave(fileName, dpi = 300)
    # print plot it could change to save the graph image somewhere
    print(plot1);


}


numberOfReactions <- function(reactionDF,geneList) {

    f <- lapply(geneList, function(x) {

        x <- paste("\\",x,"\\b",sep="")

        nbreReactions <- length(grep(x, reactionDF$ko))

        return <- nbreReactions;
    })
    f <- do.call(rbind, f)
    f <- as.vector(f)

    return <- f
}

numberOfMetabolites <- function(nodeDF,metaboliteList) {
    f <- lapply(metaboliteList, function(x) {
        nbreMetabolites <- length(grep(x, nodeDF$keggId))

        return <- nbreMetabolites;
    })
    f <- do.call(rbind, f)
    f <- as.vector(f)

    return <- f
}


removeNotInGraph <- function(df) {
    row_sub <- apply(df, 1, function(row)
        all(row != 0))
    df <- df[row_sub,]

    return <- df;

}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <-
    function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
        library(grid)

        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)

        numPlots = length(plots)

        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
            # Make the panel
            # ncol: Number of columns of plots
            # nrow: Number of rows needed, calculated from # of cols
            layout <-
                matrix(seq(1, cols * ceiling(numPlots / cols)),
                       ncol = cols, nrow = ceiling(numPlots / cols))
        }

        if (numPlots == 1) {
            print(plots[[1]])

        } else {
            # Set up the page
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(nrow(layout),
                                                       ncol(layout))))

            # Make each plot, in the correct location
            for (i in 1:numPlots) {
                # Get the i,j matrix positions of the regions that contain this subplot
                matchidx <-
                    as.data.frame(which(layout == i, arr.ind = TRUE))

                print(
                    plots[[i]], vp = viewport(
                        layout.pos.row = matchidx$row,
                        layout.pos.col = matchidx$col
                    )
                )
            }
        }
    }
