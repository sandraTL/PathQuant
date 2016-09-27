#' Function that output a heatmap to visualize distance calculated between every
#' gene-metabolite associations in input.
#'
#' Function calculting shortest distance between every genes and metabolites in
#' a gene-metabolite pairs of your association parameter on a graph model of
#' KEGG map selected, where nodes are metabolites and reactions are edges.
#'
#' If a gene or a metabolite is present on multiple edges or nodes, then the
#' shortest distance is selected.
#'
#' Output: Heatmap of distances calculated between associated genes-metabolites.
#' Columns represent genes and rows represent metabolites. The calculated
#' distance is shown in each cell with the corresponding color code
#' (from red - closest; to yellow - farthest).
#'
#' @param pathwayId KEGG Id of selected pathway.
#' @param association Dataframe with 2 columns, where each line reprensents an
#'        association. First column are genes and the sencond column are
#'        metabolites. Only use KEGG Ids.
#' @param display common names of KEGG ids input or display KEGG Ids for axis
#'        description.
#'
#' @keywords graph, heatmap, shortestDistance, KEGG
#' @export
#' @examples heatmapAsso("hsa01100", shinAndAlDF, TRUE)

heatmapAsso <- function(pathwayId, association, commonNames = FALSE){
    title <- pathwayId;

    pathwayId <- gsub("hsa:", "hsa", pathwayId)

    mError2 <-"Sorry, for each pairs of gene/metabolite entered, either the gene
             or the metabolite or both weren't mapped on the selected pathway.
             Thus, no distance was calculated"
    fileName <- paste("heatmapAsso", pathwayId, ".png", sep = "")
    #test_heatmap(pathwayId, association);
    graphe <-  createGraphFromPathway(pathwayId);
    rGeneList<-numberOfReactions(graphe@edgeDF,association[,1])
    rMetaboliteList <- numberOfMetabolites(graphe@nodeDF, association[,2])
    tempDf1 <- data.frame(cbind(g1 = as.vector(association[,1]),
                                g2 = as.vector(as.numeric(rGeneList)),
                                m1 = as.vector(association[,2]),
                                m2 = as.vector(as.numeric(rMetaboliteList))))

     tempDf1 <- removeNotInGraph(tempDf1)

     if(nrow(tempDf1) == 0 ){
         print("Allo")
         stop(mError2, call. = FALSE)
     }
    association <- subset(tempDf1[,c(1,3)])
    gene<- data.frame(c(association[1]))
    metabolite <- data.frame(c(association[2]))

    AllSP <- getDistanceAll(pathwayId, gene, metabolite);

    if(commonNames){
        geneCommonName <-
            getCommonNames(as.vector(unlist(rownames(AllSP))), "gene")
        geneCommonName <- as.vector(unlist(geneCommonName))
        metaboliteCommonName <-
            getCommonNames(as.vector(unlist(colnames(AllSP))), "metabolite")
        metaboliteCommonName <- as.vector(unlist(metaboliteCommonName))
        tempDat <- data.frame(Row = rep(rownames(AllSP), each= ncol(AllSP)),
                          Col = rep(colnames(AllSP), times= nrow(AllSP)))
        isAsso <- getAssociationForHeatmap(association, tempDat)

         dat <- data.frame(Row = rep(geneCommonName, each= ncol(AllSP)),
                       Col = rep(metaboliteCommonName, times= nrow(AllSP)),
                       Distance = c(t(AllSP)),
                       Associations = isAsso)

    }else{
        tempDat <- data.frame(Row = rep(rownames(AllSP), each= ncol(AllSP)),
                              Col = rep(colnames(AllSP), times= nrow(AllSP)))

        isAsso <- getAssociationForHeatmap(association, tempDat)
        dat <- data.frame(
            Row = as.vector(rep(rownames(AllSP), each= ncol(AllSP))),
            Col = as.vector(rep(colnames(AllSP), times= nrow(AllSP))),
                          Distance = c(t(AllSP)),
                          Associations = isAsso)
    }

    # create frame to color edges of associated genes and metabolites
    frames = dat[dat$Associations, c("Row","Col")]
    frames$Row = as.integer(frames$Row)
    frames$Col = as.integer(frames$Col)

    # from discrete to continuous values...
     if(!is.infinite(min(dat$Distance))){
          dat$Distance <- as.numeric(as.character(dat$Distance))
      }else{
          dat$Distance <- as.character(dat$Distance)
      }

    colors <- c("#BD0026","#E31A1C","#FC4E2A","#FD8D3C","#FEB24C",
                "#FED976", "#FFEDA0", "#FFFFCC")
    p2 = ggplot2::ggplot(data=dat) +
        ggplot2::stat_bin2d(ggplot2::aes(x=Row, y=Col, fill=Distance),
                            binwidth = c(1,1))+
        ggplot2::theme(
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=2),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),

            text = ggplot2::element_text(size=12),
            axis.text=ggplot2::element_text(colour="black"),
            axis.text.x = ggplot2::element_text(angle=315,vjust=1,hjust=0))+
        ggplot2::xlab("Genes")+
        ggplot2::ylab("Metabolites")+
        ggplot2::coord_fixed(ratio=0.5)+
        ggplot2::geom_rect(data=frames, size=1, fill=NA, colour="black",
        ggplot2::aes(xmin=Row-0.5, xmax=Row+0.5, ymin=Col-0.5, ymax=Col + 0.5))+
        ggplot2::geom_text(label = as.numeric(dat$Distance, 1),
                 size = 2,ggplot2::aes(x = Row, y = Col)) +

        ggplot2::ggtitle(title)

        ##Fill gradient only if all values aren't Inf
        if(typeof(dat$Distance)=="double"){
        p2 <-p2 + ggplot2::scale_fill_gradientn(colours = colors)}
        else{ p2 <-p2 + ggplot2::scale_fill_grey(start = 0.5, end = 0.5) }

        p2 <- p2 + ggplot2::ggsave(fileName, dpi=300)
    print(p2);

}

getAssociationForHeatmap<- function(data, heatmapDf){
    boolIsAsso <- FALSE;
    f <- apply(heatmapDf, 1, function(x){
             isAssociation <-  length(data[data[1] == x[1] & data[2] == x[2]])
                   if(isAssociation >= 2){
                       boolIsAsso <- TRUE;
                   }else{
                       boolIsAsso<- FALSE
                   }
               })
}




## new style heatmap ......

heatmapAsso_1 <- function(pathwayId){
    title <- pathwayId;

    pathwayId <- gsub("hsa:", "hsa", pathwayId)

    mError2 <-"Sorry, for each pairs of gene/metabolite entered, either the gene
    or the metabolite or both weren't mapped on the selected pathway.
    Thus, no distance was calculated"
    fileName <- paste("heatmapAsso", pathwayId, ".png", sep = "")
    #test_heatmap(pathwayId, association);
    graphe <-  createGraphFromPathway(pathwayId);
    # rGeneList<-numberOfReactions(graphe@edgeDF,association[,1])
    # rMetaboliteList <- numberOfMetabolites(graphe@nodeDF, association[,2])
    # tempDf1 <- data.frame(cbind(g1 = as.vector(association[,1]),
    #                             g2 = as.vector(as.numeric(rGeneList)),
    #                             m1 = as.vector(association[,2]),
    #                             m2 = as.vector(as.numeric(rMetaboliteList))))
    #
    # tempDf1 <- removeNotInGraph(tempDf1)


    iems_res = gdata::read.xls("/Users/sandra/Desktop/Metabolomique/DataSets/Thiele/PQ_input_NewHeatmap.xlsx")
    iems_resNoNA <- iems_res[complete.cases(iems_res),]


    if(nrow(iems_resNoNA) == 0 ){
        print("Allo")
        stop(mError2, call. = FALSE)
    }


    association <- subset(iems_resNoNA[,c(2,4)])
    gene <- data.frame(c(iems_resNoNA[2]))
    metabolite <- data.frame(c(iems_resNoNA[4]))

    # gene<- data.frame(c(association[1]))
    # metabolite <- data.frame(c(association[2]))


    #AllSP <- getDistanceAll(pathwayId, gene, metabolite);

    AllSP1 <- getDistanceAll("hsa01100", gene, metabolite)

    AllSP1[mapply(is.infinite, AllSP1)] <- NA

    maxDistance <- max(AllSP1, na.rm=TRUE)

    AllSP1[mapply(is.na, AllSP1)] <- Inf


    AllSP <- data.frame()

    for(row in 1:nrow(AllSP1)){

        rowName <- as.vector(rownames(AllSP1))
        print(rowName[row])
        row <- as.vector(table(factor(t(AllSP1[row,]), levels=c(0:25,Inf))))
        AllSP  <- rbind(AllSP , rowName = row)


    }

    colnames(AllSP) <- c(0:25,Inf)
    rownames(AllSP) <- rownames(AllSP1)



     AllSP[AllSP == 0] <- Inf;

    # geneCommonName <- getCommonNames(as.vector(unlist(rownames(AllSP))), "gene")
    #
    # geneCommonName <- as.vector(unlist(geneCommonName))
    #
    # metaboliteCommonName <- getCommonNames(as.vector(unlist(colnames(AllSP))),
    #                                        "metabolite")
    # metaboliteCommonName <- as.vector(unlist(metaboliteCommonName))


     RtoWrange<-colorRampPalette(c("red", "white"))
     WtoGrange<-colorRampPalette(c("white", "green"))

    tempDat <- data.frame(Row = rep(rownames(AllSP), each= ncol(AllSP)),
                          Col = rep(colnames(AllSP), times= nrow(AllSP)))
  #  isAsso <- getAssociationForHeatmap(association, tempDat)
    vec <- list();
    for(i in 1:nrow(AllSP)){vec <- c(vec,as.vector(AllSP[i,]))}
    distance <- t(data.frame("distance" = as.vector(vec)))


    dat <- data.frame( Row = rep(rownames(AllSP), each= ncol(AllSP)),
                       Col = rep(colnames(AllSP), times= nrow(AllSP)),
                       Distance = c(t(distance))
                      # Associations = isAsso

    );
print(dat)
    # create frame to color edges of associated genes and metabolites
    frames = dat[dat$Associations, c("Row","Col")]
    frames$Row = as.integer(frames$Row)
    frames$Col = as.integer(frames$Col)

    # from discrete to continuous values...
    if(!is.infinite(min(dat$Distance))){

        dat$Distance <- as.numeric(as.character(dat$Distance))

    }else{
        dat$Distance <- as.character(dat$Distance)

    }
    colors <- c("#BD0026","#E31A1C","#FC4E2A","#FD8D3C","#FEB24C",
                "#FED976", "#FFEDA0", "#FFFFCC")
    p2 = ggplot2::ggplot(data=dat) +
         ggplot2::geom_raster(ggplot2::aes(x=Row, y=Col, fill=Distance)) +
        ggplot2::stat_bin2d(ggplot2::aes(x=Row, y=Col, fill=Distance),binwidth = c(1,1))+
        ggplot2::theme(
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=2),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),

            text = ggplot2::element_text(size=12),
            axis.text=ggplot2::element_text(colour="black"),
            axis.text.x = ggplot2::element_text(angle=315,vjust=1,hjust=0))+


        ggplot2::xlab("Genes")+
        ggplot2::ylab("Metabolites")+
       # ggplot2::coord_fixed(ratio=0.5)+
        # ggplot2::geom_rect(data=frames, size=1, fill=NA, colour="black",
        #                    ggplot2::aes(xmin=Row-0.5, xmax=Row+0.5, ymin=Col-0.5, ymax=Col + 0.5))+
         ggplot2::geom_text(label = as.numeric(dat$Distance, 1),
                            size = 2,ggplot2::aes(x = Row, y = Col)) +

        ggplot2::ggtitle(title)

    ##Fill gradient only if all values aren't Inf
    if(typeof(dat$Distance)=="double"){
        p2 <-p2 + ggplot2::scale_fill_gradient2(low=RtoWrange(100), mid=WtoGrange(100), high="gray")}
    else{ p2 <-p2 + ggplot2::scale_fill_grey(start = 0.5, end = 0.5) }

    p2 <- p2 + ggplot2::ggsave(fileName, dpi=300)
    print(p2);

}





