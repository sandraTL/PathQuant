#' Function that output a heatmap to visualize distance computed between every
#' gene-metabolite associations in input using the overview map only
#'
#' Function computing shortest distance between every genes and metabolites in
#' a gene-metabolite pairs of your association dataset on a graph model of
#' KEGG map selected, where nodes are metabolites and reactions/genes are edges.
#'
#' If a gene or a metabolite is present on multiple edges or nodes, then the
#' shortest distance is selected within the group

#'
#' Output: Heatmap of distances calculated between associated genes-metabolites.
#' Columns represent genes and rows represent metabolites. The calculated
#' distance is shown in each cell with the corresponding color code
#' (from red - closest; to yellow - farthest).
#'
#' @param association dataframe with 2 columns, where each line reprensents a
#'        uniq association. First column are genes and sencond column are
#'        metabolites. Only use KEGG Ids.
#' @param display common names of KEGG ids input (TRUE) or
#'        display KEGG Ids for axis description.
#'
#' @keywords graph, heatmap, shortest.distance, KEGG
#' @export
#' @examples heatmap(shinAndAlDF, TRUE)

heatmap <- function(association, commonNames = FALSE){

    # print("heatmapAsso")

    # mError2 <-"Sorry, for each pairs of gene/metabolite entered, either the gene
    #          or the metabolite or both weren't mapped on the selected pathway.
    #          Thus, no distance was calculated"
     fileName <- paste("heatmap", "hsa01100", ".png", sep = "")


    # test_heatmap(pathwayId, association);


    # Get the distance for the associated set
      data.result.ob <- create.data.restul.ob(association,
                                pathway = c("hsa01100"),
                                F,
                                F,
                                F)

      # Keep only results with a finite numerical value
      df.result <- df.sub.data.result.numeric(data.result.ob, T)[,c(2,3)]

      # Get a uniq list of gene and metabolite to build heatmap
      gene.uniq <- sort(unique(df.result[,1]))
      metabolite.uniq <- sort(unique(df.result[,2]))

      # Build heatmap dataset to compute distance for every possible
      # pair between both list
      df.heatmap <- data.frame(
                        Row = rep(gene.uniq, each= length(metabolite.uniq)),
                        Col = rep(metabolite.uniq, times= length(gene.uniq)))

      # Define the pairs from the original dataset (TRUE)
      is.asso <- getAssociationForHeatmap(df.result, df.heatmap)

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

# getAssociationForHeatmap<- function(data, heatmapDf){
#     boolIsAsso <- FALSE;
#     f <- apply(heatmapDf, 1, function(x){
#         isAssociation <-  length(data[data[1] == x[1] & data[2] == x[2]])
#         if(isAssociation >= 2){
#             boolIsAsso <- TRUE;
#         }else{
#             boolIsAsso<- FALSE
#         }
#     })
# }
#



## new style heatmap ......
#
# heatmapAsso_1 <- function(pathwayId){
#     title <- pathwayId;
#
#     pathwayId <- gsub("hsa:", "hsa", pathwayId)
#
#     mError2 <-"Sorry, for each pairs of gene/metabolite entered, either the gene
#     or the metabolite or both weren't mapped on the selected pathway.
#     Thus, no distance was calculated"
#     fileName <- paste("heatmapAsso", pathwayId, ".png", sep = "")
#     #test_heatmap(pathwayId, association);
#     graphe <-  createGraphFromPathway(pathwayId);
#     # rGeneList<-numberOfReactions(graphe@edgeDF,association[,1])
#     # rMetaboliteList <- numberOfMetabolites(graphe@nodeDF, association[,2])
#     # tempDf1 <- data.frame(cbind(g1 = as.vector(association[,1]),
#     #                             g2 = as.vector(as.numeric(rGeneList)),
#     #                             m1 = as.vector(association[,2]),
#     #                             m2 = as.vector(as.numeric(rMetaboliteList))))
#     #
#     # tempDf1 <- removeNotInGraph(tempDf1)
#
#
#     iems_res = gdata::read.xls("/Users/sandra/Desktop/Metabolomique/DataSets/Thiele/PQ_input_NewHeatmap.xlsx")
#     iems_resNoNA <- iems_res[complete.cases(iems_res),]
#
#
#     if(nrow(iems_resNoNA) == 0 ){
#         print("Allo")
#         stop(mError2, call. = FALSE)
#     }
#
#
#     association <- subset(iems_resNoNA[,c(2,4)])
#     gene <- data.frame(c(iems_resNoNA[2]))
#     metabolite <- data.frame(c(iems_resNoNA[4]))
#
#     # gene<- data.frame(c(association[1]))
#     # metabolite <- data.frame(c(association[2]))
#
#
#     #AllSP <- getDistanceAll(pathwayId, gene, metabolite);
#
#     AllSP1 <- getDistanceAll("hsa01100", gene, metabolite)
#
#     AllSP1[mapply(is.infinite, AllSP1)] <- NA
#
#     maxDistance <- max(AllSP1, na.rm=TRUE)
#
#     AllSP1[mapply(is.na, AllSP1)] <- Inf
#
#
#     AllSP <- data.frame()
#
#     for(row in 1:nrow(AllSP1)){
#
#         rowName <- as.vector(rownames(AllSP1))
#         print(rowName[row])
#         row <- as.vector(table(factor(t(AllSP1[row,]), levels=c(0:25,Inf))))
#         AllSP  <- rbind(AllSP , rowName = row)
#
#
#     }
#
#     colnames(AllSP) <- c(0:25,Inf)
#     rownames(AllSP) <- rownames(AllSP1)
#
#
#
#     AllSP[AllSP == 0] <- Inf;
#
#     # geneCommonName <- getCommonNames(as.vector(unlist(rownames(AllSP))), "gene")
#     #
#     # geneCommonName <- as.vector(unlist(geneCommonName))
#     #
#     # metaboliteCommonName <- getCommonNames(as.vector(unlist(colnames(AllSP))),
#     #                                        "metabolite")
#     # metaboliteCommonName <- as.vector(unlist(metaboliteCommonName))
#
#
#     RtoWrange<-colorRampPalette(c("red", "white"))
#     WtoGrange<-colorRampPalette(c("white", "green"))
#
#     tempDat <- data.frame(Row = rep(rownames(AllSP), each= ncol(AllSP)),
#                           Col = rep(colnames(AllSP), times= nrow(AllSP)))
#     #  isAsso <- getAssociationForHeatmap(association, tempDat)
#     vec <- list();
#     for(i in 1:nrow(AllSP)){vec <- c(vec,as.vector(AllSP[i,]))}
#     distance <- t(data.frame("distance" = as.vector(vec)))
#
#       # Bind with heatmap dataset
#       df.heatmap <- cbind(df.heatmap, "Associations" = as.vector(is.asso))
#
#       #Compute distance for heatmap dataset
#       heatmap.result.ob <- create.data.restul.ob(df.heatmap,
#                                                 pathway = c("hsa01100"),
#                                                 F,
#                                                 F,
#                                                 F)
#
#        # Bind results
#       df.heatmap <- cbind(df.heatmap,
#                    "Distance" = as.vector(heatmap.result.ob@distance))
#       colnames(df.heatmap) <- c("Row", "Col", "Associations", "Distance")
#
#       if(commonNames){
#
#           # Get commun names for gene and metabolite
#           gene.names <- getCommonNames(as.vector(unlist(gene.uniq)), "gene")
#           gene.names <- as.vector(unlist(gene.names))
#
#           metabolite.names <-
#               getCommonNames(as.vector(unlist(metabolite.uniq)), "metabolite")
#           metabolite.names <- as.vector(unlist(metabolite.names))
# >>>>>>> sandra
#
#           #Change ids to common names
#           df.heatmap$Row <- rep(gene.names, each = length(metabolite.names))
#           df.heatmap$Col <- rep(metabolite.names, times = length(gene.names))
#
# <<<<<<< HEAD
#     dat <- data.frame( Row = rep(rownames(AllSP), each= ncol(AllSP)),
#                        Col = rep(colnames(AllSP), times= nrow(AllSP)),
#                        Distance = c(t(distance))
#                        # Associations = isAsso
#
#     );
#    # print(dat)
#     # create frame to color edges of associated genes and metabolites
#     frames = dat[dat$Associations, c("Row","Col")]
#     frames$Row = as.integer(frames$Row)
#     frames$Col = as.integer(frames$Col)
#
#     # from discrete to continuous values...
#     if(!is.infinite(min(dat$Distance))){
#
#         dat$Distance <- as.numeric(as.character(dat$Distance))
#
#     }else{
#         dat$Distance <- as.character(dat$Distance)
#
#     }
#     colors <- c("#BD0026","#E31A1C","#FC4E2A","#FD8D3C","#FEB24C",
#                 "#FED976", "#FFEDA0", "#FFFFCC")
#     p2 = ggplot2::ggplot(data=dat) +
#         ggplot2::geom_raster(ggplot2::aes(x=Row, y=Col, fill=Distance)) +
#         ggplot2::stat_bin2d(ggplot2::aes(x=Row, y=Col, fill=Distance),binwidth = c(1,1))+
#         ggplot2::theme(
#             panel.border = ggplot2::element_rect(colour="black",fill=NA,size=2),
#             panel.grid.major = ggplot2::element_blank(),
#             panel.grid.minor = ggplot2::element_blank(),
#
#             text = ggplot2::element_text(size=12),
#             axis.text=ggplot2::element_text(colour="black"),
#             axis.text.x = ggplot2::element_text(angle=315,vjust=1,hjust=0))+
#
#
#         ggplot2::xlab("Genes")+
#         ggplot2::ylab("Metabolites")+
#         # ggplot2::coord_fixed(ratio=0.5)+
#         # ggplot2::geom_rect(data=frames, size=1, fill=NA, colour="black",
#         #                    ggplot2::aes(xmin=Row-0.5, xmax=Row+0.5, ymin=Col-0.5, ymax=Col + 0.5))+
#         ggplot2::geom_text(label = as.numeric(dat$Distance, 1),
#                            size = 2,ggplot2::aes(x = Row, y = Col)) +
#
#         ggplot2::ggtitle(title)
#
#     ##Fill gradient only if all values aren't Inf
#     if(typeof(dat$Distance)=="double"){
#         p2 <-p2 + ggplot2::scale_fill_gradient2(low=RtoWrange(100), mid=WtoGrange(100), high="gray")}
#     else{ p2 <-p2 + ggplot2::scale_fill_grey(start = 0.5, end = 0.5) }
#
#     p2 <- p2 + ggplot2::ggsave(fileName, dpi=300)
#     print(p2);
# =======
#           # Order in alphabetical order before defining the frame ids
#           # because the heatmap function orders automatically
#           df.heatmap <- df.heatmap[order(df.heatmap[1], df.heatmap[2]),]
#
#           # Set the number of rows in increasing order to built frame id
#           # for associations
#           row.names(df.heatmap) <- c(1:nrow(df.heatmap))
#
#       }
#
#     # Define frame ids ready for black highlight of associations
#     frame = df.heatmap[df.heatmap$Associations, c("Row","Col")]
#     asso <- as.numeric(row.names(frame))
#     row.asso.frame <- lapply(asso, function(x) ceiling(x/length(metabolite.uniq)))
#     col.asso.frame <- lapply(asso, function(x) {
#         if(x%%(length(metabolite.uniq)) == 0){
#             length(metabolite.uniq)
#         } else {
#             x%%(length(metabolite.uniq))
#         }
#       }
#     )
#
#     frame <- data.frame("Row" = as.vector(unlist(row.asso.frame)),
#                         "Col" = as.vector(unlist(col.asso.frame)))
#     frame$Row = as.integer(frame$Row)
#     frame$Col = as.integer(frame$Col)
#
#     # Set heat colors
#     colors <- c("#BD0026","#E31A1C","#FC4E2A","#FD8D3C","#FEB24C",
#                 "#FED976", "#FFEDA0", "#FFFFCC")
# >>>>>>> sandra
#
#     # ggplot2 heatmap preparation
#     p2 = ggplot2::ggplot(df.heatmap, ggplot2::aes(x = Row, y = Col, fill = Distance),
#                                                          binwidth = c(1,1)) +
#     ggplot2::geom_tile() +
#     ggplot2::scale_fill_gradientn(colours = colors) +
#     ggplot2::theme(
#        panel.border = ggplot2::element_rect(colour="black",fill=NA,size=2),
#                      panel.grid.major = ggplot2::element_blank(),
#                      panel.grid.minor = ggplot2::element_blank(),
#
#                      text = ggplot2::element_text(size=12),
#                      axis.text = ggplot2::element_text(colour="black"),
#         axis.text.x = ggplot2::element_text(angle=315,vjust=1,hjust=0))+
#
#      ggplot2::coord_equal() +
#      ggplot2::xlab("Genes")+
#      ggplot2::ylab("Metabolites")+
#      ggplot2::ggtitle("Heatmap")+
#      ggplot2::guides(fill=ggplot2::guide_legend(title="srd"))+
#      ggplot2::coord_fixed(ratio=0.5)+
#      ggplot2::geom_rect(data=frame, size=1, fill=NA, colour="black",
#      ggplot2::aes(xmin=Row-0.5, xmax=Row+0.5, ymin=Col-0.5, ymax=Col + 0.5))+
#      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
#      ggplot2::geom_text(label = as.numeric(df.heatmap$Distance, 1),
#                        size = 2,ggplot2::aes(x = Row, y = Col))
#      print(p2)
#
#      return <- df.heatmap
# }
#
# <<<<<<< HEAD
#
#
# =======
getAssociationForHeatmap<- function(data, heatmapDf){

    # print("getAssociationForHeatmap")

    boolIsAsso <- FALSE;
    f <- apply(heatmapDf, 1, function(x){
             isAssociation <-  length(data[data[1] == x[1] & data[2] == x[2]])

                   if(isAssociation >= 2){
                       boolIsAsso <- TRUE;
                   }else{
                       boolIsAsso<- FALSE
                   }
               })
    return <- f
}








