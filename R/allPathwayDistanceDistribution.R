
#' Function that output a heatmap to visualize distance calculated between every
#' gene-metabolite associations in input.
#'
#' Function calculting shortest distance between every genes and metabolites in
#' a gene-metabolite pairs of your association parameter on a graph model of KEGG map
#' selected, where nodes are metabolites and reactions are edges.
#'
#' If a gene or a metabolite is present on multiple edges or nodes, then the shortest distance is selected.
#'
#' Output : Heatmap of distances calculated between associated genes-metabolites.
#' Columns represent genes and rows represent metabolites. The calculated distance is shown in each cell with the corresponding color code (from red - closest; to yellow - farthest).
#'
#' @param pathwayId KEGG Id of selected pathway.
#' @param association Dataframe with 2 columns, where each line reprensents an
#'        association. First column are genes and the sencond column are
#'        metabolites. Only use KEGG Ids.
#' @keywords graph, heatmap, shortestDistance, KEGG
#' @examples heatmapAsso("hsa01100", shinAndAlDF)

heatmapAssoAllMaps <- function(association, distanceDF){

    distanceDF <- getDistanceAssoAllMaps(shinAndAlDF,1)[,c(1,3,5,6)]


    mError2 <-"Sorry, for each pairs of gene/metabolite entered, either the gene
    or the metabolite or both weren't mapped on the selected pathway.
    Thus, no distance was calculated"
    fileName <- "allMapsDistribution.png"

    if(nrow(distanceDF) == 0 ){
        stop(mError2, call. = FALSE)
    }


    distanceDF$Concat_c1_c2 <- do.call(paste, c(distanceDF[c(1,2)], sep = "/"))

    associations <- data.frame("association" = distanceDF$Concat_c1_c2)
    associations <- data.frame("association" =associations
                                            [!duplicated(associations),])
    pathways <- data.frame("pathway" = distanceDF$pathwayId)
    pathways <- data.frame("pathway" = pathways
                               [!duplicated(pathways),])

    dat <- data.frame( Row = rep(as.vector(pathways[,1]), each= nrow(associations)),
                       Col = rep(as.vector(associations[,1]), times= nrow(pathways))

    );
    distances<- as.character();
    isAsso <- as.logical();
    for(i in 1:nrow(dat)){

     asso <- as.character(dat[i,1])
     path <- as.character(dat[i,2])

     dist <- distanceDF[distanceDF[,4]==asso & distanceDF[,5]==path,"distance"]

     if(length(dist)==0) dist <- "NaN";
     if(dist == Inf){
         asso <- TRUE
         isAsso <- append(isAsso, asso)
     }else{
         asso <- FALSE
         isAsso <- append(isAsso, asso)
         }
      distances <- append(distances, dist)

    }

    dat$distance <- as.vector(distances)
    dat$Associations <- as.vector(isAsso)

    frames = dat[dat$Associations, c("Row","Col")]

    frames$Row = as.integer(frames$Row)
    frames$Col = as.integer(frames$Col)

    dat$distance <- as.numeric(as.character(dat$distance))



    colors <- c("#BD0026","#E31A1C","#FC4E2A","#FD8D3C","#FEB24C",
                "#FED976", "#FFEDA0", "#FFFFCC")

    p2 <- ggplot2::ggplot(dat, ggplot2::aes(x=Row, y=Col, fill=distance)) +

        ggplot2::theme(
          panel.border = ggplot2::element_rect(colour="black",fill=NA,size=2),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle=315,vjust=1,hjust=0)
         )+

        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label=paste(distance)),na.rm = T) +
        ggplot2::geom_rect(data=frames, size=1, fill=NA, colour="black",
        ggplot2::aes(xmin=Row-0.5, xmax=Row+0.5, ymin=Col-0.5, ymax=Col + 0.5))+

        ggplot2::scale_fill_gradientn(colours = colors, na.value = "white")

    p2 <- p2 + ggplot2::ggsave(fileName, dpi=300)
       # ggplot2::ggtitle(title)
    print(p2);

}

exportDFAllMaps <- function(association, distanceDF){

    # distanceDF <- getDistanceAssoAllMaps(shinAndAlDF,1)[,c(1,3,5,6)]



  #  distanceDF$geneName <- as.vector(geneCommonName)
  #  distanceDF$metaboliteName <- as.vector(metaboliteCommonName)

    print(distanceDF)

    distanceDF_Hsa01100 <- distanceDF;
    distanceDF <- distanceDF_Hsa01100[!distanceDF_Hsa01100$pathwayId == "hsa01100",]

    distanceDF$Concat_c1_c2 <- do.call(paste, c(distanceDF[c(1,2)], sep = "/"))
    distanceDF_Hsa01100$Concat_c1_c2 <-
        do.call(paste, c(distanceDF_Hsa01100[c(1,2)], sep = "/"))



    associations <- data.frame("associations" = distanceDF$Concat_c1_c2)
    associations <- data.frame("associations" =associations
                               [!duplicated(associations),])


    # pathways <- data.frame("pathway" = distanceDF$pathwayId)
    # pathways <- data.frame("pathway" = pathways
    #                        [!duplicated(pathways),])
    #
    # dat <- data.frame( Row = rep(as.vector(pathways[,1]), each= nrow(associations)),
    #                    Col = rep(as.vector(associations[,1]), times= nrow(pathways))

   # );

    #
# print(distanceDF_Hsa01100)
#     distances <- as.character();
#     pathways <- as.character();
#     gene <- as.character()
#     metabo <- as.character()
#
#     for(i in 1:nrow(associations)){
#
#         asso <- as.character(associations[i,1])
#
#
#         dist <- distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"distance"]
#         if(all(dist == dist[1])){
#             dist <- dist[1]
#         }else{
#             dist <- paste(dist, collapse = ',')
#         }
#
#         path <- distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"pathwayId"]
#         ge <- as.character(distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"geneKEGGId"][1])
#
#         meta <- as.character(distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"metaboliteKEGGId"][1])
#
#
#         path <- paste(path, collapse = ',')
#         if(length(dist)==0) dist <- "NaN";
#
#         pathways <- append(pathways, path)
#         distances <- append(distances, dist)
#         gene <- append(gene, ge)
#         metabo <- append(metabo, meta)
#
#     }
#
#
print(associations)
    print(distanceDF_Hsa01100)
    distances <- as.character();
    pathways <- as.character();
    gene <- as.character();
    metabo <- as.character();
    assos_re <- as.character();

    for(i in 1:nrow(associations)){

        asso <- as.character(associations[i,1])


        dist <- distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"distance"]

        if(all(dist == dist[1])){
            dist <- dist[1]
        }else{
            dist <- paste(dist, collapse = ',')


        }
        assos_temp <-asso;
        path <- distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"pathwayId"]
        ge <- as.character(distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"geneKEGGId"][1])

        meta <- as.character(distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"metaboliteKEGGId"][1])
        path <- distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"pathwayId"]
        ge <- as.character(distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"geneKEGGId"][1])

        meta <- as.character(distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"metaboliteKEGGId"][1])


        path <- paste(path, collapse = ',')
        if(length(dist)==0) dist <- "NaN";
        assos_re <- append(assos_re, assos_temp)
        pathways <- append(pathways, path)
        distances <- append(distances, dist)
        gene <- append(gene, ge)
        metabo <- append(metabo, meta)
        # path <- distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"pathwayId"]
        # ge <- as.character(distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"geneKEGGId"][1])
        #
        # meta <- as.character(distanceDF_Hsa01100[distanceDF_Hsa01100[,5]==asso,"metaboliteKEGGId"][1])
        #
        #
        # path <- paste(path, collapse = ',')
        # if(length(dist)==0) dist <- "NaN";
        #
        # pathways <- append(pathways, path)
        # distances <- append(distances, dist)
        # gene <- append(gene, ge)
        # metabo <- append(metabo, meta)

    }





#
      geneCommonName <- getCommonNames(gene, "gene")
# #
     geneCommonName <- as.vector(unlist(geneCommonName))
#     #
#     # print("llll")
#     # print(geneCommonName)
#
#
#
      metaboliteCommonName <- getCommonNames(metabo, "metabolite")
      metaboliteCommonName <- as.vector(unlist(metaboliteCommonName))
    #
    # print("hhhhh")
    # print(metaboliteCommonName)



   # print(associations)
    associat <- data.frame("Associations" =  as.vector(assos_re), "Distances" =  as.vector(distances),
                           "Pathways" =  as.vector(pathways));
   # associat$asso <- as.vector(assos_re)
  #   associat$distances <- as.vector(distances)
  #   associat$pathways <- as.vector(pathways)
  #
  #
  #
  # #  print(gene)
  #  # print(metabo)
     associat$gene <- as.vector(geneCommonName)
      associat$metabolite <- as.vector(metaboliteCommonName)
  #  # print(associations)
  #
      associat$Associations<- do.call(paste, c(associat[c(4,5)], sep = " / "))
      print(associat)

#
 selectedRows <-grep("hsa01100", associat$Pathways)
# print(selectedRows)
 associat <- associat[selectedRows,]
# #rownames(associations_F) <- 1:nrow(associations_F)
 print(associat)
colors <- c("#BD0026","#E31A1C","#FC4E2A","#FD8D3C","#FEB24C",
            "#FED976", "#FFEDA0", "#FFFFCC")

     png("allMapsDF.png", width=10,height=10,units="in", res=1000)

     p<-gridExtra::tableGrob(associat[,c(1,2,3)])
#
       # ind1 <- find_cell(p, 7, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#BD0026",col="#BD0026")


      ind1 <- find_cell(p, 8, 2, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 8, 3, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 8, 4, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")

      ind1 <- find_cell(p, 13, 2, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 13, 3, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 13, 4, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")

      ind1 <- find_cell(p, 15, 2, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 15, 3, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 15, 4, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")

      ind1 <- find_cell(p, 16, 2, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 16, 3, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 16, 4, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")

      ind1 <- find_cell(p, 19, 2, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 19, 3, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 19, 4, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")

      ind1 <- find_cell(p, 20, 2, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 20, 3, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 20, 4, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")

      ind1 <- find_cell(p, 21, 2, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 21, 3, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 21, 4, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 22, 2, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 22, 3, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")
      ind1 <- find_cell(p, 22, 4, "core-bg")
      p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FFA500",col="#FFA500")


       # ind1 <- find_cell(p, 9, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#BD0026",col="#BD0026")
       #
       #
       # ind1 <- find_cell(p, 2, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#E31A1C",col="#E31A1C")
       # ind1 <- find_cell(p, 3, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#E31A1C",col="#E31A1C")
       # ind1 <- find_cell(p, 8, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#E31A1C",col="#E31A1C")
       # ind1 <- find_cell(p, 4, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "#FC4E2A",col="#FC4E2A")
       #
       # ind1 <- find_cell(p, 5, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "grey",col="grey")
       # ind1 <- find_cell(p, 6, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "grey",col="grey")
       # ind1 <- find_cell(p, 10, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "grey",col="grey")
       # ind1 <- find_cell(p, 11, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "grey",col="grey")
       # ind1 <- find_cell(p, 12, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "grey",col="grey")
       # ind1 <- find_cell(p, 13, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "grey",col="grey")
       # ind1 <- find_cell(p, 14, 3, "core-bg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "grey",col="grey")
       #
       # ind1 <- find_cell(p, 4, 4, "core-fg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fontface="bold")
       # ind1 <- find_cell(p, 5, 4, "core-fg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar( fontface="bold")
       # ind1 <- find_cell(p, 8, 4, "core-fg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fontface="bold")
       # ind1 <- find_cell(p, 9, 4, "core-fg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fontface="bold")
       # ind1 <- find_cell(p, 10, 4, "core-fg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fontface="bold")
       # ind1 <- find_cell(p, 11, 4, "core-fg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fontface="bold")
       # ind1 <- find_cell(p, 12, 4, "core-fg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fontface="bold")
       # ind1 <- find_cell(p, 13, 4, "core-fg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fontface="bold")
       # ind1 <- find_cell(p, 14, 4, "core-fg")
       # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fontface="bold")

     # ind1 <- find_cell(p, 5, 2, "core-bg")
     #
     # ind2 <- find_cell(p, 5, 3, "core-bg")
     # ind3 <- find_cell(p, 5, 4, "core-bg")
     # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "darkgoldenrod1",col="darkgoldenrod1")
     # p$grobs[ind2][[1]][["gp"]] <-grid::gpar(fill = "darkgoldenrod1",col="darkgoldenrod1")
     # p$grobs[ind3][[1]][["gp"]] <-grid::gpar(fill = "darkgoldenrod1",col="darkgoldenrod1")
     #
     #
     # ind1 <- find_cell(p, 7, 2, "core-bg")
     # ind2 <- find_cell(p, 7, 3, "core-bg")
     # ind3 <- find_cell(p, 7, 4, "core-bg")
     # p$grobs[ind1][[1]][["gp"]] <-grid::gpar(fill = "darkgoldenrod1",col="darkgoldenrod1")
     # p$grobs[ind2][[1]][["gp"]] <-grid::gpar(fill = "darkgoldenrod1",col="darkgoldenrod1")
     # p$grobs[ind3][[1]][["gp"]] <-grid::gpar(fill = "darkgoldenrod1",col="darkgoldenrod1")
     #


     grid::grid.draw(p)
     dev.off()



}

#
# g <-gridExtra::tableGrob(iris[1:4, 1:3])
find_cell <- function(table, row, col, name="core-fg"){
    l <- table$layout
    which(l$t==row & l$l==col & l$name==name)
}

#
#
# ind <- find_cell(g, 3, 2, "core-fg")
# ind2 <- find_cell(g, 2, 3, "core-bg")
# g$grobs[ind][[1]][["gp"]] <-grid::gpar(fontsize=15, fontface="bold")
# g$grobs[ind2][[1]][["gp"]] <-grid::gpar(fill="darkolivegreen1", col = "darkolivegreen4", lwd=5)
# grid::grid.draw(g)

