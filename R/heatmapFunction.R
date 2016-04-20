
#' Function that output a heatmap to visualize distance calculated between every
#' gene-metabolite associations in input.
#'
#' Function calculting shortest distance between every genes and metabolites in
#' a gene-metabolite pairs of your data parameter on a graph model of KEGG map
#' selected, where nodes are metabolites and reactions are edges.
#'
#' If a gene or a metabolite is present on multiple edges or nodes, then
#' shortest distance are calculated for every combinaison possible and the
#' shortest distance is selected.
#'
#' Output : Heatmap graphic where black countoured distances represent
#'           an association.
#'
#' @param pathwayId KEGG Id of selected pathway.
#' @param association Dataframe with 2 columns, where each line reprensents an
#'        associations. First column are the genes and the sencond column as the
#'        metabolites. Only use KEGG Ids.
#' @keywords graph, heatmap, shortestDistance, KEGG
#' @export
#' @examples heatmapAsso(metabolismOverviewMapKEGGId, shinAndAlDF)

heatmapAsso <- function(pathwayId, data){
    pathwayId <- gsub("hsa:", "hsa", pathwayId)
    mError1 <-"Error in input associatedGeneMetaboDF, please enter you data
             where colnames(df) <- c(gene,metabolite) frame with
             KEGG ids of genes (ex : hsa:00001) in first
             column and associated KEGG ids metabolites (ex: C00001)
             in second column"
    mError2 <-"Sorry, for each pairs of gene/metabolite entered, either the gene
             or the metabolite or both weren't mapped on the selected pathway.
             Thus, no distance was calculated"

    if(length(data) == 0){
        stop(mError1, call. = FALSE);
    }
    if(!is.data.frame(data) ||
       length(data[1,])< 2 ||
       length(data[1,])> 3){
        stop(mError1, call. = FALSE);
    }
    for(row in 1:nrow(data)){

        if(substr(data[row,1],1,4)!="hsa:")
            stop(mError1, call. = FALSE);
        if(substr(data[row,2],0,1) != "C"
           && length(data[row,2]) != 5)
            stop(mError1, call. = FALSE);
    }

    graphe <-  createGraphFromPathway(pathwayId);
    rGeneList<-numberOfReactions(graphe@edgeDF,data[,1])
    rMetaboliteList <- numberOfMetabolites(graphe@nodeDF, data[,2])
    tempDf1 <- data.frame(cbind(g1 = as.vector(data[,1]),
                                g2 = as.vector(as.numeric(rGeneList)),
                                m1 = as.vector(data[,2]),
                                m2 = as.vector(as.numeric(rMetaboliteList))))

    tempDf1 <- removeNotInGraph(tempDf1)

    if(nrow(tempDf1) == 0 ){
        stop(mError2, call. = FALSE)
    }
    data <- subset(tempDf1[,c(1,3)])

    data1 <- data.frame(c(data[2]))

    AllSP <- getDistanceAll(pathwayId, data, data1);

    geneCommonName <- getCommonNames(as.vector(unlist(rownames(AllSP))), "gene")
    geneCommonName <- as.vector(unlist(geneCommonName))

    metaboliteCommonName <- getCommonNames(as.vector(unlist(colnames(AllSP))),
                                           "metabolite")
    metaboliteCommonName <- as.vector(unlist(metaboliteCommonName))
    #print(cbind(AllSP, geneCommonName, metaboliteCommonName))
    tempDat <- data.frame(Row = rep(rownames(AllSP), each= ncol(AllSP)),
                          Col = rep(colnames(AllSP), times= nrow(AllSP)))
    isAsso <- getAssociationForHeatmap(data, tempDat)



    dat <- data.frame( Row = rep(geneCommonName, each= ncol(AllSP)),
                       Col = rep(metaboliteCommonName, times= nrow(AllSP)),
                       Distance = c(t(AllSP)),
                       Associations = isAsso

    );
    # create frame to color edges of associated genes and metabolites
    frames = dat[dat$Associations, c("Row","Col")]
    frames$Row = as.integer(frames$Row)
    frames$Col = as.integer(frames$Col)

    # from discrete to continuous values...
    dat$Distance <- as.numeric(as.character(dat$Distance))


    colors <- c("#800026","#BD0026","#E31A1C","#FC4E2A","#FD8D3C","#FEB24C",
                "#FED976", "#FFEDA0", "#FFFFCC")
    p2 = ggplot2::ggplot(data=dat) +
        ggplot2::geom_raster(ggplot2::aes(x=Row, y=Col, fill=Distance)) +
        ggplot2::theme(
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=2),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            text = ggplot2::element_text(size=12),
            axis.text=ggplot2::element_text(colour="black"),
            axis.text.x = ggplot2::element_text(angle=315,vjust=1,hjust=0))+
        ggplot2::scale_fill_gradientn(colours = colors)+
        ggplot2::xlab("Genes")+
        ggplot2::ylab("Metabolites")+
        ggplot2::geom_rect(data=frames, size=1, fill=NA, colour="black",
        ggplot2::aes(xmin=Row-0.5, xmax=Row+0.5, ymin=Col-0.5, ymax=Col + 0.5))+
        ggplot2::geom_text(label = as.numeric(dat$Distance, 1),
    size = 2,ggplot2::aes(x = Row, y = Col)) +



    ggplot2::ggsave("heatmapAsso1.png", width = 9, height = 7, dpi=300)
    print(p2);

}

getAssociationForHeatmap<- function(data, heatmapDf){
    boolIsAsso <- FALSE;
    f <- apply(heatmapDf, 1, function
               (x){
                   isAssociation <-  length(data[data[1] == x[1] & data[2] == x[2]])
                   if(isAssociation >= 2){
                       boolIsAsso <- TRUE;
                   }else{
                       boolIsAsso<- FALSE
                   }
               })
}


