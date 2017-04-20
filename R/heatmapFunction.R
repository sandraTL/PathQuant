

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

    # Bind with heatmap dataset
    df.heatmap <- cbind(df.heatmap, "Associations" = as.vector(is.asso))

    #Compute distance for heatmap dataset
    heatmap.result.ob <- create.data.restul.ob(df.heatmap,
                                               pathway = c("hsa01100"),
                                               F,
                                               F,
                                               F)

    # Bind results
    df.heatmap <- cbind(df.heatmap,
                        "Distance" = as.vector(heatmap.result.ob@distance))
    colnames(df.heatmap) <- c("Row", "Col", "Associations", "Distance")

    if(commonNames){

        # Get commun names for gene and metabolite
        gene.names <- getCommonNames(as.vector(unlist(gene.uniq)), "gene")
        gene.names <- as.vector(unlist(gene.names))

        metabolite.names <-
            getCommonNames(as.vector(unlist(metabolite.uniq)), "metabolite")
        metabolite.names <- as.vector(unlist(metabolite.names))

        #Change ids to common names
        df.heatmap$Row <- rep(gene.names, each = length(metabolite.names))
        df.heatmap$Col <- rep(metabolite.names, times = length(gene.names))

        # Order in alphabetical order before defining the frame ids
        # because the heatmap function orders automatically
        df.heatmap <- df.heatmap[order(df.heatmap[1], df.heatmap[2]),]

        # Set the number of rows in increasing order to built frame id
        # for associations
        row.names(df.heatmap) <- c(1:nrow(df.heatmap))

    }

    # Define frame ids ready for black highlight of associations
    frame = df.heatmap[df.heatmap$Associations, c("Row","Col")]
    asso <- as.numeric(row.names(frame))
    row.asso.frame <- lapply(asso, function(x) ceiling(x/length(metabolite.uniq)))
    col.asso.frame <- lapply(asso, function(x) {
        if(x%%(length(metabolite.uniq)) == 0){
            length(metabolite.uniq)
        } else {
            x%%(length(metabolite.uniq))
        }
    }
    )

    frame <- data.frame("Row" = as.vector(unlist(row.asso.frame)),
                        "Col" = as.vector(unlist(col.asso.frame)))
    frame$Row = as.integer(frame$Row)
    frame$Col = as.integer(frame$Col)

    # Set heat colors
    colors <- c("#BD0026","#E31A1C","#FC4E2A","#FD8D3C","#FEB24C",
                "#FED976", "#FFEDA0", "#FFFFCC")

    # ggplot2 heatmap preparation
    p2 = ggplot2::ggplot(df.heatmap, ggplot2::aes(x = Row, y = Col, fill = Distance),
                         binwidth = c(1,1)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradientn(colours = colors) +
        ggplot2::theme(
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=2),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),

            text = ggplot2::element_text(size=12),
            axis.text = ggplot2::element_text(colour="black"),
            axis.text.x = ggplot2::element_text(angle=315,vjust=1,hjust=0))+

        ggplot2::coord_equal() +
        ggplot2::xlab("Genes")+
        ggplot2::ylab("Metabolites")+
        ggplot2::ggtitle("Heatmap")+
        ggplot2::guides(fill=ggplot2::guide_legend(title="srd"))+
        ggplot2::coord_fixed(ratio=0.5)+
        ggplot2::geom_rect(data=frame, size=1, fill=NA, colour="black",
                           ggplot2::aes(xmin=Row-0.5, xmax=Row+0.5, ymin=Col-0.5, ymax=Col + 0.5))+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
        ggplot2::geom_text(label = as.numeric(df.heatmap$Distance, 1),
                           size = 2,ggplot2::aes(x = Row, y = Col))
    print(p2)

    return <- df.heatmap
}

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
