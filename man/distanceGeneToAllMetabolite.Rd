
\name{distributionGene}
\alias{distributionGene}
\title{distributionGene - Distance distribution plots for single gene}
\usage{
distributionGene(pathwayId, associations, metabolite, gene)
}
\arguments{
\item{pathwayId}{KEGG Id of selected pathway}

\item{associations}{Dataframe with 2 columns, where each line reprensents an associations. First column are the genes and the sencond column as the metabolites. Only use KEGG Ids.}

\item{metabolite}{Dataframe of 1 column, representing all the measured metabolites. KEGG Ids}

\item{gene}{Selected gene. KEGG Id}
}
\description{
Function ploting the distribution of distances between a gene and all measured metabolite, highlithing the distance of its associated metabolites.
}
\details{
Function calculating shortest distance between a selected gene and all measured metabolites.The plot is depicted as frequency bars, which represent the number of metabolites at a given distance for the selected gene. Frequency bars are shown in grey for metabolites that are not associated with the selected gene and in red if there is at least one metabolite associated with this gene.

If a gene or a metabolite is present on multiple edges or nodes, then the shortest distance is selected.
}
\examples{
distributionGene("hsa01100",
shinAndAlDF, completeMetaboDF, "hsa:1373")
}
\keyword{KEGG}
\keyword{graph,}
\keyword{shortestDistance,}

