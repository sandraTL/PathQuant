
\name{getDistanceAll}
\alias{getDistanceAll}
\title{Function calculating shortest distance between every gene in a
gene-metabolite pairs and all metabolite.}
\usage{
getDistanceAll(pathwayId,gene, metabolite)
}
\arguments{
\item{pathwayId}{KEGG Id of selected map}

\item{gene}{Dataframe of 1 column, representing all the genes reported. KEGG Ids.}

\item{metabolite}{Dataframe of 1 column, representing all the measured metabolites. KEGG Ids}
}
\description{
Function calculating shortest distance between every gene in a
gene-metabolite pairs of your data parameter and every metabolite
(in a pair or not) on a graph model of KEGG map selected, where nodes are
metabolites and reactions are edges.
}
\details{
If a gene or a metabolite is present on multiple edges or nodes, then
shortest distance are calculated for every combinaison possible and the
shortest distance is selected.

Output : dataframe with metabolite in columns and gene in rows with
         shortest distance values.
}
\examples{
getDistanceAll("hsa01100",completeGeneDF,completeMetaboDF)
}
\keyword{KEGG}
\keyword{graph,}
\keyword{shortestDistance,}

