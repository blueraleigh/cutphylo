\name{cut.phylo}
\alias{cut.phylo}
\title{Sums of squares partitioning on a phylogeny}
\usage{
cut(x,y,k,option)
}
\arguments{
\item{x}{A phylogenetic tree. An object of class "phylo".}
\item{y}{A data matrix. An object of class "matrix". Row names 
should be present and are expected to correspond to the terminal 
nodes of \code{x}.}
\item{k}{An integer. The maximum number of subtrees to partition \code{x} into.}
\item{option}{The global mean to use for partitioning the sums of squares. 
Either "arithmetic" or "phylogenetic". Default: "arithmetic"}
}
\details{
The function tries to find a set of \code{k-1} branches 
that, when cut, partition the phylogeny into \code{k} subtrees
that best explain the variance of \code{y}. Currently, a greedy
search algorithm is used that incrementally builds a solution
from the previous best solution. I.e., the solution for k=3 uses 
the solution for k=2 and so on.

When \code{option="phylogenetic"}, the function partitions a
weighted sum of squares by computing a set of tip weights, \code{w},
such that \code{sum(w*y)/sum(w)} is equal to the phylogenetic mean
of \code{y} under a Brownian motion model. These weights are returned
as part of the result.
}
\value{
A list with components:
\describe{
\item{RSS}{A numeric vector of length \code{k}. \code{RSS[i]} is the smallest 
residual sum of squares that can be obtained by partitioning \code{x} 
into \code{i} subtrees.}
\item{RSS.profile}{A numeric vector of length \code{k-1} formed by
\code{(RSS[i-1] - RSS[i]) / (RSS[i] - RSS[i+1])} for i=2,...,k-1.}
\item{r.squared}{A numeric vector of length \code{k}. The variance
explained by each partition.}
\item{subtree.r.squared}{A numeric \code{k}-by-\code{k} matrix. The variance
explained by each partition subtree in each partition.}
\item{weights}{A numeric vector of length \code{ape::Ntip(x)}. The weights
used to partition the sums of squares.}
\item{cuts}{An integer vector with the indices of the branches that 
were cut.}
\item{cutset}{An integer matrix with \code{k} columns and 
\code{ape::Ntip(x) + ape::Nnode(x)} rows. The value returned by 
\code{cutset[i,j]} is the index of the subtree that node \code{i} 
belongs to when \code{x} is partitioned into the \code{j} subtrees 
that yield the smallest residual sum of squares.}
}
}
