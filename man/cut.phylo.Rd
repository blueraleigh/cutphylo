\name{cut.phylo}
\alias{cut.phylo}
\title{Globally optimal sums of squares partitioning on phylogeny}
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
\value{
A list with three components:
\describe{
\item{RSS}{A numeric vector of length \code{k}. \code{RSS[i]} is the smallest 
residual sum of squares that can be obtained by partitioning \code{x} 
into \code{i} subtrees.}
\item{cutset}{An integer matrix with \code{k} columns and 
\code{ape::Ntip(x) + ape::Nnode(x)} rows. The value returned by 
\code{cutset[i,j]} is the index of the subtree that node \code{i} 
belongs to when \code{x} is partitioned into the \code{j} subtrees 
that yield the smallest residual sum of squares.}
\item{weight}{A numeric vector of length \code{ape::Ntip(x)}. The weights
used to partition the sums of squares.}
}
}
