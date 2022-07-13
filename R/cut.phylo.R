#' Sums of squares partitioning on a phylogeny
#'
#' @param x A phylogenetic tree. An object of class 'phylo'.
#' @param y A data matrix. An object of class 'matrix'. 
#' Row names should be present and are expected to correspond 
#' to the terminal nodes of 'x'.
#' @param k The maximum number of subtrees to partition 'x' into.
#' @param option The global mean to use for partitioning the sums
#' of squares. Either 'arithmetic' or 'phylogenetic'.
#' @details The function tries to find a set of \code{k-1} branches 
#' that, when cut, partition the phylogeny into \code{k} subtrees
#' that best explain the variance of \code{y}. Currently, a greedy
#' search algorithm is used that incrementally builds a solution
#' from the previous best solution. I.e., the solution for k=3 uses 
#' the solution for k=2 and so on.
#' @value A list with components:
#' \describe{
#' \item{RSS}{A numeric vector of length k. RSS[i] is the
#' smallest residual sum of squares that can be obtained by
#' partitioning x into i subtrees.}
#' \item{RSS.profile}{A numeric vector of length k-1 formed by
#' (RSS[i-1] - RSS[i]) / (RSS[i] - RSS[i+1]) for i=2,...,k-1.}
#' \item{r.squared}{A numeric vector of length k. The variance
#' explained by each partition.}
#' \item{subtree.r.squared}{A numeric k-by-k matrix. The variance
#' explained by each partition subtree in each partition.}
#' \item{weights}{A numeric vector of length ape::Ntip(x). The weights
#' used to partition the sums of squares.}
#' \item{cuts}{An integer vector with the indices of the branches that 
#' were cut.}
#' \item{cutset}{An integer matrix with k columns and 
#' ape::Ntip(x) + ape::Nnode(x) rows. cutset[i,j] returns
#' the index of the subtree that node i belongs to when
#' x is partitioned into the j subtrees that yield the
#' smallest residual sum of squares.}
#' }
cut.phylo = function(x, y, k, option=c("arithmetic","phylogenetic"), ...)
{
    stopifnot(inherits(x, "phylo"))
    stopifnot(inherits(y, "matrix"))
    stopifnot(!is.null(rownames(y)))
    stopifnot(ape::is.binary(x))
    stopifnot(all(rownames(y) %in% x$tip.label))
    stopifnot(all(x$tip.label %in% rownames(y)))
    ntip = ape::Ntip(x)
    nnode = ape::Nnode(x)
    k = max(1L, min(ntip, as.integer(k)))
    x = reorder(x, "pruningwise")
    y = y[x$tip.label, ]
    parent = x$edge[,1]-1L
    child = x$edge[,2]-1L
    edgelen = x$edge.length
    nedge = length(parent)

    if (is.null(edgelen))
        edgelen = rep(1, nedge)

    ord = match(1:(ntip+nnode), x$edge[,2], nomatch=0)
    brlen = numeric(ntip+nnode)
    brlen[-(ntip+1)] = edgelen[ord]

    w = numeric(ntip+nnode)
    w[1:ntip] = 1

    option = match.arg(option)

    if (option == "phylogenetic")
        w[1] = -1

    Y = t(y)
    Y = cbind(Y, matrix(0, nrow(Y), nnode))
    
    ans = .Call(C_cutphylo
        ,k,ntip,nnode,nedge,parent,child,brlen,w,Y)

    SSE = ans[[1]]
    cuts = ans[[2]]
    cutset = ans[[3]]
    SSE.profile=setNames(
            (SSE[1:(k-2)]-SSE[2:(k-1)])/(SSE[2:(k-1)]-SSE[3:k]), 2:(k-1))

    w = w[1:ntip]
    GM = colMeans(sweep(y, 1, w, "*"))
    R2 = matrix(0, k, k, dimnames=list(cuts,NULL))
    for (i in 2:k)
    {
        for (j in 1:i)
        {
            idx = which(cutset[1:ntip,i] == j)
            m = colSums(sweep(y[idx,,drop=FALSE], 1, w[idx], "*")) / sum(w[idx])
            R2[j, i] = sum(w[idx]) * sum((m - GM)^2) / SSE[1]
        }
    }

    structure(list(
        RSS=SSE
        , RSS.profile=SSE.profile
        , r.squared=1-SSE/SSE[1]
        , subtree.r.squared=R2
        , weights=w
        , cuts=cuts
        , cutset=cutset
    ), class="cutphylo")

}
