#' Globally optimal sums of squares partitioning on phylogeny
#'
#' @param x A phylogenetic tree. An object of class 'phylo'.
#' @param y A data matrix. An object of class 'matrix'. 
#' Row names should be present and are expected to correspond 
#' to the terminal nodes of 'x'.
#' @param k The maximum number of subtrees to partition 'x' into.
#' @param option The global mean to use for partitioning the sums
#' of squares. Either 'arithmetic' or 'phylogenetic'.
#' @value A list with three components:
#' \describe{
#' \item{RSS}{A numeric vector of length k. RSS[i] is the
#' smallest residual sum of squares that can be obtained by
#' partitioning x into i subtrees.}
#' \item{cutset}{An integer matrix with k columns and 
#' ape::Ntip(x) + ape::Nnode(x) rows. cutset[i,j] returns
#' the index of the subtree that node i belongs to when
#' x is partitioned into the j subtrees that yield the
#' smallest residual sum of squares.}
#' \item{weight}{A numeric vector of length ape::Ntip(x). The weights
#' used to partition the sums of squares.}
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
    brlen = x$edge.length
    nedge = length(parent)

    if (is.null(brlen))
        brlen = rep(1, nedge)

    brlen = brlen[match(1:(ntip+nnode), x$edge[,2], nomatch=0)]

    w = numeric(ntip+nnode)
    w[1:ntip] = 1

    option = match.arg(option)

    if (option == "phylogenetic")
        w[1] = -1
    
    ans = .Call(C_cutphylo
        ,k,ntip,nnode,nedge,parent,child,brlen,w,t(y))
    ans[[3]] = w[1:ntip]

    min_index = apply(ans[[1]], 2, which.min)
    ans[[1]] = ans[[1]][cbind(min_index, 1:k)]
    tmp = matrix(0L, ntip+nnode, k)
    for (i in 1:k)
        tmp[,i] = ans[[2]][, min_index[i], i]
    ans[[2]] = tmp
    structure(setNames(ans, c("RSS", "cutset", "weight")), class="cutphylo")
}
