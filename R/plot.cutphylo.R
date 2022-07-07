plot.cutphylo = function(x, ...)
{
    stopifnot(inherits(x, "cutphylo"))
    op = par(no.readonly=TRUE)
    on.exit(par(op))
    par(mfrow=c(1, 3))
    Q = x[[1]]
    K = length(Q)
    k = 2L:(K-1)
    plot(1:K, Q, las=1, xlab="K", ylab="RSS", type="b", ...)
    plot(1:K, 1-Q/Q[1], las=1, xlab="K", ylab="Variance explained", type="b", 
        ylim=c(0,1), ...)
    plot(k, (Q[k-1]-Q[k]) / (Q[k]-Q[k+1]), las=1, xlab="K", 
        ylab="RSS profile", type="b", ...)
}
