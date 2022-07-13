plot.cutphylo = function(x, ...)
{
    stopifnot(inherits(x, "cutphylo"))
    op = par(no.readonly=TRUE)
    on.exit(par(op))
    par(mfrow=c(1, 3))
    K = length(x[["RSS"]])
    k = 2L:(K-1)
    plot(1:K, x[["RSS"]], las=1, xlab="K", ylab="RSS", type="b", ...)
    plot(1:K, x[["r.squared"]], las=1, xlab="K", ylab="Variance explained", 
        type="b", ylim=c(0,1), ...)
    plot(k, x[["RSS.profile"]], las=1, xlab="K", 
        ylab="RSS profile", type="b", ...)
}
