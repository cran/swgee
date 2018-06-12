print.summary.swgee <-
function(x, ...)
{
    cat("Call: beta \n")
    printCoefmat(x$beta, P.values=TRUE, has.Pvalue=TRUE)
    cat("\n")
    cat("Call: alpha \n")
    printCoefmat(x$alpha, P.values=TRUE, has.Pvalue=TRUE)
}
