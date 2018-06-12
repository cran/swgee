summary.swgee <-
function(object, ...)
{
    beta <- matrix(object$beta, ncol = 2)
    est.beta <- unlist(beta[,1])
    names(est.beta) <- names(object$beta)[1:(length(object$beta)/2)]
    std.beta <- sqrt(unlist(beta[,2]))
    tstat.beta <- est.beta/std.beta
    pvalue.beta <- 1 - pchisq(tstat.beta*tstat.beta, 1)
    
    beta <- cbind(Estimate = est.beta,
                  StdErr = std.beta,
                  t.value = tstat.beta,
                  p.value = pvalue.beta)  
    
    alpha <- matrix(object$alpha, ncol = 2)
    est.alpha <- unlist(alpha[,1])
    names(est.alpha) <- names(object$alpha)[1:(length(object$alpha)/2)]
    std.alpha <- sqrt(unlist(alpha[,2]))
    tstat.alpha <- est.alpha/std.alpha
    pvalue.alpha <- 1 - pchisq(tstat.alpha*tstat.alpha, 1)
    
    alpha <- cbind(Estimate = est.alpha,
                  StdErr = std.alpha,
                  t.value = tstat.alpha,
                  p.value = pvalue.alpha) 
                   
    out <- list(beta=beta, alpha=alpha)
    class(out) <- "summary.swgee"
    out
}
