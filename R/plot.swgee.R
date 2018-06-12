plot.swgee <-
function(x, covariate, ...)
{  
    if(is.na(match(as.character(covariate), all.vars(x$formula))))
        stop("covariate not specified in the formula")
    
    plot(x$simex.plot[,covariate] ~ x$simex.plot$lambda, type="b", xlab = expression(lambda), ylab=as.character(covariate))

}
