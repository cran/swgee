getsimexest <-
function(indata) 
{
    lambda <- sort(unique(indata$lambda))
    length.beta<- (dim(indata)[2]-1)/2

    locdata <- lapply(lambda, function(lamk, indata){
        tmpdata <- indata[indata$lambda == lamk,]
        beta <- apply(tmpdata, 2, mean)[2:(length.beta+1)]
        tau <- apply(tmpdata, 2, mean)[(length.beta+2):(dim(indata)[2])] - apply(tmpdata, 2, var)[2:(length.beta+1)]
        return(c(lambda=lamk, beta, tau))
        }, indata=indata)
    
    locdata <- do.call("rbind", locdata)
    locdata=as.data.frame(locdata)
    locdata <- locdata[order(locdata$lambda),]

    evalest <- lapply(2:dim(locdata)[2], function(ind){
                    lambda2=lambda*lambda
                    fit <- lm(locdata[,ind] ~ lambda + lambda2)
                    fitcoef <- sum(fit$coef*c(1, -1, 1))
                    return(fitcoef)
                    })
    
    evalest <- do.call("cbind", evalest)
    colnames(evalest)=colnames(indata)[-1]
    evalest=as.data.frame(evalest)
    return(list(locdata,evalest))
}
