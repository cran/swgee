swgee <-
function(formula, data = parent.frame(), id, family=family, corstr = "independence", missingmodel,
 SIMEXvariable, SIMEX.err, repeated = FALSE, repind = NULL, B = 50, lambda = seq(0, 2, 0.5))
{
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$link <- m$family <- m$corstr <- m$missingmodel <- m$SIMEXvariable <- m$SIMEX.err <- m$repeated <- m$repind <- m$B <- m$lambda <- NULL
    if(is.null(m$id))
        m$id <- as.name("id")
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    id <- model.extract(m, id)
    iid <- as.numeric(names(id))
    namesx <- colnames(model.matrix(Terms, m, contrasts))
    vars <- all.vars(formula)
    mmvars <- all.vars(missingmodel)[-1]
    names.vars.lag <- paste(mmvars,".lag", sep="")
    if (is.null(id) )
        stop( "id variable not found.", call. = FALSE )
    if (is.null(missingmodel)) {
        warning( "Missing model should be specified", call. = FALSE )
    }
    vv <- NULL
    nSIMEXvariable = length(SIMEXvariable)
    for ( i in seq_along(SIMEXvariable) ) {
         if((!is.na(match(SIMEXvariable[i], vars))))
         vv <- c(vv,i)
         }
    if ( is.null(vv)|length(vv) != nSIMEXvariable )
        stop( "SIMEXvariable must be character and specified in the formula", call. = FALSE )
    if (repeated == FALSE) {
        SIMEX.err <- as.matrix(SIMEX.err)
        if ( nrow(SIMEX.err) != ncol(SIMEX.err) ) 
            stop( "SIMEX.err must be a square matrix", call. = FALSE )
        if ( length(SIMEXvariable) != nrow(SIMEX.err) ) 
            stop( "SIMEXvariable and SIMEX.err have non-conforming size", call. = FALSE )
    }
    else if (repeated == TRUE) {
        if (length(SIMEXvariable) != length(repind)) {
            stop("SIMEXvariable and repind have non-conforming size")
        }
    }
    if ( !is.numeric(B) | B <= 0 ) {
        stop( "B must be positive integer", call. = FALSE )
    }
    else {
        B <- ceiling(B)
    }
    if ( !is.vector(lambda) | !is.numeric(lambda) ) 
        stop( ":Invalide lambda value", call. = FALSE )
    if ( any(lambda < 0) ) {
        warning( "Lambda should be positive values. Negative values will be ignored", call. = FALSE )
        lambda <- lambda[lambda >= 0]
    }
    
    locdata <- data   
    R <- rep(0,max(iid))
    R[iid] <- 1
    max.visit <- (max(iid)/length(unique(id)))-1
    locdata$R <- R
    
    alphainfo <- betainfo <- NULL
    nlambda  <- length(lambda)
    for (b in 1:B) {
        alphainfo.B <- c()
        betainfo.B <- c()
        wb <- list()
        for (k in 1:nlambda) {
            if (repeated == FALSE) {
                wb[[k]] <- locdata[SIMEXvariable]
                if ( lambda[k] > 0 ) {
                    wb[[k]] <- locdata[SIMEXvariable] + sqrt(lambda[k]) * rmvnorm(nrow(locdata), rep(0, nSIMEXvariable), SIMEX.err)
                }
            }
            else if (repeated == TRUE) {
                sudo.simvar <- locdata[SIMEXvariable]
                for (i in 1:nSIMEXvariable) {                  
                    n.i <- length(repind[[i]])
                    mat <- lapply(1:nrow(locdata), function(obs, inx) {
                        z.i = rnorm(n.i, 0, 1)
                        constrast.i = (z.i - mean(z.i))/(sd(z.i)*sqrt(n.i-1))
                        sudo.i = mean(as.numeric(inx[obs,])) + sqrt(lambda[k]/n.i)*as.numeric(inx[obs,])%*%constrast.i
                        return(sudo.i)
                        }, inx=locdata[,repind[[i]]])
                   sudo.simvar[SIMEXvariable[i]]=as.vector(do.call("rbind", mat))
                }
                wb[[k]] <- sudo.simvar
            }
        }
        for (k in 1:nlambda){
            inlocdata <- locdata
            inlocdata[SIMEXvariable] <- wb[[k]]
            #names(inlocdata.lag) <- names.vars.lag
            varlag <- function(id, inx, lag)
            {
                lag.var <- lapply(unique(id), function(obs) {
                id.min <- min(iid[which(id == obs)])
                lag.inx <- inx[id.min:(id.min+max.visit)]
                c(rep(NA, lag), lag.inx[1:(length(lag.inx) - lag)])
                })
                lag.var <- unlist(lag.var)
                return(lag.var)
            }
            for( i in 1: length(mmvars)){
                inlocdata[names.vars.lag[i]] <- varlag(id=id, inx=inlocdata[,mmvars[i]], lag=1)
            }

            ######### gee coef Robust S.E. #####################
            miss.formula <- paste("R ~ ", paste(names.vars.lag, collapse = "+"), sep = "")
            miss.fit <- glm(miss.formula, data = inlocdata, family=quasibinomial(link="logit"))
            alpha <- as.numeric(summary(miss.fit)$coefficients[ ,1])
            alphase <- as.numeric(summary(miss.fit)$coefficients[ ,2])
            alphainfo.k <- c(lambda[k], alpha, alphase^2)
            lambdaR <- predict(miss.fit, newdata = inlocdata, type = "response")
            inlocdata$wgtprob <- inlocdata$R
            for(obs in 1:length(unique(id))) {
                index <- as.numeric(names(which(id==unique(id)[obs])))
                prob <- lambdaR[index]
                prob[1] <- 1
                pprob <- unlist(lapply(1:length(index), function(m) {
                    prod(prob[1:m])
                }))
                inlocdata$wgtprob[index] <- R[index]/pprob
            }
            inlocdata <- inlocdata[inlocdata$R==1,]
            inlocdata[vars] <- lapply(inlocdata[vars], as.numeric)
            inlocdata <- inlocdata[, !(names(inlocdata) %in% c(names.vars.lag,"R"))]
            inlocdata$id <- id
            inlocdata <- as.data.frame(inlocdata)
            wgtprob <- NULL
            fit <- suppressWarnings(geeglm(formula, id=id, weights=wgtprob, data=inlocdata, family=family, corstr=corstr))
            beta <- as.numeric(fit$geese$beta)
            varbeta <- diag(fit$geese$vbeta)
            betainfo.k <- c(lambda[k], beta, varbeta)
            betainfo.B <- rbind(betainfo.B, betainfo.k)
            alphainfo.B <- rbind(alphainfo.B, alphainfo.k)
            }
        alphainfo <- rbind(alphainfo, alphainfo.B)
        betainfo <- rbind(betainfo, betainfo.B)
    }    
    alphainfo <- data.frame(alphainfo, row.names = NULL)
    betainfo <- data.frame(betainfo, row.names = NULL)
    colnames(betainfo) <- c("lambda", namesx, paste("var.", namesx, sep=""))
    colnames(alphainfo) <- c("lambda", paste("alpha", 1:length(alpha), sep=""), paste("var.alpha", 1:length(alpha), sep=""))
    alpha.est <- getsimexest(indata=alphainfo)
    beta.est <- getsimexest(indata=betainfo)
    plot.est <- rbind(c(-1, unlist(beta.est[[2]][ , namesx[-1]])), beta.est[[1]][,c("lambda", namesx[-1])])
    plot.est <- data.frame(plot.est)
    value <- list()
    value$call <- call
    value$family <- family
    value$corstr <- corstr
    value$formula <- formula
    value$SIMEXvariable <- SIMEXvariable
    value$B <- B
    value$beta <- unlist(beta.est[[2]])
    value$alpha <- unlist(alpha.est[[2]])
    value$simex.plot <- plot.est
    class(value) <- "swgee"
    return(value)
}
