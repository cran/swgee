print.swgee <-
function(x, ...)
{  
    cat("Formula: \n")
    print(x$formula)
    cat("\n")
    cat("SIMEX variables: \n")
    print(x$SIMEXvariable)
    cat("\n")
    cat("Number of simulations: \n")
    print(x$B)
    cat("\n")
    cat("This is response process coefficients: \n")
    print(x$beta)
    cat("\n")
    cat("This is missing process coefficients:\n")
    print(x$alpha)
}
