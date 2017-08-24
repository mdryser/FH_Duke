##################################################
# Full Likelihood with 4 parameters
# Use full likelihood to estimate 4 parameters
# psi:    proportion of non-progressive cancers
# lambda: exponential rate of progression to clinical
#         state for progressive cancers
# w:      incidence of preclinical disease (N->P/I)
# beta:   screening sensitivity
# 
# Hardcoded parameter inputs:
#   Time of onset of susceptibility to P/I 
#   Screening times
##################################################
EstPsiBMucstab4Full.fun <- function(x, y, t, z, clin.FU) {
    # x: vector of parameters (psi, lambda, w, beta)
    # y: rows are: y[1:3, ]=matrix(s[1], r[1], n[1],
    #                              s[2], r[2], n[2],
    #                              ..., ..., ...)
    #    where s is the number of screen detections and
    #          r is the number of interval cases
    #          n is the number of participants
    # z: if z[1] > 0: it is the number of parameters to be fixed
    #       z[2]: value of the parameter
    if(sum(z[[1]]) > 0){
        theta <- rep(-1, 4)
        theta[z[[1]]] <- z[[2]]
        theta[theta == -1] <- x
    } else {
        theta <- x
    }
    y40 <- y[1:5, ]
    y50 <- y[6:10, ]
    n.y <- nrow(y)/2
    t40 <- t[1:7]
    t50 <- t[8:14]
    ll <- length(clin.FU)/2
    clin.FU40 <- clin.FU[1:ll]
    clin.FU50 <- clin.FU[(ll+1):(2*ll)]

    ##################################################
    # Calculate probability of screen detection at
    # each exam
    ##################################################
    v <- t(0:(n.y-1))
    d40 <- apply(v, 2, function(x) D.j(theta, t40, x))
    d50 <- apply(v, 2, function(x) D.j(theta, t50, x))

    ##################################################
    # Calculate probability of clinical diagnosis in
    # each interval
    ##################################################
    q40 <- apply(v, 2, function(x) I.j(theta, t40, x))
    q50 <- apply(v, 2, function(x) I.j(theta, t50, x))

    ##################################################
    # Calculate and return full log likelihood
    ##################################################
    #n.f40<- CF.fun(theta, t40[2]-t40[1]) # need to correct for the fact that we assumed there was no incidence in first 20 years
    #n.f50<- CF.fun(theta, t50[2]-t50[1]) # need to correct for the fact that we assumed there was no incidence in first 20 years
    
    # NEW version with exact normalization terms: 40
    n.f40 <- rep(CF.fun(theta, t40[2]-t40[1]), n.y)
    QQ40 <- cumsum(d40+ q40 )
    n.f40[2:n.y] <- n.f40[2:n.y]-QQ40[1:(n.y-1)]
    
    # NEW version with exact normalization terms: 50
    n.f50 <- rep(CF.fun(theta, t50[2]-t50[1]), n.y)
    QQ50 <- cumsum(d50+ q50 )
    n.f50[2:n.y] <- n.f50[2:n.y]-QQ50[1:(n.y-1)]
    
    
    fin <- 0

    # step 1: all but the last round
    for(jj in 1:(n.y-1)){
        fin<- fin -(y40[jj, 1]*log(d40[jj]/n.f40[jj])+y40[jj, 2]*log(q40[jj]/n.f40[jj])+(y40[jj, 3]-y40[jj, 1]-y40[jj, 2])*log(1-d40[jj]/n.f40[jj]-q40[jj]/n.f40[jj]))
        fin<- fin -(y50[jj, 1]*log(d50[jj]/n.f50[jj])+y50[jj, 2]*log(q50[jj]/n.f50[jj])+(y50[jj, 3]-y50[jj, 1]-y50[jj, 2])*log(1-d50[jj]/n.f50[jj]-q50[jj]/n.f50[jj]))
    }

    # step 2: the last round
    # The screen-detection contribution
    fin <- fin - y40[n.y, 1]*log(d40[n.y]/n.f40[n.y])
    fin <- fin - y50[n.y, 1]*log(d50[n.y]/n.f50[n.y])

    # The annual clinical incidence contributions
    p40 <- apply(t(1:ll), 2, function(x) Ifin.l(theta, t40, x))
    p50 <- apply(t(1:ll), 2, function(x) Ifin.l(theta, t50, x))
    fin <- fin - clin.FU40 %*% log(p40/n.f40[n.y])
    fin <- fin - clin.FU50 %*% log(p50/n.f50[n.y])

    # The remaining contribution
    fin <- fin-(y40[n.y, 3]-y40[n.y, 1]-sum(clin.FU40))*log(1-d40[n.y]/n.f40[n.y]-sum(p40)/n.f40[n.y])
    fin <- fin-(y50[n.y, 3]-y50[n.y, 1]-sum(clin.FU50))*log(1-d50[n.y]/n.f50[n.y]-sum(p50)/n.f50[n.y])
    return(fin)
}

