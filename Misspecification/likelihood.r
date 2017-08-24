####################################################
####################################################
# Survival function for U (mixture progression time)
#
# t is the time
# psi: fraction of indolent cancers
# lambda: rate P->C
####################################################
####################################################

Q.surv <-function(t,psi,lambda){
  return(psi+(1-psi)*exp(-lambda*t))
}


####################################################
####################################################
# Screen-detected cancers: Likelihood contribution
#
# theta: model parameters
# t: vector of times
# j: is the screen number [0,1,2,3,...]
#
####################################################
####################################################

D.j <-function(theta, t,j){
  
  psi<-theta[1]
  lambda<-theta[2]
  w<- theta[3]
  beta<-theta[4]
  
  out<-0
  j<-j+2
  
  for(k in 2:j){
    
    out<-out+beta*(1-beta)^(j-k)*(psi*(exp(-w*t[k-1])-exp(-w*t[k]))+w*(1-psi)/(lambda-w)*(exp(t[k]*(lambda-w)-lambda*t[j])-exp(t[k-1]*(lambda-w)-lambda*t[j])))
    
  }
  
  # correction 6/14
  out=out*exp(w*t[1])
  
  return(out)
  
}


####################################################
####################################################
# Interval-detected cancers: Likelihood contribution
# Prior to last screen
#
# theta: model parameters
# t: vector of times
# j: is the screen number [0,1,2,3,...]
#
####################################################
####################################################

I.j <-function(theta, t, j){
  
  psi<-theta[1]
  lambda<-theta[2]
  w<-theta[3]
  beta<-theta[4]
  
  j<-j+2
  out<-0
  
  for(k in 2:j){
    
    a<-exp(-lambda*t[j]+t[k]*(lambda-w))-exp(-lambda*t[j+1]+t[k]*(lambda-w))
    b<-exp(-lambda*t[j]+t[k-1]*(lambda-w))-exp(-lambda*t[j+1]+t[k-1]*(lambda-w))
    
    out<-out+ (1-beta)^(j-k+1)*(a-b)
    
  }
  
  out<- (1-psi)*w/(lambda-w)*out 
  
  out<-out+(1-psi)*(exp(-w*t[j])-exp(-w*t[j+1])-w/(lambda-w)*(exp(-t[j+1]*w)-exp(-lambda*t[j+1]+(lambda-w)*t[j])))
  
  
  # correction 6/14
  out=out*exp(w*t[1])
  
  return(out)
  
  
}

####################################################
####################################################
# Interval-detected cancers: Likelihood contribution
# After last screen
#
# theta: model parameters
# t: vector of times 
# l: is the year of follow up (l=1,...,L)
#
####################################################
####################################################

Ifin.l<-function(theta, t, l){
  
  psi<-theta[1]
  lambda<-theta[2]
  w<-theta[3]
  beta<-theta[4]
  
  J<-length(t)-1 # index of last screen
  tJ<-t[J] # time of last screen
  L<-t[J+1]-t[J] # number of follow-up years after last screen
  
  out<-0
  
  
  out<-exp(-lambda-w*(tJ+l-1))-exp(-l*lambda-w*tJ)  
  
  
  for(k in 2:J){
    
    out<-out + (1-beta)^(J-k+1)*(exp(-w*t[k]-lambda*(tJ+l-t[k]))-exp(-w*t[k-1]-lambda*(tJ+l-t[k-1])))

    
  }
  
  out<-out*w/(lambda-w)*(exp(lambda)-1)
  
  out<-out+exp(-w*(tJ+l-1))-exp(-w*(tJ+l))
  
  out<- out - w/(lambda-w)*(exp(-w*(tJ+l))-exp(-lambda-w*(tJ+l-1)))
  
  out<-(1-psi)*out
  
  
  # correction 6/14
  out=out*exp(w*t[1])
  
  return(out)
  
  
}



####################################################
####################################################
# Conditioning on T_P+T_C>t_0
#
# theta: model parameters
# t0: time of first screen
# j: is the screen number [0,1,2,3,...]
#
####################################################
####################################################

CF.fun <-function(theta, t0){
  
  psi<-theta[1]
  lambda<-theta[2]
  w<-theta[3]

  out<-psi+(1-psi)*(exp(-w*t0)+w/(lambda-w)*(exp(-t0*w)-exp(-t0*lambda)))
  
  return(out)
  
}

##################################################
##################################################
# Full Likelihood: 4 parameters
# Use full likelihood to estimate 4 parameters:
# psi:  proportion of non-progressive cancers
# lambda:  rate of progression (P->C)
# w: incidence of preclinical disease (N->P/I)
# beta: screening sensitivity
#
# 
# Hardcoded parameter inputs:
# Time of onset of susceptibility to P/I 
# Screening times
##################################################
##################################################

EstPsiBMucstab4Full.fun <- function(x,y,t, z, clin.FU) {
    # x: vector of parameters (psi, lambda, w, beta)
    # y: rows are: y[1:3, ]=matrix(s[1], r[1], n[1],
    #                                   s[2], r[2], n[2],
    #                                   ..., ..., ...)
    #    where s is the number of screen detections and
    #          r is the number of interval cases
    #          n is the number of participants
    # z: if z[1]>0: it is the number of parameter to be fixed
    #       z[2]: value of the parameter
    

    if (sum(z[[1]])>0){
      theta<-rep(-1,4)
      theta[z[[1]]]<-z[[2]]
      theta[theta==-1]<-x
      
    } else {
      
      theta<-x
      
      }
    
    n.y<-nrow(y)
  
    ##################################################
    # Calculate probability of screen detection at
    # each exam
    ##################################################
  
    v<-t(0:(n.y-1))
    
    d<-apply(v, 2, function(x) D.j(theta, t, x))
    

    ##################################################
    # Calculate probability of clinical diagnosis in
    # each interval
    ##################################################
    
    q<-apply(v, 2, function(x) I.j(theta, t, x))
  
    ##################################################
    # Calculate and return full log likelihood
    ##################################################
    
    # old version with only one normalization term (approximation)
    #n.f<-rep(CF.fun(theta, t[2]-t[1]),n.y)
    #n.f<-rep(1,n.y)
    
    # NEW version with exact normalization terms
    
    n.f<-rep(CF.fun(theta, t[2]-t[1]),n.y)
    QQ <- cumsum(d+ q )
    n.f[2:n.y]<-n.f[2:n.y]-QQ[1:(n.y-1)]
     

    fin<-0
    
    # step 1: all but the last round
    
    for(jj in 1:(n.y-1))  {

      fin<- fin -(y[jj, 1]*log(d[jj]/n.f[jj])+y[jj, 2]*log(q[jj]/n.f[jj])+(y[jj, 3]-y[jj,1]-y[jj,2])*log(1-d[jj]/n.f[jj]-q[jj]/n.f[jj]))

    }
        
    # step 2: the last round
    # The screen-detection contribution
    fin <- fin - y[n.y,1]*log(d[n.y]/n.f[n.y])
    
    # The annual clinical incidence contributions
    p<-apply(t(1:length(clin.FU)), 2, function(x) Ifin.l(theta, t, x))
    
    fin <- fin - clin.FU %*% log(p/n.f[n.y])
    
    fin <- fin - (y[n.y,3]-y[n.y,1]-sum(clin.FU))*log(1-d[n.y]/n.f[n.y]-sum(p)/n.f[n.y])    

    return(fin)
}



