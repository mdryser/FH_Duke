##################################################
##################################################
# Simulation study 
##################################################
##################################################


##################################################
# Parameters: Synthetic Data
##################################################

MMM<-1000 # number of runs

##################################################
# Parameters: Synthetic Data
##################################################
N<-50000   # total number of women 
psi<- 0.25  # fraction of indolent cancers
lambda<-1/2.5 # mean sojourn time P->C
w<-.0025;    # hazard rate for onset of P/I
beta<-.8  # screening sensitivity
t.vec<-c(20, 50:54,60) # time vector:
# t.vec[1]<- age where w beomes positive
# t.vec[2]<- age at first screen
# t.vec[end-1]<- last screen
# t.vec[end]<- end of follow-up

final<-data.frame( psi=numeric(MMM), psi.se=numeric(MMM),
                   lambda=numeric(MMM), lambda.se=numeric(MMM),
                   w=numeric(MMM), w.se=numeric(MMM), beta=numeric(MMM),
                   beta.se=numeric(MMM), T2E=numeric(MMM), problem=numeric(MMM))


##################################################
# Prepare likelihood functions
##################################################

source('likelihood_2017-06-26.r')


##################################################
# Initial values and limits for optimization algo
##################################################
# time unit: years

x0 <- c(psi=0.5, lambda=1/10, w=1/100, beta=0.5)
lb <- c(psi=0, lambda=0, w=0, beta=0)
ub <- c(psi=1, lambda=5, w=5,beta=1)


##################################################
# Main function to generate the data
##################################################

# z is a vector of size 2;
# z[1] = 0: all 4 parameters
#      = k>0: kth parameter is fixed
# z[2] = value of fixed parameter (if any)

tablerow4 <- function(x0, lb, ub, yname, t.vec, z, clin.FU, likelihood=EstPsiBMucstab4Full.fun){
  y <- get(yname)
  # estimate parameters
  est <- nlminb(start=x0,
                obj=likelihood,
                lower=lb,
                upper=ub,
                y=y,
                t=t.vec,
                z=z,
                clin.FU=clin.FU)
  
  # calculate standard errors
  hess <- numDeriv::hessian(func=likelihood , x=est$par, y=y, t=t.vec, z=z, clin.FU=clin.FU)
  se <- sqrt(diag(fBasics::inv(hess)))
  
  # re-estimate parameters disallowing indolent cancers
  estm1 <- nlminb(start=x0[-1],
                  obj=likelihood,
                  lower=lb[-1],
                  upper=ub[-1],
                  y=y,
                  t=t.vec,
                  z=c(1,0),
                  clin.FU=clin.FU)
  # test for null model (no indolent model)
  pval <- pchisq(2*(estm1$obj-est$obj), 1, lower.tail=FALSE)/2
  
  
  # summarize results
  row <- data.frame(psi=est$par[1], psi.se=se[1],
                    lambda=est$par[2], lambda.se=se[2],
                    w=est$par[3], w.se=se[3],
                    beta=est$par[4], beta.se=se[4],
                    T2E=pval,problem=est$convergence, 
                    row.names=NULL,
                    check.names=FALSE)
  return(row)
}




for(mmm in 1:MMM){
  
  print(mmm)
  source('synthetic_data_2017-06-26.r')
  
  final[mmm,]<-tablerow4(x0=x0, lb=lb, ub=ub, yname='trial.data', t.vec=t.vec, z=c(0,0), clin.FU=last.screen, likelihood = EstPsiBMucstab4Full.fun)
  
}



synth<-data.frame(psi=psi, psi.se=NA,
                  lambda=lambda, lambda.se=NA,
                  w=w, w.se=NA,
                  beta=beta, beta.se=NA,
                  T2E=0, 
                  row.names=NULL,
                  check.names=FALSE)


if(sum(final$problem)>0){
  
  print("watch out")
}

final$T2E<-1-1*(final$T2E<.05)
fff<-rbind(synth,c(mean(final$psi),sd(final$psi) , mean(final$lambda),sd(final$lambda) , mean(final$w),sd(final$w), mean(final$beta), sd(final$beta), mean(final$T2E)))
fff<-round(fff,5)

save(fff,file="Bias_SE_stats_7_19_2017.RData")