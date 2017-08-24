##################################################
##################################################
# Simulation study 
##################################################
##################################################


##################################################
# Parameters: Synthetic Data
##################################################

MMM<-200# number of runs

##################################################
# Parameters: Synthetic Data
##################################################
N<-50000   # total number of women 
psi.bound<- seq(from=0, to= 0.8, by=0.2)  # fraction of indolent cancers
lambda<-1/2.5 # mean sojourn time P->C
w<-.0025;    # hazard rate for onset of P/I
beta<-.4  # screening sensitivity
t.vec<-c(20, 50:54,60) # time vector:
# t.vec[1]<- age where w beomes positive
# t.vec[2]<- age at first screen
# t.vec[end-1]<- last screen
# t.vec[end]<- end of follow-up


result_mat <- expand.grid( N=c(1:MMM), psi = psi.bound)

result_mat$lambda<-0
result_mat$w<-0
result_mat$beta <-0
result_mat$problem<-0


##################################################
# Prepare likelihood functions
##################################################

source('likelihood.r')


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
  est <- nlminb(start=x0[-1],
                obj=likelihood,
                lower=lb[-1],
                upper=ub[-1],
                y=y,
                t=t.vec,
                z=z,
                clin.FU=clin.FU)
  
  # calculate standard errors
  hess <- numDeriv::hessian(func=likelihood , x=est$par, y=y, t=t.vec, z=z,clin.FU=clin.FU )
  se <- sqrt(diag(fBasics::inv(hess)))
  

  
  # summarize results
  row <- data.frame(lambda=est$par[1], lambda.se=se[1],
                    w=est$par[2], w.se=se[2],
                    beta=est$par[3], beta.se=se[3], problem=est$convergence,
                    row.names=NULL,
                    check.names=FALSE)
  return(row)
}


for(psi in psi.bound){
  
  
  IM<-matrix(0, ncol = 4, nrow = MMM)

  for(mmm in 1:MMM){
    print(c(mmm,psi))
    
    source('Synthetic_Data.R')
    clin.FU=last.screen
    
    a<-tablerow4(x0, lb, ub, 'trial.data', t.vec, c(1,0), clin.FU)
    IM[mmm,1]<-a$lambda
    IM[mmm,2]<-a$w
    IM[mmm,3]<-a$beta
    IM[mmm,4]<-a$problem
    
    
  }
  
  result_mat[which(result_mat$psi==psi),3:6]<-IM
  
}




result_mat<-as.data.frame(result_mat)

result_mat$psi<-as.factor(result_mat$psi)


save(result_mat,w, beta, lambda, N, mmm, t.vec, file="Misspecification_beta_p4.RData")

 p.1 <- ggplot(result_mat, aes(x=psi, y=lambda)) +   geom_violin() + geom_hline(yintercept = lambda, linetype="dotted") +geom_boxplot(width=0.1)
 p.2 <- ggplot(result_mat, aes(x=psi, y=w)) +   geom_violin() + geom_hline(yintercept = w, linetype="dotted") +geom_boxplot(width=0.1)
 p.3 <- ggplot(result_mat, aes(x=psi, y=beta)) +   geom_violin() + geom_hline(yintercept = beta, linetype="dotted") +geom_boxplot(width=0.1)

 grid.arrange(p.1,p.2,p.3, top = "Model Misspecification",  layout_matrix = matrix(c(1,2,3), ncol=1, byrow=TRUE))





