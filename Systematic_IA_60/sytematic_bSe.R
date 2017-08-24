

##############################################################
##############################################################
## Systematic Simulation Study - Mixture model
## PART 2: BIAS AND STANDARD ERROR OF THE ESTIMATOR
##############################################################
##############################################################


##################################################
# Parameters: Synthetic Data
##################################################

#rm(list=ls())

n.n=10 # number of grid points in psi and lambda directions
MMM<-100 # number of samples per parameter combination

N<-50000   # total number of women 
psi.bound<- seq(from=0, to=.75, length.out = n.n)  # fraction of indolent cancers
lambda.bound<-10^seq(from=-1, to=log10(2), length.out = n.n)  # mean sojourn time P->C
w<-.0025;    # hazard rate for onset of P/I
beta<-.8  # screening sensitivity
t.vec<-c(0, 50:54, 60) # time vector:
# t.vec[1]<- age where w beomes positive
# t.vec[2]<- age at first screen
# t.vec[end-1]<- last screen
# t.vec[end]<- end of follow-up


##################################################
# Prepare the output
##################################################

result_mat <- expand.grid(psi = psi.bound, lambda=lambda.bound)

# collect the bias
result_mat$bias_psi=rep(0, length(psi.bound))
result_mat$bias_lambda=rep(0, length(psi.bound))
result_mat$bias_w=rep(0, length(psi.bound))
result_mat$bias_beta=rep(0, length(psi.bound))

# collect the standard error
result_mat$SE_psi=rep(0, length(psi.bound))
result_mat$SE_lambda=rep(0, length(psi.bound))
result_mat$SE_w=rep(0, length(psi.bound))
result_mat$SE_beta=rep(0, length(psi.bound))

# final<-data.frame( psi=numeric(MMM), psi.se=numeric(MMM),
#                    lambda=numeric(MMM), lambda.se=numeric(MMM),
#                    w=numeric(MMM), w.se=numeric(MMM), beta=numeric(MMM),
#                    beta.se=numeric(MMM), p=numeric(MMM), problem=numeric(MMM))


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
                    p=pval, 
                    row.names=NULL,
                    check.names=FALSE)
  return(row)
}


##################################################
# Loop over parameter space
##################################################

for(jj in 1:dim(result_mat)[1]){
  
  inter_mat=matrix(0, nrow = MMM, ncol=4) # zwischenschritt
  psi<-result_mat$psi[jj] # fraction of indolent cancers
  lambda=result_mat$lambda[jj] # mean sojourn time P->C
  
  # number of runs per parameter set
  for(mmm in 1:MMM){
    
    print(c(jj,mmm))
    x.vec <- c(psi=psi, lambda=lambda, w=w, beta=beta)
    
    source('synthetic_data_2017-06-26.r') # generate the data
    clin.FU=last.screen
    out<-tablerow4(x0, lb, ub, 'trial.data', t.vec, c(0,0), clin.FU) # compute output
    inter_mat[mmm, ]=c(out$psi, out$lambda, out$w, out$beta) # store intermediary data
    
  }
  
  # update the output
  result_mat$bias_psi[jj]<-psi-mean(inter_mat[,1])
  result_mat$bias_lambda[jj]<-lambda-mean(inter_mat[,2])
  result_mat$bias_w[jj]<-w-mean(inter_mat[,3])
  result_mat$bias_beta[jj]<-beta-mean(inter_mat[,4])
  
  result_mat$SE_psi[jj]<-sd(inter_mat[,1])
  result_mat$SE_lambda[jj]<-sd(inter_mat[,2])
  result_mat$SE_w[jj]<-sd(inter_mat[,3])
  result_mat$SE_beta[jj]<-sd(inter_mat[,4])
  
  
  
}


save(result_mat, t.vec, w, beta, MMM, file="Exp_60_bias_se.RData")


#### PLOTS: BIAS #########

# m=c(0,max(abs(result_mat$bias_psi)))
# p.bias.psi<-ggplot(data = result_mat, aes(x = psi, y = log10(lambda))) +  geom_raster(aes(fill = abs(bias_psi)), interpolate=FALSE)+  ggtitle("Bias for psi") +theme(aspect.ratio=1)+
#   scale_fill_continuous(limits = m, guide = guide_colorbar(title="|Bias|",nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# m=c(0,max(abs(result_mat$bias_lambda)))
# p.bias.lambda<-ggplot(data = result_mat, aes(x = psi, y = log10(lambda))) +  geom_raster(aes(fill = abs(bias_lambda)), interpolate=FALSE)+  ggtitle("Bias for lambda") +theme(aspect.ratio=1)+
#   scale_fill_continuous(limits = m, guide = guide_colorbar(title="|Bias|",nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# m=c(0,max(abs(result_mat$bias_w)))
# p.bias.w<-ggplot(data = result_mat, aes(x = psi, y = log10(lambda))) +  geom_raster(aes(fill = abs(bias_w)), interpolate=FALSE)+  ggtitle("Bias for w") +theme(aspect.ratio=1)+
#   scale_fill_continuous(limits = m, guide = guide_colorbar(title="|Bias|",nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# m=c(0,max(abs(result_mat$bias_beta)))
# p.bias.beta<-ggplot(data = result_mat, aes(x = psi, y = log10(lambda))) +  geom_raster(aes(fill = abs(bias_beta)), interpolate=FALSE)+  ggtitle("Bias for beta") +theme(aspect.ratio=1)+
#   scale_fill_continuous(limits = m, guide = guide_colorbar(title="|Bias|",nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# multiplot(p.bias.psi,  p.bias.lambda, p.bias.w,  p.bias.beta, cols=2)
# 
# 
# #### PLOTS: SE #########
# 
# m=c(0,max(abs(result_mat$SE_psi)))
# p.se.psi<-ggplot(data = result_mat, aes(x = psi, y = log10(lambda))) +  geom_raster(aes(fill = abs(SE_psi)), interpolate=FALSE)+  ggtitle("Stadard error for psi") +theme(aspect.ratio=1)+
#   scale_fill_continuous(limits = m, guide = guide_colorbar(title="SE",nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# m=c(0,max(abs(result_mat$SE_lambda)))
# p.se.lambda<-ggplot(data = result_mat, aes(x = psi, y = log10(lambda))) +  geom_raster(aes(fill = abs(SE_lambda)), interpolate=FALSE)+  ggtitle("Standard error for lambda") +theme(aspect.ratio=1)+
#   scale_fill_continuous(limits = m, guide = guide_colorbar(title="SE",nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# m=c(0,max(abs(result_mat$SE_w)))
# p.se.w<-ggplot(data = result_mat, aes(x = psi, y = log10(lambda))) +  geom_raster(aes(fill = abs(SE_w)), interpolate=FALSE)+  ggtitle("Standard error for w") +theme(aspect.ratio=1)+
#   scale_fill_continuous(limits = m, guide = guide_colorbar(title="SE",nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# m=c(0,max(abs(result_mat$SE_beta)))
# p.se.beta<-ggplot(data = result_mat, aes(x = psi, y = log10(lambda))) +  geom_raster(aes(fill = abs(SE_beta)), interpolate=FALSE)+  ggtitle("Standard error for beta") +theme(aspect.ratio=1)+
#   scale_fill_continuous(limits = m, guide = guide_colorbar(title="SE",nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# multiplot(p.se.psi,  p.se.lambda, p.se.w,  p.se.beta, cols=2)

