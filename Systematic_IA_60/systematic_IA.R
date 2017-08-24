##############################################################
##############################################################
## Systematic Simulation Study - Mixture model
## PART 1: IDENTIFIABILITY AND LIKELIHOOD RATIO TEST
##############################################################
##############################################################

library(ggplot2)
# Multiple plot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

set.seed(44)
##################################################
# Parameters: Operational
##################################################
scale.factor= 3 # this is to scale rate window
per.factor=.2  # this is to scale the probability window
n.n<-10 # number of grid points in the psi-lambda plane
ZZ <-100 # number of MCMC runs

##################################################
# Parameters: Synthetic Data
##################################################
N<-50000   # total number of women 
psi.bound<- seq(from=0, to=.75, length.out = n.n)  # fraction of indolent cancers
lambda.bound<-10^seq(from=-1, to=log10(2), length.out = n.n)  # mean sojourn time P->C
w<-.0025;    # hazard rate for onset of P/I
w.bound=c(w/scale.factor, w*scale.factor)
beta<-.8  # screening sensitivity
beta.bound=c(max(0, beta-per.factor), min(1, beta+per.factor))
t.vec<-c(0, 50:54, 60) # time vector:
# t.vec[1]<- age where w beomes positive
# t.vec[2]<- age at first screen
# t.vec[end-1]<- last screen
# t.vec[end]<- end of follow-up

##################################################
# Prepare likelihood functions
##################################################
source('likelihood_2017-06-26.r')

##################################################
# Initial values and limits for optimization algo
##################################################
# time unit: years

x0 <- c(psi=0.5, lambda=1, w=1/100, beta=0.5)
lb <- c(psi=0, lambda=0, w=0, beta=0)
ub <- c(psi=1, lambda=20, w=5,beta=1)


##################################################
# Give initial values and limits for optimization
##################################################

# dataframe that contains the psi/lambda values
result_mat <- expand.grid(psi = psi.bound, lambda=lambda.bound)

# add the outputs
# IdM=Identifiable (1=yes; 0=no)
# LRM=likelihood ratio test (1=reject H0(psi>0); 0=fail to reject H0)
# FailM= # of convergence failures during sampling
result_mat$IdM=rep(0, length(psi.bound))
result_mat$LRM=rep(0, length(psi.bound))
result_mat$FailM=rep(0, length(psi.bound))


  
  # outer loop: run over parameter space
  for(kk in 1:dim(result_mat)[1]){
    
    
    # keep track of the number of failed convergence status
    fail=0
    
    psi<-result_mat$psi[kk] # fraction of indolent cancers
    lambda=result_mat$lambda[kk] # mean sojourn time P->C
    
    
    # inner loop: run ZZ samples for each parameter set
    for(zum in 1:ZZ){
      
      print(c(zum,kk,fail))
      
    x.vec <- c(psi=psi, lambda=lambda, w=w, beta=beta)
      
    source('synthetic_data_2017-06-26.r') # generate the data
    clin.FU=last.screen
    
    # compute MLE
    est<-nlminb(start=x0,
                obj=EstPsiBMucstab4Full.fun,
                lower=lb,
                upper=ub,
                y=trial.data,
                t=t.vec,
                z=c(0,0),
                clin.FU=clin.FU)
    
    # extract MLE
    peak<-est$objective
    params<-est$par
    fail=fail+est$convergence # update the failure status
    
    
    # prepare the profile likelihood boundaries
    pl.int=matrix(0,nrow=4, ncol=2)
    pl.int[1,]=c(max(0, psi-per.factor), min(1, psi+per.factor))
    pl.int[2,]=c(lambda/scale.factor, lambda*scale.factor)
    pl.int[3,]=w.bound
    pl.int[4,]=beta.bound
    
    # Compute for each paramter
    NLL.track=matrix(0, nrow=4,ncol=2)
    
    for(qqq in 1:4){
      
      # Fix qqq-th entry
      # order: 1) psi, 2) lambda, 3) w, 4) beta
      
      for(j in 1:2){
        
        ak <- nlminb(start=x0[-qqq],
                     obj=EstPsiBMucstab4Full.fun,
                     lower=lb[-qqq],
                     upper=ub[-qqq],
                     y=trial.data,
                     t=t.vec,
                     z=c(qqq,pl.int[qqq,j]),
                     clin.FU=clin.FU
        )
        
        NLL.track[qqq,j]<-ak$objective-peak         # MLE for fixed qqq-the parameter at position j
        fail=fail+est$convergence   # update failure count
        
        
      }
    }
    
    
    # compute the value at the boundary for LRT
    ak <- nlminb(start=x0[-1],
                 obj=EstPsiBMucstab4Full.fun,
                 lower=lb[-1],
                 upper=ub[-1],
                 y=trial.data,
                 t=t.vec,
                 z=c(1,0),
                 clin.FU=clin.FU
    )
    
    LRtestNLL=ak$objective-peak
    
    # ChiSquare test for LRT and CI
    CI.line<- qchisq(0.95, 1)/2
    
    #### ascertain if identifiable
    if(result_mat$psi[kk]<per.factor){  # if 0 is within per.factor of the psi value
      # in this case make sure that the right the right boundary point is above CI.line for psi, and both sides for other params
      result_mat$IdM[kk]= result_mat$IdM[kk]+ifelse(NLL.track[1,2]>CI.line & all(NLL.track[2:4,]>CI.line) , 1, 0  )
      
    } else{ # if 0 is not within per.factor of psi value
      # in this case all boundary points have to lie above the cutoff
      result_mat$IdM[kk]= result_mat$IdM[kk]+ifelse(all(NLL.track>CI.line) , 1, 0  )
      
    }
    #### if the null is rejected (only if 0 NLL exceed CI.line)
    result_mat$LRM[kk]=result_mat$LRM[kk]+ (LRtestNLL>CI.line)
    ##### keep track of the number of fails for this parameter set
    result_mat$FailM[kk]= result_mat$FailM[kk]+fail
  }
  
}


save(result_mat, t.vec, w, beta, ZZ, file="Exp_60_identifiability.RData")


#load("Exp_60_identifiability.RData")
###############################
# Plot the identifiability 
# and the LRtest maps
################################

# p.IdM<-ggplot(data = result_mat, aes(x = psi, y = log10(lambda))) +  geom_raster(aes(fill = IdM/ZZ), interpolate=FALSE)+  ggtitle("Fraction of identifiable runs") +theme(aspect.ratio=1)+
#   scale_fill_continuous(limits = c(0,1), guide = guide_colorbar(title="Fraction",nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# p.LRM<-ggplot(data = result_mat, aes(x = psi, y = log10(lambda))) + geom_raster(aes(fill = LRM/ZZ), interpolate=FALSE)+  ggtitle("Fraction with positve LR test") +theme(aspect.ratio=1)+
#   scale_fill_continuous(limits = c(0,1), guide = guide_colorbar(title="Fraction", nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# multiplot(p.IdM, p.LRM, cols=2)





