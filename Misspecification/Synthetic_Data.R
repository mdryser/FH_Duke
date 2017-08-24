#########################################
#########################################
### Generate a synthetic data set
#########################################
#########################################

len<-length(t.vec)
#########################################
## STEP 1
## Age at onset of P/I disease
#########################################
time.P<-t.vec[1]+rexp(N,rate=w)

#########################################
## STEP 2
## Progression to invasive
#########################################
prog.ind<-(runif(N)>psi)*1

#########################################
## STEP 3
## Time to onset of clinical disease
#########################################

time.C<-rep(Inf,N) # default is Inf
a<-which(prog.ind>0) # label those with prog dis
time.C[a]<-time.P[a]+rexp(length(a),rate=lambda)


############################################
## STEP 4
## Screening
############################################
# only include those without clinical before 1st 
# screen

a<-which(time.C>t.vec[2])
time.C<-time.C[a]
time.P<-time.P[a]

# vector to keep track of those not yet diagnosed
# 1 progressive patient
# 2 indolent patient
at.risk.ind<-rep(1,length(a))
at.risk.ind[(time.C==Inf)]<-2

# allocate final output
trial.data<-matrix(0, ncol=3, nrow = (length(t.vec)-2)) 
count.prog.ind<-data.frame(prog=numeric(length(t.vec)-2),indolent=numeric(length(t.vec)-2) )

# Loop over the first k-1 screens
for(k in 2:(len-2)){
  
  # number still at risk
  at.risk<-sum(1*(at.risk.ind>0))
  
  # still at risk, P before screen, C after screen: can be picked up by screen
  dum<-which((time.P<t.vec[k]) & (time.C>t.vec[k]) & (at.risk.ind))
  
  # choose those that are picked up by screen
  dum<-sample(dum, rbinom(1,length(dum), beta), replace = FALSE)   # contains indices of those that are picked up
  
  count.prog.ind[k-1,]<-c(sum(at.risk.ind[dum]==1), sum(at.risk.ind[dum]==2))
  
  at.risk.ind[dum]<-0 # remove from at risk set
  
  # identify those still at risk, become clinical in kth interval: will become clinical
  dam<-which((time.C>t.vec[k]) & (time.C<t.vec[k+1]) & (at.risk.ind))

  at.risk.ind[dam]<-0; # remove those from risk set
  
  # update the output matrix
  trial.data[k-1,]<-c(length(dum), length(dam), at.risk)
  
  
}



### Last (kth) screen, and clinical follow up

# number still at risk
at.risk<-sum(1*(at.risk.ind>0))

# still at risk, P before screen, C after screen: can be picked up by screen
dum<-which((time.P<t.vec[len-1]) & (time.C>t.vec[len-1]) & (at.risk.ind))

# choose those that are picked up by screen
dum<-sample(dum, rbinom(1,length(dum), beta), replace = FALSE)   # contains indices of those that are picked up

count.prog.ind[k-1,]<-c(sum(at.risk.ind[dum]==1), sum(at.risk.ind[dum]==2))

at.risk.ind[dum]<-0 # remove from at risk set

# number of years for clinical follow up
l.FU<-t.vec[len]-t.vec[len-1]

# vector to store number of clinical cases each year
last.screen<-rep(0, l.FU)

for(q in 1:l.FU){

# identify those still at risk, become clinical in qth interval: will become clinical
dam<-which((time.C>(t.vec[len-1]+q-1)) & (time.C<(t.vec[len-1]+q)) & (at.risk.ind))

at.risk.ind[dam]<-0; # remove those from risk set

last.screen[q]<-length(dam)

}

# update the output matrix
trial.data[len-2,]<-c(length(dum), sum(last.screen), at.risk)




