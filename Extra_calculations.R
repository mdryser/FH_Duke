## Script for:
## 1. program sensitivity
## 2. overdiagnosis
## 3. goodness of fit

#########################
## 1. Program Sensitivity
#########################

N=1000000


# The estimates from the combined CNBSS
x.vec<-c(psi=0.3588, lambda=0.4554,w= 0.0014,beta=0.2782)
source('trials_FU.r')
t.vec<-c(20,40,42:45)
set.seed(33)
source('synthetic_data_2017-06-26.r')
trial.data
vec40<-trial.data[1,]

# The estimates from the combined CNBSS
x.vec<-c(psi=0.3588, lambda=0.4554,w= 0.0014,beta=0.2782)
source('trials_FU.r')
t.vec<-c(20,50,52:56)
set.seed(33)
source('synthetic_data_2017-06-26.r')
trial.data
vec50<-trial.data[1,]

sens<-c(vec40[1]/sum(vec40[1:2]),vec50[1]/sum(vec50[1:2]), (vec40[1]+vec50[1])/(sum(vec40[1:2])+sum(vec50[1:2])))

# 1: senstivity among 40yo, 2: sensitivity among 50yo; 3: combined sensitivity
print(sens)

#########################
## 2. Overdiagnosis
#########################

N=100000.
x.vec<-c(psi=0.3588, lambda=0.4554,w= 0.0014,beta=0.2782)
#x.vec<-c(psi=0.1, lambda=0.4554,w= 0.0014,beta=1)


t.vec<-c(20,45:50)

source('synthetic_data_2017-06-26.r')

#trial.data
#count.prog.ind
#rowSums(count.prog.ind)

# overdiagnosis by round
count.prog.ind[,2]/rowSums(count.prog.ind)
# total overdiagnosis
sum(count.prog.ind$indolent)/sum(count.prog.ind)


t.vec<-c(20,55:60)

source('synthetic_data_2017-06-26.r')

#trial.data
#count.prog.ind
#rowSums(count.prog.ind)

# overdiagnosis by round
count.prog.ind[,2]/rowSums(count.prog.ind)
# total overdiagnosis
sum(count.prog.ind$indolent)/sum(count.prog.ind)


#########################
## 3. Goodness of Fit
#########################

source('trials_FU.R')

