##################################################
# Summarize analyses of CNBSS-I, and CNBSS-II
##################################################


##################################################
# CNBSS-I data:
##################################################
cnbss40 <- matrix(c(98, 19, 25214,
                    39, 16, 22424,
                    44,  8, 22066,
                    52, 10, 21839,
                    26, NA , 14146),
                  byrow=TRUE,
                  ncol=3)


## from Baines 2016 Table 2A; years of FU: 1,2,3,4,5, 10, 15, 20
FU.40=c(9,59,112,168,221) # cumulative numbers of observed cases
a<-FU.40[2:5]-FU.40[1]  # cumulative starting in year two of follow up
a<-round(a*14111/24742) # rescale back to the number of people who remained in screening until 5th scree
FU.40[2:5]<-FU.40[1]+a # cumulative numbers adjusted

cnbss40[5,2]<-FU.40[ind.FU] # complete the matrix

if(ind.FU>1){
FU.40 <- c(FU.40[1],FU.40[2:ind.FU]-FU.40[1:(ind.FU-1)]) # annual rates (instead of cumulative numbers)
} else {
  
  FU.40 <- c(FU.40[1])
  }

##################################################
# CNBSS-II data:
##################################################
cnbss50 <- matrix(c(142, 15, 19711,
                    66, 10, 17669,
                    43,  9, 17347,
                    54,  9, 17193,
                    28, NA, 9876),
                  byrow=TRUE,
                  ncol=3)

## from Baines 2016 Table 2B; years of FU: 1,2,3,4,5, 10, 15, 20
FU.50=c(5, 48, 78, 123, 181) # observed numbers
a<-FU.50[2:5]-FU.50[1]  # cumulative starting in year two of follow up
a<-round(a*9843/19111) # rescale back to the number of people who remained in screening until 5th screen
FU.50[2:5]<-FU.50[1]+a # cumulative numbers adjusted

cnbss50[5,2]<-FU.50[ind.FU] # complete the matrix

if(ind.FU>1){
FU.50 <- c(FU.50[1],FU.50[2:ind.FU]-FU.50[1:(ind.FU-1)]) # annual rates (instead of cumulative numbers)
} else {
  
  FU.50 <- c(FU.50[1])
}
