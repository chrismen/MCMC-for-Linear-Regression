
## generate artifical data for a linear regresssion model
## data will be called by MCMC codes written in R and C, respectively.

rm(list=ls())
gc()
n=10000
set.seed(666)
sigma2=0.2
x1 = rnorm(n)         
x2 = rnorm(n)
x3 = rnorm(n)
x=rnorm(n)
y = 1 + 1.5*x1 + 2.2*x2 +3*x3 + sigma2*x
    
dy <- data.frame(y)               
write.table(dy, "R-simulated-logist-data-y.csv", row.names=FALSE, na="",col.names=FALSE, sep=",")
df = data.frame(y=y,x1=x1,x2=x2,x3=x3)
dx=data.frame(x1,x2,x3)
write.table(dx, "C-simulated-logist-data-x.txt",  row.names = FALSE, col.names = FALSE)
write.table(dx, "R-simulated-logist-data-x.csv", row.names=FALSE, na="",col.names=FALSE, sep=",")
write.table(dy, "C-simulated-logist-data-y.txt",row.names = FALSE, col.names = FALSE)





