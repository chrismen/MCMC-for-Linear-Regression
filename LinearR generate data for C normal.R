
# an example for dummy variables from the Michigan state university

# clean everything first 

ls()
rm(list=ls())
gc()



# ----------- an R example from R   ----------------------


m0 <- matrix(NA, 4, 0)
rownames(m0)

m2 <- cbind(1, 1:4)

colnames(m2, do.NULL = FALSE)

colnames(m2) <- c("x","Y")

rownames(m2) <- rownames(m2, do.NULL = FALSE, prefix = "Obs.")
m2


m2 <- cbind(1, 1:4)



dimnames(m2) = list(   c("row1", "row2", "row3", "row4"),         # row names 
                       c("col1", "col2")) # column names 

m2["row2", "col2"] # element at 2nd row, 2nd column 

m2["row4", ] # element at 2nd row

aa="row3"

m2[aa, ] # element at 2nd row, 3rd column


# --- end of the example 



# ---------------generate a mixture of two normal distribution------------
par(mar=c(1,1,1,1))

par(mar=c(5,4,4,5)+  2)

n =4000 # the sample size

nn=n+1000


b = seq(1:nn)
#z = seq(1:nn)


beta1 = 2

beta2 = 3

beta3= 4



mu1 = -0.75

#mu2 = 0.75

sigma1 = 0.2

#sigma2 = 0.6

#rho1 = 0.25

likelihood=0





for (i in 1:nn)
{
  
  {b[i]=rnorm(1,mu1,sigma1)}
  
  
}


par(mar=c(1,1,1,1))

hist(b, breaks=50,col="red", xlab="",ylab="",yaxt="n",xaxt='n'  ,main=""  ) 

title(main="Histogram with density Curve for a mixture of two normal distributions", font = 3, cex.main = 2, col.main="blue")




axis(2, col.axis="red", las=1,cex.axis=2)

par(mar=c(1,1,1,1))
par(new=TRUE)
# Kernel Density Plot
d <- density(b) # returns the density data 


plot(d, xlab="",ylab="", main = "",type="l",col="blue", lwd=2 ,xaxt='n',yaxt='n') # plots the results

axis(4,  col.axis="blue", las=2, cex.axis=2, tck=-.01)

axis(1,  col.axis="blue", las=2, cex.axis=2, tck=-.01)





# ---------------generate a mixture of two Gamma distribution ------------



x1 = seq(1:nn)
x2 = seq(1:nn)

x3 = seq(1:nn)



#x0=seq(1:n)

y = seq(1:nn)



# generate an artifical data set for simulation studies

i=1

x1[i] = (2 * i - 1) / n - 1

x2[i] = x1[i] ^ 2

x3[i] = x1[i]^3


y[i] = b[i] + beta1 * x1[i] + beta2 * x2[i] +  beta3 * x3[i] 

dx <- data.frame( x1[i], x2[i], x3[i])                ## example data
write.table(dx, "simulated-data-normal x.txt", row.names = FALSE, col.names = FALSE)


dy <- data.frame(y[i])                ## example data
write.table(dy, "simulated-data-normal y.txt", row.names = FALSE, col.names = FALSE)





for (i in 2:nn)
  
  
{
  # x0[i]=1
  
  x1[i] = (2 * i - 1) / n - 1
  
  x2[i] = x1[i] ^ 2
  
  x3[i] = x1[i]^3
  
  
  y[i] = beta1 * x1[i] + beta2 * x2[i] +  beta3 * x3[i] + b[i]
  
  dx <- data.frame( x1[i], x2[i], x3[i])                ## example data
  write.table(dx, "simulated-data-normal x.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  
  
  dy <- data.frame(y[i])                ## example data
  write.table(dy, "simulated-data-normal y.txt", row.names = FALSE, col.names = FALSE, append = TRUE)
  
  
  
}

plot(y)



h <- read.table("simulated-data-normal x.txt")


te=dim(h)

rows=te[1]

cols=te[2]

a=t(h[,1])

for (i in 2:cols)
  
{
  
  a=cbind(a, t(h[,i]))
  
  
}


write.table(t(a), file="simulated-data-normal x colum.txt", row.names=FALSE, col.names=FALSE)

