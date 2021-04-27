#--------------------------------------------------------------------------
# # Design operator calX (as a n-list of d1Xd2-matrices) and noise
# Uniform distribution. Indexing sample count uses "columnwise stacking" 
#--------------------------------------------------------------------------
data_genXe <- function(n,d1,d2){
  p = d1*d2
  elm = matrix(rep(0,p), nrow = d1)
  calX = rep(list(elm), n)
  noise = rep(0,n)
  aux = rdunif(n, a = 1, b = p)
  for (i in 1:n) {
    column = ceiling(aux[i]/d1)
    row = aux[i]%%d1
    if(row == 0){
      row = d1
    }
    calX[[i]][row,column] = 1
    noise[i] = rnorm(1, mean = 0, sd = 1)
  }
  list = list(calX = calX, noise = noise)
  return(list)
}

#-----------------------------------------------
# Function to produce true low-rank parameter
#-----------------------------------------------
data_genP <- function(n,d1,d2,r,B.mag){
  library(matlib)
  library(purrr)
  library(mclust)
  
  U = randomOrthogonalMatrix(d1, r)
  V = randomOrthogonalMatrix(d2, r)
  D = diag(B.mag,r)
  B.true = U%*%D%*%t(V)
  
  return(B.true)
}

#---------------------------------------------------------------------
# Function to produce corrupted response y given (calX, noise, parameter)
#---------------------------------------------------------------------
data_genY <- function(n,step,scale,calX,noise,B.true,theta.mag){
  
  y = rep(0,n)
  for (i in 1:n) {
    y[i] = sum(diag( t(calX[[i]]) %*% B.true )) + noise[i]
  } 
  
  #Produce corrupted response for grid of corruptions
  aux = (n/step)/scale
  Y = matrix(rep(0,n*aux), nrow = n)
  for (j in 1:aux) {
    Y[,j] = y
  }
  for (j in 1:aux) {
    Y[1:(step*j),j] =  Y[1:(step*j),j] + theta.mag
  }
  #Adjoint clean response sample in first column
  Y = cbind(as.matrix(y),Y)
  
  return(Y)
}

#----------------------------------------------------------------------------------
# Square of Euclidean norm function
#----------------------------------------------------------------------------------
norm.sq <- function(x){sum(x^2)}

#----------------------------------------------------------------------------------
# Soft-thresholding function
#----------------------------------------------------------------------------------
soft <- function(x,lambda){ sign(x)*max(abs(x)-lambda,0) }

#----------------------------------------------------------------------------------
# Soft-thresholding function contrained on the l-infty ball of radius a
#----------------------------------------------------------------------------------
soft.linfty <- function(x,lambda,a){ sign(x)*min( max(abs(x)-lambda,0), a) }

#----------------------------------------------------------------------------------
# Function to do low-rank trace regression solver given ( calX, y, B.max ). 
#It will solve for the normalization B/sqrt(p) and design sqrt(p)*calX
#----------------------------------------------------------------------------------
trace.reg <- function(n,d1,d2,B.max,SlopeT,A,calX,y,L,t.final){
  
  #call for prox_sorted_L1()
  library(glmSLOPE)
  
  #Parameter regularization
  d = d1+d2
  p = d1*d2
  m = min(d1,d2)
  aux = 2*log(d)
  lambda = max(B.max,1)*sqrt(aux)*sqrt(p/(m*n))
  lambda.aux = max(B.max,1)*aux*sqrt(log(m))*(sqrt(p)/n)
  lambda = 0.5*max(lambda,lambda.aux)
  
  #Corruption regularization
  if (SlopeT == T){
    tau = 1*max(B.max,1)*sqrt(1/n)
    omega = rep(0,n)
    for (i in 1:n) {
      omega[i] = log(A*n/i)
    }
  } 
  else{
    tau = 1*max(B.max,1)*sqrt(log(n)/n)
    omega = rep(1,n)
  }
  
  #Iinitialization
  B = matrix(0, nrow = d1, ncol = d2)
  calX.B = rep(0, n)
  theta = as.matrix(rep(0,n))
  t = 1
  value = rep(0,t.final)
  
  #Proximal iterations
  for (t in 1:t.final) {
    for (i in 1:n) {
      calX.B[i] = sum(diag(t( sqrt(p)*calX[[i]] )%*%B)) #Scaled design is sqrt(p)*calX
    }
    aux = calX.B + sqrt(n)*theta - sqrt(p)*y #Scaled response is sqrt(p)*y
    value[t] = (1/(2*n))*norm.sq(aux)
    
    #Parameter proximal iteration
    gradB = matrix(rep(0,d1*d2), nrow = d1)
    for (i in 1:n) {
      gradB = gradB + (1/n)*aux[i]*sqrt(p)*calX[[i]] #Scaled design is sqrt(p)*calX
      i = i + 1
    }
    W = B - (1/L)*gradB
    B.aux = svd(W, nu = min(d1, d2), nv = min(d1,d2))
    D = B.aux$d
    for (j in 1:min(d1, d2)) {
      D[j] = soft.linfty(D[j], lambda*sqrt(p)/L, B.max)
      #D[j] = soft(D[j], lambda/L)
    }
    B = (B.aux$u)%*%diag(D)%*%t(B.aux$v)
    
    #Corruption proximal iteration
    gradtheta = (1/sqrt(n))*aux
    w = theta - (1/L)*gradtheta
    aux1 = (tau*sqrt(p)/L)*omega
    theta =  as.matrix(prox_sorted_L1(w,aux1))
    
    t = t + 1 
  }
  
  ls.trace.reg = list(B = B, value = value)
  return(ls.trace.reg)
}

#------------------------------------------------------------------------------------------------------------------
# Function to call data, do trace regression for fixed low-rank parameter with rank r and grid of corruption 
# then compute MSE averaged over N.rep repetitions. Given repetition, the same data is used across all corruptions
#------------------------------------------------------------------------------------------------------------------
trace.reg.MSE <- function(n,d1,d2,r,SlopeT,A,step,scale,B.mag,theta.mag,L,t.final,N.rep){
  
  # Generate low-rank parameter
  B.true = data_genP(n,d1,d2,r,B.mag)
  B.max = max(abs(B.true))
  
  # Initialization of repetition count
  aux = (n/step)/scale
  MSE = rep(0,aux + 1)
  Mvalue = matrix( rep(0, t.final*(aux + 1) ) , nrow = t.final)
  
  for (index in 1:N.rep) {
    # calX and noise generation
    ls0 <- data_genXe(n,d1,d2)
    calX = ls0$calX
    noise = ls0$noise
    
    for (j in 1:(aux+1)) {
      # Parameter and response generation
      Y = data_genY(n,step,scale,calX,noise,B.true,theta.mag)
      y = Y[,j]
      
      #Regression
      ls2 = trace.reg(n,d1,d2,B.max,SlopeT,A,calX,y,L,t.final)
      MSE[j] = MSE[j] + norm.sq( ls2$B - B.true ) 
      Mvalue[,j] = Mvalue[,j] + as.matrix(ls2$value)
    }
    index = index + 1
  }
  
  MSE = MSE/N.rep
  Mvalue = Mvalue/N.rep
  return(list(MSE = MSE, Mvalue = Mvalue))
}

#-------------------------------------------------------------------------
# For specific parameters plot MSE together for rank r = 1, 2, 3 
#-------------------------------------------------------------------------
n = 80
d1 = 10
d2 = 10
p=d1*d2
step = 1
scale = 2
SlopeT = T
A = 10
B.mag = 1
theta.mag = 10
L = 2
t.final = 30
N.rep = 100

aux = (n/step)/scale #must be an integer

# Rank r = 1
r = 1
reg = trace.reg.MSE(n,d1,d2,r,SlopeT,A,step,scale,B.mag,theta.mag,L,t.final,N.rep)
dfA = reg$MSE

# Rank r = 5
r = 5
reg = trace.reg.MSE(n,d1,d2,r,SlopeT,A,step,scale,B.mag,theta.mag,L,t.final,N.rep)
dfB = reg$MSE

# Rank r = 10
r = 10
reg = trace.reg.MSE(n,d1,d2,r,SlopeT,A,step,scale,B.mag,theta.mag,L,t.final,N.rep)
dfC = reg$MSE

#All 3 plots
epsilon = rep(0,aux+1)
for (j in  1:aux) {
  epsilon[j+1] = step*j/n
}
epsilon.sel = data.frame( epsilon = epsilon )


# Aggregating data points in single data frame factored by the 'rank'
dfA.sel = dfA[ (epsilon.sel > 0.04) & (epsilon.sel <= 0.5) ]
dfB.sel = dfB[ (epsilon.sel > 0.04) & (epsilon.sel <= 0.5) ]
dfC.sel = dfC[ (epsilon.sel > 0.04) & (epsilon.sel <= 0.5) ]
epsilon.sel = epsilon.sel[ (epsilon.sel > 0.04) & (epsilon.sel <= 0.5) ]
SEcomb = c( dfA.sel, dfB.sel, dfC.sel )
epscomb = c( epsilon.sel, epsilon.sel, epsilon.sel )
aux = dim(as.array(epsilon.sel))
rank = as.factor( c( rep(1, aux ), rep( 5, aux ), rep( 10, aux )  ) )
df <- data.frame( SE = SEcomb, epsilon = epscomb, rank = rank)

library(ggplot2)
library(extrafont)
loadfonts(device = "win")
pl <- ggplot( data = df, aes(x = epsilon, y = SE, color = rank, shape = rank ))
pl <- pl + geom_point(size = 4, alpha = 1) 
pl <- pl + geom_smooth(method = 'lm',se=FALSE, linetype="dashed")
pl <- pl + theme_bw()
pl <- pl + theme(text=element_text(family="serif", size =30)) + labs( x =  expression(epsilon), y = 'MSE')
