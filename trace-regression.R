#-------------------------------------------------------------------------------------
# Function to generate design operator calX (as a n-list of d1Xd2-matrices) and noise
#-------------------------------------------------------------------------------------
data_genXe <- function(n,d1,d2,sd){
  p = d1*d2
  elm = matrix(rep(0,p), nrow = d1)
  calX = rep(list(elm), n)
  noise = rep(0,n)
  for (i in 1:n) {
    calX[[i]] = matrix(rnorm(p, mean = 0, sd = 1), nrow = d1, ncol = d2)
    noise[i] = rnorm(1, mean = 0, sd = sd)
  }
  list = list(calX = calX, noise = noise)
  return(list)
}

#-----------------------------------------------------------
# Function to produce true low-rank parameter
#-----------------------------------------------------------
data_genB <- function(d1,d2,r,B.mag){
  library(matlib)
  library(mclust)
  
  #Low-rank parameter
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
  
  #Produce clean response
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
soft <- function(x,lambda){ sign(x)*max( abs(x) - lambda, 0) }

#----------------------------------------------------------------------------------
# Low-rank trace regression solver given ( calX, y )
#----------------------------------------------------------------------------------
trace.reg <- function(n,d1,d2,SlopeT,A,sd,calX,y,L,t.final){
  
  #call for prox_sorted_L1()
  library(glmSLOPE)
  
  #Parameter regularization tuning
  lambda = 4*sd*(sqrt(d1/n)+sqrt(d2/n))
  
  #Corruption regularization tuning
  if (SlopeT == T){
    tau = 10*sd*sqrt(1/n)
    omega = rep(0,n)
    for (i in 1:n) {
      omega[i] = log(A*n/i)
    }
  } 
  else{
    tau = 10*sd*sqrt(log(n)/n)
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
      calX.B[i] = sum(diag(t(calX[[i]])%*%B))
    }
    aux = calX.B + sqrt(n)*theta-y
    value[t] = (1/(2*n))*norm.sq(aux)
    
    #Parameter proximal iteration
    gradB = matrix(rep(0,d1*d2), nrow = d1)
    for (i in 1:n) {
      gradB = gradB + (1/n)*aux[i]*calX[[i]]
      i = i + 1
    }
    W = B - (1/L)*gradB
    B.aux = svd(W, nu = min(d1, d2), nv = min(d1,d2))
    D = B.aux$d
    for (j in 1:min(d1, d2)) {
      D[j] = soft(D[j], lambda/L)
    }
    B = (B.aux$u)%*%diag(D)%*%t(B.aux$v)
    
    #Corruption proximal iteration
    gradtheta = (1/sqrt(n))*aux
    w = theta - (1/L)*gradtheta
    aux1 = (tau/L)*omega
    theta =  as.matrix(prox_sorted_L1(w,aux1))
    
    t = t + 1 
  }
  
  ls.trace.reg = list(B = B, value = value)
  return(ls.trace.reg)
}

#------------------------------------------------------------------------------------------------------------------
# Function to call data, do trace regression for each rank on vector vr over grid of corruption fraction.
# Compute and plot sqrt of MSE versus corruption fraction averaged over N.rep repetitions for each rank in vr. 
# Given repetition, the same data is used across all corruptions. 
# For each corruption fraction, the same low-rank parameter is used across repetitions. 
#------------------------------------------------------------------------------------------------------------------
trace.reg.MSE <- function(n,d1,d2,r,SlopeT,A,step,scale,B.mag,theta.mag,L,t.final,sd,N.rep){
  
  # Generate low-rank parameter
  B.true = data_genB(d1,d2,r,B.mag)
  
  # Initialization of repetition count
  aux = (n/step)/scale
  MSE = rep(0,aux + 1)
  Mvalue = matrix( rep(0, t.final*(aux + 1) ) , nrow = t.final)
  
  for (index in 1:N.rep) {
    # calX and noise generation
    ls0 <- data_genXe(n,d1,d2,sd)
    calX = ls0$calX
    noise = ls0$noise
    
    for (j in 1:(aux+1)) {
      Y = data_genY(n,step,scale,calX,noise,B.true,theta.mag)
      y = Y[,j]
      
      #Regression
      ls2 = trace.reg(n,d1,d2,SlopeT,A,sd,calX,y,L,t.final)
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
# For specific parameters plot MSE together for rank r = 2,5,8
#-------------------------------------------------------------------------
n = 1000
d1 = 10
d2 = 10
step = 50
scale = 2
SlopeT = T
A = 10
B.mag = 10
theta.mag = 10
sd = 1
L = 5
t.final = 50
N.rep = 50


aux = (n/step)/scale #must be an integer

# Rank r = 1
r = 1
reg = trace.reg.MSE(n,d1,d2,r,SlopeT,A,step,scale,B.mag,theta.mag,L,t.final,sd,N.rep)
dfA = sqrt(reg$MSE)

# Rank r = 2
r = 2
reg = trace.reg.MSE(n,d1,d2,r,SlopeT,A,step,scale,B.mag,theta.mag,L,t.final,sd,N.rep)
dfB = sqrt(reg$MSE)

# Rank r = 3
r = 3
reg = trace.reg.MSE(n,d1,d2,r,SlopeT,A,step,scale,B.mag,theta.mag,L,t.final,sd,N.rep)
dfC = sqrt(reg$MSE)

#Contamination variable 
epsilon = rep(0,aux+1)
for (j in  1:aux) {
  epsilon[j+1] = step*j/n
}
epsilon.sel = data.frame( epsilon = epsilon )


# Aggregating data points in single data frame factored by the 'rank'
SEcomb = c( dfA[epsilon.sel <= 0.5], dfB[epsilon.sel <= 0.5], dfC[epsilon.sel <= 0.5] )
epscomb = c( epsilon.sel[epsilon.sel <= 0.5], epsilon.sel[epsilon.sel <= 0.5], epsilon.sel[epsilon.sel <= 0.5] )
aux = dim(as.array(epsilon.sel[epsilon.sel <= 0.5]))
rank = as.factor( c( rep(1, aux ), rep( 2, aux ), rep( 3, aux )  ) )
df <- data.frame( SE = SEcomb, epsilon = epscomb, rank = rank)

library(ggplot2)
library(extrafont)
loadfonts(device = "win")
pl <- ggplot( data = df, aes(x = epsilon, y = SE, color = rank, shape = rank ))
pl <- pl + geom_point(size = 4, alpha = 1) 
pl <- pl + geom_smooth(method = 'lm',se=FALSE, linetype="dashed")
pl <- pl + theme_bw()
pl <- pl + theme(text=element_text(family="serif", size = 30)) + labs( x =  expression(epsilon), y = 'Square root of MSE')



