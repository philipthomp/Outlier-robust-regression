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
data_genB <- function(n,d1,d2,r,B.mag){
  library(matlib)
  library(mclust)
  
  #Low-rank parameter
  U = randomOrthogonalMatrix(d1, r)
  V = randomOrthogonalMatrix(d2, r)
  D = diag(B.mag,r)
  B.true = U%*%D%*%t(V)
  B.max = max(abs(B.true))
  scale = B.mag/sqrt(n) 
  B.true = (scale/B.max)*B.true # r-rank parameter will have sup norm less than B.mag/sqrt(n)
  
  return(B.true)
}

#--------------------------------------------------------------------------------
# Function to assign matrix entry from a index vector (columnwise) and magnitude
# This will be used to assign randomly nonzero entries to sparse matrix parameter
#--------------------------------------------------------------------------------
matrix.entry <- function(d1,d2,indexv, Gamma.mag){
  Gamma = matrix(0, nrow = d1, ncol = d2)
  for (i in indexv) {
    column = ceiling(i/d1)
    row = i%%d1
    if(row == 0){
      row = d1
    }
    Gamma[row,column] = Gamma.mag
  }
  return(Gamma)
}
#--------------------------------------------------------------------------------

#-----------------------------------------------------------
# Function to produce true sparse parameter
#-----------------------------------------------------------
data_genGamma <- function(d1,d2,s,Gamma.mag,random){
  p = d1*d2
  #Sparse parameter
  if(random == T){
    indexv = sample(1:p, size = s, replace = F)
    Gamma.true = matrix.entry(d1,d2,indexv,Gamma.mag)
  }
  else{
    indexv = 1:s
    Gamma.true = matrix.entry(d1,d2,indexv,Gamma.mag)
  }
  
  return(Gamma.true)
}

#---------------------------------------------------------------------
# Function to produce clean response y given (calX, noise, parameter)
#---------------------------------------------------------------------
data_genY <- function(n,calX,noise,B.true,Gamma.true){
  
  y = rep(0,n)
  for (i in 1:n) {
    y[i] = sum(diag( t(calX[[i]]) %*% (B.true + Gamma.true) )) + noise[i]
  } 
  
  return(y)
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
# Soft-thresholding function contrained on the l-infty ball of radius a
#----------------------------------------------------------------------------------
soft.linfty <- function(x,lambda,a){ sign(x)*min( max(abs(x)-lambda,0), a) }

#----------------------------------------------------------------------------------
# Low-rank trace regression solver given (calX,y) and (sd,B.mag)
#----------------------------------------------------------------------------------
trace.reg.MD <- function(n,d1,d2,B.mag,sd,calX,y,L,t.final){
  
  #Low-rank parameter regularization tuning
  lambda = 1*sd*(sqrt(d1/n)+sqrt(d2/n))
  
  #Sparse parameter regularization tuning
  p=d1*d2
  tau = 1*sd*sqrt(log(p)/n) + 1*(B.mag/sqrt(n))
  
  #Iinitialization
  B = matrix(0, nrow = d1, ncol = d2)
  Gamma = matrix(0, nrow = d1, ncol = d2)
  calXBGamma = rep(0, n)
  t = 1
  value = rep(0,t.final)
  
  #Proximal iterations
  for (t in 1:t.final) {
    for (i in 1:n) {
      calXBGamma[i] = sum(diag(t(calX[[i]])%*%(B + Gamma))) #trace inner product for ith data point 
    }
    aux = calXBGamma - y
    value[t] = (1/(2*n))*norm.sq(aux)
    
    #Common gradient computation
    grad = matrix(rep(0,p), nrow = d1)
    for (i in 1:n) {
      grad = grad + (1/n)*aux[i]*calX[[i]]
      i = i + 1
    }
    
    #Low-rank parameter proximal step
    W = B - (1/L)*grad
    B.aux = svd(W, nu = min(d1, d2), nv = min(d1,d2))
    D = B.aux$d
    for (j in 1:min(d1, d2)) {
      D[j] = soft.linfty( D[j], lambda/L, B.mag/sqrt(n) )
    }
    B = (B.aux$u)%*%diag(D)%*%t(B.aux$v)
    
    #Sparse parameter proximal step
    V = Gamma - (1/L)*grad
    for (j in 1:p) {
      Gamma[j] = soft( V[j], tau/L )
    }
    
    t = t + 1 
  }
  
  ls.trace.reg = list(B = B, Gamma = Gamma, value = value)
  return(ls.trace.reg)
}

#------------------------------------------------------------------------------------------------------------------
# Function to call data, do trace regression with matrix decomposition for each rank on vector vr 
# over grid of sparsity for the sparse parameter. 
# Compute and plot MSE versus Sparsity averaged over N.rep repetitions for each rank in vr. 
# Given repetition, the same data is used across all corruptions. 
# For each sparsity s, the same sparse parameter is used across repetitions. 
#------------------------------------------------------------------------------------------------------------------
trace.reg.MD.MSE.sparse <- function(n,d1,d2,vr,S,step,random,sd,B.mag,Gamma.mag,L,t.final,N.rep){
  
  p = d1*d2
  R = dim(as.array(vr))
  
  # Save low-parameter
  elm = matrix(rep(0,p), nrow = d1)
  B.true = rep(list(elm), R)
  for (j in 1:R) {
    B.true[[j]] = data_genB(n,d1,d2,vr[j],B.mag)
  }
  
  # Save grid of sparse parameters
  elm = matrix(rep(0,p), nrow = d1)
  K = (S-1)/step # number of steps in grid of sparsity initializing at 1. Total number of points is K+1
  Gamma.true = rep(list(elm), K+1)
  Gamma.true[[1]] = data_genGamma(d1,d2,1,Gamma.mag,random)
  for (k in 1:K) {
    s = 1 + k*step
    Gamma.true[[k+1]] = data_genGamma(d1,d2,s,Gamma.mag,random)
  }
  
  # Initialization of repetition count
  elm = rep(0,K+1)
  MSE = rep(list(elm), R)
  elm = matrix( rep(0, t.final*(K+1) ) , nrow = t.final )
  Mvalue = rep(list(elm), R)
  
  for (index in 1:N.rep) {
    # calX and noise generation
    ls0 <- data_genXe(n,d1,d2,sd)
    calX = ls0$calX
    noise = ls0$noise
    
    for (j in 1:R){
      for (k in 1:(K+1)) {
        # Parameter call and response generation
        y = data_genY(n,calX,noise,B.true[[j]],Gamma.true[[k]])
        
        #Regression
        ls1 = trace.reg.MD(n,d1,d2,B.mag,sd,calX,y,L,t.final)
        MSE[[j]][k] = MSE[[j]][k] + norm.sq( ls1$B - B.true[[j]] ) + norm.sq( ls1$Gamma - Gamma.true[[k]] )
        Mvalue[[j]][,k] = Mvalue[[j]][,k] + as.matrix(ls1$value)
      }
    }  
    
    index = index + 1
  }
  
  for (j in 1:R) {
    MSE[[j]] = MSE[[j]]/N.rep
    Mvalue[[j]] = Mvalue[[j]]/N.rep
  }
  
  return(list(MSE = MSE, Mvalue = Mvalue))
}

#-----------------------------------------------------------------------------------
# For specific parameters plot MSE together for sparsity s over grid of rank 1:R
#-----------------------------------------------------------------------------------
n = 1000
d1 = 10
d2 = 10
step = 5
S = 51
K = (S-1)/step #must be natural number
B.mag = 1
Gamma.mag = 10
random = T
sd = 1
L = 5
t.final = 50
N.rep = 20

#---------------------------------------------------------------------------------------
# 1 PLOT:
#---------------------------------------------------------------------------------------
vr = c(5)
ls1 <- trace.reg.MD.MSE.sparse(n,d1,d2,vr[1],S,step,random,sd,B.mag,Gamma.mag,L,t.final,N.rep)

#Aggregating data points in single data frame factored by the 'sparsity'
MSEcomb = c( ls1$MSE[[1]])
v = seq(1, S, step)
Scomb = c( v )
Rank = as.factor( c( rep(vr[1], K+1 )  ) )
df <- data.frame( MSE = MSEcomb, Sparsity = Scomb, Rank = Rank)


#---------------------------------------------------------------------------------------
# 2 PLOTS:
#---------------------------------------------------------------------------------------
#vr = c(1,5)
#ls1 <- trace.reg.MD.MSE.sparse(n,d1,d2,vr[1],S,step,random,sd,B.mag,Gamma.mag,L,t.final,N.rep)
#ls2 <- trace.reg.MD.MSE.sparse(n,d1,d2,vr[2],S,step,random,sd,B.mag,Gamma.mag,L,t.final,N.rep)

#Aggregating data points in single data frame factored by the 'sparsity'
#MSEcomb = c( ls1$MSE[[1]], ls2$MSE[[1]])
#v = seq(1, S, step)
#Scomb = c( v, v )
#Rank = as.factor( c( rep(vr[1], K+1 ), rep( vr[2], K+1 )  ) )
#df <- data.frame( MSE = MSEcomb, Sparsity = Scomb, Rank = Rank)


#---------------------------------------------------------------------------------------
# 3 PLOTS:
#---------------------------------------------------------------------------------------
#vr = c(1,5,10)
#ls1 <- trace.reg.MD.MSE.sparse(n,d1,d2,vr[1],S,step,random,sd,B.mag,Gamma.mag,L,t.final,N.rep)
#ls2 <- trace.reg.MD.MSE.sparse(n,d1,d2,vr[2],S,step,random,sd,B.mag,Gamma.mag,L,t.final,N.rep)
#ls3 <- trace.reg.MD.MSE.sparse(n,d1,d2,vr[3],S,step,random,sd,B.mag,Gamma.mag,L,t.final,N.rep)

#Aggregating data points in single data frame factored by the 'sparsity'
#MSEcomb = c( ls1$MSE[[1]], ls2$MSE[[1]], ls3$MSE[[1]] )
#v = seq(1, S, step)
#Scomb = c( v, v, v )
#Rank = as.factor( c( rep(vr[1], K+1 ), rep( vr[2], K+1 ), rep( vr[3], K+1 )  ) )
#df <- data.frame( MSE = MSEcomb, Sparsity = Scomb, Rank = Rank)


library(ggplot2)
library(extrafont)
loadfonts(device = "win")
pl <- ggplot( data = df, aes(x = Sparsity, y = MSE, color = Rank, shape = Rank ))
pl <- pl + geom_point(size = 4, alpha = 1) 
pl <- pl + geom_smooth(method = 'lm',se=FALSE, linetype="dashed")
pl <- pl + theme_bw()
pl <- pl + theme(text=element_text(family="serif",size=30)) + labs( x =  'Sparsity', y = 'MSE')



