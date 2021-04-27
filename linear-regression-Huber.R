#-------------------------------------------------------------------------------------
# Function to generate design matrix X and noise
#-------------------------------------------------------------------------------------
data_genXe <- function(n,p,sd){
  
  X = matrix(rep(0,n*p), nrow = n)
  noise = rep(0,n)
  for (i in 1:n) {
    X[i,] = rnorm(p, mean = 0, sd = 1)
    noise[i] = rnorm(1, mean = 0, sd = sd)
  }
  list = list(X = X, noise = noise)
  return(list)
}

#-------------------------------------------------------------------------------------
# Function to generate sparse parameter
#-------------------------------------------------------------------------------------
data_genP <- function(p,s,b.mag){
  
  b.true = rep(0,p)
  for (j in 1:s) {
    b.true[j] = b.mag
  }
  
  return(b.true)
}

#---------------------------------------------------------------------
# Function to produce corrupted response y given (X, noise, parameter)
#---------------------------------------------------------------------
data_genY <- function(n,step,scale,X,noise,b.true,theta.mag){
  
  #Produce clean response
  y = rep(0,n)
  for (i in 1:n) {
    y[i] = t(X[i,]) %*% b.true + noise[i]
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
# Sparse linear regression solver given (X,y)
#----------------------------------------------------------------------------------
linear.reg <- function(n,p,SlopeB,s,SlopeT,A1,A,X,y,L,sd,t.final){
  
  #call for prox_sorted_L1()
  library(glmSLOPE)
  
  #Parameter regularization
  if (SlopeB == T){
    lambda = 4*sd*sqrt(1/n)
    omegab = rep(0,p)
    for (j in 1:p) {
      omegab[j] = log(A1*p/j)
    }
  } 
  else{
    lambda = 4*sd*sqrt(log(p/s)/n)
    omegab = rep(1,p)
  } 
  
  #Corruption regularization
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
  b = as.matrix(rep(0,p))
  theta = as.matrix(rep(0,n))
  t = 1
  value = rep(0,t.final)
  
  #Proximal iterations
  for (t in 1:t.final) {
    aux = X%*%b+sqrt(n)*theta-y
    value[t] = (1/(2*n))*norm.sq(aux)
    
    gradb = (1/n)*t(X)%*%aux
    u = b - (1/L)*gradb
    aux1 = (lambda/L)*omegab
    b =  as.matrix(prox_sorted_L1(u,aux1))
    
    gradtheta = (1/sqrt(n))*aux
    v = theta - (1/L)*gradtheta
    aux1 = (tau/L)*omega
    theta =  as.matrix(prox_sorted_L1(v,aux1))
    
    t = t + 1 
  }
  
  ls.linear.reg = list(b = b, value = value)
  return(ls.linear.reg)
}

#------------------------------------------------------------------------------------------------------------------
# Function to call data, do linear regression for each for sparsity s grid of corruption fraction. 
# The vector vt contains two possibilities in {T,F}, corresponding to the SlopeT variable. 
# Compute and plot sqrt of MSE versus corruption fraction averaged over N.rep repetitions for each sparsity in vs. 
# Given repetition, the same data is used across all corruptions. 
# For each corruption fraction, the same sparse parameter is used across repetitions. 
#------------------------------------------------------------------------------------------------------------------
linear.reg.MSE <- function(n,p,SlopeB,s,vt,A1,A,step,scale,b.mag,theta.mag,L,t.final,sd,N.rep){
  
  # Generate sparse parameter
  b.true = data_genP(p,s,b.mag)
  
  # Initialization of repetition count
  aux = (n/step)/scale
  VT = dim(as.array(vt))
  elm = rep(0,aux + 1)
  MSE = rep(list(elm), VT)
  elm = matrix( rep(0, t.final*(aux + 1) ) , nrow = t.final)
  Mvalue = rep(list(elm), VT)
  
  for (index in 1:N.rep) {
    
    # calX and noise generation
    ls0 <- data_genXe(n,p,sd)
    X = ls0$X
    noise = ls0$noise
    
    # Generate corrupted response
    Y = data_genY(n,step,scale,X,noise,b.true,theta.mag)
    
    for (l in 1:2){
      for (j in 1:(aux+1)) {
        y = Y[,j]
      
        #Regression
        ls2 = linear.reg(n,p,SlopeB,s,vt[l],A1,A,X,y,L,sd,t.final)
        MSE[[l]][j] = MSE[[l]][j] + norm.sq( ls2$b - b.true ) 
        Mvalue[[l]][,j] = Mvalue[[l]][,j] + as.matrix(ls2$value)
      }
    }  
    
    index = index + 1
  }
  
  for (l in 1:2) {
    MSE[[l]] = MSE[[l]]/N.rep
    Mvalue[[l]] = Mvalue[[l]]/N.rep
  }
  
  
  return(list(MSE = MSE, Mvalue = Mvalue))
}

#----------------------------------------------------------------------------------------------------
# For specific parameters plot MSE together for s = 25 for Huber regression and Sorted Huber loss
#----------------------------------------------------------------------------------------------------
n = 1000
p = 100
step = 50
scale = 2
SlopeB = T
A1 = 10
A = 10
b.mag = 50
theta.mag = 50
L = 5
t.final = 50
sd = 1
N.rep = 100

aux = (n/step)/scale #must be an integer

s = 25
vt = c(T,F)

# Regression
ls <- linear.reg.MSE(n,p,SlopeB,s,vt,A1,A,step,scale,b.mag,theta.mag,L,t.final,sd,N.rep)
dfA = sqrt(ls$MSE[[1]])
dfB = sqrt(ls$MSE[[2]])

# Contamination variable
epsilon = rep(0,aux+1)
for (j in  1:aux) {
  epsilon[j+1] = step*j/n
}
epsilon.sel = data.frame( epsilon = epsilon )


# Aggregating data points in single data frame factored by the type of loss
SEcomb = c( dfA[epsilon.sel <= 0.5], dfB[epsilon.sel <= 0.5] )
epscomb = c( epsilon.sel[epsilon.sel <= 0.5], epsilon.sel[epsilon.sel <= 0.5])
aux = dim(as.array(epsilon.sel[epsilon.sel <= 0.5]))
loss = as.factor( c( rep('Sorted Huber', aux ), rep( 'Huber', aux )  ) )
df <- data.frame( SE = SEcomb, epsilon = epscomb, loss = loss)

library(ggplot2)
library(extrafont)
loadfonts(device = "win")
pl <- ggplot( data = df, aes(x = epsilon, y = SE, color = loss, shape = loss))
pl <- pl + geom_point(size = 4, alpha = 1)
pl <- pl + geom_smooth(method = 'lm',se=FALSE, linetype="dashed")
pl <- pl + theme_bw()
pl <- pl + theme(text=element_text(family="serif",size=30)) + labs( x =  expression(epsilon), y = 'Square root of MSE')



