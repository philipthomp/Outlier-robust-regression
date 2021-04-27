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
# for different levels of corruption in vector vl. 
# Then compute MSE averaged over N.rep repetitions. Given repetition, the same data is used across all corruptions
#------------------------------------------------------------------------------------------------------------------
trace.reg.MSE <- function(n,d1,d2,r,SlopeT,A,step,scale,B.mag,vl,L,t.final,N.rep){
  
  VL = dim(as.array(vl))
  
  # Generate low-rank parameter
  B.true = data_genP(n,d1,d2,r,B.mag)
  B.max = max(abs(B.true))
  
  # Initialization of repetition count
  aux = (n/step)/scale
  elm = rep(0,aux + 1)
  MSE = rep(list(elm), VL)
  elm = matrix( rep(0, t.final*(aux + 1) ) , nrow = t.final)
  Mvalue = rep(list(elm), VL)
  
  for (index in 1:N.rep) {
    # calX and noise generation
    ls0 <- data_genXe(n,d1,d2)
    calX = ls0$calX
    noise = ls0$noise
    
    for (l in 1:VL){
      for (j in 1:(aux+1)) {
        # Parameter and response generation
        Y = data_genY(n,step,scale,calX,noise,B.true,vl[l])
        y = Y[,j]
      
        #Regression
        ls2 = trace.reg(n,d1,d2,B.max,SlopeT,A,calX,y,L,t.final)
        MSE[[l]][j] = MSE[[l]][j] + norm.sq( ls2$B - B.true ) 
        Mvalue[[l]][,j] = Mvalue[[l]][,j] + as.matrix(ls2$value)
      }
    }
    
    index = index + 1
  }
  
  for (l in 1:VL) {
    MSE[[l]] = MSE[[l]]/N.rep
    Mvalue[[l]] = Mvalue[[l]]/N.rep
  }
  
  return(list(MSE = MSE, Mvalue = Mvalue))
}

#-----------------------------------------------------------------------------------------
# For specific parameters plot MSE together for theta.mag = 10, 100, 1000 and rank = 3
#-----------------------------------------------------------------------------------------
n = 80
d1 = 10
d2 = 10
p=d1*d2
step = 4
scale = 2
SlopeT = T
A = 10
B.mag = 1
L = 2
t.final = 30
N.rep = 100
aux = (n/step)/scale 

r = 3

#------------------------------------------------------------------------------------------
# 3 different MSE for different levels 10, 100, 1000 using independent data sets
#------------------------------------------------------------------------------------------
# theta.mag = 10
#theta.mag = 10
#reg = trace.reg.MSE(n,d1,d2,r,SlopeT,A,step,scale,B.mag,theta.mag,L,t.final,N.rep)
#MSE = reg$MSE
#SE = sqrt(MSE[2:aux])
#SE.sel = SE[3:aux]
#dfA = data.frame(SE.A = SE.sel)

# theta.mag = 100
#theta.mag = 100
#reg = trace.reg.MSE(n,d1,d2,r,SlopeT,A,step,scale,B.mag,theta.mag,L,t.final,N.rep)
#MSE = reg$MSE
#SE = sqrt(MSE[2:aux])
#SE.sel = SE[3:aux]
#dfB = data.frame(SE.B = SE.sel)

# theta.mag = 1000
#theta.mag = 1000
#reg = trace.reg.MSE(n,d1,d2,r,SlopeT,A,step,scale,B.mag,theta.mag,L,t.final,N.rep)
#MSE = reg$MSE
#SE = sqrt(MSE[2:aux])
#SE.sel = SE[3:aux]
#dfC = data.frame(SE.C = SE.sel)
#-------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# 3 different MSE for different levels 10, 100, 1000 using same data set
#------------------------------------------------------------------------------------------
vl = c(10, 100, 1000)
ls <- trace.reg.MSE(n,d1,d2,r,SlopeT,A,step,scale,B.mag,vl,L,t.final,N.rep)
SE.A = sqrt(ls$MSE[[1]][1:aux+1])
SE.B = sqrt(ls$MSE[[2]][1:aux+1])
SE.C = sqrt(ls$MSE[[3]][1:aux+1])


#Continuing aggregation
#Corruption grid
epsilon = rep(0,aux)
for (j in  1:aux) {
  epsilon[j] = step*j/n
}

SEcomb = c( SE.A[epsilon <= 0.5],  SE.B[epsilon <= 0.5], SE.C[epsilon <= 0.5] )
epscomb = c( epsilon[epsilon <= 0.5], epsilon[epsilon <= 0.5], epsilon[epsilon <= 0.5] )
aux = dim(as.array(epsilon[epsilon <= 0.5]))
level = as.factor( c( rep(10, aux ), rep( 100, aux ), rep( 1000, aux )  ) )
df <- data.frame( SE = SEcomb, epsilon = epscomb, level = level)

library(ggplot2)
library(extrafont)
loadfonts(device = "win")
pl <- ggplot( data = df, aes(x = epsilon, y = SE, color = level, shape = level ))
pl <- pl + geom_point(size = 4, alpha = 1) + ylim(1,4)
pl <- pl + geom_smooth(method = 'lm',se=FALSE, linetype="dashed")
pl <- pl + theme_bw()
pl <- pl + theme(text=element_text(family="serif")) + labs( x =  expression(epsilon), y = 'Square root of MSE', title = 'Robust matrix completion: error versus contamination')
