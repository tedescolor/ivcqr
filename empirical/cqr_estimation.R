path = getwd() #set path to the folder simulations/
setwd(path)
test = FALSE # In this script, this variable is defined only for consistency with the replication package, but it is not used since the script finishes in short time
{
  # Package names
  packages <-
    c("survival",
      "foreach",
      "doParallel",
      "quantreg",
      "dplyr")
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages],repos='http://cran.us.r-project.org')
  }
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
} #LOAD PACKAGES

tauhat.func <- function(y0, x0, z, x, delta,h){
  # tau0(y0, x0) = F(T<y0|x0); so y0 is the C_i, and x0 is the xi in the paper
  # z is observed vector of response variable
  # x is the observed covariate
  # delta is the censoring indicator function
  # h is the bandwidth
  n<-length(z)
  ###kernel weights#########################
  Bn = Bnk.func(x0, x, h)
  if (y0<max(z))
  {
    # sort the data z, and the delta, Bn correspondingly to the order of sorted y
    z2 = sort(z)
    Order = order(z) # so z[Order] = z2
    Bn2 = Bn[Order]
    delta2 = delta[Order]
    eta = which(delta2==1 & z2<=y0) # the index of those observations satisfying delta2==1 & z2<=y0
    Bn3 = Bn2[n:1]  # change the order of Bn2, make the first obs of Bn2 to be the last of Bn3
    tmp = 1- Bn2 /cumsum(Bn3)[n:1]  
    out = 1-prod(tmp[eta], na.rm=T) # na.rm=T, as some of those tmp=NA as the denom =0
  } 
  else out<-1 
  return(out)
}
# Wang and Wang (2009) (Code from webpage Huixia Wang)
Bnk.func    <- function(x0, x, h){
  # the kernel weight function Bnk(x0, x), where x0 is a scalar, and x is a vector
  # returns a vector
  # h is the bandwidth
  xx<-(x-x0)/h  
  xx[abs(xx)>=1]<-1
  w<-15*(1-xx^2)^2/16  #biquadratic kernel 
  w<-w/sum(w)
  return(w)
} 
# MM version of Wang and Wang's estimator 
MMLQR.WW    <- function(y,X,delta,tau,beta,h,toler=1e-8,maxit=5000){
  # MM version of Wang and Wang's estimator 
  iteration <- 0
  n.obs     <- length(y)
  df.mat    <- cbind(rep(1,n.obs),X)
  
  # Calculation of epsilon
  tn        <- toler/n.obs
  e0        <- -tn/log(tn)
  eps       <- (e0-tn)/(1+log(e0))
  
  # Calculation of weights \hat{w}
  ind <- which(delta==0)
  w   <- rep(1, n.obs)  # the weight
  
  for(i in 1:length(ind))
  {
    x0 = X[ind[i]]
    y0 = y[ind[i]]
    tau.star = tauhat.func(y0,x0, y, X, delta,h=h)
    if (tau>tau.star) w[ind[i]] = (tau-tau.star)/(1-tau.star)
  }
  
  # pseudo observations
  ind2       <- which(w!=1)
  y.pse      <- rep(max(y)+100, length(ind2))
  x.pse      <- X[ind2,]
  df.mat.pse <- cbind(rep(1,length(ind2)),x.pse)
  
  yy    <- c(y, y.pse)
  xx    <- c(X, x.pse)
  ww    <- c(w, 1-w[ind2])
  
  # Initialization of condition for break
  cond <- T
  while(cond)
  { 
    beta.prev <- beta
    r.vec     <- y-df.mat%*%as.matrix(beta)
    r.vec.pse <- y.pse-df.mat.pse%*%as.matrix(beta)
    
    A.entries     <- c(w/(eps+abs(r.vec)))/2
    A.mat         <- diag(A.entries)
    A.entries.pse <- c((1-w[ind2])/(eps+abs(r.vec.pse)))/2
    A.mat.pse     <- diag(A.entries.pse)
    
    B.vec     <- as.matrix(rep(tau-.5,n.obs))
    
    beta      <- solve(t(df.mat)%*%A.mat%*%df.mat+t(df.mat.pse)%*%A.mat.pse%*%df.mat.pse)%*%(t(df.mat)%*%(A.mat%*%as.matrix(y)+B.vec)+t(df.mat.pse)%*%A.mat.pse%*%as.matrix(y.pse))
    cond      <- max(abs(beta-beta.prev)) > toler
    iteration <- iteration + 1
    if(iteration > maxit){warning("WARNING: Algorithm did not converge"); break} 
  }
  return(list("beta"=c(beta),"IterN"=iteration))
} 

df = read.csv("data.csv")

df = df[complete.cases(df),]  #remove rows with missing data
print(paste(nrow(df), "rows"))

df = df[ df$white == 0 & df$male == 0 & df$single==1 & df$children==1,]
table(df$jtpa, df$treatment)/nrow(df)
us =  seq(0.1,0.9,by=0.1)
#portnoy
res = crq(formula = Surv(log(df$days),df$delta,type = "right")~treatment+age,
          taus =us,data = df,method = "Portnoy")
coef(res,taus = us)
summary(res,taus = us)

# wang wang
beta= c()
for(u in us){
beta = rbind(beta, 
             c(u,MMLQR.WW(as.numeric(log(df$days)), cbind(df$treatment, df$age), as.numeric(df$delta), 
                          beta = c(0,0,0),tau=u, h=.05)$beta)) 
}
beta = as.data.frame(beta)
names(beta) = c("u","intercept", "treatment", "age")
write.csv(beta, "cqr.csv",row.names = FALSE)
