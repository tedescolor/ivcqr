path = getwd() #set path to the folder simulations/
setwd(path)
test = FALSE
if(test){
  ##### test code setting
  ns = c(100)#sample size 
  Nsimulations = 2#number of simulations
  us = c(0.3)#u of interest
}else{
  ### simulation paper setting
  ns = c(500,1000)#sample size
  Nsimulations = 500#number of simulations
  us = c(0.3,0.5,0.7)#u of interest
}

{
  # Package names
  packages <-
    c('survival',
      'foreach',
      'doParallel',
      'quantreg',
      'dplyr')
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages],repos='http://cran.us.r-project.org')
  }
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
} #LOAD PACKAGES


createSimulationDataFunction = #create the simulationData(n,beta,ucoeffs), given the distribution of W,X1,X2
  function(Wdensity, #density of W (as string)
           X1density,#density of X1 (as string)
           X2density, #density of X2 (as string)
           Cdensity #density of C (as string)
  ){
    return(
      function(n,
               #definition of ucoeffs: beta(u) = (beta1,beta2) + ucoeffs * u
               verbose = FALSE #print intermidiate results
      ){ # simulate data 2 covariate
        Ws = eval(parse(text=Wdensity)) # IV
        Us = runif(n) # U variable
        Xs = cbind(1,eval(parse(text=X1density)), eval(parse(text=X2density))) #Covariate
        # OSS: we consider beta(u) = (beta1 + ucoeff1* u, beta2 + ucoeff2 * u)
        Ts <- exp(  apply( cbind(Us, Us, Us) * Xs,1,sum) ) #Times variable
        Cs = eval(parse(text=Cdensity))
        Ds = as.numeric(Ts <= Cs) # Delta
        # OSS: the supp of C consists of all real numbers, therefore the model is identifiable.
        Ys = pmin(Ts,Cs) # Observed Times
        if(verbose){
          par(mfrow=c(3,1)); hist(Ts);hist(Cs);hist(Ys);
          print(paste("censor = ", sum(Ds == 0)/n, sep = ""))}
        data = as.data.frame(cbind(Ys,Ds,as.data.frame(Xs),Ws))
        names(data) = c("y", "delta", "x0","x1", "x2", "w")
        return(data)
      }
    )
  }


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


mat = matrix(nrow = 0, ncol = 10)
alldensities = 
  rbind(
    c("rexp(n)", "0.5*Us + Ws -1 > 0", "runif(n)","rexp(n, 0.176)"), 
    c("rexp(n)", "0.5*Us + Ws -1 > 0", "runif(n)","rexp(n, 0.069)"),
    c("as.numeric(runif(n))>0.5", "0.5*Us + Ws -1 > 0", "runif(n)","rexp(n, 0.175)"),
    c("as.numeric(runif(n))>0.5", "0.5*Us + Ws -1 > 0", "runif(n)","rexp(n, 0.07)"),
    c("exp(rnorm(n))", "0.5*Us + Ws + 0.2*runif(n)", "rexp(n)","rexp(n, 0.065)"),
    c("exp(rnorm(n))", "0.5*Us + Ws + 0.2*runif(n)", "rexp(n)","rexp(n, 0.0173)"))



covNames = c("x0","x1","x2"); IVNames = c("w","x2");
allresults = c()
for(u in us){
  for(i in 1:nrow(alldensities)){
    for(n in ns){
      densities = alldensities[i,]
      current = matrix(nrow = 0, ncol = 3)
      betau = c(u,u,u)
      print(paste(u,i,n))
      Cdensity= densities[4];
      cores=detectCores()
      cl <- makeCluster(cores-1) #not to overload your computer
      registerDoParallel(cl)
      results <- foreach(j=1:Nsimulations, .combine=rbind, .packages = packages) %dopar% {
        set.seed(j)
        Wdensity = densities[1]; X1density = densities[2];X2density = densities[3]; Cdensity= densities[4];
        simulationData = createSimulationDataFunction(Wdensity = Wdensity,
                                                      X1density = X1density,
                                                      X2density = X2density,
                                                      Cdensity = Cdensity)
      
      data = simulationData(n)
      beta = MMLQR.WW(as.numeric(log(data$y)), cbind(data$x1, data$x2), as.numeric(data$delta), beta = c(0,0,0),tau=u, h=.05)$beta 
      bdata = data[sample.int(n,n,replace = TRUE),]
      bbeta = MMLQR.WW(as.numeric(log(bdata$y)), cbind(bdata$x1, bdata$x2), as.numeric(bdata$delta), beta = c(0,0,0),tau=u, h=.05)$beta
      return(c(beta,bbeta))
      }
      stopCluster(cl)
      allresults = rbind(allresults, cbind(Cdensity,n,u,results) )
      betau= c(u,u,u)
      rmse = sqrt(mean(apply(results[,1:3],1,function(x){sum((x-betau)^2)})))
      CP95 = array(NA,dim = 3)
      for(comp in 1:3){
        CP95[comp] = sum( (results[,comp] - u) >quantile( results[,comp+3] - results[,comp],0.025) & 
                            (results[,comp] - u) <quantile( results[,comp+3] - results[,comp],0.975) )/Nsimulations
      }
      mat = rbind(mat, c(colMeans(results[,1:3]-betau), Cdensity,n,u,rmse,CP95 ))
    }
  
  }
}
mat = as.data.frame(mat)
betanames = c("beta0","beta1","beta2")
names(mat) = c(sapply(betanames, function(x){paste("crq_bias_", x, sep = "")}), "Cdensity", "n","u", "rmse",sapply(betanames, function(x){paste("CP95_", x, sep = "")}))
mat

ids = sapply(mat$Cdensity, function(c){which(alldensities[,4] == c)}) #CDesity used as unique identifier
mat$Wdensity = alldensities[ids,1]
mat$X1density = alldensities[ids,2]
mat$X2density = alldensities[ids,3]
mat$censor = c(40,20,40,20,40,20)[ids]
head(mat)
#construct columns referring to the alternative estimator of table 2 of the paper
table2 = c()
for(Wdensity in c("rexp(n)","exp(rnorm(n))","as.numeric(runif(n))>0.5")){
  current = mat[mat$Wdensity == Wdensity,]
  current = as.data.frame(current)
  table2 = rbind(table2, current[order(current$u,as.numeric(current$n),current$censor),])
}
table2 = table2[,c("Wdensity","u","n","censor",
                   sapply(betanames, function(x){paste("crq_bias_", x, sep = "")}),
                   "rmse")]
table2 <- cbind(table2$Wdensity,as.data.frame(lapply(table2[,2:ncol(table2)], as.numeric)))
options("digits" = 3)
table2
write.csv(table2, file = "cqr_results.csv", row.names = FALSE)
