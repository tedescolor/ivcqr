path = getwd() #set path to the folder simulations/
setwd(path)
test = FALSE
if(test){
##### test code setting
  ns = c(100)#sample size
  Nsimulations = 10#number of simulations
  us = c(0.1)#u of interest
  attempts = 2#attempts
  battempts = 2 #number bootstrap attempts
}else{
  ##### simulation paper setting
  ns = c(500,1000)#sample size
  Nsimulations = 500#number of simulations
  us = c(0.3,0.5,0.7)#u of interest
  attempts = 100#attempts
  battempts = 100 #number bootstrap attempts
}

{
  # Package names
  packages <-
    c("survival",
      "foreach",
      "doParallel",
      "dplyr")
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages],repos='http://cran.us.r-project.org')
  }
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
} #LOAD PACKAGES


#################################
#### ESTIMATION FUNCTIONS
#################################
KMfunction = function(KMsurvfit){ # from survfit KM to function KM
  return(function(t){
    #if(t < min(KMsurvfit$time)){return(1.0)}
    index = sum( KMsurvfit$time -rep(t,length(KMsurvfit$time)) <= 0 )
    return(KMsurvfit$surv[index]) # return the $surv value
  })
}


createf = function(data,u,deltaName,timeName, covNames, IVNames, G){ return( 
  function(b){
    n = nrow(data)
    ys = data[,timeName]
    Zs = data.matrix(data[,covNames])
    deltas = data$delta
    Gs = sapply(ys, G)
    Gs[Gs == 0] = 1 
    w1 = array(NA,dim = c(n,n,length(IVNames)))
    w2 = array(NA,dim = c(n,n,length(IVNames)))
    for(i in 1:n){
      w1[,i,] = data.matrix(data[,IVNames])
      w2[i,,] = data.matrix(data[,IVNames])
    }
    W = array(as.numeric(apply(w1 <= w2,c(1,2),all)),dim = c(n,n))
    dG = deltas/Gs
    ab = t(
      matrix(
        as.numeric(ys <= exp(Zs %*% b))*dG - u,
        nrow = n, ncol = 1)) %*% W
    return( ab %*% t(ab))# / (n^3))# useless in terms of minimization to do the last division
  })
}






estimate = function(data, #dataset
                    u, #quantile
                    method = "Brent", # estimation method (see optim function).
                    lower = -5, upper = 5,
                    epsilon = NA,# if start is given, the starting point for the estimation is taken in start +- epsilon
                    start = NA, #start value
                    # 1 dimension: (lower,upper) = interval where beta is searched (see optim)
                    # more dimensions: interval for uniform random starting point for beta, in each component
                    deltaName = c("delta"),
                    timeName = c("y"),
                    covNames = c("x"),# covariate names
                    IVNames=c("w"), # Instrumental variables names
                    verbose = FALSE, # print during estimation,
                    maxiter = 5,
                    control = list(trace =1) #see optim
){
  formulaSurv = formula(paste("Surv(",timeName[1],", 1 -", deltaName[1],") ~ 1", sep = "")) # 1- deltaName because we want to KM of Censor variable
  KMsurvfit = survfit(formulaSurv,data = data) #Kaplan Maier estimator for censor variable
  G = KMfunction(KMsurvfit) # Kaplan Maier estimator for censor variable as function
  f = createf(data = data,u=u,deltaName=deltaName,timeName=timeName, covNames=covNames, IVNames=IVNames , G=G)
  if(is.na(start) || is.na(epsilon) ){
    currentStart =  runif(length(covNames),min = lower, max = upper)
  }else{
    currentStart = start + runif(length(covNames),min = -epsilon, max = epsilon)
  }
  #find minimum
  optimRes = optim(par = currentStart,
                   fn = f,method = method #, lower = lower, upper = upper, 
                   #control = control#just for print more intermediate results
  )
  if(min(optimRes$par) > lower & max(optimRes$par) < upper){ #if minim is inside the compact
    return(c(optimRes$par, optimRes$value)) #return minimum
  }
  return(c(rep(NA,length(covNames)), Inf))# return alternative 
}

#################################
#### SIMULATE DATA FUNCTION STUDY
#################################
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
#simulation design 
designs = rbind(
    c("rexp(n)", "0.5*Us + Ws -1 > 0", "runif(n)","rexp(n, 0.176)"), 
    c("rexp(n)", "0.5*Us + Ws -1 > 0", "runif(n)","rexp(n, 0.069)"), 
    c("as.numeric(runif(n))>0.5", "0.5*Us + Ws -1 > 0", "runif(n)","rexp(n, 0.175)"), 
    c("as.numeric(runif(n))>0.5", "0.5*Us + Ws -1 > 0", "runif(n)","rexp(n, 0.07)"), 
    c("exp(rnorm(n))", "0.5*Us + Ws + 0.2*runif(n)", "rexp(n)","rexp(n, 0.065)"), 
    c("exp(rnorm(n))", "0.5*Us + Ws + 0.2*runif(n)", "rexp(n)","rexp(n, 0.0173)"))


inputs = c() #inputs for parallel foor loop 
for(u in us){
  for(i in 1:nrow(designs)){
    for(n in ns){
      for(seed in 1:Nsimulations){
        inputs =rbind(inputs, 
                      c(seed,u,n,designs[i,])
                      )
      }
    }
  }
}

inputs = as.data.frame(inputs); #convert to dataframe
names(inputs) = c("seed", "u", "n", "Wdensity", "X1density", "X2density","Cdensity") #label columns
inputs$seed = as.numeric(inputs$seed);inputs$u = as.numeric(inputs$u);inputs$n = as.numeric(inputs$n); #convert to numeric seed, u, n 


method = "Nelder-Mead" #estimation method 
lower = 0 # see estimate
upper = 1# see estimate
bepsilon = .1# bbeta searched in beta +- bepsilon
covNames = c("x0","x1","x2"); IVNames = c("w","x2");
betau = c(u,  u , u) 

#################### create data

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores-1) #not to overload your computer
registerDoParallel(cl)
allResults <- foreach(i=1:nrow(inputs), .combine=rbind, .packages = packages) %dopar% {
  set.seed(inputs$seed[i])
  n = inputs$n[i]
  #current densities 
  Wdensity = inputs$Wdensity[i]; Cdensity= inputs$Cdensity[i];
  X1density = inputs$X1density[i]; X2density = inputs$X2density[i];
  
  ## generate simulation data 
  simulationData = createSimulationDataFunction(Wdensity = Wdensity,
                                              X1density = X1density,
                                              X2density = X2density,
                                              Cdensity = Cdensity)
  data = simulationData(n)
  percCensor = 1 - sum(data$delta)/n #percentage of censoring
  bs = matrix(NA, nrow = attempts, ncol = length(covNames) + 1) #initialize matrix 
  for(i in 1:attempts){
    bs[i, ] = estimate(data = data,
                       u = u,
                       method = method,
                       covNames = covNames,
                       IVNames = IVNames,
                       lower = lower,
                       upper = upper
                       #, control = list(trace =1, maxit = maxit)
    )
  }
  b = bs[which.min(bs[,ncol(bs)]), 1:length(covNames) ] #select the global minimum 
  bdata = data[sample(1:nrow(data),nrow(data),replace = T),] #boostrap data
  bootbs = matrix(NA, nrow = battempts, ncol = length(covNames) + 1) #initialize matrix 
  for(i in 1:battempts){
    bootbs[i, ] = estimate(data = bdata,
                       u = u,
                       start = b,
                       epsilon = bepsilon, 
                       method = method,
                       covNames = covNames,
                       IVNames = IVNames,
                       lower = lower,
                       upper = upper
                       #, control = list(trace =1, maxit = maxit)
    )
  }
  bboot = bootbs[which.min(bootbs[,ncol(bootbs)]),1:length(covNames)]#select the global minimum 
  configuration = c(u,
                    Wdensity,X1density,X2density,Cdensity, 
                    percCensor,
                    betau,
                    n,
                    attempts,
                    battempts,method)
  return(   rbind(c(FALSE,b,configuration),
                  c(TRUE,bboot,configuration) )
  )
}
stopCluster(cl)

betanames = c("beta0","beta1","beta2")
truebetanames = c("truebeta0","truebeta1","truebeta2")
namesResults = c("isBootstrap",
                 betanames,
                 "u",
                 c("Wdensity", "X1density", "X2density", "Cdensity"),
                 "censor",
                 truebetanames,
                 "n",
                 "attempts",
                 "battempts","method")
allResults = as.data.frame(allResults)
names(allResults) = namesResults
allResults
#save partial results
write.csv(allResults,"raw_simulations.csv", row.names = F)

#### compute statistics
allResults = read.csv("raw_simulations.csv")
Cdensities = unique(allResults$Cdensity) #it will work as identifier of the design
res = c() #initialize result matrix
for(Cdensity in sort(Cdensities)){
  for(n in unique(allResults$n)){
    for(u in unique(allResults$u)){
      print(paste(Cdensity,n))
      results  = allResults[allResults$Cdensity == Cdensity & allResults$n == n & allResults$u == u
                            ,]
      if(nrow(results)){
        current = results[results$isBootstrap == F,]; 
        bcurrent = results[results$isBootstrap == T,];
        ## remove NA values
        noNas = apply(bcurrent, 1, function(x){!any(is.na(x))}) & apply(current, 1, function(x){!any(is.na(x))})
        current = current[noNas ,]
        bcurrent = bcurrent[noNas,]
        n = as.numeric(current[1,c("n")])
        
        biases = colMeans(current[,betanames ] - current[,truebetanames ]) #compute bias
        sds = apply(current[,betanames ] - current[,truebetanames ],2,sd) #compute sd
        bootsds = apply(bcurrent[,betanames ] - current[,betanames ],2,sd)  #compute bootstrap sd
        RMSE = sqrt(mean(apply(current[,betanames ] - current[,truebetanames ],1,function(x){sum(x^2)}))) #compute MSE
        matrixCP = matrix(NA, ncol = length(betanames), nrow =nrow(current)) #initialize coverage probability matrix
        ##compute coverage probability
        for(i in 1:length(betanames)){
          for(j in  1:nrow(current)){
            matrixCP[j,i] = 
              (current[j,c(betanames[i])] - current[j,c(truebetanames[i])] ) > quantile( ( bcurrent[,c(betanames[i])] - current[,c(betanames[i])]),0.025)   &
              (current[j,c(betanames[i])] - current[j,c(truebetanames[i])] ) < quantile( ( bcurrent[,c(betanames[i])] - current[,c(betanames[i])]),0.975)
          }
        }
      
        res = rbind(res,
                    c(nrow(current),current[1,c("n")],current[1,c("u")], current[1,c("Wdensity")], current[1,c("X1density")], current[1,c("X2density")], Cdensity, mean(current$censor), 
                      biases,RMSE,
                      sds,
                      bootsds,
                      colSums(matrixCP)/nrow(current) )   # merge results
        )
      }
      
    }
  }
}
## covert to dataframe
res = as.data.frame(res)
names(res) = c("nsimulation",
               "n",
               "u", "Wdensity", "X1density", "X2density", "Cdensity", "censor",
               sapply(truebetanames, function(x){paste("bias_", x, sep = "")}),
               "RMSE",
               sapply(truebetanames, function(x){paste("sd_", x, sep = "")}),
               sapply(truebetanames, function(x){paste("boot_sd_", x, sep = "")}),
               sapply(truebetanames, function(x){paste("CP95_", x, sep = "")})
)
#construct columns referring to the proposed estimator of table 2 of the paper
table2 = c()
for(Wdensity in c("rexp(n)","exp(rnorm(n))","as.numeric(runif(n))>0.5")){
  current = res[res$Wdensity == Wdensity,]
  current = as.data.frame(current)
  table2 = rbind(table2, current[order(current$u,as.numeric(current$n),current$censor),])
}
table2 = table2[,c("Wdensity","u","n","censor",
          sapply(truebetanames, function(x){paste("bias_", x, sep = "")}),
          "RMSE",
          sapply(truebetanames, function(x){paste("CP95_", x, sep = "")})
          )]
table2 <- cbind(table2$Wdensity,as.data.frame(lapply(table2[,2:ncol(table2)], as.numeric)))
options("digits" = 3)
table2
write.csv(table2,"simulations.csv",row.names = FALSE)
