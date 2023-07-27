path = getwd() #set path to the folder empirical/
setwd(path)
test = FALSE
#### test setting
if(test){
  attempts = 5
  N = 2 #boostrap resampling
  us = c(0.1)
}else{
  ### paper setting
  attempts = 100
  N = 250 #boostrap resampling
  us = c(0.1,0.2,0.3,0.4,0.5,0.6)
}

{
  # Package names
  packages <-
    c("survival",
      "foreach",
      "doParallel")
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
    index = which.max( # index of the max
      KMsurvfit$time[ # of the times
        KMsurvfit$time -rep(t,length(KMsurvfit$time)) <= 0  # among the times which are lower or equal then t
      ]
    )
    return(KMsurvfit$surv[ # return the $surv value
      ifelse(index == 0, 1, index) #of the index, if not 0, otherwise 1
    ])
  })
}
#given the data, returns the function data has to be minimized
createf = function(data,u,deltaName,timeName, covNames, IVNames, G, lower, upper){ return( 
  function(beta){
    #print(beta)
    n = nrow(data)
    res = 0
    # OSS: weights can be computed only once:
    weigths = sapply(data[,timeName],FUN = G)
    # OSS: sapply(data$y,FUN = G) can return zero.
    #So we set them to NA and when we will sum (see for loop), we set na.rm = TRUE
    weigths[weigths==0] <- NA
    deltas = data[,deltaName] #extract delta column
    deltas[deltas == 0] = NA # sub with NA, so that we remove in the next sum that terms by na.rm = TRUE
    weigths = deltas / weigths
    for(j in 1:n){
      #for optimization, we compute the indicator function for wi<=wj just once, and save it in the following wsj variable:
      wsj = as.numeric(apply( data.matrix(data[,IVNames]) <= data.matrix(data[rep(j, nrow(data)),IVNames]),1,all))
      # compute the statistic for the value wj
      res = res + ( (sum(
        data[,timeName]  <= exp(matrix(beta, nrow = 1) %*% t(data.matrix(data[,covNames]))) *
          wsj * weigths, na.rm = TRUE) - u * sum(wsj) )/n )^2 
    }
    return(res/n)  # it should be res/n, but might be it is < than computer tolerance, and not proper minimization occurs.
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
                    maxiter = 1,
                    control = list(trace =1) #see optim
){
  formulaSurv = formula(paste("Surv(",timeName[1],", 1 -", deltaName[1],") ~ 1", sep = ""))
  KMsurvfit = survfit(formulaSurv,data = data) #Kaplan Maier estimator for censor variable
  G = KMfunction(KMsurvfit) # Kaplan Maier estimator for censor variable as function
  # to avoid local minima, we minimize the function multiple times (attempts)
  # and return the value which obtains the minimum
  f = createf(data = data,u=u,deltaName=deltaName,timeName=timeName, covNames=covNames, IVNames=IVNames , G=G,lower = lower, upper = upper)
  iter = 0
  while(iter<maxiter){
      iter = iter+1; 
      if(verbose){print(paste(i, " - ", iter, sep = ""))}
      if(is.na(start) || is.na(epsilon) ){
        currentStart = runif(length(covNames), min= lower, max = upper)
      }else{
        currentStart = start + runif(length(covNames),min = -epsilon, max = epsilon)
      }
      optimRes = optim(par = currentStart,
                       fn  = f,
                       method = method,
                       #control = control#just for print more intermediate results
      )
      if(min(optimRes$par) > lower & max(optimRes$par) < upper){
        return(c(optimRes$par, f(optimRes$par)))
      } 
    }
    return(c(rep(NA, length(covNames)), Inf))

}
df = read.csv("data.csv")
head(df)
df = df[complete.cases(df),]  #remove rows with missing data
print(paste(nrow(df), "rows"))
df = df[ df$white == 0 & df$male == 0 & df$single==1 & df$children==1,] # interesting subset
bepsilon = .5;# bootstrap beta search in estimated beta + - bepsilon
estimations = read.csv("estimation.csv")
method = "Nelder-Mead";# method minimization 
lower = -5; upper = 10; #beta search in [lower, upper] interval (each component)
df$intercept = 1;
deltaName = c("delta");
timeName = c("days");#df$intxrel.days = 30* df$intxrel; timeName = c("intxrel.days") 
covNames = c("intercept","treatment","age")
IVNames = c("jtpa","age")
final_results = c()
for(u in us){
  print(u)
  for(k in 1:N){
      print(k)
      set.seed(k)
      bdf = df[sample(1:nrow(df), nrow(df), replace=T), ]
      start = estimations[estimations$u == u,covNames ]
      #setup parallel backend to use many processors
      cores=detectCores()
      cl <- makeCluster(cores[1]) #not to overload your computer
      registerDoParallel(cl)
      results <- foreach(i=1:attempts, .combine=rbind, .packages = packages) %dopar% {
        set.seed(i)
        b = 
          estimate(data = bdf,
                       u = u,
                       start = start + runif(3,min=-bepsilon,max = bepsilon),
                       epsilon = 0,
                       method = method,
                       lower = lower,
                       upper = upper,
                       deltaName = deltaName,covNames = covNames, timeName = timeName,IVNames = IVNames,
                        verbose = F)
        return(c(T,b)
        )
        
      }
      #stop cluster
      stopCluster(cl)
      results = as.data.frame(results)
      names(results) = c("isBootstrap", covNames, "f")
      results$u = u; results$method = method; results$lower = lower; results$upper = upper;
      results = results[which.min(results$f),]
      final_results = rbind(final_results, results)
  }
}
names(final_results) = names(results)
final_results
#cnstruct confidence intervals
CIs = c()
for(u in us){
  current = final_results[final_results$u==u,covNames] 
  #current = current - estimations[estimations$u == u,covNames ][rep(1,nrow(current)),]
  CIs = rbind(CIs, cbind(u, as.numeric(estimations[estimations$u == u,covNames ]),
                         apply(current,2,quantile,probs = 0.025),
                         apply(current,2,quantile,probs = 0.975)))
}
CIs = as.data.frame(CIs)
names(CIs) = c("u","Estimation", "CI025","CI975")
CIs
write.csv(CIs,"bootstrap_results.csv", row.names = F)


