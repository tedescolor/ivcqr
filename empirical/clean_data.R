path = getwd()#/Data subfolder
setwd(path)
#LOAD PACKAGES
{
  # Package names
  packages <-
    c("survival",
      "quantreg",
      "foreach",
      "foreign",
      "haven",
      "doParallel")
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages],
                     repos='http://cran.us.r-project.org')
  }
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
}

# spells 
spells = read_dta("jtpa/Data/expmon/expmon_stata/empevr.dta")
# filter time data: only ID and evempl-- variables
spells = cbind(spells$recid, 
              spells[,names(spells)[sapply(names(spells),function(cov_name){startsWith(cov_name, "evempl")}) ]])

# only people starting unemployed (evempl00 == 0)
spells = spells[spells$evempl00 == 0, ]
#initialize array for event-time and deltas
deltas = array(0, dim = nrow(spells)) 
events = array(0, dim = nrow(spells))

for(i in 1:nrow(spells)){ #for every subject
  row = spells[i,-c(1)] #take the spell 
  if(sum(row!=0)> 0){ #check if something occurs
    t = which(row != 0)[1] # check when it occurs
    events[i] = t #save 
    if(row[t] == 1){ #check if it is uncensored
      deltas[i] = 1;
    }
  }else{
    events[i] = length(row) #nothing occurred so censored time at the end of the spell
  }
}

#create data wih id, time, deltas
events = as.numeric(events)
data = as.data.frame(cbind(spells[,c(1)],events,as.numeric(deltas)))
names(data) = c("recid", "time","delta")

# jtpa programm variable
jtpa = as.data.frame(read_dta(file = "jtpa/Data/expfup1/expfup1_stata/main.dta"))
names(jtpa)
jtpa = jtpa[,c("recid","f1d7")]
table(jtpa$f1d7)
jtpa = jtpa[jtpa$f1d7 %in% c(1,2),] # only known variable, 2 = not partecipating

# covariates
covariates = read_dta("jtpa/Data/expbif/expbif.dta")
names(covariates)
# idname + covariate to save
cov_names = c("recid",# id
              "white", # 1= white, 0=not-white
              "sex", #1 = male, 2=female
              "age", 
              "yearearn",#year salary
              "ra_stat", #1 treatment, 2 control
              "hsged", #high school completed = 1, unknown = 9
              "bfmarsta" # married status, 2 = true, 1=single, then other, 9= unknown
              ) 
children = apply(covariates[,c("chunder4","ch4or5", "ch6to12","chover18","ch13to18")],1,function(x){
  count = sum(as.numeric(x)) #sum of the children indicators 
  return(ifelse(count>=9,9,#if there is a 9 in the row ==> count>=9 ==> return 9 (unknown),
                ifelse(count==0,0,1)))  #otherwise check if there count is 0 (no children)
})
covariates$children = children
cov_names = c(cov_names, "children")
#filter the covariate
covariates = as.data.frame(covariates[,names(covariates) %in% cov_names])
head(covariates)

# add covariate to data
data = merge(data,covariates)
# add jtpa partecipation
data = merge(data,jtpa)
#transform to numeric all the data
data[] <- lapply(data, function(x) as.numeric(as.character(x)))


head(data)
#only complete rows
data = data[complete.cases(data),]

#remove rows undefined values
data = data[ data$children != 9  &
               data$bfmarsta != 9 &
              data$yearearn != 99999 &
               data$hsged != 9
            ,]
data = data[complete.cases(data),]
head(data)
#pecentage censoring
1-sum(as.numeric(data$delta))/nrow(data)

#transform variable to Boolean
data$jtpa = as.numeric(data$f1d7 == 1)
data$treatment = as.numeric(data$ra_stat == 1)
data$single = as.numeric(data$bfmarsta==1)
data$hsged = as.numeric(data$hsged == 1)
data$male = as.numeric(data$sex == 1)
data = data[,-which(names(data) %in% c("bfmarsta","ra_stat","f1d7", "sex"))]
head(data)

headf1<-read_dta("jtpa/Data/expfup1/expfup1_stata/head_fp1.dta")
headf2<-read_dta("jtpa/Data/expfup2/expfup2_stata/head_fp2.dta")

jobf1<-read_dta("jtpa/Data/expfup1/expfup1_stata/bjob.dta")
jobf2<-read_dta("jtpa/Data/expfup2/expfup2_stata/bjob.dta")


beginendfu1<-cbind.data.frame(as.numeric(headf1$recid),headf1$f1beg,headf1$f1ref) # beginning of study and first follow up date
colnames(beginendfu1)<-c("recid","beg","end1")
datejobfu1<-cbind.data.frame(jobf1$recid,jobf1$f1b2bdd,jobf1$f1b2bmm,jobf1$f1b2byy) #check for 92s date of job and multiple dates!!!!!!

endfu2<-cbind.data.frame(as.numeric(headf2$recid), headf2$fdat_fin) # day of second follow upinterview
colnames(endfu2)<-c("recid","end2")
datejobfu2<-cbind.data.frame(jobf2$recid,jobf2$f2b2bdd,jobf2$f2b2bmm,jobf2$f2b2byy) #check for 92s date of job and multiple dates!!!!!!

datejobfu1[datejobfu1$`jobf1$f1b2bdd` == 92,c("jobf1$f1b2bdd")] = "01"# 92 == start of the month
datejobfu1[datejobfu1$`jobf1$f1b2bdd` == 93,c("jobf1$f1b2bdd")] = "14"# 93 == middle of the month
datejobfu1[datejobfu1$`jobf1$f1b2bdd` == 94,c("jobf1$f1b2bdd")] = "28"# 93 == end of the month

datejobfu1 = datejobfu1[!(datejobfu1$`jobf1$f1b2bdd` %in% c(98,99)),]#remove 98 = dont known, 99 =missing
datejobfu1 = datejobfu1[(datejobfu1$`jobf1$f1b2byy` != 99),] #99  unknown


datejobfu1 = datejobfu1[datejobfu1$`jobf1$recid` %in% data$recid, ] #only the sibjects in our data 
head(datejobfu1)



datejobfu1$date = sapply(1:nrow(datejobfu1),function(i){
  year = datejobfu1$`jobf1$f1b2byy`[i]
  month = datejobfu1$`jobf1$f1b2bmm`[i]
  day = datejobfu1$`jobf1$f1b2bdd`[i]
  return(paste("19",year,"-",month,"-",day,sep = ""))
})

datejobfu1_clean = c()
done_recid = c()
for(i in 1:nrow(datejobfu1)){
  id = datejobfu1$`jobf1$recid`[i]
  if(!(id%in%done_recid)){
    done_recid = c(done_recid, id)
    current = datejobfu1[datejobfu1$`jobf1$recid` == id, c("date")]
    datejobfu1_clean = rbind(datejobfu1_clean, c(id,min(current)))
  }

}
datejobfu1_clean = as.data.frame(datejobfu1_clean)
names(datejobfu1_clean) = c("recid", "date")
head(datejobfu1_clean)

data = merge(x = data, y = datejobfu1_clean, by = "recid", all.x = TRUE)
sum(is.na(data$date))/nrow(data)
sum(datejobfu1$`jobf1$recid` %in% beginendfu1$recid)/length(unique(beginendfu1$recid))



### second follow up

datejobfu2[datejobfu2$`jobf2$f2b2bdd` == 92,c("jobf2$f2b2bdd")] = "01"# 92 == start of the month
datejobfu2[datejobfu2$`jobf2$f2b2bdd` == 93,c("jobf2$f2b2bdd")] = "14"# 93 == middle of the month
datejobfu2[datejobfu2$`jobf2$f2b2bdd` == 94,c("jobf2$f2b2bdd")] = "28"# 93 == end of the month

datejobfu2 = datejobfu2[!(datejobfu2$`jobf2$f2b2bdd` %in% c(98,99)),]#remove 98 = dont known, 99 =missing
datejobfu2 = datejobfu2[(datejobfu2$`jobf2$f2b2byy` != 99),] #99  unknown


#only the sibjects in our data that have NA as date
done_recid = c()
datejobfu2_cleaned = c()
for(i in 1:nrow(datejobfu2)){
  id = datejobfu2$`jobf2$recid`[i]
  if( !(id %in% done_recid) ){
    if( (id %in% data$recid) && is.na(data[data$recid == id, c("date")])){
      datejobfu2_cleaned = rbind(datejobfu2_cleaned, datejobfu2[i,])
    }
  }
}
datejobfu2_cleaned = as.data.frame(datejobfu2_cleaned)
names(datejobfu2_cleaned) = names(datejobfu2)


datejobfu2_cleaned$date = sapply(1:nrow(datejobfu2_cleaned),function(i){
  year = datejobfu2_cleaned$`jobf2$f2b2byy`[i]
  month = datejobfu2_cleaned$`jobf2$f2b2bmm`[i]
  day = datejobfu2_cleaned$`jobf2$f2b2bdd`[i]
  return(paste("19",year,"-",month,"-",day,sep = ""))
})


datejobfu2_cleaned2 = c()
done_recid = c()
for(i in 1:nrow(datejobfu2_cleaned)){
  id = datejobfu2_cleaned$`jobf2$recid`[i]
  if(!(id%in%done_recid)){
    done_recid = c(done_recid, id)
    current = datejobfu2_cleaned[datejobfu2_cleaned$`jobf2$recid` == id, c("date")]
    datejobfu2_cleaned2 = rbind(datejobfu2_cleaned2, c(id,min(current)))
  }
  
}
datejobfu2_cleaned2 = as.data.frame(datejobfu2_cleaned2)
names(datejobfu2_cleaned2) = c("recid", "date")
head(datejobfu2_cleaned2)

count = 0
for(i in 1:nrow(data)){
  id = data$recid[i]
  if(is.na(data$date[i])){
    if(id %in% datejobfu2_cleaned2$recid){
      count = count +1
      data[i,c("date")] = datejobfu2_cleaned2[datejobfu2_cleaned2$recid == id, c("date")]
    }
  }
}
sum(is.na(data$date))/nrow(data)

data$has_date = !is.na(data$date)

table(data$has_date, data$delta)
#remove inconstency with emprf
data = data[!(data$has_date == TRUE & data$delta == 0), ]
data = data[!(data$has_date == FALSE & data$delta == 1), ]

months = as.data.frame(rbind(c("JAN","01"),
                             c("FEB","02"),
                             c("MAR","03"),
                             c("APR","04"),
                             c("MAY","05"),
                             c("JUN","06"),
                             c("JUL","07"),
                             c("AUG","08"),
                             c("SEP","09"),
                             c("OCT","10"),
                             c("NOV","11"),
                             c("DEC","12")))
names(months) = c("name", "number")
endfu2$date = sapply(1:nrow(endfu2),function(i){
  str =  endfu2$end2[i]
  if(str == ""){
    return(NA)
  }
  day = substr(str,1,2)
  month = months[months$name == substr(str, 3,5), c("number") ]
  year = substr(str, 6,7)
  return(paste("19",year,"-",month,"-",day,sep = ""))})

  # assign censoring date, for follow up 2
for(i in 1:nrow(data)){
  if(data$delta[i] == 0){
    id = data$recid[i]
    if(id %in% endfu2$recid){
      data[i,c("date")] = endfu2$date[endfu2$recid == id]
    }
  }
}

# assign censoring date, for follow up 1

beginendfu1$date_end = sapply(1:nrow(beginendfu1),function(i){
  str = beginendfu1$end1[i]
  return(paste("19",substr(str,1,2),"-",substr(str,3,4),"-",substr(str,5,6),sep = ""))
})

beginendfu1$date_beg = sapply(1:nrow(beginendfu1),function(i){
  str = beginendfu1$beg[i]
  return(paste("19",substr(str,1,2),"-",substr(str,3,4),"-",substr(str,5,6),sep = ""))
})

sum(is.na(data))
for(i in 1:nrow(data)){
  if(is.na(data$date[i])){
    id = data$recid[i]
    if(id %in% beginendfu1$recid){
      data[i,c("date")] = beginendfu1$date_end[beginendfu1$recid == id]
    }
  }
}
data$date_start = NA
for(i in 1:nrow(data)){
  id = data$recid[i]
  data[i,c("date_start")] = beginendfu1[beginendfu1$recid == id, c("date_beg")]
}

data$days =  sapply(1:nrow(data), function(i){
  return(as.integer(difftime(data$date[i],data$date_start[i], units = "days")))
})
data = data[data$days>0,]
head(data)
par(mfrow=c(2,1))
hist(data$time)
hist(data$day)

par(mfrow=c(2,1))

hist(data[data$delta == 1, c("days")])
hist(data[data$delta == 0, c("days")])
data = data[,-c(2,13,14,15)]
nrow(data)

write.csv(data, "data.csv",row.names = FALSE)

