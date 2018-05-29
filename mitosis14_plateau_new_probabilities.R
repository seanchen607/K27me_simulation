rm(list=ls())
options(scipen=999)
chromlength=1000
prc2=1
populationSize=10
life=20
Tp01=0.9
Tp02=0
Tp03=0
Tp12=0.5
Tp13=0
Tp23=0.2
Tp00=0
Tp11=0
Tp22=0

mitop=0.1
depop=0.9
newdir=paste("~/Desktop/Histone_Mark_Simulation/Steadystate/",gsub(":","_",gsub(" ","-",date())),"-","chrmlngth_",chromlength,"-prc2_",prc2,"-pop_",populationSize,"-add01-12-23_",Tp01,"-",Tp12,"-",Tp23,"-mitop_",mitop,"/",sep = "")
system(paste("mkdir -p ",newdir,sep = ""))
system(paste("mkdir -p ",newdir,"cell1/",sep = ""))
setwd(newdir)

#### Defining normal distribution with min and max and mean and sd
mysamp <- function(n, m, s, lwr, upr, nnorm) {
  samp <- rnorm(nnorm, m, s)
  samp <- samp[samp >= lwr & samp <= upr]
  if (length(samp) >= n) {
    return(sample(samp, n))
  }  
  stop(simpleError("Not enough values to sample from. Try increasing nnorm."))
}

#### Defining Chromatin 0 (with chromSize nuecleosomes all me0)
createNakedChromatin <- function(chromSize)
{
  chromatin<-data.frame(integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize))
  colnames(chromatin)<-c("N","S","me0","me1","me2","me3")
  for (s in 1:chromSize)
  {
    chromatin$N[s]=s
    chromatin$S[s]=0
    chromatin$me0[s]=1
    chromatin$me1[s]=0
    chromatin$me2[s]=0
    chromatin$me3[s]=0
  }
  results<-list("add1"=0,"add2"=0,"add3"=0,"add0"=0,"rando"=0,"chr"=chromatin)
  return(results)
}

chromatin<-createNakedChromatin(chromlength)


#### Defining the probability of PRC2 falling off chromatin based on nucleosomeNumber and length of the chromosome
prc2hl <- function(nucleosomeNumber,chromSize)
{
  # PRC2 half-life
   if (1-(nucleosomeNumber/chromSize)-(2/10) < 0){lower=0} else {lower=1-(nucleosomeNumber/chromSize)-(2/10)}
   if (1-(nucleosomeNumber/chromSize)+(2/10) > 1){upper=1} else {upper=1-(nucleosomeNumber/chromSize)+(2/10)}
   meaner=(upper+lower)/2
   hl=mysamp(n=1, m=meaner, s=0.1, lwr=lower, upr=upper, nnorm=100)
  return(hl)  
}

deposit <- function(chr,chromSize)
{
  for (nucleosomeNumber in 1:chromSize)
  {
    falloff=prc2hl(nucleosomeNumber,chromSize)
    rando=runif(n=1,min=0,max=100)/100
    rando=rando
    if (chr$S[nucleosomeNumber] == 0)
    {
      if (rando <= Tp00*falloff)
      {
        chr$S[nucleosomeNumber]=0
        chr$me0[nucleosomeNumber]=1
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=0
      } else if(rando < Tp01*falloff) {
        chr$S[nucleosomeNumber]=1
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=1
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=0
      } else  if(rando < Tp02*falloff ) {
        chr$S[nucleosomeNumber]=2
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=1
        chr$me3[nucleosomeNumber]=0
      } else  if(rando < Tp03*falloff ) {
        chr$S[nucleosomeNumber]=3
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=1
      }
    } else if (chr$S[nucleosomeNumber] == 1)  {
      if(rando < Tp11*falloff)
      {
        chr$S[nucleosomeNumber]=1
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=1
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=0
      } else if(rando < Tp12*falloff) {
        chr$S[nucleosomeNumber]=2
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=1
        chr$me3[nucleosomeNumber]=0
      } else  if(rando < Tp13*falloff) {
        chr$S[nucleosomeNumber]=3
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=1
      } 
    } else if (chr$S[nucleosomeNumber] == 2)  {
      if(rando < Tp22*falloff)
      {
        chr$S[nucleosomeNumber]=2
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=1
        chr$me3[nucleosomeNumber]=0
      } else if(rando < Tp23*falloff)  {
        chr$S[nucleosomeNumber]=3
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=1
      }
    }
  } 


  results<-list("Tp01"=Tp01,"Tp12"=Tp12,"Tp23"=Tp12,"rando"=rando,"chr"=chr)
  return(results)
}


#### Defining mitosis function
mitosis <- function (chromatin){
  d <- chromatin
  for (i in 1:nrow(chromatin)) {
    check=runif(n=1,min = 0,max=1)
    if (check > 0.5)
    {
      d$S[i]=0
      d$me0[i]=1
      d$me1[i]=0
      d$me2[i]=0
      d$me3[i]=0
    }
  }
  return (d)
}

# creating the initial population 
population<-list()
currentpop<-list()  
currentpop<-list("lifecycle"=0,"chr"=createNakedChromatin(chromlength)[["chr"]])  
currentpop[["chr"]]$me0=0
mitonum=0
lifecycle=1
me0avg<-c()
me1avg<-c()
me2avg<-c()
me3avg<-c()

while(lifecycle <= life)
{

  temp<-data.frame(integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength))
  colnames(temp)<-c("N","S","me0","me1","me2","me3")
  for (s in 1:chromlength)
  {
    temp$N[s]=s
    temp$S[s]=0
    temp$me0[s]=0
    temp$me1[s]=0
    temp$me2[s]=0
    temp$me3[s]=0
  }
  
  
  print(lifecycle)
  for (p in 1:populationSize)
  {
#  print (paste("individual: ",p,sep=""))
  if (lifecycle==1){population[[p]]<-createNakedChromatin(chromlength)}  
  for (round in 1:prc2)
    {
    if (runif(n=1,min=0,max=1) >= 1-depop)
      {
      population[[p]]<-deposit(population[[p]][["chr"]],chromlength)
      }
    }
    if (runif(n=1,min=0,max=1) >= 1-mitop)
      {
      population[[p]][["chr"]]<-mitosis(population[[p]][["chr"]])
      mitonum=mitonum+1
      }
    # if(populationSize>1){
      temp<-temp+population[[p]][["chr"]]
      # }
    # else(temp<-population[[p]][["chr"]])
  }
  
# if (lifecycle==1)
#   {
currentpop=list("lifecycle"=c(currentpop[["lifecycle"]],lifecycle),"chr"=temp)
# } else {
# currentpop=list("lifecycle"=c(currentpop[["lifecycle"]],lifecycle),"chr"=temp+currentpop[["chr"]])
# }
  
 png(paste("cell1/",lifecycle,".ind.","1","-cycles",lifecycle,"size",chromlength,".png",sep=""),width = 800,height = 600)
 par(mfrow=c(4,1))
 plot(population[[1]][["chr"]]$me0,type="h",ylim=c(0,1),ylab = "K36me0", main = paste("Cycle: ",lifecycle," Individual: 1",sep = ""),xlab = "Chromatin")
 plot(population[[1]][["chr"]]$me1,type="h",ylim=c(0,1),ylab = "K36me1",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$me2,type="h",ylim=c(0,1),ylab = "K36me2",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$me3,type="h",ylim=c(0,1),ylab = "K36me3",xlab = "Chromatin")
 dev.off()
# 

png(paste(populationSize,"_pop.","-cycles_",lifecycle,"-chromsize_",chromlength,".png",sep=""),width = 800,height = 600)
par(mfrow=c(4,1))
plot((currentpop[["chr"]]$me0)/(populationSize),type="h",ylab = "K36me0", main = paste("Cycle: ",lifecycle,sep = ""),xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$me1)/(populationSize),type="h",ylab = "K36me1", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$me2)/(populationSize),type="h",ylab = "K36me2", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$me3)/(populationSize),type="h",ylab = "K36me3", xlab = "Chromatin",ylim=c(0,1))

dev.off()

me0avg=c(me0avg,mean((currentpop[["chr"]]$me0)/populationSize))
me1avg=c(me1avg,mean((currentpop[["chr"]]$me1)/populationSize))
me2avg=c(me2avg,mean((currentpop[["chr"]]$me2)/populationSize))
me3avg=c(me3avg,mean((currentpop[["chr"]]$me3)/populationSize))

if (lifecycle > 10){
var0=var(me0avg[lifecycle:(lifecycle-10)])
var1=var(me1avg[lifecycle:(lifecycle-10)])
var2=var(me2avg[lifecycle:(lifecycle-10)])
var3=var(me3avg[lifecycle:(lifecycle-10)])

if (var0 < 0.0001 && var1 < 0.0001 && var2 < 0.0001 && var3 < 0.0001 )
{
  lastlap=lifecycle
  lifecycle=life
}
}
lastlap=lifecycle

lifecycle=lifecycle+1
currentpop[["chr"]]$N=population[[p]][["chr"]]$N
}
print(paste(mitonum," mitosis events happened during this simulation.",sep=""))
png(paste("4plots-",populationSize,"_pop.","-cycles_",lastlap,"-chromsize_",chromlength,".png",sep=""),width = 800,height = 600)
par(mfrow=c(2,2))
plot(me0avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K36me0",xlab = "Cycle")
text(x = (lastlap/10)+5,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+5),y = 0.85,labels = paste("Average mark: ",me0avg[lastlap],sep = ""),pos = 1,offset = 0)
plot(me1avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K36me1",xlab = "Cycle")
text(x = (lastlap/10)+5,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+5),y = 0.85,labels = paste("Average mark: ",me1avg[lastlap],sep = ""),pos = 1,offset = 0)
plot(me2avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K36me2",xlab = "Cycle")
text(x = (lastlap/10)+5,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+5),y = 0.85,labels = paste("Average mark: ",me2avg[lastlap],sep = ""),pos = 1,offset = 0)
plot(me3avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K36me3",xlab = "Cycle")
text(x = (lastlap/10)+5,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+5),y = 0.85,labels = paste("Average mark: ",me3avg[lastlap],sep = ""),pos = 1,offset = 0)
dev.off()
