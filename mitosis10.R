rm(list=ls())
newdir=paste("~/Desktop/Histone_Mark_Simulation/",gsub(":","_",gsub(" ","-",date())),"/",sep = "")
system(paste("mkdir -p ",newdir,sep = ""))
system(paste("mkdir -p ",newdir,"cell1/",sep = ""))
setwd(newdir)
options(scipen=999)
chromlength=1000
prc2=5
populationSize=100
burden12=0.3
burden23=0.6
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
return(chromatin)
}

chromatin<-createNakedChromatin(chromlength)


#### Defining the probability of PRC2 falling off chromatin based on nucleosomeNumber and length of the chromosome
prc2hl <- function(nucleosomeNumber,chromSize)
{
  # PRC2 half-life
  if (0.3-(nucleosomeNumber/chromSize) < 0){lower=0} else {lower=0.3-(nucleosomeNumber/chromSize)}
  if (1.3-(nucleosomeNumber/chromSize) > 1){upper=1} else {upper=1.3-(nucleosomeNumber/chromSize)}
  if (0.8-(nucleosomeNumber/chromSize) < 0){meaner=0} else {meaner=0.8-(nucleosomeNumber/chromSize)}
  hl=mysamp(n=1, m=meaner, s=0.05, lwr=lower, upr=upper, nnorm=100)
  return(hl)  
}

deposit <- function(chr,chromSize)
{
  for (nucleosomeNumber in 1:chromSize)
  {
    falloff=prc2hl(nucleosomeNumber,chromSize)
   # print(falloff)
    add1=mysamp(n=1, m=0.90, s=0.01, lwr=0, upr=1, nnorm=100)*falloff
    add2=mysamp(n=1, m=0.01, s=0.001, lwr=0, upr=0.1, nnorm=100)*falloff
    add3=mysamp(n=1, m=0.00, s=0.0, lwr=0, upr=0, nnorm=100)*falloff
    add0=1-(sum(add1+add2+add3))
    rando=runif(n=1,min=0,max=100)/100
    
    if (chr$S[nucleosomeNumber] == 0)
    {
       if (rando >= (add1+add2+add3))
      {
        chr$S[nucleosomeNumber]=0
        chr$me0[nucleosomeNumber]=1
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=0
      } else if(rando < add1 && rando > 0) {
        chr$S[nucleosomeNumber]=1
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=1
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=0
      } else  if(rando < (add2+add1) && rando > add1 ) {
        chr$S[nucleosomeNumber]=2
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=1
        chr$me3[nucleosomeNumber]=0
      } else  if(rando < (add3+add2+add1) && rando > (add1+add2) ) {
        chr$S[nucleosomeNumber]=3
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=1
      }
    } else if (chr$S[nucleosomeNumber] == 1)  {
      if(rando >= add1+add2+add3+burden12)
      {
        chr$S[nucleosomeNumber]=1
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=1
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=0
      } else if(rando < add1 && rando > 0) {
        chr$S[nucleosomeNumber]=2
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=1
        chr$me3[nucleosomeNumber]=0
      } else  if(rando < (add2+add1) && rando > add1) {
        chr$S[nucleosomeNumber]=3
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=1
      } 
    } else if (chr$S[nucleosomeNumber] == 2)  {
      if(rando > add1+add2+add3+burden23)
      {
        chr$S[nucleosomeNumber]=2
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=1
        chr$me3[nucleosomeNumber]=0
      } else if(rando < add1 && rando > 0)  {
        chr$S[nucleosomeNumber]=3
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=1
      }
    }
  } 
  results<-list("add1"=add1,"add2"=add2,"add3"=add3,"add0"=add0,"rando"=rando,"chr"=chr)
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

depositedchr<-createNakedChromatin(chromlength)


for (i in 1:prc2){
png(paste("cell1/",i,".cell1.png",sep=""),width = 800,height = 600)

depositedchr<-deposit(depositedchr,chromlength)[["chr"]]
mitotic<-mitosis(depositedchr)

par(mfrow=c(4,2))
plot(depositedchr$me0,type="h",main=paste("Cell-1 K36me0 : PRC2 ",i,sep=""),ylim=c(0,1))
plot(mitotic$me0,type="h",main=paste("Mitotic Cell-1 K36me0 : PRC2 ",i,sep=""),ylim=c(0,1))
plot(depositedchr$me1,type="h",main=paste("Cell-1 K36me1 : PRC2 ",i,sep=""),ylim=c(0,1))
plot(mitotic$me1,type="h",main=paste("Mitotic Cell-1 K36me1 : PRC2 ",i,sep=""),ylim=c(0,1))
plot(depositedchr$me2,type="h",main=paste("Cell-1 K36me2 : PRC2 ",i,sep=""),ylim=c(0,1))
plot(mitotic$me2,type="h",main=paste("Mitotic Cell-1 K36me0 : PRC2 ",i,sep=""),ylim=c(0,1))
plot(depositedchr$me3,type="h",main=paste("Cell-1 K36me3 : PRC2 ",i,sep=""),ylim=c(0,1))
plot(mitotic$me3,type="h",main=paste("Mitotic Cell-1 K36me0 : PRC2 ",i,sep=""),ylim=c(0,1))
dev.off()

}

# creating the initial population 
# the number of columns depends on the number of prc2 passages
colnumber=(prc2-1)*4
population<- data.frame(matrix(nrow = chromlength,ncol = (4*prc2)))
popcolnames<-c()
for (i in 1:prc2){
  popcolnames<-c(popcolnames,paste(i,"me0",sep="_"))
  popcolnames<-c(popcolnames,paste(i,"me1",sep="_"))
  popcolnames<-c(popcolnames,paste(i,"me2",sep="_"))
  popcolnames<-c(popcolnames,paste(i,"me3",sep="_"))
}
colnames(population)<-popcolnames
for (n in 1:chromlength)
{
  population[n,]=0
}


### Writing marks on population
for (cell in 1:populationSize){
depositedchr<-createNakedChromatin(chromlength)
depositedchr<-deposit(chromatin,chromlength)[["chr"]]
for (i in 1:prc2)
  {
  print (paste("Original - ",cell,":",i,sep=""))
  population[,(i*4)-3]=population[,(i*4)-3]+depositedchr$me0
  population[,(i*4)-2]=population[,(i*4)-2]+depositedchr$me1
  population[,(i*4)-1]=population[,(i*4)-1]+depositedchr$me2
  population[,(i*4)]=population[,(i*4)]+depositedchr$me3
  depositedchr<-deposit(depositedchr,chromlength)[["chr"]]
  }
}


# for (i in 1:prc2)
# {
#   png(paste(i,".mitosis.png",sep=""),width = 800,height = 600)
#   par(mfrow=c(3,1))
#   # plot((population[,(i*4)-3])/populationSize,type="h",ylim=c(0,1),main=paste("K36me0 prc2 round: ",i,sep=""))
#   plot((population[,(i*4)-2])/populationSize,type="h",ylim=c(0,1),main=paste("K36me1 prc2 round: ",i,sep=""))
#   plot((population[,(i*4)-1])/populationSize,type="h",ylim=c(0,1),main=paste("K36me2 prc2 round: ",i,sep=""))
#   plot((population[,(i*4)-0])/populationSize,type="h",ylim=c(0,1),main=paste("K36me3 prc2 round: ",i,sep=""))
#   dev.off()
#   
# }



##################################################################
##################################################################
######################      MITOSIS    ###########################
##################################################################
##################################################################

colnumber=(prc2-1)*4
mpopulation<- data.frame(matrix(nrow = chromlength,ncol = (4*prc2)))
popcolnames<-c()
for (i in 1:prc2){
  popcolnames<-c(popcolnames,paste(i,"me0",sep="_"))
  popcolnames<-c(popcolnames,paste(i,"me1",sep="_"))
  popcolnames<-c(popcolnames,paste(i,"me2",sep="_"))
  popcolnames<-c(popcolnames,paste(i,"me3",sep="_"))
}
colnames(mpopulation)<-popcolnames
for (n in 1:chromlength)
{
  mpopulation[n,]=0
}


### Writing marks on mpopulation
for (cell in 1:populationSize){
  depositedchr<-createNakedChromatin(chromlength)
  depositedchr<-deposit(chromatin,chromlength)[["chr"]]
  i=1
  while (i <= prc2)
  {
    mpopulation[,(i*4)-3]=mpopulation[,(i*4)-3]+depositedchr$me0
    mpopulation[,(i*4)-2]=mpopulation[,(i*4)-2]+depositedchr$me1
    mpopulation[,(i*4)-1]=mpopulation[,(i*4)-1]+depositedchr$me2
    mpopulation[,(i*4)]=mpopulation[,(i*4)]+depositedchr$me3
    ### Check for mitosis
    if (runif(n=1,min=0,max=1) >= 0.9){depositedchr<-mitosis(depositedchr);i=1}
    print (paste("Mitosis - ",cell,":",i,sep=""))
    depositedchr<-deposit(depositedchr,chromlength)[["chr"]]
    i=i+1
  }
}


for (i in c(1:prc2))
{
  png(paste(i,".mitosis.png",sep=""),width = 1200,height = 600)
  par(mfrow=c(3,2))
  # plot((population[,(i*4)-3])/populationSize,type="h",ylim=c(0,1),main=paste("K36me0 prc2 round: ",i,sep=""))
  plot((population[,(i*4)-2])/populationSize,type="h",ylim=c(0,1),main=paste("K36me1 prc2 round: ",i,sep=""))
  plot((mpopulation[,(i*4)-2])/populationSize,type="h",ylim=c(0,1),main=paste("Mitosis: K36me1 prc2 round: ",i,sep=""))
  
  plot((population[,(i*4)-1])/populationSize,type="h",ylim=c(0,1),main=paste("K36me2 prc2 round: ",i,sep=""))
  plot((mpopulation[,(i*4)-1])/populationSize,type="h",ylim=c(0,1),main=paste("Mitosis: K36me2 prc2 round: ",i,sep=""))
  
  plot((population[,(i*4)-0])/populationSize,type="h",ylim=c(0,1),main=paste("K36me3 prc2 round: ",i,sep=""))
  plot((mpopulation[,(i*4)-0])/populationSize,type="h",ylim=c(0,1),main=paste("Mitosis: K36me3 prc2 round: ",i,sep=""))
  
  dev.off()
  
}

