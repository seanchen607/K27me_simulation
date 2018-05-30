rm(list=ls())
#for (Tp23 in c(0.05))
#{
options(scipen=999)
chromlength=1000
prc2=1
populationSize=1000
life=100
Tp01=0.9
Tp02=0
Tp03=0
Tp12=0.5
Tp13=0
Tp23=0.05
Tp00=0
Tp11=0
Tp22=0
mitosis_prob=0.001
depop=0.9

newdir=paste("~/Desktop/Histone_Mark_Simulation/Steadystate/",gsub(":","_",gsub(" ","-",date())),"-","chrmlngth_",chromlength,"-prc2_",prc2,"-pop_",populationSize,"-add01-12-23_",Tp01,"-",Tp12,"-",Tp23,"-mitosis_prob_",mitosis_prob,"/",sep = "")
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




#### Defining mitosis function
mitosis <- function (chrom){
  mitochromatin<-chrom
  for (i in 1:nrow(chrom)) {
    check=runif(n=1,min = 0,max=1)
    if (check > 0.5)
    {
      mitochromatin$S[i]=0
      mitochromatin$me0[i]=1
      mitochromatin$me1[i]=0
      mitochromatin$me2[i]=0
      mitochromatin$me3[i]=0
    }
  }
  return (mitochromatin)
}
## Testing mitosis
test<-chromatin
test[["chr"]]$S=1
test[["chr"]]$me1=1
test[["chr"]]$me0=0
test<-mitosis(chromatin[["chr"]])
#View(test)


#### Defining the probability of PRC2 falling off chromatin based on nucleosomeNumber and length of the chromosome
prc2hl <- function(nucleosomeNumber)
{
  # PRC2 half-life
  hl<-(1-0.002)^nucleosomeNumber
  return(hl)  
}

### Testin prc2hl
d<-c()
for (z in 1:1000)
{
d<-c(d,prc2hl(z))
}
#plot(d)



### Depositing marks
deposit <- function(chr,chromSize)
{
mitosiscount=0
  
  for (nucleosomeNumber in 1:chromSize)
  {
    rando=runif(n=1,min=0,max=1)
    falloff=prc2hl(nucleosomeNumber)

    if (chr$S[nucleosomeNumber] == 0)
    {
      if (rando < Tp00*falloff)
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
# mitosis possibility at each nucleosome, there is a chance that mitosis happens which is very slim (0.001) but considering the length of the chromatin, on average is 1
    mitosisCheck=runif(n=1,min=0,max=1)
    if (mitosisCheck <= mitosis_prob)
    {
    chr<-mitosis(chr)
    mitosiscount=mitosiscount+1
#    print(paste("mitosis in individual",p," lifecyle: ",lifecycle,sep=""))
# 50% chance that the daughter chromatin receives the prc2, if it doesn't the new prc2 will start from first nucleosome again otherwise, the prc2 continues from where it was interrupted.
      if(runif(n=1,min=0,max=1) > 0.5)
      {
      nucleosomeNumber=1
      }
    } 
    # else if (mitosisCheck > mitosis_prob) {
    #   #if the mitosis doesn't happen, the output of the deposition chromatin would be the same as what it accumulated so far without slashing it in half
    #   chr2<-chr
    # }
  } 


  results<-list("Tp01"=Tp01,"Tp12"=Tp12,"Tp23"=Tp12,"rando"=rando,"chr"=chr,"numberofmitosis"=mitosiscount)
  return(results)
}




# creating the initial population 
population<-list()
currentpop<-list()  
currentpop<-list("lifecycle"=1,"chr"=createNakedChromatin(chromlength)[["chr"]])  
currentpop[["chr"]]$me0=0

lifecycle=1
me0avg<-c()
me1avg<-c()
me2avg<-c()
me3avg<-c()

##### Cycles #######
totalmitosis=0
while(lifecycle <= life)
{
## Creating a place-holder for future plotting of the population at each life-cycle
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
  if (lifecycle==1) 
    {
    population[[p]]<-createNakedChromatin(chromlength)
    } else if (lifecycle > 1) {
    if (runif(n=1,min=0,max=1) <= depop)
      {
      population[[p]]<-deposit(population[[p]][["chr"]],chromlength)
      totalmitosis=totalmitosis+population[[p]][["numberofmitosis"]]
      }
    }
  temp<-temp+population[[p]][["chr"]]
#  totalmitosis=totalmitosis+population[["numberofmitosis"]]
  
  }
if (lifecycle > 1){
currentpop=list("lifecycle"=c(currentpop[["lifecycle"]],lifecycle),"chr"=temp)
} else {
  currentpop=list("lifecycle"=lifecycle,"chr"=temp)
}
  
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

if (lifecycle > 5){
var0=var(me0avg[lifecycle:(lifecycle-5)])
var1=var(me1avg[lifecycle:(lifecycle-5)])
var2=var(me2avg[lifecycle:(lifecycle-5)])
var3=var(me3avg[lifecycle:(lifecycle-5)])

if (var0 < 0.0001 && var1 < 0.0001 && var2 < 0.0001 && var3 < 0.0001 )
{
  lastlap=lifecycle
#  lifecycle=life
  break
}
}

lifecycle=lifecycle+1

currentpop[["chr"]]$N=population[[p]][["chr"]]$N

}

mitonum=totalmitosis
print(paste(mitonum," mitosis events happened during this simulation.",sep=""))
png(paste("4plots-",populationSize,"_pop.","-cycles_",lastlap,"-chromsize_",chromlength,".png",sep=""),width = 800,height = 600)
par(mfrow=c(2,2))
plot(me0avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K36me0",xlab = "Cycle")
text(x = (lastlap/10)+6,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+6),y = 0.85,labels = paste("Average mark: ",me0avg[lastlap],sep = ""),pos = 1,offset = 0)
plot(me1avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K36me1",xlab = "Cycle")
text(x = (lastlap/10)+6,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+6),y = 0.85,labels = paste("Average mark: ",me1avg[lastlap],sep = ""),pos = 1,offset = 0)
plot(me2avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K36me2",xlab = "Cycle")
text(x = (lastlap/10)+6,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+6),y = 0.85,labels = paste("Average mark: ",me2avg[lastlap],sep = ""),pos = 1,offset = 0)
plot(me3avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K36me3",xlab = "Cycle")
text(x = (lastlap/10)+6,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+6),y = 0.85,labels = paste("Average mark: ",me3avg[lastlap],sep = ""),pos = 1,offset = 0)
dev.off()
#}
