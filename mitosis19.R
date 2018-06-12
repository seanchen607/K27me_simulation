rm(list=ls())
# for (Tp01 in c(0.97))
# {
#   for(Tp12 in c(0.8))
#     {
#     for(Tp23 in c(0.2))
#       {
#       for (Tp02 in c(0.1))
#       {
options(scipen=999)
chromlength=1000
prc2=1
populationSize=100
life=100
Tp01=0.95
Tp03=0
Tp12=0.8
Tp13=0.05
Tp02=0.05
Tp23=0.2
Tp00=0
Tp11=0
Tp22=0
mitosis_prob=0.0001
depop=1
genic_structure=list(c(200:250),c(400:550),c(650:800))
expression_structure=list(c(400:550),c(650:800))
k36me2_structure=list(c(400:550),c(650:800))


 newdir=paste("~/Desktop/Histone_Mark_Simulation/Steadystate/",gsub(":","_",gsub(" ","-",date())),"-","chrmlngth_",chromlength,"-prc2_",prc2,"-pop_",populationSize,"-add01-12-23_0213-",Tp01,"-",Tp12,"-",Tp23,"-",Tp02,"-mitosis_prob_",mitosis_prob,"/",sep = "")
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
  chromatin<-data.frame(integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize))
  colnames(chromatin)<-c("N","S","me0","me1","me2","me3","K36me2","genic","expression")
  for (s in 1:chromSize)
  {
    chromatin$N[s]=s
    chromatin$S[s]=0
    chromatin$me0[s]=1
    chromatin$me1[s]=0
    chromatin$me2[s]=0
    chromatin$me3[s]=0
    chromatin$K36me2[s]=0
    chromatin$genic[s]=0
    chromatin$expression[s]=0
  }
  results<-list("add1"=0,"add2"=0,"add3"=0,"add0"=0,"rando"=0,"chr"=chromatin)
  return(results)
}
chromatin<-createNakedChromatin(chromlength)

### uniform distribution of k36me2
# addk36me2 <- function(chromatin,k36me2_structure)
# {
#      for (i in 1:length(k36me2_structure))
#      {
#        chromatin[["chr"]]$K36me2[k36me2_structure[i]]=1
#      }
#      return(chromatin)
# }
 
   
## Add genic structure
depositgene<-function(chrom,genic_structure)
{
  for (k in genic_structure)
  {
    start=k[1]
    end=k[length(genic_structure[k])]
    for (q in start:end)
    {
      chrom$genic[q]=1
    }
  }
  return(chrom)
}


## add uniform expression
depositexpression<-function(chrom,expression_structure)
{
  for (k in expression_structure)
  {
    start=k[1]
    end=k[length(expression_structure[k])]
    for (q in start:end)
    {
      chrom$expression[q]=1
    }
  }
  return(chrom)
}


#### Defining regions of K36me2
### Define the type of peak
distk36me2<-function(start,end)
{
  d1<-c()
  k36me2length=(end-start)/2
  for (i in 1:k36me2length)
  {
    d1<-c(d1,(1-0.02)^(i))
  }
  d2<-rev(d1)
  d<-c(d2,d1)
  return(d)
}
### deposit the mark based on the peak model
depositk36me2<-function(chrom,k36me2_structure)
{
  for (k in k36me2_structure)
  {
    start=k[1]
    end=k[length(k36me2_structure[k])]
    # print(start)
    # print(end)
    i=1
    for (q in start:end)
    {
      chrom$K36me2[q]=distk36me2(start,end)[i]
      i=i+1
    }
  }
  chrom$K36me2[is.na(chrom$K36me2)] <- 0
  return(chrom)
}





#### Defining mitosis function
mitosis <- function (chrom){
  mitochromatin<-chrom
  for (i in 1:nrow(chrom)) 
  {
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
test<-mitosis(test[["chr"]])
#View(test)


#### Defining the probability of PRC2 falling off chromatin based on nucleosomeNumber and length of the chromosome
prc2hl <- function(nucleosomeNumber)
{
  # PRC2 half-life
  hl<-(1-0.002)^nucleosomeNumber
  return(hl)  
}

### Testin prc2hl
# d<-c()
# for (z in 1:1000)
# {
# d<-c(d,prc2hl(z))
# }
# plot(d)





### Preparing the naked chromosome for K27m123 deposition by adding K36me2 marks as well as defining genic regions and expression in different regions defined above

chromatin[["chr"]]<-depositgene(chromatin[["chr"]],genic_structure)
chromatin[["chr"]]<-depositexpression(chromatin[["chr"]],expression_structure)
chromatin[["chr"]]<-depositk36me2(chromatin[["chr"]],k36me2_structure)


### Depositing marks
deposit <- function(chr,chromSize)
{
mitosiscount=0
  
  for (nucleosomeNumber in 1:chromSize)
  {
    rando=runif(n=1,min=0,max=1)
    falloff=prc2hl(nucleosomeNumber)

    # if(chr$K36me2[nucleosomeNumber] == 1 )
    # {
    #   #penaltyk36me2=mysamp(n=1, m=0.5, s=0.1, lwr=0, upr=1, nnorm=100)
    #   penaltyk36me2=0
    # } else { penaltyk36me2=1 }
    penaltyk36me2=1-(chr$K36me2[nucleosomeNumber])
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
      } else if(rando < Tp12*falloff*penaltyk36me2) {
        chr$S[nucleosomeNumber]=2
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=1
        chr$me3[nucleosomeNumber]=0
      } else  if(rando < Tp13*falloff*penaltyk36me2) {
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
      } else if(rando < Tp23*falloff*penaltyk36me2)  {
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
temp<-data.frame(integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength))
colnames(temp)<-c("N","S","me0","me1","me2","me3","K36me2","genic","expression")

  for (s in 1:chromlength)
  {
    temp$N[s]=s
    temp$S[s]=0
    temp$me0[s]=0
    temp$me1[s]=0
    temp$me2[s]=0
    temp$me3[s]=0
    temp$K36me2[s]=0
    temp$genic[s]=0
    temp$expression[s]=0
  }



# temp$K36me2<-depositk36me2(temp,k36me2_structure)$K36me2
# temp$genic<-depositgene(temp,genic_structure)$genic
# temp$expression<-depositexpression(temp,expression_structure)$expression
# for (s in 1:length(k36me2_structure))
# {
#   temp$K36me2[k36me2_structure[s]]=1
#   temp$genic[genic_structure[s]]=1
#   temp$expression[expression_structure[s]]=1
# }

# for (i in K36me2start:K36me2end)
# {
#   temp$K36me2[s]=1
# }
  
  print(lifecycle)
  for (p in 1:populationSize)
  {
  if (lifecycle==1) 
    {
    population[[p]]<-createNakedChromatin(chromlength)
    population[[p]][["chr"]]<-depositk36me2(population[[p]][["chr"]],k36me2_structure)
    population[[p]][["chr"]]<-depositgene(population[[p]][["chr"]],genic_structure)
    population[[p]][["chr"]]<-depositexpression(population[[p]][["chr"]],expression_structure)
    } else if (lifecycle > 1) {
    if (runif(n=1,min=0,max=1) <= depop)
      {
      #population[[p]][["chr"]]$K36me2<-population[[p]][["chr"]]$K36me2+depositk36me2(population[[p]][["chr"]],k36me2_structure)$K36me2
      population[[p]][["chr"]]$K36me2<-(population[[p]][["chr"]]$K36me2+((depositk36me2(population[[p]][["chr"]],k36me2_structure)$K36me2)*(lifecycle-1)))/(lifecycle)
      
      #population[[p]][["chr"]]$genic<-population[[p]][["chr"]]$genic+depositgene(population[[p]][["chr"]],genic_structure)$genic
      population[[p]][["chr"]]$genic<-(population[[p]][["chr"]]$genic+((depositgene(population[[p]][["chr"]],genic_structure)$genic)*(lifecycle-1)))/(lifecycle)
      
      #population[[p]][["chr"]]$expression<-population[[p]][["chr"]]$expression+depositexpression(population[[p]][["chr"]],expression_structure)$expression
      population[[p]][["chr"]]$expression<-(population[[p]][["chr"]]$expression+((depositexpression(population[[p]][["chr"]],expression_structure)$expression)*(lifecycle-1)))/(lifecycle)
      
      population[[p]]<-deposit(population[[p]][["chr"]],chromlength)
      
      totalmitosis=totalmitosis+population[[p]][["numberofmitosis"]]
      }
    }
  temp<-temp+population[[p]][["chr"]]
  }
  
if (lifecycle > 1){
currentpop=list("lifecycle"=c(currentpop[["lifecycle"]],lifecycle),"chr"=temp)
} else {
  currentpop=list("lifecycle"=lifecycle,"chr"=temp)
}
  
 png(paste("cell1/",lifecycle,".ind.","1","-cycles",lifecycle,"size",chromlength,".png",sep=""),width = 800,height = 1000)
 par(mfrow=c(7,1))
 plot(population[[1]][["chr"]]$me0,type="h",ylim=c(0,1),ylab = "K27me0", main = paste("Cycle: ",lifecycle," Individual: 1",sep = ""),xlab = "Chromatin")
 plot(population[[1]][["chr"]]$me1,type="h",ylim=c(0,1),ylab = "K27me1",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$me2,type="h",ylim=c(0,1),ylab = "K27me2",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$me3,type="h",ylim=c(0,1),ylab = "K27me3",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$K36me2,type="h",ylim=c(0,1),ylab = "K36me2",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$expression,type="h",ylim=c(0,1),ylab = "Expression",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$genic,type="h",ylim=c(0,1),ylab = "Genic/Intergenic",xlab = "Chromatin")
 
 dev.off()
# 

png(paste(populationSize,"_pop.","-cycles_",lifecycle,"-chromsize_",chromlength,".png",sep=""),width = 800,height = 1000)
par(mfrow=c(7,1))
plot((currentpop[["chr"]]$me0)/(populationSize),type="h",ylab = "K27me0", main = paste("Cycle: ",lifecycle,sep = ""),xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$me1)/(populationSize),type="h",ylab = "K27me1", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$me2)/(populationSize),type="h",ylab = "K27me2", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$me3)/(populationSize),type="h",ylab = "K27me3", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$K36me2)/(populationSize),type="h",ylab = "K36me2", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$expression)/(populationSize),type="h",ylab = "Expression", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$genic)/(populationSize),type="h",ylab = "Genic/Intergenic", xlab = "Chromatin",ylim=c(0,1))
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

if (var0 < 0.001 && var1 < 0.001 && var2 < 0.001 && var3 < 0.001 )
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
plot(me0avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K27me0",xlab = "Cycle")
text(x = (lastlap/10)+6,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+6),y = 0.85,labels = paste("Average mark: ",me0avg[lastlap],sep = ""),pos = 1,offset = 0)
plot(me1avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K27me1",xlab = "Cycle")
text(x = (lastlap/10)+6,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+6),y = 0.85,labels = paste("Average mark: ",me1avg[lastlap],sep = ""),pos = 1,offset = 0)
plot(me2avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K27me2",xlab = "Cycle")
text(x = (lastlap/10)+6,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+6),y = 0.85,labels = paste("Average mark: ",me2avg[lastlap],sep = ""),pos = 1,offset = 0)
plot(me3avg,type = "l",ylim=c(0,1),xlim = c(1,lastlap), main="Average H3K27me3",xlab = "Cycle")
text(x = (lastlap/10)+6,y = 1,labels = paste("Stopped at cycle: ",lastlap," - #Mitosis: ",mitonum,sep = ""),pos = 1,offset = 0)
text(x = round((lastlap/10)+6),y = 0.85,labels = paste("Average mark: ",me3avg[lastlap],sep = ""),pos = 1,offset = 0)
dev.off()
#       }
#     }
#   }
# }
print("Don't forget about replinishing K36me2 and also its distribution in an area shouldn't be uniform. they also have peaks no matter how broad. also check if the prc2 should fall of chromatin once seeing k36me2 (and thus not printing anything afterwards or just not printing on nucleosomes where there is k36me2")
