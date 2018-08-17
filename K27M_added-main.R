#This version includes addition of K27M and including these mutations.
rm(list=ls())
this_version="A3.3-K27M"

## Home directory on your computer, please make sure you have permission to write
homedir="~/Desktop/Desktop_June_2018/Histone_Mark_Simulation/Steadystate/"


## Please do not change this, this is to avoid the scietific connotation in the results
options(scipen=999)

## Length of the chromatin (i.e. number of uncleosomes)
chromlength=1000

## Population size (number of chromatins)
populationSize=100

### Number of cycles (maximum) before the simulation stops. This will reach only if we don't get to a steady state before
life=100

## Transition probabilities, the first number is the current state the second number os the next state, e.g. Tp01 means the probability of transitioning from K27me0 to K27me1 and so on.
Tp01=0.99
Tp03=0
Tp12=0.99
Tp13=0.0
Tp02=0.0
Tp23=0.05
Tp00=0
Tp11=0
Tp22=0

### The probability of mitosis happening at any move of PRC2
mitosis_prob=0.0001

## Probability of PRC2 does its function (depositin); 0-1, 1 means everytime if he sits on a nucleosome and if everything else (like other marks etc) look suitable, it does its job. Zero means it won't even if everything else looks okay.
depop=1

### Distribution of the genes, and other mark domains. It can be a list of several domains. Be careful with overlaps and going over the chromosome length
genic_structure=list(c(100:300),c(800:1000))
expression_structure=list(c(100:300))
k36me2_structure=list(c(400:600))
k36me3_structure=list(c(100:300))
K27Mprobability=0.05
## The value below shows how long K27M in scale of mitosis clock will stall the PRC2. In other words, if the probability of mitosis is 1/1000 and the value below is 1/10, that means K27M will stall PRC2 100 times more than its regular speed which is 1000 rounds per mitosis.
how_long_K27M_stalls=0.1
### How hard K36me2 affects the deposition of K27me3 (zero= completely prevents it, 1= not affecting at all)
howhardk36me2=0.1

### How hard K36me3 affects the deposition of K27me2 and K27me3 (zero= completely prevents it, 1= not affecting at all)
howhardk36me3=0

#### Peakiness of K36me2 and K36me3 [0,0.1] (0 is a uniform block, 0.1 is a very sharp peak)
peakiness=0.0005

#### Slope of PRC2 falling off [0,0.1] (0 is a non-falling off, 0.1 is a very sharp slope)
prc2slop=0.0025
#### Defining the probability of PRC2 falling off chromatin based on nucleosomeNumber and length of the chromosome
prc2hl <- function(nucleosomeNumber)
{
  # PRC2 half-life
  #  hl<-(1-0.004)^nucleosomeNumber
  hl<-(1-prc2slop)^nucleosomeNumber
  
  return(hl)  
}

### Testin prc2hl
d<-c()
for (z in 1:1000)
{
  d<-c(d,prc2hl(z))
}
plot(d, main = paste("slope= ",prc2slop,sep = ""))







### This factor indictes how big the previous nucleosome's status affects the deposition of K27me3. i.e. 2 means it makes it twice as likely etc.
neighborK27me3bonus=10

### Set this to zero if you want to avoid random walk movement for PRC2
randomwk=0
## The maximum number of steps for PRC2 in case of random walk before it falls off the chromatin
maximumsteps=10000
## The maximum size of each step in random walk 
stepsize=10





newdir=paste(homedir,this_version,"-",gsub(":","_",gsub(" ","-",date())),"-","chrmlngth_",chromlength,"-pop_",populationSize,"-add01-12-23_0213-",Tp01,"-",Tp12,"-",Tp23,"-",Tp02,"-mitosis_prob_",mitosis_prob,"_Randomwalk_",randomwk,"_neighbor_",neighborK27me3bonus,"x","/",sep = "")
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
  chromatin<-data.frame(integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize))
  colnames(chromatin)<-c("N","S","me0","me1","me2","me3","K36me2","K36me3","genic","expression","K27M")
  for (s in 1:chromSize)
  {
    chromatin$N[s]=s
    chromatin$S[s]=0
    chromatin$me0[s]=1
    chromatin$me1[s]=0
    chromatin$me2[s]=0
    chromatin$me3[s]=0
    chromatin$K36me2[s]=0
    chromatin$K36me3[s]=0
    chromatin$genic[s]=0
    chromatin$expression[s]=0
    chromatin$K27M=0
  }
  results<-list("add1"=0,"add2"=0,"add3"=0,"add0"=0,"rando"=0,"chr"=chromatin)
  return(results)
}
chromatin<-createNakedChromatin(chromlength)

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


#### Defining regions of K36me2 and K36me3
### Define the type of peak
shapeOfK36mark<-function(start,end)
{
  d1<-c()
  k36memarklength=(end-start)/2
  for (i in 1:k36memarklength)
  {
    d1<-c(d1,(1-peakiness)^(i))
  }
  d2<-rev(d1)
  d<-c(d2,d1)
  return(d)
}
### deposit the k36me2 mark based on the peak model
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
      #chrom$K36me2[q]=shapeOfK36mark(start,end)[i]
      if (!is.na(shapeOfK36mark(start,end)[i]) && shapeOfK36mark(start,end)[i] > runif(n=1,min=0,max=1))
      {
        chrom$K36me2[q]=1
      } else { chrom$K36me2[q]=0 }
      i=i+1
    }
  }
  chrom$K36me2[is.na(chrom$K36me2)] <- 0
  return(chrom)
}


### deposit the k36me3 mark based on the peak model
depositk36me3<-function(chrom,k36me_structure)
{
  for (k in k36me3_structure)
  {
    start=k[1]
    end=k[length(k36me3_structure[k])]
    # print(start)
    # print(end)
    i=1
    for (q in start:end)
    {
      #chrom$K36me3[q]=shapeOfK36mark(start,end)[i]
      if (!is.na(shapeOfK36mark(start,end)[i]) && shapeOfK36mark(start,end)[i] > runif(n=1,min=0,max=1))
      {
        chrom$K36me3[q]=1
      } else { chrom$K36me3[q]=0 }
      i=i+1
    }
  }
  chrom$K36me3[is.na(chrom$K36me3)] <- 0
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
# test<-chromatin
# test[["chr"]]$S=1
# test[["chr"]]$me1=1
# test[["chr"]]$me0=0
# test<-mitosis(test[["chr"]])
#View(test)


### Add K27M ###
### Note, that the first K27M (since it's on H3.3 and H3.3 usually doesn't exist in K27me3 domains), should be far down from the PRC2 nucleation site). I use the arbitrary distance of 50 nucleosomes 
nucleation_distance=50
addk27m <- function(chrom,K27Mprobability)
{
  for (i in nucleation_distance:nrow(chrom)) 
  {
    check=runif(n=1,min = 0,max=1)
    if (check < K27Mprobability)
    {
      chrom$K27M[i]=1
    }
  }
  return (chrom)
}




### Preparing the naked chromosome for K27m123 deposition by adding K36me2 marks as well as defining genic regions and expression in different regions defined above

chromatin[["chr"]]<-depositgene(chromatin[["chr"]],genic_structure)
chromatin[["chr"]]<-depositexpression(chromatin[["chr"]],expression_structure)
chromatin[["chr"]]<-depositk36me2(chromatin[["chr"]],k36me2_structure)
chromatin[["chr"]]<-depositk36me3(chromatin[["chr"]],k36me3_structure)
chromatin[["chr"]]<-addk27m(chromatin[["chr"]],K27Mprobability)

#sample(c(1,chromlength))
### Depositing marks

randomwalkme<-function(chromlength,maxstepsize,maxsteps)
{
  d<-c()
  distance=0
  chromlength=chromlength
  i=1
  while(abs(distance) <= chromlength && i < maxsteps){
    right=round(runif(n=1,min=1,max=maxstepsize),digits = 0)
    left=round(runif(n=1,min=-1*(maxstepsize),max=-1),digits = 0)
    distance=distance+sample(x = c(left,right),size = 1,replace = T,prob = c(0.5,0.5))
#    print (i)
#    print (distance)
    d<-c(d,distance)
    i=i+1
  }
  d<-subset(d,d[]>0 & d[] <=chromlength)
  if(length(d)<1) {d=1}
  return(d)
}
 # test<-randomwalkme(1000,10,10000)
 # plot(test,type="l",ylab = paste("Distance; max dist = ",max(test),sep=""),xlab = paste("Number of Steps in the positive range; last step = ",length(test),sep=""))

deposit <- function(chr,chromSize)
{
mitosiscount=0
## If randomwalk is chosen for PRC2 move, the "howPRC2moves" variable which is an vector made of randomwalked variables will be the order of PRC2 moves, otherwise, if randomwalk is set to zero, it'll be a vector of 0 to size of the chromatin.
if(randomwk==0){
  howPRC2moves=c(1:chromSize)
}
else {howPRC2moves=randomwalkme(1000,stepsize,maximumsteps)}

## Based on what vector we have from previous part, the deposition starts from the first element of the vector to the last element. note the vector's name is howPRC2moves
    for (nucleosomeNumber in howPRC2moves)
  {
    ### First, checking if the nucleosome contains K27M, if so, we give the cell the chance to go through mitosis with a chance of 10%
      
      if (chr$K27M[nucleosomeNumber]==1)
      {
        mitosisCheck=how_long_K27M_stalls
        ## The random number below is chosen to be compared to 0.5 on each nucleosome that has a H3.3K27M. If the random number is below 0.5, that means PRC2 still can methylate the other H3.3 from the two in each nucleosome. Basically it writes the methyl mark half of the times on nucleosomes with K27M.
        k27m=runif(n=1,min=0,max=1)
      } else {
        mitosisCheck=runif(n=1,min=0,max=1)
        k27m=0
        }
      
      
      
    #choosing a random number between 0 and 1 for testing with the probability of each scenario
    rando=runif(n=1,min=0,max=1)
    ## Calculating the probability of PRC2 falling off chromatin based on its position 
    falloff=prc2hl(nucleosomeNumber)
    ## Checking the status of the previous nucleosome's K27me3. If we are at the first nucleosome, we put 1 for the state otherwise if it's larger than one, we check the previous one.
    if (nucleosomeNumber > 1){
    prevnucK27me3state=chr$me3[(nucleosomeNumber)-1]
    } else {prevnucK27me3state=1}

    neighborK27me3bonus=1
    # If the previous nucleosome is K27me3 (or if we are at the first uncleosome), the bonus is as defined in the parameters on the header of the script otherwise, it will be equal to one (no bous)
     if (prevnucK27me3state == 1) 
     {
       neighborK27me3bonus=neighborK27me3bonus
     }
     else {
       neighborK27me3bonus=1
       }
    # If the current nucleosome has K36me2 and K36me3, the penalties will be applied as defined in the parameters on the header of the script otherwise, it will be equal to one (no penalty)    
    if (chr$K36me2[nucleosomeNumber]==1)
    {penaltyk36me2=howhardk36me2}
    else {penaltyk36me2=1}
    if (chr$K36me3[nucleosomeNumber]==1)
    {penaltyk36me3=howhardk36me3}
    else {penaltyk36me3=1}

    #When PRC2 makes it to each nucleosome, checks the current state (S). If the current state shows no K27me (1,2, or 3)   
    if (chr$S[nucleosomeNumber] == 0)
    {
      #If the random generate number falls between zero and Tp00 (times the chance of falling off) nothing changes
      if (rando < Tp00*falloff && k27m<0.5)
      {
        chr$S[nucleosomeNumber]=0
        chr$me0[nucleosomeNumber]=1
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=0
      } 
      # If it falls between zero and Tp01(*fall off chance), changes the S from 0 to 1
      else if(rando < Tp01*falloff && k27m<0.5) {
        chr$S[nucleosomeNumber]=1
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=1
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=0
      } 
      # If it falls between zero and Tp02(*fall off chance), changes the S from 0 to 2
    
      else  if(rando < Tp02*falloff  && k27m<0.5) {
        chr$S[nucleosomeNumber]=2
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=1
        chr$me3[nucleosomeNumber]=0
      } 
      # If it falls between zero and Tp03(*fall off chance), changes the S from 0 to 3
      else  if(rando < Tp03*falloff  && k27m<0.5) {
        chr$S[nucleosomeNumber]=3
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=1
      }
    }
    ## Same story as above if the current state is 1
    else if (chr$S[nucleosomeNumber] == 1)  {
      if(rando < Tp11*falloff && k27m<0.5)
      {
        chr$S[nucleosomeNumber]=1
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=1
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=0
      } else if(rando < Tp12*falloff*penaltyk36me3 && k27m<0.5) {
        chr$S[nucleosomeNumber]=2
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=1
        chr$me3[nucleosomeNumber]=0
      } else  if(rando < Tp13*falloff*penaltyk36me3 && k27m<0.5) {
        chr$S[nucleosomeNumber]=3
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=1
      } 
    } 
    ## Same story as above if the current state is 2, except here we have the bonus for the previous nucleosome having K27me3 mark
    else if (chr$S[nucleosomeNumber] == 2)  {
      if(rando < Tp22*falloff && k27m<0.5)
      {
        chr$S[nucleosomeNumber]=2
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=1
        chr$me3[nucleosomeNumber]=0
      } else if(rando < Tp23*falloff*penaltyk36me2*penaltyk36me3*neighborK27me3bonus && k27m<0.5)  {
        chr$S[nucleosomeNumber]=3
        chr$me0[nucleosomeNumber]=0
        chr$me1[nucleosomeNumber]=0
        chr$me2[nucleosomeNumber]=0
        chr$me3[nucleosomeNumber]=1
      }
    }
# mitosis possibility at each nucleosome, there is a chance that mitosis happens which is very slim (0.001) but considering the length of the chromatin, on average is 1
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
k36me2avg<-c()
k36me3avg<-c()
k27mavg<-c()

##### Cycles #######
# Deciding where the K27M mutated nucleosomes will appear on the chromatin
thispopulationk27mdistribution<-chromatin[["chr"]]$K27M

totalmitosis=0

while(lifecycle <= life)
{
## Creating a place-holder for future plotting of the population at each life-cycle
temp<-data.frame(integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength),integer(chromlength))
colnames(temp)<-c("N","S","me0","me1","me2","me3","K36me2","K36me3","genic","expression","K27M")

  for (s in 1:chromlength)
  {
    temp$N[s]=s
    temp$S[s]=0
    temp$me0[s]=0
    temp$me1[s]=0
    temp$me2[s]=0
    temp$me3[s]=0
    temp$K36me2[s]=0
    temp$K36me3[s]=0
    temp$genic[s]=0
    temp$expression[s]=0
    temp$K27M[s]=0
  }



  
  print(lifecycle)
  for (p in 1:populationSize)
  {
  if (lifecycle==1) 
    {
    population[[p]]<-createNakedChromatin(chromlength)
    population[[p]][["chr"]]<-depositk36me2(population[[p]][["chr"]],k36me2_structure)
    population[[p]][["chr"]]<-depositk36me3(population[[p]][["chr"]],k36me3_structure)
    population[[p]][["chr"]]<-depositgene(population[[p]][["chr"]],genic_structure)
    population[[p]][["chr"]]<-depositexpression(population[[p]][["chr"]],expression_structure)
    population[[p]][["chr"]]$K27M<-thispopulationk27mdistribution
    
    } else if (lifecycle > 1) {
    if (runif(n=1,min=0,max=1) <= depop)
      {
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
 par(mfrow=c(9,1))
 plot(population[[1]][["chr"]]$me0,type="h",ylim=c(0,1),ylab = "K27me0", main = paste("Cycle: ",lifecycle," Individual: 1",sep = ""),xlab = "Chromatin")
 plot(population[[1]][["chr"]]$me1,type="h",ylim=c(0,1),ylab = "K27me1",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$me2,type="h",ylim=c(0,1),ylab = "K27me2",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$me3,type="h",ylim=c(0,1),ylab = "K27me3",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$K36me2,type="h",ylim=c(0,1),ylab = "K36me2",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$K36me3,type="h",ylim=c(0,1),ylab = "K36me3",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$K27M,type="h",ylim=c(0,1),ylab = "K27M",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$expression,type="h",ylim=c(0,1),ylab = "Expression",xlab = "Chromatin")
 plot(population[[1]][["chr"]]$genic,type="h",ylim=c(0,1),ylab = "Genic/Intergenic",xlab = "Chromatin")
 
 dev.off()
# 

png(paste(populationSize,"_pop.","-cycles_",lifecycle,"-chromsize_",chromlength,".png",sep=""),width = 800,height = 1000)
par(mfrow=c(9,1))
plot((currentpop[["chr"]]$me0)/(populationSize),type="h",ylab = "K27me0", main = paste("Cycle: ",lifecycle,sep = ""),xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$me1)/(populationSize),type="h",ylab = "K27me1", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$me2)/(populationSize),type="h",ylab = "K27me2", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$me3)/(populationSize),type="h",ylab = "K27me3", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$K36me2)/(populationSize),type="h",ylab = "K36me2", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$K36me3)/(populationSize),type="h",ylab = "K36me3", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$K27M)/(populationSize),type="h",ylab = "K27M", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$expression)/(populationSize),type="h",ylab = "Expression", xlab = "Chromatin",ylim=c(0,1))
plot((currentpop[["chr"]]$genic)/(populationSize),type="h",ylab = "Genic/Intergenic", xlab = "Chromatin",ylim=c(0,1))
dev.off()

me0avg=c(me0avg,mean((currentpop[["chr"]]$me0)/populationSize))
me1avg=c(me1avg,mean((currentpop[["chr"]]$me1)/populationSize))
me2avg=c(me2avg,mean((currentpop[["chr"]]$me2)/populationSize))
me3avg=c(me3avg,mean((currentpop[["chr"]]$me3)/populationSize))
k36me2avg=c(k36me2avg,mean((currentpop[["chr"]]$K36me2)/populationSize))
k36me3avg=c(k36me3avg,mean((currentpop[["chr"]]$K36me3)/populationSize))
k27mavg=c(k27mavg,mean((currentpop[["chr"]]$K27M)/populationSize))


if (lifecycle > 10){
var0=var(me0avg[lifecycle:(lifecycle-10)])
var1=var(me1avg[lifecycle:(lifecycle-10)])
var2=var(me2avg[lifecycle:(lifecycle-10)])
var3=var(me3avg[lifecycle:(lifecycle-10)])

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


# To be added:
# 5% aporadic K27M (done)
# increasing the chance of mitosis once hitting K27M to 10% thus every K27M stalls the process as long as 10% of a cell cycle. so whenever prc2 hits k27m, subject the chromosome to mitosis (with 10% chance)
# check the sum of probability at each position, it has to come to 1
