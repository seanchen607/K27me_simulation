## 1- Preparing Simulation
rm(list=ls())
# this_version = a prefix to the results' directory which contains information on the parameters and a version name. Just to keep tracks of versions.
this_version="A6.1-Constant-Clock-H3.3-Only-in-OpenChr-Prob-K27M"

# homdir= any directory on your computer where you want to store the simulation results and plots, please make sure you have permission to write
homedir="~/Desktop/Desktop_June_2018/Histone_Mark_Simulation/Constant_clock/"

# Please do not change this, this is to avoid the scietific connotation in the results
options(scipen=999)


## 2- Defining parameters
# Length of the chromatin (i.e. number of uncleosomes)
chromlength=1000

# Population size (number of chromatins)
populationSize=100

# Number of cycles (maximum) before the simulation stops. This will reach only if we don't get to a steady state before
life=100

# Transition probabilities, the first number is the current state the second number os the next state, e.g. Tp01 means the probability of transitioning from K27me0 to K27me1 and so on.
# Note that adding me to higher methylation states is apparently harder thus the diffrence
Tp00=0.01
Tp01=0.99
Tp02=0.0
Tp03=0

Tp11=0.05
Tp12=0.95
Tp13=0.0

Tp22=0.1
Tp23=0.9

# The probability of mitosis happening at any move of PRC2
mitosis_prob=0.0001

# Probability of PRC2 does its function (depositin); 0-1, 1 means everytime if he sits on a nucleosome and if everything else (like other marks etc) look suitable, it does its job. Zero means it won't even if everything else looks okay.
depop=1

# Distribution of the genes, and other mark domains. It can be a list of several domains. Be careful with overlaps and going over the chromosome length
# genic_structure=list(c(100:300),c(800:1000))
# expression_structure=list(c(100:300))
# openchromatin=list(c(50:350),c(650:800))
# k36me2_structure=list(c(400:600))
# k36me3_structure=list(c(100:300))

genic_structure=list(c(980,1000))
expression_structure=list(c(980:1000))
openchromatin=list(c(980,1000))
k36me2_structure=list(c(980,1000))
k36me3_structure=list(c(980,1000))

# 5% of the genome but mostly found in expressed genes (open chromatin regions)
K27Mprobability=0.05*((chromlength)/(sum(lengths(openchromatin))))

# The value below shows how long K27M in scale of mitosis clock will stall the PRC2. In other words, if the probability of mitosis is 1/1000 and the value below is 1/10, that means K27M will stall PRC2 100 times more than its regular speed which is 1000 rounds per mitosis.
how_long_K27M_stalls=0.1

# How hard K36me2 affects the deposition of K27me3 (zero= completely prevents it, 1= not affecting at all)
howhardk36me2=0.1

# How hard K36me3 affects the deposition of K27me2 and K27me3 (zero= completely prevents it, 1= not affecting at all)
howhardk36me3=0

# Peakiness of K36me2 and K36me3 [0,0.1] (0 is a uniform block, 0.1 is a very sharp peak)
K36_slope=0.0005

# Slope of PRC2 falling off [0,0.1] (0 is a non-falling off, 0.1 is a very sharp slope)
prc2_slop=0.0025

# This factor indictes how big the previous nucleosome's status affects the deposition of K27me3. i.e. 2 means it makes it twice as likely etc.
neighborK27me3bonus=10

# Set this to zero if you want to avoid random walk movement for PRC2
randomwk=0

# The maximum number of steps for PRC2 in case of random walk before it falls off the chromatin
maximumsteps=10000

# The maximum size of each step in random walk 
stepsize=10





# Creating a directory to save the plots
newdir=paste(homedir,this_version,"-",gsub(":","_",gsub(" ","-",date())),"-","chrmlngth_",chromlength,"-pop_",populationSize,"-add01-12-23_0213-",Tp01,"-",Tp12,"-",Tp23,"-",Tp02,"-mitosis_prob_",mitosis_prob,"_Randomwalk_",randomwk,"_neighbor_",neighborK27me3bonus,"x","/",sep = "")
system(paste("mkdir -p ",newdir,sep = ""))
system(paste("mkdir -p ",newdir,"cell1/",sep = ""))
setwd(newdir)





#### Defining the probability of PRC2 falling off chromatin based on nucleosomeNumber and length of the chromosome
slope <- function(slope_length,slope_factor)
{
  # this changed: prc2hl()--> slope()
  hl<-(1-slope_factor)^slope_length
  return(hl)
}





#### Testing PRC2 half life slope
# Testin prc2hl
d<-c()
for (z in 1:1000)
{
  d<-c(d,slope(z, prc2_slop))
}
plot(d, main = paste("PRC2 Slope= ",prc2_slop,sep = ""))


#### Defining normal distribution with min and max and mean and sd
mysamp <- function(n, m, s, lwr, upr, nnorm) {
  samp <- rnorm(nnorm, m, s)
  samp <- samp[samp >= lwr & samp <= upr]
  if (length(samp) >= n) {
    return(sample(samp, n))
  }  
  stop(simpleError("Not enough values to sample from. Try increasing nnorm."))
}





#### Testing Normal distribution shape
##### Drawing a normal distribution with n values, standard deviation of sd with predefined minimum, maximum and average values from a pool of nnorm values.
# Testin prc2hl
d<-c()
n=1000
m=0.5
s=0.1
lwr=0
upr=1
nnorm=10000
d<-mysamp(n,m,s,lwr,upr,nnorm)
hist(d, main = paste("n= ",n," ,mean= ",m," ,sd= ",s, " , min=",lwr," ,max= ",upr," ,tries= ",nnorm))




#### Instructions for creating chromatin zero (with predefined Size and all the marks and structures set to 0)
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


################################## Functions for modifying chromatin zero ##################################

## Adding genic structure
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


## Adding uniform expression
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


#### Adding K36me2 and K36me3

### First, defining the type of peak (how sharp or flat)
shapeOfK36mark<-function(start,end,K36_slope)
{
   d1<-c()
   k36memarklength=(end-start)/2
   K36peakiness = K36_slope
   for (i in 1:k36memarklength)
   {
     d1<-c(d1,slope(i, K36peakiness))
   }
  d2<-rev(d1)
  d<-c(d2,d1)
  return(d)
}



### depositing k36me2 mark based on the peak model
depositk36me2<-function(chrom,k36me2_structure)
{
  for (k in k36me2_structure)
  {
    start=k[1]
    end=k[length(k36me2_structure[k])]
    i=1
    for (q in start:end)
    {
      if (!is.na(shapeOfK36mark(start,end,K36_slope)[i]) && shapeOfK36mark(start,end,K36_slope)[i] > runif(n=1,min=0,max=1))
      {
        chrom$K36me2[q]=1
      } else { chrom$K36me2[q]=0 }
      i=i+1
    }
  }
  chrom$K36me2[is.na(chrom$K36me2)] <- 0
  return(chrom)
}


### depositing k36me3 mark based on the peak model
depositk36me3<-function(chrom,k36me_structure)
{
  for (k in k36me3_structure)
  {
    start=k[1]
    end=k[length(k36me3_structure[k])]
    i=1
    for (q in start:end)
    {
      if (!is.na(shapeOfK36mark(start,end,K36_slope)[i]) && shapeOfK36mark(start,end,K36_slope)[i] > runif(n=1,min=0,max=1))
      {
        chrom$K36me3[q]=1
      } else { chrom$K36me3[q]=0 }
      i=i+1
    }
  }
  chrom$K36me3[is.na(chrom$K36me3)] <- 0
  return(chrom)
}


#### Adding K27M mutation
## below is the changed version based on the presence of H3.3 only in open chromatin i.e. expression_structure
addk27m <- function(chrom,openchromatin)
{
  for (k in openchromatin)
  {
    start=k[1]
    end=k[length(openchromatin[k])]
    for (q in start:end)
    {
      check=runif(n=1,min = 0,max=1)
      if (check < K27Mprobability)
      {
        chrom$K27M[q]=1
      }
    }
  }
  return(chrom)
}

#### PRC2 random-walking on the chromatin
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
    if (distance < 0)
    {
      distance=0
    }
    d<-c(d,distance)
    i=i+1
  }
  d<-subset(d,d[]>0 & d[] <=chromlength)
  if(length(d)<1) {d=1}
  return(d)
}

##### Testing Randomwalk function
test<-randomwalkme(1000,10,10000)
plot(test,type="l",ylab = paste("Distance; max dist = ",max(test),sep=""),xlab = paste("Number of Steps in the positive range; last step = ",length(test),sep=""))



#### Creating chromatin zero
chromatin<-createNakedChromatin(chromlength)
### Preparing the naked chromosome for K27me1,2,3 simulation by deposition of predefined K36me2,3 marks as well as defining genic regions, K27M mutation in open chromatin, and expression in different regions

chromatin[["chr"]]<-depositgene(chromatin[["chr"]],genic_structure)
chromatin[["chr"]]<-depositexpression(chromatin[["chr"]],expression_structure)
chromatin[["chr"]]<-depositk36me2(chromatin[["chr"]],k36me2_structure)
chromatin[["chr"]]<-depositk36me3(chromatin[["chr"]],k36me3_structure)
chromatin[["chr"]]<-addk27m(chromatin[["chr"]],openchromatin)

###################################### Depositing K27 me marks ######################################




