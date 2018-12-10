#### 1- Preparing Simulation ####
rm(list=ls())
par(mfrow=c(1,1))

# this_version = a prefix to the results' directory which contains information on the parameters and a version name. Just to keep tracks of versions.
this_version="A6.1-Constant-prc2location-H3.3-Only-in-OpenChr-Prob-K27M"

# homdir= any directory on your computer where you want to store the simulation results and plots, please make sure you have permission to write
homedir="~/Desktop/Desktop_June_2018/Histone_Mark_Simulation/Constant_clock/"

# Please do not change this, this is to avoid the scietific connotation in the results
options(scipen=999)

#### 2- Defining parameters ####
# Length of the chromatin (i.e. number of uncleosomes)
chromlength=1000

# For the purpose of plotting
snapshot_interval=1000


# Population size (number of chromatins)
populationSize=100

# Keeping record of populations' different marks throughout the life-time.
# This is a list that at each timer tick, a dataframe will be added to. These data frames contain the summation of the entire populations' 4 marks.
# This list will be used for plotting the average values of the entire population for chromosome plots at each interval (set by user) and the average plots for 4 different marks
# as well as where the simulation should stop based on the variance changes in the average values between consecutive timer ticks
population_history<-list()

# Maximum number of cycles (maximum) before the simulation stops. This will reach only if we don't get to a steady state before
life=50
var_threshold=0.00001

# Transition probabilities, the first number is the current state the second number os the next state, e.g. Tp01 means the probability of transitioning from K27me0 to K27me1 and so on.
# Note that adding me to higher methylation states is apparently harder thus the diffrence
Tp00=0.08
Tp01=0.92
Tp02=0
Tp03=0

Tp11=0.15
Tp12=0.85
Tp13=0

Tp22=0.44
Tp23=0.56

# The probability of mitosis happening at any move of PRC2
mitosis_prob=0.000001

# Distribution of the genes, and other mark domains. It can be a list of several domains. Be careful with overlaps and going over the chromosome length
genic_structure=list(c(981,1000))
expression_structure=list(c(981:1000))
openchromatin=list(c(0:99),c(400:499))
k36me2_structure=list(c(701,850))
k36me3_structure=list(c(701,850))

# 5% of the genome but mostly found in expressed genes (open chromatin regions)
K27Mprobability=0.5

# likelihood of H3.3 to appear on a nucleosome in open chromatin
h33probInOpenChr=0.8

#Length of open chromatin
openChrlength=sum(sapply(openchromatin, length))

# Probability of H3 outside of oepn chromatin
h33probInChr=(openChrlength-(openChrlength*h33probInOpenChr))/(chromlength-openChrlength)

# The value below shows how long K27M in scale of timer ticks prc2location will stall the PRC2. In other words, if the probability of mitosis is 1/1000 and the value below is 1/10, that means K27M will stall PRC2 100 times more than its regular speed which is 1000 rounds per mitosis.
stallingforK27M=50

# How much K36me2 affects the deposition of K27me3 (zero= completely prevents it, 1= not affecting at all) This value will be multiplied by the probability of the deposition, thus closer to 0 means lowering the probability more.
k36me2Effect=1

# How much K36me3 affects the deposition of K27me3 (zero= completely prevents it, 1= not affecting at all) This value will be multiplied by the probability of the deposition, thus closer to 0 means lowering the probability more.
k36me3Effect=1

# Peakiness of K36me2 and K36me3 [0,0.1] (0 is a uniform block, 0.1 is a very sharp peak)
K36_slope=0.0005

# This factor indicates how big the previous nucleosome's status affects the deposition of K27me3. i.e. 2 means it makes it twice as likely, 10 means times likely etc.
neighborK27me3bonus=10

# Prc2 movement model (0=sequenctial from nucleosome 1, 1=random walk, 2=random association)
prc2MovementModel=2

# Slope of PRC2 falling off [0,0.1] (0 is a non-falling off, 0.1 is a very sharp slope)
prc2_slop=0.0005
if(prc2MovementModel==2)
{prc2_slop=0.0}


# If the random walk model is set, this value sets the maximum number of steps for PRC2 before it falls off the chromatin
maximumsteps=10000

# The maximum size of each step in random walk model
stepsize=10





# Creating a directory to save the plots
#newdir=paste(homedir,this_version,"-",gsub(":","_",gsub(" ","-",date())),"-","chrmlngth_",chromlength,"-pop_",populationSize,"-add01-12-23_0213-",Tp01,"-",Tp12,"-",Tp23,"-",Tp02,"-mitosis_prob_",mitosis_prob,"_Randomwalk_",prc2MovementModel,"_neighbor_",neighborK27me3bonus,"x","/",sep = "")
#system(paste("mkdir -p ",newdir,sep = ""))
#system(paste("mkdir -p ",newdir,"cell1/",sep = ""))
#setwd(newdir)





#### Defining the probability of PRC2 falling off chromatin based on nucleosomeNumber and length of the chromosome
slope <- function(current_position,slope_factor)
{
  # this changed: prc2hl()--> slope()
  hl<-(1-slope_factor)^current_position
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




#### 3- Instructions for creating chromatin zero (with predefined Size and all the marks and structures set to 0) ####
createNakedChromatin <- function(chromSize)
{
  chromatin<-data.frame(integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),double(chromSize))
  colnames(chromatin)<-c("N","S","me0","me1","me2","me3","K36me2","K36me3","genic","expression","h33","K27M","prc2_falling_chance")
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
    chromatin$h33[s]=0
    chromatin$K27M=0
    chromatin$prc2_falling_chance=0
  }
  results<-list("add1"=0,"add2"=0,"add3"=0,"add0"=0,"rando"=0,"prc2_attached"=1,"chr"=chromatin)
  return(results)
}


#### 4- Functions for modifying chromatin zero ####

## 4-1: Adding genic structure
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


## 4-2: Adding uniform expression
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


#### 4-3: Adding K36me2 and K36me3

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



### 4-4: depositing k36me2 mark based on the peak model
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


### 4-5: depositing k36me3 mark based on the peak model
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


####4-6-0: Distribution of H3.3 in Open chromatin
## around 10% of the nucleosomes have H3.3 which are the histones that get mutated (with ~50% likelihood)
## These histones are mainly concentrated in open chromatin (likelihoods can change). 
## Thiss function defines the distribution of H3.3 (in open chromatin as long as rest of the genome with lower likelihood)



addH33 <- function(chrom,openchromatin,h33probInChr,h33probInOpenChr)
{
  for (k in 1:nrow(chrom))
  {
    j=0
    for (i in 1:length(openchromatin))
    { 
      if (k %in% openchromatin[[i]])
      {
      j=1
      } 
    }
    
    if(j==1)
    {
    
      if (runif(n=1,min = 0,max=1) <= h33probInOpenChr)
      {
      chrom$h33[k]=1
      }
    }
    else if(j==0)
    {
    if (runif(n=1,min = 0,max=1) <= h33probInChr)
      {
      chrom$h33[k]=1
      }
    }
  }
return(chrom)
}


#### 4-6: Adding K27M mutation
## below is the changed version based on the presence of H3.3 only in open chromatin i.e. expression_structure
## around 10% of H3s are H3.3 and around half of them are mutated
addk27m <- function(chrom)
{
  for (k in 1:nrow(chrom))
  {
    check=runif(n=1,min = 0,max=1)
    if (check < K27Mprobability && chrom$h33[k]==1)
    {
      chrom$K27M[k]=1
    }
  }
  return(chrom)
}

#### 4-7: PRC2 random-walking on the chromatin
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



## If randomwalk is chosen for PRC2 move, the "PRC2moveVector" variable which is an vector made of randomwalked variables will be the order of PRC2 moves, otherwise, if randomwalk is set to zero, it'll be a vector of 0 to size of the chromatin.
if (prc2MovementModel==0)
{
  PRC2moveVector=c(1:chromlength)
} else if(prc2MovementModel==1)
{
  PRC2moveVector=randomwalkme(chromlength,stepsize,maximumsteps)
} else if (prc2MovementModel==2)
{
  PRC2moveVector=sample(c(1:chromlength),size = life*snapshot_interval,replace = TRUE)
}

#### Creating chromatin zero
chromatin<-createNakedChromatin(chromlength)
### Preparing the naked chromosome for K27me1,2,3 simulation by deposition of predefined K36me2,3 marks as well as defining genic regions, K27M mutation in open chromatin, and expression in different regions

chromatin[["chr"]]<-depositgene(chromatin[["chr"]],genic_structure)
chromatin[["chr"]]<-depositexpression(chromatin[["chr"]],expression_structure)
chromatin[["chr"]]<-depositk36me2(chromatin[["chr"]],k36me2_structure)
chromatin[["chr"]]<-depositk36me3(chromatin[["chr"]],k36me3_structure)



#### 4-8: Defining mitosis function
mitosis <- function (chromatin){
  for (i in 1:nrow(chromatin)) 
  {
    check=runif(n=1,min = 0,max=1)
    if (check > 0.5)
    {
      chromatin$S[i]=0
      chromatin$me0[i]=1
      chromatin$me1[i]=0
      chromatin$me2[i]=0
      chromatin$me3[i]=0
    }
  }
  return (chromatin)
}
# Testing mitosis
# test<-chromatin[["chr"]]
# test$S=1
# test$me1=1
# test$me0=0
# test2<-mitosis(test)
# View(test)
# View(test2)
# 



#### 4-9: Depositing K27 me marks on each nucleosome
deposit<-function(chromatin,nucleosomeNumber)
{
  
  #choosing a random number between 0 and 1 for testing with the probability of each scenario
  rando=runif(n=1,min=0,max=1) 
  
  k36me2penalty=k36me3penalty=0
  if(chromatin[["chr"]]$K36me2[nucleosomeNumber]==1)
  {k36me2penalty=k36me2Effect}
  if(chromatin[["chr"]]$K36me3[nucleosomeNumber]==1)
  {k36me3penalty=k36me3Effect}
  
  #When PRC2 makes it to each nucleosome, checks the current state (S). If the current state shows no K27me (1,2, or 3)   
  if (chromatin[["chr"]]$S[nucleosomeNumber] == 0)
  {
    #If the random generate number falls between zero and Tp00 (times the chance of falling off) nothing changes
    if (rando < Tp00)
    {
      chromatin[["chr"]]$S[nucleosomeNumber]=0
      chromatin[["chr"]]$me0[nucleosomeNumber]=1
      chromatin[["chr"]]$me1[nucleosomeNumber]=0
      chromatin[["chr"]]$me2[nucleosomeNumber]=0
      chromatin[["chr"]]$me3[nucleosomeNumber]=0
    } 
    # If it falls between zero and Tp01(*fall off chance), changes the S from 0 to 1
    else if(rando < Tp01) {
      chromatin[["chr"]]$S[nucleosomeNumber]=1
      chromatin[["chr"]]$me0[nucleosomeNumber]=0
      chromatin[["chr"]]$me1[nucleosomeNumber]=1
      chromatin[["chr"]]$me2[nucleosomeNumber]=0
      chromatin[["chr"]]$me3[nucleosomeNumber]=0
    } 
    # If it falls between zero and Tp02(*fall off chance), changes the S from 0 to 2
    else  if(rando < Tp02) {
      chromatin[["chr"]]$S[nucleosomeNumber]=2
      chromatin[["chr"]]$me0[nucleosomeNumber]=0
      chromatin[["chr"]]$me1[nucleosomeNumber]=0
      chromatin[["chr"]]$me2[nucleosomeNumber]=1
      chromatin[["chr"]]$me3[nucleosomeNumber]=0
    } 
    # If it falls between zero and Tp03(*fall off chance), changes the S from 0 to 3
    else  if(rando < Tp03) {
      chromatin[["chr"]]$S[nucleosomeNumber]=3
      chromatin[["chr"]]$me0[nucleosomeNumber]=0
      chromatin[["chr"]]$me1[nucleosomeNumber]=0
      chromatin[["chr"]]$me2[nucleosomeNumber]=0
      chromatin[["chr"]]$me3[nucleosomeNumber]=1
    }
  }
  ## Same story as above if the current state is 1
  else if (chromatin[["chr"]]$S[nucleosomeNumber] == 1)  {
    if(rando < Tp11)
    {
      chromatin[["chr"]]$S[nucleosomeNumber]=1
      chromatin[["chr"]]$me0[nucleosomeNumber]=0
      chromatin[["chr"]]$me1[nucleosomeNumber]=1
      chromatin[["chr"]]$me2[nucleosomeNumber]=0
      chromatin[["chr"]]$me3[nucleosomeNumber]=0
    } else if(rando < Tp12) {
      chromatin[["chr"]]$S[nucleosomeNumber]=2
      chromatin[["chr"]]$me0[nucleosomeNumber]=0
      chromatin[["chr"]]$me1[nucleosomeNumber]=0
      chromatin[["chr"]]$me2[nucleosomeNumber]=1
      chromatin[["chr"]]$me3[nucleosomeNumber]=0
    } else  if(rando < Tp13) {
      chromatin[["chr"]]$S[nucleosomeNumber]=3
      chromatin[["chr"]]$me0[nucleosomeNumber]=0
      chromatin[["chr"]]$me1[nucleosomeNumber]=0
      chromatin[["chr"]]$me2[nucleosomeNumber]=0
      chromatin[["chr"]]$me3[nucleosomeNumber]=1
    } 
  } 
  ## Same story as above if the current state is 2, except here we have the bonus for the previous nucleosome having K27me3 mark
  else if (chromatin[["chr"]]$S[nucleosomeNumber] == 2)  {
    if(rando < Tp22)
    {
      chromatin[["chr"]]$S[nucleosomeNumber]=2
      chromatin[["chr"]]$me0[nucleosomeNumber]=0
      chromatin[["chr"]]$me1[nucleosomeNumber]=0
      chromatin[["chr"]]$me2[nucleosomeNumber]=1
      chromatin[["chr"]]$me3[nucleosomeNumber]=0
    } else if(rando < Tp23)  {
      chromatin[["chr"]]$S[nucleosomeNumber]=3
      chromatin[["chr"]]$me0[nucleosomeNumber]=0
      chromatin[["chr"]]$me1[nucleosomeNumber]=0
      chromatin[["chr"]]$me2[nucleosomeNumber]=0
      chromatin[["chr"]]$me3[nucleosomeNumber]=1
    }
  }
  #  print (paste("rando: ",rando,sep=""))
  return(chromatin)
}




## Creating the initial population 
population<-list()
#currentpop<-list()  
currentpop<-list("prc2location"=1,"chr"=createNakedChromatin(chromlength)[["chr"]])  
currentpop[["chr"]]$me0=0


me0avg<-c()
me1avg<-c()
me2avg<-c()
me3avg<-c()
k36me2avg<-c()
k36me3avg<-c()
k27mavg<-c()


# Deciding where the K27M mutated nucleosomes will appear on the chromatin
thispopulationk27mdistribution<-chromatin[["chr"]]$K27M


#### 5- RUNNING POPULATION SIMULATION  ####
# Setting up the population abd the simulation parameters

stall=0
prc2location=1
timer=1
mitosis_count=0
prc2_round_counter=0

#### Create the naked chromatin for the population simulation
for (p in 1:populationSize)
{
  population[[p]]<-createNakedChromatin(chromlength)
  population[[p]][["chr"]]<-depositk36me2(population[[p]][["chr"]],k36me2_structure)
  population[[p]][["chr"]]<-depositk36me3(population[[p]][["chr"]],k36me3_structure)
  population[[p]][["chr"]]<-depositgene(population[[p]][["chr"]],genic_structure)
  population[[p]][["chr"]]<-depositexpression(population[[p]][["chr"]],expression_structure)
  population[[p]][["chr"]]<-addH33(population[[p]][["chr"]],openchromatin,h33probInChr,h33probInOpenChr)
  population[[p]][["chr"]]<-addk27m(population[[p]][["chr"]])
  
}  

## Testing the distribution of H3.3 and K27M
populationTotal<-data.frame(integer(chromlength),integer(chromlength))
colnames(populationTotal)<-c("h33","k27m")
for (i in 1:populationSize)
{
  populationTotal$h33<-populationTotal$h33+population[[i]][["chr"]]$h33
  populationTotal$k27m<-populationTotal$k27m+population[[i]][["chr"]]$K27M
}
## End of test


## Initializing the plotting events
population_Marks_sum<-data.frame(me0=integer(chromlength),me1=integer(chromlength),me2=integer(chromlength),me3=integer(chromlength))
population_Marks_sum$me0=0
population_Marks_sum$me1=0
population_Marks_sum$me2=0
population_Marks_sum$me3=0


average_values0<-c()
average_values1<-c()
average_values2<-c()
average_values3<-c()



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
# Begginning of timer's life
prc2pathindex=1
breaking=0
while (timer <= life*snapshot_interval)
{
  
  if(prc2location%%chromlength==0 || prc2location==1 && prc2MovementModel<2){print (paste("Timer: ",timer," > prc2location: ",prc2location,sep = ""))}
  
  ### Check if the prc2 has reached the end of chromosome and if so, reset it to first nucleosome  
  if ((prc2location-1) %% chromlength == 0 && prc2MovementModel<2) 
  {
    if(prc2MovementModel==1){
      prc2path<-randomwalkme(chromlength,stepsize,life)
    } else if (prc2MovementModel==0) {prc2path<-c(1:chromlength)}
    
    prc2pathindex=1
    prc2location=1
    prc2_round_counter=prc2_round_counter+1
    
  }
  ###NEW
  if (prc2MovementModel==2) 
  {
  prc2path<-PRC2moveVector
  prc2_round_counter=prc2_round_counter+1
  }
  ### Checking out if PRC2 is still attached
  
  
  #### If the simulation has been going on for at least one round of time, then start the deposition
  if (prc2location > 0 && timer > 0) ### Start of the actual deposition 
  {
    prc2_threshold<-slope((prc2location),prc2_slop)
    
    for (p in 1:populationSize) ## We go position by position across all individuals in the population
    {
      if(population[[p]][["chr"]]$prc2_falling_chance[prc2location] <= prc2_threshold || prc2location==1)
      {
        population[[p]]$prc2_attached=1
      }
      else {
        population[[p]]$prc2_attached=0
        population[[p]][["chr"]]$prc2_falling_chance[prc2location]=1
      }
      
      ####### Now if PRC2 is attached:
      if (population[[p]]$prc2_attached==1)
      {
        population[[p]]<-deposit(population[[p]],prc2location)
        population[[p]][["chr"]]$prc2_falling_chance[prc2location] <- runif(1,min = 0,max=1)
      }
      mitosis_chance=runif(n=1,min=0,max=1)
      if (mitosis_chance < mitosis_prob )
      {
        population[[p]][["chr"]]<-mitosis(population[[p]][["chr"]])
        mitosis_count=mitosis_count+1
        mitosis_chance=0
        print(paste("individual: ",p," > timer: ",timer," > prc2 attachment: ",prc2_round_counter,sep=""))
      }
      if(prc2location%%snapshot_interval==0 || timer==1)
      {      
        population_Marks_sum$me0<-population_Marks_sum$me0+population[[p]][["chr"]]$me0
        population_Marks_sum$me1<-population_Marks_sum$me1+population[[p]][["chr"]]$me1
        population_Marks_sum$me2<-population_Marks_sum$me2+population[[p]][["chr"]]$me2
        population_Marks_sum$me3<-population_Marks_sum$me3+population[[p]][["chr"]]$me3
        if(p==populationSize)
        {
          print (paste("prc2location=",prc2location,"~",chromlength," > p=",p,"~",populationSize))
          average_values0<-c(average_values0,sum(population_Marks_sum$me0)/(populationSize*chromlength))
          average_values1<-c(average_values1,sum(population_Marks_sum$me1)/(populationSize*chromlength))
          average_values2<-c(average_values2,sum(population_Marks_sum$me2)/(populationSize*chromlength))
          average_values3<-c(average_values3,sum(population_Marks_sum$me3)/(populationSize*chromlength))
          population_Marks_sum$me0=0
          population_Marks_sum$me1=0
          population_Marks_sum$me2=0
          population_Marks_sum$me3=0
          
          # Check if K27me3 has reached plateau (variance of the last 5 average values < 0.001)
          if(length(average_values3) > 5)
          {
            var3<-var(average_values3[(length(average_values3)-5):length(average_values3)])
            var2<-var(average_values2[(length(average_values2)-5):length(average_values2)])
            var1<-var(average_values1[(length(average_values1)-5):length(average_values1)])
            var0<-var(average_values0[(length(average_values0)-5):length(average_values0)])
            
            if((var3 < var_threshold && var2 < var_threshold && var1 < var_threshold && var0 < var_threshold) || timer==life)
            {
              #print (paste("breaking time:",timer,sep=""))
              breaking=1;
            }
          }
        }
      }
    }
    
    
    # Stall starts at zero. At the end of each timer round, stall is checked, if zero, 
    # that means there was no K27M at the position where prc2 is landed at the time=timer and it can move on, 
    # if stall was larger than zero, then it means there was a K27M and the prc2location cannot change 
    # unless the stall surpasses a threshold (for example 50 which equals to timer cycle going on for 50 rounds or 
    # in other words the prc2 has been stalled for 50 round)
    if (population[[p]][["chr"]]$K27M[(prc2location+1)]==1 && prc2location<1000)
    {  
      stall=stall+1
    }
    
  } ### End of actions if prc2 is located farther than the start point and timer > 0
  
  
  timer=timer+1
  if(stall > stallingforK27M || stall ==0)
  {
    stall=0
    prc2pathindex=prc2pathindex+1
    if(prc2pathindex > chromlength){prc2pathindex=1}
    prc2location=prc2path[prc2pathindex]
    
    ## Changing the prc2 attachment chance only if the prc2 is moving
    population[[p]][["chr"]]$prc2_falling_chance<-runif(1,min = 0,max=1)
  } else {prc2location=prc2location}
  
  
  if(breaking==1)
  {
    print (paste("breaking time:",timer,sep=""))
    break; 
  }
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
} ### End of timer's life


print(paste("Number of mitosis: ",mitosis_count,sep=""))
print(paste("Prc2 Rounds: ",prc2_round_counter,sep=""))

#### 6- Plotting the average values ####

# par(mfrow=c(2,2))
# yaxis_marks<-seq(from = 0, to = timer, by = snapshot_interval)
# plot(average_values0,type = "l",xaxt="n",xlab="Timer",ylab="me0",ylim=c(0,1.2))
# axis(1, at=(yaxis_marks/snapshot_interval), labels=yaxis_marks) 
# plot(average_values1,type = "l",xaxt="n",xlab="Timer",ylab="me1")
# axis(1, at=(yaxis_marks/snapshot_interval), labels=yaxis_marks) 
# plot(average_values2,type = "l",xaxt="n",xlab="Timer",ylab="me2")
# axis(1, at=(yaxis_marks/snapshot_interval), labels=yaxis_marks) 
# plot(average_values3,type = "l",xaxt="n",xlab="Timer",ylab="me3")
# axis(1, at=(yaxis_marks/snapshot_interval), labels=yaxis_marks) 
# 

#### 7- Plotting the end chromatin situation ####

plot_df<-data.frame(me0=integer(chromlength),me1=integer(chromlength),me2=integer(chromlength),me3=integer(chromlength))
plot_df$me0=0
plot_df$me1=0
plot_df$me2=0
plot_df$me3=0
plot_df$k36me2=0
plot_df$k36me3=0


for (i in 1:populationSize)
{
  plot_df$me0<-plot_df$me0+population[[i]][["chr"]]$me0
  plot_df$me1<-plot_df$me1+population[[i]][["chr"]]$me1
  plot_df$me2<-plot_df$me2+population[[i]][["chr"]]$me2
  plot_df$me3<-plot_df$me3+population[[i]][["chr"]]$me3
  plot_df$k36me2<-plot_df$k36me2+population[[i]][["chr"]]$K36me2
  plot_df$k36me3<-plot_df$k36me3+population[[i]][["chr"]]$K36me3
}

plot_df<-plot_df/populationSize


par(mfrow=c(4,1))
plot(populationTotal$h33/populationSize,ylim=c(0,1),type="h", ylab = "h3.3")
plot(populationTotal$k27m/populationSize,ylim=c(0,1),type="h",ylab="k27m")
plot(plot_df$k36me2,ylim=c(0,1),ylab = "k36me2",type = "h")
plot(plot_df$k36me3,ylim=c(0,1),ylab = "k36me3",type = "h")


par(mfrow=c(4,1))
plot(plot_df$me0,ylim=c(0,1),ylab = "me0",type = "h")
plot(plot_df$me1,ylim=c(0,1),ylab = "me1",type = "h")
plot(plot_df$me2,ylim=c(0,1),ylab = "me2",type = "h")
plot(plot_df$me3,ylim=c(0,1),ylab = "me3",type = "h")





### ToDo: To put back the effectors on the probability of depositing marks such as K36me2,3 
### ToDo: To put back the effect of neighboring nucleosomes (last step)
