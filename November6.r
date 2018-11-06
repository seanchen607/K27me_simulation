## 1- Preparing Simulation
rm(list=ls())
par(mfrow=c(1,1))

# this_version = a prefix to the results' directory which contains information on the parameters and a version name. Just to keep tracks of versions.
this_version="A6.1-Constant-prc2location-H3.3-Only-in-OpenChr-Prob-K27M"

# homdir= any directory on your computer where you want to store the simulation results and plots, please make sure you have permission to write
homedir="~/Desktop/Desktop_June_2018/Histone_Mark_Simulation/Constant_clock/"

# Please do not change this, this is to avoid the scietific connotation in the results
options(scipen=999)


## 2- Defining parameters
# Length of the chromatin (i.e. number of uncleosomes)
chromlength=1000

# Population size (number of chromatins)
populationSize=20

# Number of cycles (maximum) before the simulation stops. This will reach only if we don't get to a steady state before
life=50000

# Transition probabilities, the first number is the current state the second number os the next state, e.g. Tp01 means the probability of transitioning from K27me0 to K27me1 and so on.
# Note that adding me to higher methylation states is apparently harder thus the diffrence
Tp00=0.1
Tp01=0.9
Tp02=0.0
Tp03=0

Tp11=0.3
Tp12=0.7
Tp13=0.0

Tp22=0.4
Tp23=0.6

# The probability of mitosis happening at any move of PRC2
mitosis_prob=0.00001

# Probability of PRC2 does its function (depositin); 0-1, 1 means everytime if he sits on a nucleosome and if everything else (like other marks etc) look suitable, it does its job. Zero means it won't even if everything else looks okay.
#depop=1

# Distribution of the genes, and other mark domains. It can be a list of several domains. Be careful with overlaps and going over the chromosome length
# genic_structure=list(c(100:300),c(800:1000))
# expression_structure=list(c(100:300))
# openchromatin=list(c(50:350),c(650:800))
# k36me2_structure=list(c(400:600))
# k36me3_structure=list(c(100:300))

genic_structure=list(c(980,1000))
expression_structure=list(c(980:1000))
openchromatin=list(c(400:460))
k36me2_structure=list(c(980,1000))
k36me3_structure=list(c(980,1000))

# 5% of the genome but mostly found in expressed genes (open chromatin regions)
K27Mprobability=0.05
#*((chromlength)/(sum(lengths(openchromatin))))

# The value below shows how long K27M in scale of mitosis prc2location will stall the PRC2. In other words, if the probability of mitosis is 1/1000 and the value below is 1/10, that means K27M will stall PRC2 100 times more than its regular speed which is 1000 rounds per mitosis.
stallingforK27M=5

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
#system(paste("mkdir -p ",newdir,sep = ""))
#system(paste("mkdir -p ",newdir,"cell1/",sep = ""))
setwd(newdir)





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




#### Instructions for creating chromatin zero (with predefined Size and all the marks and structures set to 0)
createNakedChromatin <- function(chromSize)
{
  chromatin<-data.frame(integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),integer(chromSize),double(chromSize))
  colnames(chromatin)<-c("N","S","me0","me1","me2","me3","K36me2","K36me3","genic","expression","K27M","prc2_falling_chance")
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
    chromatin$prc2_falling_chance=0
  }
  results<-list("add1"=0,"add2"=0,"add3"=0,"add0"=0,"rando"=0,"prc2_attached"=1,"chr"=chromatin)
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



## If randomwalk is chosen for PRC2 move, the "PRC2moveVector" variable which is an vector made of randomwalked variables will be the order of PRC2 moves, otherwise, if randomwalk is set to zero, it'll be a vector of 0 to size of the chromatin.
if (randomwk==0)
{
  PRC2moveVector=c(1:chromlength)
} else 
{
  PRC2moveVector=randomwalkme(chromlength,stepsize,maximumsteps)
}


#### Creating chromatin zero
chromatin<-createNakedChromatin(chromlength)
### Preparing the naked chromosome for K27me1,2,3 simulation by deposition of predefined K36me2,3 marks as well as defining genic regions, K27M mutation in open chromatin, and expression in different regions

chromatin[["chr"]]<-depositgene(chromatin[["chr"]],genic_structure)
chromatin[["chr"]]<-depositexpression(chromatin[["chr"]],expression_structure)
chromatin[["chr"]]<-depositk36me2(chromatin[["chr"]],k36me2_structure)
chromatin[["chr"]]<-depositk36me3(chromatin[["chr"]],k36me3_structure)
chromatin[["chr"]]<-addk27m(chromatin[["chr"]],openchromatin)






#### Defining mitosis function
mitosis <- function (chromatin){
  mitochromatin<-chromatin
  for (i in 1:nrow(chromatin[["chr"]])) 
  {
    check=runif(n=1,min = 0,max=1)
    if (check > 0.5)
    {
      mitochromatin[["chr"]]$S[i]=0
      mitochromatin[["chr"]]$me0[i]=1
      mitochromatin[["chr"]]$me1[i]=0
      mitochromatin[["chr"]]$me2[i]=0
      mitochromatin[["chr"]]$me3[i]=0
    }
  }
  return (mitochromatin)
}
## Testing mitosis
#test<-chromatin
#test[["chr"]]$S=1
#test[["chr"]]$me1=1
#test[["chr"]]$me0=0
#test2<-mitosis(test)
#View(test[["chr"]])
#View(test2[["chr"]])










###################################### Depositing K27 me marks on each nucleosome ######################################

deposit<-function(chromatin,nucleosomeNumber)
{
  
  #choosing a random number between 0 and 1 for testing with the probability of each scenario
  rando=runif(n=1,min=0,max=1)  
  
  
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


#############################################################
############### RUNNING POPULATION SIMULATION  ##############
#############################################################

stall=0
totalmitosis=0
prc2location=0
timer=0


## Below for plotting the average total marks

average0total=c()
average1total=c()
average2total=c()
average3total=c()
total_me0=total_me1=total_me2=total_me3=0


mitosis_count=0
prc2_round_counter=0


while (timer <= life)
{
  #mitosis_count=0
#  print (paste("Timer: ",timer,sep = ""))
  
  if ((prc2location-1) %% chromlength == 0 && prc2location > 1)
  {
    prc2location=1
    prc2_round_counter=prc2_round_counter+1
#    print("#################")
  }
  
  ### Checking out if PRC2 is still attached
  
  
  
  
  #### If the simulation just started, create the naked chromatin for the population
  if (prc2location==0 && timer==0) 
  {
    for (p in 1:populationSize)
    {
      population[[p]]<-createNakedChromatin(chromlength)
      population[[p]][["chr"]]<-depositk36me2(population[[p]][["chr"]],k36me2_structure)
      population[[p]][["chr"]]<-depositk36me3(population[[p]][["chr"]],k36me3_structure)
      population[[p]][["chr"]]<-depositgene(population[[p]][["chr"]],genic_structure)
      population[[p]][["chr"]]<-depositexpression(population[[p]][["chr"]],expression_structure)
      population[[p]][["chr"]]$K27M<-thispopulationk27mdistribution
    }  
  } 
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

      if(prc2location %% 100 ==0 && p==5){#print(paste("Round: ",prc2_round_counter," >> P= ",p," > prc2 location: ",prc2location," > prc2_falling_chance: ",population[[p]][["chr"]]$prc2_falling_chance," > prc2 attached: ",population[[p]]$prc2_attached,sep = ""))}
      #      print(paste("P: ",p," >>> prc2 location: ",prc2location," >>> prc2 attached: ",population[[p]]$prc2_attached,sep=""))
      ## A random number between 0 and 1 to be compared to the value calculated based on the prc2location or timer and the prc2slop (i.e. how fast the halflife of prc2 declines)
####### Now if PRC2 is attached:
      if (population[[p]]$prc2_attached==1)
      {
        population[[p]]<-deposit(population[[p]],prc2location)
        if(prc2location %% 50 ==0)
      #  print(paste("threshold: ",prc2_threshold," > chance: ",population[[p]][["chr"]]$prc2_falling_chance[prc2location],sep = ""))
        population[[p]][["chr"]]$prc2_falling_chance[prc2location] <- runif(1,min = 0,max=1)
        }
      mitosis_chance=runif(n=1,min=0,max=1)
      if (mitosis_chance < mitosis_prob )
      {
        population[[p]]<-mitosis(population[[p]])
        mitosis_count=mitosis_count+1
      }
      mitosis_chance=0
    }
    
    if (population[[p]][["chr"]]$K27M[(prc2location+1)]==1 && prc2location<1000)
    {  
      #      print(paste(prc2location,stall,sep = "-"))
      # Stall starts at zero. At the end of each timer cycle, stall is checked, if zero, that means there was no K27M and it can move on, if stall was larger than zero, then it means there was a K27M and the prc2location cannot change unless the stall surpasses a threshold (for example 50 which equals to timer cycle going on for 50 rounds or in other words the prc2 has been stalled for 50 round)
      stall=stall+1
    }
  }
  
  
  
  
  #if (prc2location == chromlength || timer == life || timer==1 || timer==0){
  
  
  average0=c()
  average1=c()
  average2=c()
  average3=c()
  
  for (p in 1:populationSize)
  {
    total_me0_timer=sum(population[[p]][["chr"]]$me0)
    total_me1_timer=sum(population[[p]][["chr"]]$me1)
    total_me2_timer=sum(population[[p]][["chr"]]$me2)
    total_me3_timer=sum(population[[p]][["chr"]]$me3)
  }
  
  if (prc2location%%1000==0){
    ### Average plots
    for (nucleosome in 1:chromlength)
    {
      summation0=0
      summation1=0
      summation2=0
      summation3=0
      
      for (p in 1:populationSize)
      {
        summation0<-summation0+population[[p]][["chr"]]$me0[nucleosome] 
        summation1<-summation1+population[[p]][["chr"]]$me1[nucleosome] 
        summation2<-summation2+population[[p]][["chr"]]$me2[nucleosome] 
        summation3<-summation3+population[[p]][["chr"]]$me3[nucleosome] 
      }
      
      average0<-c(average0,(summation0/populationSize))
      average1<-c(average1,(summation1/populationSize))
      average2<-c(average2,summation2/populationSize)
      average3<-c(average3,summation3/populationSize)
    }
    
    
    
    # average0total<-c(average0total,(total_me0_timer/chromlength))
    # average1total<-c(average1total,(total_me1_timer/chromlength))
    # average2total<-c(average2total,(total_me2_timer/chromlength))
    # average3total<-c(average3total,(total_me3_timer/chromlength))
    
    
    
    #### Snapshot ####
    #if (timer%%1000==0){
    par(mfrow=c(5,1))
    plot(average0, ylim = c(0,1),type = "h", main = paste("Snapshot at time: ",timer," PRC2 Location: ",prc2location,sep=""), ylab="me0")
    plot(average1, ylim = c(0,1),type = "h", ylab = "me1")
    plot(average2, ylim = c(0,1),type = "h", ylab = "me2")
    plot(average3, ylim = c(0,1),type = "h", ylab = "me3")
    plot(population[[1]][["chr"]]$K27M, ylim = c(0,1),type = "h", ylab = "K27M")
  }
  ###NEW
  average0total<-c(average0total,(total_me0_timer/chromlength))
  average1total<-c(average1total,(total_me1_timer/chromlength))
  average2total<-c(average2total,(total_me2_timer/chromlength))
  average3total<-c(average3total,(total_me3_timer/chromlength))
  
  
  timer=timer+1
  if(stall > stallingforK27M || stall ==0)
  {
    stall=0
    prc2location=prc2location+1
    ## Changing the prc2 attachment chance only if the prc2 is moving
    population[[p]][["chr"]]$prc2_falling_chance<-runif(1,min = 0,max=1)
  } else {prc2location=prc2location}
}
print(paste("Number of mitosis: ",mitosis_count,sep=""))
print(paste("Prc2 Rounds: ",prc2_round_counter,sep=""))

par(mfrow=c(2,2))
plot(average0total, type = "l", main=(paste("mitosis prob= ",mitosis_prob,sep = "")))
plot(average1total, type = "l")
plot(average2total, type = "l")
plot(average3total, type = "l")


##### To add a reset for mitosis probability checker and a vector that holds 0 or 1 for if the prc2 is attached for each nucleosome
### ToDo: To implement effectors on the probability of depositing marks such as K36me2,3 
### ToDo: To implement the effect of neighboring nucleosomes (last step)

