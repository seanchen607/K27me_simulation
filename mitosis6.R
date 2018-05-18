setwd("~/Desktop/Histone_Mark_Simulation/Images/new/")
options(scipen=999)
mysamp <- function(n, m, s, lwr, upr, nnorm) {
  samp <- rnorm(nnorm, m, s)
  samp <- samp[samp >= lwr & samp <= upr]
  if (length(samp) >= n) {
    return(sample(samp, n))
  }  
  stop(simpleError("Not enough values to sample from. Try increasing nnorm."))
}

nnuc=500
### prc2 initiation
prc2=0
prc2rounds=20
#prc2rounds=prc2rounds
burden12=0.00
burden23=0.55
# population initiation
finalz<-data.frame(integer(nnuc),integer(nnuc),integer(nnuc),integer(nnuc),integer(nnuc))
colnames(finalz)<-c("N","me0","me1","me2","me3")
for (n in 1:nnuc)
{
  finalz$N[n]=n
  finalz$me0[n]=0
  finalz$me1[n]=0
  finalz$me2[n]=0
  finalz$me3[n]=0
}


passages<- data.frame(matrix(nrow = nnuc,ncol = (3*prc2rounds)))
passcolnames<-c()
#passcolnames<-c("N")
for (i in 1:prc2rounds){
  #  passcolnames<-c(passcolnames,paste(i,"me0",sep="_"))
  passcolnames<-c(passcolnames,paste(i,"me1",sep="_"))
  passcolnames<-c(passcolnames,paste(i,"me2",sep="_"))
  passcolnames<-c(passcolnames,paste(i,"me3",sep="_"))
}
colnames(passages)<-passcolnames
for (n in 1:nnuc)
{
  passages[n,]=0
}

initial_cell_count=50
final_cell_count=initial_cell_count

for (cell in 1:initial_cell_count)  {
  
  go_to_next_cell=0
  
  ### Creating Naked Chromatin
  chromatin<-data.frame(integer(nnuc),integer(nnuc),integer(nnuc),integer(nnuc),integer(nnuc),integer(nnuc))
  colnames(chromatin)<-c("N","S","me0","me1","me2","me3")
  
  
  for (n in 1:nnuc)
  {
    chromatin$N[n]=n
    chromatin$S[n]=0
    chromatin$me0[n]=1
    chromatin$me1[n]=0
    chromatin$me2[n]=0
    chromatin$me3[n]=0
  }
  if (cell ==1)
  {
    #png(paste("~/Desktop/Histone_Mark_Simulation/Images/000-3marks-average.png",sep=""),width = 800,height = 600)
    
    par(mfrow=c(4,1))
    plot(chromatin$me0,type = "h", ylim = c(0,1), main = paste("State 0 : K27me0","  -- PRC2 passage: ",prc2,sep=""),ylab="")
    plot(chromatin$me1,type = "h", ylim = c(0,1), main = paste("State 0 : K27me1","  -- PRC2 passage: ",prc2,sep=""),ylab="")
    plot(chromatin$me2,type = "h", ylim = c(0,1), main = paste("State 0 : K27me2","  -- PRC2 passage: ",prc2,sep=""),ylab="")
    plot(chromatin$me3,type = "h", ylim = c(0,1), main = paste("State 0 : K27me3","  -- PRC2 passage: ",prc2,sep=""),ylab="")
    #dev.off()
  }
  
  padd0=c()
  padd1=c()
  padd2=c()
  padd3=c()
  
#  prc2rounds=prc2rounds 
  print(prc2rounds)  
  
  ### Runs of PRC2
  for (prc2 in 1:prc2rounds) {
    ### Depositing marks
    for (n in 1:nnuc)
    {
      
      padd1[n]=mysamp(n=1, m=0.69, s=0.1, lwr=0, upr=1, nnorm=100)
      padd2[n]=mysamp(n=1, m=0.01, s=0.1, lwr=0, upr=1, nnorm=100)
      padd3[n]=mysamp(n=1, m=0.00, s=0.1, lwr=0, upr=1, nnorm=100)
      padd0[n]=1-(sum(padd1[n],padd2[n],padd3[n]))
      
      
      # PRC2 halflife
      if (0.3-(n/nnuc) < 0 ){lower=0} else {lower=0.3-(n/nnuc)}
      if (1.3-(n/nnuc) > 1){upper=1} else {upper=1.3-(n/nnuc)}
      if (0.8-(n/nnuc) < 0){meaner=0} else {meaner=0.8-(n/nnuc)}
      
      chromatin$pE[n]=mysamp(n=1, m=meaner, s=0.05, lwr=lower, upr=upper, nnorm=100)
      
      chromatin$pDep0[n]=padd0[n]
      chromatin$pDep1[n]=padd1[n]*chromatin$pE[n]
      chromatin$pDep2[n]=padd2[n]*chromatin$pE[n]
      chromatin$pDep3[n]=padd3[n]*chromatin$pE[n]
      chromatin$pmax[n]=max(c(chromatin$pDep0[n],chromatin$pDep1[n],chromatin$pDep2[n],chromatin$pDep3[n]))
      #     print (paste ("Cell: ",cell," --- PRC2 passage: ",prc2))
      
      if (chromatin$S[n] == 0)
      {
        if(chromatin$pDep0[n] == chromatin$pmax[n])
        {
          chromatin$me0[n]=1
          chromatin$me1[n]=0
          chromatin$me2[n]=0
          chromatin$me3[n]=0
        } else if(chromatin$pDep1[n] == chromatin$pmax[n]) {
          chromatin$S[n]=1
          chromatin$me0[n]=0
          chromatin$me1[n]=1
          chromatin$me2[n]=0
          chromatin$me3[n]=0
        } else  if(chromatin$pDep2[n] == chromatin$pmax[n]) {
          chromatin$S[n]=2
          chromatin$me0[n]=0
          chromatin$me1[n]=0
          chromatin$me2[n]=1
          chromatin$me3[n]=0
        } else  if(chromatin$pDep3[n] == chromatin$pmax[n]) {
          chromatin$S[n]=3
          chromatin$me0[n]=0
          chromatin$me1[n]=0
          chromatin$me2[n]=0
          chromatin$me3[n]=1
        }
      } else if (chromatin$S[n] == 1)  {
        if(chromatin$pDep0[n] == chromatin$pmax[n])
        {
          chromatin$S[n]=1
          chromatin$me0[n]=0
          chromatin$me1[n]=1
          chromatin$me2[n]=0
          chromatin$me3[n]=0
        } else if(chromatin$pDep1[n] == chromatin$pmax[n] & rnorm(1,mean=0.5,sd=0.1) > burden12) {
          chromatin$S[n]=2
          chromatin$me0[n]=0
          chromatin$me1[n]=0
          chromatin$me2[n]=1
          chromatin$me3[n]=0
        } else  if(chromatin$pDep2[n] == chromatin$pmax[n]) {
          chromatin$S[n]=3
          chromatin$me0[n]=0
          chromatin$me1[n]=0
          chromatin$me2[n]=0
          chromatin$me3[n]=1
        } 
      } else if (chromatin$S[n] == 2)  {
        if(chromatin$pDep0[n] == chromatin$pmax[n])
        {
          chromatin$S[n]=2
          chromatin$me0[n]=0
          chromatin$me1[n]=0
          chromatin$me2[n]=1
          chromatin$me3[n]=0
        } else if(chromatin$pDep1[n] == chromatin$pmax[n]  & rnorm(1,mean=0.5,sd=0.1) > burden23)  {
          chromatin$S[n]=3
          chromatin$me0[n]=0
          chromatin$me1[n]=0
          chromatin$me2[n]=0
          chromatin$me3[n]=1
        }
      }
    }
    
    
    

    if (cell ==1 && prc2>0){
      #png(paste("~/Desktop/Histone_Mark_Simulation/Images/",prc2,"-3marks-average.png",sep=""),width = 800,height = 600)
      
      par(mfrow=c(4,1))
      plot(chromatin$me0,type = "h", ylim = c(0,1), main = paste("Cell 1: K27me0","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      plot(chromatin$me1,type = "h", ylim = c(0,1), main = paste("Cell 1: K27me1","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      plot(chromatin$me2,type = "h", ylim = c(0,1), main = paste("Cell 1: K27me2","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      plot(chromatin$me3,type = "h", ylim = c(0,1), main = paste("Cell 1: K27me3","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      #dev.off()
    }
    
    
    colnumber=(prc2-1)*3
    #      passages[,colnumber]=passages[,colnumber]+chromatin$me0
    passages[,(colnumber+1)]=passages[,(colnumber+1)]+chromatin$me1
    passages[,(colnumber+2)]=passages[,(colnumber+2)]+chromatin$me2
    passages[,(colnumber+3)]=passages[,(colnumber+3)]+chromatin$me3      
    
    
    
    
  }
  
  
  finalz$me0<-finalz$me0+chromatin$me0
  finalz$me1<-finalz$me1+chromatin$me1
  finalz$me2<-finalz$me2+chromatin$me2
  finalz$me3<-finalz$me3+chromatin$me3
  print ("cell #: ")
  print (cell)
  #print ("chromatin3")
  #print(chromatin$me3)
  
  
  # Mitosis  
  if (runif(n=1,min = 0,max = 100) >= 50  && prc2 > 5)
  {
    #Adding one to the number of final cell population
    final_cell_count=final_cell_count+1
    go_to_next_cell=1
    ## Removing the chromatin values already added to the pool (passages) so when we add the values of the daughter chromatins, it won't be doubled
    
    passages[,(colnumber+1)]=passages[,(colnumber+1)]-chromatin$me1
    passages[,(colnumber+2)]=passages[,(colnumber+2)]-chromatin$me2
    passages[,(colnumber+3)]=passages[,(colnumber+3)]-chromatin$me3      
    finalz$me0<-finalz$me0-chromatin$me0
    finalz$me1<-finalz$me1-chromatin$me1
    finalz$me2<-finalz$me2-chromatin$me2
    finalz$me3<-finalz$me3-chromatin$me3
    
    # Creating two empty daughter chromatins
    daughter1<-data.frame(integer(nnuc),integer(nnuc),integer(nnuc),integer(nnuc),integer(nnuc),integer(nnuc))
    colnames(daughter1)<-c("N","S","me0","me1","me2","me3")
    
    
    for (n in 1:nnuc)
    {
      daughter1$N[n]=n
      daughter1$S[n]=0
      daughter1$me0[n]=1
      daughter1$me1[n]=0
      daughter1$me2[n]=0
      daughter1$me3[n]=0
    }
    
    daughter2<-data.frame(integer(nnuc),integer(nnuc),integer(nnuc),integer(nnuc),integer(nnuc),integer(nnuc))
    colnames(daughter2)<-c("N","S","me0","me1","me2","me3")
    
    
    for (n in 1:nnuc)
    {
      daughter2$N[n]=n
      daughter2$S[n]=0
      daughter2$me0[n]=1
      daughter2$me1[n]=0
      daughter2$me2[n]=0
      daughter2$me3[n]=0
    }
    
    
    #randomly assigning marks from the mother chromatin to its daughter chromatins (likelihood of 50% on each nucleosome to be assigned to a different daughter chromatin)
    
    for (m in 1:nnuc){
      if (runif(n=1,min = 0,max = 100) > 50)
      {
        daughter1$S[m]=chromatin$S[m]
        daughter1$me0[m]=chromatin$me0[m]
        daughter1$me1[m]=chromatin$me1[m]
        daughter1$me2[m]=chromatin$me2[m]
        daughter1$me3[m]=chromatin$me3[m]
      }
      else {
        daughter2$S[m]=chromatin$S[m]
        daughter2$me0[m]=chromatin$me0[m]
        daughter2$me1[m]=chromatin$me1[m]
        daughter2$me2[m]=chromatin$me2[m]
        daughter2$me3[m]=chromatin$me3[m]
      }
    }
    
    
    
    # print daughter chromatins
    if (cell ==1 && prc2==prc2rounds){
      #png(paste("~/Desktop/Histone_Mark_Simulation/Images/",prc2,"-3marks-average.png",sep=""),width = 800,height = 600)
      par(mfrow=c(4,1))
      plot(daughter1$me0,type = "h", ylim = c(0,1), main = paste("Original Daughter1 Cell 1: K27me0","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      plot(daughter1$me1,type = "h", ylim = c(0,1), main = paste("Original Daughter1 Cell 1: K27me1","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      plot(daughter1$me2,type = "h", ylim = c(0,1), main = paste("Original Daughter1 Cell 1: K27me2","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      plot(daughter1$me3,type = "h", ylim = c(0,1), main = paste("Original Daughter1 Cell 1: K27me3","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      
      plot(daughter2$me0,type = "h", ylim = c(0,1), main = paste("Original Daughter2 Cell 2: K27me0","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      plot(daughter2$me1,type = "h", ylim = c(0,1), main = paste("Original Daughter2 Cell 2: K27me1","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      plot(daughter2$me2,type = "h", ylim = c(0,1), main = paste("Original Daughter2 Cell 2: K27me2","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      plot(daughter2$me3,type = "h", ylim = c(0,1), main = paste("Original Daughter2 Cell 2: K27me3","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
      
      
      #dev.off()
    }
    
currentprc2round=prc2    
    ####### Now putting daughter chromatins through another round of Lysine methylation  ################
    ################## Daughter 1 #######################
    ### Runs of PRC2
    for (prc2 in 1:prc2rounds) {
      ### Depositing marks on daughter1
      for (n in 1:nnuc)
      {
        
        padd1[n]=mysamp(n=1, m=0.69, s=0.1, lwr=0, upr=1, nnorm=100)
        padd2[n]=mysamp(n=1, m=0.01, s=0.1, lwr=0, upr=1, nnorm=100)
        padd3[n]=mysamp(n=1, m=0.00, s=0.1, lwr=0, upr=1, nnorm=100)
        padd0[n]=1-(sum(padd1[n],padd2[n],padd3[n]))
        
        
        # PRC2 halflife
        if (0.3-(n/nnuc) < 0 ){lower=0} else {lower=0.3-(n/nnuc)}
        if (1.3-(n/nnuc) > 1){upper=1} else {upper=1.3-(n/nnuc)}
        if (0.8-(n/nnuc) < 0){meaner=0} else {meaner=0.8-(n/nnuc)}
        
        daughter1$pE[n]=mysamp(n=1, m=meaner, s=0.05, lwr=lower, upr=upper, nnorm=100)
        
        daughter1$pDep0[n]=padd0[n]
        daughter1$pDep1[n]=padd1[n]*daughter1$pE[n]
        daughter1$pDep2[n]=padd2[n]*daughter1$pE[n]
        daughter1$pDep3[n]=padd3[n]*daughter1$pE[n]
        daughter1$pmax[n]=max(c(daughter1$pDep0[n],daughter1$pDep1[n],daughter1$pDep2[n],daughter1$pDep3[n]))
        #     print (paste ("Cell: ",cell," --- PRC2 passage: ",prc2))
        
        if (daughter1$S[n] == 0)
        {
          if(daughter1$pDep0[n] == daughter1$pmax[n])
          {
            daughter1$me0[n]=1
            daughter1$me1[n]=0
            daughter1$me2[n]=0
            daughter1$me3[n]=0
          } else if(daughter1$pDep1[n] == daughter1$pmax[n]) {
            daughter1$S[n]=1
            daughter1$me0[n]=0
            daughter1$me1[n]=1
            daughter1$me2[n]=0
            daughter1$me3[n]=0
          } else  if(daughter1$pDep2[n] == daughter1$pmax[n]) {
            daughter1$S[n]=2
            daughter1$me0[n]=0
            daughter1$me1[n]=0
            daughter1$me2[n]=1
            daughter1$me3[n]=0
          } else  if(daughter1$pDep3[n] == daughter1$pmax[n]) {
            daughter1$S[n]=3
            daughter1$me0[n]=0
            daughter1$me1[n]=0
            daughter1$me2[n]=0
            daughter1$me3[n]=1
          }
        } else if (daughter1$S[n] == 1)  {
          if(daughter1$pDep0[n] == daughter1$pmax[n])
          {
            daughter1$S[n]=1
            daughter1$me0[n]=0
            daughter1$me1[n]=1
            daughter1$me2[n]=0
            daughter1$me3[n]=0
          } else if(daughter1$pDep1[n] == daughter1$pmax[n] & rnorm(1,mean=0.5,sd=0.1) > burden12) {
            daughter1$S[n]=2
            daughter1$me0[n]=0
            daughter1$me1[n]=0
            daughter1$me2[n]=1
            daughter1$me3[n]=0
          } else  if(daughter1$pDep2[n] == daughter1$pmax[n]) {
            daughter1$S[n]=3
            daughter1$me0[n]=0
            daughter1$me1[n]=0
            daughter1$me2[n]=0
            daughter1$me3[n]=1
          } 
        } else if (daughter1$S[n] == 2)  {
          if(daughter1$pDep0[n] == daughter1$pmax[n])
          {
            daughter1$S[n]=2
            daughter1$me0[n]=0
            daughter1$me1[n]=0
            daughter1$me2[n]=1
            daughter1$me3[n]=0
          } else if(daughter1$pDep1[n] == daughter1$pmax[n]  & rnorm(1,mean=0.5,sd=0.1) > burden23)  {
            daughter1$S[n]=3
            daughter1$me0[n]=0
            daughter1$me1[n]=0
            daughter1$me2[n]=0
            daughter1$me3[n]=1
          }
        }
      }
      
      
      
      #prc2=(prc2)+1
      #print (prc2)
      
      if (cell ==1 && prc2>0){
        #png(paste("~/Desktop/Histone_Mark_Simulation/Images/",prc2,"-3marks-average.png",sep=""),width = 800,height = 600)
        
        par(mfrow=c(4,1))
        plot(daughter1$me0,type = "h", ylim = c(0,1), main = paste("Daughter1 Cell 1: K27me0","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
        plot(daughter1$me1,type = "h", ylim = c(0,1), main = paste("Daughter1 Cell 1: K27me1","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
        plot(daughter1$me2,type = "h", ylim = c(0,1), main = paste("Daughter1 Cell 1: K27me2","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
        plot(daughter1$me3,type = "h", ylim = c(0,1), main = paste("Daughter1 Cell 1: K27me3","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
        #dev.off()
      }
      
      
      colnumber=(prc2-1)*3
      #passages[,colnumber]=passages[,colnumber]+daughter1$me0
      passages[,(colnumber+1)]=passages[,(colnumber+1)]+daughter1$me1
      passages[,(colnumber+2)]=passages[,(colnumber+2)]+daughter1$me2
      passages[,(colnumber+3)]=passages[,(colnumber+3)]+daughter1$me3
      finalz$me0<-finalz$me0+daughter1$me0
      finalz$me1<-finalz$me1+daughter1$me1
      finalz$me2<-finalz$me2+daughter1$me2
      finalz$me3<-finalz$me3+daughter1$me3
      
      
      
      
    }
    ################## Daughter 2 #######################
    
    ### Runs of PRC2
    for (prc2 in 1:prc2rounds) {
      ### Depositing marks on daughter 2
      for (n in 1:nnuc)
      {
        
        padd1[n]=mysamp(n=1, m=0.69, s=0.1, lwr=0, upr=1, nnorm=100)
        padd2[n]=mysamp(n=1, m=0.01, s=0.1, lwr=0, upr=1, nnorm=100)
        padd3[n]=mysamp(n=1, m=0.00, s=0.1, lwr=0, upr=1, nnorm=100)
        padd0[n]=1-(sum(padd1[n],padd2[n],padd3[n]))
        
        
        # PRC2 halflife
        if (0.3-(n/nnuc) < 0 ){lower=0} else {lower=0.3-(n/nnuc)}
        if (1.3-(n/nnuc) > 1){upper=1} else {upper=1.3-(n/nnuc)}
        if (0.8-(n/nnuc) < 0){meaner=0} else {meaner=0.8-(n/nnuc)}
        
        daughter2$pE[n]=mysamp(n=1, m=meaner, s=0.05, lwr=lower, upr=upper, nnorm=100)
        
        daughter2$pDep0[n]=padd0[n]
        daughter2$pDep1[n]=padd1[n]*daughter2$pE[n]
        daughter2$pDep2[n]=padd2[n]*daughter2$pE[n]
        daughter2$pDep3[n]=padd3[n]*daughter2$pE[n]
        daughter2$pmax[n]=max(c(daughter2$pDep0[n],daughter2$pDep1[n],daughter2$pDep2[n],daughter2$pDep3[n]))
        #     print (paste ("Cell: ",cell," --- PRC2 passage: ",prc2))
        
        if (daughter2$S[n] == 0)
        {
          if(daughter2$pDep0[n] == daughter2$pmax[n])
          {
            daughter2$me0[n]=1
            daughter2$me1[n]=0
            daughter2$me2[n]=0
            daughter2$me3[n]=0
          } else if(daughter2$pDep1[n] == daughter2$pmax[n]) {
            daughter2$S[n]=1
            daughter2$me0[n]=0
            daughter2$me1[n]=1
            daughter2$me2[n]=0
            daughter2$me3[n]=0
          } else  if(daughter2$pDep2[n] == daughter2$pmax[n]) {
            daughter2$S[n]=2
            daughter2$me0[n]=0
            daughter2$me1[n]=0
            daughter2$me2[n]=1
            daughter2$me3[n]=0
          } else  if(daughter2$pDep3[n] == daughter2$pmax[n]) {
            daughter2$S[n]=3
            daughter2$me0[n]=0
            daughter2$me1[n]=0
            daughter2$me2[n]=0
            daughter2$me3[n]=1
          }
        } else if (daughter2$S[n] == 1)  {
          if(daughter2$pDep0[n] == daughter2$pmax[n])
          {
            daughter2$S[n]=1
            daughter2$me0[n]=0
            daughter2$me1[n]=1
            daughter2$me2[n]=0
            daughter2$me3[n]=0
          } else if(daughter2$pDep1[n] == daughter2$pmax[n] & rnorm(1,mean=0.5,sd=0.1) > burden12) {
            daughter2$S[n]=2
            daughter2$me0[n]=0
            daughter2$me1[n]=0
            daughter2$me2[n]=1
            daughter2$me3[n]=0
          } else  if(daughter2$pDep2[n] == daughter2$pmax[n]) {
            daughter2$S[n]=3
            daughter2$me0[n]=0
            daughter2$me1[n]=0
            daughter2$me2[n]=0
            daughter2$me3[n]=1
          } 
        } else if (daughter2$S[n] == 2)  {
          if(daughter2$pDep0[n] == daughter2$pmax[n])
          {
            daughter2$S[n]=2
            daughter2$me0[n]=0
            daughter2$me1[n]=0
            daughter2$me2[n]=1
            daughter2$me3[n]=0
          } else if(daughter2$pDep1[n] == daughter2$pmax[n]  & rnorm(1,mean=0.5,sd=0.1) > burden23)  {
            daughter2$S[n]=3
            daughter2$me0[n]=0
            daughter2$me1[n]=0
            daughter2$me2[n]=0
            daughter2$me3[n]=1
          }
        }
      }
      
      
      
      #prc2=(prc2)+1
      #print (prc2)
      
      if (cell ==1 && prc2>0){
        #png(paste("~/Desktop/Histone_Mark_Simulation/Images/",prc2,"-3marks-average.png",sep=""),width = 800,height = 600)
        
        par(mfrow=c(4,1))
        plot(daughter2$me0,type = "h", ylim = c(0,1), main = paste("Daughter2 Cell 1: K27me0","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
        plot(daughter2$me1,type = "h", ylim = c(0,1), main = paste("Daughter2 Cell 1: K27me1","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
        plot(daughter2$me2,type = "h", ylim = c(0,1), main = paste("Daughter2 Cell 1: K27me2","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
        plot(daughter2$me3,type = "h", ylim = c(0,1), main = paste("Daughter2 Cell 1: K27me3","  -- PRC2 passage: ",prc2," Cell: ",cell,sep=""),ylab="")
        #dev.off()
      }
      
      
      colnumber=(prc2-1)*3
      #      passages[,colnumber]=passages[,colnumber]+daughter2$me0
      passages[,(colnumber+1)]=passages[,(colnumber+1)]+daughter2$me1
      passages[,(colnumber+2)]=passages[,(colnumber+2)]+daughter2$me2
      passages[,(colnumber+3)]=passages[,(colnumber+3)]+daughter2$me3      
      finalz$me0<-finalz$me0+daughter2$me0
      finalz$me1<-finalz$me1+daughter2$me1
      finalz$me2<-finalz$me2+daughter2$me2
      finalz$me3<-finalz$me3+daughter2$me3
      
      
      
    }    
    
    

    ####### End of putting daughter chromatins through another round of Lysine methylation  #############
    
    
  }
  ### End of mitosis
  if (go_to_next_cell==1) break
  
  
  
}




finalz$me0av<-finalz$me0/cell
finalz$me1av<-finalz$me1/cell
finalz$me2av<-finalz$me2/cell
finalz$me3av<-finalz$me3/cell




for (state in 0:(prc2rounds-1))
{
  #png(paste("~/Desktop/Histone_Mark_Simulation/Images/Average_State-",state,"-3marks-average.png",sep=""),width = 800,height = 600)
  par(mfrow=c(3,1))
  
  plot((passages[,(3*state)+1]/final_cell_count),type = "h", ylim = c(0,1), main = paste("K27me1","  -- Average over : ",cell," Cells"," PRC2 passage: ",state,sep=""),ylab="")
  plot((passages[,(3*state)+2]/final_cell_count),type = "h", ylim = c(0,1), main = paste("K27me2","  -- Average over : ",cell," Cells"," PRC2 passage: ",state,sep=""),ylab="")
  plot((passages[,(3*state)+3]/final_cell_count),type = "h", ylim = c(0,1), main = paste("K27me3","  -- Average over : ",cell," Cells"," PRC2 passage: ",state,sep=""),ylab="")
  #dev.off()
  
}



# for (state in (prc2rounds-1):(prc2rounds-1))
# {
#   #png(paste("~/Desktop/Histone_Mark_Simulation/Images/Average_State-",state,"-3marks-average.png",sep=""),width = 800,height = 600)
#   par(mfrow=c(3,1))
#   
#   plot((passages[,(3*state)+1]/initial_cell_count),type = "h", ylim = c(0,1), main = paste("K27me1","  -- Average over : ",cell," Cells"," PRC2 passage: ",state,sep=""),ylab="")
#   plot((passages[,(3*state)+2]/initial_cell_count),type = "h", ylim = c(0,1), main = paste("K27me2","  -- Average over : ",cell," Cells"," PRC2 passage: ",state,sep=""),ylab="")
#   plot((passages[,(3*state)+3]/initial_cell_count),type = "h", ylim = c(0,1), main = paste("K27me3","  -- Average over : ",cell," Cells"," PRC2 passage: ",state,sep=""),ylab="")
#   #dev.off()
#   
# }
