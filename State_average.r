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
mitosis=10
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


passages<- data.frame(matrix(nrow = nnuc,ncol = (3*mitosis)))
passcolnames<-c()
#passcolnames<-c("N")
for (i in 1:mitosis){
#  passcolnames<-c(passcolnames,paste(i,"me0",sep="_"))
  passcolnames<-c(passcolnames,paste(i,"me1",sep="_"))
  passcolnames<-c(passcolnames,paste(i,"me2",sep="_"))
  passcolnames<-c(passcolnames,paste(i,"me3",sep="_"))
}
colnames(passages)<-passcolnames
 for (n in 1:nnuc)
 {
#   passages$N[n]=n
  passages[n,]=0
 }




 for (cell in 1:50)  {
   set.seed(cell)
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
     png(paste("~/Desktop/Histone_Mark_Simulation/Images/000-3marks-average.png",sep=""),width = 800,height = 600)
   
     par(mfrow=c(4,1))
     plot(chromatin$me0,type = "h", ylim = c(0,1), main = paste("State 0 : K27me0","  -- PRC2 passage: ",prc2,sep=""),ylab="")
     plot(chromatin$me1,type = "h", ylim = c(0,1), main = paste("State 0 : K27me1","  -- PRC2 passage: ",prc2,sep=""),ylab="")
     plot(chromatin$me2,type = "h", ylim = c(0,1), main = paste("State 0 : K27me2","  -- PRC2 passage: ",prc2,sep=""),ylab="")
     plot(chromatin$me3,type = "h", ylim = c(0,1), main = paste("State 0 : K27me3","  -- PRC2 passage: ",prc2,sep=""),ylab="")
   dev.off()
    }

    padd0=c()
    padd1=c()
    padd2=c()
    padd3=c()


    ### Runs of PRC2
    for (prc2 in 1:mitosis) {
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
  
  
  
#prc2=(prc2)+1
#print (prc2)

if (cell ==1 && prc2>0){
png(paste("~/Desktop/Histone_Mark_Simulation/Images/",prc2,"-3marks-average.png",sep=""),width = 800,height = 600)
  
 par(mfrow=c(4,1))
 plot(chromatin$me0,type = "h", ylim = c(0,1), main = paste("Cell 1: K27me0","  -- PRC2 passage: ",prc2,sep=""),ylab="")
 plot(chromatin$me1,type = "h", ylim = c(0,1), main = paste("Cell 1: K27me1","  -- PRC2 passage: ",prc2,sep=""),ylab="")
 plot(chromatin$me2,type = "h", ylim = c(0,1), main = paste("Cell 1: K27me2","  -- PRC2 passage: ",prc2,sep=""),ylab="")
 plot(chromatin$me3,type = "h", ylim = c(0,1), main = paste("Cell 1: K27me3","  -- PRC2 passage: ",prc2,sep=""),ylab="")
dev.off()
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
print ("iteration")
print (cell)
#print ("chromatin3")
#print(chromatin$me3)


}




finalz$me0av<-finalz$me0/cell
finalz$me1av<-finalz$me1/cell
finalz$me2av<-finalz$me2/cell
finalz$me3av<-finalz$me3/cell

# 
#png(paste("~/Desktop/Histone_Mark_Simulation/Images/999-Average-over ",cell,"-3marks-average.png",sep=""),width = 800,height = 600)

par(mfrow=c(3,1))
#plot(finalz$me0av,type = "h", ylim = c(0,1), main = paste("K27me0","  -- Average over : ",cell," Cells",sep=""),ylab="")
plot(finalz$me1av,type = "h", ylim = c(0,1), main = paste("K27me1","  -- Average over : ",cell," Cells",sep=""),ylab="")
plot(finalz$me2av,type = "h", ylim = c(0,1), main = paste("K27me2","  -- Average over : ",cell," Cells",sep=""),ylab="")
plot(finalz$me3av,type = "h", ylim = c(0,1), main = paste("K27me3","  -- Average over : ",cell," Cells",sep=""),ylab="")
#dev.off()


 for (state in 0:(mitosis-1))
 {
   png(paste("~/Desktop/Histone_Mark_Simulation/Images/Average_State-",state,"-3marks-average.png",sep=""),width = 800,height = 600)
   par(mfrow=c(3,1))
   
   plot((passages[,(3*state)+1]/cell),type = "h", ylim = c(0,1), main = paste("K27me1","  -- Average over : ",cell," Cells",sep=""),ylab="")
   plot((passages[,(3*state)+2]/cell),type = "h", ylim = c(0,1), main = paste("K27me2","  -- Average over : ",cell," Cells",sep=""),ylab="")
   plot((passages[,(3*state)+3]/cell),type = "h", ylim = c(0,1), main = paste("K27me3","  -- Average over : ",cell," Cells",sep=""),ylab="")
   dev.off()
   
}

# 3 4 5
# 7 8 9















# 
# x<-c() 
# nnuc=1000
# for (n in 1:nnuc){if (0.5-(n/nnuc) < 0 ){lower=0} else {lower=0.5-(n/nnuc)};if (1.5-(n/nnuc) > 1){upper=1} else {upper=1.5-(n/nnuc)};x<-c(x,mysamp(n=1, m=(1-(n/nnuc)), s=0.1, lwr=lower, upr=upper, nnorm=10000))}
# par(mfrow=c(1,0)) 
# plot(x)
# summary(x)
