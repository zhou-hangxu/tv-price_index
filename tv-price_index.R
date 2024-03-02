setwd("/Users/ZHOUHANGXU/Desktop/期末report")
p <- read.csv("./pdata.csv", header = TRUE)
q <- read.csv("./qdata.csv", header = TRUE)

#View(p)
#View(q)

### (1)指数算式
## unweighted indexes
# carli
PC <- NULL
for(t in 1:84){
  tmp <- 0
  for(i in 1:402){
    tmp <- tmp + (p[t,i]/p[1,i])
  }
  PC <- c(PC,tmp/402)
}

# jevons
PJ <- NULL
for(t in 1:84){
  tmp <- 1
  for(i in 1:402){
    tmp <- tmp * ((p[t,i]/p[1,i])^(1/402))
  }
  PJ <- c(PJ,tmp)
}

# dutot
PD <- NULL
for(t in 1:84){
  tmp <- sum(p[t,])/sum(p[1,])
  PD <- c(PD,tmp)
}

# graph
x11(width = 8, height = 4)
plot(PC, type='l',col='blue', xlab='time', ylab='indexes', main = 'unweighted indexes')
lines(PD, col='red')
lines(PJ, col='green')
legend('top', legend=c('Carli','Dutot','Jevons'),col=c('blue','red','green'),lty=1, ncol=3)
##weighted 
# laspeyres
PL <- NULL
for(t in 1:84){
  PL <- c(PL,sum(p[t,]*q[1,])/sum(p[1,]*q[1,]))
}

# paasche
PP <- NULL
for(t in 1:84){
  PP <- c(PP,sum(p[t,]*q[t,])/sum(p[1,]*q[t,]))
}

# fisher
PF <- (PP*PL)^(1/2)

# walsh
PW <- NULL
for(t in 1:84){
  PW <- c(PW,sum(p[t,]*(q[t,]*q[1,])^(1/2))/sum(p[1,]*(q[t,]*q[1,])^(1/2)))
}

# tornquist
PT <- NULL
for(t in 1:84){
  tmp <- 1
  for(i in 1:402){
    w1 <- p[1,i]*q[1,i]/sum(p[1,]*q[1,])
    wt <- p[t,i]*q[t,i]/sum(p[t,]*q[t,])
    tmp <- tmp*(p[t,i]*p[1,i])^(0.5*(w1+wt))
  }
 PT <- c(PT, tmp)
}

# graph
x11(width = 8, height = 4)
plot(PL, type='l',col='blue', xlab='time', ylab='indexes', main = 'weighted indexes')
lines(PP, col='red')
lines(PF, col='green')
lines(PW, col='grey')
lines(PT, col='purple')
legend('top', legend=c('laspeyres','paasche','fisher','walsh','tornquist'),col=c('blue','red','green','grey','purple'),lty=1, ncol=4)


## weighted(chain)
# laspeyres(chain)
PLC <- 1
for(t in 2:84){
  PLC <- c(PLC,PLC[t-1]*sum(p[t,]*q[t-1,])/sum(p[t-1,]*q[t-1,]))
}

# paasche(chain)
PPC <- 1
for(t in 2:84){
  PPC <- c(PPC,PPC[t-1]*sum(p[t,]*q[t,])/sum(p[t-1,]*q[t,]))
}

# fisher(chain)
PFC <- (PLC*PPC)^(1/2)

# graph
x11(width = 8, height = 4)
plot(PLC, type='l',col='blue', xlab='time', ylab='indexes', main = 'weighted indexes(chain)')
lines(PPC, col='red')
lines(PFC, col='green')
legend('top', legend=c('laspeyres','paasche','fisher'),col=c('blue','red','green'),lty=1, ncol=3)


### (2)品質調整済み価格指数

## data process
jan <- rep(1:ncol(p),each=nrow(p))
time <- rep(rownames(p),ncol(p))
p1<- as.vector(as.matrix(p))
q1 <- as.vector(as.matrix(q))
df <- data.frame(JAN = jan, TIME = time, P = p1, Q = q1)

# time dummy
tlist <- unique(sort(df$TIME))
df_D <- NULL
for(t in 1:84){
  tmp <- ifelse(df$TIME==tlist[t],1,0)
  if(t==1){
    df_D <- data.frame(D01=tmp)
  }else{
    df_D <- data.frame(df_D,tmp)
    tmp <- paste("D",formatC(t,width=2,flag=0),sep="")
    names(df_D)[dim(df_D)[2]] <- tmp
  }
}

# jan dummy
JANlist <- unique(sort(df$JAN))
df_JAN <- NULL
for(n in 1:402){
  tmp <- ifelse(df$JAN==JANlist[n],1,0)
  if(n==1){
    df_JAN <- data.frame(J001=tmp)
  }else{
    df_JAN <- data.frame(df_JAN,tmp)
    tmp <- paste("J",formatC(n,width=3,flag=0),sep="")
    names(df_JAN)[dim(df_JAN)[2]] <- tmp
  }
}

## unweighted time product dummy (PTPD)
res <- lm("log(P)~.+0",data=data.frame(P=df$P,df_D[,-1],df_JAN))
summary(res)
summary(res)$coef
PTPD <- c(1,exp(summary(res)$coef[1:83]))

## weighted time product dummy (WPTPD)
df$E = df$P * df$Q

# calculate weight
EA <- tapply(df$E,df$TIME,sum)
df_EAI <- data.frame(TIME=names(EA),EAI=1/EA)
df <- merge(df,df_EAI,by="TIME")
df$S <- df$E*df$EAI
df$SSR <- df$S**0.5
df$SLP <- df$SSR*log(df$P)
df_DS <- df_D*df$SSR
df_JANS <- df_JAN*df$SSR

# regression
res <- lm("SLP~.+0",data=data.frame(SLP=df$SLP,df_DS[,-1],df_JANS))
summary(res)
summary(res)$coef
WPTPD <- c(1,exp(summary(res)$coef[1:83]))

## price index with adjacent periods(PAP)
## weighted price index with adjacent periods(WPAP)
Tlist <- unique(sort(df$TIME))
Rselect <- which(df$TIME==Tlist[1]|df$TIME==Tlist[2])
res <- lm("log(P)~.+0",data=data.frame(P=df$P[Rselect],ONE=1,df_D[Rselect,-1],df_JAN[Rselect,]))
summary(res)
summary(res)$coef
PAP <- c(1,exp(summary(res)$coef[2]))

res <- lm("SLP~.+0",data=data.frame(SLP=df$SLP[Rselect],SSR=df$SSR[Rselect],df_DS[Rselect,-1],df_JANS[Rselect,]))
PWAP <- c(1,exp(summary(res)$coef[2]))

for(t in 3:84){
  Rselect <- which(df$TIME==Tlist[t-1]|df$TIME==Tlist[t])
  res <- lm("log(P)~.+0",data=data.frame(P=df$P[Rselect],ONE=1,df_D[Rselect,-c(1:(t-1))],df_JAN[Rselect,]))
  PAP <- c(PAP,PAP[t-1]*exp(summary(res)$coef[2]))  
  res <- lm("SLP~.+0",data=data.frame(SLP=df$SLP[Rselect],SSR=df$SSR[Rselect],df_DS[Rselect,-c(1:(t-1))],df_JANS[Rselect,]))
  PWAP <- c(PWAP,PWAP[t-1]*exp(summary(res)$coef[2]))
}          

## average price(PA)
PAT <- NULL
for(t in 1:84){
  N = sum(q[t,])
  PAT <- c(PAT,sum(p[t,])/N)
}

PA <- PAT/PAT[1]

## unit value price(PUV)
PUVT <- NULL
for(t in 1:84){
  Q = sum(q[t,])
  PQ = sum(p[t,]*q[t,])
  PUVT <- c(PUVT, PQ/Q)
}
PUV <- PUVT/PUVT[1]

# graph
x11(width = 8, height = 4)
plot(PTPD, type='l',col='blue', xlab='time', ylab='indexes', main = 'quality_adjusted_indexs',ylim=c(0,2))
lines(WPTPD, col='red')
lines(PAP, col='green')
lines(PWAP, col='grey')
lines(PA, col='black')
lines(PUV, col='purple')
legend('bottom', legend=c('PTPD','WPTPD','PAP','PWAP','PA','PUV'),col=c('blue','red','green','grey','black','purple'),lty=1, ncol=3)


### (3)季節財

## 欠損値
p$month <- rep(1:12, times = nrow(p)/12)
q$month <- rep(1:12, times = nrow(q)/12)
mlist <- unique(p$month)
ssec2idx <- matrix(1,ncol=10,nrow=1)

for(mm in 1:12){
  pp <- p[p$month==mlist[mm],c(1:402)]
  qq <- q[q$month==mlist[mm],c(1:402)]
  for(nn in 1:402){
    if(length(which(pp[,nn]==0))!=0&length(which(pp[,nn]==0))!=7){
      for(yy in 1:7){
        if(yy==1&pp[yy,nn]==0){
          pp[yy,nn] <- pp[min(which(pp[,nn]!=0)),nn]
        }else if(pp[yy,nn]==0){
          tmp <- which(pp[,nn]!=0)
          pp[yy,nn] <- pp[max(tmp[tmp<yy]),nn]
        }
      }
    }
  }
}

##3.Monthly Indexes using Carry Forward prices
# Laspeyres、Paasche、Fisher
PL <- NULL
PP <- NULL
for(ii in 1:7){
  PL <- c(PL,sum(pp[ii,]*qq[1,])/sum(pp[1,]*qq[1,]))
  PP <- c(PP,sum(pp[ii,]*qq[ii,])/sum(pp[1,]*qq[ii,]))
}
PF <- sqrt(PP*PL)

# tornqvist
PT <- NULL
for(ii in 1:7){
  tmp <- 0
  for(nn in 1:402){
    if(pp[ii,nn]!=0&pp[1,nn]!=0){
      s1 <- pp[1,nn]*qq[1,nn]/sum(pp[1,]*qq[1,])
      st <- pp[ii,nn]*qq[ii,nn]/sum(pp[ii,]*qq[ii,])
      tmp <- tmp + 0.5*(s1+st)*log(pp[ii,nn]/pp[1,nn])
    }
  }
  PT <- c(PT,exp(tmp))
}

# chain Index
PLC <- 1
for(ii in 2:7){
  PLC <- c(PLC,PLC[ii-1]*(sum(pp[ii,]*qq[ii-1,])/sum(pp[ii-1,]*qq[ii-1,])))
}

PPC <- 1
for(ii in 2:7){
  PPC <- c(PPC,PPC[ii-1]*(sum(pp[ii,]*qq[ii,])/sum(pp[ii-1,]*qq[ii,])))
}

PFC <- sqrt(PLC*PPC)

PTC <- 1
for(ii in 2:7){
  tmp <- 0
  for(nn in 1:402){
    if(pp[ii,nn]!=0&pp[ii-1,nn]!=0){
      s1 <- pp[ii-1,nn]*qq[ii-1,nn]/sum(pp[ii-1,]*qq[ii-1,])
      st <- pp[ii,nn]*qq[ii,nn]/sum(pp[ii,]*qq[ii,])
      tmp <- tmp + 0.5*(s1+st)*log(pp[ii,nn]/pp[ii-1,nn])
    }
  }
  PTC <- c(PTC,PTC[ii-1]*exp(tmp))
}

# GEKS
PGEKS <- NULL
for(yy in 1:7){
  tmp <- 1
  for(zz in 1:7){
    PLyz <- sum(pp[yy,]*qq[zz,])/sum(pp[zz,]*qq[zz,])
    PPyz <- sum(pp[yy,]*qq[yy,])/sum(pp[zz,]*qq[yy,])
    PFyz <- (PLyz*PPyz)**(1/2)
    tmp <- tmp * PFyz
  }
  PGEKS <- c(PGEKS,(tmp)^(1/7))
}
PGEKS <- PGEKS/PGEKS[1]

# Dissimilarity
Delta <- matrix(0,nrow=7,ncol=7)
for(y in 1:7){
  for(z in 1:7){
    yy <- sum(pp[y,]*qq[y,])
    zy <- sum(pp[z,]*qq[y,])
    zz <- sum(pp[z,]*qq[z,])
    yz <- sum(pp[y,]*qq[z,])
    Delta1 <- 0
    Delta2 <- 0
    for(nn in 1:402){
      Delta1 <- Delta1 + (pp[y,nn]*qq[y,nn]/yy-pp[z,nn]*qq[y,nn]/zy)**2
      Delta2 <- Delta2 + (pp[z,nn]*qq[z,nn]/zz-pp[y,nn]*qq[z,nn]/yz)**2
    }
    Delta[y,z] <- Delta1+Delta2
  }
}

# Similarity Linked Index
PS <- 1
PS <- c(PS,PF[2])
for(yy in 3:7){
  tt <- which(Delta[yy,1:yy-1]==min(Delta[yy,1:yy-1]))
  PLCtmp <- (sum(pp[yy,]*qq[tt,])/sum(pp[tt,]*qq[tt,]))
  PPCtmp <- (sum(pp[yy,]*qq[yy,])/sum(pp[tt,]*qq[yy,]))
  PS <- c(PS,PS[tt]*(PLCtmp*PPCtmp)^(1/2))
}
ssec2idx <- rbind(ssec2idx,cbind(ssec2idx[dim(ssec2idx)[1],1]*PL,
                                 ssec2idx[dim(ssec2idx)[1],2]*PP,
                                 ssec2idx[dim(ssec2idx)[1],3]*PF,
                                 ssec2idx[dim(ssec2idx)[1],4]*PT,
                                 ssec2idx[dim(ssec2idx)[1],5]*PLC,
                                 ssec2idx[dim(ssec2idx)[1],6]*PPC,
                                 ssec2idx[dim(ssec2idx)[1],7]*PFC,
                                 ssec2idx[dim(ssec2idx)[1],8]*PTC,
                                 ssec2idx[dim(ssec2idx)[1],9]*PGEKS,
                                 ssec2idx[dim(ssec2idx)[1],10]*PS))

write.csv(ssec2idx,"ssec2idx.csv")

x11(width = 8, height = 4)
plot(PL, type='l',col='blue', xlab='time', ylab='indexes', main = 'monthly indexes using carry forward prices')
lines(PP, col='red')
lines(PF, col='green')
lines(PT, col='grey')
lines(PLC, col='black')
lines(PPC, col='pink')
lines(PFC, COL='yellow')
lines(PTC, COL='brown')
lines(PGEKS, COL='orange')
lines(PS, COL='#FF0090')
legend('topleft', legend=c('PL','PP','PF','PT','PLC','PPC','PFC','PTC','PGEKS','PS'),col=c('blue','red','green','grey','black','pink','yellow','brown','orange','#FF0090'),lty=1, ncol=5)

## 4.maximum overlap year over year monthly indexes sec3

pp <- p[p$month=="1",c(1:402)]
qq <- q[q$month=="1",c(1:402)]
ssec3idx <- matrix(1,ncol=10,nrow=1)
for(mm in 1:12){
pp <- p[p$month==mlist[mm],c(1:402)]
qq <- q[q$month==mlist[mm],c(1:402)]
}
PL <- NULL
PP <- NULL
for(ii in 1:7){
    tmp1 <- ifelse(pp[ii,]>0,1,0)
    tmp2 <- ifelse(pp[1,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    PL <- c(PL,sum(pp[ii,overlap]*qq[1,overlap])/sum(pp[1,overlap]*qq[1,overlap]))
    PP <- c(PP,sum(pp[ii,overlap]*qq[ii,overlap])/sum(pp[1,overlap]*qq[ii,overlap]))
}
print(cbind(PL,PP))
PF <- sqrt(PP*PL)
print(PF)
PLC <- 1
PPC <- 1
for(ii in 2:7){
    tmp1 <- ifelse(pp[ii,]>0,1,0)
    tmp2 <- ifelse(pp[ii-1,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    PLC <- c(PLC,PLC[ii-1]*(sum(pp[ii,overlap]*qq[ii-1,overlap])/sum(pp[ii-1,overlap]*qq[ii-1,overlap])))
    PPC <- c(PPC,PPC[ii-1]*(sum(pp[ii,overlap]*qq[ii,overlap])/sum(pp[ii-1,overlap]*qq[ii,overlap])))
}
PFC <- sqrt(PLC*PPC)
cbind(PLC,PPC,PFC)
PT <- NULL
for(ii in 1:7){
    tmp <- 0
    tmp1 <- ifelse(pp[ii,]>0,1,0)
    tmp2 <- ifelse(pp[1,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    for(nn in overlap){
      if(pp[ii,nn]!=0&pp[1,nn]!=0){
        s1 <- pp[1,nn]*qq[1,nn]/sum(pp[1,overlap]*qq[1,overlap])
        st <- pp[ii,nn]*qq[ii,nn]/sum(pp[ii,overlap]*qq[ii,overlap])
        tmp <- tmp + 0.5*(s1+st)*log(pp[ii,nn]/pp[1,nn])
      }
    }
    PT <- c(PT,exp(tmp))
}
PT
PTC <- 1
for(ii in 2:7){
    tmp <- 0
    tmp1 <- ifelse(pp[ii,]>0,1,0)
    tmp2 <- ifelse(pp[ii-1,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    for(nn in overlap){
      if(pp[ii,nn]!=0&pp[ii-1,nn]!=0){
        s1 <- pp[ii-1,nn]*qq[ii-1,nn]/sum(pp[ii-1,overlap]*qq[ii-1,overlap])
        st <- pp[ii,nn]*qq[ii,nn]/sum(pp[ii,overlap]*qq[ii,overlap])
        tmp <- tmp + 0.5*(s1+st)*log(pp[ii,nn]/pp[ii-1,nn])
      }
    }
    PTC <- c(PTC,PTC[ii-1]*exp(tmp))
}
PTC
PGEKS <- NULL
for(yy in 1:7){
    tmp <- 1
    tmp1 <- ifelse(pp[yy,]>0,1,0)
    for(zz in 1:7){
      tmp2 <- ifelse(pp[zz,]>0,1,0)
      overlap <- which(tmp1*tmp2==1)
      PLyz <- sum(pp[yy,overlap]*qq[zz,overlap])/sum(pp[zz,overlap]*qq[zz,overlap])
      PPyz <- sum(pp[yy,overlap]*qq[yy,overlap])/sum(pp[zz,overlap]*qq[yy,overlap])
      PFyz <- (PLyz*PPyz)**(1/2)
      tmp <- tmp * PFyz
    }
    PGEKS <- c(PGEKS,(tmp)^(1/7))
}
PGEKS <- PGEKS/PGEKS[1]
  #Dissimilarity
Delta <- matrix(0,nrow=7,ncol=7)
for(y in 1:7){
    for(z in 1:7){
      yy <- sum(pp[y,]*qq[y,])
      zy <- sum(pp[z,]*qq[y,])
      zz <- sum(pp[z,]*qq[z,])
      yz <- sum(pp[y,]*qq[z,])
      Delta1 <- 0
      Delta2 <- 0
      for(nn in 1:402){
        Delta1 <- Delta1 + (pp[y,nn]*qq[y,nn]/yy-pp[z,nn]*qq[y,nn]/zy)**2
        Delta2 <- Delta2 + (pp[z,nn]*qq[z,nn]/zz-pp[y,nn]*qq[z,nn]/yz)**2
      }
      Delta[y,z] <- Delta1+Delta2
    }
}
  #Similarity Linked Index
PS <- 1
PS <- c(PS,PF[2])
for(yy in 3:7){
    tt <- which(Delta[yy,1:yy-1]==min(Delta[yy,1:yy-1]))
    tmp1 <- ifelse(pp[yy,]>0,1,0)
    tmp2 <- ifelse(pp[tt,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    PLCtmp <- (sum(pp[yy,overlap]*qq[tt,overlap])/sum(pp[tt,overlap]*qq[tt,overlap]))
    PPCtmp <- (sum(pp[yy,overlap]*qq[yy,overlap])/sum(pp[tt,overlap]*qq[yy,overlap]))
    PS <- c(PS,PS[tt]*(PLCtmp*PPCtmp)^(1/2))
}
ssec3idx <- rbind(ssec3idx,cbind(ssec3idx[dim(ssec3idx)[1],1]*PL,
                                 ssec3idx[dim(ssec3idx)[1],2]*PP,
                                 ssec3idx[dim(ssec3idx)[1],3]*PF,
                                 ssec3idx[dim(ssec3idx)[1],4]*PT,
                                 ssec3idx[dim(ssec3idx)[1],5]*PLC,
                                 ssec3idx[dim(ssec3idx)[1],6]*PPC,
                                 ssec3idx[dim(ssec3idx)[1],7]*PFC,
                                 ssec3idx[dim(ssec3idx)[1],8]*PTC,
                                 ssec3idx[dim(ssec3idx)[1],9]*PGEKS,
                                 ssec3idx[dim(ssec3idx)[1],10]*PS))

write.csv(ssec3idx,"ssec3idx.csv")
x11(width = 8, height = 4)
plot(PL, type='l',col='blue', xlab='time', ylab='indexes', main = 'maximum overlap year over year monthly indexes')
lines(PP, col='red')
lines(PF, col='green')
lines(PT, col='grey')
lines(PLC, col='black')
lines(PPC, col='pink')
lines(PFC, col='yellow')
lines(PTC, col='brown')
lines(PGEKS, col='orange')
lines(PS, col='#FF0090')
legend('topleft', legend=c('PL','PP','PF','PT','PLC','PPC','PFC','PTC','PGEKS','PS'),col=c('blue','red','green','grey','black','pink','yellow','brown','orange','#FF0090'),lty=1, ncol=5)



## 5.The Construction of Annual Indexes using Carry Forward Prices P14 SEC4
pp <- p[,-c(403)]
qq <- q[,-c(403)]
#imputeprice
mlist <- unique(p$month)
ppp <- NULL
qqq <- NULL
for(mm in 1:12){
  pp <- p[p$month==mlist[mm],-c(403)]
  qq <- q[q$month==mlist[mm],-c(403)]
  for(nn in 1:402){
    if(length(which(pp[,nn]==0))!=0&length(which(pp[,nn]==0))!=7){
      for(yy in 1:7){
        if(yy==1&pp[yy,nn]==0){
          pp[yy,nn] <- pp[min(which(pp[,nn]!=0)),nn]
        }else if(pp[yy,nn]==0){
          tmp <- which(pp[,nn]!=0)
          pp[yy,nn] <- pp[max(tmp[tmp<yy]),nn]
        }
      }
    }
  }
  ppp <- rbind(ppp,pp)
  qqq <- rbind(qqq,qq)
}
#Laspeyres,Paasche
PL <- NULL
PP <- NULL
for(yy in 1:7){
  PLN <- 0
  PLD <- 0
  PPN <- 0
  PPD <- 0
  for(mm in 1:12){
    PLN <- PLN + sum(ppp[(yy-1)+(mm-1)*7+1,] * qqq[(mm-1)*7+1,])
    PLD <- PLD + sum(ppp[(mm-1)*7+1,] * qqq[(mm-1)*7+1,])
    PPN <- PPN + sum(ppp[(mm-1)*7+1,] * qqq[(yy-1)+(mm-1)*7+1,])
    PPD <- PPD + sum(ppp[(yy-1)+(mm-1)*7+1,] * qqq[(yy-1)+(mm-1)*7+1,])
  }
  PL <- c(PL,PLN/PLD)
  PP <- c(PP,PPD/PPN)
}
PF <- (PL*PP)^(1/2)
cbind(PL,PP,PF)
PLC <- 1
PPC <- 1
for(yy in 2:7){
  PLN <- 0
  PLD <- 0
  PPN <- 0
  PPD <- 0
  for(mm in 1:12){
    PLN <- PLN + sum(ppp[(yy-1)+(mm-1)*7+1,] * qqq[(yy-2)+(mm-1)*7+1,])
    PLD <- PLD + sum(ppp[(yy-2)+(mm-1)*7+1,] * qqq[(yy-2)+(mm-1)*7+1,])
    PPN <- PPN + sum(ppp[(yy-2)+(mm-1)*7+1,] * qqq[(yy-1)+(mm-1)*7+1,])
    PPD <- PPD + sum(ppp[(yy-1)+(mm-1)*7+1,] * qqq[(yy-1)+(mm-1)*7+1,])
  }
  PLC <- c(PLC,PLC[yy-1]*PLN/PLD)
  PPC <- c(PPC,PPC[yy-1]*PPD/PPN)
}
PFC <- (PLC*PPC)^(1/2)
cbind(PLC,PPC,PFC)
#Tornqvist
PT <- NULL
for(yy in 1:7){
  tmp <- 0
  SyD <- 0
  S1D <- 0
  for(mm in 1:12){
    SyD <- SyD + sum(ppp[(yy-1)+(mm-1)*7+1,]*qqq[(yy-1)+(mm-1)*7+1,])
    S1D <- S1D + sum(ppp[(mm-1)*7+1,]*qqq[(mm-1)*7+1,])
  }
  for(mm in 1:12){
    SyN <- sum(ppp[(yy-1)+(mm-1)*7+1,]*qqq[(yy-1)+(mm-1)*7+1,])
    S1N <- sum(ppp[(mm-1)*7+1,]*qqq[(mm-1)*7+1,])
    for(nn in 1:402){
      if(ppp[(yy-1)+(mm-1)*7+1,nn]!=0&ppp[(mm-1)*7+1,nn]!=0){
        sy <- ppp[(yy-1)+(mm-1)*7+1,nn]*qqq[(yy-1)+(mm-1)*7+1,nn]/sum(ppp[(yy-1)+(mm-1)*7+1,]*qqq[(yy-1)+(mm-1)*7+1,])
        s1 <- ppp[(mm-1)*7+1,nn]*qqq[(mm-1)*7+1,nn]/sum(ppp[(mm-1)*7+1,]*qqq[(mm-1)*7+1,])
        tmp <- tmp + 0.5*(s1*S1N/S1D+sy*SyN/SyD)*log(ppp[(yy-1)+(mm-1)*7+1,nn]/ppp[(mm-1)*7+1,nn])
      }
    }
  }
  PT <- c(PT,exp(tmp))
}
PTC <- 1
for(yy in 2:7){
  tmp <- 0
  SyD <- 0
  S1D <- 0
  for(mm in 1:12){
    SyD <- SyD + sum(ppp[(yy-1)+(mm-1)*7+1,]*qqq[(yy-1)+(mm-1)*7+1,])
    S1D <- S1D + sum(ppp[(yy-2)+(mm-1)*7+1,]*qqq[(yy-2)+(mm-1)*7+1,])
  }
  for(mm in 1:12){
    SyN <- sum(ppp[(yy-1)+(mm-1)*7+1,]*qqq[(yy-1)+(mm-1)*7+1,])
    S1N <- sum(ppp[(yy-2)+(mm-1)*7+1,]*qqq[(yy-2)+(mm-1)*7+1,])
    for(nn in 1:402){
      if(ppp[(yy-1)+(mm-1)*7+1,nn]!=0&ppp[(yy-2)+(mm-1)*7+1,nn]!=0){
        sy <- ppp[(yy-1)+(mm-1)*7+1,nn]*qqq[(yy-1)+(mm-1)*7+1,nn]/sum(ppp[(yy-1)+(mm-1)*7+1,]*qqq[(yy-1)+(mm-1)*7+1,])
        s1 <- ppp[(yy-2)+(mm-1)*7+1,nn]*qqq[(yy-2)+(mm-1)*7+1,nn]/sum(ppp[(yy-2)+(mm-1)*7+1,]*qqq[(yy-2)+(mm-1)*7+1,])
        tmp <- tmp + 0.5*(s1*S1N/S1D+sy*SyN/SyD)*log(ppp[(yy-1)+(mm-1)*7+1,nn]/ppp[(yy-2)+(mm-1)*7+1,nn])
      }
    }
  }
  PTC <- c(PTC,PTC[yy-1]*exp(tmp))
}
#GEKS
PGEKS <- NULL
for(yy in 1:7){
  tmp <- 1
  for(zz in 1:7){
    PLN <- 0
    PLD <- 0
    PPN <- 0
    PPD <- 0
    for(mm in 1:12){
      PLN <- PLN + sum(ppp[(yy-1)+(mm-1)*7+1,] * qqq[(zz-1)+(mm-1)*7+1,])
      PLD <- PLD + sum(ppp[(zz-1)+(mm-1)*7+1,] * qqq[(zz-1)+(mm-1)*7+1,])
      PPN <- PPN + sum(ppp[(zz-1)+(mm-1)*7+1,] * qqq[(yy-1)+(mm-1)*7+1,])
      PPD <- PPD + sum(ppp[(yy-1)+(mm-1)*7+1,] * qqq[(yy-1)+(mm-1)*7+1,])
    }
    tmpPL <- PLN/PLD
    tmpPP <- PPD/PPN
    tmp <- tmp * (tmpPL*tmpPP)^(1/2)
  }
  PGEKS <- c(PGEKS,tmp^(1/7))
}
PGEKS <- PGEKS / PGEKS[1]
#CCDI
PCCDI <- c(1,1,1,1,1,1,1)
for(yy in 1:7){
  for(zz in 1:7){
    tmp <- 0
    SyD <- 0
    SzD <- 0
    for(mm in 1:12){
      SyD <- SyD + sum(ppp[(yy-1)+(mm-1)*7+1,]*qqq[(yy-1)+(mm-1)*7+1,])
      SzD <- SzD + sum(ppp[(zz-1)+(mm-1)*7+1,]*qqq[(zz-1)+(mm-1)*7+1,])
    }
    for(mm in 1:12){
      SyN <- sum(ppp[(yy-1)+(mm-1)*7+1,]*qqq[(yy-1)+(mm-1)*7+1,])
      SzN <- sum(ppp[(zz-1)+(mm-1)*7+1,]*qqq[(zz-1)+(mm-1)*7+1,])
      for(nn in 1:402){
        if(ppp[(yy-1)+(mm-1)*7+1,nn]!=0&ppp[(zz-1)+(mm-1)*7+1,nn]!=0){
          sy <- ppp[(yy-1)+(mm-1)*7+1,nn]*qqq[(yy-1)+(mm-1)*7+1,nn]/sum(ppp[(yy-1)+(mm-1)*7+1,]*qqq[(yy-1)+(mm-1)*7+1,])
          s1 <- ppp[(zz-1)+(mm-1)*7+1,nn]*qqq[(zz-1)+(mm-1)*7+1,nn]/sum(ppp[(zz-1)+(mm-1)*7+1,]*qqq[(zz-1)+(mm-1)*7+1,])
          tmp <- tmp + 0.5*(s1*SzN/SzD+sy*SyN/SyD)*log(ppp[(yy-1)+(mm-1)*7+1,nn]/ppp[(zz-1)+(mm-1)*7+1,nn])
        }
      }
    }
    PCCDI[yy] <- PCCDI[yy]*exp(tmp)
  }
  PCCDI[yy] <- PCCDI[yy]^(1/7)
}
PCCDI <- PCCDI/PCCDI[1]
#Dissimilarity
Delta <- matrix(0,nrow=7,ncol=7)

for(y in 1:7){
  for(z in 1:7){
    Delta1 <- 0
    Delta2 <- 0
    for(mm in 1:12){
      yy <- sum(ppp[(y-1)+(mm-1)*7+1,]*qqq[(y-1)+(mm-1)*7+1,])
      zy <- sum(ppp[(z-1)+(mm-1)*7+1,]*qqq[(y-1)+(mm-1)*7+1,])
      zz <- sum(ppp[(z-1)+(mm-1)*7+1,]*qqq[(z-1)+(mm-1)*7+1,])
      yz <- sum(ppp[(y-1)+(mm-1)*7+1,]*qqq[(z-1)+(mm-1)*7+1,])
      for(nn in 1:402){
        Delta1 <- Delta1 + (ppp[(y-1)+(mm-1)*7+1,nn]*qqq[(y-1)+(mm-1)*7+1,nn]/yy-ppp[(z-1)+(mm-1)*7+1,nn]*qqq[(y-1)+(mm-1)*7+1,nn]/zy)**2
        Delta2 <- Delta2 + (ppp[(z-1)+(mm-1)*7+1,nn]*qqq[(z-1)+(mm-1)*7+1,nn]/zz-ppp[(y-1)+(mm-1)*7+1,nn]*qqq[(z-1)+(mm-1)*7+1,nn]/yz)**2
      }
    }
    Delta[y,z] <- Delta1+Delta2
  }
}
PS <- c(1,PF[2])
for(yy in 3:7){
  PLN <- 0
  PLD <- 0
  PPN <- 0
  PPD <- 0
  tt <- which(Delta[yy,1:yy-1]==min(Delta[yy,1:yy-1]))
  for(mm in 1:12){
    PLN <- PLN + sum(ppp[(yy-1)+(mm-1)*7+1,] * qqq[(tt-1)+(mm-1)*7+1,])
    PLD <- PLD + sum(ppp[(tt-1)+(mm-1)*7+1,] * qqq[(tt-1)+(mm-1)*7+1,])
    PPN <- PPN + sum(ppp[(tt-1)+(mm-1)*7+1,] * qqq[(yy-1)+(mm-1)*7+1,])
    PPD <- PPD + sum(ppp[(yy-1)+(mm-1)*7+1,] * qqq[(yy-1)+(mm-1)*7+1,])
  }
  PLtmp <- PLN/PLD
  PPtmp <- PPD/PPN
  PFtmp <- (PLtmp*PPtmp)^(1/2)
  PS <- c(PS,PS[tt]*PFtmp)
}
ssec4idx <- cbind(PL,PLC,PP,PPC,PF,PFC,PT,PTC,PGEKS,PCCDI,PS)
write.csv(ssec4idx,"./ssec4idx.csv")
x11(width = 8, height = 4)

plot(PL, type='l',col='blue', xlab='time', ylab='indexes', main = 'annual indexes using carry forward prices')
lines(PP, col='red')
lines(PF, col='green')
lines(PT, col='grey')
lines(PLC, col='black')
lines(PPC, col='pink')
lines(PFC, COL='yellow')
lines(PTC, COL='brown')
lines(PGEKS, COL='orange')
lines(PS, COL='#FF0090')
lines(PCCDI, COL='#F000F0')
legend('topleft', legend=c('PL','PP','PF','PT','PLC','PPC','PFC','PTC','PGEKS','PS','PCCDI'),col=c('blue','red','green','grey','black','pink','yellow','brown','orange','#FF0090','#F000F0'),lty=1, ncol=5)


#6.The Construction of Annnual indexes using Maximum overlap Bilateral Indexes P18
pp <- p[,-c(403)]
qq <- q[,-c(403)]

ppp <- NULL
qqq <- NULL
for(mm in 1:12){
  for(yy in 1:7){
    ppp <- rbind(ppp,pp[(yy-1)*12+(mm-1)+1,])
    qqq <- rbind(qqq,qq[(yy-1)*12+(mm-1)+1,])
  }
}

# Laspeyres,Paasche,Fisher
PL <- NULL
PP <- NULL
for(yy in 1:7){
  PLN <- 0
  PLD <- 0
  PPN <- 0
  PPD <- 0
  for(mm in 1:12){
    tmp1 <- ifelse(ppp[(yy-1)+(mm-1)*6+1,]>0,1,0)
    tmp2 <- ifelse(ppp[(mm-1)*6+1,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    PLN <- PLN + sum(ppp[(yy-1)+(mm-1)*6+1,overlap] * qqq[(mm-1)*6+1,overlap])
    PLD <- PLD + sum(ppp[(mm-1)*6+1,overlap] * qqq[(mm-1)*6+1,overlap])
    PPN <- PPN + sum(ppp[(mm-1)*6+1,overlap] * qqq[(yy-1)+(mm-1)*6+1,overlap])
    PPD <- PPD + sum(ppp[(yy-1)+(mm-1)*6+1,overlap] * qqq[(yy-1)+(mm-1)*6+1,overlap])
  }
  PL <- c(PL,PLN/PLD)
  PP <- c(PP,PPD/PPN)
}
PF <- (PL*PP)^(1/2)
cbind(PL,PP,PF)
# chain Laspeyres, chain Paasche
PLC <- 1
PPC <- 1
for(yy in 2:7){
  PLN <- 0
  PLD <- 0
  PPN <- 0
  PPD <- 0
  for(mm in 1:12){
    tmp1 <- ifelse(ppp[(yy-1)+(mm-1)*6+1,]>0,1,0)
    tmp2 <- ifelse(ppp[(yy-2)+(mm-1)*6+1,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    PLN <- PLN + sum(ppp[(yy-1)+(mm-1)*6+1,overlap] * qqq[(yy-2)+(mm-1)*6+1,overlap])
    PLD <- PLD + sum(ppp[(yy-2)+(mm-1)*6+1,overlap] * qqq[(yy-2)+(mm-1)*6+1,overlap])
    PPN <- PPN + sum(ppp[(yy-2)+(mm-1)*6+1,overlap] * qqq[(yy-1)+(mm-1)*6+1,overlap])
    PPD <- PPD + sum(ppp[(yy-1)+(mm-1)*6+1,overlap] * qqq[(yy-1)+(mm-1)*6+1,overlap])
  }
  PLC <- c(PLC,PLC[yy-1]*PLN/PLD)
  PPC <- c(PPC,PPC[yy-1]*PPD/PPN)
}
PFC <- (PLC*PPC)^(1/2)
cbind(PLC,PPC,PFC)
#Tornqvust
PT <- NULL
for(yy in 1:6){
  tmp <- 0
  SyD <- 0
  S1D <- 0
  for(mm in 1:12){
    tmp1 <- ifelse(ppp[(yy-1)+(mm-1)*6+1,]>0,1,0)
    tmp2 <- ifelse(ppp[(mm-1)*6+1,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    SyD <- SyD + sum(ppp[(yy-1)+(mm-1)*6+1,overlap]*qqq[(yy-1)+(mm-1)*6+1,overlap])
    S1D <- S1D + sum(ppp[(mm-1)*6+1,overlap]*qqq[(mm-1)*6+1,overlap])
  }
  for(mm in 1:12){
    tmp1 <- ifelse(ppp[(yy-1)+(mm-1)*6+1,]>0,1,0)
    tmp2 <- ifelse(ppp[(mm-1)*6+1,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    SyN <- sum(ppp[(yy-1)+(mm-1)*6+1,overlap]*qqq[(yy-1)+(mm-1)*6+1,overlap])
    S1N <- sum(ppp[(mm-1)*6+1,overlap]*qqq[(mm-1)*6+1,overlap])
    for(nn in overlap){
      if(ppp[(yy-1)+(mm-1)*6+1,nn]!=0&ppp[(mm-1)*6+1,nn]!=0){
        sy <- ppp[(yy-1)+(mm-1)*6+1,nn]*qqq[(yy-1)+(mm-1)*6+1,nn]/sum(ppp[(yy-1)+(mm-1)*6+1,overlap]*qqq[(yy-1)+(mm-1)*6+1,overlap])
        s1 <- ppp[(mm-1)*6+1,nn]*qqq[(mm-1)*6+1,nn]/sum(ppp[(mm-1)*6+1,overlap]*qqq[(mm-1)*6+1,overlap])
        tmp <- tmp + 0.5*(s1*S1N/S1D+sy*SyN/SyD)*log(ppp[(yy-1)+(mm-1)*6+1,nn]/ppp[(mm-1)*6+1,nn])
      }
    }
  }
  PT <- c(PT,exp(tmp))
}
#ChainTornqvust
PTC <- 1
for(yy in 2:6){
  tmp <- 0
  SyD <- 0
  S1D <- 0
  for(mm in 1:12){
    tmp1 <- ifelse(ppp[(yy-1)+(mm-1)*6+1,]>0,1,0)
    tmp2 <- ifelse(ppp[(yy-2)+(mm-1)*6+1,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    SyD <- SyD + sum(ppp[(yy-1)+(mm-1)*6+1,overlap]*qqq[(yy-1)+(mm-1)*6+1,overlap])
    S1D <- S1D + sum(ppp[(yy-2)+(mm-1)*6+1,overlap]*qqq[(yy-2)+(mm-1)*6+1,overlap])
  }
  for(mm in 1:12){
    tmp1 <- ifelse(ppp[(yy-1)+(mm-1)*6+1,]>0,1,0)
    tmp2 <- ifelse(ppp[(yy-2)+(mm-1)*6+1,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    SyN <- sum(ppp[(yy-1)+(mm-1)*6+1,overlap]*qqq[(yy-1)+(mm-1)*6+1,overlap])
    S1N <- sum(ppp[(yy-2)+(mm-1)*6+1,overlap]*qqq[(yy-2)+(mm-1)*6+1,overlap])
    for(nn in overlap){
      if(ppp[(yy-1)+(mm-1)*6+1,nn]!=0&ppp[(yy-2)+(mm-1)*6+1,nn]!=0){
        sy <- ppp[(yy-1)+(mm-1)*6+1,nn]*qqq[(yy-1)+(mm-1)*6+1,nn]/sum(ppp[(yy-1)+(mm-1)*6+1,overlap]*qqq[(yy-1)+(mm-1)*6+1,overlap])
        s1 <- ppp[(yy-2)+(mm-1)*6+1,nn]*qqq[(yy-2)+(mm-1)*6+1,nn]/sum(ppp[(yy-2)+(mm-1)*6+1,overlap]*qqq[(yy-2)+(mm-1)*6+1,overlap])
        tmp <- tmp + 0.5*(s1*S1N/S1D+sy*SyN/SyD)*log(ppp[(yy-1)+(mm-1)*6+1,nn]/ppp[(yy-2)+(mm-1)*6+1,nn])
      }
    }
  }
  PTC <- c(PTC,PTC[yy-1]*exp(tmp))
}
#PGEKS
PGEKS <- NULL
for(yy in 1:6){
  tmp <- 1
  for(zz in 1:6){
    PLN <- 0
    PLD <- 0
    PPN <- 0
    PPD <- 0
    for(mm in 1:12){
      tmp1 <- ifelse(ppp[(yy-1)+(mm-1)*6+1,]>0,1,0)
      tmp2 <- ifelse(ppp[(zz-1)+(mm-1)*6+1,]>0,1,0)
      overlap <- which(tmp1*tmp2==1)
      PLN <- PLN + sum(ppp[(yy-1)+(mm-1)*6+1,overlap] * qqq[(zz-1)+(mm-1)*6+1,overlap])
      PLD <- PLD + sum(ppp[(zz-1)+(mm-1)*6+1,overlap] * qqq[(zz-1)+(mm-1)*6+1,overlap])
      PPN <- PPN + sum(ppp[(zz-1)+(mm-1)*6+1,overlap] * qqq[(yy-1)+(mm-1)*6+1,overlap])
      PPD <- PPD + sum(ppp[(yy-1)+(mm-1)*6+1,overlap] * qqq[(yy-1)+(mm-1)*6+1,overlap])
    }
    tmpPL <- PLN/PLD
    tmpPP <- PPD/PPN
    tmp <- tmp * (tmpPL*tmpPP)^(1/2)
  }
  PGEKS <- c(PGEKS,tmp^(1/6))
}
PGEKS <- PGEKS / PGEKS[1]
#CCDI
PCCDI <- c(1,1,1,1,1,1)
for(yy in 1:6){
  for(zz in 1:6){
    tmp <- 0
    SyD <- 0
    S1D <- 0
    for(mm in 1:12){
      tmp1 <- ifelse(ppp[(yy-1)+(mm-1)*6+1,]>0,1,0)
      tmp2 <- ifelse(ppp[(zz-1)+(mm-1)*6+1,]>0,1,0)
      overlap <- which(tmp1*tmp2==1)
      SyD <- SyD + sum(ppp[(yy-1)+(mm-1)*6+1,overlap]*qqq[(yy-1)+(mm-1)*6+1,overlap])
      S1D <- S1D + sum(ppp[(zz-1)+(mm-1)*6+1,overlap]*qqq[(zz-1)+(mm-1)*6+1,overlap])
    }
    for(mm in 1:12){
      tmp1 <- ifelse(ppp[(yy-1)+(mm-1)*6+1,]>0,1,0)
      tmp2 <- ifelse(ppp[(zz-1)+(mm-1)*6+1,]>0,1,0)
      overlap <- which(tmp1*tmp2==1)
      SyN <- sum(ppp[(yy-1)+(mm-1)*6+1,overlap]*qqq[(yy-1)+(mm-1)*6+1,overlap])
      S1N <- sum(ppp[(zz-1)+(mm-1)*6+1,overlap]*qqq[(zz-1)+(mm-1)*6+1,overlap])
      for(nn in overlap){
        if(ppp[(yy-1)+(mm-1)*6+1,nn]!=0&ppp[(zz-1)+(mm-1)*6+1,nn]!=0){
          sy <- ppp[(yy-1)+(mm-1)*6+1,nn]*qqq[(yy-1)+(mm-1)*6+1,nn]/sum(ppp[(yy-1)+(mm-1)*6+1,overlap]*qqq[(yy-1)+(mm-1)*6+1,overlap])
          s1 <- ppp[(zz-1)+(mm-1)*6+1,nn]*qqq[(zz-1)+(mm-1)*6+1,nn]/sum(ppp[(zz-1)+(mm-1)*6+1,overlap]*qqq[(zz-1)+(mm-1)*6+1,overlap])
          tmp <- tmp + 0.5*(s1*S1N/S1D+sy*SyN/SyD)*log(ppp[(yy-1)+(mm-1)*6+1,nn]/ppp[(zz-1)+(mm-1)*6+1,nn])
        }
      }
    }
    PCCDI[yy] <- PCCDI[yy]*exp(tmp)
  }
  PCCDI[yy] <- PCCDI[yy]^(1/6)
}
PCCDI <- PCCDI / PCCDI[1]
#Dissimilarity
Delta <- matrix(0,nrow=6,ncol=6)
for(y in 1:6){
  for(z in 1:6){
    Delta1 <- 0
    Delta2 <- 0
    for(mm in 1:12){
      yy <- sum(ppp[(y-1)+(mm-1)*6+1,]*qqq[(y-1)+(mm-1)*6+1,])
      zy <- sum(ppp[(z-1)+(mm-1)*6+1,]*qqq[(y-1)+(mm-1)*6+1,])
      zz <- sum(ppp[(z-1)+(mm-1)*6+1,]*qqq[(z-1)+(mm-1)*6+1,])
      yz <- sum(ppp[(y-1)+(mm-1)*6+1,]*qqq[(z-1)+(mm-1)*6+1,])
      for(nn in 1:14){
        Delta1 <- Delta1 + (ppp[(y-1)+(mm-1)*6+1,nn]*qqq[(y-1)+(mm-1)*6+1,nn]/yy-ppp[(z-1)+(mm-1)*6+1,nn]*qqq[(y-1)+(mm-1)*6+1,nn]/zy)**2
        Delta2 <- Delta2 + (ppp[(z-1)+(mm-1)*6+1,nn]*qqq[(z-1)+(mm-1)*6+1,nn]/zz-ppp[(y-1)+(mm-1)*6+1,nn]*qqq[(z-1)+(mm-1)*6+1,nn]/yz)**2
      }
    }
    Delta[y,z] <- Delta1+Delta2
  }
}
PS <- c(1,PF[2])
for(yy in 3:6){
  PLN <- 0
  PLD <- 0
  PPN <- 0
  PPD <- 0
  tt <- which(Delta[yy,1:yy-1]==min(Delta[yy,1:yy-1]))
  for(mm in 1:12){
    tmp1 <- ifelse(ppp[(yy-1)+(mm-1)*6+1,]>0,1,0)
    tmp2 <- ifelse(ppp[(tt-1)+(mm-1)*6+1,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    PLN <- PLN + sum(ppp[(yy-1)+(mm-1)*6+1,overlap] * qqq[(tt-1)+(mm-1)*6+1,overlap])
    PLD <- PLD + sum(ppp[(tt-1)+(mm-1)*6+1,overlap] * qqq[(tt-1)+(mm-1)*6+1,overlap])
    PPN <- PPN + sum(ppp[(tt-1)+(mm-1)*6+1,overlap] * qqq[(yy-1)+(mm-1)*6+1,overlap])
    PPD <- PPD + sum(ppp[(yy-1)+(mm-1)*6+1,overlap] * qqq[(yy-1)+(mm-1)*6+1,overlap])
  }
  PLtmp <- PLN/PLD
  PPtmp <- PPD/PPN
  PFtmp <- (PLtmp*PPtmp)^(1/2)
  PS <- c(PS,PS[tt]*PFtmp)
}
ssec5idx <- cbind(PL,PLC,PP,PPC,PF,PFC,PT,PTC,PGEKS,PCCDI,PS)
write.csv(ssec5idx,"./ssec5idx.csv")

x11(width = 8, height = 4)
plot(PL, type='l',col='blue', xlab='time', ylab='indexes', main = 'annnual indexes using maximum overlap bilateral indexes')
lines(PP, col='red')
lines(PF, col='green')
lines(PT, col='grey')
lines(PLC, col='black')
lines(PPC, col='pink')
lines(PFC, COL='yellow')
lines(PTC, COL='brown')
lines(PGEKS, COL='orange')
lines(PS, COL='#FF0090')
lines(PCCDI, COL='#F000F0')
legend('topleft', legend=c('PL','PP','PF','PT','PLC','PPC','PFC','PTC','PGEKS','PS','PCCDI'),col=c('blue','red','green','grey','black','pink','yellow','brown','orange','#FF0090','#F000F0'),lty=1, ncol=5)


## 7.Month to Month Indexers using Carry Forward prices P22 Sec6
pp <- p[,c(1:402)]
qq <- q[,c(1:402)]
#impute
for(nn in 1:402){
  if(length(which(pp[,nn]==0))!=0&length(which(pp[,nn]==0))!=7){
    for(yy in 1:84){
      if(yy==1&pp[yy,nn]==0){
        pp[yy,nn] <- pp[min(which(pp[,nn]!=0)),nn]
      }else if(pp[yy,nn]==0){
        tmp <- which(pp[,nn]!=0)
        pp[yy,nn] <- pp[max(tmp[tmp<yy]),nn]
      }
    }
  }
}
# Laspeyres,Paasche
PL <- NULL
PP <- NULL
for(yy in 1:84){
  PL <- c(PL,sum(pp[yy,]*qq[1,])/sum(pp[1,]*qq[1,]))
  PP <- c(PP,sum(pp[yy,]*qq[yy,])/sum(pp[1,]*qq[yy,]))
}

# Fisher
PF <- (PL*PP)^(1/2)
cbind(PL,PP,PF)

# chain Laspeyres, chain Paasche, chain Fisher
PLC <- 1
PPC <- 1
for(yy in 2:84){
  tmpPL <- sum(pp[yy,]*qq[(yy-1),])/sum(pp[(yy-1),]*qq[(yy-1),])
  tmpPP <- sum(pp[yy,]*qq[yy,])/sum(pp[(yy-1),]*qq[yy,])
  PLC <- c(PLC,PLC[yy-1]*tmpPL)
  PPC <- c(PPC,PPC[yy-1]*tmpPP)
}
PFC <- (PLC*PPC)^(1/2)
cbind(PLC,PPC,PFC)

# GEKS
PGEKS <- NULL
for(yy in 1:84){
  tmp <- 1
  for(rr in 1:84){
    tmpPL <- sum(pp[yy,]*qq[rr,])/sum(pp[rr,]*qq[rr,])
    tmpPP <- sum(pp[yy,]*qq[yy,])/sum(pp[rr,]*qq[yy,])
    tmp <- tmp * (tmpPL*tmpPP)^(1/2)
  }
  PGEKS <- c(PGEKS,tmp^(1/84))
}
PGEKS <- PGEKS / PGEKS[1]

# Dissimilarity
Delta <- matrix(0,nrow=84,ncol=84)
for(t in 1:84){
  for(r in 1:84){
    tt <- sum(pp[t,]*qq[t,])
    rt <- sum(pp[r,]*qq[t,])
    rr <- sum(pp[r,]*qq[r,])
    tr <- sum(pp[t,]*qq[r,])
    Delta1 <- 0
    Delta2 <- 0
    for(nn in 1:402){
      Delta1 <- Delta1 + (pp[t,nn]*qq[t,nn]/tt-pp[r,nn]*qq[t,nn]/rt)**2
      Delta2 <- Delta2 + (pp[r,nn]*qq[r,nn]/rr-pp[t,nn]*qq[r,nn]/tr)**2
    }
    Delta[t,r] <- Delta1 + Delta2    
  }
}
PS <- c(1,PF[2])
for(yy in 3:84){
  tt <- which(Delta[yy,1:yy-1]==min(Delta[yy,1:yy-1]))
  tmpPL <- sum(pp[yy,]*qq[tt,])/sum(pp[tt,]*qq[tt,])
  tmpPP <- sum(pp[yy,]*qq[yy,])/sum(pp[tt,]*qq[yy,])
  PS <- c(PS,PS[tt]*(tmpPL*tmpPP)^(1/2))
}
ssec6idx <- cbind(PL,PLC,PP,PPC,PFC,PF,PGEKS,PS)
write.csv(ssec6idx,"./ssec6idx.csv")

x11(width = 8, height = 4)
plot(PLC, type='l',col='blue', xlab='time', ylab='indexes', main = 'month to month indexes using carry forward prices')
lines(PP, col='red')
lines(PF, col='green')
lines(PL, col='black')
lines(PPC, col='pink')
lines(PFC, COL='yellow')
lines(PGEKS, COL='orange')
lines(PS, COL='#FF0090')
legend('topleft', legend=c('PLC','PP','PF','PL','PPC','PFC','PGEKS','PS'),col=c('blue','red','green','grey','black','pink','yellow','orange','#FF0090'),lty=1, ncol=5)

## 8. Month to Month Indexes using Maximum Overlap Bilateral Indexes P24
pp <- p[,c(1:402)]
qq <- q[,c(1:402)]
PL <- NULL
PP <- NULL
for(ii in 1:84){
  tmp1 <- ifelse(pp[ii,]>0,1,0)
  tmp2 <- ifelse(pp[1,]>0,1,0)
  overlap <- which(tmp1*tmp2==1)
  PL <- c(PL,sum(pp[ii,overlap]*qq[1,overlap])/sum(pp[1,overlap]*qq[1,overlap]))
  PP <- c(PP,sum(pp[ii,overlap]*qq[ii,overlap])/sum(pp[1,overlap]*qq[ii,overlap]))
}
print(cbind(PL,PP))
PF <- sqrt(PP*PL)
print(PF)
PLC <- 1
PPC <- 1
for(ii in 2:84){
  tmp1 <- ifelse(pp[ii,]>0,1,0)
  tmp2 <- ifelse(pp[ii-1,]>0,1,0)
  overlap <- which(tmp1*tmp2==1)
  PLC <- c(PLC,PLC[ii-1]*(sum(pp[ii,overlap]*qq[ii-1,overlap])/sum(pp[ii-1,overlap]*qq[ii-1,overlap])))
  PPC <- c(PPC,PPC[ii-1]*(sum(pp[ii,overlap]*qq[ii,overlap])/sum(pp[ii-1,overlap]*qq[ii,overlap])))
}
PFC <- sqrt(PLC*PPC)
cbind(PLC,PPC,PFC)

# PGEKS
PGEKS <- NULL
for(yy in 1:84){
  tmp <- 1
  for(rr in 1:84){
    tmp1 <- ifelse(pp[yy,]>0,1,0)
    tmp2 <- ifelse(pp[rr,]>0,1,0)
    overlap <- which(tmp1*tmp2==1)
    tmpPL <- sum(pp[yy,overlap]*qq[rr,overlap])/sum(pp[rr,overlap]*qq[rr,overlap])
    tmpPP <- sum(pp[yy,overlap]*qq[yy,overlap])/sum(pp[rr,overlap]*qq[yy,overlap])
    tmp <- tmp * (tmpPL*tmpPP)^(1/2)
  }
  PGEKS <- c(PGEKS,tmp^(1/84))
}
PGEKS <- PGEKS/PGEKS[1]

# Dissimilarity
Delta <- matrix(0,nrow=84,ncol=84)
for(t in 1:84){
  for(r in 1:84){
    tt <- sum(pp[t,]*qq[t,])
    rt <- sum(pp[r,]*qq[t,])
    rr <- sum(pp[r,]*qq[r,])
    tr <- sum(pp[t,]*qq[r,])
    Delta1 <- 0
    Delta2 <- 0
    for(nn in 1:402){
      Delta1 <- Delta1 + (pp[t,nn]*qq[t,nn]/tt-pp[r,nn]*qq[t,nn]/rt)**2
      Delta2 <- Delta2 + (pp[r,nn]*qq[r,nn]/rr-pp[t,nn]*qq[r,nn]/tr)**2
    }
    Delta[t,r] <- Delta1 + Delta2    
  }
}
PS <- c(1,PF[2])
for(yy in 3:84){
  tt <- which(Delta[yy,1:yy-1]==min(Delta[yy,1:yy-1]))
  tmp1 <- ifelse(pp[yy,]>0,1,0)
  tmp2 <- ifelse(pp[tt,]>0,1,0)
  overlap <- which(tmp1*tmp2==1)
  tmpPL <- sum(pp[yy,overlap]*qq[tt,overlap])/sum(pp[tt,overlap]*qq[tt,overlap])
  tmpPP <- sum(pp[yy,overlap]*qq[yy,overlap])/sum(pp[tt,overlap]*qq[yy,overlap])
  PS <- c(PS,PS[tt]*(tmpPL*tmpPP)^(1/2))
}
ssec7idx <- cbind(PL,PLC,PP,PPC,PFC,PF,PGEKS,PS)
write.csv(ssec7idx,"./ssec7idx.csv")

x11(width = 8, height = 4)
plot(PL, type='l',col='blue', xlab='time', ylab='indexes', main = 'month to month indexes using maximum overlap bilateral indexes')
lines(PP, col='red')
lines(PF, col='green')
lines(PLC, col='black')
lines(PPC, col='pink')
lines(PFC, COL='yellow')
lines(PGEKS, COL='orange')
lines(PS, COL='#FF0090')
legend('topleft', legend=c('PL','PP','PF','PLC','PPC','PFC','PGEKS','PS'),col=c('blue','red','green','grey','black','pink','yellow','orange','#FF0090'),lty=1, ncol=5)

## 9.Month to Month Unweighted Prices Indexes Using Carry Forward prices P26 sec8
#impute
for(nn in 1:402){
  if(length(which(pp[,nn]==0))!=0&length(which(pp[,nn]==0))!=7){
    for(yy in 1:84){
      if(yy==1&pp[yy,nn]==0){
        pp[yy,nn] <- pp[min(which(pp[,nn]!=0)),nn]
      }else if(pp[yy,nn]==0){
        tmp <- which(pp[,nn]!=0)
        pp[yy,nn] <- pp[max(tmp[tmp<yy]),nn]
      }
    }
  }
}
#Dutot,Carli,Jevons
PD <- NULL
PC <- NULL
PJ <- NULL
for(tt in 1:84){
  PD <- c(PD,sum(pp[tt,])/sum(pp[1,]))
  PC <- c(PC,sum(pp[tt,]/pp[1,])/402)
  PJ <- c(PJ,prod(pp[tt,]/pp[1,])^(1/402))
}
#chain carli
PCC <- 1
for(tt in 2:84){
  PCC <- c(PCC,PCC[tt-1]*sum(pp[tt,]/pp[tt-1,])/402)
}
ssec8idx <- cbind(PJ,PD,PC,PCC,PGEKS)
write.csv(ssec8idx,"ssec8idx.csv")

x11(width = 8, height = 4)
plot(PJ, type='l',col='blue', xlab='time', ylab='indexes', main = 'month to month unweighted prices indexes using carry forward prices')
lines(PD, col='red')
lines(PC, col='green')
lines(PCC, col='black')
lines(PGEKS, col='pink')
legend('topleft', legend=c('PL','PD','PD','PCC','PGEKS'),col=c('blue','red','green','grey','black','pink'),lty=1, ncol=3)

## 10.  Month to Month Unweighted Prices Indexes Using Maximum Overlap Bilateral Indexes. P27 SEC9
#Dutot,Carli,Jevons
PD <- NULL
PC <- NULL
PJ <- NULL
for(tt in 1:84){
  tmp1 <- ifelse(pp[tt,]>0,1,0)
  tmp2 <- ifelse(pp[1,]>0,1,0)
  overlap <- which(tmp1*tmp2==1)
  PD <- c(PD,sum(pp[tt,overlap])/sum(pp[1,overlap]))
  PC <- c(PC,sum(pp[tt,overlap]/pp[1,overlap])/length(overlap))
  PJ <- c(PJ,prod(pp[tt,overlap]/pp[1,overlap])^(1/length(overlap)))
}
PDC <- 1
PCC <- 1
PJC <- 1
for(tt in 2:84){
  tmp1 <- ifelse(pp[tt,]>0,1,0)
  tmp2 <- ifelse(pp[tt-1,]>0,1,0)
  overlap <- which(tmp1*tmp2==1)
  PDC <- c(PDC,PDC[tt-1]*sum(pp[tt,overlap])/sum(pp[tt-1,overlap]))
  PCC <- c(PCC,PCC[tt-1]*sum(pp[tt,overlap]/pp[tt-1,overlap])/length(overlap))
  PJC <- c(PJC,prod(pp[tt,overlap]/pp[tt-1,overlap])^(1/length(overlap)))
}
ssec9idx <- cbind(ssec8idx,PD,PC,PJ,PDC,PCC,PJC)
write.csv(ssec9idx,"./ssec9idx.csv")

x11(width = 8, height = 4)
plot(PD, type='l',col='blue', xlab='time', ylab='indexes', main = 'month to month unweighted indexes Using maximum overlap bilateral indexes')
lines(PC, col='red')
lines(PJ, col='green')
lines(PDC, col='black')
lines(PCC, col='pink')
lines(PJC, COL='yellow')
legend('topleft', legend=c('PD','PC','PJ','PDC','PCC','PJC'),col=c('blue','red','green','grey','black','pink','yellow'),lty=1, ncol=3)

