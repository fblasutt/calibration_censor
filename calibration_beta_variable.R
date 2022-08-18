################################################
## Calibration of the Church Censorship Model
## with time varying beta
################################################

###PART-1: LOAD PACKAGES########
## remove (almost) everything in the working environment.
rm(list = ls())

#increase memory size
memory.limit(size=50000)

# to import xl files
library(readxl)

#to write excel files
require(writexl)

# to output towards latex
require(stargazer)

# for generating random variables
library(evd)

# for latex symbols in the graphs
library(latex2exp)

# genetic algorithm for estimation
library(GA)

# for a bootstraped CI
library(boot)

#extract from a word
library(tidyverse)

#For fancy graphs
library(tikzDevice)
library(ggplot2)

#For list sving
library(rlist)

#Set seed
set.seed(2)
###PART 0: WHICH SAMPLE/MODEL DO YOU WANT?#############################################


impbeta <<-1    #imperfectSqzr application of censorship-max (5./6.)
maxper<<-5      #max number of periods to consider




###PART 1: GET THE DATA#############################################

  #Load results
  setwd("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\data_work\\Fabio")
  mbau<-list.load("data.Rdata")
  
  # for (j in 1:30){print(max(momentss[[j]]-mbau[[j]]))}
  
  
  
  
  
  assign(names(mbau)[1],mbau[[1]])
  assign(names(mbau)[2],mbau[[2]])
  assign(names(mbau)[3],mbau[[3]])
  assign(names(mbau)[4],mbau[[4]])
  assign(names(mbau)[5],mbau[[5]])
  assign(names(mbau)[6],mbau[[6]])
  assign(names(mbau)[7],mbau[[7]])
  assign(names(mbau)[8],mbau[[8]])
  assign(names(mbau)[9],mbau[[9]])
  assign(names(mbau)[10],mbau[[10]])
  assign(names(mbau)[11],mbau[[11]])
  assign(names(mbau)[12],mbau[[12]])
  assign(names(mbau)[13],mbau[[13]])
  assign(names(mbau)[14],mbau[[14]])
  assign(names(mbau)[15],mbau[[15]])
  assign(names(mbau)[16],mbau[[16]])
  assign(names(mbau)[17],mbau[[17]])
  assign(names(mbau)[18],mbau[[18]])
  assign(names(mbau)[19],mbau[[19]])
  assign(names(mbau)[20],mbau[[20]])
  assign(names(mbau)[21],mbau[[21]])
  assign(names(mbau)[22],mbau[[22]])
  assign(names(mbau)[23],mbau[[23]])
  assign(names(mbau)[24],mbau[[24]])
  assign(names(mbau)[25],mbau[[25]])
  assign(names(mbau)[26],mbau[[26]])
  assign(names(mbau)[27],mbau[[27]])
  assign(names(mbau)[28],mbau[[28]])
  assign(names(mbau)[29],mbau[[29]])
  assign(names(mbau)[30],mbau[[30]])
  



###PART 2: ESTIMATION#################################################


t<-numeric()
v<-numeric()
pe<-numeric()
z0<-numeric()
up<-numeric()
Sqc<-numeric(length=5)
Sqr<-numeric(length=5)
SNC<-numeric(length=5)
SNC75<-numeric(length=5)
Sq<-numeric(length=5)
pre<-numeric()
pre1<-numeric()
pre2<-numeric()
i<-numeric()
thet<-numeric()
eqrt<-numeric()
eqt<-numeric()
eqrtt<-numeric()
eqct<-numeric()
# match initial conditions


# some functions to get the simulated variables

f<-function(x,the,k,pr,eqt,eqrt) {
  
  -eqt+(((eqrt/x)^(1/the))/((eqrt/x)^(1/the)+pr))*eqrt+
    (1-(((eqrt/x)^(1/the))/((eqrt/x)^(1/the)+pr)))*x}

fm<-function(x,the,k,pr,eqt,eqrt) {
  
  -eqt+(((eqrt/x)^(1/the))/((eqrt/x)^(1/the)+pr))*eqrt+
    (1-(((eqrt/x)^(1/the))/((eqrt/x)^(1/the)+pr)))*x}

fq<-function(the,pr,eqct,eqrt) {
  
  zeta=(eqrt/eqct)^(1/the)
  emma=zeta/(zeta+pr)
  
  (emma)*eqrt+(1-emma)*eqct}


zt<-function(beta,pe,t,up,thet,eqct,eqrtt) {
  
  
  z0=(eqrtt/eqct)^(1/thet)
  if(t>2){
    z0=pe/(1)*(z0*(1)/pe)^(2)
    pe/(1-beta)*(z0*(1-beta)/pe)^(2^(t-2))}else if(t==2){pe/(1)*(z0*(1)/pe)^(2^(t-1))}else{z0}}


mt<-function(beta,pr,t,up,thet,eqct,eqrtt) zt(beta,pr,t,up,thet,eqct,eqrtt)/(pr+zt(beta,pr,t,up,thet,eqct,eqrtt))

mbetat<-function(beta,pr,t,up,thet,eqct,eqrtt) beta*impbeta*mt(beta,pr,t,up,thet,eqct,eqrtt)



mu<-function(t)ifelse(t==1,200,ifelse(t==2,87.8,ifelse(t==3,78.7,ifelse(t==4,82.8,ifelse(t==5,85.1,0)))))/100





qrs<-function(beta,pe,t,up,pre,nu,thet,eqct,eqrtt) if(t<=2){gamma(1-thet)*((1+nu)*mu(t)*((1)*(pre/gamma(1-thet))^(1/thet)*
                                                                                           (mt(0,pe,t-1,up,thet,eqct,eqrtt))))^thet}else
                                                                                           {gamma(1-thet)*((1+nu)*mu(t)*((1-beta)*(pre/gamma(1-thet))^(1/thet)*
                                                                                                                           (mt(beta,pe,t-1,up,thet,eqct,eqrtt))))^thet} 


qcs<-function(beta,pe,t,up,pre,nu,thet,eqct,eqrtt) if(t<=2){gamma(1-thet)*(((1+nu)*mu(t)*((pre/gamma(1-thet))^(1/thet)*
                                                                                            (1-mt(0,pe,t-1,up,thet,eqct,eqrtt)))))^thet}else
                                                                                            {gamma(1-thet)*(((1+nu)*mu(t)*((pre/gamma(1-thet))^(1/thet)*
                                                                                                                             (1-mt(beta,pe,t-1,up,thet,eqct,eqrtt)))))^thet}             


qs<-function(beta,pe,t,up,pre1,pre2,nu,thet,eqct,eqrtt) if(t<=2){qrs(0,pe,t,up,pre1,nu,thet,eqct,eqrtt)*mt(0,pe,t,up,thet,eqct,eqrtt)+
    (1-mt(0,pe,t,up,thet,eqct,eqrtt))*qcs(0,pe,t,up,pre2,nu,thet,eqct,eqrtt)}else
    {qrs(beta,pe,t,up,pre1,nu,thet,eqct,eqrtt)*mt(beta,pe,t,up,thet,eqct,eqrtt)+
        (1-mt(beta,pe,t,up,thet,eqct,eqrtt))*qcs(beta,pe,t,up,pre2,nu,thet,eqct,eqrtt)}

nc<-function(beta,pe,t,up,pre1,pre2,nu,thet,eqct,eqrtt) if(t<=2){qrs(0,pe,t,up,pre1,nu,thet,eqct,eqrtt)*
    (((1)*mt(0,pe,t,up,thet,eqct,eqrtt))/(1))+
    ((1-mt(0,pe,t,up,thet,eqct,eqrtt))/(1))*
    qcs(0,pe,t,up,pre2,nu,thet,eqct,eqrtt)}else
    {qrs(beta,pe,t,up,pre1,nu,thet,eqct,eqrtt)*
        (((1-beta)*mt(beta,pe,t,up,thet,eqct,eqrtt))/(1-beta*mt(beta,pe,t,up,thet,eqct,eqrtt)))+
        ((1-mt(beta,pe,t,up,thet,eqct,eqrtt))/(1-beta*mt(beta,pe,t,up,thet,eqct,eqrtt)))*
        qcs(beta,pe,t,up,pre2,nu,thet,eqct,eqrtt)}



#GUARDA A INITIAL CONDITION CIU C

# define the function that we want to minimize taking beta
maxg<-100
max<-numeric()
maxx<-numeric()
nuu<-numeric()
betaa<-numeric()
sum<-numeric(length=1)
MT<-numeric()
MR<-numeric()
MT75<-numeric()
MR75<-numeric()
upper<-seq(0.01,200,length.out=maxg)
qu<-numeric()
qur<-numeric()

# compute empirical and non empirical
n<-1000#how many draws from distribution
B=numeric(length=n)
B1=numeric(length=n)
Bc=numeric(length=n)
Br=numeric(length=n)
set.seed(2)
B1=runif(n, min = 0, max = 1)  # shocks for deciding later from which distribution to draw



#parameters for the simulation

mini<-function(b1,b2,b3,b4,pee,nuu,thet,qu,qur){
  
  sum=0     
  pe<-exp(pee)
  Sqr[1]=qur
  Sqc[1]=qu
  Sq[1]=fq(thet,pe,Sqc[1],Sqr[1])
  betaa<-c(0,b1,b2,b3,b4)
  
  
  for(i in 1:maxper){
    if(i>=2){
      
      Sqc[i]<-qcs(betaa[i],pe,i,maxx,Sqc[i-1],nuu,thet,Sqc[1],Sqr[1])
      Sqr[i]<-qrs(betaa[i],pe,i,maxx,Sqr[i-1],nuu,thet,Sqc[1],Sqr[1])
      #Sq[i] <-qs(betaa[i],pe,i,maxx,Sqr[i-1],Sqc[i-1],nuu,thet,Sqc[1],Sqr[1])          
      
    }}
  
  
  
  
  
  for (i in 1:maxper){
    mcur<-mt(betaa[i],pe,i,maxx,thet,Sqc[1],Sqr[1])        
    
    
    
    
    if(i>=2){
      
      sum=sum+((mcur*betaa[i]*impbeta-Embeta[i])/Embeta[i])^2
      #Moments with median
      
      MR[i]<-(((qrs(betaa[i],pe,i,maxx,Sqr[i-1],nuu,thet,Sqc[1],Sqr[1])/(gamma(1-thet)))^(1/1))/
                ((log(2))^thet))
      
      #Model with 75p
      
      MR75[i]<-(((qrs(betaa[i],pe,i,maxx,Sqr[i-1],nuu,thet,Sqc[1],Sqr[1])/(gamma(1-thet)))^(1/1))/
                  ((log(4/3))^thet))
      
      # draw random realization from the two frechet distibutions
      set.seed(2)
      Br=rfrechet(n, loc=0, scale=Sqr[i]/(gamma(1-thet)),shape=1/thet)
      set.seed(2)
      Bc=rfrechet(n, loc=0, scale=Sqc[i]/(gamma(1-thet)),shape=1/thet)
      pr<-mcur#*(1-betaa[i])/(1-betaa[i]*mcur)
      for (j in 1:n) {if (B1[j]>pr) {B[j]=Bc[j]} else {B[j]=Br[j]}}
      MT[i]=median(B)
      MT75[i]=quantile(B,0.75)
      
      #Moments with average
      #MT[i]<-qs(betaa[i],pe,i,maxx,nuu,Sqr[i-1],Sqc[i-1],thet,Sq[1],Sqr[1])
      #MR[i]<-qrs(betaa[i],pe,i,maxx,nuu,Sqr[i-1],thet,Sq[1],Sqr[1])
      
      
      sum=sum+((MT[i]-EqM[i])/EqM[i])^2
      sum=sum+((MT75[i]-EqQ75[i])/EqQ75[i])^2
      #=sum+((MR[i]-EqrM[i])/EqrM[i])^2
      #sum=sum+((MR75[i]-EqrQ75[i])/EqrQ75[i])^2
      
      
      
    }else{
      
      #Median 
      
      MR[1]<-(((Sqr[1]/(gamma(1-thet)))^(1/1))/
                ((log(2))^thet))
      
      #75 percentile
      
      MR75[1]<-(((Sqr[1]/(gamma(1-thet)))^(1/1))/
                  ((log(4/3))^thet))
      
      # draw random realization from the two frechet distibutions
      set.seed(2)
      Br=rfrechet(n, loc=0, scale=Sqr[i]/(gamma(1-thet)),shape=1/thet)
      set.seed(2)
      Bc=rfrechet(n, loc=0, scale=Sqc[i]/(gamma(1-thet)),shape=1/thet)
      pr<-mcur#*(1-betaa[i])/(1-betaa[i]*mcur)
      for (j in 1:n) {if (B1[j]>pr) {B[j]=Bc[j]} else {B[j]=Br[j]}}
      MT[i]=median(B)
      MT75[i]=quantile(B,0.75)
      
      #Average
      #MT[1]<-Sq[1]
      #MR[1]<-Sqr[1]
      
      sum=sum+((MT[1]-EqM[1])/EqM[1])^2
      sum=sum+((MT75[1]-EqQ75[1])/EqQ75[1])^2
      #sum=sum+((MR[1]-EqrM[1])/EqrM[1])^2
      #sum=sum+((MR75[1]-EqrQ75[1])/EqrQ75[1])^2
    }}
  
  # if(mt(betaa,pe,1,maxx,thet)<0.5){sum=10000}      
  
  return(sum)}


ff<-function(x) -mini(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9])



# call the genetic alogorithm for the estimation
GA <- ga(type = "real-valued", 
         fitness = ff, 
         lower = c(0.12,0.12,0.12,0.12,1.7,1.0,0.15,3.5-0.35,Eqr[1]-0.78),
         upper = c(0.25,0.25,0.25,0.25,2.7,2.5,0.4,3.5+0.4 ,Eqr[1]-0.1),
         popSize = 500, maxiter = 180, run = 50,
         elitism = max(1, round(200*0.15)),
         optim=FALSE,
         seed=1,
         optimArgs = list(method = "L-BFGS-B",poptim = 0.2,pressel = 0.5,
                          control = list(fnscale = -1, maxit = 20)))

# visualize the solution 
summary(GA)
plot(GA)

beta<-c(0,GA@solution[1,1],GA@solution[1,2],GA@solution[1,3],GA@solution[1,4])
p<-exp(GA@solution[1,5])
nu<-GA@solution[1,6]
theta<-GA@solution[1,7]
Sqr[1]<-GA@solution[1,9]
Sq[1]<-fq(theta,p,GA@solution[1,8],Sqr[1])
Sqc[1]<-GA@solution[1,8]



###PART 3: Simulation + graphs########################################################

betaT<-numeric(length = )





Smbeta=numeric(length=5)
Sz=numeric(length=5)
Sm=numeric(length=5)
Ez=numeric(length=5)
Em=numeric(length=5)
SqrM<-numeric(length=5)
SNCM<-numeric(length=5)
SqM<-numeric(length=5)
SqrQ75<-numeric(length=5)
SNCQ75<-numeric(length=5)
SqQ75<-numeric(length=5)
max<-1
betat<-beta



for(i in 1:5){
  
  
  Sm[i]=mt(betat[i],p,i,max,theta,Sqc[1],Sqr[1])
  Sz[i]=zt(betat[i],p,i,max,theta,Sqc[1],Sqr[1])
  Eqc[1]=Sqc[1]
  Sqc[1]=Eqc[1]
  Ez[i]=(Eqr[i]/Eqc[i])^(1/theta)
  Em[i]=Ez[i]/(p+Ez[i])
  Smbeta[i]=betat[i]*mt(betat[i],p,i,max,theta,Sqc[1],Sqr[1])
  
  
  
  
  if(i>=2 ){
    
    Sqc[i]=qcs(betat[i],p,i,max,Sqc[i-1],nu,theta,Sqc[1],Sqr[1])
    Sqr[i]=qrs(betat[i],p,i,max,Sqr[i-1],nu,theta,Sqc[1],Sqr[1])
    Sq[i]=qs(betat[i],p,i,max,Sqr[i-1],Sqc[i-1],nu,theta,Sqc[1],Sqr[1])
    SNC[i]=nc(betat[i],p,i,max,Sqr[i-1],Sqc[i-1],nu,theta,Sqc[1],Sqr[1])
    
    #Below ONLY IF WE CONSIDER SHORT ESTIMATION PERIOD#
    if(maxper<5 & i==5){
      Smbeta[5]=betat[i]*mt(betat[i],p,i,max,theta,Sqc[4],Sqr[4])

      Sm[i]=mt(0,p,2,max,theta,Sqc[4],Sqr[4])
      Sqc[i]=qcs(0,p,2,max,Sqc[4],nu,theta,Sqc[4],Sqr[4])
      Sqr[i]=qrs(0,p,2,max,Sqr[4],nu,theta,Sqc[4],Sqr[4])
      Sq[i]=qs(0,p,2,max,Sqr[4],Sqc[4],nu,theta,Sqc[4],Sqr[4])
      SNC[i]=nc(0,p,2,max,Sqr[4],Sqc[4],nu,theta,Sqc[4],Sqr[4])
      
      
      
    }
    
    #median
    SqM[i]=(((Sq[i]/(gamma(1-theta)))^(1/1))/
              ((log(2))^theta))
    SqrM[i]=(((Sqr[i]/(gamma(1-theta)))^(1/1))/
               ((log(2))^theta))
    SqQ75[i]=(((Sq[i]/(gamma(1-theta)))^(1/1))/
                ((log(4/3))^theta))
    SqrQ75[i]=(((Sqr[i]/(gamma(1-theta)))^(1/1))/
                 ((log(4/3))^theta))
    
    # draw random realization from the two frechet distibutions
    Br=rfrechet(n, loc=0, 
                scale=(qrs(betat[i],p,i,max,Sqr[i-1],nu,theta,Sqc[1],Sqr[1])/(gamma(1-theta))),
                shape=1/theta)
    
    Bc=rfrechet(n, loc=0, 
                scale=qcs(betat[i],p,i,max,Sqc[i-1],nu,theta,Sqc[1],Sqr[1])/(gamma(1-theta)),
                shape=1/theta)
    
    pr<-Sm[i]*(1-betat[i])/(1-betat[i]*Sm[i])
    for (j in 1:n) {if (B1[j]>pr) {B[j]=Bc[j]} else {B[j]=Br[j]}}
    SNCM[i]=median(B)
    SNCQ75[i]=quantile(B,0.75)
    
    pr<-Sm[i]#*(1-betat)/(1-betat*Sm[i])
    for (j in 1:n) {if (B1[j]>pr) {B[j]=Bc[j]} else {B[j]=Br[j]}}
    SqM[i]=median(B)
    SqQ75[i]=quantile(B,0.75)
    
    
    
  }else{
    
    
    
    
    
    #SqQ75[1]=(((Sq[1]/(gamma(1-theta)))^(1/1))/
    #           ((log(4/3))^theta))
    SqrQ75[1]=(((Sqr[1]/(gamma(1-theta)))^(1/1))/
                 ((log(4/3))^theta))
    #SqM[1]=(((Sq[1]/(gamma(1-theta)))^(1/1))/
    #         ((log(2))^theta))
    SqrM[1]=(((Sqr[1]/(gamma(1-theta)))^(1/1))/
               ((log(2))^theta))
    SNC[1]=Sqr[1]*
      (((1)*mt(betat[i],p,i,max,theta,Sqc[1],Sqr[1]))/(1))+
      ((1-mt(betat[i],p,i,max,theta,Sqc[1],Sqr[1]))/(1))*
      Sqc[1]
    
    # draw random realization from the two frechet distibutions
    Br=rfrechet(n, loc=0, 
                scale=Sqr[1]/(gamma(1-theta)),
                shape=1/theta)
    
    Bc=rfrechet(n, loc=0, 
                scale=Sqc[1]/gamma(1-theta),
                shape=1/theta)
    
    pr<-Sm[i]
    for (j in 1:n) {if (B1[j]>pr) {B[j]=Bc[j]} else {B[j]=Br[j]}}
    
    SNCM[1]=median(B)
    SNCQ75[1]=quantile(B,0.75)
    SqM[1]=median(B)
    SqQ75[1]=quantile(B,0.75)
    
  }
  
}



# graph to visualize the reults
setwd("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\paper")

  
  ###############################
  #Graph Q50+Q75
  ################################
  time=seq(1,5,1) 
  q.df<- data.frame(time, EqM,SqM,EqQ75,SqQ75,CIq[1,],CIq[2,],CIq75[1,],CIq75[2,])
  
  tikz(file = "q.tex", width = 5, height = 5)
  
  #Simple plot of the dummy data using LaTeX elements
  plot <- ggplot(q.df, aes(x = time)) + 
    geom_ribbon(aes(ymin = CIq.1..., ymax = CIq.2...), fill = "red",alpha=0.15)+
    geom_line(aes(y = EqM), color = "red",size=2) + 
    geom_line(aes(y = SqM), color="red", linetype="dashed",size=2) +
    geom_ribbon(aes(ymin = CIq75.1..., ymax = CIq75.2...), fill = "blue",alpha=0.15)+
    geom_line(aes(y = EqQ75), color = "blue",size=2) + 
    geom_line(aes(y = SqQ75), color="blue", linetype="dashed",size=2) +
    #Space does not appear after Latex 
    
    labs( x = "Time", y = "Knowledge Quality") +
    
    ylim(2.9, 9)+
    
    theme_classic(
      base_family = "",
      base_line_size = 1/22,
      base_rect_size = 1/22)+
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      #axis.title.y = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0)),
      axis.text.x  = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      axis.text.y  = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0))
    )
  
  #Necessary to close or the tikxDevice .tex file will not be written
  print(plot)
  dev.off()
  
  ###############################
  #Graph Qr50+Qr75
  ################################
  qr.df<- data.frame(time, EqrM,SqrM,EqrQ75,SqrQ75,CIqr[1,],CIqr[2,],CIqr75[1,],CIqr75[2,])
  
  tikz(file = "qr.tex", width = 5, height = 5)
  
  #Simple plot of the dummy data using LaTeX elements
  plot <- ggplot(qr.df, aes(x = time)) + 
    geom_ribbon(aes(ymin = CIqr.1..., ymax = CIqr.2...), fill = "red",alpha=0.15)+
    geom_line(aes(y = EqrM), color = "red",size=2) + 
    geom_line(aes(y = SqrM), color="red", linetype="dashed",size=2) +
    geom_ribbon(aes(ymin = CIqr75.1..., ymax = CIqr75.2...), fill = "blue",alpha=0.15)+
    geom_line(aes(y = EqrQ75), color = "blue",size=2) + 
    geom_line(aes(y = SqrQ75), color="blue", linetype="dashed",size=2) +
    geom_rect(aes(xmin=1,
                  xmax =2,
                  ymin = -Inf,
                  ymax = Inf), fill = 'white') +
    #Space does not appear after Latex 
    
    labs( x = "Time", y = "Knowledge Quality") +
    
    ylim(2.9, 9)+
    
    theme_classic(
      base_family = "",
      base_line_size = 1/22,
      base_rect_size = 1/22)+
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      #axis.title.y = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0)),
      axis.text.x  = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      axis.text.y  = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0))
    )
  
  #Necessary to close or the tikxDevice .tex file will not be written
  print(plot)
  dev.off()
  
  
  ###############################
  #Graph Qr50+Qr75
  ################################
  qnc.df<- data.frame(time, ENCM,SNCM,ENCQ75,SNCQ75,CINC[1,],CINC[2,],CINC75[1,],CINC75[2,])
  
  tikz(file = "qnc.tex", width = 5, height = 5)
  
  #Simple plot of the dummy data using LaTeX elements
  plot <- ggplot(qnc.df, aes(x = time)) + 
    geom_ribbon(aes(ymin = CINC.1..., ymax = CINC.2...), fill = "red",alpha=0.15)+
    geom_line(aes(y = ENCM), color = "red",size=2) + 
    geom_line(aes(y = SNCM), color="red", linetype="dashed",size=2) +
    geom_ribbon(aes(ymin = CINC75.1..., ymax = CINC75.2...), fill = "blue",alpha=0.15)+
    geom_line(aes(y = ENCQ75), color = "blue",size=2) + 
    geom_line(aes(y = SNCQ75), color="blue", linetype="dashed",size=2) +
    geom_rect(aes(xmin=1,
                  xmax =2,
                  ymin = -Inf,
                  ymax = Inf), fill = 'white') +
    #Space does not appear after Latex 
    
    labs( x = "Time", y = "Knowledge Quality") +
    
    ylim(2.9, 9)+
    
    theme_classic(
      base_family = "",
      base_line_size = 1/22,
      base_rect_size = 1/22)+
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      #axis.title.y = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0)),
      axis.text.x  = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      axis.text.y  = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0))
    )
  
  #Necessary to close or the tikxDevice .tex file will not be written
  print(plot)
  dev.off()
  
  ###############################
  #Graph Qr50+Qr75
  ################################
  Embeta<-Embeta*100
  Smbeta<-Smbeta*100*impbeta
  CImbeta<-CImbeta*100
  b.df<- data.frame(time, Embeta,Smbeta,CImbeta[1,],CImbeta[2,])
  
  tikz(file = "b.tex", width = 5, height = 5)
  
  #Simple plot of the dummy data using LaTeX elements
  plot <- ggplot(b.df, aes(x = time)) + 
    geom_ribbon(aes(ymin = CImbeta.1..., ymax = CImbeta.2...), fill = "darkgrey",alpha=0.2)+
    geom_line(aes(y = Embeta), color = "darkgrey",size=2) + 
    geom_line(aes(y = Smbeta), color="darkgrey", linetype="dashed",size=2) +
    geom_rect(aes(xmin=1,
                  xmax =2,
                  ymin = -Inf,
                  ymax = Inf), fill = 'white') +
    #Space does not appear after Latex 
    
    labs( x = "Time", y = "\\% censored authors") +
    
    
    theme_classic(
      base_family = "",
      base_line_size = 1/22,
      base_rect_size = 1/22)+
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      #axis.title.y = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0)),
      axis.text.x  = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      axis.text.y  = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0))
    )
  
  #Necessary to close or the tikxDevice .tex file will not be written
  print(plot)
  dev.off()
  




###PART 4: counterfactual experiments##################

#No censorship
SzC=numeric(length=5)
SmC=numeric(length=5)
SmbetaC=numeric(length=5)
SqcC=numeric(length=5)
SqrC=numeric(length=5)
SqC=numeric(length=5)
for(i in 1:5){
  
  SmC[i]=mt(0,p,i,max,theta,Sqc[1],Sqr[1])
  SzC[i]=zt(0,p,i,max,theta,Sqc[1],Sqr[1])
  SmbetaC[i]=mbetat(0,p,i,max,theta,Sqc[1],Sqr[1])
  
  
  if(i>=2){
    
    SqcC[i]=qcs(0,p,i,max,SqcC[i-1],nu,theta,Sqc[1],Sqr[1])
    SqrC[i]=qrs(0,p,i,max,SqrC[i-1],nu,theta,Sqc[1],Sqr[1])
    SqC[i]=qs(0,p,i,max,SqrC[i-1],SqcC[i-1],nu,theta,Sqc[1],Sqr[1])
  }else{
    
    SqcC[1]=Sqc[1]
    SqC[1]=Sq[1]
    SqrC[1]=Sqr[1]}}



  ###############################
  #Graph bm
  ################################
  SmbetaC<-SmbetaC*100
  Smbeta<-Smbeta*0
  bmc.df<- data.frame(time, Smbeta,SmbetaC)
  
  tikz(file = "bmc.tex", width = 5, height = 5)
  
  #Simple plot of the dummy data using LaTeX elements
  plot <- ggplot(b.df, aes(x = time)) + 
    geom_line(aes(y = Smbeta), color = "darkgrey",size=2) + 
    geom_line(aes(y = SmbetaC), color="darkgrey", linetype="dashed",size=2) +
    geom_rect(aes(xmin=1,
                  xmax =2,
                  ymin = -Inf,
                  ymax = Inf), fill = 'white') +
    #Space does not appear after Latex 
    
    labs( x = "Time") +
    
    
    theme_classic(
      base_family = "",
      base_line_size = 1/22,
      base_rect_size = 1/22)+
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      #axis.title.y = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0)),
      axis.text.x  = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      axis.text.y  = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0))
    )
  
  #Necessary to close or the tikxDevice .tex file will not be written
  print(plot)
  dev.off()
  
  
  ###############################
  #Graph m
  ################################
  SmC<-SmC*100
  Sm<-Sm*100
  bc.df<- data.frame(time, Sm,SmC)
  
  tikz(file = "bc.tex", width = 5, height = 5)
  
  #Simple plot of the dummy data using LaTeX elements
  plot <- ggplot(bc.df, aes(x = time)) + 
    geom_line(aes(y = Sm), color = "darkgrey",size=2) + 
    geom_line(aes(y = SmC), color="darkgrey", linetype="dashed",size=2) +
    #Space does not appear after Latex 
    
    labs( x = "Time") +
    
    
    theme_classic(
      base_family = "",
      base_line_size = 1/22,
      base_rect_size = 1/22)+
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      #axis.title.y = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0)),
      axis.text.x  = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      axis.text.y  = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0))
    )
  
  #Necessary to close or the tikxDevice .tex file will not be written
  print(plot)
  dev.off()
  
  ###############################
  #Graph sq
  ################################
  sq.df<- data.frame(time, Sq,SqC,Sqr,SqrC,Sqc,SqcC)
  
  tikz(file = "sq.tex", width = 5, height = 5)
  
  #Simple plot of the dummy data using LaTeX elements
  plot <- ggplot(b.df, aes(x = time)) + 
    geom_line(aes(y = Sqr), color ="red",size=2) + 
    geom_line(aes(y = SqrC), color="red", linetype="dashed",size=2) +
    geom_line(aes(y = Sq), color ="blue",size=2) + 
    geom_line(aes(y = SqC), color="blue", linetype="dashed",size=2) +
    geom_line(aes(y = Sqc), color ="orange",size=2) + 
    geom_line(aes(y = SqcC), color="orange", linetype="dashed",size=2) +
    #Space does not appear after Latex 
    
    labs( x = "Time") +
    
    
    theme_classic(
      base_family = "",
      base_line_size = 1/22,
      base_rect_size = 1/22)+
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      #axis.title.y = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0)),
      axis.text.x  = element_text(size=18,margin = margin(t = 10, r =  0, b = 0, l = 0)),
      axis.text.y  = element_text(size=18,margin = margin(t = 0,  r = 10, b = 0, l = 0))
    )
  
  #Necessary to close or the tikxDevice .tex file will not be written
  print(plot)
  dev.off()

 
###PART 5: A final computation of the Impact on knowledge#########


#overall drop

print("The overall % drop is:")
print((Sq[5]-SqC[5])/Sq[5])

####
#decomposition of the drop
###
delta=Sq[5]-SqC[5]
deltas=Sqr[5]-Sqc[5]
deltasC=SqrC[5]-SqcC[5]
ch_within<-((SmC[5]/100)*(Sqr[5]-SqrC[5])+(1-SmC[5]/100)*(Sqc[5]-SqcC[5]))#change within
ch_int<-((Sm[5]-SmC[5])/100*(deltas-deltasC))#interaction
ch_across<-((Sm[5]-SmC[5])/100*(SqrC[5])+(SmC[5]-Sm[5])/100*(SqcC[5]))#change in prevalence revol.


cat("The overall absolute drop is: ",delta)
cat("The within absolute drop is: ",ch_within,"absolute is",-ch_within/delta)
cat("The across absolute drop is: ",ch_across,"absolute is",-ch_across/delta)
cat("The int    absolute drop is: ",ch_int,"absolute is",-ch_int/delta)





#How much of the drop?
(SqC[5]-Sq[5])/(Sq[2]-Sq[5])


print(c((Sq[5]-SqC[5])/SqC[5],
        (Sm[5]-SmC[5])/SmC[5],
        beta,
        Sm[5]))


#LINE REGARDING WITH THE RESULTS
betaa<-(beta[2]+beta[3]+beta[4]+beta[5])/4
setwd("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\paper")
nome<-"betav"
name <- paste(nome,".Rnw", sep = "")
results<-c((Sq[5]-SqC[5])/Sq[5]*100,(Sm[5]-SmC[5])/Sm[5]*100,betaa*100,Sm[5])


sink(name)
cat(paste('&  ',round(results[1], digits=0),'\\hspace{-0.1cm}\\% & ',round(results[2], digits=0),'\\hspace{-0.1cm}\\% 
           &  ',round(results[3], digits=0),'\\hspace{-0.1cm}\\%$^*$ &' ,round(results[4], digits=0),'\\hspace{-0.1cm}\\%\\\\'))
sink()
Sweave(name)
