################################################
## Calibration of the Church Censorship Model
## Authors: Fabio Fabio and David de la Croix
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
set.seed(1)



###PART 1: GET THE DATA####
setwd("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\data_work\\Fabio")
mbau<-list.load("data.Rdata")

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

setwd("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\calibration")
###PART 2: ESTIMATION##########


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
pr1<-numeric()
p1<-numeric()
eqrt<-numeric()
eqt<-numeric()
eqrtt<-numeric()
eqtt<-numeric()
# match initial conditions


# some functions to get the simulated variables

f<-function(x,the,k,pr,eqt,eqrt) {
  
  -eqt+(((eqrt/x)^(1/the))/((eqrt/x)^(1/the)+pr))*eqrt+
    (1-(((eqrt/x)^(1/the))/((eqrt/x)^(1/the)+pr)))*x}

fm<-function(x,the,k,pr,eqt,eqrt) {
  
  -eqt+(((eqrt/x)^(1/the))/((eqrt/x)^(1/the)+pr))*eqrt+
    (1-(((eqrt/x)^(1/the))/((eqrt/x)^(1/the)+pr)))*x}

zt<-function(beta,pe,t,up,thet,eqtt,eqrtt,p1) {
  
  # get z0 first
  v=uniroot(f,c(0.00001,up),tol=0.00000001,maxiter=10000,the=thet,k=1,
            pr=p1,eqt=eqtt,eqrt=eqrtt)
  z0=(eqrtt/v$root)^(1/thet)
  if(t>2){
    z0=p1*(z0/p1)^(2)
    pe/(1-beta)*(z0*(1-beta)/pe)^(2^(t-2))}else if(t==2){p1*(z0/p1)^(2)}else{z0}}

mue<-function(t)ifelse(t==1,200,ifelse(t==2,87.8,ifelse(t==3,78.7,ifelse(t==4,82.8,ifelse(t==5,85.1,0)))))/100

mt<-function(beta,pr,t,up,thet,eqtt,eqrtt,p1) if(t<=2){zt(beta,pr,t,up,thet,eqtt,eqrtt,p1)/(p1+zt(beta,pr,t,up,thet,eqtt,eqrtt,p1))}else
{zt(beta,pr,t,up,thet,eqtt,eqrtt,p1)/(pr+zt(beta,pr,t,up,thet,eqtt,eqrtt,p1))}  

mbetat<-function(beta,pr,t,up,thet,eqtt,eqrtt,p1) beta*mt(beta,pr,t,up,thet,eqtt,eqrtt,p1)

qrs<-function(beta,pe,t,up,mu,pre,thet,eqtt,eqrtt,p1) if(t<=2){gamma(1-thet)*(mue(t)*mu*((1)*(pre/gamma(1-thet))^(1/thet)*
                                                                                    (mt(0,p1,t-1,up,thet,eqtt,eqrtt,p1))))^thet}else
                                                                                    {gamma(1-thet)*(mue(t)*mu*((1-beta)*(pre/gamma(1-thet))^(1/thet)*
                                                                                                          (mt(beta,pe,t-1,up,thet,eqtt,eqrtt,p1))))^thet} 


qcs<-function(beta,pe,t,up,mu,pre,thet,eqtt,eqrtt,p1) if(t<=2){gamma(1-thet)*((mue(t)*mu*((pre/gamma(1-thet))^(1/thet)*
                                                                                     (1-mt(0,p1,t-1,up,thet,eqtt,eqrtt,p1)))))^thet}else
                                                                                     {gamma(1-thet)*((mue(t)*mu*((pre/gamma(1-thet))^(1/thet)*
                                                                                                            (1-mt(beta,pe,t-1,up,thet,eqtt,eqrtt,p1)))))^thet}             


qs<-function(beta,pe,t,up,mu,pre1,pre2,thet,eqtt,eqrtt,p1) if(t<=2){qrs(0,p1,t,up,mu,pre1,thet,eqtt,eqrtt,p1)*mt(0,p1,t,up,thet,eqtt,eqrtt,p1)+
    (1-mt(0,p1,t,up,thet,eqtt,eqrtt,p1))*qcs(0,p1,t,up,mu,pre2,thet,eqtt,eqrtt,p1)}else
    {qrs(beta,pe,t,up,mu,pre1,thet,eqtt,eqrtt,p1)*mt(beta,pe,t,up,thet,eqtt,eqrtt,p1)+
        (1-mt(beta,pe,t,up,thet,eqtt,eqrtt,p1))*qcs(beta,pe,t,up,mu,pre2,thet,eqtt,eqrtt,p1)}

nc<-function(beta,pe,t,up,mu,pre1,pre2,thet,eqtt,eqrtt,p1) if(t<=2){qrs(0,p1,t,up,mu,pre1,thet,eqtt,eqrtt,p1)*
    (((1)*mt(0,p1,t,up,thet,eqtt,eqrtt,p1))/(1-mt(0,p1,t,up,thet,eqtt,eqrtt,p1)))+
    ((1-mt(0,p1,t,up,thet,eqtt,eqrtt,p1))/(1-mt(0,p1,t,up,thet,eqtt,eqrtt,p1)))*
    qcs(0,p1,t,up,mu,pre2,thet,eqtt,eqrtt,p1)}else
    {qrs(beta,pe,t,up,mu,pre1,thet,eqtt,eqrtt,p1)*
        (((1-beta)*mt(beta,pe,t,up,thet,eqtt,eqrtt,p1))/(1-beta*mt(beta,pe,t,up,thet,eqtt,eqrtt,p1)))+
        ((1-mt(beta,pe,t,up,thet,eqtt,eqrtt,p1))/(1-beta*mt(beta,pe,t,up,thet,eqtt,eqrtt,p1)))*
        qcs(beta,pe,t,up,mu,pre2,thet,eqtt,eqrtt,p1)}



#GUARDA A INITIAL CONDITION CIU C
# compute empirical and non empirical
n<-1000 #how many draws from distribution
B=numeric(length=n)
B1=numeric(length=n)
Bc=numeric(length=n)
Br=numeric(length=n)
B1=runif(n, min = 0, max = 1)  # shocks for deciding later from which distribution to draw

# define the function that we want to minimize taking beta
maxg<-100
max<-numeric()
maxx<-numeric()
muu<-numeric()
betaa<-numeric()
pe1<-numeric()
sum<-numeric(length=1)
MT<-numeric()
MR<-numeric()
MT75<-numeric()
MR75<-numeric()
qu<-numeric()
qur<-numeric()
upper<-seq(0.01,20,length.out=maxg)



#parameters for the simulation

mini<-function(betaa,pee,muu,thet,pee1,qu,qur){
  sum=0     
  
  pe=exp(pee)
  pe1=pe*pee1
  Sqr[1]= ((log(2))^thet)*gamma(1-thet)*EqrM[1]
  Sq[1]= ((log(2))^thet)*gamma(1-thet)*EqM[1]
  
  Sqr[1]=qur
  Sq[1]=qu
  # first check that I can find a Kc that rationalize data with
  # current parameter
  for(j in maxg:1){
    
    if(f(0.00001,thet,1,pe1,Sq[1],Sqr[1])*f(upper[j],thet,1,pe1,Sq[1],Sqr[1])<0){
      
      maxx=upper[j]
      break}
    
    maxx=400}
  
  # for(j in maxg:1){
  
  #   if(f(0.00001,thet,1,pe,Sq[1],Sqr[1])*f(upper[j],thet,1,pe,Sq[1],Sqr[1])<0){
  
  #   maxx=upper[j]
  #   break}
  
  # maxx=400}
  
  # compute tha loss function
  if(f(0.00001,thet,1,pe1,Sq[1],Sqr[1])*f(maxx,thet,1,pe1,Sq[1],Sqr[1])>0){return(10000)}
  #if(f(0.00001,thet,1,pe1,Sq[1],Sqr[1])*f(maxx,thet,1,pe,Sq[1],Sqr[1])>0){return(10000)}
  
  # get Sqc[1],Eqc[1] consistent with parameters
  v=uniroot(f,c(0.00001,maxx),tol=0.00000001,maxiter=10000,the=thet,k=1,pr=pe1,
            eqt=Sq[1],eqrt=Sqr[1])
  Eqc[1]<-v$root
  Sqc[1]=Eqc[1]
  
  
  for(i in 1:5){
    if(i>=2){
      
      Sqc[i]<-qcs(betaa,pe,i,maxx,muu,Sqc[i-1],thet,Sq[1],Sqr[1],pe1)
      Sqr[i]<-qrs(betaa,pe,i,maxx,muu,Sqr[i-1],thet,Sq[1],Sqr[1],pe1)
      Sq[i] <-qs(betaa,pe,i,maxx,muu,Sqr[i-1],Sqc[i-1],thet,Sq[1],Sqr[1],pe1)}}         
  
  
  
  
  
  
  for (i in 1:5){
    
   
      
      mcur<-mt(betaa,pe,i,maxx,thet,Sq[1],Sqr[1],pe1)        
      
      
      
      
      
      
      
      
      if(i>=2){
        
        #Moments with average
        sum=sum+((mcur*betaa-Embeta[i])/Embeta[i])^2 
      

      MR[i]<-(((qrs(betaa,pe,i,maxx,muu,Sqr[i-1],thet,Sq[1],Sqr[1],pe1)/(gamma(1-thet)))^(1/1))/
                ((log(2))^thet))
      
     
      MR75[i]<-(((qrs(betaa,pe,i,maxx,muu,Sqr[i-1],thet,Sq[1],Sqr[1],pe1)/(gamma(1-thet)))^(1/1))/
                  ((log(4/3))^thet))
      
      
      
      # draw random realization from the two frechet distibutions
      set.seed(2)
      Br=rfrechet(n, loc=0, scale=Sqr[i]/(gamma(1-thet)),shape=1/thet)
      set.seed(2)
      Bc=rfrechet(n, loc=0, scale=Sqc[i]/(gamma(1-thet)),shape=1/thet)
      pr<-mcur#*(1-betaa)/(1-betaa*mcur)
      for (j in 1:n) {if (B1[j]>pr) {B[j]=Bc[j]} else {B[j]=Br[j]}}
      MT[i]=median(B)
      MT75[i]=quantile(B,0.75)
      
      
      sum=sum+((MT[i]-EqM[i])/EqM[i])^2
      sum=sum+((MT75[i]-EqQ75[i])/EqQ75[i])^2
      #sum=sum+((MR[i]-EqrM[i])/EqrM[i])^2
      #sum=sum+((MR75[i]-EqrQ75[i])/EqrQ75[i])^2
      
      
    }else{
      
      #Median & 75 percentile
     
      MR[1]<-(((Sqr[1]/(gamma(1-thet)))^(1/1))/
                ((log(2))^thet))
      
     
      MR75[1]<-(((Sqr[1]/(gamma(1-thet)))^(1/1))/
                  ((log(4/3))^thet))
      
      # draw random realization from the two frechet distibutions
      set.seed(2)
      Br=rfrechet(n, loc=0, scale=Sqr[i]/(gamma(1-thet)),shape=1/thet)
      set.seed(2)
      Bc=rfrechet(n, loc=0, scale=Sqc[i]/(gamma(1-thet)),shape=1/thet)
      pr<-mcur#*(1-betaa)/(1-betaa*mcur)
      for (j in 1:n) {if (B1[j]>pr) {B[j]=Bc[j]} else {B[j]=Br[j]}}
      MT[i]=median(B)
      MT75[i]=quantile(B,0.75)
      
      sum=sum+((MT[1]-EqM[1])/EqM[1])^2
      sum=sum+((MT75[1]-EqQ75[1])/EqQ75[1])^2
      #sum=sum+((MR[1]-EqrM[1])/EqrM[1])^2
      #sum=sum+((MR75[1]-EqrQ75[1])/EqrQ75[1])^2
    }}
  
  # if(mt(betaa,pe,1,maxx,thet)<0.5){sum=10000}      
  
  return(sum)}


ff<-function(x) -mini(x[1],x[2],x[3],x[4],x[5],x[6],x[7])


# call the genetic alogorithm for the estimation
GA <- ga(type = "real-valued", 
         fitness = ff, 
         lower = c(0.14,2.4,2.1,0.25,0.85 ,Eq[1]+0.3    ,Eqr[1]-0.8),
         upper = c(0.20,3.6,3.1,0.35,1.05,Eq[1]+1.4,Eqr[1]-0.2),
         popSize = 500, maxiter = 200, run = 20,
         elitism = max(1, round(200*0.15)),
         pcrossover = 0.8,
         pmutation = 0.9,
         keepBest = TRUE,
         optim=TRUE,
         seed=1,
         optimArgs = list(method = "L-BFGS-B",poptim = 0.2,pressel = 0.2,
                          control = list(fnscale = -1, maxit = 20)))

# visualize the solution 
summary(GA)
plot(GA)

beta<-GA@solution[1,1]
p<-exp(GA@solution[1,2])
mu<-GA@solution[1,3]
theta<-GA@solution[1,4]
pr1<-GA@solution[1,5]*p
Sq[1]<-GA@solution[1,6]
Sqr[1]<-GA@solution[1,7]







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

#Sqr[1]= ((log(2))^theta)*gamma(1-theta)*EqrM[1]
#q[1]= ((log(2))^theta)*gamma(1-theta)*EqM[1]

#Sqr[1]=Eqr[1]-0.5
#Sq[1]=Eq[1]

for(j in maxg:1){
  
  
  
  #  if(f(0.001,theta,1,p)*f(upper[j],theta,1,p)<0 & f(0.001,theta,2,p)*f(upper[j],theta,2,p)<0 &
  # f(0.001,theta,3,p)*f(upper[j],theta,3,p)<0 & f(0.001,theta,4,p)*f(upper[j],theta,4,p)<0){
  if(f(0.00001,theta,1,pr1,Sq[1],Sqr[1])*f(upper[j],theta,1,pr1,Sq[1],Sqr[1])<0){
    max=upper[j]
    break}}


for(i in 1:5){
  
  
  Sm[i]=mt(beta,p,i,max,theta,Sq[1],Sqr[1],pr1)
  Sz[i]=zt(beta,p,i,max,theta,Sq[1],Sqr[1],pr1)
  Smbeta[i]=beta*mt(beta,p,i,max,theta,Sq[1],Sqr[1],pr1)
  
  
  
  #v=uniroot(f,c(0.00001,max),tol=0.0001,maxiter=100,the=theta,k=i,pr=p)
  #Eqc[i]=v$root
  v=uniroot(f,c(0.00001,max),tol=0.0001,maxiter=100,the=theta,k=1,pr=pr1,eqt=Sq[1],eqrt=Sqr[1])
  Eqc[1]=v$root
  Sqc[1]=Eqc[1]
  Ez[i]=(Eqr[i]/Eqc[i])^(1/theta)
  Em[i]=Ez[i]/(pr1+Ez[i])
  
  
  
  if(i>=2){
    
    Sqc[i]=qcs(beta,p,i,max,mu,Sqc[i-1],theta,Sq[1],Sqr[1],pr1)
    Sqr[i]=qrs(beta,p,i,max,mu,Sqr[i-1],theta,Sq[1],Sqr[1],pr1)
    Sq[i]=qs(beta,p,i,max,mu,Sqr[i-1],Sqc[i-1],theta,Sq[1],Sqr[1],pr1)
    SNC[i]=nc(beta,p,i,max,mu,Sqr[i-1],Sqc[i-1],theta,Sq[1],Sqr[1],pr1)
    
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
                scale=(qrs(beta,p,i,max,mu,Sqr[i-1],theta,Sq[1],Sqr[1],pr1)/(gamma(1-theta))),
                shape=1/theta)
    
    Bc=rfrechet(n, loc=0, 
                scale=qcs(beta,p,i,max,mu,Sqc[i-1],theta,Sq[1],Sqr[1],pr1)/(gamma(1-theta)),
                shape=1/theta)
    
    for (j in 1:n) {if (B1[j]>mt(beta,p,i,max,theta,Sq[1],Sqr[1],pr1)) {B[j]=Bc[j]} else {B[j]=Br[j]}}
    SNCM[i]=median(B)
    SNCQ75[i]=quantile(B,0.75)
    
    
  }else{
    
    
    
    
    
    SqQ75[1]=(((Sq[1]/(gamma(1-theta)))^(1/1))/
                ((log(4/3))^theta))
    SqrQ75[1]=(((Sqr[1]/(gamma(1-theta)))^(1/1))/
                 ((log(4/3))^theta))
    SqM[1]=(((Sq[1]/(gamma(1-theta)))^(1/1))/
              ((log(2))^theta))
    SqrM[1]=(((Sqr[1]/(gamma(1-theta)))^(1/1))/
               ((log(2))^theta))
    SNC[1]=Sqr[1]*
      (((1)*mt(beta,pr1,i,max,theta,Sq[1],Sqr[1],pr1))/(1-mt(beta,pr1,i,max,theta,Sq[1],Sqr[1],pr1)))+
      ((1-mt(beta,pr1,i,max,theta,Sq[1],Sqr[1],pr1))/(1-mt(beta,pr1,i,max,theta,Sq[1],Sqr[1],pr1)))*
      Sqc[1]
    
    # draw random realization from the two frechet distibutions
    Br=rfrechet(n, loc=0, 
                scale=Sqr[1]/(gamma(1-theta)),
                shape=1/theta)
    
    Bc=rfrechet(n, loc=0, 
                scale=Sqc[1]/gamma(1-theta),
                shape=1/theta)
    
    for (j in 1:n) {if (B1[j]>mt(beta,pr1,1,max,theta,Sq[1],Sqr[1],pr1)) {B[j]=Bc[j]} else {B[j]=Br[j]}}
    SNCM[1]=median(B)
    SNCQ75[1]=quantile(B,0.75)
    
  }
  
  
  
}




##########################################
# counterfactual experiment: no censorship
##########################################
SzC=numeric(length=5)
SmC=numeric(length=5)
SmbetaC=numeric(length=5)
SqcC=numeric(length=5)
SqrC=numeric(length=5)
SqC=numeric(length=5)
for(i in 1:5){
  
  SmC[i]=mt(0,pr1,i,max,theta,Sq[1],Sqr[1],pr1)
  SzC[i]=zt(0,pr1,i,max,theta,Sq[1],Sqr[1],pr1)
  SmbetaC[i]=mbetat(0,pr1,i,max,theta,Sq[1],Sqr[1],pr1)
  
  
  if(i>=2){
    
    SqcC[i]=qcs(0,pr1,i,max,mu,SqcC[i-1],theta,Sq[1],Sqr[1],pr1)
    SqrC[i]=qrs(0,pr1,i,max,mu,SqrC[i-1],theta,Sq[1],Sqr[1],pr1)
    SqC[i]=qs(0,pr1,i,max,mu,SqrC[i-1],SqcC[i-1],theta,Sq[1],Sqr[1],pr1)
  }else{
    
    SqcC[1]=Sqc[1]
    SqC[1]=Sq[1]
    SqrC[1]=Sqr[1]}}

##########################################
# counterfactual experiment: no self censorship
##########################################
SzCC=numeric(length=5)
SmCC=numeric(length=5)
SmbetaCC=numeric(length=5)
SqcCC=numeric(length=5)
SqrCC=numeric(length=5)
SqCC=numeric(length=5)
for(i in 1:5){
  
  SmCC[i]=mt(0,p,i,max,theta,Sq[1],Sqr[1],pr1)
  SzCC[i]=zt(0,p,i,max,theta,Sq[1],Sqr[1],pr1)
  SmbetaCC[i]=mbetat(0,p,i,max,theta,Sq[1],Sqr[1],pr1)
  
  
  if(i>=2){
    
    SqcCC[i]=qcs(0,p,i,max,mu,SqcCC[i-1],theta,Sq[1],Sqr[1],pr1)
    SqrCC[i]=qrs(0,p,i,max,mu,SqrCC[i-1],theta,Sq[1],Sqr[1],pr1)
    SqCC[i]=qs(0,p,i,max,mu,SqrCC[i-1],SqcC[i-1],theta,Sq[1],Sqr[1],pr1)
  }else{
    
    SqcCC[1]=Sqc[1]
    SqCC[1]=Sq[1]
    SqrCC[1]=Sqr[1]}}


##########################################
# counterfactual experiment: no self censorship
##########################################
SzCCC=numeric(length=5)
SmCCC=numeric(length=5)
SmbetaCCC=numeric(length=5)
SqcCCC=numeric(length=5)
SqrCCC=numeric(length=5)
SqCCC=numeric(length=5)
for(i in 1:5){
  
  SmCCC[i]=mt(beta,pr1,i,max,theta,Sq[1],Sqr[1],pr1)
  SzCCC[i]=zt(beta,pr1,i,max,theta,Sq[1],Sqr[1],pr1)
  SmbetaCCC[i]=mbetat(beta,pr1,i,max,theta,Sq[1],Sqr[1],pr1)
  
  
  if(i>=2){
    
    SqcCCC[i]=qcs(beta,pr1,i,max,mu,SqcCCC[i-1],theta,Sq[1],Sqr[1],pr1)
    SqrCCC[i]=qrs(beta,pr1,i,max,mu,SqrCCC[i-1],theta,Sq[1],Sqr[1],pr1)
    SqCCC[i]=qs(beta,pr1,i,max,mu,SqrCCC[i-1],SqcC[i-1],theta,Sq[1],Sqr[1],pr1)
  }else{
    
    SqcCCC[1]=Sqc[1]
    SqCCC[1]=Sq[1]
    SqrCCC[1]=Sqr[1]}}

#A final computation of the Impact on knowledge

#counter:no censorship, no inquistion
print(c((Sq[5]-SqC[5])/Sq[5],
        (Sm[5]-SmC[5])/Sm[5],
        beta,
        Sm[5]))

#counter:inquisition, no cenosrhip
print(c((Sq[5]-SqCC[5])/Sq[5],
        (Sm[5]-SmCC[5])/Sm[5],
        beta,
        Sm[5]))

#overall:no inquisition,censorship
print(c((Sq[5]-SqCCC[5])/Sq[5],
        (Sm[5]-SmCCC[5])/Sm[5],
        beta,
        Sm[5]))

p^(-theta)/(pr1^(-theta))

#LINE REGARDING WITH THE RESULTS
setwd("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\paper")
nome<-"selfc"
normal<-FALSE
name <- paste(nome,".Rnw", sep = "")
if(normal){results<-c((Sq[5]-SqCC[5])/Sq[5]*100,(Sm[5]-SmCC[5])/Sm[5]*100,beta*100,Sm[5])}
if(!normal){results<-c((Sq[5]-SqCC[5])/Sq[5]*100,(Sm[5]-SmCC[5])/Sm[5]*100,beta*100,Sm[5]*100)}


sink(name)
cat(paste('&  ',round(results[1], digits=0),'\\hspace{-0.1cm}\\% & ',round(results[2], digits=0),'\\hspace{-0.1cm}\\% 
           &  ',round(results[3], digits=0),'\\hspace{-0.1cm}\\% &' ,round(results[4], digits=0),'\\hspace{-0.1cm}\\%\\\\'))
sink()
Sweave(name)
