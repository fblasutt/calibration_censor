  ################################################
  ## Calibration of the Church Censorship Model
  ## Authors: Blasutto Fabio and David de la Croix
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
  
  #Parallel Computing
  require(doParallel)
  
  #Set seed
  set.seed(2)
  ###PART 0: WHICH SAMPLE/MODEL DO YOU WANT?#############################################
  cmom<-TRUE    #want to redo-data part? 
  serrors<-FALSE  #want to compute se for parameters?
  nboots<-500    #Bootstrapping nummber
  
  
  impbeta <<-1    #imperfectSqzr application of censorship-max (5./6.)
  maxper<<-5      #max number of periods to consider
  
  ita<-FALSE      #only italian scholars
  itan<-FALSE     #only northern italian scholars
  itas<-FALSE     #only southern italian scholars
  notonlyby<-FALSE#not only publications by person 
  weak<-FALSE     #do not consider weak links
  unio<-FALSE     #Only univresity, no academies
  identif<-FALSE  #this serves for understanding how the model is identified
  fivpmod<-FALSE   #5 periods model?
  robc<-FALSE
  wiki<-FALSE
  longe<-FALSE
  
  
  if(impbeta>1){imperfect<-TRUE}else{imperfect<-FALSE}
  if(maxper<5){maxt<-TRUE}else{maxt<-FALSE}
  if(ita+itan+itas+notonlyby+maxt+imperfect+weak+robc+unio+wiki+longe+fivpmod>0){normal<-FALSE}else{normal<-TRUE}
  
  #Get name to print results
  names<-c("ita","itan","itas","notonlyby","maxt","imperfect","weak","normal","robc","unio","wiki","longe")
  valll<-c(ita,itan,itas,notonlyby,maxt,imperfect,weak,normal,robc,unio,wiki,longe)
  namet<-which.max(valll)
  nome<-names[namet]
  
  #Part below if for timing
  if(fivpmod){leng=5;pcens=2}else{leng=10;pcens=4;maxper=10}
  
  
  ###PART 1: GET THE DATA#############################################
  
  if (cmom){
    
    if(fivpmod)
      
    {
      #Moments for 5 periods model
      source("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\calibration_censor\\moments_5periods.R") 
      
    }
    
    else
      
      {
        #Moments for 10 periods model
        source("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\calibration_censor\\moments_10periods.R") 
      
      }
    
    
  }else{
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
    
  }
  
  
  ###PART 2: ESTIMATION#################################################
    
    

  t<-numeric()
  v<-numeric()
  pe<-numeric()
  z0<-numeric()
  up<-numeric()
  Sqc<-numeric(length=leng)
  Sqr<-numeric(length=leng)
  SNC<-numeric(length=leng)
  SNC75<-numeric(length=leng)
  Sq<-numeric(length=leng)
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
  #Embeta=c(0.0984931+(0.0984931-0.0917645),0.0984931,(0.0984931+0.0917645)/2,0.0917645,(0.0917645+0.08468171)/2,0.0846817,
   #        (0.0846817+0.07079675)/2,0.07079675,(0.07079675+0.046189)/2,0.046189)
  
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
    if(t>pcens){
      z0=pe/(1)*(z0*(1)/pe)^(2^(pcens-1))#pe/(1)*(z0*(1)/pe)^(2)
      min(pe/(1-beta)*(z0*(1-beta)/pe)^(2^(t-pcens)),1000000000000000000000)}else if(t<=pcens & t>1){pe/(1)*(z0*(1)/pe)^(2^(t-1))}else{z0}}
  
    
    
  
  mt<-function(beta,pr,t,up,thet,eqct,eqrtt) zt(beta,pr,t,up,thet,eqct,eqrtt)/(pr+zt(beta,pr,t,up,thet,eqct,eqrtt))
  
  mbetat<-function(beta,pr,t,up,thet,eqct,eqrtt) beta*impbeta*mt(beta,pr,t,up,thet,eqct,eqrtt)
  
  #Below Macro shocks
  if(leng==5)
     {
        if (longe) {mu<-function(t) ifelse(t==1,68.26-18,ifelse(t==2,64.03-18,ifelse(t==3,65.17-18,ifelse(t==4,64.83-18,ifelse(t==5,69.86-18,0)))))/(68.18-18)}
        if (!longe){mu<-function(t)ifelse(t==1,200,ifelse(t==2,87.8,ifelse(t==3,78.7,ifelse(t==4,82.8,ifelse(t==5,85.1,0)))))/100}
  }else{
  
        mu<-function(t) ifelse(t==1,100,ifelse(t==2,101.1,ifelse(t==3,88.4,ifelse(t==4,88.2,ifelse(t==5,81.5,ifelse(t==6,76.8,ifelse(t==7,83.1,ifelse(t==8,83.5,ifelse(t==9,81.5,ifelse(t==10,92.2,0))))))))))/100}
    
  

  
  
  
  qrs<-function(beta,pe,t,up,pre,nu,thet,eqct,eqrtt) if(t<=pcens){gamma(1-thet)*((1+nu)*mu(t)*((1)*(pre/gamma(1-thet))^(1/thet)*
                                                                                   (mt(0,pe,t-1,up,thet,eqct,eqrtt))))^thet}else
                                                                                   {gamma(1-thet)*((1+nu)*mu(t)*((1-beta)*(pre/gamma(1-thet))^(1/thet)*
                                                                                                         (mt(beta,pe,t-1,up,thet,eqct,eqrtt))))^thet} 
  
  
  qcs<-function(beta,pe,t,up,pre,nu,thet,eqct,eqrtt) if(t<=pcens){gamma(1-thet)*(((1+nu)*mu(t)*((pre/gamma(1-thet))^(1/thet)*
                                                                                    (1-mt(0,pe,t-1,up,thet,eqct,eqrtt)))))^thet}else
                                                                                    {gamma(1-thet)*(((1+nu)*mu(t)*((pre/gamma(1-thet))^(1/thet)*
                                                                                                           (1-mt(beta,pe,t-1,up,thet,eqct,eqrtt)))))^thet}             
  
  
  qs<-function(beta,pe,t,up,pre1,pre2,nu,thet,eqct,eqrtt) if(t<=pcens){qrs(0,pe,t,up,pre1,nu,thet,eqct,eqrtt)*mt(0,pe,t,up,thet,eqct,eqrtt)+
      (1-mt(0,pe,t,up,thet,eqct,eqrtt))*qcs(0,pe,t,up,pre2,nu,thet,eqct,eqrtt)}else
      {qrs(beta,pe,t,up,pre1,nu,thet,eqct,eqrtt)*mt(beta,pe,t,up,thet,eqct,eqrtt)+
          (1-mt(beta,pe,t,up,thet,eqct,eqrtt))*qcs(beta,pe,t,up,pre2,nu,thet,eqct,eqrtt)}
  
  nc<-function(beta,pe,t,up,pre1,pre2,nu,thet,eqct,eqrtt) if(t<=pcens){qrs(0,pe,t,up,pre1,nu,thet,eqct,eqrtt)*
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
  
  mini<-function(betaa,pee,nuu,thet,qu,qur){
    
    sum=0     
    pe<-exp(pee)
    Sqr[1]=qur
    Sqc[1]=qu
    Sq[1]=fq(thet,pe,Sqc[1],Sqr[1])
    
    
    
    for(i in 1:maxper){
      if(i>=2){
        
        Sqc[i]<-qcs(betaa,pe,i,maxx,Sqc[i-1],nuu,thet,Sqc[1],Sqr[1])
        Sqr[i]<-qrs(betaa,pe,i,maxx,Sqr[i-1],nuu,thet,Sqc[1],Sqr[1])
        #Sq[i] <-qs(betaa,pe,i,maxx,Sqr[i-1],Sqc[i-1],nuu,thet,Sqc[1],Sqr[1])          
        
      }}
    
    
    
    
    
    for (i in 1:maxper){
      mcur<-mt(betaa,pe,i,maxx,thet,Sqc[1],Sqr[1])        
      
      
      
      
      if(i>=2){
        
        #Share censored is target only from period before censorship pcens
        if(i>=pcens){sum=sum+((mcur*betaa*impbeta-Embeta[i])/Embeta[i])^2}
        if(i==4 & mcur>=0.507){sum=sum+10000}
        #Moments with median
        
        MR[i]<-(((qrs(betaa,pe,i,maxx,Sqr[i-1],nuu,thet,Sqc[1],Sqr[1])/(gamma(1-thet)))^(1/1))/
                  ((log(2))^thet))
        
        #Model with 75p
        
        MR75[i]<-(((qrs(betaa,pe,i,maxx,Sqr[i-1],nuu,thet,Sqc[1],Sqr[1])/(gamma(1-thet)))^(1/1))/
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
        
        #Moments with average
        #MT[i]<-qs(betaa,pe,i,maxx,nuu,Sqr[i-1],Sqc[i-1],thet,Sq[1],Sqr[1])
        #MR[i]<-qrs(betaa,pe,i,maxx,nuu,Sqr[i-1],thet,Sq[1],Sqr[1])
        
        
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
        pr<-mcur#*(1-betaa)/(1-betaa*mcur)
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
  
  
  ff<-function(x) -mini(x[1],x[2],x[3],x[4],x[5],x[6])
  
  
  
  # call the genetic alogorithm for the estimation
  GA <- ga(type = "real-valued", 
           fitness = ff, 
           #lower = c(0.12,1.7,0.0,0.05,3.5-0.35,Eqr[1]-0.68),
           #upper = c(0.25,2.5,3  ,0.45,3.5+0.4 ,Eqr[1]-0.2),
           lower = c(0.05,1.5,0.0,0.15,3.5-0.65,Eqr[1]-1.78),
           upper = c(0.20,2.8,2.3  ,0.45,3.5+0.6 ,Eqr[1]-0.8),
           popSize = 1000, maxiter = 200, run = 50,
           #elitism = max(1, round(20*0.15)),
           optim=FALSE,
           seed=1,
           parallel=TRUE,
           optimArgs = list(method = "L-BFGS-B",poptim = 0.1,pressel = 0.1,
                            control = list(fnscale = -1, maxit = 20)))
  
  # visualize the solution 
  summary(GA)
  plot(GA)
  
  beta<-GA@solution[1,1]
  p<-exp(GA@solution[1,2])
  nu<-GA@solution[1,3]
  theta<-GA@solution[1,4]
  Sqr[1]<-GA@solution[1,6]
  Sq[1]<-fq(theta,p,GA@solution[1,5],Sqr[1])
  Sqc[1]<-GA@solution[1,5]
  
  
  
  ###PART 3: Simulation + graphs########################################################
  
  betaT<-numeric(length = )
  
  
  
  
  
  Smbeta=numeric(length=leng)
  Sz=numeric(length=leng)
  Sm=numeric(length=leng)
  Ez=numeric(length=leng)
  Em=numeric(length=leng)
  SqrM<-numeric(length=leng)
  SNCM<-numeric(length=leng)
  SqM<-numeric(length=leng)
  SqrQ75<-numeric(length=leng)
  SNCQ75<-numeric(length=leng)
  SqQ75<-numeric(length=leng)
  max<-1
  betat<-beta
  
  
  
  for(i in 1:leng){
    
    
    Sm[i]=mt(betat,p,i,max,theta,Sqc[1],Sqr[1])
    Sz[i]=zt(betat,p,i,max,theta,Sqc[1],Sqr[1])
    Eqc[1]=Sqc[1]
    Sqc[1]=Eqc[1]
    Ez[i]=(Eqr[i]/Eqc[i])^(1/theta)
    Em[i]=Ez[i]/(p+Ez[i])
    Smbeta[i]=betat*mt(betat,p,i,max,theta,Sqc[1],Sqr[1])
    
    
    
    
    if(i>=2 ){
      
      Sqc[i]=qcs(betat,p,i,max,Sqc[i-1],nu,theta,Sqc[1],Sqr[1])
      Sqr[i]=qrs(betat,p,i,max,Sqr[i-1],nu,theta,Sqc[1],Sqr[1])
      Sq[i]=qs(betat,p,i,max,Sqr[i-1],Sqc[i-1],nu,theta,Sqc[1],Sqr[1])
      SNC[i]=nc(betat,p,i,max,Sqr[i-1],Sqc[i-1],nu,theta,Sqc[1],Sqr[1])
      
      #Below ONLY IF WE CONSIDER SHORT ESTIMATION PERIOD#
      if(maxper<leng & i==leng){
        Smbeta[5]=betat*mt(betat,p,i,max,theta,Sqc[4],Sqr[4])
        betat<-0
        Sm[i]=mt(betat,p,2,max,theta,Sqc[4],Sqr[4])
        Sqc[i]=qcs(betat,p,2,max,Sqc[4],nu,theta,Sqc[4],Sqr[4])
        Sqr[i]=qrs(betat,p,2,max,Sqr[4],nu,theta,Sqc[4],Sqr[4])
        Sq[i]=qs(betat,p,2,max,Sqr[4],Sqc[4],nu,theta,Sqc[4],Sqr[4])
        SNC[i]=nc(betat,p,2,max,Sqr[4],Sqc[4],nu,theta,Sqc[4],Sqr[4])
        
        
        
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
                  scale=(qrs(betat,p,i,max,Sqr[i-1],nu,theta,Sqc[1],Sqr[1])/(gamma(1-theta))),
                  shape=1/theta)
      
      Bc=rfrechet(n, loc=0, 
                  scale=qcs(betat,p,i,max,Sqc[i-1],nu,theta,Sqc[1],Sqr[1])/(gamma(1-theta)),
                  shape=1/theta)
      
      pr<-Sm[i]*(1-betat)/(1-betat*Sm[i])
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
        (((1)*mt(betat,p,i,max,theta,Sqc[1],Sqr[1]))/(1))+
        ((1-mt(betat,p,i,max,theta,Sqc[1],Sqr[1]))/(1))*
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
  if(normal){
    
    ###############################
    #Graph Q50+Q75
    ################################
    time=seq(1,leng,1) 
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
    
  }
  
  
  
  ###PART 4: counterfactual experiments##################
  
  #No censorship
  SzC=numeric(length=leng)
  SmC=numeric(length=leng)
  SmbetaC=numeric(length=leng)
  SqcC=numeric(length=leng)
  SqrC=numeric(length=leng)
  SqC=numeric(length=leng)
  for(i in 1:leng){
    
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
  
  
  #No Economic decline
  SzCd=numeric(length=leng)
  SmCd=numeric(length=leng)
  SmbetaCd=numeric(length=leng)
  SqcCd=numeric(length=leng)
  SqrCd=numeric(length=leng)
  SqCd=numeric(length=leng)
  
  #Update: NO decline of Italy
  mu<-function(t) 1
  for(i in 1:leng){
    
    SmCd[i]=mt(betat,p,i,max,theta,Sqc[1],Sqr[1])
    SzCd[i]=zt(betat,p,i,max,theta,Sqc[1],Sqr[1])
    SmbetaCd[i]=mbetat(betat,p,i,max,theta,Sqc[1],Sqr[1])
    
    
    if(i>=2){
      
      SqcCd[i]=qcs(betat,p,i,max,SqcCd[i-1],nu,theta,Sqc[1],Sqr[1])
      SqrCd[i]=qrs(betat,p,i,max,SqrCd[i-1],nu,theta,Sqc[1],Sqr[1])
      SqCd[i]=qs(betat,p,i,max,SqrCd[i-1],SqcC[i-1],nu,theta,Sqc[1],Sqr[1])
    }else{
      
      SqcCd[1]=Sqc[1]
      SqCd[1]=Sq[1]
      SqrCd[1]=Sqr[1]}}
  
  #No censorship, no macro shocks
  SzCdp=numeric(length=leng)
  SmCdp=numeric(length=leng)
  SmbetaCdp=numeric(length=leng)
  SqcCdp=numeric(length=leng)
  SqrCdp=numeric(length=leng)
  SqCdp=numeric(length=leng)
  for(i in 1:leng){
    
    SmCdp[i]=mt(0,p,i,max,theta,Sqc[1],Sqr[1])
    SzCdp[i]=zt(0,p,i,max,theta,Sqc[1],Sqr[1])
    SmbetaCdp[i]=mbetat(0,p,i,max,theta,Sqc[1],Sqr[1])
    
    
    if(i>=2){
      
      SqcCdp[i]=qcs(0,p,i,max,SqcCdp[i-1],nu,theta,Sqc[1],Sqr[1])
      SqrCdp[i]=qrs(0,p,i,max,SqrCdp[i-1],nu,theta,Sqc[1],Sqr[1])
      SqCdp[i]=qs(0,p,i,max,SqrCdp[i-1],SqcC[i-1],nu,theta,Sqc[1],Sqr[1])
    }else{
      
      SqcCdp[1]=Sqc[1]
      SqCdp[1]=Sq[1]
      SqrCdp[1]=Sqr[1]}}
  
  
  if(normal){
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
    
    ######################################
    #Graph sq under the two counterfacual
    #######################################
    sq.df<- data.frame(time, Sq,SqC,SqCd)
    
    tikz(file = "sqCd.tex", width = 10, height = 5)
    
    #Simple plot of the dummy data using LaTeX elements
    plot <- ggplot(b.df, aes(x = time)) + 
      geom_line(aes(y = (SqC-Sq)/Sq), color ="blue",size=2) + 
      geom_line(aes(y = (SqCd-Sq)/Sq), color ="red",size=2) + 
      #Space does not appear after Latex 
      
      labs( x = "Time") +
      
      ylab("Gains in average quality (in %)") +
      
      
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
  }
  ###PART 5: A final computation of the Impact on knowledge#########
  
  
  #overall drop
  
  print("The overall % drop is:")
  print((Sq[leng]-SqC[leng])/Sq[leng])
  
  ####
  #decomposition of the drop
  ###
  delta=Sq[leng]-SqC[leng]
  deltas=Sqr[leng]-Sqc[leng]
  deltasC=SqrC[leng]-SqcC[leng]
  ch_within<-((SmC[leng]/100)*(Sqr[leng]-SqrC[leng])+(1-SmC[leng]/100)*(Sqc[leng]-SqcC[leng]))#change within
  ch_int<-((Sm[leng]-SmC[leng])/100*(deltas-deltasC))#interaction
  ch_across<-((Sm[leng]-SmC[leng])/100*(SqrC[leng])+(SmC[leng]-Sm[leng])/100*(SqcC[leng]))#change in prevalence revol.
  
  
  cat("The overall absolute drop is: ",delta)
  cat("The within absolute drop is: ",ch_within,"absolute is",-ch_within/delta)
  cat("The across absolute drop is: ",ch_across,"absolute is",-ch_across/delta)
  cat("The int    absolute drop is: ",ch_int,"absolute is",-ch_int/delta)
  
  
  
  
  
  #How much of the drop?
  (SqC[leng]-Sq[leng])/(Sq[pcens]-Sq[leng])
  
  
  print(c((Sq[leng]-SqC[leng])/SqC[leng],
          (Sm[leng]-SmC[leng])/SmC[leng],
          beta,
          Sm[leng]))
  
  
  ###PART 6: Get value functions#########
  betax<-beta
  source("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\calibration_censor\\vf.R") 
  
  #Do you want the robustness graph?
  if(robc){
    
    #how many versions?
    seqq<-rev(seq(from = 0.8, to = 0.9, by=0.1))
    gpsi1<-numeric(length=length(seqq))
    gpsi2<-numeric(length=length(seqq))
    gmm<-numeric(length=length(seqq))
    for(j in 1:length(seqq)){
      
      
      #Re-set imperfect censorship
      impbeta<-1/seqq[j]
      
      #Re-estimate the model
      source("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\calibration_censor\\vf1.R") 
      betax<-GA@solution[1,1]
      Sqr[1]<-GA@solution[1,6]
      Sqc[1]<-GA@solution[1,5]
      source("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\calibration_censor\\vf.R") 
      gpsi1[j]<-psi1
      gpsi2[j]<-psi2
      gmm[j]<-mt(0,exp(GA@solution[1,2]),1,1,GA@solution[1,4],Sqc[1],Sqr[1])
    }
    
    #Cool graph with results
    lim.df<- data.frame(seqq,gpsi1,gpsi2,gmm)
    
    tikz(file = "lim.tex", width = 5, height = 5)
    
    #Simple plot of the dummy data using LaTeX elements
    plot <- ggplot(lim.df, aes(x = seqq)) + 
      #geom_ribbon(aes(ymin = CIq.1..., ymax = CIq.2...), fill = "red",alpha=0.15)+
      geom_line(aes(y = gpsi1), color = "blue",size=1) + 
      geom_line(aes(y = gpsi2), color="red",size=1) +
      #geom_ribbon(aes(ymin = CIq75.1..., ymax = CIq75.2...), fill = "blue",alpha=0.15)+
      #geom_line(aes(y = gmm), color = "black",size=2) + 
      #Space does not appear after Latex 
      
      labs( x = "Time", y = "Knowledge Quality") +
      
      ylim(min(gpsi1,gpsi2), max(gpsi1,gpsi2))+
      
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
    
    
    
  }
  ###PART 7: Bootrstrapped standard errors#########
  betaB<-numeric(length=nboots)
  pB<-numeric(length=nboots)
  nuB<-numeric(length=nboots)
  thetaB<-numeric(length=nboots)
  SqB<-numeric(length=nboots)
  SqrB<-numeric(length=nboots)
  SqcB<-numeric(length=nboots)
  
  
  betaS<-numeric(length=1)
  pS<-numeric(length=1)
  nuS<-numeric(length=1)
  thetaS<-numeric(length=1)
  SqS<-numeric(length=1)
  SqrS<-numeric(length=1)
  
  
  if(serrors & normal){
    mu<-function(t)ifelse(t==1,200,ifelse(t==2,87.8,ifelse(t==3,78.7,ifelse(t==4,82.8,ifelse(t==5,85.1,0)))))/100
    for(j in 1:nboots){
      print(j)
      
      #Get bootstraped moments  
      for(i in 1:leng){
        Embeta[i]=Bmbeta[i,j]
        Eqr[i]=Bqr[i,j]
        Eq[i]=Bq[i,j]
        ENC[i]=BNC[i,j]
        EqrM[i]=BqrM[i,j]
        EqM[i]=BqM[i,j]
        ENCM[i]=BNCM[i,j]
        EqrQ75[i]=Bqr75[i,j]
        EqQ75[i]=Bq75[i,j]
        ENCQ75[i]=BNC75[i,j]
      }

      #find the parameters for this draw  
      GA <- ga(type = "real-valued", 
               fitness = ff, 
               lower = c(0.1,1.3,0.0,0.05,3.5-0.65,Eqr[1]-0.8),
               upper = c(0.3,2.9,4  ,0.45,3.5+0.6 ,Eqr[1]+0.5),
               popSize = 100, maxiter = 180, run = 50,
               elitism = max(1, round(200*0.15)),
               optim=FALSE,
               monitor=FALSE,
               keepBest = TRUE,
               seed=1,
               optimArgs = list(method = "L-BFGS-B",poptim = 0.2,pressel = 0.2,
                                control = list(fnscale = -1, maxit = 20)))
      
      betaB[j]<-GA@solution[1,1]
      pB[j]<-exp(GA@solution[1,2])
      nuB[j]<-GA@solution[1,3]
      thetaB[j]<-GA@solution[1,4]
      #SqB[j]<-GA@solution[1,5]
      #SqrB[j]<-GA@solution[1,6]
      
      SqrB[j]<-GA@solution[1,6]
      Sq[j]<-fq(thetaB[j],pB[j],GA@solution[1,5],SqrB[j])
      SqcB[j]<-GA@solution[1,5]
      
      pB[j]=pB[j]^(-thetaB[j])
      
      
      
      
    }
    
    #Compute standard errors
    betaS<-sd(betaB)
    pS<-sd(pB)
    nuS<-sd(nuB)
    thetaS<-sd(thetaB)
    SqrS<-sd(SqrB)
    SqcS<-sd(SqcB)
    
    #Standard Errors+other nice things for table
    kc<-(Sqc[1]/(gamma(1-theta)))^(1/theta)
    kr<-(Sqr[1]/(gamma(1-theta)))^(1/theta)
    krs<-sd((SqrB/(gamma(1-theta)))^(1/theta))
    
    
    
    
    #SqcB<-(SqB-nuB*SqrB)/(1-nuB)
    kcs<-sd((SqcB/(gamma(1-theta)))^(1/theta))
    
    #for the table
    dkc<-nchar(round(kc,digits=0))
    dkr<-nchar(round(kr,digits=0))
    dkcs<-nchar(round(kcs,digits=0))
    dkrs<-nchar(round(krs,digits=0))
    
    #kc<-kc/(10^(dkc-1))
    #kr<-kr/(10^(dkr-1))
    #kcs<-kcs/(10^(dkcs-1))
    #krs<-krs/(10^(dkrs-1))
    
    psi1s<-psi1#*1000
    psi2s<-psi2#*1000
    #Also for qc and qr
    print(c(Sqc[1],sd(SqcB)))
    print(c(Sqr[1],sd(SqrB)))
    
    #Build the tables for the paper
    sink("tempp.Rnw")
    
    cat(paste('\\begin{table}[htpb]
      \\caption{Identification of Parameters}
      \\centering % used for centering table
      \\begin{tabular}{@{\\extracolsep{5pt}}l c c c c} 
      \\hline\\hline%inserts double horizontal lines
      \\rule{-4pt}{2.5ex}
       Calibrated Parameters &  & Value & Standard Errors & Target \\\\ [0.05ex] % inserts table
       \\hline
      \\rule{-4pt}{2.5ex}
       Discount Factor  & $\\delta$   &0.06& -  &RBC literature\\\\[0.15ex]
       Fixed Cost of Censorship  & $\\hat{\\psi}$   &(\\hspace{-0.1cm}',round(psi1s, digits=3),'-',round(psi2s, digits=3),'\\hspace{-0.1cm})& - &Index set-up \\\\[0.15ex]
       \\hline % inserts single horizontal line
       \\rule{-4pt}{2.5ex}
       Estimated Parameters &  & Value & Standard Errors & Target  \\\\ [0.05ex] % inserts table
        %heading
      \\hline % inserts single horizontal line
      \\rule{-4pt}{2.5ex}
      Compliant knowledge in 1  & $k^C_1$   &',round(kc, digits=1),'&',round(kcs, digits=2),' &  $\\Omega(\\vartheta)$  \\\\[0.15ex]
      Rev. knowledge in 1  & $k^R_1$   &',round(kr, digits=1),'&',round(krs, digits=2),' & $\\Omega(\\vartheta)$ \\\\[0.15ex]
      Productivity of books  & $\\theta$   &',round(theta, digits=2),'& ',round(thetaS, digits=3),' & $\\Omega(\\vartheta)$\\\\[0.15ex]
      Max Censorship  & $\\overline{\\beta}$   &',round(beta, digits=2),'& ',round(betaS, digits=3),'& $\\Omega(\\vartheta)$\\\\[0.15ex]
      Knowledge Growth   & $\\nu$   &',round(nu, digits=2),'& ',round(nuS, digits=3),' & $\\Omega(\\vartheta)$\\\\[0.15ex]
      Price of rev. books   & $p$   &',round(p^(-theta), digits=2),'& ',round(pS, digits=3),' & $\\Omega(\\vartheta)$\\\\[0.15ex]
      \\hline\\hline
      \\end{tabular}
      \\label{table:param}
      \\end{table}'))
    
    sink()
    Sweave("tempp.Rnw")
    
  }
  
  
  ###PART 8: Final Tables#########
  #Embeta<-Embeta*100
  
  
  
  
  #LINE REGARDING WITH THE RESULTS
  
  name <- paste(nome,".Rnw", sep = "")
  if(normal){results<-c((Sq[leng]-SqC[leng])/Sq[leng]*100,(Sm[leng]-SmC[leng])/Sm[leng]*100,beta*100*impbeta,Sm[leng])}
  if(!normal){results<-c((Sq[leng]-SqC[leng])/Sq[leng]*100,(Sm[leng]-SmC[leng])/Sm[leng]*100,beta*100*impbeta,Sm[leng]*100)}
  
  
  sink(name)
  cat(paste('&  ',round(results[1], digits=0),'\\hspace{-0.1cm}\\% & ',round(results[2], digits=0),'\\hspace{-0.1cm}\\% 
             &  ',round(results[3], digits=0),'\\hspace{-0.1cm}\\% &' ,round(results[4], digits=0),'\\hspace{-0.1cm}\\%\\\\'))
  sink()
  Sweave(name)
  
  #for trying the identification
  if(identif){
    setwd("C:\\Users\\Fabio\\Dropbox\\Roman_Church_censorship_growth\\data_work\\Fabio")
    momentss<-list("Embeta"=Sm*beta/100,"Eqr"=Sqr,
                   "ENC"=SNC,"Eqc"=Sqc,
                   "Em"=Sm,"EqrQ75"=SqrQ75,
                   "ENCQ75"=SNCQ75,"EqQ75"=SqQ75,
                   "EqrM"=SqrM,"ENCM"=SNCM,
                   "EqM"=SqM,"CImbeta"=CImbeta,
                   "CIq"=CIq,"CIqr"=CIqr,
                   "CINC"=CINC,"CIq75"=CIq75,
                   "CIqr75"=CIqr75,"CINC75"=CINC75,
                   "Eq"=Eq,
                   "GAs"=summary(GA))
    
    list.save(momentss, 'data.Rdata') 
  }
  