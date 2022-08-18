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
  
  # get z0 first
  #v=uniroot(f,c(0.00001,up),tol=0.00000001,maxiter=10000,the=thet,k=1,
  #          pr=pe,eqt=eqtt,eqrt=eqrtt)
  #z0=(eqrtt/v$root)^(1/thet)
  z0=(eqrtt/eqct)^(1/thet)
  if(t>2){
    z0=pe/(1)*(z0*(1)/pe)^(2)
    pe/(1-beta)*(z0*(1-beta)/pe)^(2^(t-2))}else if(t==2){pe/(1)*(z0*(1)/pe)^(2^(t-1))}else{z0}}


mt<-function(beta,pr,t,up,thet,eqct,eqrtt) zt(beta,pr,t,up,thet,eqct,eqrtt)/(pr+zt(beta,pr,t,up,thet,eqct,eqrtt))

mbetat<-function(beta,pr,t,up,thet,eqct,eqrtt) beta*impbeta*mt(beta,pr,t,up,thet,eqct,eqrtt)

qrs<-function(beta,pe,t,up,mu,pre,thet,eqct,eqrtt) if(t<=2){gamma(1-thet)*(mu*((1)*(pre/gamma(1-thet))^(1/thet)*
                                                                                 (mt(0,pe,t-1,up,thet,eqct,eqrtt))))^thet}else
                                                                                 {gamma(1-thet)*(mu*((1-beta)*(pre/gamma(1-thet))^(1/thet)*
                                                                                                       (mt(beta,pe,t-1,up,thet,eqct,eqrtt))))^thet} 


qcs<-function(beta,pe,t,up,mu,pre,thet,eqct,eqrtt) if(t<=2){gamma(1-thet)*((mu*((pre/gamma(1-thet))^(1/thet)*
                                                                                  (1-mt(0,pe,t-1,up,thet,eqct,eqrtt)))))^thet}else
                                                                                  {gamma(1-thet)*((mu*((pre/gamma(1-thet))^(1/thet)*
                                                                                                         (1-mt(beta,pe,t-1,up,thet,eqct,eqrtt)))))^thet}             


qs<-function(beta,pe,t,up,mu,pre1,pre2,thet,eqct,eqrtt) if(t<=2){qrs(0,pe,t,up,mu,pre1,thet,eqct,eqrtt)*mt(0,pe,t,up,thet,eqct,eqrtt)+
    (1-mt(0,pe,t,up,thet,eqct,eqrtt))*qcs(0,pe,t,up,mu,pre2,thet,eqct,eqrtt)}else
    {qrs(beta,pe,t,up,mu,pre1,thet,eqct,eqrtt)*mt(beta,pe,t,up,thet,eqct,eqrtt)+
        (1-mt(beta,pe,t,up,thet,eqct,eqrtt))*qcs(beta,pe,t,up,mu,pre2,thet,eqct,eqrtt)}

nc<-function(beta,pe,t,up,mu,pre1,pre2,thet,eqct,eqrtt) if(t<=2){qrs(0,pe,t,up,mu,pre1,thet,eqct,eqrtt)*
    (((1)*mt(0,pe,t,up,thet,eqct,eqrtt))/(1))+
    ((1-mt(0,pe,t,up,thet,eqct,eqrtt))/(1))*
    qcs(0,pe,t,up,mu,pre2,thet,eqct,eqrtt)}else
    {qrs(beta,pe,t,up,mu,pre1,thet,eqct,eqrtt)*
        (((1-beta)*mt(beta,pe,t,up,thet,eqct,eqrtt))/(1-beta*mt(beta,pe,t,up,thet,eqct,eqrtt)))+
        ((1-mt(beta,pe,t,up,thet,eqct,eqrtt))/(1-beta*mt(beta,pe,t,up,thet,eqct,eqrtt)))*
        qcs(beta,pe,t,up,mu,pre2,thet,eqct,eqrtt)}




GA <- ga(type = "real-valued", 
         fitness = ff, 
         lower = c(0.12,1.7,1.0,0.05,3.5-0.35,Eqr[1]-0.4),
         upper = c(0.25,2.5,4  ,0.45,3.5+0.4 ,Eqr[1]+0.4),
         popSize = 100, maxiter = 180, run = 50,
         elitism = max(1, round(200*0.15)),
         monitor=FALSE,
         optim=FALSE,
         seed=1,
         optimArgs = list(method = "L-BFGS-B",poptim = 0.2,pressel = 0.5,
                          control = list(fnscale = -1, maxit = 20)))
