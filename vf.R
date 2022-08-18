#####################################
#Compute value function of censoring
######################################

#Initial setup
itermax<-10000                       #max inter
psi<-0.003420013#0.003428373#
sig<-1e-10                           #max error
delta<-0.06                          #discount factor
u<-function(c){1-c}                  #church's utility
mg<-seq(from = 0, to = 1, by =.0001) #grid of m
vr0<-u(mg)/(1-delta) 
vc0<-u(mg)/(1-delta) #Initial guess
mpb<-function(x){(1-betax)*x^2/
                (1-x*((betax-2)*x+2))}#evolution of m with censorship
mp0<-function(x){         x^2/
                (1-x*((-2)*x+2))}    #evolution of m w/o censorship
vr1<-vr0
vc1<-vc0
v0<-vr1
v1<-v0
mgb<-mpb(mg)
mg0<-mp0(mg)

#Get the value function vr
for(i in 1:itermax){
  
  #Get new guess
  vr1<-u(mg)+delta*approx(mg,vr0,mgb)$y 
  vr1[length(mg)]<-u(mg[length(mg)])/(1-delta)
  if(max(abs(vr1-vr0)/max(abs(vr0), 1e-10))<sig){break}else{vr0<-vr1}
  
  
}

for(i in 1:itermax){
  
  #Get new guess
  vc1<-u(mg)+delta*approx(mg,v0,mg0)$y
  v1<-pmax(vr1-psi,vc1)
  #Get tolerance
  if(max(abs(v1-v0)/max(abs(v0), 1e-10))<sig){break}else{
    vc0<-vc1
    v0<-v1}
}

#Now get the thrat points
SmC2<-mt(0,exp(GA@solution[1,2]),2,1,GA@solution[1,4],Sqc[1],Sqr[1]) 
SmC3<-mt(0,exp(GA@solution[1,2]),3,1,GA@solution[1,4],Sqc[1],Sqr[1])
SmC4<-mt(0,exp(GA@solution[1,2]),4,1,GA@solution[1,4],Sqc[1],Sqr[1])
psi1<-delta/(1-delta)*(approx(mg,vr1,c(mpb(SmC2)))$y-approx(mg,vr1,c(SmC3))$y)/(approx(mg,vr1,1/(2-betax))$y -approx(mg,vc1,1/(2-betax))$y)
psi2<-delta/(1-delta)*(approx(mg,vr1,c(mpb(SmC3)))$y-approx(mg,vr1,c(SmC4))$y)/(approx(mg,vr1,1/(2-betax))$y -approx(mg,vc1,1/(2-betax))$y)
