!##############################################################################
! PROGRAM Church Optimizing Behavior
!
! copyright: Fabio Blasutto and David de la Croix
!            
!##############################################################################
include "param.f90"
!include "ToolboxKind.f90"



program Church

    ! modules
    use globals
    use toolbox


    ! set the fixed cost of censorship
    psig=linspace(0.1d0,0.2d0,NP) 


    do jj=1,NP

    !ge the psi
    psi=psig(jj)
    
    !solve the model
    call solve()

    enddo


   

    contains


    !solve the model
    subroutine solve()
 
       implicit none

        ! initialize grid
    l=0d0
    r=1d0
    m=linspace(l,r,NM)
    chi=1d0

    ! initial Guess on the Value Functions
    do i=1,NM
       
      Vc0(i)=utility(1d0-m(i))/(1-delta)
      Vr0(i)=utility(1d0-m(i))/(1-delta)
      
    enddo
    
    VrF0=Vr0-psi
    V0=max(Vc0,VrF0)
    
    ! iterate until policy function converges
    do iter = 1, itermax



        ! linear interpolation
        do i=1,NM

          ! get the share of revolutionary books in t+1
          mr(i)=(-1d0+Bbeta)*m(i)**2/(-1+m(i)*(2d0+m(i)*(-2d0+Bbeta)))
          mc(i)=(-1d0)*m(i)**2/(-1+m(i)*(2d0+m(i)*(-2d0)))
      
          
          ! compliant
          call linint_Equi_1(mc(i),l,r,NM-1,a,b,var)
          Vc_c(i)=Vc0(a+1)*var+Vc0(b+1)*(1d0-var)
          Vc_r(i)=Vr0(a+1)*var+Vr0(b+1)*(1d0-var)-psi
          Vc1(i)=max(Vc_c(i),Vc_r(i))

          ! revolutionary         
          call linint_Equi_1(mr(i),l,r,NM-1,a,b,var)
          Vr1(i)=Vr0(a+1)*var+Vr0(b+1)*(1d0-var)

        enddo

        ! get the new value functions
        do i=1,NM

            
            Vc(i)=utility(1d0-m(i))+delta*Vc1(i)
            Vr(i)=utility(1d0-m(i))+delta*Vr1(i)
    
        enddo
        VrF=Vr-psi


        ! get convergence level
        V=max(VrF,Vc)

        con_lev = maxval(abs(V - V0)/max(abs(V0), 1d-10))
        !write(*,'(i5,2x,f20.7)')iter, con_lev        
        
        ! check for convergence
        if(con_lev < sig)then  
         
 
         !check whether consistency with the data
         call linint_Equi_1(m1,l,r,NM-1,a,b,var)
          Vc_c(1)=Vc0(a+1)*var+Vc0(b+1)*(1d0-var)
          Vc_r(1)=Vr0(a+1)*var+Vr0(b+1)*(1d0-var)-psi
          
         if(Vc_c(1)>Vc_r(1))then
         call linint_Equi_1(m2,l,r,NM-1,a,b,var)

          Vc_c(1)=Vc0(a+1)*var+Vc0(b+1)*(1d0-var)
          Vc_r(1)=Vr0(a+1)*var+Vr0(b+1)*(1d0-var)-psi


          if(Vc_c(1)<Vc_r(1))then
 
          write(*,'(i5,2x,f20.7)')jj,psi  
             

          endif


         endif

         !go to next step
         return
        endif        

        ! update the value function
        V0=V
        VrF0=VrF
        Vr0=Vr
        Vc0=Vc
        
    enddo

   ! write(*,*)'No Convergence'

    end subroutine


    


    ! compute interesting and useful thresholds
    subroutine threshold()
 
       implicit none


        ! threshold over which m->1 even with censorship
        m_tR=1d0/(2d0-Bbeta)

        ! threshold over which m->1 with censorship
        m_tC=0.5d0


        ! lower threshold for switching to censorship
        trL=0d0
        do i=1,NM

          if(VrF(i)>Vc(i))then

          trL=m(i)
          goto 1
          
          endif

        enddo


        ! higher threshold for switching to censorship
1       trH=1d0
        do i=NM,1,-1
          
          if(VrF(i)>Vc(i))then

          trH=m(i)
          goto 2

          endif
 
        enddo

        
        ! threshold such that m(t+1)=trH
!2       m_trH=(2d0*trH-sqrt(4d0*trH**2d0-4d0*(2d0*trH-1d0)*trH))&
       !             /(2d0*(2d0*trH-1d0))


        ! threshold such that m(t+1)=m_tr
      !  m_trM=(2d0*m_tR-sqrt(4d0*m_tR**2d0-4d0*(2d0*m_tR-1d0)*m_tR))&
        !            /(2d0*(2d0*m_tR-1d0))

2 write(*,*)trL
write(*,*)trH
    end subroutine



    ! compute dynamics starting from different intial m
    subroutine dynamics()
 
       implicit none

       index=0  !index that says whether censorship is aready set up
       
       do i=1,NM
         
         ! get the staring point
         dyn(i,1)=m(i)

         if(i.eq.NM)then
             
           dyn(i,1)=trL-0.0001d0
                
         endif

         
         if(i.eq.NM-1)then
             
           dyn(i,1)=trH-0.0001d0
                
         endif
         
         do t=2,TT


           ! get m(t+1) accordin to the policy function
           if((dyn(i,t-1)>trL.and.dyn(i,t-1)<trH).or.(index(i).eq.1))then

             dyn(i,t)=(1d0-Bbeta)*chi*dyn(i,t-1)**2/&
                 ((1d0-Bbeta)*chi*dyn(i,t-1)**2+(1d0-dyn(i,t-1))**2)

          !   dyn(i,t)=chi*dyn(i,t-1)**2*(dyn(i,t-1)**Bbeta)/&
           !      (chi*dyn(i,t-1)**2*(dyn(i,t-1)**Bbeta)+(1d0-dyn(i,t-1))**2)
                 
             index(i)=1

           else

             dyn(i,t)=chi*dyn(i,t-1)**2/(chi*dyn(i,t-1)**2+(1d0-dyn(i,t-1))**2)

           endif


         enddo
       enddo



    end subroutine

    ! calibrate the model to the data
    subroutine calibration

        implicit none
        real*8:: KR_t(1:TF),KC_t(1:TF), beta_t(1:TF), m_t(1:TF), m_beta(1:TF), z_t(1:TF)
        real*8:: s_KR_t(1:TF),s_KC_t(1:TF), s_m_beta(1:TF), s_z_t(1:TF),s_m_t(1:TF),s_beta(1:TF)
        integer:: ts,IOUT
        
        ! get the initial data
        s_KC_t(1)= 3.260d0
        s_KR_t(1)= 5.970d0
        s_m_beta(1)=0.102d0
        
        s_KC_t(2)= 2.760d0
        s_KR_t(2)= 4.790d0
        s_m_beta(2)=0.062d0
        
        s_KC_t(3)= 3.050d0
        s_KR_t(3)= 4.990d0
        s_m_beta(3)=0.037d0

        s_z_t=s_KR_t/s_KC_t
        s_m_t=s_KR_t/(s_KR_t+s_KC_t)
        s_beta=s_m_beta/s_m_t


        ! get initial conditions
        KC_t(1)=s_KC_t(1)
        KR_t(1)=s_KR_t(1)
        m_beta(1)=s_m_beta(1)
        m_t(1)=KR_t(1)/(KR_t(1)+KC_t(1))
        do t=1,TF
        beta_t(t)=m_beta(1)/m_t(1)
        enddo

        z_t(1)=KR_t(1)/KC_t(1)

        !call shock()
        

        ! simulation of the model

        ! censorship is zero at the beginnig
        index(1)=1
        
        do t=2,TF

           ! get m(t+1) accordin to the policy function
           if((m_t(t-1)>trL.and.m_t(t-1)<trH).or.(index(1).eq.1))then

             m_t(t)=(1d0-beta_t(t))*chi*m_t(t-1)**2/&
                 ((1d0-beta_t(t))*chi*m_t(t-1)**2+(1d0-m_t(t-1))**2)

             z_t(t)=(1d0-beta_t(t))*chi*(z_t(t-1))**2


            !m_t(t)=chi*m_t(t-1)**2*(m_t(t-1)**Bbeta)/&
             !    (chi*m_t(t-1)**2*(m_t(t-1)**Bbeta)+(1d0-m_t(t-1))**2)
                 
             index(1)=1

           else

             m_t(t)=chi*m_t(t-1)**2/(chi*m_t(t-1)**2+(1d0-m_t(t-1))**2)
             z_t(t)=chi*(z_t(t-1))**2

           endif

           ! m*beta+z
           m_beta(t)=m_t(t)*beta_t(t)


        enddo

        ! compare the two creating a txt file
        IOUT=3
        OPEN(IOUT,FILE='OUTPUT.TXT'  ,STATUS='UNKNOWN')

        WRITE(IOUT,1)
   1    FORMAT(&
         '###############################################################',/,&
         '###     This file compares data with model dynamics ###########',/,&
         '###############################################################',/,&
         '                                                               ')
          write(IOUT,'(a,a,a,a,a,a,a,a,a)')&
          'T','     ','z(sim),z(data)',&
                              '       ','m(sim),m(data)',&
                                                   '          ','m*b(sim),m*b(data)'&
                                                   ,'         ','beta(sim),b(data)'
          write(IOUT,'(a)')&
         '                                                                 '

         do t=1,TF
           write(IOUT,'(I1,a,f6.4,a,f6.4,a,f6.4,a,f6.4,a,f6.4,a,f6.4,a,f6.4,a,f6.4)')&
          t,&
              '     ',z_t(t),',',s_z_t(t),&
                              '        ',m_t(t),',',s_m_t(t),&
                                                    '            ',m_beta(t),',',s_m_beta(t)&
                                                    ,'               ',beta_t(t),',',s_beta(t)


         enddo

            
        CLOSE( IOUT )
        


    end subroutine calibration


    ! for creating output plots.
    subroutine output()

        implicit none


        ! end timer
        call toc()

        write(*,*)'Convergence!'


        ! plots here
        indice=1d0




      ! plot the dyncamis
      time=linspace(1d0,real(TT,2),TT)
      
      do i=1,NM,(NM/100)-1
 
        call plot(time(1:TT), dyn(i,1:TT),legend='')

      enddo
      call plot(time(1:TT), dyn(NM,1:TT),legend='trL',marker=2)
      call plot(time(1:TT), dyn(NM-1,1:TT),legend='trH',marker=3)
      call execplot(ylabel='m', xlabel='time',&
      title='Dynamics of m',&
      filename='dyn', filetype='pdf', output='testdata')

!$$$$$$       !Plot \beta over time
!$$$$$$         call plot(time(1:TT), (1d0-dyn(1,1:TT)**Bbeta)*dyn(1,1:TT),legend='1')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(200,1:TT)**Bbeta)*dyn(200,1:TT),legend='2')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(300,1:TT)**Bbeta)*dyn(300,1:TT),legend='3')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(400,1:TT)**Bbeta)*dyn(400,1:TT),legend='4')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(500,1:TT)**Bbeta)*dyn(500,1:TT),legend='5')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(600,1:TT)**Bbeta)*dyn(600,1:TT),legend='6')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(700,1:TT)**Bbeta)*dyn(700,1:TT),legend='7')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(800,1:TT)**Bbeta)*dyn(800,1:TT),legend='8')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(900,1:TT)**Bbeta)*dyn(900,1:TT),legend='9')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(1000,1:TT)**Bbeta)*dyn(1000,1:TT),legend='trL',marker=2)
!$$$$$$         call plot(time(1:TT), (1d0-dyn(999,1:TT)**Bbeta)*dyn(999,1:TT),legend='trH',marker=3)
!$$$$$$ 
!$$$$$$    
!$$$$$$ 
!$$$$$$       call execplot(ylabel='m beta', xlabel='time',&
!$$$$$$       title='Dynamics of m beta',&
!$$$$$$       filename='dyn', filetype='pdf', output='testdata')
!$$$$$$ 
!$$$$$$ 
!$$$$$$             !Plot \beta over time
!$$$$$$         call plot(time(1:TT), (1d0-dyn(1,1:TT)**Bbeta),legend='1')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(200,1:TT)**Bbeta),legend='2')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(300,1:TT)**Bbeta),legend='3')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(400,1:TT)**Bbeta),legend='4')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(500,1:TT)**Bbeta),legend='5')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(600,1:TT)**Bbeta),legend='6')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(700,1:TT)**Bbeta),legend='7')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(800,1:TT)**Bbeta),legend='8')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(900,1:TT)**Bbeta),legend='9')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(1000,1:TT)**Bbeta),legend='trL')
!$$$$$$         call plot(time(1:TT), (1d0-dyn(999,1:TT)**Bbeta),legend='trH',marker=3)
!$$$$$$ 
!$$$$$$    
!$$$$$$ 
!$$$$$$       call execplot(ylabel='beta', xlabel='time',&
!$$$$$$       title='Dynamics of beta',&
!$$$$$$       filename='dyn', filetype='pdf', output='testdata')


      ! value funtion
      call plot(m, Vc,legend='Vc') 
      call plot(m, VrF,legend='VrF')
      call execplot(ylabel='Utility', xlabel='Initial m',&
      title='Value function and Initial Conditions',&
      filename='valuef', filetype='pdf', output='testdata')

      ! for the proof
      diff=Vc-VrF

      call plot(m, diff,legend='Vc-VrF')
      call plot(indice*m_trM, Vc,legend='crucial threshold') 
      call plot(indice*m_trH, Vc,legend='crucial threshold 2') 
      call execplot(ylabel='Vc-VrF', xlabel='Initial m',&
      title='Value function and Initial Conditions',&
      filename='study', filetype='pdf', output='testdata')

      
              
    
        ! quit program
        stop 
    
    end subroutine

    ! utility function
    function utility(c)
       
         implicit none
         
         real*8, intent(in):: c
         real*8 :: utility
       
         
         !utility=log(max(c,0.000000000000000000000001d0))
         utility=c
       
       end function utility


end program