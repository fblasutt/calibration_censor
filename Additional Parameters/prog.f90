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
    psig=linspace(0.00356d0,0.00359d0,NP) !linspace(0.000356d0,0.00359d0,NP) !


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
         !write(*,*)1, Vc0(a+1)*var+Vc0(b+1)*(1d0-var)- Vr0(a+1)*var+Vr0(b+1)*(1d0-var)
         !if(Vc_c(1)>Vc_r(1))then
         call linint_Equi_1(m2,l,r,NM-1,a,b,var)

          Vc_c(1)=Vc0(a+1)*var+Vc0(b+1)*(1d0-var)
          Vc_r(1)=Vr0(a+1)*var+Vr0(b+1)*(1d0-var)-psi
         
         ! write(*,*)2, Vc0(a+1)*var+Vc0(b+1)*(1d0-var)- Vr0(a+1)*var+Vr0(b+1)*(1d0-var)

         if(Vc_c(1)<Vc_r(1))then
 
         write(*,'(i5,2x,f20.7)')jj,psi
        
         endif


         !endif

         !go to next step
         return
        endif        

        ! update the value function
        V0=V
        VrF0=VrF
        Vr0=Vr
        Vc0=Vc
        
    enddo

    write(*,*)'No Convergence'

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