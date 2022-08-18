!##############################################################################
! MODULE globals
!
! copyright: Fabio Blasutto and David de la Croix

!##############################################################################
include "ToolboxKind.f90"

module globals


    use toolbox


    implicit none


    
    ! model parameters    
    real*8, parameter :: delta = 0.06d0   !The discount factor  (0.99**30)  
    real*8, parameter :: Bbeta = 0.18d0     !Max Censorship
    real*8, parameter :: theta   = 0.13d0   !Productivity of ideas
    real*8, parameter :: mu   = 2.47d0       !TFP growth
    real*8, parameter :: p   = 0.41d0       !relative price of revolutionary books
    real*8, parameter :: Q1   = 4.851495d0       !initial condition
    real*8, parameter :: Qr1   = 6.883995d0       !initial condition
    real*8, parameter :: m1   = 0.501d0!0.5056186d0       !initial condition
    real*8, parameter :: m2   = 0.52d0!0.5112359d0       !initial condition
    
    real*8:: chi  
    real*8 :: alpha
    real*8 :: gammap 

    ! numerical parameters
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000
    integer, parameter :: TF = 5




    ! counter variables
    integer :: i,iter

    ! for interopolation
    integer :: a,b
    real*8 :: var
    
    ! grids
    integer, parameter :: NM = 1000
    real*8 :: m(1:NM),mc(1:NM),mr(1:NM),l,r
    real*8:: psi 

    ! all the psi to be checksd
    integer, parameter :: NP = 1000
    real*8:: psig(1:NP)
    integer::jj

    ! value Functions
    real*8 :: Vc(1:NM),Vc1(1:NM),Vc_c(1:NM),Vc_r(1:NM),Vc0(1:NM)
    real*8 :: Vr(1:NM),Vr1(1:NM),VrF(1:NM),Vr0(1:NM),VrF0(1:NM),Vr_c(1:NM),Vr_r(1:NM)
    real*8 :: V(1:NM),V0(1:NM)

    ! other arrays for the plots
    real*8 :: m_tR, m_tC, m_trM, trH, trL, m_trH
    real*8 :: indice(1:NM), diff(1:NM)

    ! variables to numerically determine policy function
    real*8 :: con_lev

    ! arrays for the dynamics over time
    integer :: t,index(1:NM)
    integer, parameter :: TT=10
    real*8 :: dyn(1:NM,1:TT),time(1:TT)


    contains


    !##############################################################################
    !EQUALLY SPACED GRID
    !##############################################################################    
    function linspace(xmin,xmax,n)
    
        real*8:: xmin, xmax, linspace(1:n)
        integer:: n, i
        
        linspace = (/ (xmin + (xmax-xmin)*(i-1)/(n-1),i=1,n) /)
        
    end function linspace



  ! subroutine that compute the variance of the simulated distribution
  subroutine frec_var_s(KKR,KKC,Theta_s,var)
   
     implicit none
     
     real*8,intent(in) ::KKR,KKC, Theta_s
     real*8,intent(out) ::var
     real*8,allocatable :: shockm(:),shocks(:), smp(:)
     real*8 :: ms,mean
     integer ::qq,NB

NB=1000000
allocate(shocks(1:NB),shockm(1:NB),smp(1:NB) ) 
     !Open the file
  !   open (unit=12,file='shocks.txt',form="formatted",action='write',status="replace")
     
     !Read the shock
   !  do qq = 1 , NB
!
   !     read(12,*) (shocks(qq))
     
   !  enddo

  !   close(12)

 call simulate_uniform_n(shocks(1:NB),0d0,1d0)



     !Open the file
  !   open (unit=11,file='shockm.txt',form="formatted",action='write',status="replace")

     
     !Read the shock
   !  do qq = 1 , NB

   !     read(11,*) (shockm(qq))
     
  !   enddo

   !  close(11)

 call simulate_uniform_n(shockm(1:NB),0d0,1d0)
     ! get the realizations

     ms=KKR/(KKR+KKC)
     shockm=1d0
     do qq =1,NB

        if(shockm(qq).le.ms)then

           call frechet_cdf_inv ( shocks(qq), 1d0/Theta_s, KKC**Theta_s , smp(qq) )

        else

           call frechet_cdf_inv ( shocks(qq), 1d0/Theta_s, KKR**Theta_s, smp(qq) )

        endif




      enddo


     ! compute the variance
   mean = 0d0                           ! compute mean
   DO qq = 1, NB
      mean = mean + smp(qq)
   END DO

   mean = mean / NB
    var = 0d0                       ! compute variance
   DO qq = 1, NB
      var = var + (smp(qq) - mean)**2
   END DO
   var = var / (NB - 1d0)


  end subroutine

  subroutine frechet_cdf_inv ( cdf, alpha,s, x )
  
  !*****************************************************************************80
  !
  !! FRECHET_CDF_INV inverts the Frechet CDF.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 September 2008
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) CDF, the value of the CDF.
  !    0.0 <= CDF <= 1.0.
  !
  !    Input, real ( kind = 8 ) ALPHA, the parameter.
  !    It is required that 0.0 < ALPHA.
  !
  !    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
  !
    implicit none
  
    real*8,intent(in)  :: alpha,s
    real*8,intent(in)  :: cdf
    real*8,intent(out) :: x
  
    if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FRECHET_CDF_INV - Fatal error!'
      write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
      stop 1
    end if
  
    if ( alpha <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FRECHET_CDF_INV - Fatal error!'
      write ( *, '(a)' ) '  ALPHA <= 0.0.'
      stop 1
    end if
  
   ! if ( cdf == 0.0D+00 ) then
   !   x = 0.0D+00
    !else
     ! x =  s*(( - 1.0D+00 / log ( cdf ) ) ** ( 1.0D+00 / alpha ))
     x =  s*(-log(cdf)) ** (-1.0D+00 / alpha )
   ! end if
  
    return
  end subroutine

   

end module
