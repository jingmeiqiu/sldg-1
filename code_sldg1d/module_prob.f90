    module module_prob
    implicit none

    integer :: kkkk
    integer :: nx
    integer :: nk
    integer :: nghost
    
    integer :: iexample
    
    real :: xleft,xright
    real :: gau(1:6,1:2)
 
    !********************
    ! evolution related variables
    real :: time_final
    real :: cfl
    !********************
    
    contains

    subroutine setup
    implicit none
    real :: pi
    
    pi = 4.*atan(1.)
    !*******************************************
    ! the number of elements
    nx = 10*2**kkkk
    !*******************************************
    ! computational domain
    xleft = 0.
    xright = 2.*pi
    !*******************************************
    ! the number of ghost elements
    nghost = 11
    !*******************************************
    ! please set up here
    ! cfl is for CFL number;
    ! time_final is the final print time;
    ! nk is the number of p^{nk} DG solution space;
    cfl = 0.18
    time_final = 1.
    nk = 2
    
    iexample = 2
 

    !*******************************************
    if(nk .eq. 0)then
        gau(1,1) = 0.
        gau(1,2) = 1.
    endif
    if(nk .eq. 1) then
        !  the points of 4th order Gauss quadrature
        gau(1,1)=-sqrt(1./3.)*0.5
        gau(2,1)=sqrt(1./3.)*0.5

        gau(1,2)=0.5
        gau(2,2)=0.5
    endif
    if(nk .eq. 2) then
        !  the points of 6th order Gauss quadrature
        gau(1,1)=-sqrt(0.6)*0.5
        gau(2,1)=0.
        gau(3,1)=sqrt(0.6)*0.5
        !   coefficients of 6th order Gauss quadrature
        gau(1,2)=5./18.
        gau(2,2)=4./9.
        gau(3,2)=5./18.

    endif
    if( nk .eq. 3) then
        ! 8th order Gaussian nodes in [-1/2,1/2]
        gau(1,1)=-sqrt( 3./7.+2./7.*sqrt(1.2) )*0.5
        gau(2,1)=-sqrt( 3./7.-2./7.*sqrt(1.2) )*0.5
        gau(3,1)= sqrt( 3./7.-2./7.*sqrt(1.2) )*0.5
        gau(4,1)=sqrt( 3./7.+2./7.*sqrt(1.2) )*0.5
        !  8th order Gaussian weights in [-1/2,1/2]
        gau(1,2)=( 18.-sqrt(30.) )/36.*0.5
        gau(2,2)=( 18.+sqrt(30.) )/36.*0.5
        gau(3,2)=( 18.+sqrt(30.) )/36.*0.5
        gau(4,2)=( 18.-sqrt(30.) )/36.*0.5
    endif
    if( nk .eq. 4) then
        !  10th order Gaussian nodes in [-1/2,1/2]
        gau(1,1)= -1./3.*sqrt( 5.+2.*sqrt(10./7.) )*0.5

        gau(2,1)= -1./3.*sqrt( 5.-2.*sqrt(10./7.) )*0.5
        gau(3,1)= 0.
        gau(4,1)= 1./3.*sqrt( 5.-2.*sqrt(10./7.) )*0.5

        gau(5,1)= 1./3.*sqrt( 5.+2.*sqrt(10./7.) )*0.5
        ! 10th order Gaussian weights in [-1/2,1/2]
        gau(1,2)= ( 322.-13.*sqrt(70.) )/900.*0.5

        gau(2,2)= ( 322.+13.*sqrt(70.) )/900.*0.5
        gau(3,2)= 128./225.*0.5
        gau(4,2)= ( 322.+13.*sqrt(70.) )/900.*0.5

        gau(5,2)= ( 322.-13.*sqrt(70.) )/900.*0.5
    endif

    end subroutine setup
    !************************************
    real function fun_init(x)
    implicit none
    real,intent(in) :: x

    if(iexample==1)then
        fun_init = sin( x )
    elseif( iexample==2 )then
        fun_init = 1.
    endif
    end function fun_init
    
    !************************************
    real function ax(x)
    implicit none
    real,intent(in) :: x

    if(iexample==1)then
        ax = 1.
    elseif(iexample ==2)then
        ax = sin(x)
    endif

    end function ax
    !*******************************
    real function exact(x,t)
    implicit none
    real,intent(in) :: x,t

    if(iexample==1)then
        exact = sin( x-t )
    elseif( iexample==2 )then
        exact = sin(  2.*atan( exp(-t)*tan(x/2.)  ) )/sin(x)
    endif

    end function exact
    
    
    end module module_prob