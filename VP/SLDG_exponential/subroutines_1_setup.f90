!<1><set up>
subroutine setup
    implicit none
    real :: alpha,beta,gamma
    real :: phi,sigma
 
    nx = 32*kkkk

    ny = nx
    nghost = 4
    n_moment = 6
    iqc = 1

    time_final = 0.2

    if(kkkk==0)then
        cfl = 1.
    else
        cfl = (kkkk)*5.
    endif

    cfl=kkkk*0.1
    cfl = 10.
    cfl = 0.5
    nghost = int(cfl) +6

    !iexprk =1

    !iexample = 2
    irk = 3

    !*****************************************************
    !VP setup

    iexp_case = 7

    if( iexp_case == 1 )then
        bt_con(2,1) = 1.
        iadditive(0) = 0
        iexprk = 1
    elseif(iexp_case == 2)then
        bt_con(2,1) = 0.5
        iadditive(0) = 0
        bt_con(3,1) = 0. ;        bt_con(3,2) = 1.
        iadditive(1) = 0
        iexprk = 2
    elseif( iexp_case == 3 )then
        bt_con(2,1) = (2.-sqrt(2.) )/2.
        iadditive(0) = 0
        bt_con(3,1) = (-2.*sqrt(2.) )/3. ;      bt_con(3,2) = 1 - bt_con(3,1)
        iadditive(1) = 0
        bt_con(4,1) =   0. ;                    bt_con(4,2) = 1- bt_con(2,1) ; bt_con(4,3) = bt_con(2,1)
        iadditive(2) = 0
        iexprk = 3
    elseif( iexp_case == 4 )then
        bt_con(2,1) = 1./2.
        iadditive(0) = 0
        bt_con(3,1) = -1.     ;      bt_con(3,2) = 2.;
        iadditive(1) = 0
        bt_con(4,1) = 1./12.  ;      bt_con(4,2) = 1./3. ;   bt_con(4,3) = -1./4.;
        iadditive(2) = 0
        bt_con(5,1) = 1./12.  ;      bt_con(5,2) = 1./3. ;   bt_con(5,3) = 5./12.;
        iadditive(3) = 1
        iexprk = 4
    elseif( iexp_case == 5 )then
        gamma = (3.+sqrt(3.) )/6.
        phi = 1./(6.*(2.*gamma-1.) )
        bt_con(2,1) = gamma
        iadditive(0) = 0
        bt_con(3,1) = gamma-1. ;      bt_con(3,2) = 2.*(1.-gamma) ;
        iadditive(1) = 0
        bt_con(4,1) = 0.; bt_con(4,2) = 1./2.-phi ;   bt_con(4,3) = 1./2.+phi ;
        iadditive(2) = 0
        bt_con(5,1) = 0.;      bt_con(5,2) = phi ;   bt_con(5,3) = -phi;
        iadditive(3) = 1
        iexprk = 4
    elseif( iexp_case == 6 )then
        alpha = 0.5
        beta = 1./6.
        gamma = ( 3.+sqrt(3.) )/6.
        sigma = (  alpha+beta*(1-2.*gamma)-1./3.    )/( 1-2.*gamma )
        bt_con(2,1) = gamma
        iadditive(0) = 0
        bt_con(3,1) = gamma-1.     ;      bt_con(3,2) = 2.*(1.-gamma) ;
        iadditive(1) = 0
        bt_con(4,1) = alpha  ;      bt_con(4,2) = beta ;   bt_con(4,3) = sigma ;
        iadditive(2) = 0
        bt_con(5,1) = -alpha  ;      bt_con(5,2) = 0.5-beta ;   bt_con(5,3) = 0.5-sigma ;
        iadditive(3) = 1
        iexprk = 4
    elseif(iexp_case == 7)then
        bt_con(2,1) = 1./3.
        iadditive(0) = 0
        bt_con(3,1) = 0. ;        bt_con(3,2) = 2./3.
        iadditive(1) = 0
        bt_con(4,1) = -1./12. ;   bt_con(4,2) = 0. ;      bt_con(4,3) = 3./4.
        iadditive(2) = 0
        iexprk = 3
    elseif( iexp_case == 8 )then
        bt_con(2,1) = 0.5;
        iadditive(0) = 0
        bt_con(3,1) = 0.0;   bt_con(3,2) = 0.5;
        iadditive(1) = 0
        bt_con(4,1) = -0.5;  bt_con(4,2) = 0.;  bt_con(4,3) = 1.;
        iadditive(2) = 0
        bt_con(5,1) = 0.25;  bt_con(5,2) = 1./6.; bt_con(5,3) = 1./6.; bt_con(5,4) = -1./12.;
        iadditive(3) = 0
        bt_con(6,1) = -1./12.; bt_con(6,2)=1./6.; bt_con(6,3) = 1./6.; bt_con(6,4) = 0.25
        iadditive(4) = 1
        iexprk = 5
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif

    kdg = 3
    call SLDG_poisson1d_paras

    ! write(*,*) 'i_case: i_case: 0 weak lan; 1 strong lan; 2 two strea; 3 two strea2'
    i_case = 0

    ! compatational domain
    if(i_case .eq. 0)then
        xleft = 0.
        xright = 4.*pi
        ybottom = -2.*pi
        ytop = 2.*pi
    elseif(i_case .eq. 1)then
        xleft = 0.
        xright = 4.*pi
        ybottom = -2.*pi
        ytop = 2.*pi
    elseif(i_case .eq. 2)then
        xleft = 0.
        xright = 4.*pi
        ybottom = -2.*pi
        ytop = 2.*pi

        ybottom = -10.
        ytop = 10.
    elseif(i_case .eq. 3)then
        xleft = 0.
        xright = 13.*pi
        ybottom = -2.*pi
        ytop = 2.*pi

    elseif(i_case .eq. 4)then
        xleft = 0.
        xright = 20.*pi/3.
        ybottom = -13.
        ytop = 13.
    endif


    end subroutine setup
!--------------------------------------------------~  
    

!<2><set dt>
subroutine setdt
    implicit none
    real :: temp_max2
 
    temp_max2 = 0.0

    do i = 1 , nx
        temp_max2 = max(temp_max2,abs(eE(i,0) ) )
    enddo

    dt = cfl/( ytop/dx + temp_max2/dy )

    if(time+dt>time_final) dt = time_final- time

    time = time +dt

end subroutine setdt
!--------------------------------------------------~     
    
    