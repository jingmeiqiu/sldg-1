    subroutine setup
    implicit none


    nx = 20*2**(kkkk-1)
    !nx = 4
    ny = nx
    nghost = 5
    n_moment = 6
    iqc = 1

    time_final = 1.5
    cfl = 2 

    iexample = 1
    irk = 3
    irelative = 1

    idebug = 0

    if(iexample==1)then
        xleft = -pi
        xright = pi
        ybottom = -pi
        ytop = pi
    elseif(iexample==2)then
        xleft = -pi
        xright = pi
        ybottom = -pi
        ytop = pi
    elseif(iexample ==3)then
        xleft = -pi*2.
        xright = pi*2.
        ybottom = -pi*2.
        ytop = pi*2.
    endif



    



    dx = (xright-xleft)/nx
    dy = (ytop-ybottom)/ny

    end subroutine setup

    ! problem
    real function  vel_x( x,y,t )
    implicit none
    real,intent(in) :: x,y,t

    if( iexample == 1 )then
        vel_x = 1.
    elseif( iexample == 2 )then
        vel_x = -cos(x/2)**2 * sin(y) * cos(pi*t/time_final )*pi
    elseif( iexample == 3 )then
        vel_x = -y
    endif

    end function vel_x

    real function  vel_y( x,y,t )
    implicit none

    real,intent(in) :: x,y,t

    if( iexample == 1 )then
        vel_y = 1.
    elseif( iexample == 2 )then
        vel_y = sin(x)*cos(y/2.)**2 * cos(pi*t/time_final )*pi
    elseif(iexample == 3)then
        vel_y = x
    endif

    end function vel_y