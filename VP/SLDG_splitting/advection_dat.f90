    subroutine advection_dat
    implicit none
    real :: temp_p
    integer :: lx,ly

    call smooth_sin

    do j = 1 , nel_y
        do i = 1 , nel_x
            elem(i,j)%psi(1:n_g,1:n_g) = fel(1:n_g,1:n_g,i,j)
        enddo
    enddo
    
    if(i_case .eq. 0)then
        !pint = 0.
        !do i = 1, nx
        !    do j = 1 , ny
        !        temp_p = 0.
        !        do lx = 1 , mmp
        !            do ly = 1 ,mmp
        !                temp_p = temp_p + dx*wq(lx)*wq(ly)*exact( x(i)+dx*vgau(lx,1),y(j)+dy*vgau(ly,1),0. )*dy
        !            enddo
        !        enddo
        !        pint = pint + temp_p
        !    enddo
        !enddo
        !pint = -pint/(xright-xleft)

        pint = -1.
    elseif(i_case .eq. 1)then
        !pint = 0.
        !do i = 1, nx
        !    do j = 1 , ny
        !        temp_p = 0.
        !        do lx = 1 , mmp
        !            do ly = 1 ,mmp
        !                temp_p = temp_p + dx*wq(lx)*wq(ly)*exact( x(i)+dx*vgau(lx,1),y(j)+dy*vgau(ly,1),0. )*dy
        !            enddo
        !        enddo
        !        pint = pint + temp_p
        !    enddo
        !enddo
        !pint = -pint/(xright-xleft)

        pint = -1.
    elseif(i_case > 1)then
        pint = 0.
        do i = 1, nel_x
            do j = 1 , nel_y
                temp_p = 0.
                do lx = 1 , n_g
                    do ly = 1 , n_g
                        temp_p = temp_p + dx*w_g(lx)*w_g(ly)* u0( x(i)+dx*x_g(lx) ,y(j)+dy* x_g(ly)  )*dy
                    enddo
                enddo
                pint = pint + temp_p
            enddo
        enddo
        pint = -pint/(xright-xleft)

    endif

    end subroutine advection_dat
    !**************************************************************************************
    subroutine smooth_sin
    implicit none
    integer :: n1,n2

    do j = 1 , nel_y
        do i = 1 , nel_x
            !
            !do n2 = 1 , n_gl
            !    do n1 = 1 ,n_gl
            !        !
            !        fel( n1,n2,i,j ) = sin( glx(n1,n2,i,j) + gly(n1,n2,i,j)  )
            !        fel( n1,n2,i,j ) = sin( glx(n1,n2,i,j)  )
            !
            !    enddo
            !enddo

            do n2 = 1 , n_g
                do n1 = 1 ,n_g
                    !
                    !fel( n1,n2,i,j ) = sin( gx(n1,n2,i,j) + gy(n1,n2,i,j)  )

                    fel( n1,n2,i,j ) = u0(  gx(n1,n2,i,j),gy(n1,n2,i,j) )
                    !fel( n1,n2,i,j ) = sin( gx(n1,n2,i,j)  )
                    !fel( n1,n2,i,j ) = sin(  gy(n1,n2,i,j)  )

                enddo
            enddo
        enddo
    enddo


    end subroutine smooth_sin
    !***************************************************
    !  initial condition
    real function u0(x,y)
    !common /const/ pi, cfl, dt, v, T, dxmin, dymin
    implicit none
    real,intent(in) :: x,y
    real :: rb,rb0
    real :: pk
    real :: alpha

    real :: pu,vth,vth2
    real :: fbt

    if(i_case .eq. 0)then
        alpha = 0.01
        pk = 0.5
        u0 = 1./sqrt(2*pi) * (1. +alpha*( cos(0.5*x) ) )*exp(- y**2/2. )
    elseif(i_case .eq. 1)then
        alpha = 0.5
        pk = 0.5
        u0 = 1./sqrt(2*pi) * (1. +alpha*( cos(0.5*x) ) )*exp(- y**2/2. )
    elseif(i_case .eq. 2)then
        alpha = 0.01
        pk = 0.5
        u0 = 2./7./sqrt(2.*pi) *(1+5.*y**2 ) * ( 1+ alpha*( ( cos(2.*pk*x ) + cos(3.*pk*x) )/1.2 +cos(pk*x)  ) )*exp(-y**2/2.)

    elseif(i_case .eq. 3)then
        pu = 0.99
        vth = 0.3
        pk = 2./13.
        vth2 = 2*vth*vth
        u0 = 0.5/vth/sqrt(2.*pi)*(exp(-( y -pu)**2./vth2)+exp(-( y +pu)**2./vth2))*(1.+0.05*cos(pk*x ))
    elseif(i_case .eq. 4)then
        pu = 4.5
        vth = 0.5
        pk = 0.3
        vth2 = 2.*vth*vth
        fbt = 9./(10.*sqrt(2.*pi) ) *exp( -y**2/2. ) + 2./(10.*sqrt(2.*pi) ) * exp( -(y-pu)**2/vth2 )
        u0 = fbt * (1+0.04*cos(pk*x) )
    endif
    !u0 = 1./sqrt(2.*pi)*y*y*exp(-0.5*y*y)*(1+0.05*cos(0.5*x))


    end function u0
    !*****************************************************
    !***********************************************
    real function exact(x,y,t)
    implicit none
    real, intent(in) :: x,y,t
    ! real, dimension(1:4) :: exact
    real :: rb,rb0

    ! exact=  sin(x+y-2.*t)
    !exact =  exp(- x**2-y**2 )


    ! cosine

    rb = sqrt( (x-0.15/0.5*pi)**2 + (y)**2  )
    rb0 = 0.3*pi
    if(rb < rb0)then
        exact =  rb0 * cos( pi*rb/(2.*rb0) )**6
    else
        exact = 0.
    endif

    return
    end function exact
    !***************************************************************************************
 