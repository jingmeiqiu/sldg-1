
    subroutine init
    implicit none
    real :: utemp
    integer :: kk,lx,ly
    !******************************************************
    !subroutine mesh_generator
    !implicit none
    ! x(i) ----> center
    ! y(j) ----> center
    DO I = 1 - nghost, nx + nghost
        X(I) = XLEFT + (I-0.5) * DX
    ENDDO

    DO J = 1 - nghost , NY + nghost
        Y(J) = Ybottom + (J-0.5) * DY
    ENDDO

    DO I = 1 - nghost, nx + 1  + nghost
        Xgrid(I) = XLEFT + (I-1.) * DX
    ENDDO

    DO J =  1-nghost, NY+ 1 + nghost
        ygrid(J) = Ybottom + (J-1.) * DY
    ENDDO

    do i = 1,nx+1
        do j = 1,ny+1
            vertex(i,j)%coor(1) = xleft + (i-1) * dx
            vertex(i,j)%coor(2) = ybottom + (j-1)*dy
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny+1
            nodex(i,j)%coor(1) = xleft + (i-0.5) * dx
            nodex(i,j)%coor(2) = ybottom + (j-1)*dy
        enddo
    enddo

    do i = 1,nx+1
        do j = 1,ny
            nodey(i,j)%coor(1) = xleft + (i-1) * dx
            nodey(i,j)%coor(2) = ybottom + (j-0.5)*dy
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny
            nodec(i,j)%coor(1) = xleft + (i-0.5) * dx
            nodec(i,j)%coor(2) = ybottom + (j-0.5)*dy
        enddo
    enddo

    do i = 1 , nx
        do j = 1 ,ny
            element(i,j)%vertex1 => vertex(i,j)
            element(i,j)%vertex2 => vertex(i+1,j)
            element(i,j)%vertex3 => vertex(i+1,j+1)
            element(i,j)%vertex4 => vertex(i,j+1)

            if(n_moment==6)then
                element(i,j)%vertex5 => nodex(i,j)
                element(i,j)%vertex6 => nodey(i+1,j)
                element(i,j)%vertex7 => nodex(i,j+1)
                element(i,j)%vertex8 => nodey(i,j)
                element(i,j)%vertex9 => nodec(i,j)
            endif

        enddo
    enddo
    !end subroutine mesh_generator
    !******************************************************

    do i = 1 , nx
        do j = 1 , ny
            do kk = 1 , n_moment
                utemp =0.0
                do lx=1,6
                    do ly=1,6
                        utemp = utemp +exact(x(i)+xg(lx)*dx ,y(j)+xg(ly)*dy,0.) *fphi(kk,xg(lx),xg(ly) )*wg(lx)*wg(ly)
                    enddo
                enddo
                element(i,j)%umodal(kk)= utemp * ai(kk)

            enddo
        enddo
    enddo


    end subroutine init

    real function exact(x,y,t)
    implicit none
    real, intent(in) :: x,y,t
    real :: rb,rb0

    if( iexample == 1 )then
        exact = sin(x+y-2*t)
    elseif( iexample == 2 )then
        rb = sqrt( (x-0.15/0.5*pi)**2 + (y)**2  )
        rb0 = 0.3*pi
        if(rb < rb0)then
            exact =  rb0 * cos( pi*rb/(2.*rb0) )**6
        else
            exact = 0.
        endif
    elseif( iexample == 3 )then
        exact = exp(-x**2-y**2)
    endif

    return
    end function exact