!<><Initialization>
subroutine init
    implicit none
    real :: utemp
    integer :: kk,lx,ly
    integer :: k


    dx = (xright-xleft)/nx
    dy = (ytop-ybottom)/ny
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
            do k = 0 , iexprk
                element(i,j,k)%vertex1 => vertex(i,j)
                element(i,j,k)%vertex2 => vertex(i+1,j)
                element(i,j,k)%vertex3 => vertex(i+1,j+1)
                element(i,j,k)%vertex4 => vertex(i,j+1)

                if(n_moment==6)then
                    element(i,j,k)%vertex5 => nodex(i,j)
                    element(i,j,k)%vertex6 => nodey(i+1,j)
                    element(i,j,k)%vertex7 => nodex(i,j+1)
                    element(i,j,k)%vertex8 => nodey(i,j)
                    element(i,j,k)%vertex9 => nodec(i,j)
                endif

            enddo

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
                element(i,j,0)%umodal(kk)= utemp * ai(kk)

            enddo
        enddo
    enddo


    end subroutine init
!--------------------------------------------------~  
real function exact(x,y,t)
    implicit none
    real, intent(in) :: x,y,t
    ! real, dimension(1:4) :: exact
    real :: rb,rb0
    real :: pk
    real :: alpha


    real :: pu,vth,vth2
    real :: fbt


    if(i_case .eq. 0)then
        alpha = 0.01
        pk = 0.5
        exact = 1./sqrt(2*pi) * (1. +alpha*( cos(0.5*x) ) )*exp(- y**2/2. )
    elseif(i_case .eq. 1)then
        alpha = 0.5
        pk = 0.5
        exact = 1./sqrt(2*pi) * (1. +alpha*( cos(0.5*x) ) )*exp(- y**2/2. )
    elseif(i_case .eq. 2)then
        alpha = 0.01
        pk = 0.5
        exact = 2./7./sqrt(2.*pi) *(1+5.*y**2 ) * ( 1+ alpha*( ( cos(2.*pk*x ) + cos(3.*pk*x) )/1.2 +cos(pk*x)  ) )*exp(-y**2/2.)

    elseif(i_case .eq. 3)then
        pu = 0.99
        vth = 0.3
        pk = 2./13.
        vth2 = 2*vth*vth
        exact = 0.5/vth/sqrt(2.*pi)*(exp(-( y -pu)**2./vth2)+exp(-( y +pu)**2./vth2))*(1.+0.05*cos(pk*x ))
    elseif(i_case .eq. 4)then
        pu = 4.5
        vth = 0.5
        pk = 0.3
        vth2 = 2.*vth*vth
        fbt = 9./(10.*sqrt(2.*pi) ) *exp( -y**2/2. ) + 2./(10.*sqrt(2.*pi) ) * exp( -(y-pu)**2/vth2 )
        exact = fbt * (1+0.04*cos(pk*x) )
    endif

    return
    end function exact
!--------------------------------------------------~  
    

