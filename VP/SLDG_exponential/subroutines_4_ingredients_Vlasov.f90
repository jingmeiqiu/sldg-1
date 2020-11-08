!<1><All ingredients>
subroutine get_rk_stage_velocity_VP(io1, phi, vel)
    implicit none
    integer,intent(in) :: io1
    real,intent(in) :: phi(0:io1)
    real,intent(out) ::  vel  
    integer :: ii
 
    vel =0.
    do ii = 0,io1
        vel = vel + bt_con(io+2,ii+1)*phi(ii)
    enddo
 
    end subroutine get_rk_stage_velocity_VP
!--------------------------------------------------~  
subroutine get_velocity_node_VP(idx, io1, vx, phi)
    implicit none
    integer,intent(in) :: idx,io1
    real,intent(in) :: vx
    real,intent(out) :: phi(0:io1)

    integer :: ii
    real :: eps1
    real :: t1,t2

    eps1 = eps
    if( abs( vx-(x(idx)-0.5*dx) )<=eps1 )then
        do ii = 0,io1
            t1 = polynomial_1d( ele_dg( (idx-2)*kdg+1:(idx-1)*kdg ,ii ) ,vx ,x(idx-1),dx,kdg-1)
            t2 = polynomial_1d( ele_dg( (idx-1)*kdg+1:idx*kdg ,ii ) ,vx ,x(idx),dx,kdg-1)
            phi(ii) = 0.5*t1 + 0.5*t2            
        enddo
    elseif( abs( vx-(x(idx)+0.5*dx) )<=eps1 )then
        do ii = 0,io1
            t1 = polynomial_1d( ele_dg( (idx-1)*kdg+1:(idx)*kdg ,ii ) ,vx ,x(idx),dx,kdg-1)
            t2 = polynomial_1d( ele_dg( (idx)*kdg+1:(idx+1)*kdg ,ii ) ,vx ,x(idx+1),dx,kdg-1)
            phi(ii) = 0.5*t1 + 0.5*t2            
        enddo   
    else
        do ii = 0,io1
            phi(ii) = polynomial_1d( ele_dg( (idx-1)*kdg+1:(idx)*kdg ,ii ) ,vx ,x(idx),dx,kdg-1)
        enddo
    endif

    end subroutine get_velocity_node_VP
!--------------------------------------------------~  
subroutine local_rk4_VP
    implicit none
    real :: vel_x, vel_y
    real :: vx1(1:2),vx2(1:2),vx3(1:2)
    integer :: idx,idy
    real :: t1,t2,t3

    integer :: ii,io1

    !real :: phix(0:3),phiy(0:3)
    !real :: phi(0:3)
    
    real :: phix(0:5),phiy(0:5)
    real :: phi(0:5)
    !**********************************
    ! not module well!

    io1=io
    if(io==3) io1= io-1
    
    if(iexp_case == 8)then
    if(io>=3) io1= 3 ! for cf4
    endif

    do i = 1,nx+1
        do j = 1,ny+1
            ! 
            phi(0:io1) = vertex( i,j )%coor(2)           
            call get_rk_stage_velocity_VP( io1, phi(0:io1) ,vel_x )           
            
            idx = ceiling( (vertex( i,j )%coor(1)-xleft)/dx )
            call get_velocity_node_VP( idx, io1, vertex( i,j )%coor(1), phi(0:io1) )
            call get_rk_stage_velocity_VP( io1, phi(0:io1) ,vel_y )            

            vx1(1) = vertex( i,j )%coor(1) - 0.5* vel_x * dt
            vx1(2) = vertex( i,j )%coor(2) - 0.5* vel_y * dt


            !***********************************************           
            phi(0:io1) = vx1(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx1(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx1(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )

            vx2(1) = vertex( i,j )%coor(1) - 0.5*vel_x*dt  
            vx2(2) = vertex( i,j )%coor(2) - 0.5*vel_y* dt  

            !***********************************************           
            phi(0:io1) = vx2(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx2(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx2(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )

            vx3(1) = vertex( i,j )%coor(1) - vel_x*dt  
            vx3(2) = vertex( i,j )%coor(2) - vel_y* dt  
            
            !********************************************************
            phi(0:io1) = vx3(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx3(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx3(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )
            !
            vertex_star( i,j )%coor(1) = 1./3.*( -vertex( i,j )%coor(1)+vx1(1)+2*vx2(1)+vx3(1) ) &
                &  -1./6.*vel_x*dt
            vertex_star( i,j )%coor(2) = 1./3.*( -vertex( i,j )%coor(2)+vx1(2)+2*vx2(2)+vx3(2) ) &
                &  -1./6.*vel_y*dt

            vertex_star( i,j )%id(1) = ceiling( (vertex_star( i,j )%coor(1)-xleft)/dx )
            vertex_star( i,j )%id(2) = ceiling( (vertex_star( i,j )%coor(2)-ybottom)/dy )
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny+1
            ! 
            phi(0:io1) = nodex( i,j )%coor(2)           
            call get_rk_stage_velocity_VP( io1, phi(0:io1) ,vel_x )           
            
            idx = ceiling( (nodex( i,j )%coor(1)-xleft)/dx )
            call get_velocity_node_VP( idx, io1, nodex( i,j )%coor(1), phi(0:io1) )
            call get_rk_stage_velocity_VP( io1, phi(0:io1) ,vel_y )            

            vx1(1) = nodex( i,j )%coor(1) - 0.5*vel_x * dt
            vx1(2) = nodex( i,j )%coor(2) - 0.5*vel_y * dt


            !***********************************************           
            phi(0:io1) = vx1(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx1(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx1(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )

            vx2(1) = nodex( i,j )%coor(1) - 0.5*vel_x*dt  
            vx2(2) = nodex( i,j )%coor(2) - 0.5*vel_y* dt  

            !***********************************************           
            phi(0:io1) = vx2(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx2(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx2(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )

            vx3(1) = nodex( i,j )%coor(1) - vel_x*dt  
            vx3(2) = nodex( i,j )%coor(2) - vel_y* dt  
            
            !********************************************************
            phi(0:io1) = vx3(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx3(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx3(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )
            !
            nodex_star( i,j )%coor(1) = 1./3.*( -nodex( i,j )%coor(1)+vx1(1)+2*vx2(1)+vx3(1) ) &
                &  -1./6.*vel_x*dt
            nodex_star( i,j )%coor(2) = 1./3.*( -nodex( i,j )%coor(2)+vx1(2)+2*vx2(2)+vx3(2) ) &
                &  -1./6.*vel_y*dt

        enddo
    enddo

    do i = 1,nx+1
        do j = 1,ny
            ! 
            phi(0:io1) = nodey( i,j )%coor(2)           
            call get_rk_stage_velocity_VP( io1, phi(0:io1) ,vel_x )           
            
            idx = ceiling( (nodey( i,j )%coor(1)-xleft)/dx )
            call get_velocity_node_VP( idx, io1, nodey( i,j )%coor(1), phi(0:io1) )
            call get_rk_stage_velocity_VP( io1, phi(0:io1) ,vel_y )            

            vx1(1) = nodey( i,j )%coor(1) - 0.5*vel_x * dt
            vx1(2) = nodey( i,j )%coor(2) - 0.5*vel_y * dt


            !***********************************************           
            phi(0:io1) = vx1(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx1(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx1(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )

            vx2(1) = nodey( i,j )%coor(1) - 0.5*vel_x*dt  
            vx2(2) = nodey( i,j )%coor(2) - 0.5*vel_y* dt  

            !***********************************************           
            phi(0:io1) = vx2(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx2(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx2(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )

            vx3(1) = nodey( i,j )%coor(1) - vel_x*dt  
            vx3(2) = nodey( i,j )%coor(2) - vel_y* dt  
            
            !********************************************************
            phi(0:io1) = vx3(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx3(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx3(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )
            !
            nodey_star( i,j )%coor(1) = 1./3.*( -nodey( i,j )%coor(1)+vx1(1)+2*vx2(1)+vx3(1) ) &
                &  -1./6.*vel_x*dt
            nodey_star( i,j )%coor(2) = 1./3.*( -nodey( i,j )%coor(2)+vx1(2)+2*vx2(2)+vx3(2) ) &
                &  -1./6.*vel_y*dt
        enddo
    enddo

    do i = 1,nx
        do j = 1,ny
            ! 
            phi(0:io1) = nodec( i,j )%coor(2)           
            call get_rk_stage_velocity_VP( io1, phi(0:io1) ,vel_x )           
            
            idx = ceiling( (nodec( i,j )%coor(1)-xleft)/dx )
            call get_velocity_node_VP( idx, io1, nodec( i,j )%coor(1), phi(0:io1) )
            call get_rk_stage_velocity_VP( io1, phi(0:io1) ,vel_y )            

            vx1(1) = nodec( i,j )%coor(1) - 0.5*vel_x * dt
            vx1(2) = nodec( i,j )%coor(2) - 0.5*vel_y * dt


            !***********************************************           
            phi(0:io1) = vx1(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx1(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx1(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )

            vx2(1) = nodec( i,j )%coor(1) - 0.5*vel_x*dt  
            vx2(2) = nodec( i,j )%coor(2) - 0.5*vel_y* dt  

            !***********************************************           
            phi(0:io1) = vx2(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx2(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx2(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )

            vx3(1) = nodec( i,j )%coor(1) - vel_x*dt  
            vx3(2) = nodec( i,j )%coor(2) - vel_y* dt  
            
            !********************************************************
            phi(0:io1) = vx3(2)
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_x )
            
            idx = ceiling( (vx3(1)-xleft)/dx )
            call get_velocity_node_VP( idx,io1,vx3(1),phi(0:io1) )
            call get_rk_stage_velocity_VP( io1,phi(0:io1),vel_y )
            !
            nodec_star( i,j )%coor(1) = 1./3.*( -nodec( i,j )%coor(1)+vx1(1)+2*vx2(1)+vx3(1) ) &
                &  -1./6.*vel_x*dt
            nodec_star( i,j )%coor(2) = 1./3.*( -nodec( i,j )%coor(2)+vx1(2)+2*vx2(2)+vx3(2) ) &
                &  -1./6.*vel_y*dt

        enddo
    enddo

    ! get the vertexes of face_lr
    do j = 1 , ny+1
        do i = 1 , nx
            face_lr(i,j)%point_origin = vertex_star(i,j)
            face_lr(i,j)%point_end = vertex_star(i+1,j)
            ! quadratic-curved
            face_lr(i,j)%point_midt = nodex_star(i,j)
        enddo
    enddo

    ! get the vertexes of face_bt
    do i = 1 , nx+1
        do j = 1 , ny
            face_bt(i,j)%point_origin = vertex_star(i,j)
            face_bt(i,j)%point_end = vertex_star(i,j+1)
            ! quadratic-curved
            face_bt(i,j)%point_midt = nodey_star(i,j)
        enddo
    enddo
    
    end subroutine local_rk4_VP 
!--------------------------------------------------~  
subroutine SLDG_RKEI
    implicit none

    do i = 1 , nx
        do j = 1 , ny
            temp11(i,j,1:n_moment) = element(i,j,io)%umodal(1:n_moment)
        enddo
    enddo

    if(iadditive(io)==0 )then
        call SLDG_Poisson1d( temp11(1:nx,1:ny,1:n_moment) )
    endif

    if(io==0) call setdt

    call local_rk4_VP

    call get_intersections_outersegments
    call get_innersegments

    call get_integral

    end subroutine SLDG_RKEI
!--------------------------------------------------~  
subroutine SLDG_poisson1d( temp16  )
    implicit none

    real,intent(in) :: temp16(nx,ny,n_moment)
    !integer,intent(in) :: kdg

    integer :: l,kk,k
    real :: egamma
    integer :: ly

    real,allocatable :: ff(:)
    real,allocatable :: rho(:,:)

    allocate( ff(kdg*nx) )

    allocate( rho(0:nx,5)  )

    !call LDG_parameters
    !pint = -1.

    rho = pint
    do i = 0, nx-1
        do l = 1, kdg
            do j  = 1 , ny
                do k = 1, n_moment
                    do ly = 1,kdg
                        rho(i,l) = rho(i,l) + temp16(i+1,j,k) * fphi(k,vgau(l,1),vgau(ly,1) ) * wq(ly) * dy
                    enddo
                enddo

            enddo
        enddo
        !
    enddo
 
    ff = 0.
    do i = 0 , nx -1
        do k = 0 , kdg -1
            do l = 1 , kdg
                ff(i*kdg + k +1 ) =  ff(i*kdg + k + 1)  + dx * wq(l) *rho(i,l) * vgau(l,k)
            enddo
        enddo
    enddo

    egamma = 0.
    do i = 0 , nx - 1
        do l = 1,kdg
            egamma = egamma + dx * wq(l) *rho(i,l) * ( x(i+1)+vgau(l,1)*dx )

        enddo
    enddo

    egamma = egamma/(xright-xleft)

    do k = 0, kdg - 1
        
        ff(1+k) = ff(1+k) + egamma*vp(k)

    enddo

    ele_dg(:,io) = 0.
    do k = 0 , kdg -1
        do kk = 1 , kdg
            ele_dg(1+k,io) = ele_dg(1+k,io) + ain(k+1,kk) * ff(kk)

        enddo
    enddo

    do i = 1, nx -1
        do k = 0 , kdg -1
            do kk = 1 , kdg
                ff(i*kdg+ k+1) = ff(i*kdg+1+k) - BL(k+1,kk) * ele_dg( (i-1)*kdg+ kk, io)
            enddo
        enddo
        !
        do k = 0 , kdg -1
            do kk = 1 , kdg
                ele_dg(i*kdg+ k+1 ,io ) = ele_dg( i*kdg+1+k  ,io) + ain(k+1,kk) * ff(i*kdg+ kk)
            enddo
        enddo
    enddo

    EE(0,io) = 0.5* polynomial_1d(ele_dg( (0)*kdg+1:1*kdg ,io ) ,x(1)-0.5*dx,x(1),dx,kdg-1) &
        +0.5*polynomial_1d(ele_dg( (nx-1)*kdg+1:nx*kdg ,io ) ,x(nx)+0.5*dx,x(nx),dx,kdg-1)
    do i = 1,nx-1  ! 1---1/2dx;
        EE(i,io) = 0.5* polynomial_1d(ele_dg( (i-1)*kdg+1:i*kdg ,io ) ,x(i)+0.5*dx,x(i),dx,kdg-1)  &
            + 0.5* polynomial_1d(ele_dg( (i)*kdg+1:(i+1)*kdg ,io ) ,x(i+1)-0.5*dx,x(i+1),dx,kdg-1)

    enddo
    EE(nx,io) = EE(0,io)
 
    do i = 1 ,nx
        ee_c(i,io) = polynomial_1d(ele_dg( (i-1)*kdg+1:i*kdg ,io ) ,x(i),x(i),dx,kdg-1)

    enddo


    !***********************periodic boundary
    ele_dg(-kdg*nghost:0 ,io) = ele_dg( (nx*kdg-kdg*nghost ):nx*kdg, io)
    ele_dg(nx*kdg+1:nx*kdg+kdg*nghost, io ) = ele_dg(1:kdg*nghost ,io )

    deallocate( ff )
    deallocate( rho )

    end subroutine SLDG_poisson1d
!--------------------------------------------------~  
subroutine SLDG_poisson1d_paras
    implicit none
    real :: temp,temp1,temp2


    vm(0) = 1.
    vm(1) = 0.5
    vm(2) = 1./6.
    vm(3) = 1./20.

    vp(0) = 1.
    vp(1) = -0.5
    vp(2) = 1./6.
    vp(3) = -1./20.
    

    !   guassian weight and point for basis
    if(kdg .eq. 1)then
        !
        wq(1) = 1.
        ! v0 4 orders of gaussian qudrature rule( 2 points)
        vxgau(1,0) = 0.  ! l denotes the gaussian point, k denotes the number of basis

        vgau(1,0) = 1.
        
        vgau(1,1) = 0.

        ain(1,1) = 1.
        BL(1,1) = -1.
    elseif(kdg.eq.2) then
        !
        !cfl1 = 1./4.

        wq(1) = 0.5
        wq(2) = 0.5
        ! v0 4 orders of gaussian qudrature rule( 2 points)
        vxgau(1,0) = 0.  ! l denotes the gaussian point, k denotes the number of basis
        vxgau(2,0) = 0.
        vgau(1,0) = 1.
        vgau(2,0) = 1.
        ! v1
        vxgau(1,1) = 1.
        vxgau(2,1) = 1.
        vgau(1,1) = -sqrt(3.)/6.
        vgau(2,1) = -vgau(1,1) 

        ain(1,1) = 0.5; ain(1,2)=-1.; 
        ain(2,1) = 1.; ain(2,2) = 2.; 


        BL(1,1) = -1.; BL(1,2)= -0.5; 
        BL(2,1) = 0.5; BL(2,2) = 0.25 ; 


    else if(kdg.eq.3) then
        ! 
        wq(1) = 5./18.
        wq(2) = 4./9.
        wq(3) = wq(1)
        ! v0 6 orders of gaussian qudrature rule( 3 points)
        vxgau(1,0) = 0.
        vxgau(2,0) = 0.
        vxgau(3,0) = 0.
        vgau(1,0) = 1.
        vgau(2,0) = 1.
        vgau(3,0) = 1.
        ! v1
        temp = sqrt(15.)/10.
        vxgau(1,1) = 1.
        vxgau(2,1) = 1.
        vxgau(3,1) = 1.
        vgau(1,1) = -temp
        vgau(2,1) = 0.
        vgau(3,1) = temp
        ! v2
        vxgau(1,2) = -2.*temp
        vxgau(2,2) = 0.
        vxgau(3,2) = 2.*temp
        vgau(1,2) = temp**2. - 1./12.
        vgau(2,2) = - 1./12.
        vgau(3,2) = vgau(1,2)

        ain(1,1) = 0.5; ain(1,2)=-1.; ain(1,3) = 0.; 
        ain(2,1) = 1.; ain(2,2) = 0.; ain(2,3) = -6.; 
        ain(3,1) = 0.; ain(3,2) = 6.; ain(3,3) = 18.;


        BL(1,1) = -1.; BL(1,2)= -0.5; BL(1,3) = -1./6.;
        BL(2,1) = 0.5; BL(2,2) = 0.25 ; BL(2,3) = 1./12.;
        BL(3,1) = -1./6.; BL(3,2) = -1./12.; BL(3,3) = -1./36.;


    else if(kdg.eq.4) then
        !cfl1 = 1./8.
        wq(1) = 0.25 - sqrt(30.)/72.
        wq(2) = 0.25 + sqrt(30.)/72.
        wq(3) = wq(2)
        wq(4) = wq(1)

        ! v0 8 orders of gaussian qudrature rule( 4 points)
        vxgau(1,0) = 0.
        vxgau(2,0) = 0.
        vxgau(3,0) = 0.
        vxgau(4,0) = 0.
        vgau(1,0) = 1.
        vgau(2,0) = 1.
        vgau(3,0) = 1.
        vgau(4,0) = 1.
        ! v1
        temp1 = sqrt(3./7.+2./35.*sqrt(30.))/2.
        temp2 = sqrt(3./7.-2./35.*sqrt(30.))/2.
        vxgau(1,1) = 1.
        vxgau(2,1) = 1.
        vxgau(3,1) = 1.
        vxgau(4,1) = 1.
        vgau(1,1) = -temp1
        vgau(2,1) = -temp2
        vgau(3,1) = temp2
        vgau(4,1) = temp1
        ! v2
        vxgau(1,2) = 2.*vgau(1,1)
        vxgau(2,2) = 2.*vgau(2,1)
        vxgau(3,2) = 2.*vgau(3,1)
        vxgau(4,2) = 2.*vgau(4,1)
        vgau(1,2) = temp1**2.-1./12.
        vgau(2,2) = temp2**2.-1./12.
        vgau(3,2) = vgau(2,2)
        vgau(4,2) = vgau(1,2)
        ! v3

        vxgau(1,3) = 3.*temp1**2.-3./20.
        vxgau(2,3) = 3.*temp2**2.-3./20.
        vxgau(3,3) = vxgau(2,3)
        vxgau(4,3) = vxgau(1,3)
        vgau(1,3) = -temp1**3.+3./20.*temp1
        vgau(2,3) = -temp2**3.+3./20.*temp2
        vgau(3,3) = -vgau(2,3)
        vgau(4,3) = -vgau(1,3)

        ain(1,1) = 0.5; ain(1,2)=-1.; ain(1,3) = 0.; ain(1,4) = 0.;
        ain(2,1) = 1.; ain(2,2) = 0.; ain(2,3) = -6.; ain(2,4) = 0.;
        ain(3,1) = 0.; ain(3,2) = 6.; ain(3,3) = 0.; ain(3,4) = -60.;
        ain(4,1) = 0.; ain(4,2) = 0.; ain(4,3) = 60.; ain(4,4) = 200.

        BL(1,1) = -1.; BL(1,2)= -0.5; BL(1,3) = -1./6.; BL(1,4) = -1./20.;
        BL(2,1) = 0.5; BL(2,2) = 0.25 ; BL(2,3) = 1./12.; BL(2,4) = 1./40.;
        BL(3,1) = -1./6.; BL(3,2) = -1./12.; BL(3,3) = -1./36.; BL(3,4) = -1./120.;
        BL(4,1) = 1./20.; BL(4,2) = 1./40.; BL(4,3) = 1./120.; BL(4,4) = 1./400.     

    else
        !
        stop

    endif


    end subroutine SLDG_poisson1d_paras
!--------------------------------------------------~  
subroutine VP_after_init
    implicit none
    integer :: lx,ly
    real :: temp_p

    !******************************************************
    if(i_case .eq. 0)then
 
        pint = -1.
    elseif(i_case .eq. 1)then
   
        pint = -1.
    elseif(i_case .eq. 2)then
        pint = 0.
        do i = 1, nx
            do j = 1 , ny
                temp_p = 0.
                do lx = 1 , mmp
                    do ly = 1 ,mmp
                        temp_p = temp_p + dx*wq(lx)*wq(ly)*exact( x(i)+dx*vgau(lx,1),y(j)+dy*vgau(ly,1),0. )*dy
                    enddo
                enddo
                pint = pint + temp_p
            enddo
        enddo
        pint = -pint/(xright-xleft)

    elseif(i_case .eq. 3)then
        pint = 0.
        do i = 1, nx
            do j = 1 , ny
                temp_p = 0.
                do lx = 1 , mmp
                    do ly = 1 ,mmp
                        temp_p = temp_p + dx*wq(lx)*wq(ly)*exact( x(i)+dx*vgau(lx,1),y(j)+dy*vgau(ly,1),0. )*dy
                    enddo
                enddo
                pint = pint + temp_p
            enddo
        enddo
        pint = -pint/(xright-xleft)
    elseif(i_case .eq. 4)then
        pint = 0.
        do i = 1, nx
            do j = 1 , ny
                temp_p = 0.
                do lx = 1 , mmp
                    do ly = 1 ,mmp
                        temp_p = temp_p + dx*wq(lx)*wq(ly)*exact( x(i)+dx*vgau(lx,1),y(j)+dy*vgau(ly,1),0. )*dy
                    enddo
                enddo
                pint = pint + temp_p
            enddo
        enddo
        pint = -pint/(xright-xleft)
    endif

    end subroutine VP_after_init
!--------------------------------------------------~  
subroutine reverse
    implicit none
    real,allocatable :: temp_rev(:,:,:,:),temp_rev1(:,:,:,:)
    !real :: temp_rev(nx,ny,kdg,kdg),temp_rev1(nx,ny,kdg,kdg)
    integer :: lx,ly,k


    allocate( temp_rev(1:nx,1:ny,1:kdg,1:kdg) )
    allocate( temp_rev1(1:nx,1:ny,1:kdg,1:kdg) )

    temp_rev = 0.
    do i = 1 , nx
        do j = 1 ,ny
            !
            do lx = 1,kdg
                do ly = 1,kdg
                    do k = 1,n_moment
                        !
                        temp_rev( i,j,lx,ly ) = temp_rev( i,j,lx,ly ) &
                        +  element(i,j,0)%umodal(k)* fphi(k,vgau(lx,1) ,vgau(ly,1) )
                    enddo
                enddo
            enddo
        enddo
    enddo



    do i = 1 , nx
        do j = 1 ,ny
            !
            do lx = 1,kdg
                do ly = 1,kdg
                    !
                    temp_rev1( i,j,lx,ly ) = temp_rev( i,ny+1-j,lx, kdg +1 -ly ) 

                enddo
            enddo
        enddo
    enddo

    do i = 1 , nx
        do j = 1 ,ny
            !
            do k = 1 , n_moment 
                element(i,j,0)%umodal(k) = 0.
            enddo
        enddo
    enddo

    do i = 1 , nx
        do j = 1 ,ny
            !
            do k = 1 , n_moment         
                do lx = 1,kdg
                    do ly = 1,kdg
                        !
                        element(i,j,0)%umodal(k) = element(i,j,0)%umodal(k)&
                        +  temp_rev1( i,j,lx,ly ) * ai(k) * fphi(k,vgau(lx,1) ,vgau(ly,1) ) *wq(lx)*wq(ly)
                    enddo
                enddo
            enddo
        enddo
    enddo


    deallocate( temp_rev )
    deallocate( temp_rev1 )


    end subroutine reverse