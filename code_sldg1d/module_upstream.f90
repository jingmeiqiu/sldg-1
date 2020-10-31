    module module_upstream

    use module_prob,only : nx
    use module_upstream_data
    implicit none

    contains

    subroutine allocate_var_upstream
    implicit none

    allocate( vertex_star(1:nx+1) )
    allocate( element_star(1:nx) )
    end subroutine allocate_var_upstream

    subroutine deallocate_var_upstream
    implicit none

    deallocate( vertex_star )
    deallocate( element_star )

    end subroutine deallocate_var_upstream
    !*************************************************************************
    subroutine get_upstream_tn(dt)
    use module_dg1d_data, only : vertex
    use module_dg1d_data, only : element1d_Eulerian
    use module_dg1d_data, only : element

    use module_dg1d, only : nx
    use module_dg1d, only : dx

    use module_prob, only : nk
    use module_prob, only : xleft

    use module_prob, only : gau


    implicit none
    real,intent(in) :: dt
    integer :: i,ii
    !******************************************************************
    !
    !******************************************************************
    type(element1d_upstream), pointer :: p
    type(element1d_Eulerian), pointer :: pe

    real :: out_umodal(1:nk+1)
    
    real :: upstream_length

    !***************************
    do i = 1 ,nx+1
        call runge_kutta( vertex(i)%coor ,  vertex_star(i)%coor   ,dt)
        !nk+1 in the
        vertex_star(i)%id = ceiling( (vertex_star(i)%coor-xleft)/dx )
    enddo
    !***************************

    ! assemble gauss lobatto of an upstream element
    ! prepare for interpolating the test function.
    do i = 1,nx
        ! consider an upstream element
        p => element_star(i)
        pe => element(i)

        p%x(1) = vertex_star(i)%coor
        p%x(2) = vertex_star(i+1)%coor
 
        do ii = 1,nk+1
      
            call runge_kutta( pe%xg(ii) , p%xg_star(ii)   ,dt)
        enddo

 
    enddo

    ! get all subinterval of an upstream element
    do i = 1,nx
        ! consider an upstream element
        p => element_star(i)
        p%point_origin = vertex_star(i)
        p%point_end = vertex_star(i+1)

        call search_segment(p)

    enddo
    do i = 1,nx
        ! consider an upstream element
        p => element_star(i)

        call get_integral_pk(p ,  out_umodal(1:nk+1) )

        element_star(i)%umodal(1:nk+1) = out_umodal(1:nk+1)


    enddo



    end subroutine get_upstream_tn
    !*************************************************************************
    subroutine get_matrix( xi,amatrix ,nk)
    implicit none
    real,intent(in) :: xi(1:nk+1)
    integer,intent(in) :: nk
    real,intent(out) :: amatrix(1:nk+1,1:nk+1)

    if( nk==1 )then
        amatrix(1,1) = 1.; amatrix(1,2) = xi(1);
        amatrix(2,1) = 1.; amatrix(2,2) = xi(2);
    elseif( nk==2 )then
        amatrix(1,1)=1.; amatrix(1,2)=xi(1); amatrix(1,3) = xi(1)**2-1./12.;
        amatrix(2,1)=1.; amatrix(2,2)=xi(2); amatrix(2,3) = xi(2)**2-1./12.;
        amatrix(3,1)=1.; amatrix(3,2)=xi(3); amatrix(3,3) = xi(3)**2-1./12.;
    endif

    end subroutine get_matrix
    !*************************************************************************
    subroutine search_segment(p)
    use module_dg1d, only : xgrid

    implicit none
    type(element1d_upstream), pointer :: p
    integer :: mx,kk
    integer :: inter

    mx = p%point_end%id - p%point_origin%id

    p%point_inter(0) = p%point_origin
    p%point_inter(1+mx) = p%point_end
    p%nsub = mx+1

    if(mx .ne. 0)then
        do kk = 1 , mx
            inter = p%point_origin%id + kk
            p%point_inter(kk)%coor = xgrid( inter )
            p%point_inter(kk)%id = inter
        enddo
    endif
    !
    do kk = 1 , 1 + mx
        p%segment(kk)%porigin = p%point_inter(kk-1)
        p%segment(kk)%pend = p%point_inter(kk)
        p%segment(kk)%id = p%point_inter(kk-1)%id
    enddo

    end subroutine search_segment
    !*******************************
    !*************************************************************************
    subroutine get_integral_pk(p,  out_umodal )
 
    use module_prob, only : gau
    use module_prob, only : nk

    use module_dg1d, only : fle
    use module_dg1d_data, only : element
    use module_dg1d, only : x
    use module_dg1d, only : dx
    use module_dg1d, only : ai

    use module_lu, only : solve17,doolittle

    use module_polynomials, only : ortho_poly1d
    implicit none

    type(element1d_upstream), pointer :: p
 
    real,intent(out) :: out_umodal(1:nk+1)
    integer :: kk
    real :: sum
    integer :: idx

    integer :: nm
    real ::  bb(1:nk+1)

    real :: amatrix(1:nk+1,1:nk+1)
    real :: store_L(1:nk+1,1:nk+1),store_U(1:nk+1,1:nk+1)
    integer :: ig

    real :: gau_t(1:nk+1),st
    real :: gau_tt(1:3)


    real :: dx_t

    real :: xs_o,xs_e,xc_star,xi(5)
    integer :: igl
 
    do nm = 1,nk+1
        !*********
        !get test function on element_star(i)

        !********************************
        ! interpolate test function
        if( nm==1 )then
            p%aa(1,1) = 1.
            p%aa(2:nk+1,1) = 0.
        else
            if(nm==2)then

                xc_star = ( p%x(1)+p%x(2) )/2.
 
                
                do ig = 1,nk+1
                    xi(ig)=(p%xg_star(ig)-xc_star)/(p%x(2)-p%x(1))
                enddo

                call get_matrix( xi(1:nk+1) ,amatrix(1:nk+1,1:nk+1),nk)

                call doolittle(Amatrix(1:nk+1,1:nk+1),store_L(1:nk+1,1:nk+1),store_U(1:nk+1,1:nk+1),nk+1)

            endif

            do kk = 1 , nk + 1
                bb(kk) = fle(nm-1, Gau(kk,1) )
            enddo

            call solve17(store_L(1:nk+1,1:nk+1),store_U(1:nk+1,1:nk+1),bb(1:nk+1),p%aa(1:nk+1,nm),nk+1)

        endif
        !********************************


        sum = 0
        do kk = 1 ,p%nsub
            idx = p%segment(kk)%id

            gau_t(1:nk+1) = (p%segment(kk)%pend%coor + p%segment(kk)%porigin%coor)/2. &
                + gau(1:nk+1,1)* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)

            dx_t = p%point_end%coor-p%point_origin%coor
            xc_star = (p%point_origin%coor+p%point_end%coor)/2.
            st =0.
            do ig = 1,nk+1
                st = st + ortho_poly1d(element(idx)%umodal(1:nk+1),gau_t(ig) ,x(idx) ,dx ,nk)&
                    *ortho_poly1d( p%aa(1:nk+1,nm),gau_t(ig),xc_star ,dx_t,nk)*gau(ig,2)
            enddo

            sum = sum + st* (p%segment(kk)%pend%coor - p%segment(kk)%porigin%coor)/dx  !/(p%point_end%coor-p%point_origin%coor)
        enddo

        out_umodal(nm)=sum*ai(nm)
    enddo

    end subroutine get_integral_pk
    !*************************************************************************

    subroutine runge_kutta( vx,vx_star, dt)
    use module_prob, only : ax
    implicit none

    real,intent(in) :: dt
    real,intent(in) :: vx
    real,intent(out) :: vx_star
    real :: vx1,vx2
    real :: vx3


    ! TVD RK3

    !vx1 = vx  - ax( vx  ) * dt
    !vx2 = 0.75*vx  + 0.25*(  vx1  - ax( vx1  ) * dt  )
    !vx_star  = 1./3.*vx  + 2./3.*( vx2  - ax( vx2  ) * dt   )


    !vx_star = vx  - ax( vx  ) * dt

    ! runge-Kutta
    vx1 = vx - 0.5*ax( vx ) * dt
    vx2 = vx - 0.5*ax( vx1 ) * dt
    vx3 = vx - ax( vx2 ) * dt
    vx_star = 1./3. * ( -vx+vx1+2.*vx2+vx3 ) - 1./6.*ax(vx3)*dt

    end subroutine runge_kutta


    subroutine update_solution
    use module_dg1d, only : nk
    use module_dg1d_data, only : element
    use module_dg1d, only : x
    implicit none
    integer :: i

    do i = 1,nx

        element(i)%umodal(1:nk+1) = element_star(i)%umodal(1:nk+1)

    enddo


    end subroutine update_solution

    end module module_upstream