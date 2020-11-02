    !*******************************************************************
    ! A well-defined SLDG code
    !             coded by 
    !                    Jingmei Qiu, Wei Guo, Xiaofeng Cai
    !*******************************************************************
    program SLDG2D_quadrilateral
    use module_vertex_face_element
    use globals2d
    use LU
    implicit none

    call parameters
    do kkkk = 1,5

        call setup
        call allocate_var 
        call init

        time = 0.
        nt = 0
        !******************** BEGIN TIME EVOLUTION ***************************
        do while(time<time_final  )

            call get_norm
            call setdt
            call boundary
            call RK_to_upstream
 
            call get_intersections_outersegments
            call get_innersegments

            call get_integral
            
            !call output
            !pause
            nt = nt + 1
            if(nt/1*1==nt) print *,time,time/time_final*100,"%"
        enddo
        
        call get_norm
        call order_DG
    !open(121,file='ave.plt')
    !write(121,*)'zone    ','i=',nx,',    j=',ny    
    !do i = 1 , nx
    !    do j = 1 , ny
    !        WRITE(121,*) X(I),Y(J),  element(i,j)%umodal(1)  
    !    enddo
    !enddo
    !
    !close(121)
    !pause  
        call output
        !pause

        call deallocate_var 

    enddo

    contains
    
    include "RK2D.f90"
    
    include "get_norm.f90"

    include "setdt.f90"
    include "setup.f90"
    include "order_DG.f90"

    include "parameters.f90"

    include "init.f90"
    include "boundary.f90"

    include "get_integral.f90"
    include "green.f90"
    include "green_p1_gauss2_line_integral.f90"
    include "green_p2_gauss3_line_integral.f90"
    include "green_gauss_line_integral.f90"

    include "get_matrix_vector.f90"
    include "get_intersections_outersegments.f90"
    include "get_intersections_qc.f90"
    include "get_intersections.f90"
    include "get_outersegments.f90"
    include "get_innersegments.f90"

    include "polynomials.f90"
    
    include "output.f90"
    
    include "get_mono_inters.f90"
    
    include "super2final_inner_x.f90"
    include "super2final_inner_y.f90"
    
    
    include "get_subfaces.f90"
    include "get_mono_inters_QC.f90"
    
    include "green_p2qc_gauss3_outer.f90"
    
    !include "achive_get_intersections.f90"
    !******************************************************
    subroutine RK_to_upstream
    implicit none
    real :: ran1,ran2

    do i = 1,nx+1
        do j = 1,ny+1
            
            !CALL RANDOM_NUMBER(ran1)
            !CALL RANDOM_NUMBER(ran2)
            !!print *,ran1*eps,ran2*eps
            !!pause
            !vertex_star( i,j )%coor(1:2) = vertex( i,j )%coor(1:2)
            !vertex_star( i,j )%coor(1) = vertex( i,j )%coor(1) - 10*ran1*eps*dt
            !vertex_star( i,j )%coor(2) = vertex( i,j )%coor(2) - 10*ran2*eps*dt

            !vertex_star( i,j )%coor(1) = vertex( i,j )%coor(1) - dt
            !vertex_star( i,j )%coor(2) = vertex( i,j )%coor(2) - dt
            call RK2D( vertex( i,j )%coor(1:2),vertex_star( i,j )%coor(1:2),time,dt,-1. )
                       
            vertex_star( i,j )%id(1) = ceiling( (vertex_star( i,j )%coor(1)-xleft)/dx )
            vertex_star( i,j )%id(2) = ceiling( (vertex_star( i,j )%coor(2)-ybottom)/dy )
        enddo
    enddo
    
    do i = 1,nx 
        do j = 1,ny+1
            call RK2D( nodex( i,j )%coor(1:2),nodex_star( i,j )%coor(1:2),time,dt,-1. )
        enddo
    enddo
    
    do i = 1,nx+1
        do j = 1,ny 
            call RK2D( nodey( i,j )%coor(1:2),nodey_star( i,j )%coor(1:2),time,dt,-1. )
        enddo
    enddo
    
    do i = 1,nx 
        do j = 1,ny 
            call RK2D( nodec( i,j )%coor(1:2),nodec_star( i,j )%coor(1:2),time,dt,-1. )
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


    end subroutine RK_to_upstream

    subroutine allocate_var 
    implicit none

    allocate( xgrid(1-nghost:nx+1+nghost) )
    allocate( ygrid(1-nghost:ny+1+nghost) )

    allocate( x(1-nghost:nx+nghost) )
    allocate( y(1-nghost:ny+nghost) )

    allocate( vertex(  1:(nx+1) , 1:(ny+1)  ) )
    allocate( vertex_star(  1:(nx+1) , 1:(ny+1)  ) )
    
    allocate( nodex( 1:nx,1:ny+1 ) )
    allocate( nodey( 1:nx+1,1:ny ) )
    allocate( nodec( 1:nx,1:ny ) )
    allocate( nodex_star( 1:nx,1:ny+1 ) )
    allocate( nodey_star( 1:nx+1,1:ny ) )
    allocate( nodec_star( 1:nx,1:ny ) )

    allocate( face_lr(1:nx,1:ny+1) )
    allocate( face_bt(1:nx+1,1:ny) )

    allocate( element(1-nghost:nx+nghost,1-nghost:ny+nghost) )
    allocate( element_star(1:nx,1:ny) )
    do i = 1-nghost , nx+nghost
        do j = 1-nghost , ny+nghost
            allocate( element(i,j)%umodal(1:n_moment) )
        enddo
    enddo    

    allocate( umod_t(1:nx,1:ny,1:n_moment) )
    
    allocate( com_mass(1:nx,1:ny) )

    end subroutine allocate_var 
    !******************************************************
    subroutine deallocate_var 
    implicit none

    deallocate( xgrid )
    deallocate( ygrid )

    deallocate( x )
    deallocate( y )

    deallocate( vertex )
    deallocate( vertex_star )
    
    deallocate( nodex )
    deallocate( nodey )
    deallocate( nodec )
    deallocate( nodex_star )
    deallocate( nodey_star )
    deallocate( nodec_star )

    deallocate( face_lr )
    deallocate( face_bt )
    
    do i = 1-nghost , nx+nghost
        do j = 1-nghost , ny+nghost
            deallocate( element(i,j)%umodal )
        enddo
    enddo
    
    deallocate( element )
    deallocate( element_star )


    deallocate( umod_t )
    
    deallocate( com_mass )

    end subroutine deallocate_var 
    !******************************************************
    end program SLDG2D_quadrilateral