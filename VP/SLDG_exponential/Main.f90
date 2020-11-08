    !*******************************************************************
    ! A well-defined SLDG code couple with the exponential integrators
    ! for the 1D1V Vlasov-Poisson system
    !
    !                             Author:
    !                                  Xiaofeng Cai
    !                                December 15th, 2018
    !*******************************************************************
    program SLDG2D
    use module_data2d
    use module_globals2d
    implicit none
    
    call Parameters

    allocate( ref(1:200,1:200,1:6 ,1:6) )
    do kkkk = 1,5

        call setup
        call allocate_var
        call init
        call VP_after_init

        !--------------------------------------------------~  
        ! positive-preserving limiter
        !if( n_moment ==3 )then
        !    call pp_limiter_p1(0)
        !elseif(n_moment == 6)then
        !    call pp_limiter_p2(0)
        !endif
        !--------------------------------------------------~  
        time = 0.
        nt = 0

        !-------------------- BEGIN TIME EVOLUTION --------------------~  
        do while(time<time_final  )
 
            do io = 0,iexprk-1
                call boundary
                call SLDG_RKEI
                !*******************************************
                ! positive-preserving limiter
                !if( n_moment ==3 )then
                !    call pp_limiter_p1(io+1)
                !elseif(n_moment == 6)then
                !    call pp_limiter_p2(io+1)
                !endif
                !********************************************
            enddo

            !update solution!
            do i = 1 ,nx
                do j = 1,ny
                    element(i,j,0)%umodal(1:n_moment) =  element(i,j,iexprk)%umodal(1:n_moment)
                enddo
            enddo

            nt = nt + 1
            if(nt/1*1==nt) print *,time,time/time_final*100,"%"
        enddo


        call reverse
        
        time_final = time_final*2
 
        do while(time<time_final  )
 
            do io = 0,iexprk-1
                call boundary
                call SLDG_RKEI
                !*******************************************
                ! positive-preserving limiter
                if( n_moment ==3 )then
                    call pp_limiter_p1(io+1)
                elseif(n_moment == 6)then
                    call pp_limiter_p2(io+1)
                endif
                !********************************************
            enddo

            !update solution!
            do i = 1 ,nx
                do j = 1,ny
                    element(i,j,0)%umodal(1:n_moment) =  element(i,j,iexprk)%umodal(1:n_moment)
                enddo
            enddo

            nt = nt + 1
            if(nt/1*1==nt) print *,time,time/time_final*100,"%"
        enddo
 
        call order_DG
  
        call deallocate_var
    enddo

    deallocate( ref )

    contains
    
    !--------------------------------------------------~ 
    include "subroutines_1_setup.f90"
    include "subroutines_2_initialization.f90"
    include "subroutines_3_boundary_condition.f90"
    
    include "subroutines_4_ingredients_Vlasov.f90"
    
    include "subroutines_5_clipping.f90"
    include "subroutines_6_intergrals_on_upstream_cells.f90"
    include "subroutines_7_least_square_for_testfunction.f90"
    include "subroutines_8_output_data.f90"
    include "subroutines_9_PP_limiters.f90"
    
    include "functions_polynomial_calculation.f90"
    !--------------------------------------------------~  
    
    !--------------------------------------------------~ 
    subroutine allocate_var
    implicit none
    integer :: k

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

    allocate( element(1-nghost:nx+nghost,1-nghost:ny+nghost,0:iexprk) )
    allocate( element_star(1:nx,1:ny,0:iexprk) )
    do i = 1-nghost , nx+nghost
        do j = 1-nghost , ny+nghost
            do k = 0,iexprk
                allocate( element(i,j,k)%umodal(1:n_moment) )
            enddo
        enddo
    enddo

    allocate( umod_t(1:nx,1:ny,1:n_moment) )
 
    !*******************************************
    allocate( temp11(1:nx,1:ny,1:n_moment) )
    allocate( ele_dg(-kdg*nghost  :nx*kdg     + kdg*nghost ,0:iexprk ) )
    allocate( ee(0:nx,0:iexprk) )
    allocate( ee_c(1:nx,0:iexprk) )


    end subroutine allocate_var
    !--------------------------------------------------~  
    subroutine deallocate_var
    implicit none
    integer :: k

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
            do k = 0,iexprk
                deallocate( element(i,j,k)%umodal )
            enddo
        enddo
    enddo

    deallocate( element )
    deallocate( element_star )


    deallocate( umod_t )
 
    !*******************************************
    deallocate( temp11 )
    deallocate( ele_dg )
    deallocate( ee )
    deallocate( ee_c )


    end subroutine deallocate_var
    !--------------------------------------------------~  
    subroutine Parameters
    implicit none
        
    pi = 4.*atan(1.)
    eps = 10* EPSILON(0.)
    
    gau2(1,1)=-sqrt(3.0)/6.0
    gau2(2,1)=sqrt(3.0)/6.0
    gau2(1,2)=0.5
    gau2(2,2)=0.5
    
    gau3(1,1) = -sqrt(15.)/10.
    gau3(2,1) = 0.
    gau3(3,1) = sqrt(15.)/10.
    gau3(1,2) = 5./18.
    gau3(2,2) =  4./9.
    gau3(3,2) = 5./18.
    
    ai(1) = 1.
    ai(2) = 12.
    ai(3) = 12.
    ai(4) = 180.
    ai(5) = 144.
    ai(6) = 180.
    !ai(4:6) = 0.
    
    xg(1)=-0.466234757101576013906150777246997304567d0
    xg(2)=-0.330604693233132256830699797509952673503d0
    xg(3)=-0.119309593041598454315250860840355967709d0
    xg(4)=-xg(3)
    xg(5)=-xg(2)
    xg(6)=-xg(1)
    wg(1)=1.71324492379170345040296142172733d-1/2d0
    wg(2)=3.60761573048138607569833513837716d-1/2d0
    wg(3)=4.67913934572691047389870343989551d-1/2d0
    wg(4)=wg(3)
    wg(5)=wg(2)
    wg(6)=wg(1)

    end subroutine Parameters
    !--------------------------------------------------~  
    end program SLDG2D