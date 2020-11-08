    !-------------------------------------------------------------------
    ! semi-Lagrangian Discontinuous Galerkin (SLDG) : 2d advection test
    ! 2D Cartesian plane a Nodal DG formulation
    !                  code by
    !                         Jingmei Qiu, Wei Guo, Xiaofeng Cai
    !                         July 10th, 2017
    !-------------------------------------------------------------------
    Program sldg2d
    use variable
    use element_mod

    implicit none

    
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

    do kkkk = 5,5

        !nel_x = 10*2**(kkkk-1)
        nel_x = 32*kkkk
        
        !nel_x = 200
               
        nel_y = nel_x
        
        norder(kkkk) = nel_x

        nel = max(nel_x,nel_y)

        nod = 2
        n_moment = nod +1
        ! for LDG
        kdg = 3
        n_moment_2d = (nod+1)*(nod + 2)/2
        ! for LDG
        n_gl = 2
        n_g = nod + 1
        i_case = 1
        irk = 3
        tprint = 60
        irelative = 1

        cfl = 20.
        ighost = int(cfl) + 4
        !ighost = int(cfl*4.*pi) + 3

        call allocate_variable
        call LDG_parameters
        ! grid, GL points & Initial (exact) data
        call grid_maker_GL
        call advection_dat
        

        call CPU_TIME(begin_time)
        
        call poisson
        ! strang splitting
        tn = 0.0d0
        nt = 0
        do while(tn<tprint)         !Time Loop

            !call poisson

            call setdt

            call splitting

        enddo ! do while
        !call output
        !pause
        goto 120
        call reverse_nodal

        tprint = tprint * 2.
        do while(tn<tprint)         !Time Loop

            !call poisson

            call setdt

            call splitting

        enddo ! do while

120 continue        
        call CPU_TIME(end_time)
        
        write(2016,*) 'total time',end_time-begin_time
        
        call order_DG
        call output
        !pause
        call deallocate_variable

    enddo

    contains
    include "grid_maker_GL.f90"
    include "allocate_variable.f90"
    include "advection_dat.f90"
    include "setdt.f90"

    include "polynomial.f90"

    include "rk_x.f90"
    include "rk_y.f90"

    include "SLDG1D_QiuGuoCai.f90"
    include "id_get.f90"
    include "green2.f90"
    include "green3.f90"

    include "splitting.f90"
    include "boundary.f90"
    include "order.f90"
    include "output.f90"

    include "SLDG_Poisson.f90"
    include "SLDG_Poisson_para.f90"
    include "poisson.f90"
    include "reverse_nodal.f90"
    !include "init_VP.f90"

    include "PP_limiter.f90"
    
    include "get_norm_DG.f90"

    end program sldg2d