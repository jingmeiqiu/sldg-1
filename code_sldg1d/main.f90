    !*************************************************************************
    ! A high order semi-Lagrangian discontinuous Galerkin method for
    ! one-dimensional linear passive transport equation:
    !
    !                  u_t + (a(x)u)_x = 0.
    !
    ! Reference:
    ! [1] Cai, Xiaofeng, Wei Guo, and Jing-Mei Qiu. "A high order conservative
    ! semi-Lagrangian discontinuous Galerkin method for two-dimensional
    ! transport simulations." Journal of Scientific Computing
    ! 73.2-3 (2017): 514-542.
    !*************************************************************************
    ! Versions:
    !
    ! A version by Xiaofeng CAI,University of Delaware,10/29/2020
    !-------------------------------------------------------------------------
    !----this code only works for pk, k=0,1,2;
    !                     (It can be easy to make a general pk, k>2)
    !            you can set up nk= 0,1,2 in setup.f90
    !----two problems:
    !    #1#  linear problem, u_t + u_x =0,                      iexample=1
    !    #2#                  u_t +(sin(x)u)_x =0, by setting up iexample=2
    !*************************************************************************
    program main
    use module_globals ! global variables of the main program

    use module_prob,only : kkkk
    use module_prob,only : setup
    use module_prob,only : time_final

    use module_dg1d,only : allocate_var,deallocate_var
    use module_dg1d,only : init
    use module_upstream,only : allocate_var_upstream,deallocate_var_upstream
    use module_upstream,only : get_upstream_tn
    use module_upstream,only : update_solution
    implicit none


    open(101,file='error_dg.txt')

    do kkkk = 1,7
        ! setup in module_prob
        call setup
        call allocate_var
        call allocate_var_upstream

        ! initial condition
        call init

        time = 0.

        do while( time<time_final )

            call setdt(time_final)

            call boundary

            call get_upstream_tn(dt)

            ! update solution on background elements
            call update_solution
        enddo

        call order_dg

        call deallocate_var
        call deallocate_var_upstream
    enddo

    contains
    include "setdt.f90"
    include "boundary.f90"

    include "order_dg.f90"

    end program main
