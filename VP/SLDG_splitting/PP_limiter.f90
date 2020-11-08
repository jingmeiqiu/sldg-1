    subroutine PP_limiter( pp_a , n_moment, pp_theta )
    implicit none
    integer,intent(in) :: n_moment
    real,intent(in) :: pp_a(1:n_moment)
    real,intent(out) :: pp_theta
    real :: pmin

    !
    !if(pp_a(1)<0.)then
    !    print *,ie,je,pp_a(1)
    !    pause
    !endif
    
    
    pp_theta = 1.
    if( n_moment == 1 )then

    elseif( n_moment == 2 )then
        pmin = 20171218.
        pmin = min( pmin , pp_a(1)-0.5*pp_a(2) )
        pmin = min( pmin , pp_a(1)+0.5*pp_a(2) )
        !if( pmin == pp_a(1) )then
        !    pp_theta = 1.
        !else
            pp_theta = min(1., abs(  (pp_a(1) -1d-12)/(pmin-pp_a(1)  )  )       )
     
        !endif
    elseif( n_moment == 3 )then
        pmin = 20171218.
        pmin = min( pmin , pp_a(1) - 0.5*pp_a(2) + pp_a(3)/6. )
        pmin = min( pmin , pp_a(1) + 0.5*pp_a(2) + pp_a(3)/6. )
        if( abs( pp_a(3) )< 1d-15 )then

        elseif( abs( pp_a(2)/(2.*pp_a(3) ) ) <0.5 )then
            pmin = min( pmin , ( 12.*pp_a(1)*pp_a(3)-3.*pp_a(2)*pp_a(2)-pp_a(3)*pp_a(3) )/( 12.*pp_a(3) ) )
        else

        endif

        if( pmin == pp_a(1) )then
            pp_theta = 1.
        else
            pp_theta = min(1., abs(  (pp_a(1)-1d-12) /(pmin-pp_a(1)  )  )       )
        endif

    else

    endif

    end subroutine PP_limiter