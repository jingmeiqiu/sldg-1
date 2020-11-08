    subroutine grid_maker_GL
    implicit none



    ! compatational domain
    if(i_case .eq. 0)then
        xleft = 0.
        xright = 4.*pi
        yleft = -2.*pi
        yright = 2.*pi
    elseif(i_case .eq. 1)then
        xleft = 0.
        xright = 4.*pi
        yleft = -2.*pi
        yright = 2.*pi
    elseif(i_case .eq. 2)then
        xleft = 0.
        xright = 4.*pi
        yleft = -2.*pi
        yright = 2.*pi
        yleft = -10.
        yright = 10.
    elseif(i_case .eq. 3)then
        xleft = 0.
        xright = 13.*pi
        yleft = -2.*pi
        yright = 2.*pi

        ! yleft = -7.
        ! yright = 7.
    elseif(i_case .eq. 4)then
        xleft = 0.
        xright = 20.*pi/3.
        yleft = -13.
        yright = 13.
    endif

    ! width of an element
    dx = (xright - xleft)/nel_x
    dy = (yright - yleft)/nel_y

    ! x(i) ----> center
    ! y(j) ----> center
    DO I = 1 - ighost, nel_x + ighost
        X(I) = XLEFT + (I-0.5) * DX
    ENDDO

    DO J =  1 -ighost , Nel_y + ighost
        Y(J) = YLEFT + (J-0.5) * DY
    ENDDO

    DO I =   1 -ighost , nel_x + 1 + ighost
        Xgrid(I) = XLEFT + (I-1.) * DX
    ENDDO

    DO J =   1 - ighost , Nel_y+ 1 + ighost
        ygrid(J) = YLEFT + (J-1.) * DY
    ENDDO


    if( n_g .eq. 2 )then
        x_GL( 1 ) = - 0.5
        x_GL( 2 ) = 0.5

        w_GL( 1 ) = 0.5
        w_GL( 2 ) = 0.5


        x_g(1) = -sqrt(3.)/6.
        x_g(2) = -x_g(1)

        w_g(1) = 0.5
        w_g(2) = 0.5

        ai(1) = 1.
        ai(2) = 12.
    elseif( n_g .eq. 3 )then
        x_GL( 1 ) = -0.5
        x_GL( 2 ) = 0.5

        x_g(1) = -sqrt(15.)/10.
        x_g(2) = 0.
        x_g(3) = sqrt(15.)/10.

        w_g(1) = 5./18.
        w_g(2) = 4./9.
        w_g(3) = 5./18.


        ai( 1 ) = 1.
        ai( 2 ) = 12.
        ai( 3 ) = 180.
    elseif(n_g .eq. 4)then
        x_GL( 1 ) = -0.5
        x_GL( 2 ) = -0.5/sqrt(5.)
        x_GL( 3 ) = 0.5/sqrt(5.)
        x_GL( 4 ) = 0.5

        w_GL( 1 ) = 1./12.
        w_GL( 2 ) = 5./12.
        w_GL( 3 ) = 5./12.
        w_GL( 4 ) = 1./12.


        ai( 1 ) = 1.
        ai( 2 ) = 12.
        ai( 3 ) = 180.
        !ai( 4 ) = 2800.
    elseif( n_g .eq. 5 )then
        x_GL( 1 ) = -0.5
        x_GL( 2 ) = -0.5*sqrt( 3./7. )
        x_GL( 3 ) = 0.
        x_GL( 4 ) = -x_GL(2)
        x_GL( 5 ) = -x_GL(1)

        w_GL( 1 ) = 1./20.
        w_GL( 2 ) = 49./180.
        w_GL( 3 ) = 16./45.
        w_GL( 4 ) = w_GL(2)
        w_GL( 5 ) = w_GL(1)


        ai( 1 ) = 1.
        ai( 2 ) = 12.
        ai( 3 ) = 180.
        ai( 4 ) = 2800.
        !ai( 5 ) = 44100.0
    endif

    do i = 1 , nel_x
        do k = 1 , n_gl
            glx(k,:,i,: ) = x_GL( k )*dx + x(i)

        enddo
        do k = 1 , n_g
            gx(k,:,i,: ) = x_G( k )*dx + x(i)

        enddo
    enddo

    do j = 1 , nel_y
        do k = 1 , n_gl
            gly( :,k,:,j ) = x_GL( k )*dy + y(j)
        enddo
        do k = 1 , n_g
            gy( :,k,:,j ) = x_G( k )*dy + y(j)
        enddo
    enddo


    ai2d(1) = 1.
    ai2d(2) = 12.
    ai2d(3) = 12.
    ai2d(4) = 180.
    ai2d(5) = 144.
    ai2d(6) = 180.

    end subroutine grid_maker_GL