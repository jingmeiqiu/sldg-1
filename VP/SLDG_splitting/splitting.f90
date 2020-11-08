    subroutine splitting
    implicit none
    real :: umod( n_moment, 1 - ighost : nel + ighost  )
    real :: xrg
    real :: ftemp,tnum
    real :: theta

    tnum = tn -0.5*dt
    do je = 1, nel_y
        !call trans_to_split_modal
        do kk= 1 ,n_g
            do ie = 1 ,nel_x
                !************************compute split modal
                do k = 1 , n_moment
                    ftemp = 0.
                    do ig = 1,n_g
                        ftemp = ftemp + elem(ie,je)%psi( ig,kk )*fle( k-1, x_g(ig) )*w_g(ig)
                        !print *, ftemp ,elem(ie,je)%psi( ig,kk ),fle( k-1, x_g(ig) ),w_g(ig),k
                        !pause
                    enddo
                    elem(ie,je)%split_modal( k ) =  ftemp * ai(k)
                enddo!k
                !************************end of computing split modal
                call PP_limiter(  elem(ie,je)%split_modal( 1:n_moment ) , n_moment, theta )

                elem(ie,je)%split_modal( 2:n_moment ) = theta * elem(ie,je)%split_modal( 2:n_moment )
            enddo
            call bc_x

            do ie = 1 ,nel_x
                ! call 1D SLDG subroutine
                call SLDG1D_QiuGuoCai( 1,dx,tnum, dt/2., glx( 1:n_gl ,1,ie,je),gx( 1:n_g ,kk,ie,je) , gy( 1  ,kk,ie,je),umod( 1:n_moment,ie ) )
                call PP_limiter( umod( 1:n_moment,ie ) , n_moment, theta )

                umod( 2:n_moment ,ie ) = theta * umod( 2:n_moment , ie )
            enddo !ie

            do ie = 1 , nel_x
                elem(ie,je)%split_modal( 1:n_moment ) = umod(1:n_moment,ie)

                !************************compute nodal
                do ig = 1,n_g
                    xrg = x(ie) + dx * x_g(ig)
                    elem(ie,je)%psi(ig,kk) = polynomial( elem(ie,je)%split_modal( 1:n_moment ), xrg, x(ie) ,dx,n_moment-1 )
                enddo
                !************************end of computing nodal
            enddo


        enddo!kk

    enddo! je


    tnum = tn

    call poisson
    call get_norm_DG

    do ie = 1, nel_x

        !call trans_to_split_modal
        do kk= 1 ,n_g
            do je = 1 ,nel_y
                !************************compute split modal
                do k = 1 , n_moment
                    ftemp = 0.
                    do ig = 1,n_g
                        ftemp = ftemp + elem(ie,je)%psi( kk,ig )*fle( k-1, x_g(ig) )*w_g(ig)
                    enddo
                    elem(ie,je)%split_modal( k ) =  ftemp * ai(k)
                enddo!k
                !************************end of computing split modal
            enddo
            call bc_y

            do je = 1 ,nel_y
                ! call 1D SLDG subroutine
                call SLDG1D_QiuGuoCai( 2,dy,tnum,dt, gly( 1 ,1:n_gl,ie,je),gy( kk,1:n_g ,ie,je),gx( kk,1  ,ie,je) ,umod( 1:n_moment,je ) )
                call PP_limiter( umod( 1:n_moment,je ) , n_moment, theta )

                umod( 2:n_moment ,je ) = theta * umod( 2:n_moment , je )
                !if( ie == 12 .and. je == 1 )then
                !    print *,umod(1:n_moment ,je )
                !    pause
                !endif

            enddo !je

            do je = 1 , nel_y
                elem(ie,je)%split_modal( 1:n_moment ) = umod(1:n_moment,je)

                !************************compute nodal
                do ig = 1,n_g
                    xrg = y(je) + dy * x_g(ig)
                    elem(ie,je)%psi(kk,ig) = polynomial( elem(ie,je)%split_modal( 1:n_moment ), xrg, y(je) ,dy,n_moment-1 )
                    !if( ie == 12 .and. je == 1 )then
                    !    print *,elem(ie,je)%split_modal( 1:n_moment )
                    !    print *,elem(ie,je)%psi(kk,ig),kk,ig
                    !    pause
                    !endif

                enddo
                !************************end of computing nodal
            enddo


        enddo!kk

    enddo! ie


    !print *,1111111111
    !***************************************************************************************************
    !***************************************************************************************************
    tnum = tn
    do je = 1, nel_y
        !call trans_to_split_modal
        do kk= 1 ,n_g
            do ie = 1 ,nel_x
                !************************compute split modal
                do k = 1 , n_moment
                    ftemp = 0.
                    do ig = 1,n_g
                        ftemp = ftemp + elem(ie,je)%psi( ig,kk )*fle( k-1, x_g(ig) )*w_g(ig)
                        !print *, ftemp ,elem(ie,je)%psi( ig,kk ),fle( k-1, x_g(ig) ),w_g(ig),k
                        !pause
                    enddo
                    elem(ie,je)%split_modal( k ) =  ftemp * ai(k)

                    !if( elem(ie,je)%split_modal( 1 ) <0. )then
                    !    print *,ie,je,elem(ie,je)%split_modal( 1 )
                    !endif

                enddo!k
                !************************end of computing split modal
            enddo
            call bc_x

            do ie = 1 ,nel_x
                ! call 1D SLDG subroutine
                call SLDG1D_QiuGuoCai( 1,dx,tnum,dt/2., glx( 1:n_gl ,1,ie,je),gx( 1:n_g ,kk,ie,je) ,gy( 1  ,kk,ie,je) ,umod( 1:n_moment,ie ) )
                call PP_limiter( umod( 1:n_moment,ie ) , n_moment, theta )

                umod( 2:n_moment ,ie ) = theta * umod( 2:n_moment , ie )
            enddo !ie

            do ie = 1 , nel_x
                elem(ie,je)%split_modal( 1:n_moment ) = umod(1:n_moment,ie)

                !************************compute nodal
                do ig = 1,n_g
                    xrg = x(ie) + dx * x_g(ig)
                    elem(ie,je)%psi(ig,kk) = polynomial( elem(ie,je)%split_modal( 1:n_moment ), xrg, x(ie) ,dx,n_moment-1 )
                enddo
                !************************end of computing nodal
            enddo
        enddo!kk

    enddo! je

    end subroutine splitting