    subroutine poisson
    implicit none
    real :: ftemp
    !integer :: ig,jg


    ! compute modal

    do je = 1 , nel_y
        do ie = 1,nel_x
            do k = 1 , n_moment_2d
                ftemp = 0.
                do ig = 1,n_g
                    do jg = 1, n_g
                        ftemp = ftemp + elem(ie,je)%psi(ig,jg)* fphi( k,x_g(ig),x_g(jg) ) * w_g(ig)*w_g(jg)
                    enddo
                enddo
                elem(ie,je)%umodal( k ) = ftemp * ai2d(k)
            enddo
        enddo
    enddo

    do i = 1 , nel_x
        do j = 1 , nel_y
            temp11(i,j,1:n_moment_2d ) = elem(i,j)%umodal(1:n_moment_2d )
        enddo
    enddo

    call SLDG_Poisson( temp11(1:nel_x, 1:nel_y ,1:n_moment_2d) )

    end subroutine poisson