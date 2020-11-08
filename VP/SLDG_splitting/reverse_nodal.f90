    subroutine reverse_nodal
    implicit none
    real,allocatable :: temp_rev(:,:,:,:),temp_rev1(:,:,:,:)
    !real :: temp_rev(nx,ny,kdg,kdg),temp_rev1(nx,ny,kdg,kdg)
    integer :: lx,ly


    allocate( temp_rev(1:nel_x,1:nel_y,1: n_g,1: n_g) )
    allocate( temp_rev1(1:nel_x,1:nel_y,1: n_g,1: n_g) )

    do i = 1 , nel_x
        do j = 1 ,nel_y
            !
            do lx = 1,n_g
                do ly = 1,n_g
                    temp_rev(i,j,lx,ly) = elem(i,j) % psi( lx, ly )
                enddo
            enddo
        enddo
    enddo



    do i = 1 , nel_x
        do j = 1 ,nel_y
            !
            do lx = 1, n_g
                do ly = 1, n_g
                    !
                    temp_rev1( i,j,lx,ly ) = temp_rev( i,nel_y+1-j,lx, n_g+1 -ly )

                enddo
            enddo
        enddo
    enddo



    do i = 1 , nel_x
        do j = 1 ,nel_y
            !
            do lx = 1, n_g
                do ly = 1, n_g
                    !
                    elem(i,j) % psi( lx, ly )= temp_rev1(i,j,lx,ly)
                enddo
            enddo

        enddo
    enddo


    deallocate( temp_rev )
    deallocate( temp_rev1 )


    end subroutine reverse_nodal