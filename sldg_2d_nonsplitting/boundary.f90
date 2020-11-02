    subroutine boundary
    implicit none
    integer :: ib

    do i = 1 , nx
        do ib = 1 , nghost
            element(i,1-ib)%umodal(1:n_moment) = element(i,ny+1-ib)%umodal(1:n_moment)
            element(i,ny+ib)%umodal(1:n_moment) = element(i,0+ib)%umodal(1:n_moment)
        enddo
    enddo

    do j = 1 - nghost , ny + nghost
        do ib = 1 , nghost
            element(1-ib,j)%umodal(1:n_moment) = element(nx+1-ib,j)%umodal(1:n_moment)
            element(nx+ib,j)%umodal(1:n_moment) = element(0+ib,j)%umodal(1:n_moment)
        enddo
    enddo


    end subroutine boundary