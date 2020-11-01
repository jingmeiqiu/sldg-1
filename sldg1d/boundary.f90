    !*************************************************************************
    subroutine boundary
    use module_prob, only : nghost
    use module_prob, only : nk
    use module_dg1d_data, only : element
    use module_dg1d, only : nx
    implicit none
    integer :: i

    do i = 1,nghost
        element(1-i)%umodal(1:nk+1) = element(nx+1-i)%umodal(1:nk+1)
        element(nx+i)%umodal(1:nk+1) = element(0+i)%umodal(1:nk+1)
    enddo

    end subroutine boundary
