    subroutine SLDG_Poisson( temp16 )
    implicit none

    !integer,intent(in) :: nx,ny
    real,intent(in) :: temp16(nel_x,nel_y,n_moment_2d)
    !real :: pint
    !real :: rho(0:500,5 )
    integer :: l,kk
    real :: egamma,xrg
    integer :: ly

    real,allocatable :: ff(:)
    real,allocatable :: rho(:,:)

    allocate( ff(kdg*nel_x) )

    allocate( rho(0:nel_x,5)  )

    !call LDG_parameters
    !pint = -1.

    rho = pint
    do i = 0, nel_x-1
        do l = 1, kdg
            do j  = 1 , nel_y
                do k = 1, n_moment_2d
                    do ly = 1,kdg
                        rho(i,l) = rho(i,l) +  temp16(i+1,j,k) *fphi(k,vgau(l,1),vgau(ly,1) ) *wq(ly) *dy
                    enddo
                enddo

            enddo
        enddo
        !
    enddo


    ff = 0.
    do i = 0 , nel_x -1
        do k = 0 , kdg -1
            do l = 1 , kdg
                ff(i*kdg + k +1 ) =  ff(i*kdg + k + 1)  + dx * wq(l) *rho(i,l) * vgau(l,k)
            enddo
        enddo
    enddo

    egamma = 0.
    do i = 0 , nel_x - 1
        do l = 1,kdg
            egamma = egamma + dx * wq(l) *rho(i,l) * ( x(i+1)+vgau(l,1)*dx )
        enddo
    enddo
    !print*,egamma
    egamma = egamma/(xright-xleft)


    do k = 0, kdg - 1
        !
        ff(1+k) = ff(1+k) + egamma*vp(k)
    enddo

    ele_dg = 0.
    do k = 0 , kdg -1
        do kk = 1 , kdg
            ele_dg(1+k) = ele_dg(1+k) + ain(k+1,kk) * ff(kk)
        enddo
    enddo

    do i = 1, nel_x -1
        do k = 0 , kdg -1
            do kk = 1 , kdg
                ff(i*kdg+ k+1) = ff(i*kdg+1+k) - BL(k+1,kk) * ele_dg( (i-1)*kdg+ kk)
            enddo
        enddo
        !
        do k = 0 , kdg -1
            do kk = 1 , kdg
                ele_dg(i*kdg+ k+1) = ele_dg( i*kdg+1+k) + ain(k+1,kk) * ff(i*kdg+ kk)
            enddo
        enddo
    enddo
    
    emax = 0.

    do i = 1 ,nel_x
        do ig = 1 , n_g
            xrg = x(i) + dx*x_g(ig)
            ee_g(ig,i) = polynomial_1d(ele_dg( (i-1)*kdg+1:i*kdg ) ,xrg,x(i),dx,kdg-1)
            emax = max(emax, abs(ee_g(ig,i) ) )
        enddo
    enddo

    deallocate( ff )
    deallocate( rho )

    end subroutine SLDG_Poisson

