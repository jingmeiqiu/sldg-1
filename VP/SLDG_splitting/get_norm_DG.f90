    !****************************************************************************************
    subroutine get_norm_DG
    implicit none

    real :: emax,e2  ! get norm
    real :: norm_1   ! L1 norm
    real :: norm_2   ! L2 norm
    real :: energy
    real :: norm_entropy
    real :: enstrophy

    integer :: lx,ly,kk0,kk1
    real :: xrg,yrg

    real :: temptemp


    emax = 0.D0
    e2 = 0.D0
    do i = 1, nel_x
        do lx = 1,6
            xrg =  x(i)+xg(lx)*dx
            e2 = e2 + polynomial_1d(ele_dg( (i-1)*kdg+1:(i)*kdg ) ,xrg,x(i),dx,kdg-1)**2 *wg(lx)

            emax = max(emax, dabs(polynomial_1d(ele_dg( (i-1)*kdg+1:(i)*kdg ) ,xrg,x(i),dx,kdg-1) ) )
        enddo
    enddo
    e2 = sqrt(e2*dx)

    write(3000,*) tn, e2, emax
    write(300,*) tn,log(e2),log(emax)
    write(30,*) tn,-0.1533*tn -3.


    !**********************************************************Mass
    norm_1 = 0.D0
    norm_2 = 0.D0
    do i = 1 , nel_x
        do j = 1 , nel_y
            !
            do lx = 1,n_g
                do ly = 1,n_g
                    norm_1 = norm_1 + elem(i,j)%psi(lx,ly)*w_g(lx)*w_g(ly)
                    norm_2 = norm_2 + elem(i,j)%psi(lx,ly)**2.D0*w_g(lx)*w_g(ly)
                enddo
            enddo

        enddo
    enddo
    norm_1 = norm_1 * dx*dy
    enstrophy = norm_2 * dx*dy
    norm_2 = sqrt( norm_2 * dx*dy )



    energy = 0.D0
    norm_entropy = 0.D0
    do i = 1 , nel_x
        do j = 1, nel_y
            !norm_entropy=norm_entropy+abs(    Dij(i,j)%umodal(1)   )  *log(  abs( Dij(i,j)%umodal(1) )  )
            do lx = 1,n_g
                do ly = 1,n_g
                    yrg =  y(j)+x_g(ly)*dy
                    temptemp = abs( elem(i,j)%psi(lx,ly) )

                    if(temptemp==0.) temptemp = 1d-14
                    norm_entropy = norm_entropy &
                        + elem(i,j)%psi(lx,ly)   *log( temptemp ) *w_g(lx)*w_g(ly)

                    energy = energy &
                        + abs( elem(i,j)%psi(lx,ly) ) &
                        * ( yrg )**2. *w_g(lx)*w_g(ly)*dx*dy
                enddo
            enddo

            !energy = energy + abs( Dij(i,j)%umodal(1) ) * Y(j)**2.D0 * dx * dy
        enddo
        do lx = 1,n_g
            xrg =  x(i)+x_g(lx)*dx
            energy = energy &
                + polynomial_1d(ele_dg( (i-1)*kdg+1:(i)*kdg ) ,xrg,x(i),dx,kdg-1)**2 *w_g(lx) *dx
        enddo
        !energy = energy + Ee(i)**2.D0 * dx
    enddo
    norm_entropy = - norm_entropy * dx * dy

    if(irelative==1)then
        r_entropy =norm_entropy
        r_energy = energy
        r_norm1 = norm_1
        r_enstrophy = enstrophy
        r_norm2 = norm_2
    endif
    irelative = 0

    write(44,*) tn, (norm_entropy-r_entropy)/r_entropy, ( energy - r_energy)/r_energy
    write(555,*) tn, (norm_1-r_norm1)/r_norm1, (norm_2-r_norm2)/r_norm2


    write(4000, *) tn, norm_entropy, energy

    write(77,*) tn,enstrophy
    write(777,*) tn,( enstrophy -r_enstrophy )/enstrophy
    write(5000, *) tn, norm_1, norm_2

    end subroutine get_norm_DG