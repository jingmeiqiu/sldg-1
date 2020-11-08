    subroutine setdt
    implicit none
    real :: vx_max,vy_max


    !vx_max = 0.001
    !vy_max = 0.001
    !!dt = cfl /(1./dx + 1./dy )
    !do j = 1, nel_y
    !    do kk = 1 , n_g
    !        do i = 1 , nel_x
    !            do k = 1 , n_g
    !                vx_max = max( vx_max, ax(gx( k ,kk,i,j),gy( k ,kk,i,j),tn)   )
    !                vy_max = max( vy_max, by(gx( k ,kk,i,j),gy( k ,kk,i,j),tn)   )
    !            enddo
    !        enddo
    !    enddo
    !enddo


    !dt = cfl*dx
    dt = cfl /( yright/dx + emax/dy )

    if(tn+dt >tprint )then
        dt = tprint - tn
    endif

    tn = tn + dt

    nt = nt + 1
    if(nt/1*1==nt) print *,tn,tn/tprint*100,"%"

    end subroutine setdt