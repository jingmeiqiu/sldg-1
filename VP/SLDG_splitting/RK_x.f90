    subroutine RK_x( y_split , vx16,vx_star16,num16,tnum,dt16 )
    implicit none
    integer,intent(in) :: num16 
    real,intent(in) :: tnum,dt16
    real,intent(in) :: vx16(num16)
    real,intent(in) :: y_split
    real :: vx1(num16),vx2(num16),vx3(num16)
    real :: rk_k1(num16),rk_k2(num16),rk_k3(num16),rk_k4(num16),rk_k5(num16),rk_k6(num16)
    real,intent(out) :: vx_star16(num16)
    integer :: ii,ik


    ! TVD RK3
    if(irk ==3)then
        do ik = 1 , num16
            vx1(ik) = vx16(ik) - ax( vx16(ik),y_split ,tnum ) * dt16
            vx2(ik) = 0.75*vx16(ik) + 0.25*(  vx1(ik) - ax( vx1( ik) ,y_split ,tnum -dt16  ) * dt16 )
            vx_star16(ik) = 1./3.*vx16(ik) + 2./3.*( vx2(ik) - ax( vx2(ik) ,y_split ,tnum -dt16*0.5 ) * dt16  )
        enddo
    elseif(irk ==4)then
        do ik = 1 , num16
            vx1(ik) = vx16(ik) - 0.5* ax( vx16(ik),y_split ,tnum ) * dt16
            vx2(ik) = vx16(ik) - 0.5* ax( vx1(ik),y_split ,tnum - dt16*0.5 ) * dt16
            vx3(ik) = vx16(ik) - ax( vx2(ik),y_split  ,tnum-dt16*0.5 ) *dt16
            vx_star16(ik) = 1./3.*( -vx16(ik) + vx1(ik) +2.*vx2(ik)+vx3(ik) ) - 1./6.*dt16* ax( vx3(ik),y_split ,tnum-dt16 )
        enddo

    elseif(irk==5)then
        do ik = 1, num16
            rk_k1(ik) = ax( vx16(ik), y_split ,tnum )

            rk_k2(ik) = ax( vx16(ik)-0.25*rk_k1(ik)*dt16 ,y_split ,tnum-0.25*dt16 )

            rk_k3(ik) = ax( vx16(ik)-rk_k1(ik)*dt16/8.-rk_k2(ik)*dt16/8. ,y_split , tnum-0.25*dt16 )

            rk_k4(ik) = ax(  vx16(ik)+rk_k2(ik)*dt16/2.-rk_k3(ik)*dt16 ,y_split ,tnum-0.5*dt16 )

            rk_k5(ik) = ax(  vx16(ik)-rk_k1(ik)*dt16*3./16.-rk_k4(ik)*dt16*9./16. ,y_split ,tnum-0.75*dt16)

            rk_k6(ik) = ax( vx16(ik)+rk_k1(ik)*dt16*3./7.-rk_k2(ik)*dt16*2./7.-rk_k3(ik)*dt16*12./7.+rk_k4(ik)*dt16*12./7. &
                - rk_k5(ik)*dt16*8./7.  ,y_split ,tnum-dt16 )

            vx_star16(ik) = vx16(ik) - (7.*rk_k1(ik)+32.*rk_k3(ik)+12.*rk_k4(ik)+32.*rk_k5(ik)+7.*rk_k6(ik) )*dt16/90.

        enddo
    endif

    end subroutine RK_x

    ! problem
    real function  ax( x16 , y16 ,t16 )
    implicit none

    real,intent(in) :: x16,y16,t16

    !ax =  -cos(x16/2)**2 * sin(y16) * cos(pi*t16/tprint )*pi
    
    ax = y16

    end function ax
    !*************************************
    !************************************************************************************************************************************

