    !*************************************************************************
    subroutine order_dg
    use module_prob, only : nx
    use module_prob, only : nk
    use module_prob, only : exact
    
    use module_dg1d_data, only : element
    
    use module_dg1d, only : x
    use module_dg1d, only : dx
    
    use module_polynomials, only : ortho_poly1d
    implicit none
    integer :: kk0,ig
    real :: exact_ave
    real :: xrg
    real :: error1,error2,error3
    real :: rr1,rr2,rr3
    real :: exe
    REAL, PARAMETER, DIMENSION(1:6):: xg=(/ &
        &-0.466234757101576013906150777246997304567e0,&
        &-0.330604693233132256830699797509952673503e0,&
        &-0.119309593041598454315250860840355967709e0,&
        &0.119309593041598454315250860840355967709e0,&
        &0.330604693233132256830699797509952673503e0,&
        &0.466234757101576013906150777246997304567e0 /)
    REAL, PARAMETER, DIMENSION(1:6):: wg=(/ &
        &1.71324492379170345040296142172733e-1/2e0,&
        &3.60761573048138607569833513837716e-1/2e0,&
        &4.67913934572691047389870343989551e-1/2e0,&
        &4.67913934572691047389870343989551e-1/2e0,&
        &3.60761573048138607569833513837716e-1/2e0,&
        &1.71324492379170345040296142172733e-1/2e0 /)

    do kk0 = 1 , nx
        write(1,*) x(kk0),element(kk0)%umodal(1)
    enddo
    close(1)

  
    !open(101,file='error_dg.txt')
    error1=0.0
    error2=0.0
    error3=0.0
    do kk0=1,nx
        exact_ave = 0.
        do ig = 1 , 6
            xrg = x(kk0) + dx * xg(ig)

            exact_ave = exact_ave + exact( xrg,time_final ) * wg(ig)

        enddo

        do ig = 1 , 6
            xrg =x(kk0) + dx * xg(ig)

            exe=exact(xrg,time_final)-ortho_poly1d(element(kk0)%umodal(1:nk+1),xrg,x(kk0),dx,nk)

            error1=error1+abs(exe) *wg(ig)
            error2=error2+(exe)**2 *wg(ig)

            error3=max(error3,abs(exe) )

        enddo

    enddo
    error1=error1/nx
    error2=sqrt(error2/nx)
    if(kkkk.eq.1) write(101,103) nx,error1,error2 ,error3
    write(*,*) 'error',error1,error2
    if(kkkk.gt.1) then
        rr1=log(er1/error1)/log(2.)
        rr2=log(er2/error2)/log(2.)
        rr3=log(er3/error3)/log(2.)
        write(101,102) nx,error1,rr1,error2, rr2 ,error3,rr3
        write(*,*) nx,rr1,rr2,rr3
    endif
    er1=error1
    er2=error2
    er3=error3

102 format(i6,1x,3('&',1x, es12.2e2,1x,'&' 1x,f8.2 ,1x),'\\',1x,'\hline')
103 format(i6,1x,3('&',1x,es12.2E2,1x,'&',1x),'\\',1x,'\hline')
123 format(4(1x,f16.6))




    open(2,file="exact.plt")		
    do kk0= 1 , nx
        exe = 0.
        do ig = 1 , 6
            xrg = x(kk0) + dx * xg(ig)
            exe = exe + exact( xrg , time_final )  * wg(ig)

        enddo
        write(2,123)  x(kk0)  , exe
    enddo
    close(2)


    end subroutine order_dg