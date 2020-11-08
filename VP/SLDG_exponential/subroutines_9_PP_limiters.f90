!<1><PP limiter for P1>
subroutine pp_limiter_p1(istage)
    implicit none
    integer,intent(in) :: istage
    real :: pmin
    real :: theta
    integer :: im
    
    real :: ve(1:4,1:2)

    do i = 1 , nx
        do j = 1 ,ny
            pmin= 1000.
            !
            ve(1,1:2) = element(i,j,istage)%vertex1%coor(1:2)
            ve(2,1:2) = element(i,j,istage)%vertex2%coor(1:2)
            ve(3,1:2) = element(i,j,istage)%vertex3%coor(1:2)
            ve(4,1:2) = element(i,j,istage)%vertex4%coor(1:2)
            do im = 1,4
                pmin = min(pmin, &
                    &polynomial( element(i,j,istage)%umodal(1:n_moment),ve(im,1),x(i),dx,ve(im,2) ,y(j),dy,n_moment ) )
            enddo

            if( abs( pmin-element(i,j,istage)%umodal(1) )<eps )then
                theta = 1.
            else
                theta = min(1., abs( element(i,j,istage)%umodal(1) /(pmin-element(i,j,istage)%umodal(1)  )  )       )
            endif

            element(i,j,istage)%umodal(2:3)  = theta* element(i,j,istage)%umodal(2:3)

        enddo
    enddo


    end subroutine pp_limiter_p1
!--------------------------------------------------~  
    
    
!<2><PP limiter for P2>
subroutine pp_limiter_p2(istage)
    implicit none
    integer,intent(in) :: istage
    real :: ve(1:4,1:2)
    real :: pmin
    real :: theta
    integer :: im

    real :: uhat1,uhat2,uhat3,uhat4,uhat5,uhat6
    real :: xmix,ymix
    real :: x111,y111,x222,y222,x333,y333,x444,y444
    ! this subroutine need more comments

    do i = 1 , nx
        do j = 1 ,ny
            pmin= 1000.
            !
            ve(1,1:2) = element(i,j,istage)%vertex1%coor(1:2)
            ve(2,1:2) = element(i,j,istage)%vertex2%coor(1:2)
            ve(3,1:2) = element(i,j,istage)%vertex3%coor(1:2)
            ve(4,1:2) = element(i,j,istage)%vertex4%coor(1:2)            
            do im = 1,4
                pmin = min(pmin,polynomial(element(i,j,istage)%umodal(1:n_moment),ve(im,1),x(i),dx,ve(im,2),y(j),dy,n_moment ) )
            enddo


            !x = -(2*dx*uhat2*uhat6-dx*uhat3*uhat5-4*uhat4*uhat6*x(i)+uhat5**2*x(i))/(4*uhat4*uhat6-uhat5**2),
            !y = (dy*uhat2*uhat5-2*dy*uhat3*uhat4+4*uhat4*uhat6*y(j)-uhat5**2*y(j))/(4*uhat4*uhat6-uhat5**2)
            uhat1 = element(i,j,istage)%umodal(1)
            uhat2 = element(i,j,istage)%umodal(2)
            uhat3 = element(i,j,istage)%umodal(3)
            uhat4 = element(i,j,istage)%umodal(4)
            uhat5 = element(i,j,istage)%umodal(5)
            uhat6 = element(i,j,istage)%umodal(6)

            if( abs(4.*uhat4*uhat6-uhat5**2 - 0.)>eps  )then
                xmix = -(2.*dx*uhat2*uhat6-dx*uhat3*uhat5-4.*uhat4*uhat6*x(i)+uhat5**2*x(i))/(4.*uhat4*uhat6-uhat5**2)
                ymix = (dy*uhat2*uhat5-2.*dy*uhat3*uhat4+4.*uhat4*uhat6*y(j)-uhat5**2*y(j))/(4.*uhat4*uhat6-uhat5**2)
                if( (xmix-xgrid(i) )*(xmix-xgrid(i+1) )<=0. .and.  (ymix-ygrid(j) )*(ymix-ygrid(j+1) )<=0.   )then
                    pmin = min(pmin,polynomial( element(i,j,istage)%umodal(1:n_moment), xmix,x(i),dx,ymix,y(j),dy,n_moment ) )
                endif
            endif

            ! x = x_{-1/2}
            ! y = (1/4)*(-2*dy*uhat3+dy*uhat5+4*uhat6*y(j))/uhat6

            ! x = x_{1/2}
            ! y = -(1/4)*(2*dy*uhat3+dy*uhat5-4*uhat6*y(j))/uhat6

            if( abs(uhat6-0.)>eps  )then
                x111 = xgrid(i)
                y111 = (1./4.)*(-2.*dy*uhat3+dy*uhat5+4.*uhat6*y(j))/uhat6
                if( (y111-ygrid(j) )*(y111-ygrid(j+1) )<=0. )then
                    pmin = min(pmin,polynomial( element(i,j,istage)%umodal(1:n_moment), x111,x(i),dx,y111,y(j),dy,n_moment ) )
                endif

                x222 =  xgrid(i+1)
                y222 = -(1./4.)*(2.*dy*uhat3+dy*uhat5-4.*uhat6*y(j))/uhat6
                if( (y222-ygrid(j) )*(y222-ygrid(j+1) )<=0. )then
                    pmin = min(pmin,polynomial( element(i,j,istage)%umodal(1:n_moment), x222,x(i),dx,y222,y(j),dy,n_moment ) )
                endif

            endif


            ! x = (1/4)*(-2*dx*uhat2+dx*uhat5+4*uhat4*x(i))/uhat4
            ! y = y_{-1/2}

            ! x = -(1/4)*(2*dx*uhat2+dx*uhat5-4*uhat4*x(i))/uhat4
            ! y = y_{1/2}

            if( abs(uhat4-0.)>eps )then
                x333 = (1./4.)*(-2*dx*uhat2+dx*uhat5+4*uhat4*x(i))/uhat4
                y333 = ygrid(j)
                if( (x333-xgrid(i) )*(x333-xgrid(i+1) )<=0. )then
                    pmin = min(pmin,polynomial( element(i,j,istage)%umodal(1:n_moment), x333,x(i),dx,y333,y(j),dy,n_moment ) )
                endif

                x444 = -(1./4.)*(2*dx*uhat2+dx*uhat5-4*uhat4*x(i))/uhat4
                y444 = ygrid(j+1)
                if( ( x444-xgrid(i) )*(x444-xgrid(i+1) )<=0.  )then
                    pmin = min(pmin,polynomial( element(i,j,istage)%umodal(1:n_moment), x444,x(i),dx,y444,y(j),dy,n_moment ) )
                endif

            endif


            if( abs( pmin -element(i,j,istage)%umodal(1) )<eps  )then
                theta = 1.
            else
                 
                theta = min(1., abs(  ( element(i,j,istage)%umodal(1) ) /(pmin-element(i,j,istage)%umodal(1)  )  )       )
            endif

            element(i,j,istage)%umodal(2:6)  = theta* element(i,j,istage)%umodal(2:6)
        enddo
    enddo


    end subroutine pp_limiter_p2
!--------------------------------------------------~  