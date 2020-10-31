    subroutine green_p2_gauss3_line_integral(v1,v2,x_base,y_base,idx,idy,s_ab)
    implicit none
    real,intent(in) :: v1(1:2),v2(1:2),x_base,y_base
    integer,intent(in) :: idx,idy
    real,intent(out) :: s_ab(1:n_moment)
    real :: x1,y1,x2,y2
    real :: xshift,yshift

    real :: xcc,ycc,xcc2,ycc2,x2m1,y2m1
    real :: xgs(3),ygs(3)
 
    real :: slope 
    real :: green_00,green_10,green_01,green_20,green_11,green_02,green_30,green_21,green_12,green_03
    real :: green_40,green_31,green_22,green_13,green_04
    real :: c1_00,c1_10,c1_01,c1_20,c1_11,c1_02
    real :: slope_x21,slope_y21
    real :: xg1_2,xg2_2,xg3_2,yg1_2,yg2_2,yg3_2
    real :: xy1,xy2,xy3
    real :: xg1_4,xg2_4,xg3_4,yg1_4,yg2_4,yg3_4
    
    real :: temp2,temp3
 
 
    x1 = (v1(1) -x(idx) )/dx ;
    y1 = (v1(2) -y(idy) )/dy ;
    x2 = (v2(1) -x(idx) )/dx ;
    y2 = (v2(2) -y(idy) )/dy ;
    xshift = ( x(idx) - x_base ) /dx
    yshift = ( y(idy) - y_base ) /dy

    c1_20 = element(idx,idy)%umodal(4)
    c1_11 = element(idx,idy)%umodal(5)
    c1_02 = element(idx,idy)%umodal(6)
    c1_00 = element(idx,idy)%umodal(1) -1./12.*( c1_20  +c1_02 )
    c1_10 = element(idx,idy)%umodal(2)
    c1_01 = element(idx,idy)%umodal(3)
    
    xcc = ( x1+x2 )*0.5
    ycc = ( y1+y2 )*0.5


    x2m1 = x2-x1
    y2m1 = y2-y1
    xgs(1) = xcc + x2m1 * gau2(1,1)
    xgs(2) = xcc + x2m1 * gau2(2,1)
    ygs(1) = ycc + y2m1 * gau2(1,1)
    ygs(2) = ycc + y2m1 * gau2(2,1)

    xgs(1) = xcc + x2m1 * gau3(1,1)
    xgs(2) = xcc + x2m1 * gau3(2,1)
    xgs(3) = xcc + x2m1 * gau3(3,1)
    ygs(1) = ycc + y2m1 * gau3(1,1)
    ygs(2) = ycc + y2m1 * gau3(2,1)
    ygs(3) = ycc + y2m1 * gau3(3,1)

    !******************************************
    !xgs(1) = xcc + x2m1 * gau2(1,1)
    !xgs(2) = xcc + x2m1 * gau2(2,1)
    !ygs(1) = ycc + y2m1 * gau2(1,1)
    !ygs(2) = ycc + y2m1 * gau2(2,1)
    !xgs(3) = 0.
    !ygs(3)= 0.
    !gau3(1,2) = 0.5
    !gau3(2,2) =  0.5
    !gau3(3,2) = 0.
    !**************************
    xg1_2 = xgs(1)*xgs(1)
    xg2_2 = xgs(2)*xgs(2)
    xg3_2 = xgs(3)*xgs(3)
    yg1_2 = ygs(1)*ygs(1)
    yg2_2 = ygs(2)*ygs(2)
    yg3_2 = ygs(3)*ygs(3)

    xy1 = xgs(1)*ygs(1)
    xy2 = xgs(2)*ygs(2)
    xy3 = xgs(3)*ygs(3)
    xg1_4 = xg1_2*xg1_2
    xg2_4 = xg2_2 *xg2_2
    xg3_4 = xg3_2 *xg3_2

    yg1_4 = yg1_2 * yg1_2
    yg2_4 = yg2_2 * yg2_2
    yg3_4 = yg3_2 * yg3_2
    !*******************
    ! c2_10 = c1_00
    ! c2_20 = c1_10
    ! c2_11 = c1_01
    !*******************
    ! c3_01 = c1_00
    ! c3_11 = c1_10
    ! c3_02 = c1_01
    !*******************
    if( abs(x2m1)< 1d-12 .and. abs(y2m1)<1d-12 )then
        s_ab(1:6) = 0.
    else
        if( abs(x2m1) > abs(y2m1) )then
            slope = (y2m1)/(x2m1)
            slope_x21 = slope*x2m1
            green_00 = ( xgs(1) *gau3(1,2)+xgs(2) *gau3(2,2) + xgs(3)*gau3(3,2) ) * slope_x21
            green_10 = ( xg1_2 *gau3(1,2) + xg2_2 *gau3(2,2) + xg3_2*gau3(3,2) )*0.5 *slope_x21
            green_01 = -( yg1_2*gau3(1,2) + yg2_2*gau3(2,2) + yg3_2*gau3(3,2) )*0.5*x2m1
            green_20 = ( xg1_2*xgs(1)*gau3(1,2)  + xgs(2)*xg2_2*gau3(2,2)  +xgs(3)*xg3_2*gau3(3,2) )/3. *slope_x21
            green_11 = ( xg1_2*ygs(1)*gau3(1,2)  + xg2_2*ygs(2)*gau3(2,2) +xg3_2*ygs(3)*gau3(3,2)  )*0.5 *slope_x21
            green_02 = -( ygs(1)*yg1_2*gau3(1,2) + ygs(2)*yg2_2*gau3(2,2)  + ygs(3)*yg3_2*gau3(3,2) )/3.*x2m1
            !
            green_30 = ( xg1_4*gau3(1,2) + xg2_4*gau3(2,2)   + xg3_4*gau3(3,2) )*0.25* slope_x21
            green_21 = ( xg1_2*xy1*gau3(1,2) + xg2_2*xy2*gau3(2,2)  + xg3_2*xy3*gau3(3,2) )/3.*slope_x21
            green_12 = -( xy1*yg1_2*gau3(1,2) + xy2*yg2_2*gau3(2,2)  + xy3*yg3_2*gau3(3,2) )/3.*x2m1
            !
            green_03 =  -( yg1_4*gau3(1,2) + yg2_4*gau3(2,2)  +yg3_4*gau3(3,2) )*0.25*x2m1
            !
            green_40 = ( xg1_4*xgs(1)*gau3(1,2) + xg2_4*xgs(2)*gau3(2,2) +xg3_4*xgs(3)*gau3(3,2) )*0.2* slope_x21
            green_31 = ( xg1_4*ygs(1)*gau3(1,2) + xg2_4*ygs(2)*gau3(2,2) + xg3_4*ygs(3)*gau3(3,2) )*0.25* slope_x21
            green_22 = ( xgs(1)*xg1_2*yg1_2*gau3(1,2) + xgs(2)*xg2_2*yg2_2*gau3(2,2)  + xgs(3)*xg3_2*yg3_2*gau3(3,2) )/3.*slope_x21
            !
            green_13 = -( yg1_4*xgs(1)*gau3(1,2) + yg2_4*xgs(2)*gau3(2,2)  + yg3_4*xgs(3)*gau3(3,2) )*0.25*x2m1
            !
            green_04 = -( ygs(1)*yg1_4*gau3(1,2) + ygs(2)*yg2_4*gau3(2,2)  + ygs(3)*yg3_4*gau3(3,2) )*0.2*x2m1
            !*****************
            !green_10 = ( xgs(1)**2 + xgs(2)**2 )*0.25 *slope_x21
            !green_01 = -( ygs(1)**2 + ygs(2)**2 )*0.25*x2m1
            !green_20 = ( xgs(1)**3 + xgs(2)**3 )*0.5*one_third *slope_x21
            !green_11 = ( xgs(1)**2*ygs(1) + xgs(2)**2*ygs(2) )*0.25 *slope_x21
            !green_02 = -( ygs(1)**3 + ygs(2)**3 )*0.5*one_third*x2m1
            !
            !green_30 = ( xgs(1)**4 + xgs(2)**4 )*0.125* slope_x21
            !green_21 = ( xgs(1)**3*ygs(1) + xgs(2)**3*ygs(2) )*0.5*one_third *slope_x21
            !green_12 = -( xgs(1)*ygs(1)**3 + xgs(2)*ygs(2)**3 )*0.5*one_third*x2m1
            !
            !green_03 =  -( ygs(1)**4 + ygs(2)**4 )*0.125*x2m1
            !
            !green_40 = ( xgs(1)**5 + xgs(2)**5 )*0.1* slope_x21
            !green_31 = ( xgs(1)**4*ygs(1) + xgs(2)**4*ygs(2) )*0.125* slope_x21
            !green_22 = ( xgs(1)**3*ygs(1)**2 + xgs(2)**3*ygs(2)**2 )*0.5*one_third *slope_x21
            !
            !green_13 = -( ygs(1)**4*xgs(1) + ygs(2)**4*xgs(2) )*0.125*x2m1
            !
            !green_04 = -( ygs(1)**5 + ygs(2)**5 )*0.1*x2m1
 
        else
            slope = (x2m1)/(y2m1)
            slope_y21 = slope*y2m1
            green_00 = ( xgs(1) *gau3(1,2)+xgs(2) *gau3(2,2) + xgs(3)*gau3(3,2) ) * y2m1
            green_10 = ( xg1_2 *gau3(1,2) + xg2_2 *gau3(2,2) + xg3_2*gau3(3,2) )*0.5 *y2m1
            green_01 = -( yg1_2*gau3(1,2) + yg2_2*gau3(2,2) + yg3_2*gau3(3,2) )*0.5*slope_y21
            green_20 = ( xg1_2*xgs(1)*gau3(1,2)  + xgs(2)*xg2_2*gau3(2,2)  +xgs(3)*xg3_2*gau3(3,2) )/3. *y2m1
            green_11 = ( xg1_2*ygs(1)*gau3(1,2)  + xg2_2*ygs(2)*gau3(2,2) +xg3_2*ygs(3)*gau3(3,2)  )*0.5 *y2m1
            green_02 = -( ygs(1)*yg1_2*gau3(1,2) + ygs(2)*yg2_2*gau3(2,2)  + ygs(3)*yg3_2*gau3(3,2) )/3.*slope_y21
            !
            green_30 = ( xg1_4*gau3(1,2) + xg2_4*gau3(2,2)   + xg3_4*gau3(3,2) )*0.25* y2m1
            green_21 = ( xg1_2*xy1*gau3(1,2) + xg2_2*xy2*gau3(2,2)  + xg3_2*xy3*gau3(3,2) )/3.*y2m1
            green_12 = -( xy1*yg1_2*gau3(1,2) + xy2*yg2_2*gau3(2,2)  + xy3*yg3_2*gau3(3,2) )/3.*slope_y21
            !
            green_03 =  -( yg1_4*gau3(1,2) + yg2_4*gau3(2,2)  +yg3_4*gau3(3,2) )*0.25*slope_y21
            !
            green_40 = ( xg1_4*xgs(1)*gau3(1,2) + xg2_4*xgs(2)*gau3(2,2) +xg3_4*xgs(3)*gau3(3,2) )*0.2* y2m1
            green_31 = ( xg1_4*ygs(1)*gau3(1,2) + xg2_4*ygs(2)*gau3(2,2) + xg3_4*ygs(3)*gau3(3,2) )*0.25* y2m1
            green_22 = ( xgs(1)*xg1_2*yg1_2*gau3(1,2) + xgs(2)*xg2_2*yg2_2*gau3(2,2)  + xgs(3)*xg3_2*yg3_2*gau3(3,2) )/3.*y2m1
            !
            green_13 = -( yg1_4*xgs(1)*gau3(1,2) + yg2_4*xgs(2)*gau3(2,2)  + yg3_4*xgs(3)*gau3(3,2) )*0.25*slope_y21
            !
            green_04 = -( ygs(1)*yg1_4*gau3(1,2) + ygs(2)*yg2_4*gau3(2,2)  + ygs(3)*yg3_4*gau3(3,2) )*0.2*slope_y21

        endif
            s_ab(1) =  c1_00 * green_00 + c1_10*green_10 + c1_01*green_01 &
                + c1_20 * green_20 + c1_11*green_11 + c1_02*green_02

            temp2 = c1_00 * green_10 + c1_10*green_20 + c1_01*green_11 &
                + c1_20 * green_30 + c1_11*green_21 + c1_02*green_12
            s_ab(2) =  xshift* s_ab(1) &
                + temp2

            temp3 = c1_00 * green_01 + c1_10*green_11 + c1_01*green_02 &
                + c1_20 * green_21 + c1_11*green_12 + c1_02*green_03
            s_ab(3) = yshift* s_ab(1) &
                +temp3

            s_ab(4) =  (xshift*xshift-1./12.)* s_ab(1) &
                + temp2*2.*xshift &
                +  c1_00 * green_20 + c1_10*green_30 + c1_01*green_21 &
                + c1_20 * green_40 + c1_11*green_31 + c1_02*green_22

            s_ab(5) = xshift*yshift*s_ab(1) &
                + yshift * temp2 + xshift * temp3 &
                + c1_00 * green_11 + c1_10*green_21 + c1_01*green_12 &
                + c1_20 * green_31 + c1_11*green_22 + c1_02*green_13

            s_ab(6) = (yshift*yshift-1./12.)* s_ab(1) &
                + 2. * yshift * temp3 &
                +  c1_00 * green_02 + c1_10*green_12 + c1_01*green_03 &
                + c1_20 * green_22 + c1_11*green_13 + c1_02*green_04
    endif

    end subroutine green_p2_gauss3_line_integral