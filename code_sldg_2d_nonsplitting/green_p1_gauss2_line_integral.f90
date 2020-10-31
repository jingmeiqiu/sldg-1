    subroutine green_p1_gauss2_line_integral(v1,v2,x_base,y_base,idx,idy,s_ab)
    implicit none
    real,intent(in) :: v1(1:2),v2(1:2),x_base,y_base
    integer,intent(in) :: idx,idy
    real,intent(out) :: s_ab(1:n_moment)
    real :: x1,y1,x2,y2
    real :: xshift,yshift
    real :: c1_00,c1_10,c1_01
    real :: xcc,ycc,xcc2,ycc2,x2m1,y2m1
    real :: xgs(2),ygs(2)
    real :: xg1_2,xg2_2,yg1_2,yg2_2
    real :: slope,slope_x21,slope_y21
    real :: green_00,green_10,green_01,green_20,green_11,green_02    
    
    x1 = (v1(1) -x(idx) )/dx ;
    y1 = (v1(2) -y(idy) )/dy ;
    x2 = (v2(1) -x(idx) )/dx ;
    y2 = (v2(2) -y(idy) )/dy ;
    
    xshift = ( x(idx) - x_base ) /dx
    yshift = ( y(idy) - y_base ) /dy
    c1_00 = element(idx,idy)%umodal(1)
    c1_10 = element(idx,idy)%umodal(2)
    c1_01 = element(idx,idy)%umodal(3)
    !
    xcc = ( x1+x2 )*0.5
    ycc = ( y1+y2 )*0.5
    xcc2 = xcc*xcc
    ycc2 = ycc*ycc
    x2m1 = x2-x1
    y2m1 = y2-y1

    !
    xgs(1) = xcc + x2m1 * gau2(1,1)
    xgs(2) = xcc + x2m1 * gau2(2,1)
    ygs(1) = ycc + y2m1 * gau2(1,1)
    ygs(2) = ycc + y2m1 * gau2(2,1)
    xg1_2 = xgs(1)*xgs(1)
    xg2_2 = xgs(2)*xgs(2)
    yg1_2 = ygs(1)*ygs(1)
    yg2_2 = ygs(2)*ygs(2)

    if( abs(x2m1)< 1d-13 .and. abs(y2m1)<1d-13 )then
        s_ab(1:3) = 0.
    else
        if( abs(x2m1) > abs(y2m1) )then
            slope = (y2m1)/(x2m1)
    
            slope_x21 = slope*x2m1
    
            green_00 = xcc* slope_x21
            green_10 = xcc2 *0.5*slope_x21
            green_01 = xcc*ycc*slope_x21
            green_20 = xcc2*xcc*slope_x21/3.
            green_11 = xcc2*ycc*0.5* slope_x21
            green_02 = -ycc2*ycc *x2m1/3.
    
            green_00 = ( xgs(1)+xgs(2) )*0.5* slope_x21
            green_10 = ( xg1_2 + xg2_2 )*0.25 *slope_x21
            green_01 = -( yg1_2 + yg2_2 )*0.25*x2m1
            green_20 = ( xg1_2*xgs(1) + xgs(2)*xg2_2 )/6. *slope_x21
            green_11 = ( xg1_2*ygs(1) + xg2_2*ygs(2) )*0.25 *slope_x21
            green_02 = -( ygs(1)*yg1_2 + ygs(2)*yg2_2 )/6. *x2m1
        
            s_ab(1) = c1_00 * green_00 + c1_10*green_10 + c1_01*green_01

            s_ab(2) = xshift * s_ab(1) + c1_00*green_10 + c1_10*green_20 + c1_01*green_11

            s_ab(3) = yshift * s_ab(1) + c1_00*green_01 + c1_10*green_11 + c1_01*green_02
        else
            slope = (x2m1)/(y2m1)
            slope_y21 = slope* y2m1

            green_00 = xcc* y2m1
            green_10 = xcc2 *0.5* y2m1
            green_01 = xcc*ycc* y2m1
            green_20 = xcc2*xcc * y2m1/3.
            green_11 = xcc2*ycc*0.5* y2m1
            green_02 = -ycc2*ycc* slope_y21/3.
       
            slope_y21 = slope*y2m1
            green_00 = ( xgs(1)+xgs(2) )*0.5* y2m1
            green_10 = ( xg1_2 + xg2_2 )*0.25 *y2m1
            green_01 = -( yg1_2 + yg2_2 )*0.25*slope_y21
            green_20 = ( xg1_2*xgs(1) + xgs(2)*xg2_2  )/6. *y2m1
            green_11 = ( xg1_2*ygs(1) + xg2_2*ygs(2) )*0.25 *y2m1
            green_02 = -( ygs(1)*yg1_2 + ygs(2)*yg2_2 )/6.*slope_y21    
    
            s_ab(1) = c1_00 * green_00 + c1_10*green_10 + c1_01*green_01
            s_ab(2) = xshift * s_ab(1) + c1_00*green_10 + c1_10*green_20 + c1_01*green_11
            s_ab(3) = yshift *s_ab(1) + c1_00*green_01 + c1_10*green_11 + c1_01*green_02
        endif
    endif
    

    end subroutine green_p1_gauss2_line_integral