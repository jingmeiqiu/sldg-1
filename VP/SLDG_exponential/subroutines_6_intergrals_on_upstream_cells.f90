!<1><a control subroutine>
subroutine get_integral
    implicit none
 
    do i = 1 , nx
        do j = 1 , ny
            element_star(i,j,io)%vertex5 => nodex_star(i,j)
            element_star(i,j,io)%vertex6 => nodey_star(i+1,j)
            element_star(i,j,io)%vertex7 => nodex_star(i,j+1)
            element_star(i,j,io)%vertex8 => nodey_star(i,j)
            element_star(i,j,io)%vertex9 => nodec_star(i,j )
            call green(umod_t(i,j,1:n_moment) )
        enddo
    enddo

    do i = 1 , nx
        do j = 1 , ny
            element(i,j,io+1)%umodal(1:n_moment) = umod_t(i,j,1:n_moment)
        enddo
    enddo

end subroutine get_integral
!--------------------------------------------------~   
    

!<2><get intervals of a certain upstream cell>
subroutine green(sum16)
    implicit none

    real,intent(out) :: sum16(1:n_moment)
    type(type_element_upstream),pointer :: pes
    type(type_element),pointer :: pe
    type(type_face),pointer :: pf
    real :: anticlock
    integer :: nm
    real :: sum
    !*********** least square ***********
    real :: x_base,y_base
    real :: aa(1:n_moment)
    real :: vert_temp(1:9,1:5)
    real :: A_temp(6,6)
    real :: temp_L(6,6),temp_U(6,6)
    real :: psi(1:9),b(1:6)


    pes => element_star(i,j,io)


    if(iadditive(io)==0 )then
        pe => element(i,j,0)
    elseif( iadditive(io)==1 )then
        pe => element(i,j,io)
    endif

    if(io==2 .and. iexp_case==7 )then
        pe => element(i,j,1)
    endif
    
    if(io==2 .and. iexp_case==8 )then
        pe => element(i,j,1)
    endif

    x_base = 0.25*( pes%vertex1%coor(1) + pes%vertex2%coor(1)  &
        + pes%vertex3%coor(1)  + pes%vertex4%coor(1)  )
    y_base = 0.25*( pes%vertex1%coor(2) + pes%vertex2%coor(2)  &
        + pes%vertex3%coor(2)  + pes%vertex4%coor(2)  )
 

    do nm = 1 ,n_moment
        sum = 0.
        !
        if(nm == 1)then
            aa(1) = 1.
            aa(2:n_moment) = 0.
        else
            if( nm ==2 )then

                if(n_moment<=3)then
                    vert_temp(1,1) = ( pes%vertex1%coor(1) - x_base )/dx
                    vert_temp(2,1) = ( pes%vertex2%coor(1) - x_base )/dx
                    vert_temp(3,1) = ( pes%vertex3%coor(1) - x_base )/dx
                    vert_temp(4,1) = ( pes%vertex4%coor(1) - x_base )/dx

                    vert_temp(1,2) = ( pes%vertex1%coor(2) - y_base )/dy
                    vert_temp(2,2) = ( pes%vertex2%coor(2) - y_base )/dy
                    vert_temp(3,2) = ( pes%vertex3%coor(2) - y_base )/dy
                    vert_temp(4,2) = ( pes%vertex4%coor(2) - y_base )/dy
                    call get_matrix_a( vert_temp(1:4,1:2),  A_temp(1:3,1:3) )

                    call doolittle(A_temp(1:3,1:3),temp_L(1:3,1:3),temp_U(1:3,1:3),3)
                elseif(n_moment==6)then
                    vert_temp(1,1) = ( pes%vertex1%coor(1) - x_base )/dx
                    vert_temp(2,1) = ( pes%vertex2%coor(1) - x_base )/dx
                    vert_temp(3,1) = ( pes%vertex3%coor(1) - x_base )/dx
                    vert_temp(4,1) = ( pes%vertex4%coor(1) - x_base )/dx
                    vert_temp(5,1) = ( pes%vertex5%coor(1) - x_base )/dx
                    vert_temp(6,1) = ( pes%vertex6%coor(1) - x_base )/dx
                    vert_temp(7,1) = ( pes%vertex7%coor(1) - x_base )/dx
                    vert_temp(8,1) = ( pes%vertex8%coor(1) - x_base )/dx
                    vert_temp(9,1) = ( pes%vertex9%coor(1) - x_base )/dx

                    vert_temp(1,2) = ( pes%vertex1%coor(2) - y_base )/dy
                    vert_temp(2,2) = ( pes%vertex2%coor(2) - y_base )/dy
                    vert_temp(3,2) = ( pes%vertex3%coor(2) - y_base )/dy
                    vert_temp(4,2) = ( pes%vertex4%coor(2) - y_base )/dy
                    vert_temp(5,2) = ( pes%vertex5%coor(2) - y_base )/dy
                    vert_temp(6,2) = ( pes%vertex6%coor(2) - y_base )/dy
                    vert_temp(7,2) = ( pes%vertex7%coor(2) - y_base )/dy
                    vert_temp(8,2) = ( pes%vertex8%coor(2) - y_base )/dy
                    vert_temp(9,2) = ( pes%vertex9%coor(2) - y_base )/dy

                    vert_temp(1:9,3) = vert_temp(1:9,1)*vert_temp(1:9,1) - 1./12.
                    vert_temp(1:9,4) = vert_temp(1:9,1)*vert_temp(1:9,2)
                    vert_temp(1:9,5) = vert_temp(1:9,2)*vert_temp(1:9,2) - 1./12.

                    call get_matrix2_a( vert_temp(1:9,1:5),  A_temp(1:6,1:6) )
                    call doolittle(A_temp(1:6,1:6),temp_L(1:6,1:6),temp_U(1:6,1:6),6)
                endif
            endif

            if(n_moment<=3)then
                psi(1) = fphi( nm, (pe%vertex1%coor(1)-x(i) )/dx,(pe%vertex1%coor(2)-y(j))/dy ) ;
                psi(2) = fphi( nm, (pe%vertex2%coor(1)-x(i) )/dx,(pe%vertex2%coor(2)-y(j))/dy ) ;
                psi(3) = fphi( nm, (pe%vertex3%coor(1)-x(i) )/dx,(pe%vertex3%coor(2)-y(j))/dy ) ;
                psi(4) = fphi( nm, (pe%vertex4%coor(1)-x(i) )/dx,(pe%vertex4%coor(2)-y(j))/dy ) ;

                call get_vector_b( vert_temp(1:4,1:2),psi(1:4),b(1:3) )
                call solve17(temp_L(1:3,1:3),temp_U(1:3,1:3),b(1:3),aa(1:3),3)
            elseif(n_moment==6)then
                psi(1) = fphi( nm, (pe%vertex1%coor(1)-x(i) )/dx,(pe%vertex1%coor(2)-y(j))/dy ) ;
                psi(2) = fphi( nm, (pe%vertex2%coor(1)-x(i) )/dx,(pe%vertex2%coor(2)-y(j))/dy ) ;
                psi(3) = fphi( nm, (pe%vertex3%coor(1)-x(i) )/dx,(pe%vertex3%coor(2)-y(j))/dy ) ;
                psi(4) = fphi( nm, (pe%vertex4%coor(1)-x(i) )/dx,(pe%vertex4%coor(2)-y(j))/dy ) ;
                psi(5) = fphi( nm, (pe%vertex5%coor(1)-x(i) )/dx,(pe%vertex5%coor(2)-y(j))/dy ) ;
                psi(6) = fphi( nm, (pe%vertex6%coor(1)-x(i) )/dx,(pe%vertex6%coor(2)-y(j))/dy ) ;
                psi(7) = fphi( nm, (pe%vertex7%coor(1)-x(i) )/dx,(pe%vertex7%coor(2)-y(j))/dy ) ;
                psi(8) = fphi( nm, (pe%vertex8%coor(1)-x(i) )/dx,(pe%vertex8%coor(2)-y(j))/dy ) ;
                psi(9) = fphi( nm, (pe%vertex9%coor(1)-x(i) )/dx,(pe%vertex9%coor(2)-y(j))/dy ) ;

                call get_vector2_b( vert_temp(1:9,1:5),psi(1:9),b(1:6) )
                call solve17(temp_L(1:6,1:6),temp_U(1:6,1:6),b(1:6),aa(1:6),6)
            endif


            !if(time>0.8)write(99,*) i,j,nm,aa(1:6)

        endif
        !******************************************************************************************

        if( nm == 1 )then
            pf => pes%edge1_lr
            call green_outersegments_face( pf,x_base,y_base  )
            pf => pes%edge2_bt
            call green_outersegments_face( pf,x_base,y_base  )
            pf => pes%edge3_lr
            call green_outersegments_face( pf,x_base,y_base  )
            pf => pes%edge4_bt
            call green_outersegments_face( pf,x_base,y_base  )

            call green_innersegments(pes,x_base,y_base)
        endif

        pf => pes%edge1_lr
        anticlock = 1.
        call sum_outersegments(nm, pf,anticlock,aa(1:n_moment),sum )
        pf => pes%edge2_bt
        anticlock = 1.
        call sum_outersegments(nm, pf,anticlock,aa(1:n_moment),sum )
        pf => pes%edge3_lr
        anticlock = -1.
        call sum_outersegments(nm, pf,anticlock,aa(1:n_moment),sum )
        pf => pes%edge4_bt
        anticlock = -1.
        call sum_outersegments(nm, pf,anticlock,aa(1:n_moment),sum )


        call sum_innersegments(nm, pes,aa(1:n_moment),sum )


        sum16(nm) = sum * ai(nm)

    enddo

end subroutine green
!--------------------------------------------------~   
subroutine sum_innersegments( nm, pes,aa,sum )
    implicit none
    integer,intent(in) :: nm
    type(type_element_upstream),pointer :: pes
    real,intent(inout) :: sum
    real,intent(in) :: aa(1:n_moment)
    integer :: ii,kk

    integer :: idx,idy,idx1,idy1
 

    do kk = 1 , pes%nsub_inner
        do ii = 1 ,n_moment
            sum = sum + aa(ii)*pes%segment_inner(kk)%c_ab(ii)
        enddo       
    enddo

end subroutine sum_innersegments
!--------------------------------------------------~   
subroutine sum_outersegments( nm,pf,anticlock,aa,sum )
    implicit none
    integer,intent(in) :: nm
    type(type_face),pointer :: pf
    real,intent(inout) :: sum
    real,intent(in) :: aa(1:n_moment)
    real,intent(in) :: anticlock
    integer :: ii,kk

    integer :: idx,idy,idx1,idy1

    do kk = 1 , pf%nsub_outer
        do ii = 1 ,n_moment
            sum = sum + anticlock * aa(ii)*pf%segment_outer(kk)%c_ab(ii)
        enddo
    enddo

    end subroutine sum_outersegments
!--------------------------------------------------~   
subroutine green_outersegments_face(pf,x_base,y_base)
    implicit none
    type(type_face),pointer :: pf
    integer :: k
    real :: v1(1:2),v2(1:2)
    real,intent(in) :: x_base,y_base
    integer :: idx,idy
    real :: s_ab(1:n_moment)
    real :: sol(1:n_moment)
    real :: dl(1:2),de(1:2),centrl(1:2)

    real :: v11(2),v22(2),v33(2)

    ! for qc
    v11(1:2) = pf%point_origin%coor(1:2)
    v22(1:2) = pf%point_midt%coor(1:2)
    v33(1:2) = pf%point_end%coor(1:2)
    ! end for qc

    do k = 1 , pf%nsub_outer
        v1(1:2) = pf%segment_outer(k)%porigin%coor(1:2)
        v2(1:2) = pf%segment_outer(k)%pend%coor(1:2)

        idx = pf%segment_outer(k)%id(1)
        idy = pf%segment_outer(k)%id(2)
 

        if(iadditive(io)==0 )then
            sol(1:n_moment) = element(idx,idy,0)%umodal(1:n_moment)
        elseif( iadditive(io)==1 )then
            sol(1:n_moment) = element(idx,idy,io)%umodal(1:n_moment)
        endif

        if(io==2 .and. iexp_case==7 )then
            sol(1:n_moment) = element(idx,idy,1)%umodal(1:n_moment)
        endif
        
        if(io==2 .and. iexp_case==8 )then
            sol(1:n_moment) = element(idx,idy,1)%umodal(1:n_moment)
        endif        

        if( iqc==0 )then
            dl(1) = dx
            dl(2) = dy
            de(1) = dx
            de(2) = dy
            centrl(1) = x(idx)
            centrl(2) = y(idy)
            call segment_int(n_moment,v1(1:2),v2(1:2),x_base,y_base,idx,idy,centrl,dl,de,sol(1:n_moment),s_ab(1:n_moment))
        elseif( iqc==1 )then
            call green_p2qc_gauss3_outer( v1,v2,x_base,y_base,idx,idy,s_ab,v11,v22,v33,sol(1:n_moment) )
        endif
        pf%segment_outer(k)%c_ab(1:n_moment) = s_ab(1:n_moment)
    enddo

end subroutine green_outersegments_face
!--------------------------------------------------~   
subroutine green_innersegments(pes,x_base,y_base)
    implicit none
    type(type_element_upstream),pointer :: pes
    real :: x_base,y_base
    integer :: k
    integer :: idx,idy
    real :: s_ab(1:n_moment)
    real :: v1(1:2),v2(1:2)
    real :: sol(1:n_moment),dl(1:2),de(1:2),centrl(1:2)

    do k = 1 , pes%nsub_inner
        v1(1:2) = pes%segment_inner(k)%porigin%coor(1:2)
        v2(1:2) = pes%segment_inner(k)%pend%coor(1:2)

        idx = pes%segment_inner(k)%id(1)
        idy = pes%segment_inner(k)%id(2)
 
        if(iadditive(io)==0 )then
            sol(1:n_moment) = element(idx,idy,0)%umodal(1:n_moment)
        elseif( iadditive(io)==1 )then
            sol(1:n_moment) = element(idx,idy,io)%umodal(1:n_moment)
        endif

        if(io==2 .and. iexp_case==7 )then
            sol(1:n_moment) = element(idx,idy,1)%umodal(1:n_moment)
        endif
        
        if(io==2 .and. iexp_case==8 )then
            sol(1:n_moment) = element(idx,idy,1)%umodal(1:n_moment)
        endif        


        dl(1) = dx
        dl(2) = dy
        de(1) = dx
        de(2) = dy
        centrl(1) = x(idx)
        centrl(2) = y(idy)
        call segment_int(n_moment,v1(1:2),v2(1:2),x_base,y_base,idx,idy,centrl,dl,de,sol(1:n_moment),s_ab(1:n_moment))
        pes%segment_inner(k)%c_ab(1:n_moment) = s_ab(1:n_moment)
    enddo

end subroutine green_innersegments
!--------------------------------------------------~   
 
 
    
    
!<3><line interval>
subroutine segment_int(nm,v1,v2,x_base,y_base,idx,idy,centrl,dl,de,sol,s_ab)
    implicit none
    integer,intent(in) :: nm
    real,intent(in) :: v1(1:2),v2(1:2),x_base,y_base,centrl(1:2),dl(1:2),de(1:2)
    integer,intent(in) :: idx,idy
    real,intent(in) :: sol(1:nm)
    real,intent(out) :: s_ab(1:nm)
    real :: x1,y1,x2,y2
    real :: xratio,yratio
    real :: xshift,yshift
    real :: c1_00
    real :: c1_10,c1_01
    real :: c1_20,c1_11,c1_02
    real :: xcc,ycc ,x2m1,y2m1
    real :: xgs(3),ygs(3)
    real :: xg1_2,xg2_2,yg1_2,yg2_2
    real :: xg3_2,yg3_2
    real :: xg1_4,xg2_4,xg3_4
    real :: yg1_4,yg2_4,yg3_4
    real :: xy1,xy2,xy3
    real :: slope,slope_x,slope_y
    real :: green00
    real :: green10,green01,green20,green11,green02
    real :: green30,green21,green12,green03
    real :: green40,green31,green22,green13,green04
    real :: temp2,temp3

    integer :: i_degree

    x1 = (v1(1) -centrl(1) )/dl(1) ;
    y1 = (v1(2) -centrl(2) )/dl(2) ;
    x2 = (v2(1) -centrl(1) )/dl(1) ;
    y2 = (v2(2) -centrl(2) )/dl(2) ;

    xratio = dl(1)/de(1)
    yratio = dl(2)/de(2)

    xshift = ( centrl(1) - x_base ) /de(1)
    yshift = ( centrl(2) - y_base ) /de(2)
    !

    x2m1 = x2-x1
    y2m1 = y2-y1

    if( abs(x2m1)< 1d-13 .and. abs(y2m1)<1d-13 )then
        s_ab(1:nm) = 0.
    else
        if( abs(x2m1) > abs(y2m1) )then
            xcc = ( x1+x2 )*0.5
            ycc = ( y1+y2 )*0.5

            slope = (y2m1)/(x2m1)

            slope_x = slope*x2m1

            green00 = xcc* slope_x

            if(nm>1)then
                xgs(1) = xcc + x2m1 * gau2(1,1)
                xgs(2) = xcc + x2m1 * gau2(2,1)
                ygs(1) = ycc + y2m1 * gau2(1,1)
                ygs(2) = ycc + y2m1 * gau2(2,1)
                xg1_2 = xgs(1)*xgs(1)
                xg2_2 = xgs(2)*xgs(2)
                yg1_2 = ygs(1)*ygs(1)
                yg2_2 = ygs(2)*ygs(2)

                green10 = ( xg1_2 + xg2_2 )*0.25 *slope_x
                green01 = -( yg1_2 + yg2_2 )*0.25*x2m1
                green20 = ( xg1_2*xgs(1) + xgs(2)*xg2_2 )/6. *slope_x
                green11 = ( xg1_2*ygs(1) + xg2_2*ygs(2) )*0.25 *slope_x
                green02 = -( ygs(1)*yg1_2 + ygs(2)*yg2_2 )/6. *x2m1
            endif

            if(nm>3)then

                xgs(1) = xcc + x2m1 * gau3(1,1)
                xgs(2) = xcc + x2m1 * gau3(2,1)
                xgs(3) = xcc + x2m1 * gau3(3,1)
                ygs(1) = ycc + y2m1 * gau3(1,1)
                ygs(2) = ycc + y2m1 * gau3(2,1)
                ygs(3) = ycc + y2m1 * gau3(3,1)

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

                green30=(xg1_4*gau3(1,2)+xg2_4*gau3(2,2)+xg3_4*gau3(3,2))/4.*slope_x
                green21=(xg1_2*xy1*gau3(1,2)+xg2_2*xy2*gau3(2,2)+xg3_2*xy3*gau3(3,2))/3.*slope_x
                green12=-(xy1*yg1_2*gau3(1,2)+xy2*yg2_2*gau3(2,2)+xy3*yg3_2*gau3(3,2))/3.*x2m1
                green03=-(yg1_4*gau3(1,2)+yg2_4*gau3(2,2)+yg3_4*gau3(3,2))/4.*x2m1
                !
                green40=(xg1_4*xgs(1)*gau3(1,2)+xg2_4*xgs(2)*gau3(2,2)+xg3_4*xgs(3)*gau3(3,2))/5.*slope_x
                green31=(xg1_4*ygs(1)*gau3(1,2)+xg2_4*ygs(2)*gau3(2,2)+xg3_4*ygs(3)*gau3(3,2))/4.* slope_x
                green22=(xgs(1)*xg1_2*yg1_2*gau3(1,2)+xgs(2)*xg2_2*yg2_2*gau3(2,2)+xgs(3)*xg3_2*yg3_2*gau3(3,2))/3.*slope_x
                green13=-(yg1_4*xgs(1)*gau3(1,2)+yg2_4*xgs(2)*gau3(2,2)+yg3_4*xgs(3)*gau3(3,2))/4.*x2m1
                green04=-(ygs(1)*yg1_4*gau3(1,2)+ygs(2)*yg2_4*gau3(2,2)+ygs(3)*yg3_4*gau3(3,2))/5.*x2m1
            endif


        else
            xcc = ( x1+x2 )*0.5
            ycc = ( y1+y2 )*0.5

            slope = (x2m1)/(y2m1)
            slope_y = slope* y2m1

            green00 = xcc* y2m1

            if(nm>1)then
                xgs(1) = xcc + x2m1 * gau2(1,1)
                xgs(2) = xcc + x2m1 * gau2(2,1)
                ygs(1) = ycc + y2m1 * gau2(1,1)
                ygs(2) = ycc + y2m1 * gau2(2,1)
                xg1_2 = xgs(1)*xgs(1)
                xg2_2 = xgs(2)*xgs(2)
                yg1_2 = ygs(1)*ygs(1)
                yg2_2 = ygs(2)*ygs(2)

                green10 = ( xg1_2 + xg2_2 )*0.25 *y2m1
                green01 = -( yg1_2 + yg2_2 )*0.25*slope_y
                green20 = ( xg1_2*xgs(1) + xgs(2)*xg2_2  )/6. *y2m1
                green11 = ( xg1_2*ygs(1) + xg2_2*ygs(2) )*0.25 *y2m1
                green02 = -( ygs(1)*yg1_2 + ygs(2)*yg2_2 )/6.*slope_y
            endif

            if(nm>3)then
                xgs(1) = xcc + x2m1 * gau3(1,1)
                xgs(2) = xcc + x2m1 * gau3(2,1)
                xgs(3) = xcc + x2m1 * gau3(3,1)
                ygs(1) = ycc + y2m1 * gau3(1,1)
                ygs(2) = ycc + y2m1 * gau3(2,1)
                ygs(3) = ycc + y2m1 * gau3(3,1)

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

                green30=(xg1_4*gau3(1,2)+xg2_4*gau3(2,2)+xg3_4*gau3(3,2))/4.*y2m1
                green21=(xg1_2*xy1*gau3(1,2)+xg2_2*xy2*gau3(2,2)+xg3_2*xy3*gau3(3,2))/3.*y2m1
                green12=-(xy1*yg1_2*gau3(1,2)+xy2*yg2_2*gau3(2,2)+xy3*yg3_2*gau3(3,2))/3.*slope_y
                green03=-(yg1_4*gau3(1,2)+yg2_4*gau3(2,2)+yg3_4*gau3(3,2))/4.*slope_y
                !
                green40=(xg1_4*xgs(1)*gau3(1,2)+xg2_4*xgs(2)*gau3(2,2)+xg3_4*xgs(3)*gau3(3,2))/5.*y2m1
                green31=(xg1_4*ygs(1)*gau3(1,2)+xg2_4*ygs(2)*gau3(2,2)+xg3_4*ygs(3)*gau3(3,2))/4.*y2m1
                green22=(xgs(1)*xg1_2*yg1_2*gau3(1,2)+xgs(2)*xg2_2*yg2_2*gau3(2,2)+xgs(3)*xg3_2*yg3_2*gau3(3,2))/3.*y2m1
                green13=-(yg1_4*xgs(1)*gau3(1,2)+yg2_4*xgs(2)*gau3(2,2)+yg3_4*xgs(3)*gau3(3,2))/4.*slope_y
                green04=-( ygs(1)*yg1_4*gau3(1,2)+ygs(2)*yg2_4*gau3(2,2)+ygs(3)*yg3_4*gau3(3,2))/5.*slope_y
            endif

        endif

        if(nm ==1 )then
            c1_00 = sol(1)
            s_ab(1) = c1_00 * green00
        elseif(nm == 3)then
            c1_00 = sol(1)
            c1_10 = sol(2)
            c1_01 = sol(3)
            s_ab(1) = c1_00 * green00 + c1_10*green10 + c1_01*green01
            s_ab(2) = xshift *s_ab(1) + xratio*(c1_00*green10 + c1_10*green20 + c1_01*green11)
            s_ab(3) = yshift *s_ab(1) + yratio*(c1_00*green01 + c1_10*green11 + c1_01*green02)
        elseif(nm==6)then
            c1_20 = sol(4)
            c1_11 = sol(5)
            c1_02 = sol(6)
            c1_00 = sol(1) - ( c1_20  +c1_02 )/12.
            c1_10 = sol(2)
            c1_01 = sol(3)

            s_ab(1) =  c1_00 * green00 + c1_10*green10 + c1_01*green01 &
                + c1_20 * green20 + c1_11*green11 + c1_02*green02

            temp2 = c1_00 * green10 + c1_10*green20 + c1_01*green11 &
                + c1_20 * green30 + c1_11*green21 + c1_02*green12
            s_ab(2) =  xshift* s_ab(1)  + xratio*temp2

            temp3 = c1_00 * green01 + c1_10*green11 + c1_01*green02 &
                + c1_20 * green21 + c1_11*green12 + c1_02*green03
            s_ab(3) = yshift* s_ab(1)  + yratio*temp3

            s_ab(4) =  (xshift*xshift-1./12.)* s_ab(1)+ temp2*2.*xshift*xratio &
                +  xratio**2*(c1_00 * green20 + c1_10*green30 + c1_01*green21 &
                + c1_20 * green40 + c1_11*green31 + c1_02*green22)

            s_ab(5) = xshift*yshift*s_ab(1)+ xratio*yshift * temp2 + xshift*yratio * temp3 &
                + xratio*yratio*(c1_00 * green11 + c1_10*green21 + c1_01*green12 &
                + c1_20 * green31 + c1_11*green22 + c1_02*green13)

            s_ab(6) = (yshift*yshift-1./12.)* s_ab(1)+ 2. * yshift*yratio * temp3 &
                + yratio**2*( c1_00 * green02 + c1_10*green12 + c1_01*green03 &
                + c1_20 * green22 + c1_11*green13 + c1_02*green04 )
        endif

    endif

    end subroutine segment_int
!--------------------------------------------------~ 
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
    c1_00 = element(idx,idy,0)%umodal(1)
    c1_10 = element(idx,idy,0)%umodal(2)
    c1_01 = element(idx,idy,0)%umodal(3)
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
!--------------------------------------------------~ 
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

    c1_20 = element(idx,idy,0)%umodal(4)
    c1_11 = element(idx,idy,0)%umodal(5)
    c1_02 = element(idx,idy,0)%umodal(6)
    c1_00 = element(idx,idy,0)%umodal(1) -1./12.*( c1_20  +c1_02 )
    c1_10 = element(idx,idy,0)%umodal(2)
    c1_01 = element(idx,idy,0)%umodal(3)
    
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
!--------------------------------------------------~ 
subroutine green_p2qc_gauss3_outer(v1,v2,x_base,y_base,idx,idy,s_ab,v11,v22,v33,sol)
    implicit none
    real,intent(in) :: v1(1:2),v2(1:2),x_base,y_base,v11(1:2),v22(1:2),v33(1:2)
    integer,intent(in) :: idx,idy
    real,intent(out) :: s_ab(1:n_moment)
    real,intent(in) :: sol(1:n_moment)
    real :: x1,y1,x2,y2
    real :: xshift,yshift

  
    real :: green_00,green_10,green_01,green_20,green_11,green_02,green_30,green_21,green_12,green_03
    real :: green_40,green_31,green_22,green_13,green_04
    real :: c1_00,c1_10,c1_01,c1_20,c1_11,c1_02
 

    real :: temp2,temp3
    !*****
    real :: v_1(1:2),v_2(1:2),v_3(1:2)
    real :: a_trans,b_trans,c_trans,d_trans
    real :: e_x2,e_y2
    real :: ac1,bc1,cc1
    real :: ac2,bc2,cc2
    real :: v(1:2)
    real :: exx1,exx2,dx_xi
    real :: xig(1:3)
    real :: xi_gauss1,xi_gauss2,xi_gauss3
    real :: dxi_gauss1,dxi_gauss2,dxi_gauss3
    real :: yi_gauss1,yi_gauss2,yi_gauss3
    real :: dyi_gauss1,dyi_gauss2,dyi_gauss3

    x1 = (v1(1) -x(idx) )/dx ;
    y1 = (v1(2) -y(idy) )/dy ;
    x2 = (v2(1) -x(idx) )/dx ;
    y2 = (v2(2) -y(idy) )/dy ;
    xshift = ( x(idx) - x_base ) /dx
    yshift = ( y(idy) - y_base ) /dy

    c1_20 = sol(4)
    c1_11 = sol(5)
    c1_02 = sol(6)
    c1_00 = sol(1) -1./12.*( c1_20  +c1_02 )
    c1_10 = sol(2)
    c1_01 = sol(3)

    !**************************************
    v_1(1) = ( v11(1) -x(idx) )/dx ;
    v_1(2) = ( v11(2) -y(idy) )/dy ;
    v_2(1) = ( v22(1) -x(idx) )/dx ;
    v_2(2) = ( v22(2) -y(idy) )/dy ;
    v_3(1) = ( v33(1) -x(idx) )/dx ;
    v_3(2) = ( v33(2) -y(idy) )/dy ;
    call trans( v_1(1:2),v_3(1:2),a_trans,b_trans,c_trans,d_trans )
    e_x2 = trans_to_ex( a_trans ,b_trans , c_trans, v_2(1:2)   )
    e_y2 = trans_to_ey( a_trans ,b_trans , d_trans, v_2(1:2)  )

    ac1 = ( v_3(2)-v_1(2) )/2.*e_y2/( e_x2**2 -1. )
    bc1 = ( v_3(1)-v_1(1) )/2.
    cc1 = ( v_1(1)+v_3(1) )/2. - ( v_3(2)-v_1(2) )/2.*e_y2/( e_x2**2 -1. )

    ac2 = -( v_3(1)-v_1(1) )/2.*e_y2/( e_x2**2 -1. )
    bc2 = ( v_3(2) - v_1(2) )/2.
    cc2 = ( v_1(2)+v_3(2) )/2. + ( v_3(1)-v_1(1) )/2.*e_y2/( e_x2**2 -1. )

    v(1)=x1; v(2) =y1
    exx1 = trans_to_ex( a_trans ,b_trans , c_trans, v(1:2) )
    v(1)=x2; v(2) =y2
    exx2 = trans_to_ex( a_trans ,b_trans , c_trans, v(1:2) )
    dx_xi = exx2 -exx1

    xig(1) = (exx1+exx2)*0.5 + dx_xi * gau3(1,1)
    xig(2) = (exx1+exx2)*0.5 + dx_xi * gau3(2,1)
    xig(3) = (exx1+exx2)*0.5 + dx_xi * gau3(3,1)

    xi_gauss1 = ac1*xig(1)*xig(1) + bc1*xig(1) + cc1
    dxi_gauss1 = 2.*ac1*xig(1) + bc1
    xi_gauss2 = ac1*xig(2)*xig(2) + bc1*xig(2) + cc1
    dxi_gauss2 = 2.*ac1*xig(2) + bc1
    xi_gauss3 = ac1*xig(3)*xig(3) + bc1*xig(3) + cc1
    dxi_gauss3 = 2.*ac1*xig(3) + bc1

    yi_gauss1 = ac2*xig(1)*xig(1) + bc2*xig(1) + cc2
    dyi_gauss1 = 2.*ac2*xig(1) + bc2
    yi_gauss2 = ac2*xig(2)*xig(2) + bc2*xig(2) + cc2
    dyi_gauss2 = 2.*ac2*xig(2) + bc2
    yi_gauss3 = ac2*xig(3)*xig(3) + bc2*xig(3) + cc2
    dyi_gauss3 = 2.*ac2*xig(3) + bc2

    green_00 =  xi_gauss1 * dyi_gauss1*gau3(1,2) + xi_gauss2 * dyi_gauss2*gau3(2,2)+ xi_gauss3* dyi_gauss3*gau3(3,2)
    green_10 = ( xi_gauss1 * xi_gauss1 * dyi_gauss1*gau3(1,2) + xi_gauss2 *xi_gauss2 * dyi_gauss2*gau3(2,2)  &
        + xi_gauss3* xi_gauss3 * dyi_gauss3*gau3(3,2)    )*0.5
    green_01 = -( yi_gauss1*yi_gauss1*dxi_gauss1*gau3(1,2) + yi_gauss2*yi_gauss2*dxi_gauss2*gau3(2,2)        &
        +  yi_gauss3*yi_gauss3*dxi_gauss3*gau3(3,2)  )*0.5
    green_20 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * dyi_gauss1*gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2 * dyi_gauss2*gau3(2,2) &
        + xi_gauss3*xi_gauss3 *xi_gauss3 * dyi_gauss3*gau3(3,2)   )/3.
    green_11 = (  xi_gauss1 * xi_gauss1*yi_gauss1 * dyi_gauss1 *gau3(1,2) &
        + xi_gauss2 *xi_gauss2*yi_gauss2 * dyi_gauss2*gau3(2,2)   &
        + xi_gauss3 *xi_gauss3*yi_gauss3 * dyi_gauss3*gau3(3,2)   )*0.5
    green_02 = -(  yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gau3(1,2) &
        + yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gau3(2,2)  &
        + yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gau3(3,2)    )/3.
    !
    green_30 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * xi_gauss1* dyi_gauss1 *gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2*xi_gauss2 * dyi_gauss2*gau3(2,2) &
        + xi_gauss3*xi_gauss3 *xi_gauss3*xi_gauss3 * dyi_gauss3*gau3(3,2)   )*0.25
    green_21 = ( xi_gauss1*xi_gauss1 * xi_gauss1*yi_gauss1 * dyi_gauss1*gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2*yi_gauss2 * dyi_gauss2*gau3(2,2)   &
        + xi_gauss3*xi_gauss3 *xi_gauss3*yi_gauss3* dyi_gauss3*gau3(3,2)   )/3.
    green_12 = -( xi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gau3(1,2) &
        + xi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gau3(2,2)   &
        + xi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gau3(3,2)    )/3.
    !
    green_03 =  -(  yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gau3(1,2) &
        + yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gau3(2,2)     &
        + yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gau3(3,2)  )*0.25
    !
    green_40 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * xi_gauss1* xi_gauss1* dyi_gauss1*gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2*xi_gauss2*xi_gauss2 * dyi_gauss2*gau3(2,2)  &
        + xi_gauss3*xi_gauss3 *xi_gauss3*xi_gauss3*xi_gauss3 * dyi_gauss3*gau3(3,2)   )*0.2
    green_31 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * xi_gauss1*yi_gauss1* dyi_gauss1*gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2*xi_gauss2*yi_gauss2* dyi_gauss2*gau3(2,2) &
        + xi_gauss3*xi_gauss3 *xi_gauss3*xi_gauss3*yi_gauss3* dyi_gauss3*gau3(3,2)   )*0.25
    green_22 = ( xi_gauss1*xi_gauss1 * xi_gauss1*yi_gauss1*yi_gauss1* dyi_gauss1*gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2*yi_gauss2*yi_gauss2* dyi_gauss2*gau3(2,2)    &
        + xi_gauss3*xi_gauss3 *xi_gauss3*yi_gauss3*yi_gauss3* dyi_gauss3*gau3(3,2)  )/3.
    !
    green_13 = -( yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*xi_gauss1*dxi_gauss1*gau3(1,2) &
        + yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*xi_gauss2*dxi_gauss2*gau3(2,2)  &
        + yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*xi_gauss3*dxi_gauss3*gau3(3,2) )*0.25
    !
    green_04 = -( yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gau3(1,2) &
        + yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gau3(2,2)        &
        + yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gau3(3,2)    )*0.2

    s_ab(1) =  (c1_00 * green_00 + c1_10*green_10 + c1_01*green_01 &
        + c1_20 * green_20 + c1_11*green_11 + c1_02*green_02)*dx_xi
 
    temp2 = (c1_00 * green_10 + c1_10*green_20 + c1_01*green_11 &
        + c1_20 * green_30 + c1_11*green_21 + c1_02*green_12)*dx_xi
    s_ab(2) =  xshift* s_ab(1) + temp2

    temp3 = (c1_00 * green_01 + c1_10*green_11 + c1_01*green_02 &
        + c1_20 * green_21 + c1_11*green_12 + c1_02*green_03 )*dx_xi
    s_ab(3) = yshift* s_ab(1)   +temp3

    s_ab(4) =  (xshift*xshift-1./12. )* s_ab(1) &
        + temp2*2.*xshift &
        +  (c1_00 * green_20 + c1_10*green_30 + c1_01*green_21 &
        + c1_20 * green_40 + c1_11*green_31 + c1_02*green_22)*dx_xi

    s_ab(5) = xshift*yshift*s_ab(1) &
        + yshift * temp2 + xshift * temp3 &
        + (c1_00 * green_11 + c1_10*green_21 + c1_01*green_12 &
        + c1_20 * green_31 + c1_11*green_22 + c1_02*green_13)*dx_xi

    s_ab(6) = (yshift*yshift-1./12.)* s_ab(1) &
        + 2. * yshift * temp3 &
        +  (c1_00 * green_02 + c1_10*green_12 + c1_01*green_03 &
        + c1_20 * green_22 + c1_11*green_13 + c1_02*green_04)*dx_xi
    !**************************************

    end subroutine green_p2qc_gauss3_outer
!--------------------------------------------------~ 

    
    
    
    
    
    
    
    
    
    
    

