!<1><get intersections (linear quadrilateral)>
subroutine get_intersections(p)
    implicit none
    real :: vorigin(1:2),vend(1:2)
    integer :: mx,my
    integer :: ia,ib,ic,id
    integer :: kstart

    type(type_face),pointer :: p

    ia = p%point_origin%id(1)
    ib = p%point_end%id(1)

    ic = p%point_origin%id(2)
    id = p%point_end%id(2)
    vorigin(1:2) = p%point_origin%coor(1:2)
    vend(1:2) = p%point_end%coor(1:2)

    mx = abs(ib-ia)

    my = abs(id-ic)

    kstart = 0
    p%nsub_outer = 1+ mx+my
    p%point_inter( 0)%coor( 1 ) = vorigin(1)
    p%point_inter( 0)%coor( 2 ) = vorigin(2)
    p%point_inter( 0)%ixy_type = 0
    p%point_inter( 0)%idx2 = 2*ia
    p%point_inter( 0)%idy2 = 2*ic
    !******
    p%point_inter( 1+ mx+my)%coor( 1 ) = vend(1)
    p%point_inter( 1+ mx+my)%coor( 2 ) = vend(2)
    p%point_inter( 1+ mx+my)%ixy_type = 0
    p%point_inter( 1+ mx+my)%idx2 = 2*ib
    p%point_inter( 1+ mx+my)%idy2 = 2*id
    !**************************************************************************

    call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)


    end subroutine get_intersections
!--------------------------------------------------~
subroutine get_mono_inters(ia,ib,ic,id,mx,my,vorigin,vend,p,kstart)
    implicit none
    integer,intent(in) :: ia,ib,ic,id,mx,my
    real,intent(in) :: vorigin(1:2),vend(1:2)
    integer,intent(in) :: kstart
    integer :: ii,jj
    integer :: k,inter,kk,ki
    integer :: ia2,ib2,ic2,id2
    type(type_face),pointer :: p
    
     
    !**************************************************************************
    ! Algorithm of getting the intersection points between face and gridlines
    ! find intersections between face and x=x_i
    
    if(mx .ne. 0)then

        do k = 1, mx
            if(ia<ib) then
                inter = ia + k
            else
                inter = ia + 1 - k
            endif
            
            kk = kstart + k
            p%point_inter( kk)%coor( 1 ) = xgrid( inter )
            p%point_inter( kk)%coor( 2 ) = yside( vorigin(1:2), vend(1:2),  xgrid(inter) )
            p%point_inter( kk)%ixy_type = 1
            p%point_inter( kk)%igrid = inter
            p%point_inter( kk)%idx2 = 2*inter-1
            p%point_inter( kk)%id_xy = ceiling( (p%point_inter( kk)%coor( 2 )-ybottom)/dy )
            p%point_inter( kk)%idy2 = 2*p%point_inter( kk)%id_xy
        enddo
    endif
    !find intersections between face and y=y_j
    
    if(my .ne. 0)then
 
        do k = 1 + mx, mx+my
            ki = k - mx
            if(ic<id)then
                inter = ic + ki
            else
                inter = ic  + 1 - ki
            endif
            
            kk = kstart + k
            p%point_inter( kk)%coor( 2 ) = ygrid( inter )
            p%point_inter( kk)%coor( 1 ) = xside( vorigin(1:2), vend(1:2),  ygrid(inter) )
            p%point_inter( kk)%ixy_type = 2
            p%point_inter( kk)%igrid = inter
            p%point_inter( kk)%idy2 = 2*inter-1
            p%point_inter( kk)%id_xy = ceiling( (p%point_inter( kk)%coor( 1 )-xleft)/dx )
            p%point_inter( kk)%idx2 = 2*p%point_inter( kk)%id_xy

        enddo
    endif
    !**************************************************************************
    ! Algorithm of sorting point_inter(1:mx+my) based on coordinates
    if(mx+my>1)then
        if(mx>=my)then
            if( vorigin(1)<vend(1) )then
                ! sort them in a increasing order
                do ii = 1 ,mx+my-1
                    do jj = kstart + 1, kstart+mx+my-ii
                        if( p%point_inter(jj)%coor(1)>p%point_inter(jj+1)%coor(1) )then
                            temp = p%point_inter(jj+1)
                            p%point_inter(jj+1) = p%point_inter(jj)
                            p%point_inter(jj) = temp
                        endif
                    enddo
                enddo
            else
                do ii = 1,mx+my-1
                    do jj = kstart+1 , kstart+mx+my-ii
                        if( p%point_inter(jj)%coor(1)<p%point_inter(jj+1)%coor(1) )then
                            temp = p%point_inter(jj+1)
                            p%point_inter(jj+1) = p%point_inter(jj)
                            p%point_inter(jj) = temp
                        endif
                    enddo
                enddo
            endif
        elseif(mx<my)then !mx>=my
            if( vorigin(2)<vend(2) )then
                do ii = 1, mx+my-1
                    do jj = kstart+1,kstart+mx+my-ii
                        if( p%point_inter(jj)%coor(2)>p%point_inter(jj+1)%coor(2) )then
                            temp = p%point_inter(jj+1)
                            p%point_inter(jj+1) = p%point_inter(jj)
                            p%point_inter(jj) = temp
                        endif
                    enddo
                enddo
            else
                do ii = 1,mx+my-1
                    do jj = kstart+1,kstart+mx+my-ii
                        if( p%point_inter(jj)%coor(2)<p%point_inter(jj+1)%coor(2) )then
                            temp = p%point_inter(jj+1)
                            p%point_inter(jj+1) = p%point_inter(jj)
                            p%point_inter(jj) = temp
                        endif
                    enddo
                enddo
            endif
        endif
    endif

    ia2 = 2*ia; ib2 = 2*ib; ic2 = 2*ic; id2 = 2*id;
    !if( mx+my>0 )then
    if( ia2<=ib2 )then
        do kk = kstart+1,kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk-1)%idx2 +1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                elseif( p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 .and. p%point_inter(kk)%idx2>p%point_inter(kk+1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                endif

            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my,kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk+1)%idx2<p%point_inter(kk)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk+1)%idx2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                elseif( p%point_inter(kk+1)%idx2>p%point_inter(kk)%idx2 .and. p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                endif

            endif

        enddo            

    elseif( ia2>ib2 )then
        do kk = kstart+1, kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk-1)%idx2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                elseif( p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 .and. p%point_inter(kk)%idx2<p%point_inter(kk+1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my, kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk+1)%idx2>p%point_inter(kk)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk+1)%idx2 +1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                elseif( p%point_inter(kk+1)%idx2<p%point_inter(kk)%idx2 .and. p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                endif
            endif

        enddo        
    endif

    !************************
    if( ic2<=id2 )then
        do kk = kstart+1, kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk-1)%idy2 +1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                elseif( p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 .and. p%point_inter(kk)%idy2>p%point_inter(kk+1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my,kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk+1)%idy2<p%point_inter(kk)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk+1)%idy2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                elseif( p%point_inter(kk+1)%idy2>p%point_inter(kk)%idy2 .and. p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                endif

            endif

        enddo  
    elseif( ic2>id2 )then
        do kk = kstart+1,kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk-1)%idy2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                elseif( p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 .and. p%point_inter(kk)%idy2<p%point_inter(kk+1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my, kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk+1)%idy2>p%point_inter(kk)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk+1)%idy2 +1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                elseif( p%point_inter(kk+1)%idy2<p%point_inter(kk)%idy2 .and. p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
                endif
            endif

        enddo 
    endif
    !endif


    !*********************************************************************************
    ! check the monotonicity
    ia2 = 2*ia; ib2 = 2*ib; ic2 = 2*ic; id2 = 2*id;

    if( mx+my>0 )then
        if(ia2<=ib2)then
            if( p%point_inter(kstart+1)%idx2<ia2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idx2<p%point_inter(k-1)%idx2 )then
                    write(911,*) 'lose monotone!'
   
                    print *,'lose monotone!3'
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idx2>ib2 )then
                write(911,*) 'lose monotone!'
            
                print *,'lose monotone!4'
                pause
            endif
        elseif( ia2>ib2 )then
            if( p%point_inter(kstart+1)%idx2>ia2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idx2>p%point_inter(k-1)%idx2 )then
                    write(911,*) 'lose monotone!'
 
                    print *,'lose monotone!2'
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idx2<ib2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!'
                pause
            endif
        endif
        !*************************************************
        if(ic2<=id2)then
            if( p%point_inter(kstart+1)%idy2<ic2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idy2<p%point_inter(k-1)%idy2 )then
                    write(911,*) 'lose monotone!'
 
                    print *,'lose monotone!1'
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idy2>id2 )then
                write(911,*) 'lose monotone!'
                pause
            endif
        elseif( ic2>id2 )then
            if( p%point_inter(kstart+1)%idy2>ic2 )then
                write(911,*) 'lose monotone!'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idy2>p%point_inter(k-1)%idy2 )then
                    write(911,*) 'lose monotone!'
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idy2<id2 )then
                write(911,*) 'lose monotone!'
                pause
            endif
        endif
    endif
    !*********************************************************************************

    end subroutine get_mono_inters
!--------------------------------------------------~
real function xside( v1, v2, y)
    implicit none
    real,intent(in) :: v1(2),v2(2)
    real,intent(in) :: y
    real :: r

    r =  ( y- v1(2) )/(v2(2)-v1(2) )

    if( r>1. )then

        xside = v2(1)

    elseif(r<0.)then
        xside = v1(1)
    else
        xside = ( v2(1)-v1(1) )*r + v1(1)
    endif

    end function xside
!--------------------------------------------------~
real function yside( v1, v2, x)
    implicit none
    real,intent(in) :: v1(2),v2(2)
    real,intent(in) :: x
    real :: r

    r = ( x- v1(1) )/(v2(1)-v1(1) )

    if(r>1)then
        yside = v2(2)
    elseif(r<0)then
        yside = v1(2)
    else
        yside = ( v2(2)-v1(2) )*r + v1(2)
    endif

    end function yside
!--------------------------------------------------~

    
!<2><get intersections (quadratic-curved quadrilateral)>  
subroutine get_intersections_qc(p)
    implicit none
    real :: vorigin(1:2),vend(1:2),v11(1:2),v33(1:2)
    integer :: mx,my
    integer :: ia,ib,ic,id
    integer :: kstart,kend
    
    real :: xi_o,xi_e

    type(type_face),pointer :: p
    
    
    v11(1:2) = p%point_origin%coor(1:2)
    v33(1:2) = p%point_end%coor(1:2)

    p%nsub_outer = 0
    kstart = p%nsub_outer
    if( p%nsubface==1 )then
        ia = p%point_origin%id(1)
        ib = p%point_end%id(1)

        ic = p%point_origin%id(2)
        id = p%point_end%id(2)
        vorigin(1:2) = p%point_origin%coor(1:2)
        vend(1:2) = p%point_end%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
        kstart = p%nsub_outer

        p%nsub_outer = 1+ mx+my
        p%point_inter( 0)%coor( 1 ) = vorigin(1)
        p%point_inter( 0)%coor( 2 ) = vorigin(2)
        p%point_inter( 0)%ixy_type = 0
        p%point_inter( 0)%idx2 = 2*ia
        p%point_inter( 0)%idy2 = 2*ic
        !******
        !print *,vorigin(1:2),vend(1:2)
        p%point_inter( 1+ mx+my)%coor( 1 ) = vend(1)
        p%point_inter( 1+ mx+my)%coor( 2 ) = vend(2)
        p%point_inter( 1+ mx+my)%ixy_type = 0
        p%point_inter( 1+ mx+my)%idx2 = 2*ib
        p%point_inter( 1+ mx+my)%idy2 = 2*id
        !**************************************************************************
        ! we will got kstart+mx+my-1 - kstart +1 =mx+my intersections
        ! for example, kstart = 1 to kend = mx+my = kstart + mx + my - 1
         
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = -1.
        xi_e = 1.
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
        
        
    elseif( p%nsubface==2 )then
        ia = p%point_origin%id(1)
        ib = p%point_mid1%id(1)

        ic = p%point_origin%id(2)
        id = p%point_mid1%id(2)
        vorigin(1:2) = p%point_origin%coor(1:2)
        vend(1:2) = p%point_mid1%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
        kstart = p%nsub_outer
        
        p%point_inter( p%nsub_outer)%coor( 1 ) = vorigin(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vorigin(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ia
        p%point_inter( p%nsub_outer)%idy2 = 2*ic
        !******
        p%nsub_outer = p%nsub_outer+ 1+ mx+my
        !******
        p%point_inter( p%nsub_outer)%coor( 1 ) = vend(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vend(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ib
        p%point_inter( p%nsub_outer)%idy2 = 2*id
        !**************************************************************************
        ! we will got kstart+mx+my-1 - kstart +1 =mx+my intersections
        ! for example, kstart = 1 to kend = mx+my = kstart + mx + my - 1
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = -1.
        xi_e = p%point_mid1%xxii
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
        
         
        
        !***************************************************************************
        !***************************************************************************
        ia = p%point_mid1%id(1)
        ib = p%point_end%id(1)

        ic = p%point_mid1%id(2)
        id = p%point_end%id(2)
        vorigin(1:2) = p%point_mid1%coor(1:2)
        vend(1:2) = p%point_end%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
        kstart = p%nsub_outer
        p%point_inter( p%nsub_outer)%coor( 1 ) = vorigin(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vorigin(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ia
        p%point_inter( p%nsub_outer)%idy2 = 2*ic
        !******
        p%nsub_outer = p%nsub_outer+ 1+ mx+my
        !******
        p%point_inter( p%nsub_outer)%coor( 1 ) = vend(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vend(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ib
        p%point_inter( p%nsub_outer)%idy2 = 2*id
        !**************************************************************************
        ! we will got kstart+mx+my-1 - kstart +1 =mx+my intersections
        ! for example, kstart = 1 to kend = mx+my = kstart + mx + my - 1
 
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = p%point_mid1%xxii
        xi_e = 1.
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
        
 
        
    elseif( p%nsubface==3 )then
        ia = p%point_origin%id(1)
        ib = p%point_mid1%id(1)

        ic = p%point_origin%id(2)
        id = p%point_mid1%id(2)
        vorigin(1:2) = p%point_origin%coor(1:2)
        vend(1:2) = p%point_mid1%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
        kstart = p%nsub_outer
        p%point_inter( p%nsub_outer)%coor( 1 ) = vorigin(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vorigin(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ia
        p%point_inter( p%nsub_outer)%idy2 = 2*ic
        !******
        p%nsub_outer = p%nsub_outer+ 1+ mx+my
        !******
        p%point_inter( p%nsub_outer)%coor( 1 ) = vend(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vend(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ib
        p%point_inter( p%nsub_outer)%idy2 = 2*id
        !**************************************************************************
        ! we will got kstart+mx+my-1 - kstart +1 =mx+my intersections
        ! for example, kstart = 1 to kend = mx+my = kstart + mx + my - 1
 
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = -1.
        xi_e = p%point_mid1%xxii
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
        
 
        !***************************************************************************
        !***************************************************************************
        ia = p%point_mid1%id(1)
        ib = p%point_mid2%id(1)

        ic = p%point_mid1%id(2)
        id = p%point_mid2%id(2)
        vorigin(1:2) = p%point_mid1%coor(1:2)
        vend(1:2) = p%point_mid2%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
       kstart = p%nsub_outer 
        p%point_inter( p%nsub_outer)%coor( 1 ) = vorigin(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vorigin(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ia
        p%point_inter( p%nsub_outer)%idy2 = 2*ic
        !******
        p%nsub_outer = p%nsub_outer+ 1+ mx+my
        !******
        p%point_inter( p%nsub_outer)%coor( 1 ) = vend(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vend(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ib
        p%point_inter( p%nsub_outer)%idy2 = 2*id
        !**************************************************************************
 
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = p%point_mid1%xxii
        xi_e = p%point_mid2%xxii
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
 
        !***************************************************************************
        !***************************************************************************
        ia = p%point_mid2%id(1)
        ib = p%point_end%id(1)

        ic = p%point_mid2%id(2)
        id = p%point_end%id(2)
        vorigin(1:2) = p%point_mid2%coor(1:2)
        vend(1:2) = p%point_end%coor(1:2)

        mx = abs(ib-ia)
        my = abs(id-ic)
        
        kstart = p%nsub_outer
        p%point_inter( p%nsub_outer)%coor( 1 ) = vorigin(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vorigin(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ia
        p%point_inter( p%nsub_outer)%idy2 = 2*ic
        !******
        p%nsub_outer = p%nsub_outer+ 1+ mx+my
        !******
        p%point_inter( p%nsub_outer)%coor( 1 ) = vend(1)
        p%point_inter( p%nsub_outer)%coor( 2 ) = vend(2)
        p%point_inter( p%nsub_outer)%ixy_type = 0
        p%point_inter( p%nsub_outer)%idx2 = 2*ib
        p%point_inter( p%nsub_outer)%idy2 = 2*id
        !**************************************************************************
        ! we will got kstart+mx+my-1 - kstart +1 =mx+my intersections
        ! for example, kstart = 1 to kend = mx+my = kstart + mx + my - 1
 
        !call get_mono_inters(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),p,kstart)
        
        xi_o = p%point_mid2%xxii
        xi_e = 1.
        call get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin(1:2),vend(1:2),xi_o,xi_e,p,kstart,v11,v33)
 
    endif

    end subroutine get_intersections_qc
!--------------------------------------------------~ 
subroutine get_mono_inters_QC(ia,ib,ic,id,mx,my,vorigin,vend,xi_o,xi_e,p,kstart,v11,v33)
    implicit none
    integer,intent(in) :: ia,ib,ic,id,mx,my
    real,intent(in) :: vorigin(1:2),vend(1:2),xi_o,xi_e
    real,intent(in) :: v11(1:2),v33(1:2)
    integer,intent(inout) :: kstart
    integer :: ii,jj
    integer :: k,inter,kk,ki
    integer :: ia2,ib2,ic2,id2
    type(type_face),pointer :: p
    integer :: kend

    real :: cxi
    real :: t,tc
    !debug
    real :: td

    integer :: istraight

    real :: t11,t22

    !**************************************************************************
    ! Algorithm of getting the intersection points between face and gridlines
    ! find intersections between face and x=x_i


    if(mx .ne. 0)then

        do k = 1 , mx
            if(ia<ib) then
                inter = ia + k
            else
                inter = ia + 1 - k
            endif

            kk = kstart +k
            p%point_inter( kk)%coor( 1 ) = xgrid( inter )
            call get_kexi_xgrid( v11(1:2), v33(1:2),p%et2(1:2), xgrid( inter ),xi_o,xi_e,cxi,istraight  )


            t = quadratic_function( p%ya,p%yb,p%yc,cxi )
            !t11 = quadratic_function( p%ya,p%yb,p%yc,cxi )
            if( istraight==1 .or. p%ist==1 )then
                t11 = t
                call get_kexi_xgrid_s(v11(1:2), v33(1:2),p%et2(1:2), xgrid( inter ),xi_o,xi_e,cxi)
 
                t22 = reverse_to_y( v11(1:2), v33(1:2),cxi,p%et2(2) )

                if( abs(t11-t22)<0.001*dy   )then
                    t = t11
      
                else
                    t = t22
                endif
                 t = t22
  
            endif



            call confine( vorigin(2), vend(2),t,tc )
            p%point_inter( kk)%coor( 2 ) = tc
 
            td = yside( vorigin(1:2), vend(1:2),  xgrid(inter) )
 
            p%point_inter( kk)%xxii = cxi
            p%point_inter( kk)%ixy_type = 1
            p%point_inter( kk)%igrid = inter
            p%point_inter( kk)%idx2 = 2*inter-1
            p%point_inter( kk)%id_xy = ceiling( (p%point_inter( kk)%coor( 2 )-ybottom)/dy )
            p%point_inter( kk)%idy2 = 2*p%point_inter( kk)%id_xy
 
        enddo
    endif
    !find intersections between face and y=y_j

    if(my .ne. 0)then
 
        do k = 1 + mx, mx+my
            ki = k - mx
            if(ic<id)then
                inter = ic + ki
            else
                inter = ic  + 1 - ki
            endif

            kk = kstart +k
            p%point_inter( kk)%coor( 2 ) = ygrid( inter )
            call get_kexi_ygrid( v11(1:2), v33(1:2),p%et2(1:2), ygrid( inter ),xi_o,xi_e,cxi,istraight  )
            t = quadratic_function( p%xa,p%xb,p%xc,cxi )
 
            if( istraight==1 .or. p%ist==1 )then
                t11 =t
                call get_kexi_ygrid_s(v11(1:2), v33(1:2),p%et2(1:2), ygrid( inter ),xi_o,xi_e,cxi)
  
                t22 = reverse_to_x( v11(1:2), v33(1:2),cxi,p%et2(2) )

                if( abs(t11-t22)<0.001*dx  )then
                    t = t11
 
                else
                    t = t22
                endif
                t = t22
 
            endif




            call confine( vorigin(1), vend(1),t,tc )
            p%point_inter( kk)%coor( 1 ) = tc
 
            td = xside( vorigin(1:2), vend(1:2),  ygrid(inter) )
 
            p%point_inter( kk)%xxii = cxi
            p%point_inter( kk)%ixy_type = 2
            p%point_inter( kk)%igrid = inter
            p%point_inter( kk)%idy2 = 2*inter-1
            p%point_inter( kk)%id_xy = ceiling( (p%point_inter( kk)%coor( 1 )-xleft)/dx )
            p%point_inter( kk)%idx2 = 2*p%point_inter( kk)%id_xy
 

        enddo
    endif
     !
    if(mx+my>1)then
         !sort them in a increasing order
        do ii = 1 ,mx+my-1
            do jj = kstart+1 , kstart+mx+my-ii
                if( p%point_inter(jj)%xxii>p%point_inter(jj+1)%xxii )then
                    temp = p%point_inter(jj+1)
                    p%point_inter(jj+1) = p%point_inter(jj)
                    p%point_inter(jj) = temp
                endif
            enddo
        enddo
    endif

    ia2 = 2*ia; ib2 = 2*ib; ic2 = 2*ic; id2 = 2*id;
    !if( mx+my>0 )then
    if( ia2<=ib2 )then
        do kk = kstart+1,kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2 +2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
 
                elseif( p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 .and. p%point_inter(kk)%idx2>p%point_inter(kk+1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
 
                endif

            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my,kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk+1)%idx2<p%point_inter(kk)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
  
                    
                elseif( p%point_inter(kk+1)%idx2>p%point_inter(kk)%idx2 .and. p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
 
                endif

            endif

        enddo

    elseif( ia2>ib2 )then
        do kk = kstart+1, kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                
                elseif( p%point_inter(kk)%idx2<p%point_inter(kk-1)%idx2 .and. p%point_inter(kk)%idx2<p%point_inter(kk+1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                 
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my, kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 2 )then
                if( p%point_inter(kk+1)%idx2>p%point_inter(kk)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2 +2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                    
                elseif( p%point_inter(kk+1)%idx2<p%point_inter(kk)%idx2 .and. p%point_inter(kk)%idx2>p%point_inter(kk-1)%idx2 )then
                    p%point_inter(kk)%idx2 = (p%point_inter(kk)%idx2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idx2/2
                  
                endif
            endif

        enddo
    endif

    !************************
    if( ic2<=id2 )then
        do kk = kstart+1, kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2 +2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
               
                elseif( p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 .and. p%point_inter(kk)%idy2>p%point_inter(kk+1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
               
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my,kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk+1)%idy2<p%point_inter(kk)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
   
                elseif( p%point_inter(kk+1)%idy2>p%point_inter(kk)%idy2 .and. p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
         
                endif

            endif

        enddo
    elseif( ic2>id2 )then
        do kk = kstart+1,kstart+mx+my
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2 -1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
           
                elseif( p%point_inter(kk)%idy2<p%point_inter(kk-1)%idy2 .and. p%point_inter(kk)%idy2<p%point_inter(kk+1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2+2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
           
                endif
            endif

        enddo
        !******************reverse correction
        do kk = kstart+mx+my, kstart+1,-1
            if( p%point_inter(kk)%ixy_type == 1 )then
                if( p%point_inter(kk+1)%idy2>p%point_inter(kk)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2 +2)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
     
                elseif( p%point_inter(kk+1)%idy2<p%point_inter(kk)%idy2 .and. p%point_inter(kk)%idy2>p%point_inter(kk-1)%idy2 )then
                    p%point_inter(kk)%idy2 = (p%point_inter(kk)%idy2-1)/2*2
                    p%point_inter(kk)%id_xy = p%point_inter(kk)%idy2/2
        
                endif
            endif

        enddo
    endif


    !endif


    !*********************************************************************************
    ! check the monotonicity
    ia2 = 2*ia; ib2 = 2*ib; ic2 = 2*ic; id2 = 2*id;

    if( mx+my>0 )then
        if(ia2<=ib2)then
            if( p%point_inter(kstart+1)%idx2<ia2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!5'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idx2<p%point_inter(k-1)%idx2 )then
                    write(911,*) 'lose monotone!'
 
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idx2>ib2 )then
                write(911,*) 'lose monotone!'
 
                print *,'lose monotone!4'
                pause
            endif
        elseif( ia2>ib2 )then
            if( p%point_inter(kstart+1)%idx2>ia2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!6'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idx2>p%point_inter(k-1)%idx2 )then
                    write(911,*) 'lose monotone!'
 
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idx2<ib2 )then
                write(911,*) 'lose monotone!'
                print *,'lose monotone!7'
                pause
            endif
        endif
        !*************************************************
        if(ic2<=id2)then
            if( p%point_inter(kstart+1)%idy2<ic2 )then
                write(911,*) 'lose monotone!'
 
                print *,'lose monotone!8'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idy2<p%point_inter(k-1)%idy2 )then
                    write(911,*) 'lose monotone!'
 
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idy2>id2 )then
                write(911,*) 'lose monotone!'
                pause
            endif
        elseif( ic2>id2 )then
            if( p%point_inter(kstart+1)%idy2>ic2 )then
                write(911,*) 'lose monotone!'
                pause
            endif
            do k = kstart+2,kstart+mx+my
                if( p%point_inter(k)%idy2>p%point_inter(k-1)%idy2 )then
                    write(911,*) 'lose monotone!'
                    pause
                endif
            enddo
            if( p%point_inter(kstart+mx+my)%idy2<id2 )then
                write(911,*) 'lose monotone!'
                pause
            endif
        endif
    endif
    !*********************************************************************************

    end subroutine get_mono_inters_QC
!--------------------------------------------------~  
subroutine get_kexi_xgrid( v1,v3,e2, x, xorigin,xend, cxi,istraight )
    implicit none
    real,intent(in) :: v1(2),v3(2),e2(2),x,xorigin,xend
    real,intent(out) :: cxi
    integer,intent(out) :: istraight

    real :: at,bt,ct
    real :: t1,t2,t3
    real :: xi(2),eta(1:2)
    !real :: xorigin,xend

    if( abs( v3(1)-v1(1) )<=abs( v3(2)-v1(2) ) )then
        at = e2(2)/( e2(1)**2-1. )
        bt = ( v3(1)-v1(1) )/( v3(2)-v1(2) )
        ct =  2.*( ( v1(1)+v3(1) )/2.-x)/( v3(2)-v1(2) )-at
        
        if( abs(at)>eps )then 
        call get_quadratic_roots( at,bt,ct,xi(1:2) )
        elseif( abs(at)<=eps .and. abs(bt)>eps )then
       
            xi(1:2) = -ct/bt 
        else
            xi(1:2) = xorigin*0.5 + xend*0.5
        endif
        !pause
    else
        t1 = e2(2)/( e2(1)**2-1. )
        t2 = ( v3(2)-v1(2) )/( v3(1)-v1(1) )
        t3 = ( x-0.5*v1(1)-0.5*v3(1) )/( v3(1)-v1(1) )
        at = t1*t2**2
        bt = -1. - 4.*t3*t2*t1
        ct = t1*( 4.*t3**2-1. )
        
        if( abs(at)>eps )then
        call get_quadratic_roots( at,bt,ct,eta(1:2) )
        xi(1) = 2*t3-t2*eta(1)
        xi(2) = 2*t3-t2*eta(2)
        elseif( abs(at)<=eps .and. abs(bt)>eps )then
            eta(1:2) = -ct/bt
            xi(1:2) = 2.*t3-t2*eta(1:2)
        else
            xi(1:2) = xorigin*0.5 + xend*0.5
        endif        

    endif

    ! pick one root!
    !xorigin = -1.
    !xend = 1.
    call pick_1root( xi(1:2),xorigin,xend,cxi,istraight )
    
    if( bt**2 - 4.*at*ct <-eps )then
        istraight = 1
        print *,221, bt**2 - 4.*at*ct,-eps, at
        pause
    endif    
    
    if(istraight==1)then
        print *,v1(1:2),v3(1:2),e2(1:2),x,xorigin,xend
        print *,xi(1:2),xorigin,xend,cxi,istraight
        print *,abs( v3(1)-v1(1) ),abs( v3(2)-v1(2) )
        print *,eta(1:2)
        print *,2222,t3,t2,eta(1)
        pause
    endif
    end subroutine get_kexi_xgrid
!--------------------------------------------------~  
subroutine get_kexi_ygrid( v1,v3,e2, y,xorigin,xend,cxi,istraight )
    implicit none
    real,intent(in) :: v1(2),v3(2),e2(2),y,xorigin,xend
    real,intent(out) :: cxi
    integer,intent(out) :: istraight

    real :: at,bt,ct
    real :: t1,t2,t3
    real :: xi(2),eta(1:2)
    !real :: xorigin,xend

    if( abs( v3(1)-v1(1) )<=abs( v3(2)-v1(2) ) )then
        t1 = e2(2)/( e2(1)**2-1. )
        t2 = ( v3(1)-v1(1) )/( v3(2)-v1(2) )
        t3 = ( y-0.5*v1(2)-0.5*v3(2) )/( v3(2)-v1(2) )
        at = t1*t2**2
        bt = -1.+4.*t3*t2*t1
        ct =  4.*t1*t3**2  -t1
        
        if( abs(at)>eps )then
        call get_quadratic_roots( at,bt,ct,eta(1:2) )
        xi(1) = 2.*t3+t2*eta(1)
        xi(2) = 2.*t3+t2*eta(2)
        elseif( abs(at)<=eps .and. abs(bt)>eps )then
            eta(1:2) = -ct/bt
            xi(1:2) = 2.*t3+t2*eta(1:2)
        else
            xi(1:2) = xorigin*0.5 + xend*0.5
        endif

    else

        at = - e2(2)/( e2(1)**2-1. )
        bt = ( v3(2)-v1(2) )/(v3(1)-v1(1) )
        ct = 2*( ( v1(2)+v3(2) )/2.-y)/(v3(1)-v1(1) )-at
        if( abs(at)>eps )then 
        call get_quadratic_roots( at,bt,ct,xi(1:2) )
        elseif( abs(at)<=eps .and. abs(bt)>eps )then
       
            xi(1:2) = -ct/bt 
        else
            xi(1:2) = xorigin*0.5 + xend*0.5
        endif
 
    endif
 
    call pick_1root( xi(1:2),xorigin,xend,cxi,istraight )
    
 
    
    if( bt**2 - 4.*at*ct <-eps )then
        istraight = 1
        print *,abs( v3(1)-v1(1) ),abs( v3(2)-v1(2) )
        print *,223, bt**2 - 4.*at*ct,-eps, at
        pause
    endif
    
    end subroutine get_kexi_ygrid
!--------------------------------------------------~  
subroutine pick_1root( xi,xorigin,xend,cxi,istraight )
    implicit none
    real,intent(in) :: xi(1:2),xorigin,xend
    real,intent(out) :: cxi
    integer :: istraight

    istraight = 0
    if( ( xi(1) - xorigin )*( xi(1) - xend ) <=0. .and. ( xi(2) - xorigin )*( xi(2) - xend )>0.  )then
        cxi = xi(1)
    elseif( ( xi(1) - xorigin )*( xi(1) - xend ) >0. .and. ( xi(2) - xorigin )*( xi(2) - xend )<=0.  )then
        cxi = xi(2)
    elseif( ( xi(1) - xorigin )*( xi(1) - xend ) <=0. .and. ( xi(2) - xorigin )*( xi(2) - xend )<=0. )then
        cxi = xi(1)*0.5 + xi(2)*0.5
  
    elseif( xi(1)>xend .and. xi(2)>xend )then
        cxi = xend
    elseif( xi(1)<xorigin .and. xi(2)<xorigin )then
        cxi = xorigin
        ! special cases
    elseif(  xi(1)<xorigin .and. xi(2)>xend )then
        if( abs( xi(1)-xorigin )< abs( xi(2)-xend) )then
        cxi = xorigin
        elseif( abs( xi(1)-xorigin )>= abs( xi(2)-xend) )then
            cxi = xend
        endif
    elseif(  xi(2)<xorigin .and. xi(1)>xend)then
        if( abs( xi(2)-xorigin )< abs( xi(1)-xend) )then
        cxi = xorigin
        elseif( abs( xi(2)-xorigin )>= abs( xi(1)-xend) )then
            cxi = xend
        endif
 
    endif

    end subroutine pick_1root
!--------------------------------------------------~  
subroutine confine( t1,t2,t,tc )
    implicit none
    real,intent(in) :: t1,t2,t
    real,intent(out) :: tc
    real :: tmin,tmax

    tmin = min( t1,t2 )
    tmax = max( t1,t2 )
    if( t<tmin )then
        tc = tmin
    elseif( t>tmax )then
        tc = tmax
    else
        tc = t
    endif

    end subroutine confine
!--------------------------------------------------~  
subroutine confine3( t1,t2,t3,t,tc )
    implicit none
    real,intent(in) :: t1,t2,t3,t
    real,intent(out) :: tc
    real :: tmin,tmax

    tmin = min( t1,t2,t3 )
    tmax = max( t1,t2,t3 )
    if( t<tmin )then
        tc = tmin
    elseif( t>tmax )then
        tc = tmax
    else
        tc = t
    endif

    end subroutine confine3
!--------------------------------------------------~  
subroutine get_quadratic_roots( a,b,c, xi )
    implicit none
    real,intent(in) :: a,b,c
    real,intent(out) :: xi(2)

    real :: gamma,delta

    delta = abs(b**2-4.*a*c)

    if(b>=0.)then
        gamma = 1.
    else
        gamma = -1.
    endif

    xi(1) = 2.*c/( -b - gamma *sqrt(delta) )
    xi(2) = ( -b-gamma*sqrt(delta)  )/a*0.5

    end subroutine get_quadratic_roots
!--------------------------------------------------~  
subroutine get_kexi_xgrid_s( v1,v3,e2, x, xorigin,xend, cxi  )
    implicit none
    real,intent(in) :: v1(1:2),v3(1:2),e2(1:2),x,xorigin,xend
    real,intent(out) :: cxi


    real :: xi

    xi = 2.*( x-(v1(1)+v3(1))/2.-(v3(2)-v1(2))/2.*e2(2)   )/( v3(1)-v1(1) )

    call confine(  xorigin,xend,xi,cxi  )

    end subroutine get_kexi_xgrid_s
!--------------------------------------------------~  
subroutine get_kexi_ygrid_s( v1,v3,e2, y, xorigin,xend, cxi  )
    implicit none
    real,intent(in) :: v1(1:2),v3(1:2),e2(1:2),y,xorigin,xend
    real,intent(out) :: cxi


    real :: xi

    xi = 2.*( y-(v1(2)+v3(2))/2.+(v3(1)-v1(1))/2.*e2(2)   )/( v3(2)-v1(2) )

    call confine(  xorigin,xend,xi,cxi  )

    end subroutine get_kexi_ygrid_s
!--------------------------------------------------~  
real function reverse_to_x(v1,v3, exi20,eta20)
    implicit none
    real,intent(in) ::  v1(1:2),v3(1:2), exi20,eta20

    reverse_to_x = ( v3(1)-v1(1) )/2.*exi20 + ( v3(2)-v1(2) )/2.*eta20 + ( v1(1)+v3(1) )/2.


    end function reverse_to_x
!--------------------------------------------------~  
real function reverse_to_y( v1,v3, exi20,eta20 )
    implicit none
    real,intent(in) :: v1(1:2),v3(1:2), exi20,eta20

    reverse_to_y = ( v3(2)-v1(2) )/2.*exi20 - (  v3(1)-v1(1) )/2.*eta20 + (  v1(2)+v3(2) )/2.

    end function reverse_to_y
!--------------------------------------------------~  
subroutine get_subfaces(p)
    implicit none
    type(type_face),pointer :: p
    real :: v1(1:2),v2(1:2),v3(1:2)
    real :: a_trans,b_trans,c_trans,d_trans

    real :: xt1,yt1,xf1,yf1
    real :: xt2,yt2,xf2,yf2
    type :: type_extreme
        sequence
        real kexi
        integer itype
    end type

    type(type_extreme) extreme(1:2),extreme_t

    real :: eh,del

    v1(1:2) = p%point_origin%coor(1:2)
    v2(1:2) = p%point_midt%coor(1:2)
    v3(1:2) = p%point_end%coor(1:2)

    call trans(v1,v3,a_trans,b_trans,c_trans,d_trans)

    p%et2(1) = trans_to_ex( a_trans,b_trans,c_trans,v2 )
    p%et2(2) = trans_to_ey( a_trans,b_trans,d_trans,v2 )

    ! for x(\xi)
    p%xa = ( v3(2)-v1(2) )/2.*p%et2(2)/(p%et2(1)**2-1)
    p%xb = ( v3(1)-v1(1) )/2.
    p%xc = ( v3(1)+v1(1) )/2.-p%xa
    ! for y(\xi)
    p%ya = -( v3(1)-v1(1) )/2.*p%et2(2)/(p%et2(1)**2-1)
    p%yb = ( v3(2)-v1(2) )/2.
    p%yc = ( v3(2)+v1(2) )/2. - p%ya


    p%nsubface = 1


    if(  abs(p%xa) <=eps )then
        ! do nothing
    elseif(  abs(p%xa) >eps )then
        eh = -p%xb/(2*p%xa)
        !del =   p%xb**2-4*p%xa*p%xc  
        !if( eh>-1. .and. eh<1. .and. del>eps .and. -del/4./p%xa<-eps  )then
        if( eh>-1. .and. eh<1.    )then
            extreme( p%nsubface )%kexi = eh
            extreme( p%nsubface )%itype = 1
            p%nsubface = p%nsubface + 1
        endif
    endif


    if(  abs(p%ya) <=eps )then
        ! do nothing
    elseif(  abs(p%ya) >eps )then
        eh = -p%yb/(2*p%ya)
        !del =   p%yb**2-4*p%ya*p%yc 
        !if( eh>-1. .and. eh<1.  .and. del>eps .and. -del/4./p%ya<-eps   )then
        if( eh>-1. .and. eh<1.       )then
            extreme( p%nsubface )%kexi = eh
            extreme( p%nsubface )%itype = 2
            p%nsubface = p%nsubface + 1
        endif
    endif
    
    !if( abs(  p%et2(2)  )<=1000*eps .and. p%nsubface>1 )then
    !    print *,io,i,j,p%et2(1:2),p%nsubface
    !    pause
    !endif
    

    
    
    p%ist = 0
    !if( abs(  p%et2(2)  )<=100*eps   )then
    !    p%nsubface = 1
    !    p%ist = 1
    !endif 
     
    !p%nsubface = 1
    !*******************************************************
    ! Assemble the extreme points of a quadratic-curved face
    !*******************************************************
    if( p%nsubface==1 )then
        ! do nothing!
    elseif( p%nsubface==2 )then
        xt1 = quadratic_function( p%xa,p%xb,p%xc,extreme(1)%kexi )
        yt1 = quadratic_function( p%ya,p%yb,p%yc,extreme(1)%kexi )

        p%point_mid1%xxii = extreme(1)%kexi
 
        if( extreme(1)%itype==1 )then   
            xf1 = xt1
            ! confine y
            call confine( v1(2),v3(2),yt1,yf1 )
        elseif( extreme(1)%itype==2 )then
            ! confine x
            call confine( v1(1),v3(1),xt1,xf1)
            yf1 = yt1
        endif
        
        p%point_mid1%coor(1) = xf1
        p%point_mid1%coor(2) = yf1

        p%point_mid1%id(1) = ceiling( (xf1-xleft)/dx )
        p%point_mid1%id(2) = ceiling( (yf1-ybottom)/dy )
  
        !print *,extreme(1)
         !pause
    elseif( p%nsubface==3 )then
        if( extreme(1)%kexi>extreme(2)%kexi )then
            extreme_t = extreme(1)
            extreme(1) = extreme(2)
            extreme(2) = extreme_t
        endif
        xt1 = quadratic_function( p%xa,p%xb,p%xc,extreme(1)%kexi )
        yt1 = quadratic_function( p%ya,p%yb,p%yc,extreme(1)%kexi )
        xt2 = quadratic_function( p%xa,p%xb,p%xc,extreme(2)%kexi )
        yt2 = quadratic_function( p%ya,p%yb,p%yc,extreme(2)%kexi )
 
        
        if( extreme(1)%itype==1 .and. extreme(2)%itype==2 )then
            xf1 = xt1
            yf2 = yt2
            
            call confine3( v1(2),v3(2),yf2, yt1,yf1 )
            call confine3( v1(1),v3(1),xf1, xt2,xf2 )
        elseif( extreme(1)%itype==2 .and. extreme(2)%itype==1 )then
            yf1 = yt1
            xf2 = xt2 
            
            call confine3( v1(2),v3(2),yf1, yt2,yf2 )
            call confine3( v1(1),v3(1),xf2, xt1,xf1 )            
        endif
        
        p%point_mid1%xxii = extreme(1)%kexi
        p%point_mid1%coor(1) = xf1
        p%point_mid1%coor(2) = yf1
        p%point_mid1%id(1) = ceiling( (xf1-xleft)/dx )
        p%point_mid1%id(2) = ceiling( (yf1-ybottom)/dy )



        p%point_mid2%xxii = extreme(2)%kexi
        p%point_mid2%coor(1) = xf2
        p%point_mid2%coor(2) = yf2
        p%point_mid2%id(1) = ceiling( (xf2-xleft)/dx )
        p%point_mid2%id(2) = ceiling( (yf2-ybottom)/dy )

    endif


    end subroutine get_subfaces
!--------------------------------------------------~  
subroutine trans( v1,v3,a_trans,b_trans,c_trans,d_trans )
    implicit none
    real,intent(in) :: v1(1:2),v3(1:2)
    real,intent(out) :: a_trans,b_trans,c_trans,d_trans
    real :: x1,y1,x2,y2

    x1 = v1(1); y1 = v1(2); x2 = v3(1); y2 = v3(2)

    if( abs(x2-x1)< abs(y2-y1) )then
        a_trans = 2.*( (x2-x1)/(y2-y1)/(y2-y1) ) / (   ( (x2-x1)/(y2-y1) )**2  + 1.    )
        b_trans = 2./( y2-y1 )  /  (   ( (x2-x1)/(y2-y1) )**2  + 1.    )
        c_trans = (  (x1-x2)/(y1-y2)*(x1+x2)/(y1-y2)  + (y1+y2)/(y1-y2)  ) /  (   ( (x2-x1)/(y2-y1) )**2  + 1.    )
        d_trans = (  (x1+x2)/(y1-y2)  +  (y1+y2)/(y2-y1)*(x2-x1)/(y2-y1)   )  /  (   ( (x2-x1)/(y2-y1) )**2  + 1.    )
    else
        a_trans = 2./(x2-x1) / (   1.+(   (y2-y1)/(x2-x1)  )**2    )
        b_trans = 2.*(y2-y1)/(x2-x1)/(x2-x1)   / (   1.+(   (y2-y1)/(x2-x1)  )**2    )
        c_trans = (  (x1+x2)/(x1-x2) + (y1-y2)/(x1-x2)*(y1+y2)/(x1-x2)  )  / (   1.+(   (y2-y1)/(x2-x1)  )**2    )
        d_trans = (  (x1+x2)/(x1-x2)*(y1-y2)/(x1-x2)  + (y1+y2)/(x2-x1)  )  / (   1.+(   (y2-y1)/(x2-x1)  )**2    )
    endif
    end subroutine trans
!--------------------------------------------------~  
real function trans_to_ex( a,b,c,v )
    implicit none
    real,intent(in) :: a,b,c,v(1:2)

    trans_to_ex = a*v(1) + b*v(2) + c

    end function trans_to_ex
!--------------------------------------------------~  
real function trans_to_ey( a,b,d,v )
    implicit none
    real,intent(in) :: a,b,d,v(1:2)

    trans_to_ey = b*v(1) - a*v(2) + d

    end function trans_to_ey
!--------------------------------------------------~  
real function quadratic_function(a,b,c,x)
    implicit none
    real,intent(in) :: a,b,c,x

    quadratic_function = a*x**2+b*x+c

    end function quadratic_function
    
    
!<3><get innersegments>
subroutine get_innersegments
    implicit none
    type(type_element_upstream),pointer :: pes
    type(type_face),pointer :: pf
    integer :: iy0,nby
    integer :: ix0,nbx
    integer :: ii,jj
    type(type_inter_point) :: superpoint1,superpoint2

    real :: xmid_t,ymid_t

    do i = 1 , nx
        do j = 1 , ny
            !Determine an upstream element's four vertexes and four faces.
            pes => element_star(i,j,io)
            pes%vertex1 => vertex_star(i,j)
            pes%vertex2 => vertex_star(i+1,j)
            pes%vertex3 => vertex_star(i+1,j+1)
            pes%vertex4 => vertex_star(i,j+1)

            pes%edge1_lr => face_lr(i,j)
            pes%edge2_bt => face_bt(i+1,j)
            pes%edge3_lr => face_lr(i,j+1)
            pes%edge4_bt => face_bt(i,j)

            !End determine an upstream element's four vertexes and four faces.

            if( iqc==0 )then
                iam = min(pes%vertex1%id(1),pes%vertex2%id(1),pes%vertex3%id(1),pes%vertex4%id(1) )
                ibm = max(pes%vertex1%id(1),pes%vertex2%id(1),pes%vertex3%id(1),pes%vertex4%id(1) )
                icm = min(pes%vertex1%id(2),pes%vertex2%id(2),pes%vertex3%id(2),pes%vertex4%id(2) )
                idm = max(pes%vertex1%id(2),pes%vertex2%id(2),pes%vertex3%id(2),pes%vertex4%id(2) )
            elseif( iqc==1 )then
                iam = min( pes%edge1_lr%point_origin%id(1),pes%edge1_lr%point_end%id(1) )
                ibm = max( pes%edge1_lr%point_origin%id(1),pes%edge1_lr%point_end%id(1) )
                icm = min( pes%edge1_lr%point_origin%id(2),pes%edge1_lr%point_end%id(2) )
                idm = max( pes%edge1_lr%point_origin%id(2),pes%edge1_lr%point_end%id(2) )
                if( pes%edge1_lr%nsubface==2 )then
                    iam = min( iam,pes%edge1_lr%point_mid1%id(1) )
                    ibm = max( ibm,pes%edge1_lr%point_mid1%id(1) )
                    icm = min( icm,pes%edge1_lr%point_mid1%id(2) )
                    idm = max( idm,pes%edge1_lr%point_mid1%id(2) )
                elseif( pes%edge1_lr%nsubface==3 )then
                    iam = min( iam,pes%edge1_lr%point_mid1%id(1),pes%edge1_lr%point_mid2%id(1) )
                    ibm = max( ibm,pes%edge1_lr%point_mid1%id(1),pes%edge1_lr%point_mid2%id(1) )
                    icm = min( icm,pes%edge1_lr%point_mid1%id(2),pes%edge1_lr%point_mid2%id(2) )
                    idm = max( idm,pes%edge1_lr%point_mid1%id(2),pes%edge1_lr%point_mid2%id(2) )
                endif
                !*******************
                iam = min( iam,pes%edge2_bt%point_origin%id(1),pes%edge2_bt%point_end%id(1) )
                ibm = max( ibm,pes%edge2_bt%point_origin%id(1),pes%edge2_bt%point_end%id(1) )
                icm = min( icm,pes%edge2_bt%point_origin%id(2),pes%edge2_bt%point_end%id(2) )
                idm = max( idm,pes%edge2_bt%point_origin%id(2),pes%edge2_bt%point_end%id(2) )
                if( pes%edge2_bt%nsubface==2 )then
                    iam = min( iam,pes%edge2_bt%point_mid1%id(1) )
                    ibm = max( ibm,pes%edge2_bt%point_mid1%id(1) )
                    icm = min( icm,pes%edge2_bt%point_mid1%id(2) )
                    idm = max( idm,pes%edge2_bt%point_mid1%id(2) )
                elseif( pes%edge2_bt%nsubface==3 )then
                    iam = min( iam,pes%edge2_bt%point_mid1%id(1),pes%edge2_bt%point_mid2%id(1) )
                    ibm = max( ibm,pes%edge2_bt%point_mid1%id(1),pes%edge2_bt%point_mid2%id(1) )
                    icm = min( icm,pes%edge2_bt%point_mid1%id(2),pes%edge2_bt%point_mid2%id(2) )
                    idm = max( idm,pes%edge2_bt%point_mid1%id(2),pes%edge2_bt%point_mid2%id(2) )
                endif
                !*******************
                iam = min( iam,pes%edge3_lr%point_origin%id(1),pes%edge3_lr%point_end%id(1) )
                ibm = max( ibm,pes%edge3_lr%point_origin%id(1),pes%edge3_lr%point_end%id(1) )
                icm = min( icm,pes%edge3_lr%point_origin%id(2),pes%edge3_lr%point_end%id(2) )
                idm = max( idm,pes%edge3_lr%point_origin%id(2),pes%edge3_lr%point_end%id(2) )
                if( pes%edge3_lr%nsubface==2 )then
                    iam = min( iam,pes%edge3_lr%point_mid1%id(1) )
                    ibm = max( ibm,pes%edge3_lr%point_mid1%id(1) )
                    icm = min( icm,pes%edge3_lr%point_mid1%id(2) )
                    idm = max( idm,pes%edge3_lr%point_mid1%id(2) )
                elseif( pes%edge3_lr%nsubface==3 )then
                    iam = min( iam,pes%edge3_lr%point_mid1%id(1),pes%edge3_lr%point_mid2%id(1) )
                    ibm = max( ibm,pes%edge3_lr%point_mid1%id(1),pes%edge3_lr%point_mid2%id(1) )
                    icm = min( icm,pes%edge3_lr%point_mid1%id(2),pes%edge3_lr%point_mid2%id(2) )
                    idm = max( idm,pes%edge3_lr%point_mid1%id(2),pes%edge3_lr%point_mid2%id(2) )
                endif
                !*******************
                iam = min( iam,pes%edge4_bt%point_origin%id(1),pes%edge4_bt%point_end%id(1) )
                ibm = max( ibm,pes%edge4_bt%point_origin%id(1),pes%edge4_bt%point_end%id(1) )
                icm = min( icm,pes%edge4_bt%point_origin%id(2),pes%edge4_bt%point_end%id(2) )
                idm = max( idm,pes%edge4_bt%point_origin%id(2),pes%edge4_bt%point_end%id(2) )
                if( pes%edge4_bt%nsubface==2 )then
                    iam = min( iam,pes%edge4_bt%point_mid1%id(1) )
                    ibm = max( ibm,pes%edge4_bt%point_mid1%id(1) )
                    icm = min( icm,pes%edge4_bt%point_mid1%id(2) )
                    idm = max( idm,pes%edge4_bt%point_mid1%id(2) )
                elseif( pes%edge4_bt%nsubface==3 )then
                    iam = min( iam,pes%edge4_bt%point_mid1%id(1),pes%edge4_bt%point_mid2%id(1) )
                    ibm = max( ibm,pes%edge4_bt%point_mid1%id(1),pes%edge4_bt%point_mid2%id(1) )
                    icm = min( icm,pes%edge4_bt%point_mid1%id(2),pes%edge4_bt%point_mid2%id(2) )
                    idm = max( idm,pes%edge4_bt%point_mid1%id(2),pes%edge4_bt%point_mid2%id(2) )
                endif
            endif

            isx(:) = 0
            isy(:) = 0
            pf =>pes%edge1_lr
            call super_inner(pf )
            pf =>pes%edge2_bt
            call super_inner(pf )
            pf =>pes%edge3_lr
            call super_inner(pf )
            pf =>pes%edge4_bt
            call super_inner(pf )
            



            pes%nsub_inner = 0
            ! get all inner segments!!
            do ix = 1,ibm-iam
            !**********sorting, if isx(ix)>2
                if( isx(ix)>2 )then
                do ii = 1,3
                    do jj = 1, 4-ii
                        if( point_inner_x( ix,jj)%coor(2)>point_inner_x( ix,jj+1)%coor(2) )then
                temp = point_inner_x( ix,jj+1)
                point_inner_x( ix,jj+1) = point_inner_x( ix,jj)
                point_inner_x( ix,jj) = temp
                        endif
                        
                    enddo
                enddo
                
                endif
                
                superpoint1 = point_inner_x( ix,1)
                superpoint2 = point_inner_x( ix,2)
                call super2final_inner_x(pes,superpoint1,superpoint2)
                if(isx(ix)>2 )then
                    superpoint1 = point_inner_x( ix,3)
                    superpoint2 = point_inner_x( ix,4)
                    call super2final_inner_x(pes,superpoint1,superpoint2)
                endif
            enddo

            !*********************************************************

            do iy = 1,idm-icm
                
                if( isy(iy)>2 )then
                do ii = 1,3
                    do jj = 1, 4-ii
                        if( point_inner_y( iy,jj)%coor(1)>point_inner_y( iy,jj+1)%coor(1) )then
                temp = point_inner_y( iy,jj+1)
                point_inner_y( iy,jj+1) = point_inner_y( iy,jj)
                point_inner_y( iy,jj) = temp
                        endif
                        
                    enddo
                enddo
                
                endif
                
                superpoint1 = point_inner_y( iy,1)
                superpoint2 = point_inner_y( iy,2)
                call super2final_inner_y(pes,superpoint1,superpoint2)
                if( isy(iy)>2 )then
                    superpoint1 = point_inner_y( iy,3)
                    superpoint2 = point_inner_y( iy,4)
                    call super2final_inner_y(pes,superpoint1,superpoint2)
                endif

            enddo

        enddo
    enddo

    end subroutine get_innersegments
!--------------------------------------------------~
subroutine super_inner(pf )
    implicit none
    type(type_face),pointer :: pf

    integer :: k


    do k = 1, pf%nsub_outer-1
        if( pf%point_inter(k)%ixy_type==1 )then
            ix = pf%point_inter(k)%igrid - iam
            if( isx(ix) == 0 )then
                isx(ix) = 1
                point_inner_x( ix,1) = pf%point_inter(k)
            elseif( isx(ix)==1 )then
                isx(ix) = 2
                point_inner_x( ix,2) = pf%point_inter(k)
            elseif( isx(ix)==2 )then
                isx(ix) = 3
                point_inner_x( ix,3) = pf%point_inter(k)
            elseif( isx(ix)==3 )then
                isx(ix) = 4
                point_inner_x( ix,4) = pf%point_inter(k)
            endif
        elseif(pf%point_inter(k)%ixy_type==2)then
            iy = pf%point_inter(k)%igrid - icm
            if( isy(iy) == 0 )then
                isy(iy) = 1
                point_inner_y( iy,1 ) = pf%point_inter(k)
            elseif( isy(iy)==1 )then
                isy(iy) = 2
                point_inner_y( iy,2 ) = pf%point_inter(k)
            elseif( isy(iy)==2 )then
                isy(iy) = 3
                point_inner_y( iy,3 ) = pf%point_inter(k)
            elseif( isy(iy)==3 )then
                isy(iy) = 4
                point_inner_y( iy,4 ) = pf%point_inter(k)
            endif

        endif
    enddo


    end subroutine super_inner
!--------------------------------------------------~
subroutine super2final_inner_x(pes,superpoint1,superpoint2)
    implicit none
    type(type_element_upstream),pointer :: pes
    type(type_inter_point) :: superpoint1,superpoint2    
    integer :: ii,jj
    integer :: nby,iy0
    real :: ymid_t


    iy0 = superpoint1%id_xy
    nby = iy0 - superpoint2%id_xy
 

    point_inner(1) =  superpoint1

    if( nby>0 )then
        do ii = 1,nby
            point_inner(ii+1)%coor(1) = point_inner(1)%coor(1)
            point_inner(ii+1)%coor(2) = ygrid(iy0+1-ii)
            point_inner(ii+1)%ixy_type = 3
            point_inner(ii+1)%idx2 =  point_inner(1)%idx2
            point_inner(ii+1)%idy2 = 2*(iy0+1-ii)-1
        enddo
    elseif( nby<0 )then
        do ii = 1 , abs(nby)
            point_inner(ii+1)%coor(1) = point_inner(1)%coor(1)
            point_inner(ii+1)%coor(2) = ygrid(iy0+ii)
            point_inner(ii+1)%ixy_type = 3
            point_inner(ii+1)%idx2 =  point_inner(1)%idx2
            point_inner(ii+1)%idy2 = 2*(iy0+ii)-1
        enddo
    else

    endif

    point_inner(2+abs(nby) ) = superpoint2

    ! final inner segments
    do ii = 1 , 1 + abs(nby)
        jj = pes%nsub_inner + ii
        pes%segment_inner(jj)%porigin = point_inner(ii)
        pes%segment_inner(jj)%pend = point_inner(ii+1)

        if( pes%segment_inner(jj)%porigin%coor(2)<pes%segment_inner(jj)%pend%coor(2) )then
            ! from bottom to top
            pes%segment_inner(jj)%id(1) = (pes%segment_inner(jj)%porigin%idx2-1)/2
        else
            pes%segment_inner(jj)%id(1) = (pes%segment_inner(jj)%porigin%idx2+1)/2
        endif

        if( pes%segment_inner(jj)%porigin%ixy_type == 3 .and. pes%segment_inner(jj)%pend%ixy_type==3 )then
            ymid_t = (pes%segment_inner(jj)%porigin%coor(2)+pes%segment_inner(jj)%pend%coor(2) )/2.
            pes%segment_inner(jj)%id(2) = ceiling( (ymid_t-ybottom)/dy )
        elseif( pes%segment_inner(jj)%porigin%ixy_type == 1  )then
            pes%segment_inner(jj)%id(2) = pes%segment_inner(jj)%porigin%id_xy
        elseif( pes%segment_inner(jj)%porigin%ixy_type == 3 .and. pes%segment_inner(jj)%pend%ixy_type==1  )then
            pes%segment_inner(jj)%id(2) = pes%segment_inner(jj)%pend%id_xy
        endif



    enddo

    pes%nsub_inner = pes%nsub_inner + abs(nby)+1

    do ii = 1 , 1+abs(nby)
        jj = pes%nsub_inner + ii
        pes%segment_inner(jj)%porigin = point_inner(ii+1)
        pes%segment_inner(jj)%pend = point_inner(ii)

        if( pes%segment_inner(jj)%porigin%coor(2)<pes%segment_inner(jj)%pend%coor(2) )then
            ! from bottom to top
            pes%segment_inner(jj)%id(1) = (pes%segment_inner(jj)%porigin%idx2-1)/2
        else
            pes%segment_inner(jj)%id(1) = (pes%segment_inner(jj)%porigin%idx2+1)/2
        endif

        if( pes%segment_inner(jj)%porigin%ixy_type == 3 .and. pes%segment_inner(jj)%pend%ixy_type==3 )then
            ymid_t = (pes%segment_inner(jj)%porigin%coor(2)+pes%segment_inner(jj)%pend%coor(2) )/2.
            pes%segment_inner(jj)%id(2) = ceiling( (ymid_t-ybottom)/dy )
        elseif( pes%segment_inner(jj)%porigin%ixy_type == 1  )then
            pes%segment_inner(jj)%id(2) = pes%segment_inner(jj)%porigin%id_xy
        elseif( pes%segment_inner(jj)%porigin%ixy_type == 3 .and. pes%segment_inner(jj)%pend%ixy_type==1  )then
            pes%segment_inner(jj)%id(2) = pes%segment_inner(jj)%pend%id_xy
        endif
    enddo

    pes%nsub_inner = pes%nsub_inner + abs(nby)+1


    end subroutine super2final_inner_x
!--------------------------------------------------~   
subroutine super2final_inner_y(pes,superpoint1,superpoint2)
    implicit none
    type(type_element_upstream),pointer :: pes
    type(type_inter_point) :: superpoint1,superpoint2
    integer :: ix0,nbx
    integer :: ii,jj
    real :: xmid_t
    
    ix0 = superpoint1%id_xy
    nbx = ix0 - superpoint2%id_xy
 
    point_inner(1) = superpoint1

    if( nbx>0 )then
        do ii = 1 , nbx
            point_inner(ii+1)%coor(1) = xgrid(ix0+1-ii)
            point_inner(ii+1)%coor(2) = point_inner(1)%coor(2)
            point_inner(ii+1)%ixy_type = 3
            point_inner(ii+1)%idx2 =  2*(ix0+1-ii)-1
            point_inner(ii+1)%idy2 = point_inner(1)%idy2

        enddo
    elseif( nbx<0 )then
        do ii = 1 , abs(nbx)
            point_inner(ii+1)%coor(1) = xgrid( ix0+ii )
            point_inner(ii+1)%coor(2) = point_inner(1)%coor(2)
            point_inner(ii+1)%ixy_type = 3
            point_inner(ii+1)%idx2 =  2*(ix0+ii)-1
            point_inner(ii+1)%idy2 = point_inner(1)%idy2
        enddo
    else

    endif

    point_inner(2+abs(nbx) ) = superpoint2

    do ii = 1 , 1 + abs(nbx)
        jj = pes%nsub_inner + ii
        pes%segment_inner(jj)%porigin = point_inner(ii)
        pes%segment_inner(jj)%pend = point_inner(ii+1)

        if( pes%segment_inner(jj)%porigin%coor(1)<pes%segment_inner(jj)%pend%coor(1) )then
            ! from left to right
            pes%segment_inner(jj)%id(2)=(pes%segment_inner(jj)%porigin%idy2+1)/2
        else
            pes%segment_inner(jj)%id(2)=(pes%segment_inner(jj)%porigin%idy2-1)/2
        endif

        if(pes%segment_inner(jj)%porigin%ixy_type==3 .and. pes%segment_inner(jj)%pend%ixy_type==3 )then
            xmid_t = (pes%segment_inner(jj)%porigin%coor(1)+pes%segment_inner(jj)%pend%coor(1) )/2.
            pes%segment_inner(jj)%id(1) = ceiling( (xmid_t-xleft)/dx )
        elseif( pes%segment_inner(jj)%porigin%ixy_type==2 )then
            pes%segment_inner(jj)%id(1) = pes%segment_inner(jj)%porigin%id_xy
        elseif( pes%segment_inner(jj)%porigin%ixy_type==3 .and. pes%segment_inner(jj)%pend%ixy_type==2 )then
            pes%segment_inner(jj)%id(1) = pes%segment_inner(jj)%pend%id_xy
        endif

    enddo

    pes%nsub_inner = pes%nsub_inner + abs(nbx) + 1

    do ii = 1 , 1 + abs(nbx)
        jj = pes%nsub_inner + ii
        
        pes%segment_inner(jj)%porigin = point_inner(ii+1)
        pes%segment_inner(jj)%pend = point_inner(ii)


        if( pes%segment_inner(jj)%porigin%coor(1)<pes%segment_inner(jj)%pend%coor(1) )then
            ! from left to right
            pes%segment_inner(jj)%id(2)=(pes%segment_inner(jj)%porigin%idy2+1)/2
        else
            pes%segment_inner(jj)%id(2)=(pes%segment_inner(jj)%porigin%idy2-1)/2
        endif

        if(pes%segment_inner(jj)%porigin%ixy_type==3 .and. pes%segment_inner(jj)%pend%ixy_type==3 )then
            xmid_t = (pes%segment_inner(jj)%porigin%coor(1)+pes%segment_inner(jj)%pend%coor(1) )/2.
            pes%segment_inner(jj)%id(1) = ceiling( (xmid_t-xleft)/dx )
        elseif( pes%segment_inner(jj)%porigin%ixy_type==2 )then
            pes%segment_inner(jj)%id(1) = pes%segment_inner(jj)%porigin%id_xy
        elseif( pes%segment_inner(jj)%porigin%ixy_type==3 .and. pes%segment_inner(jj)%pend%ixy_type==2 )then
            pes%segment_inner(jj)%id(1) = pes%segment_inner(jj)%pend%id_xy
        endif

    enddo
    pes%nsub_inner = pes%nsub_inner + abs(nbx) + 1

    end subroutine super2final_inner_y
!--------------------------------------------------~   


!<4><get outersegments>
subroutine get_outersegments(p)
    implicit none

    type(type_face),pointer :: p
    type(type_segment),pointer :: pouter
    integer :: k
    real :: xmid_t,ymid_t

    do k = 1 , p%nsub_outer
        p%segment_outer(k)%porigin = p%point_inter(k-1)
        p%segment_outer(k)%pend = p%point_inter(k)
        
        pouter => p%segment_outer(k)
        if( pouter%porigin%ixy_type==0 )then
            pouter%id(1) = pouter%porigin%idx2/2
            pouter%id(2) = pouter%porigin%idy2/2
        elseif( pouter%porigin%ixy_type==1 .and. pouter%pend%ixy_type==0 )then
            pouter%id(1) = pouter%pend%idx2/2
            pouter%id(2) = pouter%pend%idy2/2            
        elseif( pouter%porigin%ixy_type==1 .and. pouter%pend%ixy_type==1 )then
            xmid_t = ( pouter%porigin%coor(1)+pouter%pend%coor(1) )/2.
            pouter%id(1) =  ceiling( (xmid_t - xleft)/dx )
            pouter%id(2) = pouter%porigin%idy2/2
        elseif( pouter%porigin%ixy_type==1 .and. pouter%pend%ixy_type==2 )then
            pouter%id(1) = pouter%pend%idx2/2
            pouter%id(2) = pouter%porigin%idy2/2
        elseif( pouter%porigin%ixy_type==2 .and. pouter%pend%ixy_type==0 )then
            pouter%id(1) = pouter%pend%idx2/2
            pouter%id(2) = pouter%pend%idy2/2             
        elseif( pouter%porigin%ixy_type==2 .and. pouter%pend%ixy_type==1 )then
            pouter%id(1) = pouter%porigin%idx2/2
            pouter%id(2) = pouter%pend%idy2/2
        elseif( pouter%porigin%ixy_type==2 .and. pouter%pend%ixy_type==2 )then
            pouter%id(1) = pouter%porigin%idx2/2
            ymid_t = ( pouter%porigin%coor(2)+pouter%pend%coor(2) )/2.
            pouter%id(2) = ceiling( (ymid_t - ybottom)/dy )
        endif
    enddo

    end subroutine get_outersegments
!--------------------------------------------------~   
    

!<5><a control subroutine>
subroutine get_intersections_outersegments
    implicit none    

    type(type_face),pointer :: p

    do i = 1 , nx
        do j = 1 , ny + 1
            p=>face_lr(i,j)

            if( iqc==0 )then
                call get_intersections(p)
            elseif( iqc==1 )then
                call get_subfaces(p)
                call get_intersections_qc(p)
            endif
            call get_outersegments(p)
        enddo
    enddo

    do i = 1 , nx +1
        do j = 1 , ny
            p=>face_bt(i,j)
            if( iqc==0 )then
                call get_intersections(p)
            elseif( iqc==1 )then
                call get_subfaces(p)
                call get_intersections_qc(p)
            endif
            call get_outersegments(p)
        enddo
    enddo
    
    end subroutine get_intersections_outersegments
!--------------------------------------------------~ 




    