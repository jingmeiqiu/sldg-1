    subroutine SLDG1D_QiuGuoCai(    ixy,dxy, tnum,dt17,euler_gl, euler_gauss, xy_split,umod_temp )
    implicit none
    integer,intent(in) :: ixy
    real,intent(in) :: dxy,dt17,tnum
    real,intent(in) :: euler_gl(1:n_gl)
    real,intent(in) :: euler_gauss(1:n_g)
    real,intent(in) :: xy_split
    real :: upstream_gl(1:n_gl), upstream_gauss(1:n_g)
    real,intent(out) :: umod_temp(1:n_moment )
    real :: xl_star,xr_star
    integer :: ia,ib
    integer :: mx
    integer :: kc,inter
    integer :: isub , nsub
    integer :: idx
    real,allocatable :: zz(:)
    
    integer :: igg

    ! STEP1. Locate the foots of trajectory x_l^{n,star}, x_r^{n,star}.
    ! We numerically solve the following final-value problem (trajectory equation):
    !             d x(t) / dt = a( x(t) , t )
    ! with the final-value x( t^{n+1} ) = x_l^n, x_r^n by a high order numerical integrator such as a classical fourth-order
    ! numerical integrator such as a classical fourth-order Runge-Kutta method.


    ! call RK(vx(1:nx+1) ,vx_star(1:nx+1),nx+1,dt )

    ! STEP2. locate Gauss-lobatto points


    if( ixy == 1 )then
        call RK_x( xy_split, euler_gl(1:n_gl), upstream_gl(1:n_gl) ,n_gl , tnum,dt17 )
        call RK_x( xy_split, euler_gauss(1:n_g), upstream_gauss(1:n_g) ,n_g ,tnum,dt17 )
    elseif(ixy == 2)then
        do igg = 1 , n_gl
            upstream_gl( igg ) = euler_gl( igg ) - ee_g( kk, ie ) *dt17
        enddo
        do igg = 1 , n_g 
            upstream_gauss( igg ) = euler_gauss( igg ) - ee_g( kk, ie ) *dt17
        enddo        
         
    endif


    !
    xl_star =  upstream_gl( 1 )
    xr_star =  upstream_gl( n_gl )


    call id_get17( xl_star,dxy,ia,ixy )

    !   if(ixy ==2 .and. je == 10)then
    !!print *,euler_gl( 1:n_gl ), upstream_gl(1:n_gl)
    !       print *,xr_star
    !pause
    !   endif

    call id_get17( xr_star,dxy,ib,ixy )



    !ia = id_get(  (Dn_star(k)%xl_star - xleft)/dx );  ib = id_get(  (Dn_star(k)%xr_star-xleft)/dx  )

    mx = ib - ia

    allocate( zz(1:2+mx) )

    zz(1) = xl_star
    zz(2+mx ) = xr_star

    if(mx .ne. 0)then
        do kc = 1 , mx
            inter = ia + kc
            if( ixy == 1 )then
                zz( 1+kc ) = xgrid( inter )
            elseif(ixy == 2)then
                zz( 1+kc ) = ygrid(inter)
            endif
        enddo
    endif

    do isub = 1 , 1 + mx
        Dn_star_sub(isub)%vpoint(1) = zz(isub)
        Dn_star_sub(isub)%vpoint(2) = zz(isub+1)


        call id_get17( (zz(isub)+zz(isub+1) )/2.,dxy,idx ,ixy )
        !idx = id_get( ( (zz(ii)+zz(ii+1) )/2. -xleft)/dx )
        Dn_star_sub(isub)%id =  idx

    enddo

    deallocate(zz)
    nsub = 1 + mx

    if(nod==1)then
        call Green2( upstream_gl(1:n_gl ), upstream_gauss(1:n_g),umod_temp(1:n_moment ) , nsub, dxy,ixy)

    elseif(nod==2)then
        !call Green3( umod_temp(1:n_gl ) )
        call Green3( upstream_gl(1:n_gl ), upstream_gauss(1:n_g), umod_temp(1:n_moment ) , nsub, dxy,ixy)
    elseif(nod==3)then
        !call Green4( umod_temp(1:n_gl ) )
    endif


    end subroutine SLDG1D_QiuGuoCai