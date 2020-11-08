    subroutine allocate_variable
    implicit none

    allocate( glx( 1:n_gl, 1:n_gl,  1:nel_x, 1:nel_y ) )
    allocate( gly( 1:n_gl, 1:n_gl,  1:nel_x, 1:nel_y ) )
    
    allocate( gx( 1:n_g, 1:n_g,  1:nel_x, 1:nel_y ) )
    allocate( gy( 1:n_g, 1:n_g,  1:nel_x, 1:nel_y ) )
    
    allocate( fel( 1:n_g, 1:n_g,  1:nel_x, 1:nel_y ) )
    
    allocate( x(1 - ighost:nel_x + ighost) )
    allocate( y(1 - ighost:nel_y + ighost) )
    
    allocate( xgrid(1-ighost :nel_x+1 + ighost ) )
    allocate( ygrid(1-ighost :nel_y+1 + ighost ) )
    
    ! Eulerian cells
    allocate( elem( 1 -ighost :nel_x+ighost, 1-ighost:nel_y+ighost ) )
    
    allocate( Dn_star_sub(1:3+int(cfl) ) )
    
    ! LDG
    allocate( ele_dg(1:nel_x*kdg) )
    allocate( ee_g(n_g,1:nel_x) )
    !*******************************************
    allocate( temp11(1:nel_x,1:nel_y,1:n_moment_2d) )
    !*******************************************

    end subroutine allocate_variable
    !************************************************
    subroutine deallocate_variable
    implicit none

    deallocate( glx )
    deallocate( gly )
    
    deallocate( gx )
    deallocate( gy )
    
    deallocate( fel )
    
    deallocate( x )
    deallocate( y )
    
    deallocate( xgrid )
    deallocate( ygrid )
    
    deallocate( elem )
    
    deallocate( Dn_star_sub )
    
    !LDG
    deallocate( ele_dg )
    deallocate( ee_g )
    
    deallocate( temp11 )

    end subroutine deallocate_variable