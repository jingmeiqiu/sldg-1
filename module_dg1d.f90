    module module_dg1d

    use module_dg1d_data

    use module_prob, only : nx
    use module_prob, only : nk
    use module_prob, only : nghost
    implicit none

    !*********
    ! mesh
    real,allocatable :: x(:)   !x(i) ----> center
    real,allocatable :: xgrid(:)

    real :: dx
    !*********

    REAL, PARAMETER, DIMENSION(1:5):: ai=(/ 1., 12.,180.,2800., 44100. /)

    contains
    !*************************************************************************
    subroutine allocate_var
    implicit none
    !*******************************
    !allocate_variables
    allocate( vertex(1:nx+1) )

    allocate( x(1-nghost:nx+nghost) )
    allocate( xgrid(1-nghost:nx+1+nghost) )

    allocate( element(1-nghost:nx+nghost) )
    !end allocate_variables
    !*******************************
    end subroutine allocate_var


    !*********************************
    subroutine deallocate_var
    implicit none
    !*******************************
    !deallocate_variables
    deallocate( vertex )

    deallocate( x )

    deallocate( xgrid )

    deallocate( element )
    !end deallocate_variables
    !*******************************
    end subroutine deallocate_var

    !*************************************************************************
    subroutine init

    use module_prob, only : xleft,xright

    use module_prob, only : gau
    use module_prob, only : fun_init

    implicit none
    integer :: i,k

    real :: u
    integer :: lx

    ! x(i) ----> center

    dx = (xright-xleft)/nx

    DO I = 1 - nghost, nx + nghost
        X(I) = XLEFT + (I-0.5) * DX
    ENDDO

    DO I = 1 - nghost, nx + 1  + nghost
        Xgrid(I) = XLEFT + (I-1.) * DX
    ENDDO

    do i = 1,nx+1
        vertex(i)%coor = xleft + (i-1) * dx
    enddo

    do i = 1,nx
        element(i)%x(1) = x(i)-0.5*dx
        element(i)%x(2) = x(i)+0.5*dx
        !****************************
        do k = 1,nk+1
            element(i)%xg(k) = x(i) + gau(k,1)*dx
        enddo
    enddo


    do i = 1 , nx
        do k = 0,nk
            u = 0.0
            do lx = 1 , nk+1
                u=u+fun_init(x(i)+gau(lx,1)*dx )*fle(k,gau(lx,1))*gau(lx,2)
            enddo
            element(i)%umodal(k+1) = u*ai(k+1)
        enddo
    enddo

    end subroutine init
    !*************************************************************************
    !********************************
    real function fle(k,x)
    implicit none
    integer, intent(in) :: k
    real,intent(in) :: x
    ! Legendre orthogonal polynomial
    ! v_0=1.0, v_1=(x-x_j)/dx_j, v_2=( (x-x_j)/dx_j )^2-1.0/12.0,
    ! v_3=( (x-x_j)/dx_j )^3- 0.15*(x-x_j)/dx_j
    ! v_4=( ((x-x_j)/dx_j)^2-3.0/14.0 )*((x-x_j)/dx_j)^2 + 3.0/560.0

    if(k.eq.0) then
        fle=1.0
    elseif (k.eq.1) then
        fle= x
    elseif(k.eq.2) then
        fle= x*x-1./12.0
    elseif (k.eq.3) then
        fle=x*x*x-0.15*x
    elseif(k.eq.4) then
        fle=(x**2-3./14.0)*x*x+3.0/560.0
    endif

    return
    end function fle
    end module module_dg1d