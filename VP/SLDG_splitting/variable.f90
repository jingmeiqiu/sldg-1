    module variable
    implicit none

    !
    real,parameter :: pi = 4.*atan(1.)
    !
    
    integer :: kkkk
    integer :: nod,n_moment
    !integer,parameter :: nod = 2
    integer :: i,j,k,kk
    integer :: ig,jg
    integer :: ie,je
    integer :: i_case
    integer :: irk
    ! nod is the degree of the polynomial
    integer :: n_gl,nel_x,nel_y,ighost
    integer :: nel
    real :: dx,dy
    real :: xleft,xright,yleft,yright
    
    real :: tprint,tn,dt,cfl
    integer :: nt
    
    real :: x_GL(5),w_GL(5)
    real :: x_g(5),w_g(5)
    integer :: n_g
    !real :: gauss(4,1)
    real :: ai(5),ai2d(6)
    real, allocatable :: glx(:,:,:,:),gly(:,:,:,:),fel(:,:,:,:)
  real, allocatable :: gx(:,:,:,:),gy(:,:,:,:)
    real, allocatable :: x(:),y(:)
    real, allocatable :: xgrid(:),ygrid(:)

    ! test order
    real :: er11,er22,er33
    
    !**************variable for LDG
    integer :: n_moment_2d
    real :: pint
    real :: vm(0:3),vp(0:3)
    integer :: kdg
    ! Gauss
    real :: wq(4),vgau(4,0:3), vxgau(4,0:3)
    real :: ain( 4,4 ),BL(4,4)
    real,allocatable :: ele_dg(:),ee_g(:,:)
    real,allocatable :: temp11(:,:,:)
    real :: emax
    
    !******************************************
    integer :: norder(10)
    real :: begin_time,end_time
    
    !**************************************************
    real :: xg(6),wg(6)
    real :: r_entropy
    real :: r_energy
    real :: r_norm1
    real :: r_enstrophy
    real :: r_norm2
    integer :: irelative
    
    end module variable