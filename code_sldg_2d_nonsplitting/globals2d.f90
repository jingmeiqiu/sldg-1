    module globals2d

    integer :: i,j
    integer :: nx,ny
    real :: pi,eps
    real :: xleft,xright,ybottom,ytop
    real :: dx,dy
    real,allocatable :: xgrid(:),ygrid(:)
    real,allocatable :: x(:),y(:)
    integer :: nghost
    integer :: iexample,irk
    integer :: n_moment,iqc
    integer :: nt
    integer :: irelative
    real :: r_norm1
    
    real :: time_final,time,dt,cfl
    
    real,allocatable :: umod_t(:,:,:)
    
    !********************************************
    integer :: iam,ibm,icm,idm
    integer :: isx(6),isy(6),ix,iy
    !********************************************
    
    !********************************************
    ! parameters
    real :: gau2(1:2,1:2)
    real :: gau3(1:3,1:3)
    real :: ai(1:6)
    real :: xg(1:6),wg(1:6)
    !********************************************
    
    !********************************************
    ! order
    real :: er11,er22,er33
    integer :: kkkk
    !********************************************
    
    integer :: idebug
    
    real,allocatable :: com_mass(:,:)
    
    end module globals2d