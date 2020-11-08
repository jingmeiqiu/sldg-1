    real function fle(k,x)
    implicit none
    integer, intent(in) :: k
    real,intent(in) :: x
    ! 关于Legendre正交多项式的函数
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
    !***************************************************************************
    real function polynomial(a00,x16,xc16,dx16,k)
    implicit none
    real,intent(in) :: x16,xc16,dx16
    integer,intent(in) :: k
    real,intent(in) :: a00(k+1)


    if(k==0)then
        polynomial = a00(1)
    elseif(k==1)then
        polynomial = a00(1) + a00(2)*(x16 - xc16 )/dx16
    elseif(k==2)then
        polynomial = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(  ((x16 - xc16) /dx16)**2 - 1./12. )
    elseif(k==3)then
        polynomial = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(  ((x16 - xc16) /dx16)**2 - 1./12. ) &
            +a00(4) *( ((x16 - xc16) /dx16)**3 - 3./20.*((x16 - xc16) /dx16) )
    endif

    end function polynomial
    
    real function polynomial_1d(a00,x16,xc16,dx16,k)
    implicit none
    real,intent(in) :: x16,xc16,dx16
    integer,intent(in) :: k
    real,intent(in) :: a00(k+1)


    if(k==0)then
        polynomial_1d = a00(1)
    elseif(k==1)then
        polynomial_1d = a00(1) + a00(2)*(x16 - xc16 )/dx16
    elseif(k==2)then
        polynomial_1d = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(  ((x16 - xc16) /dx16)**2 - 1./12. )
    elseif(k==3)then
        polynomial_1d = a00(1) + a00(2)*(x16 - xc16 )/dx16 +  a00(3)*(  ((x16 - xc16) /dx16)**2 - 1./12. ) &
        +a00(4) *( ((x16 - xc16) /dx16)**3 - 3./20.*((x16 - xc16) /dx16) )
    endif

    end function polynomial_1d
    
    !%
    real function fphi(k,x,y)
    implicit none
    integer,intent(in) :: k
    real,intent(in) :: x,y

    if(k .eq. 1)then
        fphi = 1
    elseif(k .eq. 2) then
        fphi = x
    elseif(k .eq. 3)then
        fphi = y
    elseif(k .eq. 4)then
        fphi = x**2 - 1./12.
    elseif(k .eq. 5)then
        fphi = x*y
    elseif(k .eq. 6)then
        fphi = y**2 - 1./12.
    endif


    end function fphi