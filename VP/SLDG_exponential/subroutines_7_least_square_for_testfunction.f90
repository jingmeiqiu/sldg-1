!*******************************************************************
! would you please help me to polish the following subroutines.
!*******************************************************************
!<1><Least square for test functions>
subroutine get_matrix_a( ver,at )
    implicit none
    real,intent(in) :: ver(1:4,1:2)
    real,intent(out) :: at(1:3,1:3)

    at(1,1) = 4.;
    at(1,2) = ver(1,1)+ver(2,1)+ver(3,1)+ver(4,1); at(1,3)=ver(1,2)+ver(2,2)+ver(3,2)+ver(4,2);
    at(2,1) = at(1,2);
    at(2,2) = ver(1,1)*ver(1,1)+ver(2,1)*ver(2,1)+ver(3,1)*ver(3,1)+ver(4,1)*ver(4,1);
    at(2,3) = ver(1,1)*ver(1,2)+ver(2,1)*ver(2,2)+ver(3,1)*ver(3,2)+ver(4,1)*ver(4,2);
    at(3,1) = at(1,3);
    at(3,2) = at(2,3);
    at(3,3) = ver(1,2)*ver(1,2)+ver(2,2)*ver(2,2)+ver(3,2)*ver(3,2)+ver(4,2)*ver(4,2);


end subroutine  get_matrix_a
!--------------------------------------------------~  
subroutine get_vector_b( ver,ps,bt )
    implicit none
    real,intent(in) :: ver(1:4,1:2),ps(1:4)
    real,intent(out) ::  bt(1:3)


    bt(1) = ps(1) + ps(2)+ps(3)+ps(4);
    bt(2) = ps(1)*ver(1,1) + ps(2)*ver(2,1)+ps(3)*ver(3,1)+ps(4)*ver(4,1);
    bt(3) = ps(1)*ver(1,2) + ps(2)*ver(2,2)+ps(3)*ver(3,2)+ps(4)*ver(4,2);

    end subroutine  get_vector_b
!--------------------------------------------------~  
subroutine get_matrix2_a( ver, at  )

    implicit none
    real,intent(in) :: ver(1:9,1:5)
    real,intent(out) :: at(1:6,1:6)

    at(1, 1) = 9
    at(1, 2) = ver(1,1)+ver(2,1)+ver(3,1)+ver(4,1)+ver(5,1)+ver(6,1)+ver(7,1)+ver(8,1)+ver(9,1)
    at(1, 3) = ver(1,2)+ver(2,2)+ver(3,2)+ver(4,2)+ver(5,2)+ver(6,2)+ver(7,2)+ver(8,2)+ver(9,2)
    at(1, 4) = ver(1,3)+ver(2,3)+ver(3,3)+ver(4,3)+ver(5,3)+ver(6,3)+ver(7,3)+ver(8,3)+ver(9,3)
    at(1, 5) = ver(1,4)+ver(2,4)+ver(3,4)+ver(4,4)+ver(5,4)+ver(6,4)+ver(7,4)+ver(8,4)+ver(9,4)
    at(1, 6) = ver(1,5)+ver(2,5)+ver(3,5)+ver(4,5)+ver(5,5)+ver(6,5)+ver(7,5)+ver(8,5)+ver(9,5)
    at(2, 1) = at(1,2)
    at(2, 2) = ver(1,1)*ver(1,1)+ver(2,1)*ver(2,1)+ver(3,1)*ver(3,1)+ver(4,1)*ver(4,1) &
        +ver(5,1)*ver(5,1)+ver(6,1)*ver(6,1)+ver(7,1)*ver(7,1)+ver(8,1)*ver(8,1)+ver(9,1)*ver(9,1)
    at(2, 3) = ver(1,1)*ver(1,2)+ver(2,1)*ver(2,2)+ver(3,1)*ver(3,2)+ver(4,1)*ver(4,2) &
        +ver(5,1)*ver(5,2)+ver(6,1)*ver(6,2)+ver(7,1)*ver(7,2)+ver(8,1)*ver(8,2)+ver(9,1)*ver(9,2)
    at(2, 4) = ver(1,1)*ver(1,3)+ver(2,1)*ver(2,3)+ver(3,1)*ver(3,3)+ver(4,1)*ver(4,3) &
        +ver(5,1)*ver(5,3)+ver(6,1)*ver(6,3)+ver(7,1)*ver(7,3)+ver(8,1)*ver(8,3)+ver(9,1)*ver(9,3)
    at(2, 5) = ver(1,1)*ver(1,4)+ver(2,1)*ver(2,4)+ver(3,1)*ver(3,4)+ver(4,1)*ver(4,4) &
        +ver(5,1)*ver(5,4)+ver(6,1)*ver(6,4)+ver(7,1)*ver(7,4)+ver(8,1)*ver(8,4)+ver(9,1)*ver(9,4)
    at(2, 6) = ver(1,1)*ver(1,5)+ver(2,1)*ver(2,5)+ver(3,1)*ver(3,5)+ver(4,1)*ver(4,5) &
        +ver(5,1)*ver(5,5)+ver(6,1)*ver(6,5)+ver(7,1)*ver(7,5)+ver(8,1)*ver(8,5)+ver(9,1)*ver(9,5)
    at(3, 1) = at(1,3)
    at(3, 2) = at(2,3)
    at(3, 3) = ver(1,2)*ver(1,2)+ver(2,2)*ver(2,2)+ver(3,2)*ver(3,2)+ver(4,2)*ver(4,2) &
        +ver(5,2)*ver(5,2)+ver(6,2)*ver(6,2)+ver(7,2)*ver(7,2)+ver(8,2)*ver(8,2)+ver(9,2)*ver(9,2)
    at(3, 4) = ver(1,2)*ver(1,3)+ver(2,2)*ver(2,3)+ver(3,2)*ver(3,3)+ver(4,2)*ver(4,3) &
        +ver(5,2)*ver(5,3)+ver(6,2)*ver(6,3)+ver(7,2)*ver(7,3)+ver(8,2)*ver(8,3)+ver(9,2)*ver(9,3)
    at(3, 5) = ver(1,2)*ver(1,4)+ver(2,2)*ver(2,4)+ver(3,2)*ver(3,4)+ver(4,2)*ver(4,4) &
        +ver(5,2)*ver(5,4)+ver(6,2)*ver(6,4)+ver(7,2)*ver(7,4)+ver(8,2)*ver(8,4)+ver(9,2)*ver(9,4)
    at(3, 6) = ver(1,2)*ver(1,5)+ver(2,2)*ver(2,5)+ver(3,2)*ver(3,5)+ver(4,2)*ver(4,5) &
        +ver(5,2)*ver(5,5)+ver(6,2)*ver(6,5)+ver(7,2)*ver(7,5)+ver(8,2)*ver(8,5)+ver(9,2)*ver(9,5)
    at(4, 1) = at(1,4)
    at(4, 2) = at(2,4)
    at(4, 3) = at(3,4)
    at(4, 4) = ver(1,3)*ver(1,3)+ver(2,3)*ver(2,3)+ver(3,3)*ver(3,3)+ver(4,3)*ver(4,3) &
        +ver(5,3)*ver(5,3)+ver(6,3)*ver(6,3)+ver(7,3)*ver(7,3)+ver(8,3)*ver(8,3)+ver(9,3)*ver(9,3)
    at(4, 5) = ver(1,3)*ver(1,4)+ver(2,3)*ver(2,4)+ver(3,3)*ver(3,4)+ver(4,3)*ver(4,4) &
        +ver(5,3)*ver(5,4)+ver(6,3)*ver(6,4)+ver(7,3)*ver(7,4)+ver(8,3)*ver(8,4)+ver(9,3)*ver(9,4)
    at(4, 6) = ver(1,3)*ver(1,5)+ver(2,3)*ver(2,5)+ver(3,3)*ver(3,5)+ver(4,3)*ver(4,5) &
        +ver(5,3)*ver(5,5)+ver(6,3)*ver(6,5)+ver(7,3)*ver(7,5)+ver(8,3)*ver(8,5)+ver(9,3)*ver(9,5)
    at(5, 1) = at(1,5)
    at(5, 2) = at(2,5)
    at(5, 3) = at(3,5)
    at(5, 4) = at(4,5)
    at(5, 5) = ver(1,4)*ver(1,4)+ver(2,4)*ver(2,4)+ver(3,4)*ver(3,4)+ver(4,4)*ver(4,4) &
        +ver(5,4)*ver(5,4)+ver(6,4)*ver(6,4)+ver(7,4)*ver(7,4)+ver(8,4)*ver(8,4)+ver(9,4)*ver(9,4)
    at(5, 6) = ver(1,4)*ver(1,5)+ver(2,4)*ver(2,5)+ver(3,4)*ver(3,5)+ver(4,4)*ver(4,5) &
        +ver(5,4)*ver(5,5)+ver(6,4)*ver(6,5)+ver(7,4)*ver(7,5)+ver(8,4)*ver(8,5)+ver(9,4)*ver(9,5)
    at(6, 1) = at(1,6)
    at(6, 2) = at(2,6)
    at(6, 3) = at(3,6)
    at(6, 4) = at(4,6)
    at(6, 5) = at(5,6)
    at(6, 6) = ver(1,5)*ver(1,5)+ver(2,5)*ver(2,5)+ver(3,5)*ver(3,5)+ver(4,5)*ver(4,5) &
        +ver(5,5)*ver(5,5)+ver(6,5)*ver(6,5)+ver(7,5)*ver(7,5)+ver(8,5)*ver(8,5)+ver(9,5)*ver(9,5)

    end subroutine  get_matrix2_a
!--------------------------------------------------~  
subroutine get_vector2_b( ver,ps,bt )

    implicit none
    real,intent(in) :: ver(1:9,1:5),ps(1:9)
    real,intent(out) :: bt(1:6)

    bt(1) = ps(1)+ps(2)+ps(3)+ps(4)+ps(5)+ps(6)+ps(7)+ps(8)+ps(9)
    bt(2) = ver(1,1)*ps(1)+ver(2,1)*ps(2)+ver(3,1)*ps(3)+ver(4,1)*ps(4)+ver(5,1)*ps(5)+ver(6,1)*ps(6)+ver(7,1)*ps(7)+ver(8,1)*ps(8)+ver(9,1)*ps(9)
    bt(3) = ver(1,2)*ps(1)+ver(2,2)*ps(2)+ver(3,2)*ps(3)+ver(4,2)*ps(4)+ver(5,2)*ps(5)+ver(6,2)*ps(6)+ver(7,2)*ps(7)+ver(8,2)*ps(8)+ver(9,2)*ps(9)
    bt(4) = ps(1)*ver(1,3)+ps(2)*ver(2,3)+ps(3)*ver(3,3)+ps(4)*ver(4,3)+ps(5)*ver(5,3)+ps(6)*ver(6,3)+ps(7)*ver(7,3)+ps(8)*ver(8,3)+ps(9)*ver(9,3)
    bt(5) = ps(1)*ver(1,4)+ps(2)*ver(2,4)+ps(3)*ver(3,4)+ps(4)*ver(4,4)+ps(5)*ver(5,4)+ps(6)*ver(6,4)+ps(7)*ver(7,4)+ps(8)*ver(8,4)+ps(9)*ver(9,4)
    bt(6) = ps(1)*ver(1,5)+ps(2)*ver(2,5)+ps(3)*ver(3,5)+ps(4)*ver(4,5)+ps(5)*ver(5,5)+ps(6)*ver(6,5)+ps(7)*ver(7,5)+ps(8)*ver(8,5)+ps(9)*ver(9,5)

    end subroutine  get_vector2_b
!--------------------------------------------------~  
subroutine solve(A,b,xx,N)
    implicit none
    integer,intent(in) :: N
    real,intent(in) :: A(N,N),b(N)
    real,intent(out) :: xx(N)

    real :: L(N,N),U(N,N)

    real :: yy(N)

    call doolittle(A,L,U,N)

    call  downtri(L,b,yy,N)

    call uptri(U,yy,xx,N)

    end subroutine solve
!--------------------------------------------------~  
subroutine solve17(L,U,b,xx,N)
    implicit none
    integer,intent(in) :: N
    real,intent(in) :: L(N,N),U(N,N),b(N)
    real,intent(out) :: xx(N)

    real :: yy(N)

    call  downtri(L,b,yy,N)

    call uptri(U,yy,xx,N)

    end subroutine solve17
!--------------------------------------------------~  
subroutine doolittle(A,L,U,N)
    !-------------------------------------------subroutine  comment
    !  Purpose   :  Doolittle function
    !              A=LU
    !--------------------------------------------------------------
    !  Input  parameters  :
    !       1.    A  n*n matrix
    !       2.    N  order
    !  Output parameters  :
    !       1.   L
    !       2.   U
    !---------------------------------------------------------------
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: A(N,N)
    real,intent(out) ::  L(N,N),U(N,N)
    integer :: k,i,j,m
    real :: s

    ! the first row of U

    U(1,:)=A(1,:)

    ! the first column of L
    L(:,1)=a(:,1)/U(1,1)

    do k=2,N

        l(k,k)=1

        do j=k,n
            s=0
            do m=1,k-1
                s=s+l(k,m)*u(m,j)
            end do
            u(k,j)=a(k,j)-s
        end do


        do i=k+1,n
            s=0
            do m=1,k-1
                s=s+l(i,m)*u(m,k)
            end do
            l(i,k)=(a(i,k)-s)/u(k,k)

        end do


    end do

    end subroutine doolittle
!--------------------------------------------------~  
subroutine uptri(A,b,xx,N)
    !--------------------------------------------subroutine  comment
    !  Purpose   :  backward substitution of upper triangular matrix
    !                 Ax=b
    !---------------------------------------------------------------
    !  Input  parameters  :
    !       1.   A(N,N)
    !       2.   b(N)
    !       3.   N
    !  Output parameters  :
    !       1.  xx
    !
    !---------------------------------------------------------------
    implicit none

    integer,intent(in) :: n

    real,intent(in) :: A(N,N),b(N)
    real,intent(out) :: xx(N)
    integer :: i,j

    xx(N)=b(N)/A(N,N)

    ! backward substitution
    do i=n-1,1,-1

        xx(i)=b(i)
        do j=i+1,N
            xx(i)=xx(i)-a(i,j)*xx(j)
        end do
        xx(i)=xx(i)/A(i,i)

    end do

    end subroutine uptri
!--------------------------------------------------~  
subroutine downtri(A,b,xx,N)
    !--------------------------------------------subroutine  comment
    !  Purpose   :  backward substitution of lower triangular matrix
    !                 Ax=b
    !----------------------------------------------------------------
    !  Input  parameters  :
    !       1.   A(N,N)
    !       2.   b(N)
    !       3.   N
    !  Output parameters  :
    !       1.  xx
    !
    !---------------------------------------------------------------
    implicit none
    integer,intent(in) :: N
    real,intent(in) :: A(N,N),b(N)
    real,intent(out) :: xx(N)
    integer :: k,i,j

    xx(1)=b(1)/a(1,1)

    do k=2,N
        xx(k)=b(k)
        do i=1,k-1
            xx(k)=xx(k)-a(k,i)*xx(i)
        end do
        xx(k)=xx(k)/a(k,k)

    end do

    end subroutine downtri

