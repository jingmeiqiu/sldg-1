    subroutine output
    implicit none
    integer :: ic,jc
    real :: xrg,yrg

    !open(121,file='temp_16.plt')
    !write(121,*)'zone    ','i=',nel_x*n_gl,',    j=',nel_y*n_gl
    !
    !DO  J=1,Nel_Y
    !    do jc = 1,n_gl
    !        DO  I=1 ,Nel_X
    !            do ic = 1,n_gl
    !                xrg =  x(i)+ x_gl(ic)*dx
    !                yrg =  y(j)+ x_gl(jc)*dy
    !                WRITE(121,*) xrg,yrg,  elem(i,j)%psi( ic,jc )
    !            enddo
    !        enddo
    !    enddo
    !enddo
    !
    !
    !CLOSE(121)

    open(121,file='temp_16.plt')
    write(121,*)'zone    ','i=',nel_x*n_g,',    j=',nel_y*n_g

    DO  J=1,Nel_Y
        do jc = 1,n_g
            DO  I=1 ,Nel_X
                do ic = 1,n_g
                    xrg =  x(i)+ x_g(ic)*dx
                    yrg =  y(j)+ x_g(jc)*dy
                    WRITE(121,*) xrg,yrg,  elem(i,j)%psi( ic,jc )
                enddo
            enddo
        enddo
    enddo


    CLOSE(121)
    end subroutine output