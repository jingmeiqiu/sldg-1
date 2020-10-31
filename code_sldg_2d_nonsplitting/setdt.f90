    subroutine setdt
    implicit none

    if(iexample == 1)then
        dt = cfl/( 1./dx + 1./dy )
    elseif( iexample == 2 )then
        dt = cfl/( pi/dx + pi/dy )
    elseif( iexample == 3 )then
        dt = cfl/( 2*pi/dx + 2*pi/dy )
    endif

    !dt = 2.5*dx

    if(time+dt>time_final) dt = time_final- time

    time = time +dt
    


    end subroutine setdt