    subroutine setdt(time_f)
    use module_prob, only : cfl
    use module_dg1d, only : dx
    implicit none
    real,intent(in) :: time_f
    ! set dt and get tn+1
    
    dt = cfl*dx
    if(time+dt>time_f) dt = time_f-time

    time = time + dt

    print *,time

    end subroutine setdt