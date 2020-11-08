    subroutine LDG_parameters
    implicit none
    real :: temp,temp1,temp2


    vm(0) = 1.
    vm(1) = 0.5
    vm(2) = 1./6.
    vm(3) = 1./20.

    vp(0) = 1.
    vp(1) = -0.5
    vp(2) = 1./6.
    vp(3) = -1./20.
    
 
    

    !   guassian weight and point for basis
    if(kdg .eq. 1)then
        !
        wq(1) = 1.
        ! v0 4 orders of gaussian qudrature rule( 2 points)
        vxgau(1,0) = 0.  ! l denotes the gaussian point, k denotes the number of basis

        vgau(1,0) = 1.
        
        vgau(1,1) = 0.

        ain(1,1) = 1.
        BL(1,1) = -1.
    elseif(kdg.eq.2) then
        !
        !cfl1 = 1./4.

        wq(1) = 0.5
        wq(2) = 0.5
        ! v0 4 orders of gaussian qudrature rule( 2 points)
        vxgau(1,0) = 0.  ! l denotes the gaussian point, k denotes the number of basis
        vxgau(2,0) = 0.
        vgau(1,0) = 1.
        vgau(2,0) = 1.
        ! v1
        vxgau(1,1) = 1.
        vxgau(2,1) = 1.
        vgau(1,1) = -sqrt(3.)/6.
        vgau(2,1) = -vgau(1,1) 

        ain(1,1) = 0.5; ain(1,2)=-1.; 
        ain(2,1) = 1.; ain(2,2) = 2.; 


        BL(1,1) = -1.; BL(1,2)= -0.5; 
        BL(2,1) = 0.5; BL(2,2) = 0.25 ; 


    else if(kdg.eq.3) then
        !
        ! cfl1 = 1./6.

        wq(1) = 5./18.
        wq(2) = 4./9.
        wq(3) = wq(1)
        ! v0 6 orders of gaussian qudrature rule( 3 points)
        vxgau(1,0) = 0.
        vxgau(2,0) = 0.
        vxgau(3,0) = 0.
        vgau(1,0) = 1.
        vgau(2,0) = 1.
        vgau(3,0) = 1.
        ! v1
        temp = sqrt(15.)/10.
        vxgau(1,1) = 1.
        vxgau(2,1) = 1.
        vxgau(3,1) = 1.
        vgau(1,1) = -temp
        vgau(2,1) = 0.
        vgau(3,1) = temp
        ! v2
        vxgau(1,2) = -2.*temp
        vxgau(2,2) = 0.
        vxgau(3,2) = 2.*temp
        vgau(1,2) = temp**2. - 1./12.
        vgau(2,2) = - 1./12.
        vgau(3,2) = vgau(1,2)

        ain(1,1) = 0.5; ain(1,2)=-1.; ain(1,3) = 0.; 
        ain(2,1) = 1.; ain(2,2) = 0.; ain(2,3) = -6.; 
        ain(3,1) = 0.; ain(3,2) = 6.; ain(3,3) = 18.;


        BL(1,1) = -1.; BL(1,2)= -0.5; BL(1,3) = -1./6.;
        BL(2,1) = 0.5; BL(2,2) = 0.25 ; BL(2,3) = 1./12.;
        BL(3,1) = -1./6.; BL(3,2) = -1./12.; BL(3,3) = -1./36.;


    else if(kdg.eq.4) then
        !cfl1 = 1./8.
        wq(1) = 0.25 - sqrt(30.)/72.
        wq(2) = 0.25 + sqrt(30.)/72.
        wq(3) = wq(2)
        wq(4) = wq(1)

        ! v0 8 orders of gaussian qudrature rule( 4 points)
        vxgau(1,0) = 0.
        vxgau(2,0) = 0.
        vxgau(3,0) = 0.
        vxgau(4,0) = 0.
        vgau(1,0) = 1.
        vgau(2,0) = 1.
        vgau(3,0) = 1.
        vgau(4,0) = 1.
        ! v1
        temp1 = sqrt(3./7.+2./35.*sqrt(30.))/2.
        temp2 = sqrt(3./7.-2./35.*sqrt(30.))/2.
        vxgau(1,1) = 1.
        vxgau(2,1) = 1.
        vxgau(3,1) = 1.
        vxgau(4,1) = 1.
        vgau(1,1) = -temp1
        vgau(2,1) = -temp2
        vgau(3,1) = temp2
        vgau(4,1) = temp1
        ! v2
        vxgau(1,2) = 2.*vgau(1,1)
        vxgau(2,2) = 2.*vgau(2,1)
        vxgau(3,2) = 2.*vgau(3,1)
        vxgau(4,2) = 2.*vgau(4,1)
        vgau(1,2) = temp1**2.-1./12.
        vgau(2,2) = temp2**2.-1./12.
        vgau(3,2) = vgau(2,2)
        vgau(4,2) = vgau(1,2)
        ! v3

        vxgau(1,3) = 3.*temp1**2.-3./20.
        vxgau(2,3) = 3.*temp2**2.-3./20.
        vxgau(3,3) = vxgau(2,3)
        vxgau(4,3) = vxgau(1,3)
        vgau(1,3) = -temp1**3.+3./20.*temp1
        vgau(2,3) = -temp2**3.+3./20.*temp2
        vgau(3,3) = -vgau(2,3)
        vgau(4,3) = -vgau(1,3)

        ain(1,1) = 0.5; ain(1,2)=-1.; ain(1,3) = 0.; ain(1,4) = 0.;
        ain(2,1) = 1.; ain(2,2) = 0.; ain(2,3) = -6.; ain(2,4) = 0.;
        ain(3,1) = 0.; ain(3,2) = 6.; ain(3,3) = 0.; ain(3,4) = -60.;
        ain(4,1) = 0.; ain(4,2) = 0.; ain(4,3) = 60.; ain(4,4) = 200.

        BL(1,1) = -1.; BL(1,2)= -0.5; BL(1,3) = -1./6.; BL(1,4) = -1./20.;
        BL(2,1) = 0.5; BL(2,2) = 0.25 ; BL(2,3) = 1./12.; BL(2,4) = 1./40.;
        BL(3,1) = -1./6.; BL(3,2) = -1./12.; BL(3,3) = -1./36.; BL(3,4) = -1./120.;
        BL(4,1) = 1./20.; BL(4,2) = 1./40.; BL(4,3) = 1./120.; BL(4,4) = 1./400.     

    else
        !
        stop

    endif




    end subroutine LDG_parameters