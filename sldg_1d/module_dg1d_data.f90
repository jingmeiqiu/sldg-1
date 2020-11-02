    module module_dg1d_data

    implicit none
    ! data structure
    !*************************************************************************
    type, public :: point1d
        sequence
        real :: coor
    end type

    type(point1d),allocatable,target,public :: vertex(:)
    !*******************************

    !*************************************************************************
    type,public :: element1d_Eulerian
        sequence
        real :: umodal(6)
        real :: x(1:2)

        real :: xg(1:5)

    end type
    type(element1d_Eulerian),allocatable,target,public :: element(:)
    !*******************************

    end module module_dg1d_data