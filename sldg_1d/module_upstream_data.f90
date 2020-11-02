    module module_upstream_data

    implicit none

    
    !*******************************
    type, public :: point_upstream
        sequence
        real :: coor
        integer :: id
    end type

    type(point_upstream),allocatable,target,public :: vertex_star(:)
    
    
    !*******************************

    type, public :: type_segment
        sequence
        type(point_upstream) :: porigin,pend
        integer :: id
    end type
    
    
    !*******************************
    type, public :: element1d_upstream
        sequence
        type(point_upstream) :: point_origin,point_end
        type(point_upstream) :: point_inter(0:11)
        type(type_segment) :: segment(10)
        integer :: nsub

        real :: x(1:2)
        real :: xg_star(1:5)
        real :: aa(1:5,1:5)

        real :: umodal(6)

    end type

    type(element1d_upstream),allocatable,target,public :: element_star(:)

    end module module_upstream_data