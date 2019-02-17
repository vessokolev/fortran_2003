program test_compute_rot_m
!
! Testing the subroutine 'compute_rot_m' supporting the rotation
! of a vector about one of the following axes:
!
! x -> [1,0,0]
! y -> [0,1,0]
! z -> [0,0,1]
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019021600
! License: GPLv2
!
use iso_c_binding,only:C_INT,C_FLOAT
use mod_rotate,only:compute_rot_m
!
implicit none
!
real(C_FLOAT) :: angle
real(C_FLOAT) :: m(3,3)

angle=0.172
!
call compute_rot_m(1,angle,m)
!
print *
print *,'Angle:',angle
print *
print *,'The elements of the matrix supporting the rotation ',&
        'about x-axis:'
print *
print *,'Rx(1,1)=',m(1,1),'Rx(1,2)=',m(2,1),'Rx(1,3)=',m(3,1)
print *,'Rx(2,1)=',m(1,2),'Rx(2,2)=',m(2,2),'Rx(2,3)=',m(3,2)
print *,'Rx(3,1)=',m(1,3),'Rx(3,2)=',m(2,3),'Rx(3,3)=',m(3,3)
!
print *
!
call compute_rot_m(2,angle,m)
!
print *
print *,'Angle:',angle
print *
print *,'The elements of the matrix supporting the rotation ',&
        'about y-axis:'
print *
print *,'Ry(1,1)=',m(1,1),'Ry(1,2)=',m(2,1),'Ry(1,3)=',m(3,1)
print *,'Ry(2,1)=',m(1,2),'Ry(2,2)=',m(2,2),'Ry(2,3)=',m(3,2)
print *,'Ry(3,1)=',m(1,3),'Ry(3,2)=',m(2,3),'Ry(3,3)=',m(3,3)
!
print *
!
call compute_rot_m(3,angle,m)
!
print *
print *,'Angle:',angle
print *
print *,'The elements of the matrix supporting the rotation ',&
        'about z-axis:'
print *
print *,'Rz(1,1)=',m(1,1),'Rz(1,2)=',m(2,1),'Rz(1,3)=',m(3,1)
print *,'Rz(2,1)=',m(1,2),'Rz(2,2)=',m(2,2),'Rz(2,3)=',m(3,2)
print *,'Rz(3,1)=',m(1,3),'Rz(3,2)=',m(2,3),'Rz(3,3)=',m(3,3)
!
print *

end program test_compute_rot_m
