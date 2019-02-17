program test_calc_improp_dih_angle
!
! This program illustrates the way to compute an improper dihedral
! angle by using the subroutine 'calc_prop_dih_angle' from the module
! 'mod_rotate.f90'.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019021600
! License: GPLv2
!
use iso_c_binding,only:C_INT,C_FLOAT
use mod_rotate,only:calc_improp_dih_angle,rad2deg
!
implicit none
!
integer(C_INT) :: i
integer(C_INT) :: j
integer(C_INT) :: k
integer(C_INT) :: l
integer(C_INT) :: m
real(C_FLOAT)  :: angle
real(C_FLOAT)  :: coord(3,4)
logical        :: flag

! Define the set of input coordinates:
!
coord(:,1)=[ 0.151, -9.828,-3.229] !5
coord(:,2)=[-0.381, -6.636,-1.200] !7
coord(:,3)=[-3.823,-10.896, 0.354] !8
coord(:,4)=[ 1.090, -3.122,-1.054] !9
!
! Set the sequence of coordinates defining the improper dihedral angle:
!
i=1
j=2
k=4
l=3
!
! Preview of the construction of the improper dihedral angle:
!
!    l
!    |
!    i
!  /   \
!  j   k
!
! Do not keep the computed angle in [0,2*pi):
!
flag=.false.
!
! Call the procedure for computing the improper dihedral angle:
!
call calc_improp_dih_angle(coord,i,j,k,l,angle,flag)
!
print *
print *,'Input coordinates :'
print *
!
do m=1,4
   !
   write(*,fmt='(I2,1x,3F8.3)') m,coord(:,m)
   !
end do
!
print *
print *,'The improper dihedral is defined as the angle between '//&
        'the vecotrs normal to the planes:'
print *
write(*,fmt='(4I2)') i,j,k
print *
print *,'and'
print *
write(*,fmt='(4I2)') j,k,l
print *
print *,'The computed improper dihedral angle is:',angle*rad2deg,'degrees'
print *

end program test_calc_improp_dih_angle
