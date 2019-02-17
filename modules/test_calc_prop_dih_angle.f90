program test_calc_prop_dih_angle
!
! This program illustrates the way to compute a proper dihedral angle by
! using the subroutine 'calc_prop_dih_angle' from the module
! 'mod_rotate.f90'.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019021600
! License: GPLv2
!
use iso_c_binding,only:C_INT,C_FLOAT
use mod_rotate,only:calc_prop_dih_angle
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
real(C_FLOAT),parameter :: rad2deg=180/4.0/atan(1.0)
logical        :: flag

! Define the set of input coordinates:
!
coord(:,1)=[-1.760,-13.138,-3.084]
coord(:,2)=[-7.187,-12.315,-1.662]
coord(:,3)=[ 0.151, -9.828,-3.229]
coord(:,4)=[ 0.374,-10.474,-9.181]
!
! Set the sequence of coordinates defining the proper dihedral angle:
!
i=2
j=1
k=3
l=4
!
! Do not keep the computed angle in [0,2*pi):
!
flag=.false.
!
! Call the procedure for computing the proper dihedral angle:
!
call calc_prop_dih_angle(coord,i,j,k,l,angle,flag)
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
print *,'The proper dihedral is defined as the following '//&
        'sequence of coordinates:'
print *
write(*,fmt='(4I2)') i,j,k,l
print *
print *,'The computed proper dihedral angle is:',angle*rad2deg,'degrees'
print *

end program test_calc_prop_dih_angle
