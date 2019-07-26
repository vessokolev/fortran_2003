program test_s_bubble_sort_i

! Small program for testing the subroutine 's_bubble_sort_i'.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019072600
! License: GPLv2
!
! HOWTO COMPILE AND RUN THE CODE:
!
! 1. GNU Fortran:
!
!   gfortran -c s_bubble_sort_i.f90
!   gfortran -o test_s_bubble_sort_i test_s_bubble_sort_i.f90 s_bubble_sort_i.o
!   ./test_s_bubble_sort_i
!
! 2. Intel Fortran:
!
!   ifort -c s_bubble_sort_i.f90
!   ifort -o test_s_bubble_sort_i test_s_bubble_sort_i.f90 s_bubble_sort_i.o
!   ./test_s_bubble_sort_i
!
! 3. PGI Fortran:
!
!   pgfortran -c s_bubble_sort_i.f90
!   pgfortran -o test_s_bubble_sort_i test_s_bubble_sort_i.f90 s_bubble_sort_i.o
!   ./test_s_bubble_sort_i

use iso_c_binding,only:C_INT
!
implicit none
!
integer(C_INT),parameter :: array_s=10
integer(C_INT)           :: array_i(array_s)
integer(C_INT)           :: array_o(array_s)

! The original array elements:
!
array_o=[219, 27, 908 , 743, 801, 12, 102, 449, 528, 204 ]
!
print *
!
print *,'The array before sorting:'
!
print *,array_o
!
print *
!
! Call the subroutine:
!
call s_bubble_sort_i(array_s,array_o,array_i)
!
print *,'The array after sorting:'
!
print *,array_o
!
print *
!
print *,'The permutation of the indexes due to the sorting:'
!
print *,array_i
!
print *

end program test_s_bubble_sort_i

