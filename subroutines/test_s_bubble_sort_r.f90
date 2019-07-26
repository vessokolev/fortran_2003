program test_s_bubble_sort_r

! Small program for testing the subroutine 's_bubble_sort_r'.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019072600
! License: GPLv2
!
! HOWTO COMPILE AND RUN THE CODE:
!
! 1. GNU Fortran:
!
!   gfortran -c s_bubble_sort_r.f90
!   gfortran -o test_s_bubble_sort_r test_s_bubble_sort_r.f90 s_bubble_sort_r.o
!   ./test_s_bubble_sort_r
!
! 2. Intel Fortran:
!
!   ifort -c s_bubble_sort_r.f90
!   ifort -o test_s_bubble_sort_r test_s_bubble_sort_r.f90 s_bubble_sort_r.o
!   ./test_s_bubble_sort_r
!
! 3. PGI Fortran:
!
!   pgfortran -c s_bubble_sort_r.f90
!   pgfortran -o test_s_bubble_sort_r test_s_bubble_sort_r.f90 s_bubble_sort_r.o
!   ./test_s_bubble_sort_r


use iso_c_binding,only:C_INT,C_FLOAT
!
implicit none
!
integer(C_INT),parameter :: array_s=10
integer(C_INT)           :: array_i(array_s)
real(C_FLOAT)            :: array_o(array_s)

! The original array elements:
!
array_o=[0.94340948, 0.27256313, 0.3733053 , 0.03770947, 0.72549825, &
         0.79981257, 0.33502679, 0.51859482, 0.59218711, 0.2659781 ]
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
call s_bubble_sort_r(array_s,array_o,array_i)
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

end program test_s_bubble_sort_r

