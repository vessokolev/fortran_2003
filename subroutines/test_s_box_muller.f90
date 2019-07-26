program test_s_box_muller

! Simple program that tests the subroutine 's_box_muller'.
!
! The program defines 5000 bins and populate them by generating
! 10000000 tests and distribute them among the bins.
!
! HOWTO COMPILE AND RUN THE CODE:
!
! 1. GNU Fortran:
!
!   gfortran -c s_box_muller.f90
!   gfortran -o test_s_box_muller test_s_box_muller.f90 s_box_muller.o
!   ./test_s_box_muller > hist.txt
!
! 2. Intel Fortran:
!
!   ifort -c s_box_muller.f90
!   ifort -o test_s_box_muller test_s_box_muller.f90 s_box_muller.o
!   ./test_s_box_muller > hist.txt
!
! 3. PGI Fortran:
!
!   pgfortran -c s_box_muller.f90
!   pgfortran -o test_s_box_muller test_s_box_muller.f90 s_box_muller.o
!   ./test_s_box_muller > hist.txt
!
! HOWTO VISUALIZE THE RESULTS
!
! Use gnuplot and plot the data from 'hist.txt'


use iso_c_binding,only:C_INT,C_FLOAT
!
implicit none
!
integer(C_INT),parameter :: num_bins=5000
integer(C_INT),parameter :: num_trials=10000000
integer(C_INT)           :: i
integer(C_INT)           :: indx
integer(C_INT)           :: bin(num_bins)
real(C_FLOAT),parameter  :: hist_ll=-6.0
real(C_FLOAT),parameter  :: hist_ul= 6.0
real(C_FLOAT)            :: bin_w
real(C_FLOAT)            :: rnd

! Initialize the bin values (set 0):
!
bin(:)=0
!
! Compute the bin width:
!
bin_w=(hist_ul-hist_ll)/num_bins
!
! Generate the tests:
!
do i=1,num_trials
   !
   call s_box_muller(rnd,1.0,0.0,.false.)
   !
   indx=(rnd-hist_ll)/bin_w
   !
   bin(indx)=bin(indx)+1
   !
end do
!
! Print out the collected PMF (Probability Mass Function) of the
! Normal Distribution (direct the output to a file):
!
do i=1,num_bins
   !
   print *,hist_ll+bin_w*i,1.0*bin(i)/num_trials
   !
end do

end program test_s_box_muller
