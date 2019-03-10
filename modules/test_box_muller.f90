program test_box_muller
!
! Small program for testing the subroutine 'box_muller' provided by the
! module 'mod_suppl'. The program calls the subroutine N times, stores
! the results in a histogram, and saves the histogram as a file. That
! file contains two columns - the first one is the center of the bin
! and the second column is the normalized normal distribution.
!
! More on the normal distribution on Wikipedia:
!
! https://en.wikipedia.org/wiki/Normal_distribution
!
! NOTES: The name of the output file is 'histogram_n.txt'!
!
!        The program generates normal distribution with mean 0.0 and
!        variance 2.0.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019030900
! License: GPLv2
!
use iso_c_binding,only:C_INT,C_LONG,C_FLOAT
use mod_suppl,only:box_muller
!
implicit none
!
integer(C_INT)              :: i
integer(C_INT)              :: j
integer(C_INT)              :: i_ll
integer(C_INT)              :: i_ul
integer(C_INT)              :: num_bins
integer(C_LONG)             :: num_trials
integer(C_LONG),allocatable :: acc(:)
real(C_FLOAT)               :: hist_ll
real(C_FLOAT)               :: hist_ul
real(C_FLOAT)               :: bin_w
real(C_FLOAT)               :: sigma
real(C_FLOAT)               :: rnd

! Specify the lower and upper limits of the interval:
!
hist_ll=-10.0
hist_ul=-hist_ll
!
! the number of bins, number of trials, and the sigma parameter of the
! distribution:
!
num_bins=1000
num_trials=10000000
sigma=2.0
!
! Estimate the bin width:
!
bin_w=(hist_ul-hist_ll)/(num_bins-1)
!
!
!hist_ll=hist_ll-bin_w_h/2
!hist_ul=hist_ul-bin_w_h/2
!
! Estimate the lower and upper index of the elements of the array, that
! will store the bins:
!
i_ll=floor(hist_ll/bin_w)
i_ul=floor(hist_ul/bin_w)
!
! Allocate the array for storing the bins (note the explicit use of the
! lower and upper index - since the index of the first element of the
! array is not 1:
!
allocate(acc(i_ll:i_ul))
!
do i=1,num_trials
   !
   ! Guess a normally distributed random number:
   !
   call box_muller(rnd,sigma,0.0,.false.)
   !
   ! Compute the index of the bin:
   !
   j=floor(rnd/bin_w)
   !
   ! Add 1 to the corresponding bin:
   !
   acc(j)=acc(j)+1
   !
end do
!
! Save the histogram as file:
!
open(666,file='histogram_n.txt',status='unknown')
!
do i=i_ll,i_ul
   !
   write(666,fmt='(E14.7,1x,E14.7)') (i+0.5)*bin_w,real(acc(i))/num_trials
   !
end do
!
close(666)
!
deallocate(acc)

end program test_box_muller
