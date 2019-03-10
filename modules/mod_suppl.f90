module mod_suppl
!
! Module that provides some useful supplementary subroutines.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019030900
! License: GPLv2
!
use iso_c_binding,only:C_INT,C_FLOAT
!
implicit none


contains


subroutine box_muller(rnd,sigma,mu,flag)
!
! Box-Muller method for generating normally distributed values with
! sigma = 1.0 and mu = 0.0. If you need to generate normally distributed
! data with sigma different than 1.0, and mu not equal to 0.0, use the
! transformation:
!
! z = z * sigma + mu
! 
! That requires to set 'flag' .true.
!
! The version of the subroutine presented bellow is CPU-effective!
!
! Interface variables:
!
real(C_FLOAT),intent(out) :: rnd
real(C_FLOAT),intent(in)  :: sigma
real(C_FLOAT),intent(in)  :: mu
logical,intent(in)        :: flag
!
! Local variables:
!
integer(C_INT),save       :: counter=0
real(C_FLOAT),parameter   :: pi=4*atan(1.0)
real(C_FLOAT),save        :: x
real(C_FLOAT),save        :: y
logical,save              :: stat=.false.

! Make the subroutine CPU effective, by calling 'random_number' per
! call. Call 'random_number' only if 'counter' becomes zero.
!
if (counter==0) then
   !
   call random_number(x)
   !
   call random_number(y)
   !
   counter=2
   !
end if
!
if (stat) then
   !
   stat=.false.
   !
   rnd=sqrt(-2*log(x))*cos(2*pi*y)
   !
else
   !
   stat=.true.
   !
   rnd=sqrt(-2*log(x))*sin(2*pi*y)
   !
end if
!
! Check if there both 'sigma' and 'mu' are actual arguments.
!
if (flag) then
   !
   rnd=rnd*sigma+mu
   !
end if
!
counter=counter-1

end subroutine box_muller


end module mod_suppl
