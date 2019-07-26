subroutine s_bubble_sort_i(array_s,array_o,array_i)
!
! Bubble sort method:
!
! https://en.wikipedia.org/wiki/Bubble_sort
!
! This is the integer version of the method.
!
! The code bellow sorts the elements of the integer array
! 'array_o', by employing the bubble sort method. In parallel,
! it manager another array, named 'array_i', whose elements are
! the original indexes (positions) of the elements.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019072600
! License: GPLv2
!
use iso_c_binding,only:C_INT
!
implicit none
!
! Interface variables:
!
integer(C_INT),intent(in)    :: array_s
integer(C_INT),intent(inout) :: array_o(array_s)
integer(C_INT),intent(out)   :: array_i(array_s)
!
! Local variables:
!
integer(C_INT)               :: i
integer(C_INT)               :: shift_i
integer(C_INT)               :: buffer_i
integer(C_INT)               :: buffer_r
logical                      :: flag

! Save the orinal indexes of the elements as array elements
! of 'array_i':
!
do i=1,array_s
   !
   array_i(i)=i
   !
end do
!
flag=.true.
!
do while(flag)
   !
   flag=.false.
   !
   do i=1,array_s-1
      !
      shift_i=i+1
      !
      if (array_o(i)>array_o(shift_i)) then
         !
         ! Shift the value:
         !
         buffer_r=array_o(i)
         !
         array_o(i)=array_o(shift_i)
         !
         array_o(shift_i)=buffer_r
         !
         ! Shift the index:
         !
         buffer_i=array_i(i)
         !
         array_i(i)=array_i(shift_i)
         !
         array_i(shift_i)=buffer_i
         !
         if (.not.flag) then
            !
            flag=.true.
            !
         end if
         !
      end if
      !
   end do
   !
end do

end subroutine s_bubble_sort_i

