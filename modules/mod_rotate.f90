module mod_rotate
!
! This module supports the rotation of the sidechain branches
! around selected axis and computes the proper dihedral angles.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019020600
! License: GPLv2
!

use iso_c_binding,only:C_INT,C_FLOAT

implicit none

public :: select_branches_to_rotate
public :: rotate_branch_around_axis
public :: norm02
public :: compute_pdih
public :: rotate_vect_around_axis
public :: get_rodrigues_matrix
public :: vector_cross_product


contains


subroutine select_branches_to_rotate(coord,branches,branches_i,center,&
                                     cutoff,is_rotatable,num_rotatable)
!
! Selects the branches falling within some "cutoff" from "center"
! (both "cutoff" and "center" are defined by the user). Only the
! branches having at least one of their atoms within the "cutoff"
! distance are marked as rotatable (is_rotatable(x)=.true. for them).
!
! Interface variables:
!
real(C_FLOAT),intent(in)   :: coord(:,:)
integer(C_INT),intent(in)  :: branches(:)
integer(C_INT),intent(in)  :: branches_i(:,:)
real(C_FLOAT),intent(in)   :: center(3)
real(C_FLOAT),intent(in)   :: cutoff
logical,intent(out)        :: is_rotatable(:)
integer(C_INT),intent(out) :: num_rotatable
!
! Local variables:
!
integer(C_INT)             :: i
integer(C_INT)             :: j
real(C_FLOAT)              :: cutoff_

cutoff_=cutoff**2
!
is_rotatable=.false.
!
num_rotatable=0
!
! If one atom of a given branch is fount to be rotatable, mark the
! whole branch as rotatable too.
!
do i=1,ubound(branches_i,2)
   !
   do j=branches_i(1,i),branches_i(2,i)
      !
      ! Use the square of the cutoff. That saves CPU time for not
      ! calling the 'sqrt' intrinsic function (that function is
      ! not fast enough).
      !
      if (sum((coord(:,branches(j))-center)**2)<cutoff_) then
         !
         is_rotatable(i)=.true.
         !
         num_rotatable=num_rotatable+1
         !
         exit
         !
      end if
      !
   end do
   !
end do

end subroutine select_branches_to_rotate


subroutine rotate_branch_around_axis(coord,axes,branches,branches_i,sel,&
                                     angle,bak_c)
!
! Rotates a branch of points around an axis.
!
! Both branch members and two atoms defining the axis, have their
! coordinates listed in the array "coord" (a 3xN array of floats).
! The array "axes" contains all pairs of indexes in "coord",
! that can define an axis. The variable "sel" point to the selected
! axis. Every branch contains the list of indexes of elements of
! "coord" that should be rotated around the selected axis. The array
! "branches_i" (containing 2 integer numbers), specifies the range
! of "branches" elements that point to the indexes of the elements of
! "coord", participating in the rotation process.
!
! Interface variables:
!
real(C_FLOAT),intent(inout) :: coord(:,:)
integer(C_INT),intent(in)   :: axes(:,:)
integer(C_INT),intent(in)   :: branches(:)
integer(C_INT),intent(in)   :: branches_i(:,:)
integer(C_INT),intent(in)   :: sel
real(C_FLOAT),intent(in)    :: angle
real(C_FLOAT),intent(out)   :: bak_c(:,:)
!
! Local variables:
!
integer(C_INT)              :: i
real(C_FLOAT)               :: rodrigues_matrix(3,3)

call get_rodrigues_matrix(coord(:,axes(1,sel))-coord(:,axes(2,sel)),&
                          angle,rodrigues_matrix)
!
do i=branches_i(1,sel),branches_i(2,sel)
   !
   bak_c(:,i-branches_i(1,sel)+1)=coord(:,branches(i))
   !
   coord(:,branches(i))=coord(:,branches(i))-coord(:,axes(1,sel))
   !
   call rotate_vect_around_axis(coord(:,branches(i)),rodrigues_matrix)
   !
   coord(:,branches(i))=coord(:,branches(i))+coord(:,axes(1,sel))
   !
end do

end subroutine rotate_branch_around_axis


subroutine restore_branch(coord,branches,branches_i,sel,bak_c)
!
! The subroutine restores the coordinates of the atoms, previously
! rotated around the axis defined by "axes(:,sel)". Note that the array
! "bak_c" contains the coordinates of the atoms before the rotation
! to take place.
!
! Interface variables:
!
real(C_FLOAT),intent(inout) :: coord(:,:)
integer(C_INT),intent(in)   :: branches(:)
integer(C_INT),intent(in)   :: branches_i(:,:)
integer(C_INT),intent(in)   :: sel
real(C_FLOAT),intent(in)    :: bak_c(:,:)
!
! Local variables:
!
integer(C_INT)              :: i

forall (i=branches_i(1,sel):branches_i(2,sel)) coord(:,branches(i))=&
    bak_c(:,i-branches_i(1,sel)+1)

end subroutine restore_branch


subroutine rotate_vect_around_axis(vect,rodrigues_matrix)
!
! Rotates 3D vector around 3D axis.
!
! NOTES: "vect" is the vector to rotate. The elements of the matrix
!        named "rodrigues_matrix" can be computed by calling the
!        subroutine "get_rodrigues_matrix".
!
! Interface variables:
!
real(C_FLOAT),intent(inout) :: vect(3)
real(C_FLOAT),intent(in)    :: rodrigues_matrix(3,3)
!
! Local variables:
!
real(C_FLOAT)               :: rotated(3)
integer(C_INT)              :: i

forall (i=1:3) rotated(i)=sum(rodrigues_matrix(:,i)*vect(:))
!
vect=rotated(:)

end subroutine rotate_vect_around_axis


subroutine get_rodrigues_matrix(vect,angle,matrix)

! Calculates the rotation matrix elements that allow the
! a pplication of the Rodrigues' formula. For more details see:
! https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

! NOTES: "vect" defines the 3D equation of the axis line.
!        It MUST begin at (0,0,0). Supply the angle "angle"
!        in radians!

! IMPORTANT!!!

! The matrix indexing foolows the Fortran way, not the
! one of C. Do take this into account!
!
! Interface variables:
!
real(C_FLOAT),intent(in)  :: vect(3)
real(C_FLOAT),intent(in)  :: angle
real(C_FLOAT),intent(out) :: matrix(3,3)
!
! Local variables:
!
integer(C_INT)            :: i
integer(C_INT),parameter  :: factor_1(2,3)=reshape([1,2,1,3,2,3],[2,3])
integer(C_INT),parameter  :: factor_2(2,3)=reshape([2,1,3,2,1,3],[2,3])
integer(C_INT),parameter  :: factor_3(2,3)=reshape([1,3,3,1,2,2],[2,3])
integer(C_INT),parameter  :: factor_4(2,3)=reshape([3,1,1,2,2,3],[2,3])
integer(C_INT),parameter  :: factor_5(2,3)=reshape([2,2,1,3,3,1],[2,3])
real(C_FLOAT)             :: norm(3)
real(C_FLOAT)             :: cos_
real(C_FLOAT)             :: sin_
real(C_FLOAT)             :: omcos_
real(C_FLOAT)             :: n(3)

norm=vect/norm02(vect)
!
cos_=cos(angle)
!
sin_=sin(angle)
!
omcos_=1-cos_
!
forall (i=1:3) n(i)=norm(factor_1(1,i))*norm(factor_1(2,i))
!
do i=1,3
   !
   matrix(i,i)=cos_+omcos_*norm(i)**2
   !
   matrix(factor_2(1,i),factor_2(2,i))=n(factor_3(1,i))*omcos_-&
                                       norm(factor_3(2,i))*sin_
   !
   matrix(factor_4(1,i),factor_4(2,i))=n(factor_5(1,i))*omcos_+&
                                       norm(factor_5(2,i))*sin_
   !
end do

end subroutine get_rodrigues_matrix


function compute_pdih(coord,i,j,k,l) result(angle)
!
! Computes the proper dihedral angle by using 'atan2'.
! The result is in radians and single precision. It
! is located within [-pi,pi].
!
! Interface variables:
!
real(C_FLOAT),intent(in)  :: coord(:,:)
integer(C_INT),intent(in) :: i
integer(C_INT),intent(in) :: j
integer(C_INT),intent(in) :: k
integer(C_INT),intent(in) :: l
real(C_FLOAT)             :: angle
!
! Local variables:
!
real(C_FLOAT)             :: b1(ubound(coord,1))
real(C_FLOAT)             :: b2(ubound(coord,1))
real(C_FLOAT)             :: b3(ubound(coord,1))
real(C_FLOAT)             :: n1(ubound(coord,1))
real(C_FLOAT)             :: n2(ubound(coord,1))
real(C_FLOAT),parameter   :: pix2=8*atan(1.0)

b1=coord(:,i)-coord(:,j)
!
b2=coord(:,j)-coord(:,k)
!
b3=coord(:,k)-coord(:,l)
!
! Compute n1 =  b1 x b2
!
n1=vector_cross_product(b1,b2)
!
! Normalize n1:
!
n1=n1/norm02(n1)
!
! Note that after n1 is computed there is no need to keep b1
! anymore. That means b1 might be used as a dummy variable to help
! reducing the number of declared variables.
!
! Compute n2 =  b2 x b3
!
n2=vector_cross_product(b2,b3)
!
! Normalize n2:
!
n2=n2/norm02(n2)
!
! Normalize b2:
!
b2=b2/norm02(b2)
!
! Bellow this line b1 is used as a dummy vector variable. So it
! can be employed in the product:
!
! m1 = n1 x b2
!
b1=vector_cross_product(n1,b2)
!
! Compute the dihedral angle using atan2:
!
angle=atan2(dot_product(b1,n2),dot_product(n1,n2))

end function compute_pdih


function norm02(vector) result(res)
!
! This function computes the norm of a vector in 3D. The result is
! single precisions floating point number.
!
! NOTICE: Some compilers (GCC gfortran, Intel Fortran) support the
!         intrinsic function 'norm2'. Unfortunally, that function
!         is not implemented in other compilers, like those by PGI.
!
! Interface variables:
!
real(C_FLOAT),intent(in) :: vector(:)
real(C_FLOAT)            :: res

! Total vectorization (implicit indexing):
!
res=sqrt(sum(vector**2))

end function norm02


function vector_cross_product(v1,v2) result(prod)
!
! Computes the cross product of two 3D vectors.
! The index selection is based on the 3D Levi-Civita
! symbol.
!
! Interface varirables:
!
real(C_FLOAT),intent(in) :: v1(:)
real(C_FLOAT),intent(in) :: v2(:)
real(C_FLOAT)            :: prod(ubound(v2,1))
!
! Local variables:
!
integer(C_INT),parameter :: factor1(3)=(/2,3,1/)
integer(C_INT),parameter :: factor2(3)=(/3,1,2/)

! Total vectorization (implicit indexing):
!
prod=v1(factor1)*v2(factor2)-v1(factor2)*v2(factor1)
!
end function vector_cross_product


subroutine angle_between_vectors(v1,v2,angle,vprod)
!
! Computes the angle between two vectors in 3D.
! It returns the angle in radians (as 'angle') and
! the vector product (as 'vprod')
!
! IMPORTANT: The vectors MUST begin at (0,0,0)!
!
! Interface variables:
!
real(C_FLOAT),intent(in)  :: v1(3)
real(C_FLOAT),intent(in)  :: v2(3)
real(C_FLOAT),intent(out) :: angle
real(C_FLOAT),intent(out) :: vprod(3)

vprod=vector_cross_product(v1,v2)
!
angle=atan2(sqrt(dot_product(vprod,vprod)),dot_product(v1,v2))

end subroutine angle_between_vectors


end module mod_rotate

