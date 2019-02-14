module mod_rotate
!
! Module supporting the rotation of vectors in 3D.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019021400
! License: GPLv2
!
use iso_c_binding,only:C_INT,C_FLOAT

implicit none


contains


subroutine rotate_branch_around_axis(axis,branch,angle,coord)
!
! Rotates a set of points (branch) about an arbitrary axis in 3D.
!
! axis - the indexes of the points defining the axis (the
!        axis is a 3D line) - two point indexes are needed,
!        and therefore the size of 1D array 'axis' is 2
!
! branch - 1D array of integers, which are the indexes of the
!          points to rotate
!
! angle - the angle of rotation (IN RADIANS!!!)
!
! coord - 3 x num_points - shaped array of floating point real
!         numbers, containing the coordinates of each point of
!         the set in 3D space of coordinates
!
! NOTES: branch(i) points to coord(:,branch(i))
!        axis(1) points to coord(:,axis(1))
!        axis(2) points to coord(:,axis(2))
!
! Interface variables:
!
integer(C_INT),intent(in)   :: axis(:)
integer(C_INT),intent(in)   :: branch(:)
real(C_FLOAT),intent(in)    :: angle
real(C_FLOAT),intent(inout) :: coord(:,:)
!
! Local variables:
!
integer(C_INT)              :: i
real(C_FLOAT)               :: matrix(3,3)

call get_rodrigues_matrix(coord(:,axis(2))-&
                          coord(:,axis(1)),&
                          angle,matrix)
!
do i=1,ubound(branch,1)
   !
   coord(:,branch(i))=coord(:,branch(i))-coord(:,axis(1))
   !
   call rotate_around_axis(coord(:,branch(i)),matrix)
   !
   coord(:,branch(i))=coord(:,branch(i))+coord(:,axis(1))
   !
end do

end subroutine rotate_branch_around_axis


subroutine rotate_around_axis(vect,rodrigues_matrix)
!
! Rotates 3D vector around 3D axis.
!
! NOTES: "vect" is the vector to rotate. The elements of the matrix
!        named "rodrigues_matrix" can be computed by calling the
!        subroutine "get_rodrigues_matrix".
!
!        IT IS IMPORTANT to note that the vector 'vect' MUST
!        begin at [0,0,0]. If not - do shift it.
!
!        Details regarding the construction and implementation of the
!        Rodrigues' formula, can be found at:
!
!        https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
!
!        The matrix indexing foolows the Fortran way, not the
!        one of C. Do take this into account when reading the code. The
!        matrix 'rodrigues_matrix' used bellow is the transpose of the
!        one given in the books.
!
! Interface variables:
!
real(C_FLOAT),intent(inout) :: vect(3)
real(C_FLOAT),intent(in)    :: rodrigues_matrix(3,3)
!
! Local variables:
!
integer(C_INT)              :: i

! Using 'forall' implicitly declares intermediate storage slot for
! storing there the initial components of 'vect'. Hence, instead of
! explicitly declaring an internal (to the subroutine) storage variable,
! which in turn requires a 'do'-loop, one could simply employ 'forall'
! and make the code shorter.
!
forall (i=1:3) vect(i)=sum(rodrigues_matrix(:,i)*vect(:))

end subroutine rotate_around_axis


subroutine get_rodrigues_matrix(vect,angle,matrix)
!
! Calculates the rotation matrix elements that allow the
! application of the Rodrigues' formula. For more details see:
!
! https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
!
! NOTES: "vect" defines the 3D equation of the axis line.
!        It MUST begin at (0,0,0). Supply the angle "angle"
!        in radians!
!
!        The matrix indexing foolows the Fortran way, not the
!        one of C. Do take this into account when reading the code. The
!        matrix 'rodrigues_matrix' used bellow is the transpose of the
!        one given in the books.
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
sin_=sin(angle)
!
omcos_=1-cos_
!
forall (i=1:3) n(i)=norm(factor_1(1,i))*norm(factor_1(2,i))
!
forall (i=1:3) matrix(i,i)=cos_+omcos_*norm(i)**2
!
forall (i=1:3) matrix(factor_2(1,i),factor_2(2,i))=&
                  n(factor_3(1,i))*omcos_-&
                  norm(factor_3(2,i))*sin_
!
forall (i=1:3) matrix(factor_4(1,i),factor_4(2,i))=&
                  n(factor_5(1,i))*omcos_+&
                  norm(factor_5(2,i))*sin_

end subroutine get_rodrigues_matrix


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

call vector_cross_product(v1,v2,vprod)
!
angle=atan2(sqrt(dot_product(vprod,vprod)),dot_product(v1,v2))

end subroutine angle_between_vectors


subroutine vector_cross_product(v1,v2,prod)
!
! Computes the cross product of two vectors in 3D.
!
! Interface varirables:
!
real(C_FLOAT),intent(in)  :: v1(3)
real(C_FLOAT),intent(in)  :: v2(3)
real(C_FLOAT),intent(out) :: prod(3)
!
! Local variables:
!
integer(C_INT)            :: i
integer(C_INT),parameter  :: factor1(3)=(/2,3,1/)
integer(C_INT),parameter  :: factor2(3)=(/3,1,2/)

prod(:)=v1(factor1(:))*v2(factor2(:))-v1(factor2(:))*v2(factor1(:))

end subroutine vector_cross_product


function norm02(argument) result(res)
!
! This function computes the norm of the N-dimensional
! vector of floating point numbers.
!
! Interface variables:
!
real(C_FLOAT),intent(in)  :: argument(:)
real(C_FLOAT)             :: res

res=sqrt(sum(argument**2))

end function norm02


subroutine compute_rot_m(axis_num,angle,matrix)
!
! Computes the components of the matrix supporting the rotation
! of a vector about the axis x, y, or z:
!
! axis_num=1 -> x (the unit vector is [1,0,0])
! axis_num=2 -> y (the unit vector is [0,1,0])
! axis_num=3 -> z (the unit vector is [0,0,1])
!
! by angle 'angle' (given in radians).
!
! NOTE: The indexes of the matrix follow the Fortran way of
!       accessing the memory. That means the matrix defined
!       bellow is the transpose of the rotation matrix given
!       in the books.
!
! Interface variables:
!
integer(C_INT),intent(in) :: axis_num
real(C_FLOAT),intent(in)  :: angle
real(C_FLOAT),intent(out) :: matrix(3,3)
!
! Local variables:
!
integer(C_INT)            :: i
integer(C_INT)            :: j
real(C_FLOAT)             :: arr(5)
integer(C_INT),parameter  :: factor(3,3,3)=[&
                                  reshape([2,1,1,1,3,5,1,4,3],[3,3]),&
                                  reshape([3,1,4,1,2,1,5,1,3],[3,3]),&
                                  reshape([3,5,1,4,3,1,1,1,2],[3,3])]

arr(1)=0.0
arr(2)=1.0
!
! Compute the sin and cos in advance and save the results as
! local variables (reduces the CPU cost).
!
arr(3)=cos(angle)
arr(4)=sin(angle)
arr(5)=-arr(4)
!
! Note, that the second index of 'matrix' specifies the row number
! of the rotation matrix (see the NOTE above):
!
do j=1,3
   !
   do i=1,3
      !
      matrix(i,j)=arr(factor(i,j,axis_num))
      !
   end do
   !
end do

end subroutine compute_rot_m


end module mod_rotate
