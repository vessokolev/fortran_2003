module mod_rotate
!
! Module supporting the rotation of vectors in 3D.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019021404
! License: GPLv2
!
use iso_c_binding,only:C_INT,C_FLOAT

implicit none


contains


subroutine rotate_branch_about_axis(axis,branch,angle,coord)
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
   call rotate_about_axis(coord(:,branch(i)),matrix)
   !
   coord(:,branch(i))=coord(:,branch(i))+coord(:,axis(1))
   !
end do

end subroutine rotate_branch_about_axis


subroutine rotate_about_axis(vect,rodrigues_matrix)
!
! Rotates 3D vector about an arbitraty 3D axis.
!
! NOTES: "vect" defines the 3D equation of the axis line. It MUST have
!        the initial point (origin) of one of the vectors defining the
!        axis line (shift the vector 'vect' with respect to that
!        vector). The elements of the matrix named "rodrigues_matrix"
!        can be computed by calling the subroutine
!        "get_rodrigues_matrix".
!
!        Details regarding the construction and implementation of the
!        Rodrigues' formula, can be found at:
!
!        https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
!
!        The matrix indexing follows the Fortran way, not the
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

! The use of 'forall' implicitly declares an intermediate storage slot
! for keeping the initial components of 'vect'. Hence, instead of
! explicitly declaring an internal (to the subroutine) storage variable,
! which is a must, when having a 'do'-loop constructor, one could simply
! employ 'forall' instead.
!
forall (i=1:3) vect(i)=sum(rodrigues_matrix(:,i)*vect(:))

end subroutine rotate_about_axis


subroutine get_rodrigues_matrix(vect,angle,matrix)
!
! Calculates the rotation matrix elements that allow the application of
! the Rodrigues' formula. For more details see:
!
! https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
!
! NOTES: "vect" defines the 3D equation of the axis line. It MUST have
!        the initial point (origin) of one of the vectors defining the
!        axis line (shift the vector 'vect' with respect to that
!        vector).
!
!        Supply the angle "angle" in radians!
!
!        The matrix indexing follows the Fortran way, not the
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
! Computes the angle between two vectors in 3D. It returns the angle
! in radians (as 'angle') and the vector product (as 'vprod').
!
! IMPORTANT: The vectors MUST share the same initial point (origin)!
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
integer(C_INT),parameter  :: factor1(3)=(/2,3,1/)
integer(C_INT),parameter  :: factor2(3)=(/3,1,2/)

prod(:)=v1(factor1(:))*v2(factor2(:))-v1(factor2(:))*v2(factor1(:))

end subroutine vector_cross_product


function norm02(argument) result(res)
!
! This function computes the norm of the N-dimensional vector of
! floating point numbers.
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
integer(C_INT)            :: selection(3,3)
integer(C_INT),parameter  :: limits(2,3)=reshape([1,3,4,6,7,9],[2,3])
integer(C_INT),parameter  :: ind(3,9)=reshape([2,1,1,1,3,5,1,4,3,&
                                               3,1,4,1,2,1,5,1,3,&
                                               3,5,1,4,3,1,1,1,2],[3,9])

selection(:,:)=ind(:,limits(1,axis_num):limits(2,axis_num))
!
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
      matrix(i,j)=arr(selection(i,j))
      !
   end do
   !
end do

end subroutine compute_rot_m


subroutine get_equation_of_plane(v1,v2,v3,equ)
!
! Computes the equation of a 3D plane in the form:
!
! f(x,y,z) = equ(1)*x + equ(2)*y + equ(3)*z + equ(4) = 0
!
! The parameters are the elements of the array "equ" and x, y, z are the
! components of a vector that belongs to the plane. Three vectors are
! required to compute the parameters 'equ(1:4).' The vector normal to
! the plane, n, can be obtained by nomalizing the vector 'equ(1:3)':
!
! n = equ(1:3)/norm02(equ(1:3))
!
! NOTE: The subroutine 'vector_cross_product' requires the input vectors
!       'v1' and 'v3' to share the same initial point. One way to do so
!       it to shift 'v1' and 'v3' by 'v2' (see the code bellow). 
!
! Interface variables:
!
real(C_FLOAT),intent(in)  :: v1(3)
real(C_FLOAT),intent(in)  :: v2(3)
real(C_FLOAT),intent(in)  :: v3(3)
real(C_FLOAT),intent(out) :: equ(4)

call vector_cross_product(v1-v2,v3-v2,equ(1:3))
!
equ(4)=v1(1)*equ(1)+v1(2)*equ(2)+v1(3)*equ(3)
!
equ(4)=-equ(4)

end subroutine get_equation_of_plane


end module mod_rotate
