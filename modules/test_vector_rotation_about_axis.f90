program test_vector_rotation_about_axis
!
! Simple program for testing the effectiveness of the subroutines:
!
! 'get_rodrigues_matrix' and 'rotate_about_axis'
!
! originally defined in 'mod_rotate.f90'. The program generates
! a series of PDB files containing the positions of the benzene atoms,
! rotated by means of the subroutines. To preview the result, invoke
! VMD as:
!
! vmd *.pdb
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019030900
! License: GPLv2
!
use iso_c_binding,only:C_INT,C_FLOAT
use mod_rotate,only:get_rodrigues_matrix,rotate_about_axis

integer(C_INT)              :: i
integer(C_INT)              :: j
integer(C_INT)              :: num_atoms
real(C_FLOAT)               :: angle
real(C_FLOAT)               :: axis(3)
real(C_FLOAT)               :: center(3)
real(C_FLOAT)               :: center_(3)
real(C_FLOAT)               :: matrix(3,3)
real(C_FLOAT),allocatable   :: coord(:,:)
real(C_FLOAT),parameter     :: pi=4.0*atan(1.0)
character(len=20)           :: file_name
!
! Benzene molecule in OPLS FF representation, in PDB format:
!
character(len=54),parameter :: line(15)=[&
'COMPND    UNNAMED                                     ',&
'AUTHOR    GENERATED BY OPEN BABEL 2.2.3               ',&
'HETATM    1  C1  BEN     1      -2.145   0.973  -0.003',&
'HETATM    2  H1  BEN     1      -3.103   0.460  -0.005',&
'HETATM    3  C2  BEN     1      -2.100   2.367   0.009',&
'HETATM    4  H2  BEN     1      -3.023   2.940   0.016',&
'HETATM    5  C3  BEN     1      -0.870   3.025   0.012',&
'HETATM    6  H3  BEN     1      -0.836   4.111   0.021',&
'HETATM    7  C4  BEN     1       0.314   2.290   0.003',&
'HETATM    8  H4  BEN     1       1.273   2.803   0.005',&
'HETATM    9  C5  BEN     1       0.270   0.895  -0.009',&
'HETATM   10  H5  BEN     1       1.193   0.322  -0.016',&
'HETATM   11  C6  BEN     1      -0.959   0.237  -0.012',&
'HETATM   12  H6  BEN     1      -0.995  -0.849  -0.021',&
'END                                                   ']

! Read the number of atom described in 'line':
!
call num_atoms_in_pdb(line,num_atoms)
!
! Allocate an array to store the atom coordinates in:
!
allocate(coord(3,num_atoms))
!
! Read the atom coordinates:
!
call read_atom_coord(line,num_atoms,coord)
!
! Compute the center of the atom coordinates (required
! to keep the molecule in the center of its own coordinate
! system):
!
center=0.0
!
do i=1,num_atoms
   !
   center=center+coord(:,i)
   !
end do
!
center=center/num_atoms
!
! Write down a PDB with the initial atom coordinates:
!
write(file_name,fmt='(I3.3,A4)') 0,'.pdb'
!
call write_pdb(file_name,line,coord)
!
! Start generating normal vectors of axes and perform a
! rotation about the axis (use randomly generated angle):
!
do i=1,100
   !
   write(file_name,fmt='(I3.3,A4)') i,'.pdb'
   !
   call random_number(axis)
   !
   call random_number(angle)
   !
   angle=2*(0.5-angle)*pi
   !
   call get_rodrigues_matrix(axis,angle,matrix)
   !
   do j=1,num_atoms
      !
      call rotate_about_axis(coord(:,j),matrix)
      !
   end do
   !
   ! Compute the cordinate center of the rotated atoms:
   !
   center_=0.0
   !
   do j=1,num_atoms
      !
      center_=center_+coord(:,j)
      !
   end do
   !
   center_=center_/num_atoms
   !
   ! Correct the rotated atoms so their center to be the one
   ! computed for the original coordinates:
   !
   forall(j=1:num_atoms) coord(:,j)=coord(:,j)-center_(:)+center(:)
   !
   ! Write down a PDB with the rotated atom coordinates:
   !
   call write_pdb(file_name,line,coord)
   !
end do


contains


subroutine num_atoms_in_pdb(line,num_atoms)
!
! Estimates the number of atoms records in the PDB description
! (stored in 'line').
!
! Interface variables:
!
character(len=*),intent(in) :: line(:)
integer(C_INT),intent(out)  :: num_atoms
!
! Local variables:
!
integer(C_INT)              :: i

! Collect the coordinates from the PDB HETATM lines.
!
! First, estimate the number of HETATM records:
!
num_atoms=0
!
do i=1,ubound(line,1)
   !
   if (line(i)(1:6)=='HETATM') then
      !
      num_atoms=num_atoms+1
      !
   end if
   !
end do

end subroutine num_atoms_in_pdb


subroutine read_atom_coord(line,num_atoms,coord)
!
! Collects the atoms coordinates from the PDB description
! (stored in 'line').
!
! Interface variables:
!
character(len=*),intent(in) :: line(:)
integer(C_INT),intent(out)  :: num_atoms
real(C_FLOAT),intent(out)   :: coord(:,:)
!
! Local variables:
!
integer(C_INT)              :: i

num_atoms=0
!
do i=1,ubound(line,1)
   !
   if (line(i)(1:6)=='HETATM') then
      !
      num_atoms=num_atoms+1
      !
      read(line(i)(31:54),fmt='(3F8.3)'),coord(:,num_atoms)
      !
   end if
   !
end do

end subroutine read_atom_coord


subroutine write_pdb(file_name,line,coord)
!
! Writes down a PDB files (the lines of the PDB file are
! stored as elements of 'line').
!
! Interface variables:
!
character(len=*),intent(in) :: file_name
character(len=*),intent(in) :: line(:)
real(C_FLOAT),intent(in)    :: coord(:,:)
!
! Local variables:
!
integer(C_INT)              :: i
integer(C_INT)              :: j
integer(C_INT)              :: stat
character(len=80)           :: line_
character(len=24)           :: coord_str

open(unit=666,file=file_name,iostat=stat,status='unknown')
!
if (stat==0) then
   !
   j=0
   !
   do i=1,ubound(line,1)
      !
      if (line(i)(1:6)=='HETATM') then
         !
         j=j+1
         !
         write(coord_str,fmt='(3F8.3)') coord(:,j)
         !
         line_=line(i)(1:30)//coord_str
         !
      else
         !
         line_=adjustl(line(i))
         !
      end if
      !
      write(unit=666,fmt='(A)') line_
      !
   end do
   !
end if
!
close(unit=666)

end subroutine write_pdb


end program test_vector_rotation_about_axis
