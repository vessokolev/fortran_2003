module mod_pdb
!
! A simple module for reading and writing PDB files. Note
! that only the ATOM and HEATAM sections are handling
! during the execution of the suboutines given bellow.
!
! Author : Veselin Kolev <vesso.kolev@gmail.com>
! Version: 2019012000
! License: GPLv2
!
use iso_c_binding,only:C_INT,C_FLOAT
!
implicit none


type resid_t
   !
   logical        :: do_count
   integer(C_INT) :: restype
   integer(C_INT) :: serial_l
   integer(C_INT) :: serial_u
   integer(C_INT) :: indx(5)
   real(C_FLOAT)  :: axis_1(3)
   real(C_FLOAT)  :: axis_2(3)
   real(C_FLOAT)  :: angle_1
   real(C_FLOAT)  :: angle_2
   real(C_FLOAT)  :: matrix_1(3,3)
   real(C_FLOAT)  :: matrix_2(3,3)
   logical        :: rot_1_possible
   logical        :: rot_2_possible
   !
   ! indx(1) - N  serial
   ! indx(2) - C  serial
   ! indx(3) - CA serial
   ! indx(4) - HA serial
   ! indx(5) - CB serial
   !
end type resid_t


contains


subroutine count_pdb_lines(pdb_file,num_lines)
!
! Estimates the number of lines detected inside the PDB file.
!
! Interface variables:
!
character(len=*),intent(in) :: pdb_file
integer(C_INT),intent(out)  :: num_lines
!
! Local variables:
!
integer(C_INT)              :: stat
integer(C_INT),parameter    :: str_len=80
character(len=str_len)      :: buffer

open(unit=666,file=pdb_file,iostat=stat)
!
num_lines=0
!
do while(stat.eq.0)
   !
   read(unit=666,fmt='(A80)',iostat=stat) buffer
   !
   if (stat==0) then
      !
      num_lines=num_lines+1
      !
   end if
   !
end do
!
close(unit=666)

end subroutine count_pdb_lines


subroutine read_pdb_lines(pdb_file,pdb_line)
!
! Reads lines detected inside the PDB file.
!
! Interface variables:
!
character(len=*),intent(in)   :: pdb_file
character(len=80),intent(out) :: pdb_line(:)
!
! Local variables:
!
integer(C_INT)                :: i
integer(C_INT)                :: stat
integer(C_INT),parameter      :: str_len=80
character(len=str_len)        :: buffer

open(unit=666,file=pdb_file,iostat=stat)
!
i=0
!
do while(stat.eq.0)
   !
   read(unit=666,fmt='(A80)',iostat=stat) buffer
   !
   if (stat==0) then
      !
      i=i+1
      !
      pdb_line(i)=buffer
      !
   end if
   !
end do
!
close(unit=666)

end subroutine read_pdb_lines


subroutine count_pdb_atoms(pdb_line,num_atoms)
!
! Reads the number of ATOM or HETATM records declared
! in the PDB file. Single model hosted inside the PDB
! file is assumed.
!
! Interface variables:
!
character(len=80),intent(in) :: pdb_line(:)
integer(C_INT),intent(out)   :: num_atoms
!
! Local variables:
!
integer(C_INT)              :: i
integer(C_INT),parameter    :: str_len=54
character(len=str_len)      :: buffer

num_atoms=0
!
do i=1,ubound(pdb_line,1)
   !
   if (pdb_line(i)(1:6)=='ATOM  '.or.&
       pdb_line(i)(1:6)=='HETATM') then
      !
      num_atoms=num_atoms+1
      !
   end if
   !
end do

end subroutine count_pdb_atoms


subroutine read_pdb_coord(pdb_line,coord)
!
! Reads the atomic coordinates declared inside the ATOM and HETATM
! records in PDB structure. Note, that one need to call first the
! subroutine "s_count_pdb_atoms" to estimate the number of those
! records.
!
! Interface variables:
!
character(len=80),intent(in) :: pdb_line(:)
real(C_FLOAT),intent(out)    :: coord(:,:)
!
! Local variables:
!
integer(C_INT)               :: i
integer(C_INT),parameter     :: str_len=54
character(len=str_len)       :: buffer

do i=1,ubound(pdb_line,1)
   !
   if (pdb_line(i)(1:6)=='ATOM  '.or.&
       pdb_line(i)(1:6)=='HETATM') then
      !
      read(pdb_line(i)(31:54),*) coord(:,i)
      !
   end if
   !
end do

end subroutine read_pdb_coord


subroutine read_pdb_num_resid(pdb_line,list,num_resid)
!
! Interface variab
!
character(len=80),intent(in) :: pdb_line(:)
character(len=3),intent(in)  :: list(:)
integer(C_INT),intent(out)   :: num_resid
!
! Local variables:
!
integer(C_INT)               :: i
integer(C_INT)               :: j
integer(C_INT)               :: k

do i=1,ubound(pdb_line,1)
   !
   if (pdb_line(i)(1:6)=='ATOM  '.or.&
       pdb_line(i)(1:6)=='HETATM') then
      !
      if (i==1) then
         !
         num_resid=1
         !
         read(pdb_line(i)(23:26),*) j
         !
      else
         !
         read(pdb_line(i)(23:26),*) k
         !
         if (k/=j) then
            !
            j=k
            !
            num_resid=num_resid+1
            !
         end if
         !
      end if
      !
   end if
   !
end do

end subroutine read_pdb_num_resid


subroutine read_pdb_resid(pdb_line,list,resid)
!
! Interface variab
!
character(len=80),intent(in) :: pdb_line(:)
character(len=3),intent(in)  :: list(:)
type(resid_t),intent(out)    :: resid(:)
!
! Local variables:
!
integer(C_INT)               :: i
integer(C_INT)               :: j
integer(C_INT)               :: k
integer(C_INT)               :: l
integer(C_INT)               :: m
logical                      :: flag

m=0
!
do i=1,ubound(pdb_line,1)
   !
   if (pdb_line(i)(1:6)=='ATOM  '.or.&
       pdb_line(i)(1:6)=='HETATM') then
      !
      m=m+1
      !
      if (i==1) then
         !
         l=1
         !
         read(pdb_line(i)(23:26),*) j
         !
         flag=.true.
         !
      else
         !
         read(pdb_line(i)(23:26),*) k
         !
         if (k/=j) then
            !
            j=k
            !
            resid(l)%serial_u=m-1
            !
            l=l+1
            !
            flag=.true.
            !
         end if
         !
      end if
      !
      if (flag) then
         !
         resid(l)%serial_l=m
         !
         resid(l)%indx=0
         !
         if (resname_in_resname_list(pdb_line(i)(18:20),&
             list,resid(l)%restype)) then
            !
            resid(l)%do_count=.true.
            !
         else
            !
            resid(l)%do_count=.false.
            !
         end if
         !
         flag=.false.
         !
      end if
      !
      if (trim(adjustl(pdb_line(i)(13:16)))=='N') then
         !
         resid(l)%indx(1)=m
         !
      end if
      !
      if (trim(adjustl(pdb_line(i)(13:16)))=='C') then
         !
         resid(l)%indx(2)=m
         !
      end if
      !
      if (trim(adjustl(pdb_line(i)(13:16)))=='CA') then
         !
         resid(l)%indx(3)=m
         !
      end if
      !
      if (trim(adjustl(pdb_line(i)(13:16)))=='HA') then
         !
         resid(l)%indx(4)=m
         !
      end if
      !
      if (trim(adjustl(pdb_line(i)(13:16)))=='CB') then
         !
         resid(l)%indx(5)=m
         !
      end if
      !
   end if
   !
end do
!
resid(l)%serial_u=m

end subroutine read_pdb_resid


function resname_in_resname_list(resname,list,indx) result(res)
character(len=3),intent(in) :: resname
character(len=3),intent(in) :: list(:)
integer(C_INT),intent(out)  :: indx
logical                     :: res
!
! Local variables:
!
integer(C_INT)              :: i

res=.false.
!
indx=0
!
do i=1,ubound(list,1)
   !
   if (resname==list(i)) then
      !
      res=.true.
      !
      indx=i
      !
      exit
      !
   end if
   !
end do

end function resname_in_resname_list


subroutine write_pdb_file(pdb_file,pdb_line,coord)
!
! Writes the atomic ATOM and HETATM records to a new PDB file.
!
! Interface variables:
!
character(len=*),intent(in)  :: pdb_file
character(len=80),intent(in) :: pdb_line(:)
real(C_FLOAT),intent(in)     :: coord(:,:)
!
! Local variables:
!
integer(C_INT)               :: i
integer(C_INT)               :: stat
integer(C_INT),parameter     :: str_len=54

open(unit=666,file=pdb_file,status='unknown',iostat=stat)
!
do i=1,ubound(pdb_line,1)
   !
   if (pdb_line(i)(1:6)=='ATOM  '.or.pdb_line(i)(1:6)=='HETATM') then
      !
      write(unit=666,fmt='(A30,3F8.3)') pdb_line(i)(1:30),coord(:,i)
      !
   else
      !
      write(unit=666,fmt='(A)') trim(adjustl(pdb_line(i)))
      !
   end if
   !
end do
!
close(unit=666)

end subroutine write_pdb_file


end module mod_pdb
