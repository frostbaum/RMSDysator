module qsort_struc_m
use structure_c
implicit none
private
public :: qsort

contains

recursive subroutine qsort(struc_rep)
  type(structure), dimension(:) :: struc_rep
  integer :: iq

  if(size(struc_rep) > 1) then
     call Partition(struc_rep,iq)
     call qsort(struc_rep(:iq-1))
     call qsort(struc_rep(iq:))
  end if
end subroutine qsort

subroutine Partition(struc_rep,marker)
  type(structure), dimension(:) :: struc_rep
  type(structure) :: x, temp
  integer, intent(out) :: marker
  integer :: i, j
  
  x = struc_rep(1)
  i= 0
  j= size(struc_rep) + 1

  do
     j = j-1
     do
        if (struc_rep(j)%energy <= x%energy) exit
        j = j-1
     end do
     i = i+1
     do
        if (struc_rep(i)%energy >= x%energy) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = struc_rep(i)
        struc_rep(i) = struc_rep(j)
        struc_rep(j) = temp
     else if (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module
