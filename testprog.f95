
program test
use rmsd_m
use structure_c
implicit none
integer :: i, j, anccnt
integer, dimension(:), allocatable :: idx
type(structure) :: s1,s2
logical :: myfail
double precision :: asd, t1, t2
character(len=4) :: ic, jc
call rmsd_set(0.4d0,.true.,.5d0)
open(111,file='indices_2way.txt')
!~open(6,file='tmp.out')
!~call cpu_time(t1)
!~call struc_read_file(s2,'try.h2o/struc0002.arc','arc')
!~call struc_read_file(s1,'try.h2o/struc0012.arc','arc')
write(6,'(A)') 'Gly.h2o'
write(111,'(A)') 'Gly.h2o'
call struc_read_file(s1,'0337_end.xyz','xyz')
call struc_read_file(s2,'0479_end.xyz','xyz')
call struc_ctrans(s1)
call struc_ctrans(s2)
call rmsd_init(s1)
allocate(idx(struc_get_natoms(s1)))
call rmsd_set_anc_s1(s1)
call rmsd_set_anc_s2(s2)
call rmsd_calc(asd,anccnt,idx)
write(6,'(A,F10.4,I5)') ' RMSD: ', asd, anccnt
write(111,'(50I4)') idx
write(6,*)
write(6,'(A)') 'Hexadien'
write(111,*)
write(111,'(A)') 'Hexadien'

do i = 1, 10
  write(ic,'(I4.4)') i
  call struc_read_file(s1,'hexadien.endFiles/'//ic//'_end.xyz','xyz')
  call struc_ctrans(s1)
  call rmsd_set_anc_s1(s1)
  deallocate(idx)
  allocate(idx(struc_get_natoms(s1)))
  do j = i+1, 11
    write(jc,'(I4.4)') j
    call struc_read_file(s2,'hexadien.endFiles/'//jc//'_end.xyz','xyz')
    call struc_ctrans(s2)
    call rmsd_set_anc_s2(s2)
    call rmsd_calc(asd,anccnt,idx)
    write(6,'(2I3,A,F10.4,I5)') i, j, ' RMSD: ', asd, anccnt
    write(111,'(50I4)') idx
  end do
end do
write(6,*)
write(6,'(A)') 'Allyl.h2o'
write(111,*)
write(111,'(A)') 'Allyl.h2o'
do i = 1, 10
  write(ic,'(I4.4)') i
  call struc_read_file(s1,'allyl.h2o.endFiles/'//ic//'_end.xyz','xyz')
  call struc_ctrans(s1)
  call rmsd_set_anc_s1(s1)
  deallocate(idx)
  allocate(idx(struc_get_natoms(s1)))
  do j = i+1, 11
    write(jc,'(I4.4)') j
    call struc_read_file(s2,'allyl.h2o.endFiles/'//jc//'_end.xyz','xyz')
    call struc_ctrans(s2)
    call rmsd_set_anc_s2(s2)
    call rmsd_calc(asd,anccnt,idx)
    write(6,'(2I3,A,F10.4,I5)') i, j, ' RMSD: ', asd, anccnt
    write(111,'(50I4)') idx
  end do
end do
write(6,*)
write(6,'(A)') 'Try.h2o'
write(111,*)
write(111,'(A)') 'Try.h2o'
do i = 1, 10
  write(ic,'(I4.4)') i
  call struc_read_file(s1,'try.h2o/struc'//ic//'.arc','arc')
  call struc_ctrans(s1)
  call rmsd_set_anc_s1(s1)
  deallocate(idx)
  allocate(idx(struc_get_natoms(s1)))
  do j = i+1, 11
    write(jc,'(I4.4)') j
    call struc_read_file(s2,'try.h2o/struc'//jc//'.arc','arc')
    call struc_ctrans(s2)
    call rmsd_set_anc_s2(s2)
    call rmsd_calc(asd,anccnt,idx)
    write(6,'(2I3,A,F10.4,I5)') i, j, ' RMSD: ', asd, anccnt
    write(111,'(50I4)') idx
  end do
end do
close(111)
!~close(6)
!~do i = 1, 10
!~am_eps = i*5.d-2
!~write(*,'(A,F8.2)') 'eps:  ', am_eps
!~call rmsd_calc(s2,asd)
!~end do
!~do i = 6, 20
!~am_eps = i*1.d-1
!~write(*,'(A,F8.2)') 'eps:  ', am_eps
!~call rmsd_calc(s2,asd)
!~end do


!~write(*,*) asd,myfail
!~call cpu_time(t2)
!~write(*,*) asd
!~write(*,*) 'time', t2-t1
!~call struc_write_file(s1,'test1.xyz','xyz')
!~call struc_write_file(s2,'test2.xyz','xyz')

end program
