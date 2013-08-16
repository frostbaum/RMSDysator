program main
  use structure_c
  use qsort_struc_m
  use rmsd_m
  implicit none
  type tag_list
    integer, dimension(:), allocatable :: i
  end type
  integer :: num_files, i, j, cand, k, cnt
  character(len=20) :: filename, num_buf
  character(len=4) :: numchar
  integer, dimension(:), allocatable :: tmp_array
  type(structure), dimension(:), allocatable :: struc_rep
  type(tag_list), dimension(:), allocatable :: tl, tl_tmp
  double precision :: rmsd, rmsd_tmp, rmsd_max = .5d0, ediff_max = .1d0
  double precision :: e_ub =0.d0
  double precision, parameter :: h2kcm = 627.509469d0
  
  num_files = 200
  rmsd_max = 0.1d0
  !~open(99,file='a.txt',status='old')
  !~  read(99,*) num_files
  !~  read(99,*) e_ub
  !~  read(99,*) rmsd_max, ediff_max
  !~close(99)
  !~allocate(struc_rep(2))
  !~call io_def(current_file,'allyl.h2o.endFiles/0267_end.xyz','xyz','r')
  !~struc_rep(1) = current_file
  !~struc_rep(1)%tag = 180
  !~call struc_cut_h(struc_rep(1))
  !~call centr_trans(struc_rep(1))
  !~call io_def(current_file,'allyl.h2o.endFiles/0180_end.xyz','xyz','r')
  !~struc_rep(2) = current_file
  !~struc_rep(2)%tag = 267
  !~call struc_cut_h(struc_rep(2))
  !~call centr_trans(struc_rep(2))
  !~call def_struc_pair(current_pair,struc_rep(2),struc_rep(1))
  !~call calc_min_rmsd(current_pair,rmsd_tmp)
  !~write(*,*) rmsd_tmp
  !~deallocate(struc_rep)
  allocate(struc_rep(num_files))
  do i=1,num_files
    write(numchar,'(I4.4)') i
    call struc_read_file(struc_rep(i),'hexadien.endFiles/'//numchar//'_end.xyz','xyz')
    struc_rep(i)%tag = i
    call struc_ctrans(struc_rep(i))
  end do
  !~do i=1,num_files
  !~  write(*,*) struc_rep(i)%tag, struc_rep(i)%energy
  !~end do
  call qsort(struc_rep)
  call rmsd_set(.1d0,.true.,1.2d0)
  call rmsd_init(struc_rep(1))
  !~write(*,*) '****************************************************'
  !~do i=1,num_files
  !~  write(*,*) struc_rep(i)%tag, struc_rep(i)%energy
  !~end do
  open(95,file='rmsd.txt')
  allocate(tl(1))
  allocate(tl(1)%i(1))
  tl(1)%i(1) = 1
  do i=2,num_files
    rmsd = 1111.d0
    cand = 0
    call rmsd_set_anc_s1(struc_rep(i))
    !~write(*,*)
    do j=1,size(tl)
      if (abs(struc_get_energy(struc_rep(i)) - struc_get_energy(struc_rep(tl(j)%i(1))))*h2kcm .gt. ediff_max) cycle
      call rmsd_set_anc_s2(struc_rep(tl(j)%i(1)))
      call rmsd_calc(rmsd_tmp,cnt)
      write(95,'(A,I4.1,A,I4.1,A,F10.4,A,F10.4,A,I4)') 'compare ',struc_rep(i)%tag,' and ',struc_rep(tl(j)%i(1))%tag,&
      &'; RMSD [A]: ',rmsd_tmp,';    Delta E [kcal/mol]: ',(struc_rep(tl(j)%i(1))%energy-struc_rep(i)%energy)*h2kcm, ' cnt:', cnt
      if (rmsd_tmp .ge. rmsd) cycle
      cand = j
      rmsd = rmsd_tmp
    end do
    if (rmsd .le. rmsd_max) then
      !~write(*,'(I4.1,A,I4.1)') i, ' is higher in energy than ', tl(cand)%i(1)
      allocate(tmp_array(size(tl(cand)%i)))
      tmp_array = tl(cand)%i
      deallocate(tl(cand)%i)
      allocate(tl(cand)%i(size(tmp_array)+1))
      tl(cand)%i(:size(tmp_array)) = tmp_array
      tl(cand)%i(size(tmp_array)+1) = i
      deallocate(tmp_array)
    else
      !~write(*,'(I4.1,A)') i, ' is a new structure'
      allocate(tl_tmp(size(tl)))
      do k=1,size(tl)
        allocate(tl_tmp(k)%i(size(tl(k)%i)))
        tl_tmp(k)%i = tl(k)%i
      end do
      deallocate(tl)
      allocate(tl(size(tl_tmp)+1))
      do k=1,size(tl_tmp)
        allocate(tl(k)%i(size(tl_tmp(k)%i)))
        tl(k)%i = tl_tmp(k)%i
      end do
      allocate(tl(size(tl_tmp)+1)%i(1))
      tl(size(tl_tmp)+1)%i(1) = i
      deallocate(tl_tmp)
    end if
  end do
  close(95)
  open(111,file='out.txt',status='replace')
  !~write(111,*) '#energy,', 'occurrence,', 'matching structures with increasing energy'
  do i=1,size(tl)
    write(111,'(F16.8,I4,999I5)') (struc_rep(tl(i)%i(1))%energy-e_ub)*h2kcm, size(tl(i)%i), struc_rep(tl(i)%i(:))%tag
  end do
  close(111)
  deallocate(struc_rep)
  stop
666 write(*,*) 'bad input'

end program
