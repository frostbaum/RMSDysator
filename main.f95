program main
  use structure_c
  use qsort_struc_m
  use rmsd_m
  implicit none
  type tag_list
    integer, dimension(:), allocatable :: i
  end type
  integer :: num_files, i, j, cand, k, cnt, curnum, tmpsize
  character(len=50) :: filename
  integer, dimension(:), allocatable :: tmp_array
  integer, dimension(0:19) :: ref_array, com_array
  type(structure), dimension(:), allocatable :: struc_rep
  type(tag_list), dimension(:), allocatable :: tl, tl_tmp
  double precision :: rmsd, rmsd_tmp, ediff, ediff_tmp, rmsd_lim = .5d0, ediff_lim = 5.d0, e_ub = 0.d0, anceps = 1.d0
  double precision, dimension(:), allocatable :: maxrmsd, maxediff, maxrmsd_tmp, maxediff_tmp
  double precision, dimension(:,:), allocatable :: rmsdmat
  double precision, parameter :: h2kcm = 627.509469d0
  
  open(123,file='tmp.txt',status='old')
  read(123,*) num_files
  read(123,*) rmsd_lim
  read(123,*) anceps
  
  allocate(struc_rep(num_files),rmsdmat(num_files,num_files))
  rmsdmat = 0.d0
  read(123,*) curnum, filename
  call struc_read_file(struc_rep(1),trim(filename),'xyz')
  call struc_set(struc_rep(1),t=curnum)
  call struc_ctrans(struc_rep(1))
  do j = 0, 19
    ref_array(j) = struc_get_natoms(struc_rep(1),j)
  end do
  
  do i = 2, num_files
    read(123,*) curnum, filename
    call struc_read_file(struc_rep(i),trim(filename),'xyz')
    call struc_set(struc_rep(i),t=curnum)
    
    do j = 0, 19
      com_array(j) = struc_get_natoms(struc_rep(i),j)
    end do
    
    if (any(com_array .ne. ref_array)) then
      write(6,'(A)') 'Structure '//trim(filename)//' does not match the reference (atm fatal)'
      stop
    end if
    call struc_ctrans(struc_rep(i))
  end do
  
  close(123)

  call qsort(struc_rep)
  
  call rmsd_set(rmsd_lim,.true.,anceps)
  call rmsd_init(struc_rep(1))
  
  do i = 1, num_files-1
    call rmsd_set_anc_s1(struc_rep(i))
    do j = i+1, num_files
      call rmsd_set_anc_s2(struc_rep(j))
      call rmsd_calc(rmsd,cnt)
      rmsdmat(i,j) = rmsd
      rmsdmat(j,i) = rmsd
    end do
  end do
  
  open(444,file='new.out')
  do i = 1, num_files
    write(444,'(100F8.3)') rmsdmat(:,i)
  end do
  
  close(444)
  
  stop
666 write(*,*) 'bad input'

end program
