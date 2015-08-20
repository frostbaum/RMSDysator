program main
  use structure_c
  use qsort_struc_m
  use rmsd_m
  implicit none
  type tag_list
    integer, dimension(:), allocatable :: i
  end type
  
  integer :: num_files, i, j, cand, k, cnt, curnum, tmpsize, tmpratio, curratio
  character(len=50) :: filename
  integer, dimension(:), allocatable :: tmp_array
  integer, dimension(0:19) :: ref_array, com_array
  type(structure), dimension(:), allocatable :: struc_rep
  type(tag_list), dimension(:), allocatable :: tl, tl_tmp
  double precision :: rmsd, rmsd_tmp, ediff, ediff_tmp, rmsd_lim = .5d0, ediff_lim = 5.d0, e_ub = 0.d0, anceps = 1.d0
  double precision, dimension(:), allocatable :: maxrmsd, maxediff, maxrmsd_tmp, maxediff_tmp
  double precision :: econv = 627.509469d0
  
  open(123,file='tmp.txt',status='old')
  read(123,*) num_files
  read(123,*) rmsd_lim, ediff_lim, econv
  read(123,*) anceps
  
  write(6,'(A)') 'Reading structures ...'
  
  allocate(struc_rep(num_files))
  
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
  e_ub = struc_get_energy(struc_rep(1))
  
  open(95,file='rmsd.txt')
  
  allocate(tl(1),maxrmsd(1),maxediff(1))
  allocate(tl(1)%i(1))
  tl(1)%i(1) = 1
  maxrmsd(1) = 0.d0
  maxediff(1) = 0.d0


  write(6,'(A)') 'Starting analysis ...'
  
  tmpratio=0
  write(6,'(A)') '---|---|---|---|---|---|---|---|---|---|'
  
  do i=2,num_files
    rmsd = 1111.d0
    cand = 0
    call rmsd_set_anc_s1(struc_rep(i))
    !~write(*,*)
	curratio = i*40/num_files
	call progressbar_add(tmpratio,curratio)
	tmpratio = curratio
    do j=1,size(tl)
      ediff_tmp = abs(struc_get_energy(struc_rep(i)) - struc_get_energy(struc_rep(tl(j)%i(1))))*econv
      if (ediff_tmp .gt. ediff_lim) cycle
      
      call rmsd_set_anc_s2(struc_rep(tl(j)%i(1)))
      call rmsd_calc(rmsd_tmp,cnt)
      
      write(95,'(A,I4.1,A,I4.1,A,F10.4,A,F10.4,A,I4)') 'compare ',struc_rep(i)%tag,' and ',struc_rep(tl(j)%i(1))%tag,&
      &'; RMSD [A]: ',rmsd_tmp,';    Delta E [kcal/mol]: ',(struc_rep(tl(j)%i(1))%energy-struc_rep(i)%energy)*econv,&
      &';  Match count:', cnt
      
      if (rmsd_tmp .ge. rmsd) cycle
      cand = j
      rmsd = rmsd_tmp
      ediff = ediff_tmp
    end do
    if (rmsd .le. rmsd_lim) then
      !~write(*,'(I4.1,A,I4.1)') i, ' is higher in energy than ', tl(cand)%i(1)
      allocate(tmp_array(size(tl(cand)%i)))
      tmp_array = tl(cand)%i
      deallocate(tl(cand)%i)
      allocate(tl(cand)%i(size(tmp_array)+1))
      tl(cand)%i(:size(tmp_array)) = tmp_array
      tl(cand)%i(size(tmp_array)+1) = i
      deallocate(tmp_array)
      
      maxrmsd(cand) = max(maxrmsd(cand),rmsd)
      maxediff(cand) = max(maxediff(cand),ediff)
    else
      !~write(*,'(I4.1,A)') i, ' is a new structure'
      tmpsize = size(tl)
      allocate(tl_tmp(tmpsize),maxediff_tmp(tmpsize),maxrmsd_tmp(tmpsize))
      maxediff_tmp = maxediff
      maxrmsd_tmp = maxrmsd
      
      do k=1,tmpsize
        allocate(tl_tmp(k)%i(size(tl(k)%i)))
        tl_tmp(k)%i = tl(k)%i
      end do
      
      deallocate(tl,maxediff,maxrmsd)
      allocate(tl(tmpsize+1),maxediff(tmpsize+1),maxrmsd(tmpsize+1))
      maxediff(:tmpsize) = maxediff_tmp
      maxrmsd(:tmpsize) = maxrmsd_tmp
      
      do k=1,tmpsize
        allocate(tl(k)%i(size(tl_tmp(k)%i)))
        tl(k)%i = tl_tmp(k)%i
      end do
      
      allocate(tl(tmpsize+1)%i(1))
      tl(tmpsize+1)%i(1) = i
      maxediff(tmpsize+1) = 0.d0
      maxrmsd(tmpsize+1) = 0.d0
      
      deallocate(tl_tmp,maxediff_tmp,maxrmsd_tmp)
    end if
  end do
  close(95)
  write(6,'(A)')
  
  open(111,file='out.txt',status='replace')
  !~write(111,*) '#energy,', 'occurrence,', 'matching structures with increasing energy'
  do i=1,size(tl)
    write(111,'(F16.5,I4,F10.5,F10.5,9999I5)') (struc_rep(tl(i)%i(1))%energy-e_ub)*econv, size(tl(i)%i), maxrmsd(i), maxediff(i),&
    & struc_rep(tl(i)%i(:))%tag
  end do
  close(111)
  deallocate(struc_rep)
  stop
666 write(*,*) 'bad input'

contains
  subroutine progressbar_add(old,current)
    integer :: current, old, loopcnt
	
	do loopcnt=old+1,current
	  write(6,'(A)',advance='no') '*'
	end do
	
  end subroutine

end program
