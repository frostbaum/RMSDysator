module rmsd_m
  use v3d_func_rep
  use structure_c
  use mk
  implicit none
  save
  private
  public :: rmsd_init, rmsd_set_anc_s1, rmsd_set_anc_s2, rmsd_calc, rmsd_set
  
  type anchor
    integer, dimension(3) :: i, t
    double precision, dimension(3) :: d
  end type
  
  type(structure), pointer :: s1 => null(), s2 => null()
  integer :: natoms = 0
  double precision :: am_eps = .4d0, rmsdlim = .5d0
  logical :: mirror = .true.
  
  type(anchor) :: anc_s1, anc_s2
  double precision, dimension(:,:), allocatable :: s1_coords, s2_coords
  integer, dimension(:), allocatable :: atypes, ntypes
  
contains
  subroutine rmsd_set(eps,mirrorflag,rmsdlimit)
    double precision :: eps, rmsdlimit
    logical :: mirrorflag
    
    am_eps = eps
    mirror = mirrorflag
    rmsdlim = rmsdlimit
  end subroutine
  
  subroutine rmsd_init(samplestruc)
    type(structure) :: samplestruc
    integer :: natyp, atyptmp, i
    
    if (natoms .ne. 0) then
      call rmsd_clr()
    end if
    
    natoms = struc_get_natoms(samplestruc)
    allocate(s1_coords(3,natoms),s2_coords(3,natoms))
    
    natyp = samplestruc%lbtype(0)%i(0)
    allocate(atypes(natyp),ntypes(0:natyp))
    
    ntypes(0) = natyp
    do i = 1, natyp
      atyptmp = samplestruc%lbtype(0)%i(i)
      atypes(i) = atyptmp
      ntypes(i) = samplestruc%lbtype(atyptmp)%i(0)
    end do
  end subroutine
  
  subroutine rmsd_clr()
    deallocate(s1_coords,s2_coords,atypes,ntypes)
    natoms = 0
  end subroutine
  
  subroutine rmsd_set_anc_s1(s1_in)
    type(structure), target :: s1_in
    
    if (struc_get_natoms(s1_in) .ne. natoms) then
      call rmsd_init(s1_in)
    end if
    
    s1 => s1_in
    
    call anc_fix(s1,anc_s1,s1_coords)
    
  end subroutine
  
  subroutine rmsd_set_anc_s2(s2_in)
    type(structure), target :: s2_in
    
    if (struc_get_natoms(s2_in) .ne. natoms) then
      call rmsd_init(s2_in)
    end if
    
    s2 => s2_in
    
    call anc_fix(s2,anc_s2,s2_coords)
    
  end subroutine
  
  subroutine rmsd_calc(rmsd,anccnt,idx)
    integer :: anccnt, tmpcnt1, tmpcnt2
    integer, dimension(natoms) :: idx1, idx2, idx
    double precision :: rmsd, rmsdtmp1, rmsdtmp2
    
    call rmsd_match_anc(s2,s1,anc_s1,s1_coords,rmsdtmp1,tmpcnt1,idx1)
    
    call rmsd_match_anc(s1,s2,anc_s2,s2_coords,rmsdtmp2,tmpcnt2,idx2)
    
    anccnt = tmpcnt1 + tmpcnt2
    if (rmsdtmp1 .lt. rmsdtmp2) then
      rmsd = rmsdtmp1
      idx = idx1
    else
      rmsd = rmsdtmp2
      idx = idx2
    end if
  end subroutine
  
  subroutine rmsd_match_anc(tgt,mdl,anc_mdl,mdl_coords,rmsd,anccnt,idx)
    type(structure) :: mdl, tgt
    type(anchor) :: anc_mdl
    double precision, dimension(3,natoms) :: mdl_coords
    double precision :: rmsd, rmsdtmp
    integer, dimension(natoms) :: idx, tmpidx
    integer :: i, j, k, ii, jj, kk, anccnt
    double precision :: tmpd1, tmpd2
    
    rmsd = 777.d0
    idx = 0
    anccnt = 0
    do i = 1, tgt%lbtype(anc_mdl%t(1))%i(0)
      ii = tgt%lbtype(anc_mdl%t(1))%i(i)
      do j = 1, tgt%lbtype(anc_mdl%t(2))%i(0)
        jj = tgt%lbtype(anc_mdl%t(2))%i(j)
        
        tmpd1 = get_dist(struc_get_coords(tgt,ii),struc_get_coords(tgt,jj))
        if (abs(tmpd1-anc_mdl%d(1)) .gt. am_eps) cycle
        
        do k = 1, tgt%lbtype(anc_mdl%t(3))%i(0)
          kk = tgt%lbtype(anc_mdl%t(3))%i(k)
            
          tmpd1 = get_dist(struc_get_coords(tgt,ii),struc_get_coords(tgt,kk))
          tmpd2 = get_dist(struc_get_coords(tgt,jj),struc_get_coords(tgt,kk))
          if (abs(tmpd1-anc_mdl%d(2)) .gt. am_eps .or. abs(tmpd2-anc_mdl%d(3)) .gt. am_eps) cycle
          
          call superpose(tgt,(/ ii, jj, kk /),mdl,mdl_coords,rmsdtmp,tmpidx)
            
          if (rmsdtmp .lt. rmsd) then
            rmsd = rmsdtmp
            idx = tmpidx
          end if
          anccnt = anccnt + 1
        end do
      end do
    end do
  end subroutine

  subroutine superpose(tgt,anc_tgt_i,mdl,mdl_coords,rmsd,idx)
    type(structure) :: mdl, tgt
    integer, dimension(3) :: anc_tgt_i
    integer, dimension(natoms) :: idx1, idx2, idx
    double precision :: rmsd1, rmsd2, rmsd
    integer :: i
    double precision, dimension(3,natoms) :: mdl_coords, mdl_tcrd, tgt_tcrd
    
    call anc_superpose(tgt,anc_tgt_i,tgt_tcrd)
    call assign_atoms(mdl,mdl_coords,tgt,mdl_tcrd,tgt_tcrd,idx1)
    call rmsd_quat(mdl_tcrd,tgt_tcrd,rmsd1)
    do i = 1, natoms
      mdl_coords(:,i) = get_refl_plane_n((/0.d0,0.d0,1.d0/),mdl_coords(:,i))
    end do
    call assign_atoms(mdl,mdl_coords,tgt,mdl_tcrd,tgt_tcrd,idx2)
    call rmsd_quat(mdl_tcrd,tgt_tcrd,rmsd2)
    !~rmsd = min(rmsd1,rmsd2)
    if (rmsd1 .lt. rmsd2) then
      rmsd = rmsd1
      idx = idx1
    else
      rmsd = rmsd2
      idx = idx2
    end if
  end subroutine

  subroutine anc_superpose(tgt,anc_tgt_i,tgt_coordstmp)
    type(structure) :: tgt
    integer, dimension(3) :: anc_tgt_i
    double precision, dimension(3) :: tx, ty, vec
    double precision, dimension(3,3) :: tmat_tgt
    double precision, dimension(3,natoms) :: tgt_coordstmp
    integer :: i
    
    tx(:) = get_normal(struc_get_coords(tgt,anc_tgt_i(2)) - struc_get_coords(tgt,anc_tgt_i(1)))
    vec(:) = struc_get_coords(tgt,anc_tgt_i(3)) - struc_get_coords(tgt,anc_tgt_i(1))
    ty(:) = get_normal(vec(:) - get_sc_prod(tx,vec)*tx(:))

    tmat_tgt(1,:) = tx
    tmat_tgt(2,:) = ty
    tmat_tgt(3,:) = get_cr_prod(tx,ty)

    do i = 1, natoms
      tgt_coordstmp(:,i) = get_lin_map(tmat_tgt,struc_get_coords(tgt,i))
    end do
    
  end subroutine
  
  subroutine assign_atoms(mdl,mdl_coords,tgt,mdl_coordstmp,tgt_coordstmp,mdl_idx)
    type(structure) :: mdl, tgt
    integer :: tn, i, j, k, nttn
    integer, dimension(natoms) :: tmpsol, mdl_idx
    double precision, dimension(natoms*natoms) :: distances
    double precision, dimension(3,natoms) :: mdl_coords, mdl_coordstmp, tgt_coordstmp
    
    do tn = 1, ntypes(0)
      if (ntypes(tn) .eq. 1) then
        mdl_idx(mdl%lbtype(atypes(tn))%i(1)) = tgt%lbtype(atypes(tn))%i(1)
      else
        nttn = ntypes(tn)
        k = 0
        do i = 1, nttn
          do j = 1, nttn
            k = k + 1
            distances(k) = get_dist(mdl_coords(:,mdl%lbtype(atypes(tn))%i(j)),tgt_coordstmp(:,tgt%lbtype(atypes(tn))%i(i)))
          end do
        end do
        
        call mk_ass(nttn,reshape(distances,(/nttn,nttn/)),tmpsol(1:nttn))
        
        do i = 1, nttn
          mdl_idx(mdl%lbtype(atypes(tn))%i(i)) = tgt%lbtype(atypes(tn))%i(tmpsol(i))
        end do
        
      end if
    end do
    
    do i = 1, natoms
      mdl_coordstmp(:,mdl_idx(i)) = mdl_coords(:,i)
      !~write(*,*) mdl_idx(i)
    end do
  end subroutine
  
  subroutine rmsd_quat(crd1,crd2,rmsd)
    double precision, dimension(3,natoms) :: crd1, crd2
    double precision, dimension(3,3) :: comat
    double precision, dimension(4,4) :: fmat
    double precision, dimension(4) :: ew
    double precision, dimension(136) :: work
    double precision :: rmsd
    integer :: k, info
    external dsyev
    
    comat = calc_comat(crd1,crd2)
    fmat = calc_fmat(comat)
    
    call dsyev('N','L',4,fmat,4,ew,work,136,info)
    
    if (info .eq. 0) then
      rmsd = 0.d0
      do k = 1, natoms
        rmsd = rmsd + get_nrm_sq(crd1(:,k)) + get_nrm_sq(crd2(:,k))
      end do
      
      rmsd = dsqrt(abs(rmsd - 2*ew(4))/natoms)
    else
      rmsd = 666.d0
    end if
  end subroutine
  
  function calc_comat(crd1,crd2) result(mat)
    integer :: i, j, k
    double precision, dimension(3,3) :: mat
    double precision, dimension(3,natoms) :: crd1, crd2
    
    do i = 1, 3
      do j = 1, 3
        mat(i,j) = 0.d0
        do k = 1, natoms
          mat(i,j) = mat(i,j) + crd1(i,k)*crd2(j,k)
        end do
      end do
    end do
  end function
  
  function calc_fmat(comat) result(fmat)
    double precision, dimension(3,3) :: comat
    double precision, dimension(4,4) :: fmat
    
    fmat(1,1) = comat(1,1) + comat(2,2) + comat(3,3)
    fmat(2,1) = comat(2,3) - comat(3,2)
    fmat(3,1) = comat(3,1) - comat(1,3)
    fmat(4,1) = comat(1,2) - comat(2,1)
    fmat(2,2) = comat(1,1) - comat(2,2) - comat(3,3)
    fmat(3,2) = comat(1,2) + comat(2,1)
    fmat(4,2) = comat(1,3) + comat(3,1)
    fmat(3,3) =-comat(1,1) + comat(2,2) - comat(3,3)
    fmat(4,3) = comat(2,3) + comat(3,2)
    fmat(4,4) =-comat(1,1) - comat(2,2) + comat(3,3)
    
  end function

  subroutine anc_fix(struc,anc,coords)
    type(structure) :: struc
    type(anchor) :: anc
    double precision, dimension(3,natoms) :: coords
    
    integer :: i, j
    double precision :: dist, tmpdist
    double precision, dimension(3) :: mc1, mc2, mc3, mx, my, vec
    double precision, dimension(3,3) :: tmat
    
    dist = 0.d0
    
    do i = 1, natoms-1
      do j = i+1, natoms
        tmpdist = get_dist(struc_get_coords(struc,i),struc_get_coords(struc,j))
        
        if (tmpdist .le. dist) cycle
        
        dist = tmpdist
        anc%i(1) = i
        anc%i(2) = j
      end do
    end do
    
    dist = 0.d0
    
    do i = 1, natoms
      tmpdist = get_dist_to_line(struc_get_coords(struc,anc%i(1)),struc_get_coords(struc,anc%i(2)),struc_get_coords(struc,i))
      
      if (tmpdist .le. dist) cycle
      
      dist = tmpdist
      anc%i(3) = i
    end do
    
    do i = 1, 3
      anc%t(i) = atlabel(struc_get_atmtype(struc,anc%i(i)))
    end do
    
    mc1 = struc_get_coords(struc,anc%i(1))
    mc2 = struc_get_coords(struc,anc%i(2))
    mc3 = struc_get_coords(struc,anc%i(3))
    
    anc%d(1) = get_dist(mc1,mc2)
    anc%d(2) = get_dist(mc1,mc3)
    anc%d(3) = get_dist(mc3,mc2)
    
    mx = get_normal(mc2 - mc1)
    vec(:) = mc3 - mc1
    my = get_normal(vec(:) - get_sc_prod(mx,vec)*mx(:))
    
    tmat(1,:) = mx
    tmat(2,:) = my
    tmat(3,:) = get_cr_prod(mx,my)
    
    do i = 1, natoms
      coords(:,i) = get_lin_map(tmat,struc_get_coords(struc,i))
    end do
  end subroutine

end module


