module globalvars

  integer,  parameter    :: dp=selected_real_kind(15,307)
  integer, parameter     :: num_threads=2
  real(dp), parameter    :: v_cutoff = 1.0d-15

  ! 0 for gamess <- currently obsolete
  ! 1 (or anything else really) for psi4
!  integer,  parameter     :: readprogram=2
  ! eigenvalue threshold for linear dependence
!  real(dp), parameter     :: threshd=0.1_dp
  ! keep at 0 for symm_orth
!  integer,  parameter     :: orthomethod=0

  real(dp), allocatable   :: h(:), t(:), s(:), rij(:), x(:), y(:), z(:)
  ! intent(in)
  real(dp), allocatable   :: eigvals(:)
  real(dp), allocatable   :: hmat(:,:), tmat(:,:), smat(:,:)
  real(dp), allocatable   :: xx(:,:), h2(:,:), t2(:,:)
  real(dp), allocatable   :: xint(:,:,:), x2(:,:,:)
  real(dp), allocatable   :: veeao(:,:,:,:)
  real(dp), allocatable   :: chmveeao(:,:,:,:)
  integer,  allocatable   :: rijindex(:,:)
  integer                 :: nrprim, nint2e, numel, numpts
  integer                 :: ldep
  integer                 :: datainp, aoints
  real(dp)                :: wtime
  character(50)           :: outform='(e23.16)'
  character(50)           :: outform2='(2(e23.16,2x))'
  character(1)            :: grdnum
  integer :: nrensp
  logical :: write_all
  real(dp) :: threshd
! character(50) :: outform='(e18.12)'

end module globalvars

module write2str
  contains
    pure integer function str_ilen(i) result(sz)
      ! Returns the length of the string representation of 'i'
      !> gives the length of an integer string
      integer, intent(in) :: i
      integer, parameter :: str_max = 100
      character(str_max) :: s
      !> if you give a string that has 100 characters, it's not gonna happen
      write(s, '(i0)') i
      sz = len_trim(s)
    end function

    pure function stri(i) result(s)
!    elemental function int2str(s, i)
      ! Converts integer "i" to string
      integer, intent(in) :: i
      character(len=str_ilen(i)) :: s

      write(s,'(i0)') i

    end function
end module

program generate_hf
  use globalvars
  use write2str
  use omp_lib

  implicit none
  integer :: i,j,k,ios
  character(256) :: iom

  open(20,file='threshd.dat',status='old')
  read(20,*) threshd
  close(20)
  write(*,'("Linear dependence threshold: "es12.5)') threshd

  write(*,*) ''
  write(*,*) "Now orthogonalizing the AO integrals. Stand by."

  open(newunit=datainp, file="data.inp", status='old', iostat=ios, iomsg=iom)
  if (ios /= 0) then
    write(*,*) 'Fatal error! could not open data.inp!!'
    write(*,*) 'Iomessage: ', trim(iom)
    stop
  endif

  read(datainp,*) nrprim
  write(*,*) "Found ", nrprim, " basis functions!!"
  read(datainp,*) nint2e
  write(*,*) "Found ", nint2e, " non-zero two-electron integrals!!"
  read(datainp,*) nrensp
  write(*,*) "Found ", nrensp, " grid points!!!!!"
!  if (orthomethod == 1) read(datainp,*) numpts
  close(datainp)

  numel = nrprim*(nrprim + 1)/2
  write(*,*) numel, " elements in the one-electron operators!!!!"
  write(*,*) '***********************************************'
  write(*,*) ''
! call omp_set_num_threads(num_threads)
! write(*,*) 'Running with ',num_threads,'thread(s)'

  ! initiate array sequence
  call charge_arrays()
  wtime = omp_get_wtime()
  ! read in integrals from file then fill up integral matrices

  do i=1, nrensp
!    call int2str(grdnum, i)
    grdnum = stri(i)

    if (i == 1) then
      write_all = .true.
    else
      write_all = .false.
    endif

    ! read integrals from psi4
    call readin_psi4()

    ! use static basis set
    call symm_ortho()
    call mult()

    ! write out integrals in compact form
    call write_reduced_out()


  enddo

  ! write wall time
  wtime = omp_get_wtime() - wtime
  write(*,'("Orthogonalization complete. Elapsed time(s): "(f10.2))') wtime

  ! discharge array pulse
  call fire_arrays()

end program


  ! this routine reads in integrals from psi4 instead of gamess
  ! reads in 1- and 2-el integrals from fort.999
  ! and number of elements from data.inp
  subroutine readin_psi4()
    use globalvars
    implicit none

    integer        :: i, j, k, l, m, n
    integer        :: i1, i2, i3, i4
    integer        :: icounter, ios
    real(dp)       :: temp, hh
    character(256) :: iom

    open(newunit=aoints, file='fort.999_'//grdnum, status='old', iostat=ios, iomsg=iom)
    if (ios /= 0) then
      write(*,*) 'Fatal error, could not open fort.999_'//grdnum//'! iomsg:'
      write(*,*) trim(iom)
      stop
    endif

    do i = 1, numel
      read(aoints,*) h(i), t(i), s(i)
    enddo

      do i = 1, numel
        read(aoints,*) x(i), y(i), z(i)
      enddo

    ! if first pass read in the two-electron integrals
    if (write_all) then
      do i = 1, nint2e
        read(aoints,*) rijindex(1,i), rijindex(2,i), &
                       rijindex(3,i), rijindex(4,i), rij(i)
      enddo
    endif

    close(aoints)

    ! the output only prints the lower diagonal, these
    ! loops ensure the full matrix is reconstructed
    k = 1
    do i=1, nrprim
      do j=1, i
        hmat(i,j) = h(k)
        hmat(j,i) = h(k)

        tmat(i,j) = t(k)
        tmat(j,i) = t(k)

        smat(i,j) = s(k)
        smat(j,i) = s(k)

        xint(1,i,j) = x(k)
        xint(1,j,i) = x(k)

        xint(2,i,j) = y(k)
        xint(2,j,i) = y(k)

        xint(3,i,j) = z(k)
        xint(3,j,i) = z(k)

        k = k + 1
      enddo
    enddo

    if (write_all) then
      veeao(:,:,:,:) = 0.0_dp

      ! change to the physicists' notation
      ! integrals are spit out in chemists notation which is x2^* x2 <1/r> x1^* x1
      ! and physicists notation is x2^* x1^* <1/r> x2 x1
      ! so the second index, needs to be switched with the third index
      do i=1, nint2e
        hh = rijindex(2,i)
        rijindex(2,i) = rijindex(3,i)
        rijindex(3,i) = hh
      enddo

      ! reconstruct the symmetry
      ! we've set the full Vee matrix to zero
      ! now we take the indices of the nonzero, symmetry reduced integrals
      ! reconstructing the full veeao matrix
      do i=1, nint2e
        i1 = rijindex(1,i)
        i2 = rijindex(2,i)
        i3 = rijindex(3,i)
        i4 = rijindex(4,i)
        veeao(i1,i2,i3,i4) = rij(i)
        veeao(i1,i4,i3,i2) = rij(i)
        veeao(i3,i2,i1,i4) = rij(i)
        veeao(i3,i4,i1,i2) = rij(i)
        veeao(i2,i1,i4,i3) = rij(i)
        veeao(i2,i3,i4,i1) = rij(i)
        veeao(i4,i1,i2,i3) = rij(i)
        veeao(i4,i3,i2,i1) = rij(i)
      enddo
    endif

    write(*,*) "Integrals from Psi4 read in"
    write(*,*) ' '
  end subroutine

  subroutine symm_ortho()
    use globalvars

    implicit none

    integer               :: i,j,k,mu,nu,la,si
    real(dp), allocatable :: smat2(:,:), one(:,:)
    real(dp), allocatable :: work(:), ws(:)

    character :: jobz, uplo
    integer   :: info, lwork

    if (.not. allocated(smat2))  allocate(smat2(nrprim,nrprim))
    if (.not. allocated(one))    allocate(one(nrprim,nrprim))
    if (.not. allocated(work))   allocate(work(3*nrprim))
    if (.not. allocated(ws))     allocate(ws(nrprim))

    ! calculate s^{-1/2}
    jobz = "V"
    uplo = "U"
    lwork = 3*nrprim

    smat2(:,:) = smat(:,:)

!  open(70,file='smat_beforediag.dat')
!  do i = 1, nrprim
!    do j = 1, nrprim
!      write(70,'(es23.16)') smat(i,j)
!    enddo
!  enddo
!  close(70)

    call dsyev(jobz,uplo,nrprim,smat,nrprim,ws,work,lwork,info)

    if (write_all) write(*,*) "smallest EV smat:", ws(1)

    eigvals(:) = ws(:)
    ldep = 0

    open(24, file="smateigvals.n")
    do i=1, nrprim
      write(24,*) i, ws(i)
      if (dabs(ws(i)) < threshd) then
        ldep = ldep + 1
      endif
      ws(i) = 1.0_dp/dsqrt(ws(i))
    enddo
    close(24)

    open(60, file='smat_xmat.dat')
    do i=1, nrprim
      do j=1, nrprim
        xx(i,j) = smat(i,j)*ws(j)
!        write(60,'(2(es23.16,2x))') xx(i,j), smat(i,j)
        write(60,'(3(es23.16,2x))') xx(i,j), smat(i,j), smat2(i,j)
      enddo
    enddo
    close(60)

  ! test 1 = xx^dagger smat2 xx
    call dgemm('t','n',nrprim,nrprim,nrprim, &
               1.0_dp,xx,nrprim,smat2,nrprim,0.0_dp,smat,nrprim)
    call dgemm('n','n',nrprim,nrprim,nrprim, &
               1.0_dp,smat,nrprim,xx,nrprim,0.0_dp,one,nrprim)

    open(24, file="one")

    do mu=1,nrprim
      do nu=1,nrprim
        if (one(mu,nu) < 1.0d-8) one(mu,nu) = 0.0_dp
        write(24,*) mu, nu, one(mu,nu)
      enddo
    enddo

    close(24)
    if (write_all) then
      write(*,*) ' '
      write(*,*) "Smat orthogonalized"
      write(*, '("Linear dependence threshold: "(es9.2))') threshd
      write(*,'((i4)" basis functions fall below this threshold")') ldep
  !    write(*,'("These will be purged from the set, leaving "(i4)" basis functions in the set")') nrprim - ldep
      write(*,*) "These will be purged from the set"
      write(*,*) ' '
    end if


    deallocate(smat2)
    deallocate(one)
    deallocate(work)
    deallocate(ws)

  end subroutine


  ! most time consuming routine
  ! basis transformation using xx
  subroutine mult()
    use globalvars
    use omp_lib

    implicit none

    integer               :: i,j,k,l,mu,nu,la,si
    real(dp), allocatable :: hel2a(:,:,:,:)

    if (write_all) then
      write(*,*) ' '
      write(*,*) 'Transforming all integrals. Stand by.'
      write(*,*) ' '
    endif

    ! transform 1-electron-integrals
    !$omp parallel &
    !$omp shared (h2, t2, x2, xx, hmat, tmat, xint, nrprim) &
    !$omp private (mu, nu, la, si)
    !$omp do
    do mu=1, nrprim
      do nu=1, nrprim
        h2(mu,nu)   = 0.0_dp
        t2(mu,nu)   = 0.0_dp
        x2(:,mu,nu) = 0.0_dp
!        x2(1,mu,nu) = 0.0_dp
!        x2(2,mu,nu) = 0.0_dp
!        x2(3,mu,nu) = 0.0_dp
        do la=1, nrprim
          do si=1, nrprim
            h2(mu,nu) = h2(mu,nu) + xx(la,mu)*hmat(la,si)*xx(si,nu)
            t2(mu,nu) = t2(mu,nu) + xx(la,mu)*tmat(la,si)*xx(si,nu)
            if (write_all) then
              x2(:,mu,nu) = x2(:,mu,nu) + xx(la,mu)*xint(:,la,si)*xx(si,nu)
!              x2(1,mu,nu) = x2(1,mu,nu) + xx(la,mu)*xint(1,la,si)*xx(si,nu)
!              x2(2,mu,nu) = x2(2,mu,nu) + xx(la,mu)*xint(2,la,si)*xx(si,nu)
!              x2(3,mu,nu) = x2(3,mu,nu) + xx(la,mu)*xint(3,la,si)*xx(si,nu)
            endif
          enddo
        enddo
      enddo
    enddo
    !$omp end do
    !$omp end parallel

    if (write_all) then

!     hel2a(mu,nu,sigma,l) = hel2a(mu,nu,sigma,l) + xx(la,l)*veeao(mu,nu,sigma,la)
!     veeao(mu,nu,k,l)     = veeao(mu,nu,k,l) + xx(sigma,k)*hel2a(mu,nu,sigma,l)
!     hel2a(mu,j,k,l)      = hel2a(mu,j,k,l) + xx(nu,j)*hel2a(mu,nu,k,l)
!     veeao(i,j,k,l)       = veeao(i,j,k,l) + xx(mu,i)*hel2a(mu,j,k,l)
      if (.not. allocated(hel2a))   allocate(hel2a(nrprim,nrprim,nrprim,nrprim))
      ! transform 2-electron-integrals to new basis
      ! LA: I'm unclear what transformation scheme is going on here
      !$omp parallel &
      !$omp shared (xx, nrprim, veeao, hel2a) &
      !$omp private (mu, nu, la, si, i, j, k, l)
      !$omp do
      do mu=1, nrprim
        do nu=1, nrprim
          do si=1, nrprim
            do l=1, nrprim
              hel2a(mu,nu,si,l) = 0.0_dp
              do la=1, nrprim
                hel2a(mu,nu,si,l) = hel2a(mu,nu,si,l) + xx(la,l)*veeao(mu,nu,si,la)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end do

      !$omp do
      do mu=1, nrprim
        do nu=1, nrprim
          do k=1, nrprim
            do l=1, nrprim
              veeao(mu,nu,k,l) = 0.0_dp
              do si=1, nrprim
                veeao(mu,nu,k,l) = veeao(mu,nu,k,l) + xx(si,k)*hel2a(mu,nu,si,l)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end do

      !$omp do
      do mu=1, nrprim
        do j=1, nrprim
          do k=1, nrprim
            do l=1, nrprim
              hel2a(mu,j,k,l) = 0.0_dp
              do nu=1, nrprim
                hel2a(mu,j,k,l) = hel2a(mu,j,k,l) + xx(nu,j)*veeao(mu,nu,k,l)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end do

      !$omp do
      do i=1, nrprim
        do j=1, nrprim
          do k=1, nrprim
            do l=1, nrprim
              veeao(i,j,k,l) = 0.0_dp
               do mu=1, nrprim
                veeao(i,j,k,l) = veeao(i,j,k,l) + xx(mu,i)*hel2a(mu,j,k,l)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end do
      !$omp end parallel

      deallocate(hel2a)
    endif

  end subroutine

  ! this routine is used if no basis transformation to be done
  subroutine nonortho()
    use globalvars

    implicit none

    h2(:,:)   = hmat(:,:)
    t2(:,:)   = tmat(:,:)
    x2(:,:,:) = xint(:,:,:)
    xx(:,:)   = smat(:,:)

  end subroutine

  ! add in routine for hdf5 -> maybe


!  subroutine write_reduced_out(write_all)
  subroutine write_reduced_out()
    use write2str
    use globalvars

    implicit none

!    logical, intent(in) :: write_all
    integer  :: i, j, mu, nu, la, si
    integer  :: nrdbf, nrdbf2, basisinfo
    integer :: k, l, N
    integer :: ntwoe, ao1e, ao2e, aobasis, hmatdat
    integer :: kl, ij, nbf_n
    character(len=10) :: fn1, fn2, fn3, fn4
    character(len=40) :: s_ij
    character(len=40) :: s_kl
!    real(dp), allocatable :: hout(:,:), tout(:,:), sout(:,:), xxout(3,:,:)


    nrdbf = ldep

    nrdbf2 = 0
    nbf_n = nrprim - ldep

    open(newunit=hmatdat, file="data_hmat.asc", status='unknown', position='append')

    ! write out h and t
    do i = nrdbf+1, nrprim
      do j = nrdbf+1, i
        write(hmatdat,outform2) h2(i,j), t2(i,j)
      enddo
    enddo

    close(hmatdat, status='keep')

    ! this only need be written once
    if (write_all) then

      open(newunit=aobasis, file="data_out.asc")

      ! write out x
      do i = nrdbf+1, nrprim
        do j = nrdbf+1, i
          write(aobasis, outform) x2(1,i,j)
        enddo
      enddo
      ! write out y
      do i = nrdbf+1, nrprim
        do j = nrdbf+1, i
          write(aobasis, outform) x2(2,i,j)
        enddo
      enddo
      ! write out z
      do i = nrdbf+1, nrprim
        do j = nrdbf+1, i
          write(aobasis, outform) x2(3,i,j)
        enddo
      enddo

    ! write out s -> though this isn't actually needed in MCEND
!      do i = nrdbf+1, nrprim
!        do j = nrdbf+1, i
!          write(aobasis, outform) xx(i,j)
!        enddo
!      enddo


      ! switch to chemists notation
      do mu = 1, nrprim
        do nu = 1, nrprim
          do la = 1, nrprim
            do si = 1, nrprim
              chmveeao(mu,nu,la,si) = veeao(mu,la,nu,si)
            enddo
          enddo
        enddo
      enddo

      ! routine to write symmetry reduced two-electron integrals
      ! there may be other ways but this works just as well
      ! but only applies for chemists notation as far as I know
!      ntwoe = 0
!      do mu = nrdbf+1, nrprim
!        do nu = nrdbf+1, nrprim
!          if (mu >= nu) then
!            do la = nrdbf+1, nrprim
!              do si = nrdbf+1, nrprim
!                if (la >= si) then
!                  call int2str(fn1, mu)
!                  call int2str(fn2, nu)
!                  s_ij = trim(fn1)//trim(fn2)
!                  call str2int(trim(s_ij), ij)
!
!                  call int2str(fn3, la)
!                  call int2str(fn4, si)
!                  s_kl = trim(fn3)//trim(fn4)
!                  call str2int(trim(s_kl), kl)
!                  if (ij >= kl) then
!                    if (dabs(chmveeao(mu,nu,la,si)) >= 1.0d-15) then
!                      ntwoe = ntwoe + 1
!                      write(aobasis,'(4(i5),(es24.16))') mu-nrdbf, nu-nrdbf, la-nrdbf, si-nrdbf, chmveeao(mu,nu,la,si)
!!                      write(aobasis, '(4(i5), (es24.16))')  mu,nu, la, si, chmveeao(mu,nu,la,si)
!                    endif
!                  endif
!                endif
!              enddo
!            enddo
!          endif
!        enddo
!      enddo
      !> alternate cleaner way
!      ntwoe = 0
!      do i=1, nrprim
!        do j=1, i
!          ij = i*(i + 1)/2 + j
!          do k=1, nrprim
!            do l=1, k
!              kl = k*(k + 1)/2 + l
!              if (ij >= kl) then
!                ntwoe = ntwoe + 1
!              endif
!            enddo
!          enddo
!        enddo
!      enddo
      ntwoe = 0
      do i=nrdbf+1, nrprim
        do j=nrdbf+1, i
          ij = i*(i + 1)/2 + j
          do k=nrdbf+1, nrprim
            do l=nrdbf+1, k
              kl = k*(k + 1)/2 + l
              if (ij >= kl) then
                if (dabs(chmveeao(i,j,k,l)) >= v_cutoff) then
                  ntwoe = ntwoe + 1
                  write(aobasis,'(4(i5),(es24.16))') i-nrdbf, j-nrdbf, k-nrdbf, l-nrdbf, chmveeao(i,j,k,l)
                endif
              endif
            enddo
          enddo
        enddo
      enddo

      close(aobasis)

      open(newunit=basisinfo, file='basis-purge.dat', status='unknown')
      write(basisinfo,*) "Linear dependence threshold: ", threshd
      write(basisinfo,'((i5)" basis functions fall below this threshold")') ldep
      write(basisinfo,'((i5)" basis functions remain ")') nrprim - ldep
      write(basisinfo,'((i10)" Nonzero two-electron integrals")') ntwoe
      close(basisinfo)

    endif

  end subroutine

  subroutine charge_arrays()
    use globalvars
    implicit none

    if (.not. allocated(h))           allocate(h(numel))
    if (.not. allocated(t))           allocate(t(numel))
    if (.not. allocated(s))           allocate(s(numel))
    if (.not. allocated(x))           allocate(x(numel))
    if (.not. allocated(y))           allocate(y(numel))
    if (.not. allocated(z))           allocate(z(numel))
    if (.not. allocated(rij))         allocate(rij(nint2e))
    if (.not. allocated(rijindex))    allocate(rijindex(4,nint2e))
    if (.not. allocated(hmat))        allocate(hmat(nrprim,nrprim))
    if (.not. allocated(tmat))        allocate(tmat(nrprim,nrprim))
    if (.not. allocated(smat))        allocate(smat(nrprim,nrprim))
    if (.not. allocated(eigvals))     allocate(eigvals(nrprim))
    if (.not. allocated(xint))        allocate(xint(3,nrprim,nrprim))
    if (.not. allocated(veeao))       allocate(veeao(nrprim,nrprim,nrprim,nrprim))
    if (.not. allocated(chmveeao))       allocate(chmveeao(nrprim,nrprim,nrprim,nrprim))

    if (.not. allocated(h2))          allocate(h2(nrprim,nrprim))
    if (.not. allocated(t2))          allocate(t2(nrprim,nrprim))
    if (.not. allocated(x2))          allocate(x2(3,nrprim,nrprim))
    if (.not. allocated(xx))          allocate(xx(nrprim,nrprim))

    return

  end subroutine charge_arrays

  subroutine fire_arrays()
    use globalvars
    implicit none

    deallocate(h)
    deallocate(t)
    deallocate(s)
    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(hmat,tmat,smat)
    deallocate(eigvals)
    deallocate(xint)
    deallocate(rij)
    deallocate(rijindex)
    deallocate(h2)
    deallocate(t2)
    deallocate(x2)
    deallocate(xx)
    deallocate(chmveeao)
    deallocate(veeao)

  end subroutine fire_arrays



  
