!IU haven't fixed yet that it removes the correct basis functions for sorting
!  The linear dependent ones seem to have the largest Ven elements
!  So they usually show up first - in which case a simple nrdbf value
!  works - but to make it bulletproof it needs to account for the
!  resorting when discarding basis functions and checking the
!  actual resorted eigvals for each AO column/line
!  reads in HF integrals and converts
!  all the matrices in the new orthogonal basis
!IU degenerate vectors may also cause problems
!ifort generate_int_v2.0.f90 -o generate_int_v2.0.o -mkl -O0
!  v2 => now has option to read in psi4 integrals instead of gamess
!  LA -> now features writing only used integrals in most symmetric way possible
!        discounting molecular symmetry
module globalvars

  integer,  parameter    :: dp=selected_real_kind(15,307)
  integer, parameter     :: num_threads=2
  real(dp), parameter    :: v_cutoff = 1.0d-15

  ! 0 for gamess <- currently obsolete
  ! 1 (or anything else really) for psi4
!  integer,  parameter     :: readprogram=2
  ! eigenvalue threshold for linear dependence
!  real(dp), parameter     :: threshd=0.1_dp
  ! prototype moving basis routines are not yet functional
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

!    if (readprogram == 0) then
!      call readin_gms()
!    else
!      call readin_psi4()
!    endif

    ! moving basis part (experimental)
    ! use static or moving basis - orthogonalize accordingly
!    if (orthomethod == 0) then
!      call symm_ortho()
!      call mult()
!    else if (orthomethod == 1) then
!!      call nonortho()
!!      call symm_ortho()
!      call QR_decomp()
!      call coherent_phase()
!      call mult()
!!      call sort_all()
!    endif

!  call write_out()
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

!!!!!!!!!!!!!  prototype moving basis routines !!!!!!!!!!!!!!!

  ! choleski decomposition
  subroutine choleski
  use globalvars

  implicit none

  integer               :: i,j,k,mu,nu,la,si
  real(dp), allocatable :: smat2(:,:), one(:,:)

  character :: uplo
  integer   :: info

  if (.not. allocated(smat2))  allocate(smat2(nrprim,nrprim))
  if (.not. allocated(one))    allocate(one(nrprim,nrprim))
  uplo = "U"
  smat2(:,:) = smat(:,:)
  open(60,file='smat_beforediag.dat')
  do i = 1, nrprim
    do j = 1, nrprim
      write(60,'(es23.16)') smat(i,j)
    enddo
  enddo
  close(60)

  call dpotrf(uplo,nrprim,smat,nrprim,info)

  ! copy upper into lower triangular
  do i = 1, nrprim
    do j = i+1, nrprim
      smat(j,i)=smat(i,j)
    enddo
  enddo

  open(60,file='smat_afterdiag.dat')
  do i = 1, nrprim
    do j = 1, nrprim
      write(60,'(es23.16)') smat(i,j)
    enddo
  enddo
  close(60)

  do i=1,nrprim
    do j=1,nrprim
      xx(i,j) = smat(i,j)
    enddo
  enddo

  ! test 1 = xx^dagger smat2 xx
  call dgemm('n','t',nrprim,nrprim,nrprim, &
             1.0_dp,xx,nrprim,xx,nrprim,0.0_dp,one,nrprim)
!  call dgemm('n','n',nrprim,nrprim,nrprim, &
!             1.0_dp,smat,nrprim,xx,nrprim,0.0_dp,one,nrprim)

  open(24,file="one")
  do mu=1,nrprim
    do nu=1,nrprim
!      if (one(mu,nu) < 1.d-8) one(mu,nu) = 0.0_dp
      write(24,*) mu, nu, one(mu,nu)
    enddo
  enddo
  close(24)


  deallocate(smat2)
  deallocate(one)

  end subroutine

! Computes all eigenvalues and eigenvectors of a real symmetric N × N matrix a. On
! output, elements of a above the diagonal are destroyed. d is a vector of length N that
! returns the eigenvalues of a. xx is an N × N matrix whose columns contain, on output, the
! normalized eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
! Does not sort the eigenvalues
  subroutine jacobi()
  use globalvars

  implicit none
  integer                              :: nrot,i,j,k
  real(dp)                             :: c,hj,sj,sm,tj,tau,theta,tresh
  real(dp), allocatable                :: smat2(:,:), one(:,:)
  real(dp), parameter                  :: abserr=1.0e-09

  if (.not. allocated(smat2))  allocate(smat2(nrprim,nrprim))
  if (.not. allocated(one))    allocate(one(nrprim,nrprim))

  smat2(:,:) = smat(:,:)
  open(60,file='smat_beforediag.dat')
  do i = 1, nrprim
    do j = 1, nrprim
      write(60,*) i,j,smat(i,j)
    enddo
  enddo
  close(60)


  ! initialize xx(i,j)=0, xx(i,i)=1
  xx(:,:)=0.0_dp
  do i=1,nrprim
    xx(i,i)=1.0_dp
  end do

  nrot=0
  ! Sum off-diagonal elements.
  sm=0.0_dp
  do i=1,nrprim-1
    do j=i+1,nrprim
      sm=sm+(smat(i,j))**2
      write(*,*) smat(i,j)
    enddo
  enddo
  ! The normal return after convergence to abserr
  if (sm <= abserr) return

  do while (sm > abserr)
  write(*,*) "off-diagonal elements are ",sm
    ! writing out some information
    open(60,file='smat_afterdiag.dat')
    do i = 1, nrprim
      do j = 1, nrprim
!        xx(j,i)=xx(i,j)
        write(60,*) i,j,xx(i,j)
      enddo
    enddo
    close(60)
    open(24, file="smateigvals.n")
    do i=1,nrprim
      write(24,*) i, smat(i,i)
    enddo
    close(24)
!    test 1 = xx^dagger smat2 xx
!    call dgemm('t','n',nrprim,nrprim,nrprim, &
!               1.0_dp,xx,nrprim,smat2,nrprim,0.0_dp,smat,nrprim)
!    call dgemm('n','n',nrprim,nrprim,nrprim, &
!               1.0_dp,smat,nrprim,xx,nrprim,0.0_dp,one,nrprim)
!    open(24,file="one")
!    do i=1,nrprim
!      do j=1,nrprim
!        if (one(i,j) < 1.d-8) one(i,j) = 0.0_dp
!        write(24,*) i,j,one(i,j)
!      enddo
!    enddo
!    close(24)
   ! average for off-diagonal elements /2
   ! tresh=0.5*sm/float(nrprim*nrprim)
   ! now start the process over
    do i=1,nrprim-1
      do j=i+1,nrprim
      write(*,*) "sm before",sm,nrot,i,j
      write(*,*) smat(:,:)
      if (smat(j,i)**2 <=tresh) cycle  ! do not touch small elements
      sm=sm-2.0_dp*smat(j,i)**2
      write(*,*) "sm after",sm
      write(*,*) smat(:,:)
      write(*,*) "-----------"
      write(*,*) xx(:,:)
      tresh=0.5_dp*sm/float(nrprim*nrprim)
      ! calculate coefficient c and s for Givens matrix
      theta=(smat(j,j)-smat(i,i))/(2.0_dp*smat(j,i))
      tau=0.5_dp*theta/dsqrt(1.0_dp+theta**2)
      write(*,*) "theta,tau",theta,tau
      sj=dsqrt(max(0.5_dp+tau,0.0))
      c=dsqrt(max(0.5_dp-tau,0.0))
      ! recalculate rows i and j
      do k=1,nrprim
        hj=c*smat(i,k)+sj*smat(j,k)
        tj=-sj*smat(i,k)+c*smat(j,k)
        smat(i,k)=hj
        smat(j,k)=tj
      enddo
      ! new matrix a_{k+1} from a_{k}, and eigenvectors
      do k=1,nrprim
        hj=c*smat(k,i)+sj*smat(k,j)
        tj=-sj*smat(k,i)+c*smat(k,j)
        smat(k,i)=hj
        smat(k,j)=tj
        hj=c*xx(k,i)+sj*xx(k,j)
        tj=-sj*xx(k,i)+c*xx(k,j)
        xx(k,i)=hj
        xx(k,j)=tj
      enddo
      enddo
    enddo
    nrot=nrot+1
    write(*,*) nrot," rotations in Jacobi."
  enddo

  write(*,*) nrot," rotations in Jacobi."
  write(*,*) "off-diagonal elements are ",sm
  write(*,*) "This is more than the maximum number! Aborting run."

  deallocate(smat2)
  deallocate(one)

  end subroutine jacobi

  ! SVD to generate transformation matrix to orthogonal basis
  ! most gentle interference with the basis, but
  ! does not provide degree of linear dependency
  subroutine svdecomp()
  use globalvars

  implicit none

  integer               :: i,j,k,mu,nu,la,si
  real(dp), allocatable :: spoh(:,:), smoh(:,:), one(:,:)
  real(dp), allocatable :: smat2(:,:)
  real(dp), allocatable :: work(:), ws(:)
  integer               :: info, lwork

  if (.not. allocated(one))   allocate(one(nrprim,nrprim))
  if (.not. allocated(spoh))  allocate(spoh(nrprim,nrprim))
  if (.not. allocated(smoh))  allocate(smoh(nrprim,nrprim))
  if (.not. allocated(smat2)) allocate(smat2(nrprim,nrprim))
  if (.not. allocated(ws))    allocate(ws(nrprim))
  lwork=5*nrprim
  if (.not. allocated(work)) allocate(work(lwork))
  xx(:,:) = 0.0_dp
  smat2(:,:)=smat(:,:)

  call dgesvd('A','A',nrprim,nrprim,smat,nrprim,ws,spoh,nrprim, &
              smoh,nrprim,work,lwork,info)
  ! remember the eigenvalues after bombing out from the routine
  eigvals(:)=ws(:)

  ! a simple threshold does not work for SVD
  ! write(*,'(es23.16)') "smallest EV smat:", ws(nrprim)
  ldep = 0
  open(24, file="smateigvals.n")
  do i=1,nrprim
    write(24,*) i, ws(i)
    if (dabs(ws(i)) < threshd) then
      ldep = ldep + 1
    endif
    ws(i) = 1.0_dp/dsqrt(ws(i))
  enddo
  close(24)

  ! dgesv returns U (spoh) and V^t (smoh)
  ! transformation matrix should be X = U V^t
  ! so that xx = spoh smoh
  ! but here only (smoh)^t does the job
  ! so that X = U (V^t)^t = U V
  ! or only spoh since
  ! here spoh = smoh^t (U = (V^t)^t)
  ! (because smat is symmetric: U smat V = U smat U^t)
  call dgemm('n','t',nrprim,nrprim,nrprim, &
           1.0_dp,spoh,nrprim,smoh,nrprim,0.0_dp,xx,nrprim)
!  do mu=1,nrprim
!    do nu=1,nrprim
!      xx(mu,nu)=spoh(mu,nu)*ws(nu)
!    enddo
!  enddo

  ! test X X^t = 1 ?
 call dgemm('n','t',nrprim,nrprim,nrprim, &
            1.0_dp,xx,nrprim,xx,nrprim,0.0_dp,one,nrprim)
 open(24,file="one")
 do mu=1,nrprim
   do nu=1,nrprim
!    one(mu,nu) = 0.0_dp
!    do la=1,nrprim
!      one(mu,nu) = one(mu,nu) + xx(la,mu)*xx(la,nu)
!    enddo
     if (one(mu,nu) < 1.d-8) one(mu,nu) = 0.0_dp
     write(24,*) mu, nu, one(mu,nu)
   enddo
 enddo
 close(24)

  deallocate(one)
  deallocate(smat2)
  deallocate(spoh)
  deallocate(smoh)
  deallocate(ws)
  deallocate(work)

  end subroutine

  subroutine QR_decomp
  use globalvars

  implicit none

  integer               :: i,j,k,mu,nu,la,si
  real(dp), allocatable :: smat2(:,:), one(:,:)
  real(dp), allocatable :: work(:), ws(:)

  integer   :: info, lwork

  lwork = nrprim
  if (.not. allocated(smat2))  allocate(smat2(nrprim,nrprim))
  if (.not. allocated(one))    allocate(one(nrprim,nrprim))
  if (.not. allocated(work))   allocate(work(lwork))
  if (.not. allocated(ws))     allocate(ws(nrprim))

  smat2(:,:) = smat(:,:)
  open(60,file='smat_beforediag.dat')
  do i = 1, nrprim
    do j = 1, nrprim
      write(60,*) i,j,smat(i,j)
    enddo
  enddo
  close(60)
  one(:,:)=0.0_dp
! do i = 1, nrprim
!   one(i,i)=1.0_dp
! enddo

  call dgeqrf(nrprim,nrprim,smat,nrprim,ws,work,lwork,info)
! matrix q is represented in reflectors as saved in ws
! together with lower triangular of smat
! matrix r is stored in upper triangular of smat
! rebuilt r
  do i = 1, nrprim
    do j = i, nrprim
      one(i,j)=smat(i,j)
      one(j,i)=one(i,j)
    enddo
  enddo
! matrix q is rebuilt through:
  call dorgqr(nrprim,nrprim,nrprim,smat,nrprim,ws,work,lwork,info)
! check A=QR
!=======

! test x s x^t=1 ?
! xx*smat2*xx
! Decide from here which basis functions get discarded -
! which ones still overlap and to what extend?
!!! is this s x x^t = 1 then?
  call dgemm('n','t',nrprim,nrprim,nrprim, &
             1.0_dp,smat2,nrprim,xx,nrprim,0.0_dp,smat,nrprim)
!>>>>>>> 58d719155b3457949802077ed5e9d65f5fc59a82
  call dgemm('n','n',nrprim,nrprim,nrprim, &
             1.0_dp,smat,nrprim,one,nrprim,0.0_dp,xx,nrprim)
  do i = 1, nrprim
    do j = 1, nrprim
      write(*,*) i,j,smat2(i,j),xx(i,j)
    enddo
  enddo
! check q^T q = 1
  call dgemm('t','n',nrprim,nrprim,nrprim, &
             1.0_dp,smat,nrprim,smat,nrprim,0.0_dp,one,nrprim)
! multiply with real matrix
! call dormqr('L','T',nrprim,nrprim,nrprim,smat,nrprim,ws,one,nrprim,WORK,LWORK,INFO)
! call dtrtrs('U','N','U',nrprim,nrprim,smat,nrprim,one,nrprim,INFO)

  open(60,file='smat_afterdiag.dat')
  do i = 1, nrprim
    do j = 1, nrprim
      write(60,*) i,j,smat(i,j),one(i,j)
    enddo
  enddo
  close(60)

  xx(:,:) = smat(:,:)

! test 1 = xx^dagger smat2 xx
!!call dgemm('n','n',nrprim,nrprim,nrprim, &
!!           1.0_dp,smat2,nrprim,xx,nrprim,0.0_dp,smat,nrprim)
! call dgemm('t','n',nrprim,nrprim,nrprim, &
!            1.0_dp,xx,nrprim,smat2,nrprim,0.0_dp,smat,nrprim)
! call dgemm('n','n',nrprim,nrprim,nrprim, &
!            1.0_dp,smat,nrprim,xx,nrprim,0.0_dp,one,nrprim)
  call dgemm('n','n',nrprim,nrprim,nrprim, &
             1.0_dp,smat2,nrprim,xx,nrprim,0.0_dp,smat,nrprim)
  call dgemm('t','n',nrprim,nrprim,nrprim, &
             1.0_dp,xx,nrprim,smat,nrprim,0.0_dp,one,nrprim)

  do mu=1,nrprim
    ws(mu)=one(mu,mu)
    do nu=1,nrprim
      one(mu,nu) = 0.0_dp
      do la=1,nrprim
        do si=1,nrprim
          one(mu,nu) = one(mu,nu) + xx(la,mu)*smat2(la,si)*xx(si,nu)
        enddo
      enddo
    enddo
  enddo

  eigvals(:) = ws(:)
  ldep = 0
  open(24, file="smateigvals.n")
  do i=1,nrprim
    write(24,*) i, ws(i)
    if (dabs(ws(i)) < threshd) then
      ldep = ldep + 1
    endif
    ws(i) = 1.0_dp/dsqrt(ws(i))
  enddo
  close(24)
!  do i=1,nrprim
!    do j=1,nrprim
!      xx(i,j) = smat(i,j)*ws(j)
!    enddo
!  enddo
!
!  call dgemm('n','n',nrprim,nrprim,nrprim, &
!             1.0_dp,smat2,nrprim,xx,nrprim,0.0_dp,smat,nrprim)
!  call dgemm('t','n',nrprim,nrprim,nrprim, &
!             1.0_dp,xx,nrprim,smat,nrprim,0.0_dp,one,nrprim)
!  do mu=1,nrprim
!    do nu=1,nrprim
!      if (one(mu,nu) < 1.d-8) one(mu,nu) = 0.0_dp
!      write(*,*) mu, nu, one(mu,nu)
!    enddo
!  enddo

  deallocate(smat2)
  deallocate(one)
  deallocate(work)
  deallocate(ws)

  end subroutine

  ! ensure that the sign of eigenvectors is consistent
  subroutine coherent_phase
  use globalvars

  implicit none
  integer            :: i,large,idamax

  ! set the phase of each column of a matrix so the largest
  ! element is positive
  do i = 1,nrprim
    large = idamax(nrprim,xx(1,i),1)
    ! just making sure we're still in the array
    if (large .le. 0) large=1
    if (large .gt. nrprim) large=1
    ! now if the largest element is negative, change the
    ! sign of the vector
    if (xx(large,i) .lt. 0.0_dp) then
      xx(:,i) = -1.0_dp*xx(:,i)
    endif
  enddo

  end subroutine

  ! this sorts the AO basis according to magnitude
  ! of Ven elements and rearranges all other matrices
  ! accordingly
  ! for the moving basis to keep elements consistent
  ! between nuclear grid points
  subroutine sort_all()
  use globalvars

  implicit none

  integer          :: i,j,k,l
  integer          :: i1,i2,i3,i4
  real(dp)         :: a
  real(dp), allocatable :: vec1(:),aux1(:,:),aux2(:,:),aux3(:,:)
  real(dp), allocatable :: aux4(:,:,:),aux5(:,:,:,:)
  integer,  allocatable  :: track(:)

! vec1 contains Ven diagonal elements and track keeps
! count of the reordering
  allocate(vec1(nrprim))
  allocate(track(nrprim))
! aux holds temporary memory
  allocate(aux1(nrprim,nrprim))
  allocate(aux2(nrprim,nrprim))
  allocate(aux3(nrprim,nrprim))
  allocate(aux4(3,nrprim,nrprim))
  allocate(aux5(nrprim,nrprim,nrprim,nrprim))

! arrays that have to be resorted
! h2(mu,nu),t2(mu,nu),x2(1-3,mu,nu),smat(mu,nu),veeao(i,j,k,l)

! look at the diagonal elements of Ven only
! and sort according to their magnitude
  do i=1,nrprim
    vec1(i) = h2(i,i) - t2(i,i)
    track(i) = i
!   write(*,*) "Diag Ven :",i, vec1(i),track(i),eigvals(i)
  enddo

! Now, look at the reordered eigenvalues to determine which AO's
! need to be removed
! Write the reordered values to allow data checks
! for this to work, we need an array that stores
! the AO number

! ldep = 0
! open(24, file="smateigvals.n", status="old", position="append", action="write")
! do i=1,nrprim
!   write(24,*) i, eigvals(i)
!   if (dabs(eigvals(i)) < threshd) then
!     write(*,*) "Found linear dependency for basis function", i
!     write(*,*) "Smat eigenvalue is ",eigvals(i)
!     ldep = ldep + 1
!   endif
! enddo
! close(24)

! Now resorting
! Using the simplest "straight insertion" routine
! which scales as (nrprim)^2
! should be ok for small nrprim (<20 definitely, <100 probably)
! for larger basis sets implement shell sort or heapsort
! quicksort probably a bit overambitious for the matrix sizes
! here and can run into memory problems
  do j=2, nrprim
    a = vec1(j)
    i = j - 1
    do while (i >= 1)
      if (vec1(i) <= a) then
        vec1(i+1) = a
        track(i+1) = j
!       eigvals(i+1) = eigvals(j)
        exit
      else
        vec1(i+1) = vec1(i)
        track(i+1) = track(i)
!       eigvals(i+1) = eigvals(i)
        vec1(i) = a
        track(i) = j
!       eigvals(i) = eigvals(j)
        i = i - 1
      endif
    enddo
  enddo

! open(50, file='diag_ven_order.dat', status='unknown', position='append')
! do i=1,nrprim
!    write(*,*) "Diag Ven :", i, vec1(i), track(i), eigvals(i)
!   write(50,'((a10,1x) (i5,1x), (es23.15), (i5,1x), (es23.15))') "Diag Ven:", i, vec1(i), track(i), eigvals(i)
! enddo
! close(50)

! Now, resort: switch rows in the arrays according to track
! only lower triangular to make use of symmetry
! at least for the one-electron part
! this can certainly be done using more efficient techniques
! but for now I am only testing if this works at all
! especially for 2el-integrals!!!
  aux1(:,:) = h2(:,:)
  aux2(:,:) = t2(:,:)
  aux3(:,:) = smat(:,:)
  aux4(:,:,:) = x2(:,:,:)
  aux5(:,:,:,:) = veeao(:,:,:,:)
  veeao(:,:,:,:) = 0.0_dp
  do i=1,nrprim
    if (track(i) /= i) then
      h2(:,i) = aux1(:,track(i))
      h2(i,i) = aux1(track(i),track(i))
      t2(:,i) = aux2(:,track(i))
      t2(i,i) = aux2(track(i),track(i))
      smat(:,i) = aux3(:,track(i))
      smat(i,i) = aux3(track(i),track(i))
      x2(:,:,i) = aux4(:,:,track(i))
      x2(:,i,i) = aux4(:,track(i),track(i))
      veeao(i,i,i,i) = aux5(track(i),track(i),track(i),track(i))
! rearrange previous indices

      do j=i, 1, -1
        h2(j,i) = aux1(track(j),track(i))
        t2(j,i) = aux2(track(j),track(i))
        smat(j,i) = aux3(track(j),track(i))
        x2(:,j,i) = aux4(:,track(j),track(i))
        do k=j, 1,-1
          do l=k, 1,-1
            veeao(k,l,j,i) = aux5(track(k),track(l),track(j),track(i))
            veeao(k,j,l,i) = aux5(track(k),track(j),track(l),track(i))
            veeao(l,j,k,i) = aux5(track(l),track(j),track(k),track(i))
            veeao(j,l,k,i) = aux5(track(j),track(l),track(k),track(i))

!            veeao(l,k,j,i) = aux5(track(l),track(k),track(j),track(i))  ! symm
!            veeao(j,k,l,i) = aux5(track(j),track(k),track(l),track(i))  ! symm
!            veeao(l,i,j,k) = aux5(track(l),track(i),track(j),track(k))  ! symm
!            veeao(j,i,l,k) = aux5(track(j),track(i),track(l),track(k))  ! symm
!
!            veeao(i,l,k,j) = aux5(track(i),track(l),track(k),track(j))  ! symm
!            veeao(i,j,k,l) = aux5(track(i),track(j),track(k),track(l))  ! symm
!            veeao(k,l,i,j) = aux5(track(k),track(l),track(i),track(j))  ! symm
!            veeao(k,j,i,l) = aux5(track(k),track(j),track(i),track(l))  ! symm
!
!           veeao(l,k,j,i) = aux5(track(l),track(k),track(j),track(i))  ! symm
            veeao(j,k,l,i) = veeao(l,k,j,i)
            veeao(l,i,j,k) = veeao(l,k,j,i)
            veeao(j,i,l,k) = veeao(l,k,j,i)

            veeao(i,l,k,j) = veeao(l,k,j,i)
            veeao(i,j,k,l) = veeao(l,k,j,i)
            veeao(k,l,i,j) = veeao(l,k,j,i)
            veeao(k,j,i,l) = veeao(l,k,j,i)

            veeao(l,j,i,k) = aux5(track(l),track(j),track(i),track(k))
            veeao(l,k,i,j) = aux5(track(l),track(k),track(i),track(j))
            veeao(j,k,i,l) = aux5(track(j),track(k),track(i),track(l))
            veeao(j,l,i,k) = aux5(track(j),track(l),track(i),track(k))

            veeao(k,i,l,j) = aux5(track(k),track(i),track(l),track(j))
            veeao(k,i,j,l) = aux5(track(k),track(i),track(j),track(l))
            veeao(l,i,k,j) = aux5(track(l),track(i),track(k),track(j))
            veeao(j,i,k,l) = aux5(track(j),track(i),track(k),track(l))


            veeao(i,k,l,j) = aux5(track(i),track(k),track(l),track(j))
            veeao(i,k,j,l) = aux5(track(i),track(k),track(j),track(l))
            veeao(i,l,j,k) = aux5(track(i),track(l),track(j),track(k))
            veeao(i,j,l,k) = aux5(track(i),track(j),track(l),track(k))
          enddo
        enddo
      enddo
    endif
  enddo

! copy lower triangular
  do i=1,nrprim
    do j=i+1,nrprim
      h2(j,i)=h2(i,j)
      t2(j,i)=t2(i,j)
      smat(j,i)=smat(i,j)
      x2(:,j,i)=x2(:,i,j)
    enddo
  enddo

! check arrays
! do i=1,nrprim
!   write(*,*) (h2(j,i),j=1,nrprim)
! enddo
! write(*,*) "-------"
! do i=1,nrprim
!   write(*,*) (aux1(j,i),j=1,nrprim)
! enddo
! write(*,*) "-------"
! do i=1,nrprim
!   write(*,*) track(i),i,vec1(i)
! enddo
! do i=1,nrprim
!   write(*,*) (x2(3,j,i),j=1,nrprim)
! enddo
! write(*,*) "-------"
! do i=1,nrprim
!   write(*,*) (aux4(3,j,i),j=1,nrprim)
! enddo
! do i=1,nrprim
! do k=1,nrprim
! do l=1,nrprim
!   write(*,*) (veeao(l,k,j,i),j=1,nrprim)
! enddo
! write(*,*) "-------"
! enddo
! enddo
! do i=1,nrprim
!   write(*,*) (aux5(5,5,j,i),j=1,nrprim)
! enddo
! do i1=1,nrprim
!   do i2=1,nrprim
!     do i3=1,nrprim
!       do i4=1,nrprim
!   write(*,*) "--------------------"
!   write(*,*) i1,i2,i3,i4,veeao(i1,i2,i3,i4)
!   write(*,*) i1,i4,i3,i2,veeao(i1,i4,i3,i2)
!   write(*,*) i3,i2,i1,i4,veeao(i3,i2,i1,i4)
!   write(*,*) i3,i4,i1,i2,veeao(i3,i4,i1,i2)
!   write(*,*) i2,i1,i4,i3,veeao(i2,i1,i4,i3)
!   write(*,*) i2,i3,i4,i1,veeao(i2,i3,i4,i1)
!   write(*,*) i4,i1,i2,i3,veeao(i4,i1,i2,i3)
!   write(*,*) i4,i3,i2,i1,veeao(i4,i3,i2,i1)
!   write(*,*) "--------------------"
!       enddo
!     enddo
!   enddo
! enddo


  deallocate(vec1)
  deallocate(track)
  deallocate(aux1,aux2,aux3,aux4,aux5)

end subroutine
!!!!! OBSOLETE SUBROUTINES !!!!!!

  ! reads in 1- and 2-el integrals from fort.999
  ! and number of elements from data.inp
  !> OBSOLETE <!
  subroutine readin_gms()
  !> obsolete in latest version of MCEND
  use globalvars

  implicit none

  integer          :: i,j,k,l,m,n, ios
  integer          :: i1,i2,i3,i4,hh
  character(256)   :: iom
  real(dp)         :: temp


  open(newunit=aoints, file='fort.999', status='old', iostat=ios, iomsg=iom)
  if (ios /= 0) then
    write(*,*) 'Fatal error, could not open fort.999, iomsg='//trim(iom)
    stop
  endif

  ! read in the h, s, and t values from the gamess output
  do i=1,numel
    read(aoints,*) h(i), s(i), t(i)
  enddo
  ! read in the two-electron indices and the values
  do i=1,nint2e
    read(aoints,*) rijindex(1,i),rijindex(2,i), &
         rijindex(3,i),rijindex(4,i),temp,rij(i)
  enddo

  ! account for the screwy gamess formatting of the dipole integrals
  m = int(nrprim/5)

  if (5*m < nrprim) m = m + 1
  if (m == 0) m = 1

  do n=1, m
    do i = (n-1)*5+1, nrprim
      l = 5*n
      if (l > i) l = i
      do j = (n-1)*5+1,l
        k = i*(i - 1)/2 + j
        read(aoints,*) x(k)
        xint(1,i,j) = x(k)
        xint(1,j,i) = x(k)
      enddo
    enddo
  enddo
  do n=1,m
    do i=(n-1)*5+1, nrprim
      l = 5*n
      if (l > i) l = i
      do j=(n-1)*5+1, l
        k = i*(i - 1)/2 + j
        read(aoints,*) y(k)
        xint(2,i,j) = y(k)
        xint(2,j,i) = y(k)
      enddo
    enddo
  enddo
  do n=1, m
    do i=(n-1)*5+1, nrprim
      l = 5*n
      if (l > i) l=i
      do j=(n-1)*5+1,l
        k = i*(i - 1)/2 + j
        read(aoints,*) z(k)
        xint(3,i,j) = z(k)
        xint(3,j,i) = z(k)
      enddo
    enddo
  enddo
  write(*,*) "integrals from gamess read in"

  close(aoints)

  k = 1
  do i=1, nrprim
    do j=1,i
      hmat(i,j) = h(k)
      hmat(j,i) = h(k)

      tmat(i,j) = t(k)
      tmat(j,i) = t(k)

      smat(i,j) = s(k)
      smat(j,i) = s(k)
      k = k + 1
    enddo
  enddo


  ! unpack 2-electron-integrals and
  ! change to "physicists notation"
    veeao(:,:,:,:) = 0.0_dp

    do i=1,nint2e
      hh = rijindex(2,i)
      rijindex(2,i) = rijindex(3,i)
      rijindex(3,i) = hh
    enddo

    do i=1, nint2e
      i1 = rijindex(1,i)
      i2 = rijindex(2,i)
      i3 = rijindex(3,i)
      i4 = rijindex(4,i)
      ! this should take care of symmetry
      veeao(i1,i2,i3,i4)=rij(i)
      veeao(i1,i4,i3,i2)=rij(i)
      veeao(i3,i2,i1,i4)=rij(i)
      veeao(i3,i4,i1,i2)=rij(i)
      veeao(i2,i1,i4,i3)=rij(i)
      veeao(i2,i3,i4,i1)=rij(i)
      veeao(i4,i1,i2,i3)=rij(i)
      veeao(i4,i3,i2,i1)=rij(i)
    enddo

  end subroutine

  subroutine write_out()
  use globalvars

  implicit none

  integer  :: i, j, mu, nu, la, si
  integer  :: nrdbf, nrdbf2, basisinfo

  ! zero for all basis functions
  ! discard one basis function for moving
  nrdbf = ldep
!  nrdbf = 0
  write(*,*) "Linear dependence threshold: ", threshd
  write(*,'((i5)" basis functions fall below this threshold")') ldep
  write(*,'((i5)" basis functions remain ")') nrprim - ldep
  write(*,*) "Purging undesirable basis functions. Stand by."

  open(newunit=basisinfo, file='basis-purge.dat', status='unknown')
  write(basisinfo,*) "Linear dependence threshold: ", threshd
  write(basisinfo,'((i5)" basis functions fall below this threshold")') ldep
  write(basisinfo,'((i5)" basis functions remain ")') nrprim - ldep
  close(basisinfo)

! For resorted arrays: Need to remove correct basis function
! i.e. if initially first basis function had small eigenvalue
! of the overlap matrix, make sure that one is still the same
! after resorting

! nrdbf=2
  nrdbf2=0
! no of discarded basis functions
! LiH,STO-3G
!      nrdbf=2
! LiH, 6-311G(d,p)
!      nrdbf=1

  open(26,file="data_out.asc")

!  write(*,*) "now writing out 1el integrals"
  !      write out h

  do i=nrdbf+1,nrprim-nrdbf2
    do j=nrdbf+1,nrprim-nrdbf2
!          if ((i /= 2) .and. (j /= 2))
      write(26,outform) h2(i,j)
    enddo
  enddo
  !      write out t
  do i=nrdbf+1,nrprim-nrdbf2
    do j=nrdbf+1,nrprim-nrdbf2
!          if ((i /= 2) .and. (j /= 2))
      write(26,outform) t2(i,j)
    enddo
  enddo
  !      write out x
  do i=nrdbf+1,nrprim-nrdbf2
    do j=nrdbf+1,nrprim-nrdbf2
!          if ((i/ = 2) .and. (j/ = 2))
      write(26,outform) x2(1,i,j)
    enddo
  enddo
  !      write out y
  do i=nrdbf+1,nrprim-nrdbf2
    do j=nrdbf+1,nrprim-nrdbf2
!          if ((i /= 2) .and. (j /= 2))
      write(26,outform) x2(2,i,j)
    enddo
  enddo
  !      write out z
  do i=nrdbf+1,nrprim-nrdbf2
    do j=nrdbf+1,nrprim-nrdbf2
!          if ((i /= 2) .and. (j /= 2))
      write(26,outform) x2(3,i,j)
    enddo
  enddo
  ! write out s
  do i=nrdbf+1,nrprim-nrdbf2
    do j=nrdbf+1,nrprim-nrdbf2
!          if ((i /= 2) .and. (j /= 2))
      write(26,outform) xx(i,j)
    enddo
  enddo

  open(29,file='smat-test.dat')
  do i=nrdbf+1,nrprim-nrdbf2
    do j=nrdbf+1,nrprim-nrdbf2
      write(29,outform) xx(i,j)
    enddo
  enddo
  close(29)

   do mu=nrdbf+1,nrprim-nrdbf2
     do nu=nrdbf+1,nrprim-nrdbf2
       do la=nrdbf+1,nrprim-nrdbf2
         do si=nrdbf+1,nrprim-nrdbf2
!           if ((mu /= 2) .and. (nu /= 2)) then
!             if ((la /= 2) .and. (si /= 2))
           write(26,outform) veeao(mu,nu,la,si)
!          write(99,*) mu,nu,la,si,veeao(mu,nu,la,si)
!            endif
         enddo
       enddo
     enddo
   enddo

  close(26)

  end subroutine
