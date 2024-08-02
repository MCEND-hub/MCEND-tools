module globalconst

  integer, parameter      :: i64 = selected_int_kind(15)
  integer, parameter      :: i32 = selected_int_kind(6)
  integer, parameter      :: dp = selected_real_kind(15,307)
  real(dp), parameter     :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp)                :: au2fs = 0.0241888_dp
  real(dp)                :: w0 = 20*2.6687
  complex(dp)             :: c0 = dcmplx(0.0_dp,0.0_dp)
  complex(dp)             :: c1 = dcmplx(1.0_dp,0.0_dp)
  complex(dp)             :: ci = dcmplx(0.0_dp,1.0_dp)
  integer                 :: fftw_forward = -1
  integer                 :: fftw_backward = 1
  integer(i64)            :: plan1, plan2

end module globalconst

module globalvars

  use globalconst

  implicit none
  complex(dp), allocatable :: sigt1(:), sigt2(:)
  complex(dp), allocatable :: sigt(:)
  complex(dp), allocatable :: sigw1(:), sigw2(:)
  real(dp), allocatable    :: time(:), omega(:), timeold(:)
  real(dp), allocatable    :: allsteps(:)
  real(dp)                 :: omega0     ! ground state energy
  real(dp)                 :: lfreq      ! carrier frequency
  integer                  :: nsteps, tinit, maxj
  character(20)            :: spec
  character                :: pulsedir
  logical                  :: dospec, window
  logical                  :: dohhg

  contains

  subroutine readin()

    integer        :: i, ios
    character(300) :: iom
    character(20)  :: line,slices
    real(dp)       :: x1,x2,x3,x4,x5,x6,x7,x8,x9
    real(dp)       :: x10,x11,x12,realvar,compvar
    integer        :: inputv, sigval, grep_status, system
    integer :: iSize
    logical :: file_exists

    write(*,*) "                       ----------------------- "
    write(*,*) "                       MCEND Spectra Transform "
    write(*,*) "                       ----------------------- "
    write(*,*) " "
    write(*,*) " Options are : "
    write(*,*) " x - z - R - acf - hhgx - hhgz - zR - xR (required) "
    write(*,*) " "
    write(*,*) " followed by w (optional) "
    write(*,*) " followed by number of windows (optional, default is 10) "
    write(*,*) " "
    write(*,*) " "
    write(*,*) " "

    dospec = .false.
    dohhg = .false.
    window = .false.
    pulsedir = 'z'
    lfreq = 0.05695_dp

    call getarg(1, spec)
    !> for now, this takes whatever string you give upon execution
    !> if string is present, FT of aucofu, if not, FT of position
    !> this can be extended to account for polarization direction of pulse
    if (len(trim(spec)) > 0) then
      write(*,*) "Found: ", spec
      if (trim(spec) == 'x') then
        write(*,*) "Spectrum along x"
        pulsedir = 'x'
      else if (trim(spec) == 'z') then
        write(*,*) "Spectrum along z"
        pulsedir = 'z'
      else if (trim(spec) == 'R') then
        write(*,*) "Spectrum along R"
        pulsedir = 'R'
      else if (trim(spec) == 'acf') then
        write(*,*) "Spectrum from autocorrelation function"
        dospec = .true.
      else if (trim(spec) == 'xR') then
        write(*,*) "Performing convolution of electronic and vibrational spectrum - x"
        pulsedir = 'xR'
      else if (trim(spec) == 'zR') then
        write(*,*) "Performing convolution of electronic and vibrational spectrum - z"
        pulsedir = 'zR'
        write(*,*) pulsedir
      else if (trim(spec) == 'hhg') then
        dohhg = .true.
        pulsedir = 'z'
      else if (trim(spec) == 'hhgz') then
        dohhg = .true.
        pulsedir = 'z'
      else if (trim(spec) == 'hhgx') then
        dohhg = .true.
        pulsedir = 'x'
      endif
    else
      write(*,*) 'No parameters given. Reading in z-direction and taking HHG'
    endif

    call getarg(2, spec)
    if (len(trim(spec)) > 0) then
      if (trim(spec) == 'w') then
        write(*,*) "Running sliding DFT"
        window = .true.
        call getarg(3, slices)
        if (len(trim(slices)) > 0) then
          read(slices,*) maxj
        else
          ! default is 10 slices
          write(*,*) 'Using default no. of windows.'
          maxj = 10
        endif
        write(*,*) "Number of windows is ", maxj
      endif
    endif

    if (dospec) then
      write(*,*) " Spectrum generation from autocorrelation function was requested."
    else if (dohhg) then
      write(*,*) " Generating high harmonic spectrum"
      ! for whatever reason system call will not find mcend*out, only *out -> ???
!      call system("grep 'Laser frequency' *.out | awk '{print$4}' > lfreq.tmp")
      grep_status = system("grep -h 'Laser frequency' *.out")
      if (grep_status /= 0) then
        write(*,*) 'Cannot grep value using default value of 0.05695 E/ea0'
      else
        grep_status = system("grep -h 'Laser frequency' *.out | awk '{print$3}' > lfreq.tmp")
        inquire(file='lfreq.tmp', size=iSize, exist=file_exists)
        if (file_exists) then
          if (iSize == 1) then
            write(*,*) 'There is nothing in the lfreq.tmp file, using value of 0.05695 Eh'
          else
            open(24, file='lfreq.tmp', status='unknown', action='read', iostat=ios, iomsg=iom)
            if (ios /= 0) then
              write(*,*) 'IO error:', trim(iom)
              write(*,*) 'Cannot read value for lfreq, using value of 0.05695 Eh'
              lfreq = 0.05695_dp
            else
              read(24,*) lfreq
              write(*,*) 'Found lfreq value of ', lfreq
              close(24)
            endif
          endif
        else
          write(*,*) 'lfreq.tmp file not found, using lfreq value of 0.0595 Eh'
        endif
      endif
    else
      write(*,*) 'Generating standard FFT'
    endif
!    elseif (dospec) write(*,*) " Spectrum generation from autocorrelation function was requested."

!    call system('wc -l < expec.t > count.txt ')
    grep_status = system('wc -l < expec.t > count.txt ')
!    call system('cat expec.t | wc -l > count.txt ')

    open(23, file='count.txt')
    read(23,*) nsteps
    close(23)

    nsteps = nsteps - 1
    allocate(sigt1(nsteps))
    allocate(sigt2(nsteps))
    allocate(sigw1(nsteps))
    allocate(sigw2(nsteps))
    allocate(omega(nsteps))
    allocate(time(nsteps))

!    call dfftw_plan_dft_1d(plan1, nsteps, sigt1, sigw1, fftw_forward, 0)
!    call dfftw_plan_dft_1d(plan2, nsteps, sigt2, sigw2, fftw_forward, 0)
    call dfftw_plan_dft_1d(plan1, nsteps, sigt1, sigw1, fftw_backward, 0)
    call dfftw_plan_dft_1d(plan2, nsteps, sigt2, sigw2, fftw_backward, 0)

    open(newunit=inputv, file='expec.t')
    read(inputv,*) line
!    open(newunit=sigval, file='sigmat.dat')
!    write(sigval,'((a7,1x), 2(a24,1x))') "Step", "Real", "Imag"
    do i=1, nsteps
      read(inputv,*) time(i), x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, realvar, compvar
!      read(inputv,*) time(i), norm, x, y, z, R, x6, x7, x8, x9, x10, Htot, S2, realvar, compvar

      if ((i == 1) .and. (dospec)) omega0 = x11

      if (dospec) then
        sigt1(i) = cmplx(realvar, compvar)
      elseif (trim(pulsedir) == 'x') then
        sigt1(i) = cmplx(x2, 0.0_dp)
      elseif (trim(pulsedir) == 'R') then
        sigt1(i) = cmplx(x5, 0.0_dp)
        if (i > 1) then
          sigt1(i) = sigt1(i) - sigt1(1)
        endif
      elseif (trim(pulsedir) == 'xR') then
        sigt1(i) = cmplx(x2*x5, 0.0_dp)
        if (i > 1) then
          sigt1(i) = sigt1(i) - sigt1(1)
        endif
      elseif (trim(pulsedir) == 'zR') then
        sigt1(i) = cmplx(x4*x5, 0.0_dp)
        if (i > 1) then
          sigt1(i) = sigt1(i) - sigt1(1)
        endif
      else
        sigt1(i) = cmplx(x4, 0.0_dp)
      endif

!       if initial dipole moment is nonzero: we only need induced DM
!      if (i > 1) then
!         sigt1(2:) = sigt1(2:) - sigt1(1)
!         sigt1(i) = sigt1(i) - sigt1(1)
!         write(sigval,'((f7.3,1x),2(es24.16,1x))') time(i), real(sigt1(i)), imag(sigt1(i))
!          endif
    enddo

    close(inputv)

    if (dospec) then
      write(*,*) "Found ground state energy of ", omega0
    else
      omega0 = 0.0_dp
    endif

  end subroutine readin

end module globalvars
