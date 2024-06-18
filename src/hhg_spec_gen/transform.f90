module data_prep

  use globalvars
  use globalconst

  implicit none

  contains
    subroutine sliding_window(j)

    integer           :: i, j, jbin, j1, j2, j3
    real(dp)          :: sigma

    ! apply gaussian window function centered around t
    ! find out size of slices
    jbin = nsteps/maxj
    j1 = (j - 1)*jbin
    j2 = j*jbin
    j3 = (j2 - j1)/2 + j1
!   sigma = 2*jbin/10.0_dp
    sigma = 2*jbin/2.0_dp
!   sigma = 0.8_dp
    sigt2(:) = 0.0_dp

!   write(*,*) "Doing sliding window..."
!   write(*,*) "Window ",j," center ", j3, j1, j2, jbin
!   write(*,*) "sigma ",sigma
    do i=1, nsteps
      sigt2(i) = sigt1(i)*dexp(-1.0_dp*((i - j3)**2/(2.0_dp*sigma**2)))
!     write(99,*) i, dexp(-1.0_dp*((i - j3)**2/(2.0_dp*sigma**2))),real(sigt2(i))
    enddo
!   do i=1, nsteps
!     write(99,*) i, real(sigt2(i)), aimag(sigt2(i))
!   enddo

    return

    end subroutine sliding_window

    subroutine prepare_data()

    integer           :: i
    real(dp)          :: sigma

    ! first apply sine, then gaussian window function
    sigma = 0.8_dp

    do i=1, nsteps
!      sigt2(i) = sigt1(i)*dsin(i*pi/(nsteps - 1.0_dp))
      sigt2(i) = sigt1(i)*dsin(i*pi/(nsteps - 1.0_dp)) &
                 *dexp(-0.5_dp*((i - nsteps - 1.0_dp)/(sigma*(nsteps - 1.0_dp))))**2
               ! don't these 2.0_dp's cancel out? -> yes they do, checked in python
!                 *dexp(-0.5_dp*(((i - nsteps - 1.0_dp)/2.0_dp)/(sigma*(nsteps - 1.0_dp)/2.0_dp)))**2
!      sigt2(i) = sigt1(i)
    enddo

    return

    end subroutine prepare_data

    subroutine transform()

      integer               :: i, j, nleast, nmost
      real(dp)              :: stepw, highestw, maxt2
      real(dp)              :: norm, aux, dmu, dt, maxt
      real(dp)              :: maxtfs, dw
      complex(dp)           :: prefac

      dt = (time(2) - time(1))/au2fs
      write(*,'("time step was ",(f9.4)," fs or ",(f9.4)," atomic units")') dt*au2fs, dt
      maxtfs = time(nsteps) - time(1)
      maxt = maxtfs/au2fs
      write(*,'("simulation ran for ",(f7.2)," fs or ",(f11.1)," atomic units")') maxtfs, maxt

      ! resolution of DFT is given by 1/(Delta t) the smaller the time step, the better
      ! the resolution because the high frequencies get resolved with a small time step
      ! a large time step will not be able to resolve fast modulations in the signal
      stepw = 2.0_dp*pi/dt

      ! the step size delta omega in frequency domain is given by the maximum propagation time
      ! the longer the propagation time, the better resolved is the
      ! signal in the frequency domain because there are more "samples"
      dw = 2.0_dp*pi/maxt

      write(*,*) "Sampling step size in frequency domain ", dw, " a.u."
      write(*,*) "maximum resolved frequency ", stepw, " a.u."
!      write(*,*) "for a laser frequency of 0.005695 Eh (800nm), that corresponds to ", &
!      if (.not. dospec) then
      if (dohhg) then
        write(*,*) "for the laser frequency of ", lfreq, " Eh, that corresponds to ", stepw/lfreq, " high harmonics"
      endif

      call dfftw_execute_dft(plan1,sigt1,sigw1)
      call dfftw_execute_dft(plan2,sigt2,sigw2)

!   do i=1, nsteps
!     write(99,*) i, real(sigt2(i)), cdabs(sigw2(i))
!   enddo
!      sigw1(:) = c0
!      sigw2(:) = c0

      norm = (1.0_dp/(time(nsteps) - time(1))*au2fs)
      write(*,*) "norm", norm, time(nsteps), time(nsteps)/au2fs

      !> the sign is related to the definition of the aucofu
      !> < Psi(t) | Psi(0) > or < Psi(0) | Psi(t) >
      !> in our case, we need a positive exponent for absorption spectra
      !> how does the fftw acount for this?, Do we just add a negative sign?

      if (dospec) then
        prefac = ci
      else
        prefac = -ci
      endif

      allocate(allsteps(nsteps))

      do i=1, nsteps
        allsteps(i) = dble(i)
      enddo

      if (dospec .and. omega0 < 0.0_dp) then
!        omega(:) = omega0 + (allsteps(:) - 1)*dw
        omega(:) = (allsteps(:) - 1)*dw - 1
      else
        omega(:) = (allsteps(:) - 1)*dw
      endif

      if (dohhg) then
        sigw1(:) = (sigw1(:)*norm)**2
        sigw2(:) = (sigw2(:)*norm)**2
      else
        sigw1(:) = (sigw1(:))/(2*pi)
        sigw2(:) = (sigw2(:))/(2*pi)
      endif

!      do i=1, nsteps
!        if (dospec .and. omega0 < 0.0_dp) then
!          omega(i) = omega0 + (i - 1)*dw
!!        omega(i) =  -1 + (i - 1)*dw
!        else
!          omega(i) = (i - 1)*dw
!        endif

        ! grid for the time
!        do j=1, nsteps
!          ! convert time from fs to a.u.
!          aux = time(j)/0.0241888_dp
!          ! either use delta omega / delta t ...
!          ! for spectra, this gives peaks at relative energies with ground state at zero
!          sigw1(i) = sigw1(i) + cdexp(prefac*omega(i)*aux)*sigt1(j)
!          ! ... or use grid numbers i and j
!          ! for acf spectra, this gives absolute energies
!!          sigw2(i) = sigw2(i) + cdexp(prefac*2.0_dp*pi*(i-1)*(j-1)/dble(nsteps))*sigt2(j)
!!          sigw2(i) = sigw2(i) + cdexp(-ci*2.0_dp*pi*(i-1)*(j-1)/dble(nsteps))*sigt2(j)
!          ! For HHG sigw2 and sigw1 are pretty much equivalent
!        enddo
!        if (.not. dospec) then
!          sigw1(i) = (sigw1(i)*norm)**2
!          sigw2(i) = (sigw2(i)*norm)**2
!        else
!          sigw1(i) = (sigw1(i))/(2*pi)
!          sigw2(i) = (sigw2(i))/(2*pi)
!         endif
!      enddo

  end subroutine transform

  subroutine write_out(j)

    integer   :: i
    integer, optional   :: j
    integer   :: outfile
    character(25) :: formout,formout2
    character(25) :: formtxt
    character(25) :: filename
    character(4) :: str
    formout = '(5(e20.12,2x))'
    formout2 = '(i5,5(e20.12,2x))'
    formtxt = '(5(a20,2x))'

    if (present(j)) then
      write(str,'(i4.4)') j
      filename='ft_omega_'//str//'.dat'
      write(*,*) "Writing to file ", filename
    else
      filename='ft_omega.dat'
      write(*,*) "Writing to file ", filename
    endif
    open(newunit=outfile, file=filename)
    ! original grid was shifted by +8.D0 Hartree to remove high-energy oscillations
    ! from the autocorrelation function and to increase numerical accuracy
    if (present(j)) then
      write(outfile,*) "  "
    else
      write(outfile,formtxt) 'Omega', 'SigmaT1', 'SigmaT2', 'SigmaOmega1', 'SigmaOmega2'
    endif
    do i = 1, nsteps
      if (dospec .and. omega0 < 0.0_dp) then
        if (present(j)) then
          write(outfile,formout2) j, omega(i) - 1 - omega0, cdabs(sigt1(i)), cdabs(sigt2(i)), cdabs(sigw1(i)), cdabs(sigw2(i))
        else
          write(outfile,formout) omega(i) - 1 - omega0, cdabs(sigt1(i)), cdabs(sigt2(i)), cdabs(sigw1(i)), cdabs(sigw2(i))
        endif
      else if (dospec .and. omega0 >= 0.0_dp) then
        if (present(j)) then
          write(outfile,formout2) j, omega(i) - omega0, cdabs(sigt1(i)), cdabs(sigt2(i)), cdabs(sigw1(i)), cdabs(sigw2(i))
        else
          write(outfile,formout) omega(i) - omega0, cdabs(sigt1(i)), cdabs(sigt2(i)), cdabs(sigw1(i)), cdabs(sigw2(i))
        endif
      else if (dohhg) then
        if (present(j)) then
!         write(outfile,formout2) j, omega(i)/0.05695_dp, cdabs(sigt1(i)), cdabs(sigt2(i)), cdabs(sigw1(i)), cdabs(sigw2(i))
          write(outfile,formout2) j, omega(i)/lfreq, cdabs(sigt1(i)), cdabs(sigt2(i)), cdabs(sigw1(i)), cdabs(sigw2(i))
!         write(outfile,formout2) j, omega(i)- omega0, cdabs(sigt1(i)), cdabs(sigt2(i)), cdabs(sigw1(i)), cdabs(sigw2(i))
        else
!         write(outfile,formout) omega(i)/0.05695_dp, cdabs(sigt1(i)), cdabs(sigt2(i)), cdabs(sigw1(i)), cdabs(sigw2(i))
          write(outfile,formout) omega(i)/lfreq, cdabs(sigt1(i)), cdabs(sigt2(i)), cdabs(sigw1(i)), cdabs(sigw2(i))
!         write(outfile,formout) omega(i)- omega0, cdabs(sigt1(i)), cdabs(sigt2(i)), cdabs(sigw1(i)), cdabs(sigw2(i))
        endif
      else
        write(outfile,formout) omega(i), cdabs(sigt1(i)), cdabs(sigt2(i)), cdabs(sigw1(i)), cdabs(sigw2(i))
      endif
    enddo

    close(outfile)

    deallocate(allsteps)

  end subroutine write_out

end module data_prep
