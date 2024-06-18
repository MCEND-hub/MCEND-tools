
!subroutine reinitstuff(newgrd,x,tmat,veeao)
subroutine reinitstuff(newgrd)
  use params
  use globalvars
  use inputvars
!  use moreadin
!  use utils
!  use hdf5
!  use hdf5_stuff
  implicit none
!  real(dp)            :: x(3,nrprime,nrprime) ! position
!  real(dp)            :: hmat(nrprime,nrprime,nrensp) !hamiltonian
!  real(dp)            :: tmat(nrprime,nrprime) !kinetic
!  real(dp)            :: veeao(nrprime,nrprime,nrprime,nrprime) !tw electron repulsion
  real(dp)            :: smoh(nrprime,nrprime) !
  real(dp)            :: temph
  integer             :: newgrd,i,j,k,l, ios, p, q
  character (len=2)   :: grdnum
  character (len=50)  :: inpformat
  character (len=255) :: workdir
  character (len=255) :: path2ints
  character (len=255) :: mofile
  character (len=255) :: iom
  integer(i64) :: total_dim
!  real(dp), allocatable :: allmos(:)
  total_dim = 6_i64*(nrprime**2) + (nrprime**4)

  ! get current working directory (cwd)

  call get_environment_variable("PWD", workdir)
  ! append cwd to integral directory
!  path2ints = trim(workdir)//trim(intdir)

  path2ints = trim(workdir)//'/'//trim(intdir)//'/'//trim(cmpdname)//'-data2-'
!  path2ints = trim(workdir)//trim(intdir)//trim(cmpdname)//'-data2-'

!  inpformat='(f18.15)'
!  inpformat='(e18.15)'
  inpformat='(e23.16)'

  write(*,*) "reading from", newgrd
!  write(grdnum,'(i3.3)') newgrd
  write(grdnum,'(i2.2)') newgrd

!  open(24, status='old', file=trim(path2ints)//trim(cmpdname)//'-r'//grdnum//'.dat')
!  if (readfromh5 /= 1) then

  open(24, status='old', file=trim(path2ints)//grdnum//'.dat', iostat=ios, iomsg=iom)
  if (ios /= 0) then
    write(*,*) 'Fatal error!'
    write(*,*) 'Iomessage: ', trim(iom)
    stop
  endif
  !
  do i=1,nrprime
    do j=1,nrprime
!     read(24,inpformat) hmat(i,j,newgrd)
      read(24,inpformat) temph
    enddo
  enddo

  do i=1,nrprime
    do j=1,nrprime
      read(24,inpformat) tmat(i,j)
    enddo
  enddo

  do i=1,nrprime
    do j=1,nrprime
      read(24,inpformat) x(1,i,j)
    enddo
  enddo

  do i=1,nrprime
    do j=1,nrprime
      read(24,inpformat) x(2,i,j)
    enddo
  enddo

  do i=1,nrprime
    do j=1,nrprime
      read(24,inpformat) x(3,i,j)
    enddo
  enddo

  ! Overlap matrix
  do i=1,nrprime
    do j=1,nrprime
      read(24,inpformat) smoh(i,j)
    enddo
  enddo

  do i=1,nrprime
    do j=1,nrprime
      do k=1,nrprime
        do l=1,nrprime
          read(24,inpformat) veeao(i,j,k,l)
        enddo
      enddo
    enddo
  enddo
  close(24)

!  do iensp=1,nrensp
!    if (iensp /= newgrd) then
!      write(grdnum,'(i2.2)') iensp
!      open(24, status='old', file=trim(path2ints)//'LiH-r'//grdnum//'.dat')
!      do i=1,nrprime
!        do j=1,nrprime
!          read(24,'(e20.13)') hmat(i,j,iensp)
!        enddo
!      enddo
!      close(24)
!    endif
!  enddo

  ! shift origin of energy axis by shifte

!  do iensp=1,nrensp
!  do i=1,nrprime
!     hmat(i,i,newgrd) = hmat(i,i,newgrd) + shifte/nel
!  enddo
!  enddo

!  elseif (readfromh5 == 1) then
!
!    allocate(allmos(total_dim))
!    call read_h5_list(trim(path2ints)//grdnum//'.h5', 'data_r'//grdnum, total_dim, allmos)
!
!    k = 1
!    do i=1, nrprime
!      do j=1, nrprime
!!        hmat(i,j,inigrd) = allmos(k)
!        temph = allmos(k)
!        k = k + 1
!      enddo
!    enddo
!
!    do i=1, nrprime
!      do j=1, nrprime
!        tmat(i,j) = allmos(k)
!        k = k + 1
!      enddo
!    enddo
!
!    do i=1, nrprime
!      do j=1, nrprime
!        x(1,i,j) = allmos(k)
!        k = k + 1
!      enddo
!    enddo
!
!    do i=1, nrprime
!      do j=1, nrprime
!        x(2,i,j) = allmos(k)
!        k = k + 1
!      enddo
!    enddo
!
!    do i=1, nrprime
!      do j=1, nrprime
!        x(3,i,j) = allmos(k)
!        k = k + 1
!      enddo
!    enddo
!
!    do i=1, nrprime
!      do j=1, nrprime
!        smoh(i,j) = allmos(k)
!        k = k + 1
!      enddo
!    enddo
!
!    do i=1, nrprime
!      do j=1, nrprime
!        do p=1, nrprime
!          do q=1, nrprime
!            veeao(i,j,p,q) = allmos(k)
!            k = k + 1
!          enddo
!        enddo
!      enddo
!    enddo
!
!    deallocate(allmos)
!
!  endif
end subroutine reinitstuff
