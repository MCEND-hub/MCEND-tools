program write_h5
  use hdf5
  use hdf5_stuff

  implicit none

  character(len=4)      :: cmpd
  character(len=400)    :: path
  character(len=2)      :: grdnum
  integer(i64)          :: nrprime, total_dim
  real(dp), allocatable :: allmos(:)
  real(dp), allocatable :: onemat(:,:,:), veemat(:,:,:,:)
!  real(dp), allocatable :: veeao(:,:)
  integer               :: i, j, k, l
  integer               :: nrensp, iensp
!  nrprime = 61_i64
!  nrprime = 41_i64
!  nrprime = 29_i64
  nrprime = 17_i64
  nrensp = 1


  total_dim = basedim*(nrprime**2) + nrprime**4
!  onee_dim = 6_i64*(nrprime**2)

!  cmpd = 'H2'
!  path = '/home/lucas/MCEND-v2.0/integrals_x2H2-bf17/'
  cmpd = 'N2'
  path = ''
!  path = ''
!  path='/home/lea0105/mcend/h2_atz_runs/integrals_H2-bf56/'
!  path='/home/lea0105/mcend/MCEND-v2.0/integrals_H2-bf56/'
!  cmpd = 'H2'
!  cmpd = 'LiH'
  ! note the size of the array dimension absolutely must be
  ! set as _i64 (i.e. 64 bit integer), or else it will cause
  ! segmentation faults

!  allocate(allmos(total_dim))
!  allocate(onemat(basedim, nrprime, nrprime))
!  allocate(veemat(nrprime, nrprime, nrprime, nrprime))
!  allocate(veeao(nrprime**2, nrprime**2))
!  allocate(hmat2(nrprime, nrprime, 4))
!  open(20,file='loophmat.dat')
!  open(25,file='reshapehmat.dat')

  do iensp=1, nrensp
    write(grdnum, '(i2.2)') iensp
!    write(*,*) trim(cmpd)//'-1e-data-'//grdnum//'.h5'

    call simple_h5_io(trim(path)//trim(cmpd)//'-data2-'//grdnum//'.h5', &
               trim(path)//trim(cmpd)//'-data2-'//grdnum//'.dat', 'data_r'//grdnum, total_dim)

!    call read_h5_list(trim(path)//trim(cmpd)//'-data2-'//grdnum//'.h5', 'data_r'//grdnum, total_dim, allmos)
!
!    call write_1e_h5(trim(path)//trim(cmpd)//'-1e-data-'//grdnum//'.h5', &
!               trim(path)//trim(cmpd)//'-data2-'//grdnum//'.dat', '1e-data_r'//grdnum, nrprime)
!
!    call write_2e_h5(trim(path)//trim(cmpd)//'-2e-data-'//grdnum//'.h5', &
!               trim(path)//trim(cmpd)//'-data2-'//grdnum//'.dat', '2e-data_r'//grdnum, nrprime)


!    call read_1eints_h5(trim(path)//trim(cmpd)//'-1e-data-'//grdnum//'.h5', '1e-data_r'//grdnum, nrprime, onemat)
!    call read_2eints_h5(trim(path)//trim(cmpd)//'-2e-data-'//grdnum//'.h5', '2e-data_r'//grdnum, nrprime, veemat)

!    do k=1, 2
!      do i = 1, 2
!        do j = 1, 2
!          write(*,'(e23.16)') onemat(k,i,j)
!        enddo
!      enddo
!    enddo
!onemat(1,:,:) = hmat(:,:,inigrd)
!onemat(2,:,:) = tmat(:,:)
!onemat(3:5,:,:) = xxmat(1:3,:,:)
!onemat(6,:,:) = smoh(:,:)
!    write(*,'(e23.16)') veemat(1,2,3,6)

  enddo

!  close(20)
!  close(25)


!  deallocate(allmos)
!  deallocate(hmat,veemat,veeao)
!  deallocate(onemat,veemat)

end program write_h5
