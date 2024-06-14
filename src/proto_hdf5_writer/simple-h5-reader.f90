module hdf5_stuff

  use hdf5
  implicit none
  integer, parameter  :: rank=1   ! rank of the data set
  integer, parameter  :: rank2=2   ! rank of the data set
  integer, parameter  :: rank3=3   ! rank of the data set
  integer, parameter  :: rank4=4   ! rank of the data set
  integer, parameter  :: dp=selected_real_kind(15,307)
  integer, parameter  :: i64=selected_int_kind(15)
  integer(i64), parameter :: basedim=6_i64
  integer(i64), parameter :: indexdim=4_i64

  contains

    subroutine simple_h5_io(filename, input_file, dname, dim0)

    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: input_file
    character(len=*), intent(in) :: dname
    integer(i64), intent(in)  :: dim0

    integer(hid_t) :: file_id, dataset_id, dataspace_id
    integer(hid_t) :: plist_id
    integer :: ierr
    integer(hsize_t) :: dims(rank)
    integer(hsize_t) :: cdims(rank)
    integer :: i, j, k, l
    integer(hsize_t)  :: data_dims(rank)
    integer :: inpt
    real(dp), allocatable :: buf(:) ! write buffer
    character(len=400) :: iom
    integer :: ios


    call h5open_f(ierr)

    call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, ierr)

    dims(:) = [dim0]
    call h5screate_simple_f(rank, dims, dataspace_id, ierr)
    call h5pcreate_f(h5p_dataset_create_f, plist_id, ierr)

!    not really sure a good way to chunk a list
!    cdims(:) = dim0/rank
    cdims(:) = dim0/50
    call h5pset_chunk_f(plist_id, rank, cdims, ierr)
    call h5pset_deflate_f(plist_id, 6, ierr)

    call h5dcreate_f(file_id, dname, h5t_native_double, dataspace_id, &
                    dataset_id, ierr, dcpl_id=plist_id)

    open(newunit=inpt, file=input_file, iostat=ios, iomsg=iom)
    if (ios /= 0) then
      write(*,*) "Fatal error!!!"
      write(*,*)  trim(iom)
      stop
    endif

    allocate(buf(dim0))
    do j = 1, dim0
      read(inpt, *) buf(j)
    enddo

    close(inpt)

    data_dims(:) = [dim0]
    call h5dwrite_f(dataset_id, h5t_native_double, buf, data_dims, ierr)

    deallocate(buf)

    ! just close all down
    call h5sclose_f(dataspace_id, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5dclose_f(dataset_id, ierr)
    call h5fclose_f(file_id, ierr)

  end subroutine

  subroutine write_1e_h5(filename, input_file, dname, dim0)

    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: input_file
    character(len=*), intent(in) :: dname
    integer(i64), intent(in)  :: dim0


    integer(hid_t) :: file_id, dataset_id, dataspace_id
    integer(hid_t) :: plist_id
    integer :: ierr
    integer(hsize_t) :: dims(rank3)
    integer(hsize_t) :: cdims(rank3)
    integer :: i, j, k, l
    integer(hsize_t)  :: data_dims(rank3)
    integer :: inpt
    integer :: chunk_size_loc
    integer :: compress_lvl_loc
    real(dp), allocatable :: buf(:,:,:) ! write buffer
    character(len=400) :: iom
    integer :: ios

    chunk_size_loc = 100000
    compress_lvl_loc = 6

    call h5open_f(ierr)

    call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, ierr)

    dims(:) = [basedim, dim0, dim0]
    call h5screate_simple_f(rank3, dims, dataspace_id, ierr)
    call h5pcreate_f(h5p_dataset_create_f, plist_id, ierr)

    cdims(:) = [basedim/rank3, dim0/rank3, dim0/rank3]
    call h5pset_chunk_f(plist_id, rank3, cdims, ierr)
    call h5pset_deflate_f(plist_id, 6, ierr)

! from the hande code
! this actually  isn't as efficient at least for the sizes i'm using
!    if (product(dims) > chunk_size_loc) then
!      cdims(1:rank3-1) = dims(1:rank3-1)
!      cdims(rank3) = chunk_size_loc/product(cdims(1:rank3-1))
!      call h5pset_chunk_f(plist_id, rank3, cdims, ierr)
!      call h5pset_deflate_f(plist_id, compress_lvl_loc, ierr)
!    end if



    call h5dcreate_f(file_id, dname, h5t_native_double, dataspace_id, &
                    dataset_id, ierr, dcpl_id=plist_id)

    open(newunit=inpt, file=input_file, iostat=ios, iomsg=iom)
    if (ios /= 0) then
      write(*,*) "Fatal error!!!"
      write(*,*)  trim(iom)
      stop
    endif

    allocate(buf(basedim,dim0,dim0))
    do k = 1, basedim
      do i = 1, dim0
        do j = 1, dim0
          read(inpt, *) buf(k, i, j)
        enddo
      enddo
    enddo

!    write(*,*) 'buf shape', shape(buf)
    close(inpt)

    data_dims(:) = [basedim, dim0, dim0]
    call h5dwrite_f(dataset_id, h5t_native_double, buf, data_dims, ierr)

    deallocate(buf)

    ! shutting it up, just close all down
    call h5sclose_f(dataspace_id, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5dclose_f(dataset_id, ierr)
    call h5fclose_f(file_id, ierr)

  end subroutine


!  subroutine write_2eindex_h5(filename, input_file, dname, dim0)
!
!    implicit none
!
!    character(len=*), intent(in) :: filename
!    character(len=*), intent(in) :: input_file
!    character(len=*), intent(in) :: dname
!    integer(i64), intent(in)  :: dim0
!
!
!    integer(hid_t) :: file_id, dataset_id, dataspace_id
!    integer(hid_t) :: plist_id
!    integer :: ierr
!    integer(hsize_t) :: dims(rank3)
!    integer(hsize_t) :: cdims(rank3)
!    integer :: i, j, k, l
!    integer(hsize_t)  :: data_dims(rank2)
!    integer :: inpt
!    integer :: chunk_size_loc
!    integer :: compress_lvl_loc
!    real(dp), allocatable :: buf(:,:,:) ! write buffer
!    character(len=400) :: iom
!    integer :: ios
!
!    chunk_size_loc = 100000
!    compress_lvl_loc = 6
!
!    call h5open_f(ierr)
!
!    call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, ierr)
!
!    dims(:) = [indexdim, dim0]
!    call h5screate_simple_f(rank2, dims, dataspace_id, ierr)
!    call h5pcreate_f(h5p_dataset_create_f, plist_id, ierr)
!
!    cdims(:) = [basedim/rank1, dim0/rank2]
!    call h5pset_chunk_f(plist_id, rank2, cdims, ierr)
!    call h5pset_deflate_f(plist_id, 6, ierr)
!
!! from the hande code
!! this actually  isn't as efficient at least for the sizes i'm using
!!    if (product(dims) > chunk_size_loc) then
!!      cdims(1:rank3-1) = dims(1:rank3-1)
!!      cdims(rank3) = chunk_size_loc/product(cdims(1:rank3-1))
!!      call h5pset_chunk_f(plist_id, rank3, cdims, ierr)
!!      call h5pset_deflate_f(plist_id, compress_lvl_loc, ierr)
!!    end if
!
!
!
!!    call h5dcreate_f(file_id, dname, h5t_native_double, dataspace_id, &
!    call h5dcreate_f(file_id, dname, h5t_native_double, dataspace_id, &
!                    dataset_id, ierr, dcpl_id=plist_id)
!
!    open(newunit=inpt, file=input_file, iostat=ios, iomsg=iom)
!    if (ios /= 0) then
!      write(*,*) "Fatal error!!!"
!      write(*,*)  trim(iom)
!      stop
!    endif
!
!    allocate(buf(indexdim,dim0))
!    do k = 1, basedim
!      do i = 1, dim0
!        read(inpt, *) buf(k, i)
!      enddo
!    enddo
!
!!    write(*,*) 'buf shape', shape(buf)
!    close(inpt)
!
!    data_dims(:) = [basedim, dim0, dim0]
!    call h5dwrite_f(dataset_id, h5t_native_double, buf, data_dims, ierr)
!
!    deallocate(buf)
!
!    ! shutting it up, just close all down
!    call h5sclose_f(dataspace_id, ierr)
!    call h5pclose_f(plist_id, ierr)
!    call h5dclose_f(dataset_id, ierr)
!    call h5fclose_f(file_id, ierr)
!
!  end subroutine


  subroutine write_2e_h5(filename, input_file, dname, dim0)

    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: input_file
    character(len=*), intent(in) :: dname
    integer(i64), intent(in)  :: dim0


    integer(hid_t) :: file_id, dataset_id, dataspace_id
    integer(hid_t) :: plist_id
    integer :: ierr
    integer(hsize_t) :: dims(rank4)
    integer(hsize_t) :: cdims(rank4)
    integer :: i, j, k, l
    integer(hsize_t)  :: data_dims(rank4)
    integer :: chunk_size_loc
    integer :: compress_lvl_loc
    integer :: inpt
    real(dp) :: trash
    real(dp), allocatable :: buf(:,:,:,:) ! write buffer
    character(len=400) :: iom
    integer :: ios
    chunk_size_loc = 100000
    compress_lvl_loc = 6


    call h5open_f(ierr)

    call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, ierr)

    dims(:) = [dim0,dim0,dim0,dim0]
    call h5screate_simple_f(rank4, dims, dataspace_id, ierr)
    call h5pcreate_f(h5p_dataset_create_f, plist_id, ierr)

    cdims(:) = dim0/rank4
    call h5pset_chunk_f(plist_id, rank4, cdims, ierr)
    call h5pset_deflate_f(plist_id, 9, ierr)

!    if (product(dims) > chunk_size_loc) then
!      cdims(1:rank4-1) = dims(1:rank4-1)
!      cdims(rank4) = chunk_size_loc/product(cdims(1:rank4-1))
!      call h5pset_chunk_f(plist_id, rank4, cdims, ierr)
!      call h5pset_deflate_f(plist_id, compress_lvl_loc, ierr)
!    end if



    call h5dcreate_f(file_id, dname, h5t_native_double, dataspace_id, &
                    dataset_id, ierr, dcpl_id=plist_id)

    open(newunit=inpt, file=input_file, iostat=ios, iomsg=iom)
    if (ios /= 0) then
      write(*,*) "Fatal error!!!"
      write(*,*)  trim(iom)
      stop
    endif

    allocate(buf(dim0,dim0,dim0,dim0))
    do k = 1, basedim
      do i = 1, dim0
        do j = 1, dim0
          read(inpt, *) trash
        enddo
      enddo
    enddo

    do i = 1, dim0
      do j = 1, dim0
        do k = 1, dim0
          do l = 1, dim0
            read(inpt,*) buf(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    close(inpt)

    data_dims(:) = dims(:)
    call h5dwrite_f(dataset_id, h5t_native_double, buf, data_dims, ierr)

    deallocate(buf)

    ! close resources
    call h5sclose_f(dataspace_id, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5dclose_f(dataset_id, ierr)
    call h5fclose_f(file_id, ierr)

  end subroutine

!  reopening routines

  subroutine read_h5_list(filename, dname, dim0, rbuf)

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dname
    integer(i64), intent(in)  :: dim0
    real(dp), intent(inout) :: rbuf(dim0)  ! read buffer

    integer(hid_t) :: file_id, dataset_id, dataspace_id ! identifiers
    integer(hid_t) :: plist_id ! property list identifier
    integer :: ierr
    integer(hsize_t) :: dims(rank) ! dimensions of data
    integer(hsize_t) :: cdims(rank) ! sizes of chunked data
    integer :: i, j, numfilt, k , l
    integer(hsize_t)  :: data_dims(rank) ! dimensions of data buffers
    call h5open_f(ierr)

    call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, ierr)
    call h5dopen_f(file_id, dname, dataset_id, ierr)

    call h5dget_create_plist_f(dataset_id, plist_id, ierr)

    data_dims(:) = [dim0]

    call h5dread_f(dataset_id, h5t_native_double, rbuf, data_dims, ierr)


    call h5dclose_f(dataset_id, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)
    return
  end subroutine read_h5_list

  subroutine read_1eints_h5(filename, dname, dim0, rbuf)

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dname
    integer(i64), intent(in)  :: dim0
    real(dp), intent(inout) :: rbuf(basedim,dim0,dim0)

    integer(hid_t) :: file_id, dataset_id, dataspace_id
    integer(hid_t) :: plist_id
    integer :: ierr
    integer :: i, j, numfilt, k
    integer(hsize_t)  :: data_dims(rank3) ! dimensions of data buffers
    call h5open_f(ierr)

    call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, ierr)
    call h5dopen_f(file_id, dname, dataset_id, ierr)

    call h5dget_create_plist_f(dataset_id, plist_id, ierr)

    data_dims(:) = [basedim,dim0,dim0]

    call h5dread_f(dataset_id, h5t_native_double, rbuf, data_dims, ierr)


    call h5dclose_f(dataset_id, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)
    return
  end subroutine


  subroutine read_2eints_h5(filename, dname, dim0, rbuf)

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dname
    integer(i64), intent(in)  :: dim0
    real(dp), intent(inout) :: rbuf(dim0,dim0,dim0,dim0)  ! read buffer

    integer(hid_t) :: file_id, dataset_id, dataspace_id ! identifiers
    integer(hid_t) :: plist_id ! property list identifier
    integer :: ierr
    integer :: i, j, numfilt, k
    integer(hsize_t)  :: data_dims(rank4) ! dimensions of data buffers
    call h5open_f(ierr)

    call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, ierr)
    call h5dopen_f(file_id, dname, dataset_id, ierr)

    call h5dget_create_plist_f(dataset_id, plist_id, ierr)

    data_dims(:) = [dim0,dim0,dim0,dim0]

    call h5dread_f(dataset_id, h5t_native_double, rbuf, data_dims, ierr)


    call h5dclose_f(dataset_id, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)
    return
  end subroutine

end module hdf5_stuff
