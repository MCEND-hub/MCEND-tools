module hdf5_stuff

  use hdf5
  use params
  implicit none
  integer, parameter  :: rank=1   ! rank of the data set
!  integer, parameter  :: dp=selected_real_kind(15,307)
!  integer, parameter  :: i64=selected_int_kind(15)

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
    integer :: i, j, k 
    integer(hsize_t)  :: data_dims(rank)
    integer :: inpt
!    real(dp) :: buf(dim0)
    real(dp), allocatable :: buf(:) ! write buffer
    character(len=400) :: iom
    integer :: ios


    call h5open_f(ierr)

    call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, ierr)

    dims(:) = [dim0]
    call h5screate_simple_f(rank, dims, dataspace_id, ierr)
    call h5pcreate_f(h5p_dataset_create_f, plist_id, ierr)

!    not really sure a good way to chunk a list
    cdims(:) = dim0/10
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
    buf(:) = 0.0_dp 

    allocate(buf(dim0))
    do j = 1, dim0
      read(inpt, *) buf(j)
    enddo 

    close(inpt)

    data_dims(:) = [dim0]
    call h5dwrite_f(dataset_id, h5t_native_double, buf, data_dims, ierr)

    deallocate(buf)

    ! close resources
    call h5sclose_f(dataspace_id, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5dclose_f(dataset_id, ierr)
    call h5fclose_f(file_id, ierr)

  end subroutine
!!!!!!!!!!!!!!!REOPEN!!!!!!!!!!!!!!!!!!!!!!1

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
    integer :: i, j, numfilt, k 
    integer(hsize_t)  :: data_dims(rank) ! dimensions of data buffers
    call h5open_f(ierr)

    call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, ierr)
    call h5dopen_f(file_id, dname, dataset_id, ierr)

  ! retrieve filter information. 
    call h5dget_create_plist_f(dataset_id, plist_id, ierr)
    
    data_dims(:) = [dim0]

    call h5dread_f(dataset_id, h5t_native_double, rbuf, data_dims, ierr)

!    do i=1, 10
!      write(*,'(e23.16)') rbuf(i)
!    enddo

!    write(*,*) rbuf(3)
      
    call h5dclose_f(dataset_id, ierr)
    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)
    return 
  end subroutine read_h5_list 

end module hdf5_stuff
