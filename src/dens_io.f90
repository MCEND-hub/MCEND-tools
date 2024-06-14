! f2py fortran subroutine
! a subroutine to handle the generation of density differences via cube files cube
! does the matrix subtraction then writes the cube file, needs the cube array input
! and the base cube file to subtract from.
! No longer needs the top lines, since they get mangled into some weird unicode form when read in
! no idea why.
!subroutine write_dens(cubefile, base_cube, cube_array, cbfile_head, headsize, ndim)
!subroutine write_dens(cubefile, ndim, base_cube, cube_array)
subroutine write_dens(cubefile, base_cube, cube_array, ndx, ndy, ndz)
  implicit none

  integer, parameter :: dp=selected_real_kind(15, 307)

  character(len=*), intent(in) :: cubefile
!  integer, intent(inout) :: ndim
!  integer, intent(inout) :: ndim(3)
  integer, intent(in) :: ndx
  integer, intent(in) :: ndy
  integer, intent(in) :: ndz
!  real(dp), intent(in) :: cube_array(ndim(1),ndim(2),ndim(3))
!  real(dp), intent(in) :: base_cube(ndim(1),ndim(2),ndim(3))
  real(dp), intent(in) :: cube_array(ndx,ndy,ndz)
  real(dp), intent(in) :: base_cube(ndx,ndy,ndz)

  real(dp), allocatable :: cube_diff(:,:,:)
  integer :: ocube, i, j

  open(newunit=ocube, file=cubefile, status='old', position='append')

!  if (.not. allocated(cube_diff)) allocate(cube_diff(ndim, ndim, ndim))
!  if (.not. allocated(cube_diff)) allocate(cube_diff(ndim(1), ndim(2), ndim(3)))
  if (.not. allocated(cube_diff)) allocate(cube_diff(ndx, ndy, ndz))

  cube_diff(:,:,:) = base_cube(:,:,:) - cube_array(:,:,:)

  do i=1, ndx
    do j=1, ndy
!      write(ocube,"(6(1pe13.5))",advance="no") cube_diff(i,j,1:ndim)
      write(ocube,"(6(1pe13.5))",advance="no") cube_diff(i,j,1:ndz)
      write(ocube,*)
    enddo
  enddo

  deallocate(cube_diff)

  close(ocube)

end subroutine
