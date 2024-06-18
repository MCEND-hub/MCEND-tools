!subroutine integral_writer(H, T, S, x, y, z, veeao, nrelem, nbf, outfilename, nrtwoeri)
subroutine integral_writer(H, T, S, x, y, z, rij, rijindex, outfilename,  nrelem, nrtwoeri)

  integer, parameter :: dp=selected_real_kind(15,307)
  integer, intent(in) :: nrelem
  integer, intent(in) :: nrtwoeri
!  integer, intent(in)   :: nbf
!  character(len=50), intent(in) :: outfilename
  character(len=*), intent(in) :: outfilename
  real(dp), intent(in)  :: H(nrelem), S(nrelem), T(nrelem)
  real(dp), intent(in)  :: x(nrelem), y(nrelem), z(nrelem)
  real(dp), intent(in)  :: rij(nrtwoeri)
  integer, intent(in)  :: rijindex(4, nrtwoeri)
  integer :: i, j
  integer :: f999

  open(newunit=f999,file=trim(outfilename))

  do i=1, nrelem
    write(f999,'(3(es26.16,1x))') H(i), T(i), S(i)
  enddo

  do i=1, nrelem
    write(f999,'(3(es26.16,1x))') x(i), y(i), z(i)
  enddo

  do i=1, nrtwoeri
    write(f999,'(4(i5),(es24.16))') rijindex(1,i), rijindex(2,i), rijindex(3,i), rijindex(4,i), rij(i)
  enddo

  close(f999)

end subroutine

subroutine hm_writer(H, T, S, x, y, z, outfilename, nrelem)

  integer, parameter :: dp=selected_real_kind(15,307)
  integer, intent(in) :: nrelem
  character(len=*), intent(in) :: outfilename
  real(dp), intent(in)  :: H(nrelem), S(nrelem), T(nrelem)
  real(dp), intent(in)  :: x(nrelem), y(nrelem), z(nrelem)
  integer :: i, j, k, l
  integer :: f999

  open(newunit=f999, file=trim(outfilename))

  do i=1, nrelem
    write(f999,'(3(es26.16,1x))') H(i), T(i), S(i)
  enddo

  do i=1, nrelem
    write(f999,'(3(es26.16,1x))') x(i), y(i), z(i)
  enddo

  close(f999)

end subroutine
