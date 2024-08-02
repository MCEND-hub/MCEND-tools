program ft_transform
  use globalvars
  use globalconst
  use data_prep
  implicit none

  integer        :: j

  ! these are the calls to subroutines
  call readin()
  if (window) then
    do j=1, maxj
      call sliding_window(j)
      call transform
      call write_out(j)
    enddo
  else
    call prepare_data()
    call transform()
    call write_out()
  endif
  call dfftw_destroy_plan(plan1)
  call dfftw_destroy_plan(plan2)
  call system('rm -f lfreq.tmp count.txt')
  deallocate(sigt1)
  deallocate(sigt2)
  deallocate(sigw1)
  deallocate(sigw2)
  deallocate(omega)
  deallocate(time)

end program
