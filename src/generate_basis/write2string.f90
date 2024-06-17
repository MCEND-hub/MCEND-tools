module write2string
  implicit none

  contains

    elemental subroutine str2int(str,intval)
      character(len=*), intent(in) :: str
      integer, intent(out)         :: intval

      read(str,*) intval

    end subroutine str2int

    elemental subroutine int2str(str,intval)

      implicit none
      integer,intent(in)         :: intval
      character(len=*), intent(out) :: str
!      character(len=)

      ! LA. sorry, there isn't any other elegant way to do this
      ! (that I know of, if there is please tell me!!)
      if (1 <= intval .and. intval < 10) then
        write(str, '(i1)') intval
      else if (10 <= intval .and. intval < 100) then
        write(str, '(i2)') intval
      else if (100 <= intval .and. intval < 1000) then
        write(str, '(i3)') intval
      else if (1000 <= intval .and. intval < 10000) then
        write(str, '(i4)') intval
      else if (10000 <= intval .and. intval < 100000) then
        write(str, '(i5)') intval
      else if (100000 <= intval .and. intval < 1000000) then
        write(str, '(i6)') intval
      else if (1000000 <= intval .and. intval < 10000000) then
        write(str, '(i7)') intval
      else if (10000000 <= intval .and. intval < 100000000) then
        write(str, '(i8)') intval
      else if (100000000 <= intval .and. intval < 1000000000) then
        write(str, '(i9)') intval
      endif

    end subroutine int2str

end module
