module gauss_quad
  implicit none
  private
  public :: gauss_leg_int, gauss_leg_int_divide, gauss_init
  integer, parameter :: dp = kind(1.d0)
  logical :: inited = .false.
  real(dp), allocatable :: gauss_x(:), gauss_w(:)
  real(dp) :: r1, r2

contains

  function gauss_leg_int_divide(f, a, b, n)
    !divides into n parts and use gauss_leg_int for each part
    real(dp) :: a, b, gauss_leg_int_divide
    intent(in) :: a, b, n
    real(dp) :: x, y
    integer :: i, n
    interface
      function f(x)
        import
        real(dp) :: f, x
        intent(in) :: x
      end function
    end interface

    gauss_leg_int_divide = 0
    do i = 0, n-1
      x = (b-a) / n * i + a
      y = (b-a) / n * (i+1) + a
      gauss_leg_int_divide = gauss_leg_int_divide + gauss_leg_int(f, x, y)
    enddo
  end function

  function gauss_leg_int(f, a, b)
    real(dp), intent(in) :: a, b
    real(dp) :: gauss_leg_int
    real(dp) :: x
    integer :: i
    interface
      function f(x)
        import
        real(dp) :: f, x
        intent(in) :: x
      end function
    end interface

    if (.not. inited) call gauss_init()
    gauss_leg_int = 0
    !need to map r1 to a, r2 to b
    do i = 1, 127
      x = (b-a) / (r2-r1) * (gauss_x(i) - r1) + a
      gauss_leg_int = gauss_leg_int + gauss_w(i) * f(x)
    enddo
    gauss_leg_int = gauss_leg_int * (b-a) / (r2-r1)
  end function

  subroutine gauss_init(datadir_in)
    integer :: i
    character(*), optional :: datadir_in
    character(1024) :: datadir = "./"
    if (present(datadir_in)) then
      datadir = trim(datadir_in)
    endif
    !use order 127
    allocate(gauss_x(127), gauss_w(127))
    open(unit=1234,file=trim(datadir) // "/leg_o127_x.txt")
    do i = 1, 127
      read(1234,*) gauss_x(i)
    enddo
    close(1234)
    open(unit=1234,file=trim(datadir) // "/leg_o127_w.txt")
    do i = 1, 127
      read(1234,*) gauss_w(i)
    enddo
    close(1234)
    open(unit=1234,file=trim(datadir) // "/leg_o127_r.txt")
    read(1234,*) r1
    read(1234,*) r2
    close(1234)
    inited = .true.
  end subroutine
end module
