module PolyLog
  !this algorithm to approximate the PolyLog function encountered in fermi dirac integrals is based on Crandall
  !reference: https://www.reed.edu/physics/faculty/crandall/papers/Polylog.pdf
  !only calculates negative real arguments and n=3, 5, 7 since I only need these values, results are real
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: L
  logical :: setP = .false.
  !Bernoulli polynomial coefficient of x^n, n from 1 to 7
  real(dp), parameter :: Bcoeff(7, 7) = &
    [[1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0], &
    [-1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0], &
    [1.d0/2.d0, -3.d0/2.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0], &
    [0.d0, 1.d0, -2.d0, 1.d0, 0.d0, 0.d0, 0.d0], &
    [-1.d0/6.d0, 0.d0, 5.d0/3.d0, -5.d0/2.d0, 1.d0, 0.d0, 0.d0], &
    [0.d0, -1.d0/2.d0, 0.d0, 5.d0/2.d0, -3.d0, 1.d0, 0.d0], &
    [1.d0/6.d0, 0.d0, -7.d0/6.d0, 0.d0, 7.d0/2.d0, -7.d0/2.d0, 1.d0]] 
  real(dp), parameter :: pi = 3.14159265358979324d0
  real(dp), parameter :: twopi = 6.28318530717958648d0
  real(dp), parameter :: pio2 = 1.57079632679489662d0
  !factorial upto 7!
  real(dp), parameter :: fac(7) = [1.d0, 2.d0, 6.d0, 24.d0, 120.d0, 720.d0, 5040.d0]
  !zeta function, zeta(1) is singular but just set to 0 here
  real(dp), parameter :: zeta(-100:7) = &
    [0.d0, 2.8382249570693707d76, 0.d0, -1.15490239239635197d74, 0.d0, 4.89623270396205532d71, &
    0.d0, -2.16456348684351856d69, 0.d0, 9.98755741757275307d66, 0.d0, -4.81432188740457694d64, &
    0.d0, 2.42673403923335241d62, 0.d0, -1.28045468879395088d60, 0.d0, 7.07987744084945806d57, 0.d0, &
    -4.10670523358102125d55, 0.d0, 2.50194790415604628d53, 0.d0, -1.60293645220089654d51, 0.d0, &
    1.08136354499716547d49, 0.d0, -7.69198587595071352d46, 0.d0, 5.77753863427704318d44, 0.d0, &
    -4.58929744324543322d42, 0.d0, 3.86142798327052589d40, 0.d0, -3.44737825582780539d38, 0.d0, &
    3.27156342364787163d36, 0.d0, -3.30660898765775767d34, 0.d0, 3.56665820953755561d32, 0.d0, &
    -4.11472887925579787d30, 0.d0, 5.08906594686622897d28, 0.d0, -6.7645882379292821d26, 0.d0, &
    9.68995788746359407d24, 0.d0, -1.50017334921539287d23, 0.d0, 2.51804719214510957d21, 0.d0, &
    -4.59798883436565035d19, 0.d0, 9.16774360319533078d17, 0.d0, -2.00403106565162527d16, 0.d0, &
    4.82414483548501704d14, 0.d0, -1.28508504993050833d13, 0.d0, 3.80879311252453688d11, 0.d0, &
    -1.26357247959166667d10, 0.d0, 4.72384867721629902d8, 0.d0, -2.00526957966880789d7, 0.d0, &
    974936.823850574713d0, 0.d0, -54827.5833333333333d0, 0.d0, 3607.5105463980464d0, 0.d0, &
    -281.460144927536232d0, 0.d0, 26.4562121212121212d0, 0.d0, -3.05395433027011974d0, 0.d0, &
    0.443259803921568627d0, 0.d0, -0.0833333333333333333d0, 0.d0, 0.0210927960927960928d0, 0.d0, &
    -0.00757575757575757576d0, 0.d0, 0.00416666666666666667d0, 0.d0, -0.00396825396825396825d0, 0.d0, &
    0.00833333333333333333d0, 0.d0, -0.0833333333333333333d0, -0.5d0, 0.d0, &
    1.64493406684822644d0, 1.20205690315959429d0, 1.08232323371113819d0, &
    1.03692775514336993d0, 1.01734306198444914d0, 1.00834927738192283d0]

  private
  public :: polyLogLi, setPolyLogPrec
  contains

    subroutine setPolyLogPrec(D)
      integer, intent(in) :: D

      L = ceiling(D * log(10.d0) / log(2.d0))
      if (3-L < -100) then
        print *, "hardcoded zeta values are not enough!"
        stop
      endif
      setP = .true.
    end subroutine setPolyLogPrec

    function polyLogLi(n, z) result(val)
      integer :: n, i
      real(dp) :: z, val
      intent(in) :: n, z

      if (.not. setP) then
        print *, "polylog precision not set! stopping"
        stop
      endif
      if (z == 1.d0) then
        val = zeta(n)
      else if (z == -1.d0) then
        val = -1 * (1.d0 - 2.d0 ** (1-n)) * zeta(n)
      else if (abs(z) < 0.5d0) then
        val = Fn0(n, L, z)
      else if (abs(z) > 2.d0) then
        val = ReGn(n, L, z) - (-1)**n * Fn0(n, L, 1.d0 / z)
      else
        val = ReFn1(n, L, z)
      endif
    end function polyLogLi

    function Fn0(n, L, z)
      integer :: L, i, n
      real(dp) :: z, Fn0
      intent(in) :: L, z, n

      Fn0 = 0
      do i = 1, L
        Fn0 = Fn0 + z ** i / (i * 1.d0) ** n
      enddo
    end function Fn0

    function ReGn(n, L, z)
      integer :: L, i, n
      real(dp) :: z, ReGn
      intent(in) :: L, z, n

      ReGn = twopi ** n / fac(n) * (-1) ** ((n-1)/2) * sum(Bcoeff(:,n) * Im_logzoi(7, z))
    end function ReGn

    function Im_logzoi(k, z)
      real(dp) :: z, Im_logzoi, r, th
      dimension :: Im_logzoi(k)
      intent(in) :: z, k
      integer :: k, i

      call get_logz_rth(z, r, th)
      do i = 1, k
        Im_logzoi(i) = (r / twopi) ** i * sin(i * (th - pio2))
      enddo
    end function

    subroutine get_logz_rth(z, r, th)
      real(dp) :: z, r, th
      intent(in) :: z
      intent(out) :: r, th

      !only works for negative, real z!
      if (z>=0.d0) then
        print *, "logz error: non-negative z", z
        stop
      endif
      r = sqrt(log(-z)**2 + pi**2)
      th = acos(log(-z) / r)
    end subroutine get_logz_rth

    function ReFn1(n, L, z)
      integer :: L, i, n
      real(dp) :: z, ReFn1, r, th, ifac, ri, H
      intent(in) :: L, z, n

      call get_logz_rth(z, r, th)
      ifac = 1.d0
      ri = 1.d0
      H=0
      ReFn1 = zeta(n)
      do i = 1, n-2
        H = H + 1.d0 / i
        ri = ri * r
        ifac = ifac * i
        ReFn1 = ReFn1 + zeta(n-i) / ifac * ri * cos(i * th)
      enddo
      ri = ri * r
      ifac = ifac * (n-1)
      ReFn1 = ReFn1 + 1.d0 / ifac * (ri * cos((n-1) * th) * (H + 1.d0 / (n-1) - log(r)) + &
        ri * sin((n-1) * th) * (th - pi))
      do i = n, L
        ri = ri * r
        ifac = ifac * i
        ReFn1 = ReFn1 + zeta(n-i) / ifac * ri * cos(i * th)
      enddo
    end function ReFn1
end module PolyLog

!!a small main program for testing
!program test
!  use PolyLog
!  implicit none
!  integer :: i, j
!  real(8) :: Li
!  call setPolyLogPrec(10)
!  print *, 'test Li 3 -1', polyLogLi(3, -1.d0)
!  do j = 3, 7, 2
!    do i = -3, 3
!      Li = polyLogLi(j, -exp(1.d0 * i))
!      print *, j, i, Li
!    enddo
!  enddo
!end program test
