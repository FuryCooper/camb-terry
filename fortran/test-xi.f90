program test_xi
implicit none
type param
  real(8) :: xi_nu(3)
end type
type(param) :: CP
real(8) :: pi=3.14159265359d0
integer :: i,j
real(8), parameter :: s12sq = 0.304, s23sq = 0.51, s13sq = 0.0219 !terry
real(8), parameter :: s12 = sqrt(s12sq), s23 = sqrt(s23sq), s13 = sqrt(s13sq) !terry
real(8), parameter :: c12 = sqrt(1.d0-s12sq), c23 = sqrt(1.d0-s23sq), c13 = sqrt(1.d0-s13sq) !terry
real(8), parameter :: t12 = s12/c12, t23 = s23/c23, t13 = s13/c13 !terry
real(8), parameter :: r23 = ((1 - t13**2*(1 - t12**2))*(c23**2 - s23**2) + 2*t12*t13**2*s13*2*s23*c23)/& !terry
  ((1 - t12**2)*(c23**2 - s23**2) - 2*s13*t12*2*s23*c23) !terry
!terry: calculate CP%xi_nu(1) CP%xi_nu(2) from CP%xi_nu(3)
do j = 0, 20
  CP%xi_nu(3)=1.d0/10*j
  CP%xi_nu(2) = r23 * CP%xi_nu(3) * (CP%xi_nu(3)**2 + pi**2)/pi**2
  do i = 1, 10
    CP%xi_nu(2) = r23 * CP%xi_nu(3) * (CP%xi_nu(3)**2 + pi**2)/(CP%xi_nu(2)**2 + pi**2)
  enddo
  CP%xi_nu(1) = -(s12sq * CP%xi_nu(2) * (CP%xi_nu(2)**2 + pi**2) + t13**2 * CP%xi_nu(3) * (CP%xi_nu(3)**2 + pi**2))/c12**2/pi**2
  do i = 1, 10
    CP%xi_nu(1) = -(s12sq * CP%xi_nu(2) * (CP%xi_nu(2)**2 + pi**2) + t13**2 * CP%xi_nu(3) &
      * (CP%xi_nu(3)**2 + pi**2))/c12**2/(CP%xi_nu(1)**2 + pi**2)
  enddo
  print *, CP%xi_nu
enddo
end program
