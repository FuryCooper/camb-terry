program aaa
implicit none
real(8) :: xi(3), pi
pi=3.14159265359d0
print *, test_sum([0.133350170372944d0, 0.291022830490415d0, 0.166955800000000d0])
print *, test_sum([2.337006068546993d-002, 2.687031677738856d-002, 0.3570307d0])
contains
function test_sum(xi)
  real(8) :: xi(3),test_sum
  integer :: i
  test_sum=0
  do i=1, 3
    test_sum=test_sum+(xi(i)/pi)**2*(2.d0+(xi(i)/pi)**2)
  enddo
  test_sum=test_sum*15.d0/7.d0
end function
end program
