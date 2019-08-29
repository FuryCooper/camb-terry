program dl2cl
implicit none
real(8) :: dl(0:2508), cl(0:2508), temp(4)
integer :: i,j,k
character(80) :: dlname
call getarg(1,dlname)
open(1234,file=trim(dlname))
cl=0.d0
do i = 2, 2508
  read(1234,*) j, dl(i), temp(4)
  cl(i) = dl(i)*2*3.14159/i/(i+1)
enddo
close(1234)
open(1234,file=trim(dlname)//".cl")
do i = 0, 2508
  write(1234,*) cl(i)
enddo
close(1234)
end program dl2cl
