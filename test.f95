program test

	use moduleGaussSeidel
	
	implicit none
	
	real(kind = 8), parameter :: omega = 1.0
	real(kind = 8), dimension(4, 4) :: a = reshape((/ 10, -1, 2, 0, -1, 11, -1, 3, &
		2, -1, 10, -1, 0, 3, -1, 8 /), (/ 4, 4 /))
	real(kind = 8), dimension(4) :: b = (/ 6, 25, -11, 15 /)
	real(kind = 8), dimension(4) :: residu, x
	integer :: i, j, rc, iter = 10000	
	
	do i = 1, 4
		do j = 1, 4
			write(*, '(f8.3, t3)', advance = 'no'), a(i, j)
		end do
		write(*, *)
	end do
	
	call seidel(0, 4, a, b, omega, x, residu, iter, rc)
	
	print *, x
		
end program test