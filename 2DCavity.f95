program twoDCavity


	use moduleParameters
	use moduleFullStep
	

	implicit none	
	
	
	real(8) 				:: u(0 : n, 0 : n + 1), v(0 : n + 1, 0 : n)
	real(8)				:: p(n, n)
	real(8) 				:: h, dt	
	integer				:: i, j

	h = 1.0 / n
	dt = 0.5 * min(0.25 * h * h * Re, h)
	
	u = 0.0
	v = 0.0
	p = 0.0
		
	do i = 1, maxIt
		!~ print *, "Iteration: ", i
		!~ print *, "----------------------------------------------------" 		
		call fullStep(h, dt, drivingV, u, v, p)
		!~ print *, p(n, n)
		!printData(u, v, p)
		if (mod(i, 100) == 0) then
			print *, i, u((n-1)/2, n/2), v(n/2, (n-1)/2), p(n/2, n/2)
		end if
	end do
	
	!~ do i = 1, n
		!~ print *, i, p(i, i)	!!!!!!!!!!!!!!!!!!!!!!!!!
	!~ end do
	
		
	! print out result
	open (unit = 1, file = "u")
	open (unit = 2, file = "v")	
	open (unit = 3, file = "p")
	do i = 1, n
		write(1, "(f15.8, 5x, f15.8)") (i-0.5) * h, u(n / 2, i)
		write(2, "(f15.8, 5x, f15.8)") (i-0.5) * h, v(i, n / 2)		
		do j = 1, n
			write(3, "(f15.8, 5x, f15.8, 5x, f15.8)") (i-0.5) * h, (j-0.5) * h, p(i, j)		
		end do
	end do
	close(1)
	close(2)
	close(3)
	
end program twoDCavity


