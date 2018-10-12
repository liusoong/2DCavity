program twoDCavity


	use moduleParameters
	use moduleFullStep
	

	implicit none	
	
	
	real(8) 				:: u(0 : n, 0 : n + 1), v(0 : n + 1, 0 : n)
	real(8)				:: p(n, n)
	real(8) 				:: h, dt	
	integer				:: i

	h = 1.0 / n
	dt = 0.01 * h * h * Re
	
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
			print *, u((n-1)/2, n/2), v(n/2, (n-1)/2), p(n/2, n/2)
		end if
	end do
	
	!~ do i = 1, n
		!~ print *, i, p(i, i)	!!!!!!!!!!!!!!!!!!!!!!!!!
	!~ end do
	
end program twoDCavity


