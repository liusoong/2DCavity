program twoDCavity


	use moduleParameters
	use moduleFullStep
	

	implicit none	
	
	
	real(8) 				:: u(0 : n, 0 : n + 1), v(0 : n + 1, 0 : n)
	real(8)				:: p(n, n)
	real(8) 				:: h, dt	
	integer				:: i

	h = 1.0 / n
	dt = 0.1 * h * h * Re
	
	u = 0.0
	v = 0.0
	p = 0.0
		
	do i = 1, maxIt
		print *, "Iteration: ", i
		print *, "----------------------------------------------------" 		
		call fullStep(h, dt, drivingV, u, v, p)
		!printData(u, v, p)
	enddo
	
end program twoDCavity


