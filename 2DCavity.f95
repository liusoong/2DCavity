program twoDCavity

	use moduleParameters
	use moduleFullStep

	implicit none	
		
	real(8) :: h = 1.0 / n	
	real(8) :: dt

	dt = 0.1 * h * h * Re
	
end program twoDCavity


