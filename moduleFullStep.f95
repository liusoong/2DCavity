module moduleFullStep


	use moduleParameters
	use moduleMultiGrid
	

	implicit none
	
	
	contains	
	
	
	subroutine fullStep(h, dt, drivingV, u, v, p)
	
		real(8), 	intent(in) 	:: h, dt, drivingV
		real(8)				:: p(n, n)
		real(8) 				:: u(0 : n, 0 : n + 1), v(0 : n + 1, 0 : n)
		real(8) 				:: uTemp(0 : n, 0 : n + 1), vTemp(0 : n + 1, 0 : n)
		integer 				:: i, j
		
		call uvTemp(h, dt, drivingV, u, v, uTemp, vTemp)	
						
		call poisson(h, dt, uTemp, vTemp, p)

		!~ do i = 1, n
				!~ print *, i, p(i, i)	!!!!!!!!!!!!!!!!!!!!!!!!!
		!~ end do

		call correctUV(h, dt, uTemp, vTemp, p, u, v)

		!~ print *, u(n-1, n), v(n, n-1), p(n, n)
		
	end subroutine fullStep

	
	subroutine uvTemp(h, dt, drivingV, u, v, uTemp, vTemp)	

		real(8),	intent(in) 		:: h, drivingV, dt
		real(8), 	intent(out)		:: uTemp(0 : n, 0 : n + 1), vTemp(0 : n + 1, 0 : n)
		real(8) 					:: u(0 : n, 0 : n + 1), v(0 : n + 1, 0 : n)		
		real(8) 					:: advU(0 : n, 0 : n + 1), advV(0 : n + 1, 0 : n)
		integer 					:: i, j
		
		uTemp = 0.0
		vTemp = 0.0
		
		call advectionUV(h, drivingV, u, v, advU, advV)	
		
		do i = 1, n - 1
			do j = 1, n
				uTemp(i, j) = u(i, j) + dt * ((u(i + 1, j) + u(i - 1, j) &
					+ u(i, j + 1) + u(i, j - 1) - 4 * u(i, j)) / (h ** 2 * Re) - advU(i, j))
				!print *, i, j, uTemp(i, j)	!!!!!!!!!!!!!!!!!!!!!!!!!
			end do
		end do		
		
		do i = 1, n
			do j = 1, n - 1
				vTemp(i, j) = v(i, j) + dt * ((v(i + 1, j) + v(i - 1, j) &
					+ v(i, j + 1) + v(i, j - 1) - 4 * v(i, j)) / (h ** 2 * Re) - advV(i, j))
				!print *, i, j, vTemp(i, j)	!!!!!!!!!!!!!!!!!!!!!!!!!
			end do
		end do
		
	end subroutine uvTemp
	
	
	subroutine advectionUV(h, drivingV, u, v, advU, advV)
	
		real(8),	intent(in) 	:: h, drivingV
		real(8) 				:: u(0 : n, 0 : n + 1), v(0 : n + 1, 0 : n)
		real(8), 	intent(out) 	:: advU(0 : n, 0 : n + 1), advV(0 : n + 1, 0 : n)
		integer 				:: i, j		
		real(8)				:: uatv, vatu
		
		advU = 0.0
		advV = 0.0
		
		forall (i = 0 : n) 
			u(i, 0) = - u(i, 1)
			u(i, n + 1) = 2.0 * drivingv - u(i, 1)
		end forall	
		
		forall (j = 0 : n)
			v(0, j) = - v(1, j)
			v(n + 1, j) = - v(n + 1, j)
		end forall
		
		do i = 1, n - 1
			do j = 1, n
				vatu = 0.25 * (v(i, j - 1) + v(i + 1, j - 1) + v(i, j) + v(i + 1, j))
				if (u(i, j) > 0.0) then
					advU(i, j) = u(i, j) * ((u(i, j) - u(i - 1, j)) / h)
				else
					advU(i, j) = u(i, j) * ((u(i + 1, j) - u(i, j)) / h)
				end if
				if (vatu > 0.0) then
					advU(i, j) = advu(i, j) + vatu * ((u(i, j) - u(i, j - 1)) / h)
				else
					advU(i, j) = advu(i, j) + vatu * ((u(i, j + 1) - u(i, j)) / h)
				end if
			end do
		end do
		
		do i = 1, n
			do j = 1, n - 1
				uatv = 0.25 * (u(i - 1, j) + u(i, j) + u(i - 1, j + 1) + u(i, j + 1))
				if (uatv > 0.0) then
					advV(i, j) = uatv * ((v(i, j) - v(i - 1, j)) / h)
				else
					advV(i, j) = uatv * ((v(i + 1, j) - v(i, j)) / h)
				end if					
				if (v(i, j) > 0.0) then
					advV(i, j) = advv(i, j) + v(i, j) * ((v(i, j) - v(i, j - 1)) / h)
				else
					advV(i, j) = advv(i, j) + v(i, j) * ((v(i, j + 1) - v(i, j)) / h)
				end if
			end do
		end do
		
	end subroutine advectionUV
	
	
	subroutine poisson(h, dt, uTemp, vTemp, p)
	
		real(8),	intent(in) 			:: h, dt
		real(8), 	intent(in) 			:: uTemp(0 : n, 0 : n + 1), vTemp(0 : n + 1, 0 : n)
		real(8)						:: p(1 : n, 1 : n)
		real(8)						:: rhs(1 : n, 1 : n)
		integer						:: i, j
		
		p = 0.
		
		Do i = 1, n
			do j = 1, n
				rhs(i, j) = (uTemp(i, j) - uTemp(i - 1, j) + vTemp(i, j) - vTemp(i, j - 1)) / h
			end do
		end do
		
		rhs = rhs / dt
		
		!~ Do i = 1, n					!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!~ do j = 1, n				!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!~ print *, i, j, rhs(i, j)		!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!~ end do					!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!~ end do						!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		call multiGridV(h, p, rhs)
	
	end subroutine poisson
	
	
	subroutine correctUV(h, dt, uTemp, vTemp, p, uNew, vNew)
	
		real(8),	intent(in) 			:: h, dt
		real(8), 	intent(in) 			:: uTemp(0 : n, 0 : n + 1), vTemp(0 : n + 1, 0 : n)
		real(8), 	intent(in) 			:: p(1 : n, 1 : n)	
		real(8), 	intent(OUT) 			:: uNew(0 : n, 0 : n + 1), vNew(0 : n + 1, 0 : n)
		integer						:: i, j
	
		do i = 1, n - 1
			do j = 1, n
				uNew(i, j) = uTemp(i, j) - dt * ((p(i + 1, j) - p(i, j)) / h)
			end do
		end do
	
		do i = 1, n
			do j = 1, n - 1
				vNew(i, j) = vTemp(i, j) - dt * ((p(i, j + 1) - p(i, j)) / h)
			end do
		end do
	
	end subroutine correctUV
	
	
end module moduleFullStep


