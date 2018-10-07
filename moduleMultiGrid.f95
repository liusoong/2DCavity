module moduleMultiGrid


	use moduleParameters
	

	implicit none

	
	contains	
		
	
	! Get 1D array index for (i, j)
	integer function getIndex1D(i, j, level) result (index1D)
	
		integer, intent(in)	:: i, j, level
		integer				:: thisLevel, nTemp
		Integer				:: before
		
		before = 0
		
		Do thisLevel = 0, level - 1
			nTemp = n / (2 ** thisLevel)
			before = before +nTemp ** 2
		end do
		
		thisLevel = level
		nTemp = n / (2 ** thisLevel)
		
		!print *, i, j, level, before,  (i - 1) * nTemp + j + before	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
		
		index1D = (i - 1) * nTemp + j + before
	
	end function getIndex1D
	
	! Get 1D array index for (i-1, j)
	integer function getIndex1DUp(i, j, level) result (Index1DForUp)
	
		integer, intent(in)	:: i, j, level
		integer				:: nTemp
		
		nTemp = n / (2 ** level)
		
		if(j == nTemp) then
			Index1DForUp = 0
		else
			Index1DForUp = getIndex1D(i, j + 1, level)
		end if
	
	end function getIndex1DUp
	
	! Get 1D array index for (i+1, j)
	integer function getIndex1DDown(i, j, level) result (index1DForDown)
	
		integer, intent(in)	:: i, j, level	
		
		if(j == 1) then
			index1DForDown = 0
		else
			index1DForDown = getIndex1D(i, j - 1, level)
		end if
	
	end function getIndex1DDown
	
	! Get 1D array index for (i, j-1)
	integer function getIndex1DLeft(i, j, level) result (index1DForLeft)
	
		integer, intent(in)	:: i, j, level
		
		if(i == 1) then
			index1DForLeft = 0
		else
			index1DForLeft = getIndex1D(i - 1, j, level)
		end if
	
	end function getIndex1DLeft
	
	! Get 1D array index for (i, j+1)
	integer function getIndex1DRight(i, j, level) result (index1DForRight)
	
		integer, intent(in)	:: i, j, level
		integer				:: nTemp
	
		nTemp = n / (2 ** level)		
		
		if(i == nTemp) then
			index1DForRight = 0
		else
			index1DForRight = getIndex1D(i + 1, j, level)
		end if
	
	end function getIndex1DRight
	
	
	subroutine relaxGS(lhs1D, rhs1D, n1D, level, nIter)
			
		integer, intent(in)	:: n1D
		real(8)				:: lhs1D(1 : n1D)
		real(8),	intent(in)	:: rhs1D(1 :  n1D)
		integer, intent(in)	:: level, nIter
		real(8) 				:: h = 1.0 / n	
		integer				:: nTemp
		integer				:: iIter, i, j, index1D, index1DUp, index1DDown, index1DLeft, index1DRight
		real(8)				:: up, down, left, right
		integer				:: coef
		
		
		nTemp = n / (2 ** level)
		print *, "nTemp = ", nTemp	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		do iIter = 1 , nIter
			do i = 1 , nTemp
				Do j = 1 , nTemp
					coef = 0
					! Get 1D array indices
					index1D = getIndex1D(i, j, level)					
					index1DUp =getIndex1DUp(i, j, level)					
					index1DDown =getIndex1DDown(i, j, level)
					index1DLeft =getIndex1DLeft(i, j, level)
					index1DRight =getIndex1DRight(i, j, level)
					! Get (i-1,j), (i+1,j), (i, j-1) and (i, j+1) values
					if (index1DUp == 0) then
						up = 0.
					else
						up = lhs1D(index1DUp)
						coef = coef + 1
					end if
					if (index1DDown == 0) then
						down = 0.
					else
						down = lhs1D(index1DDown)						
						coef = coef + 1
					end if
					if (index1DLeft == 0) then
						left = 0.
					else
						left = lhs1D(index1DLeft)
						coef = coef + 1
					end if
					if (index1DRight == 0) then
						right = 0.
					else
						right = lhs1D(index1DRight)
						coef = coef + 1
					end if
					! Calculate new (i,j) value
					lhs1D(index1D) = ((up + down + left + right) - rhs1D(index1D) * h**2) / coef
					!~ if (i == nTemp .and. j == nTemp) then
						!~ print *, i, j, index1D, h, lhs1D(index1D), rhs1D(index1D)	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					!~ end if
				end do
			end do		
		end do
	
	end subroutine relaxGS
	
	
	subroutine residue(lhs1D, rhs1d, res1d, n1D, level)
						
		integer, intent(in)	:: n1D
		real(8),	intent(in)	:: lhs1D(1 :  n1D), rhs1D(1 :  n1D)
		real(8)				:: res1D(1 : n1D)
		integer, intent(in)	:: level
		real(8) 				:: h = 1.0 / n	
		integer				:: nTemp
		integer				:: i, j, index1D, index1DUp, index1DDown, index1DLeft, index1DRight
		real(8)				:: up, down, left, right
		
		nTemp = n / (2 ** level)
		
		do i = 1 , nTemp
			Do j = 1 , nTemp
				! Get 1D array indices
				index1D = getIndex1D(i, j, level)					
				index1DUp =getIndex1DUp(i, j, level)					
				index1DDown =getIndex1DDown(i, j, level)
				index1DLeft =getIndex1DLeft(i, j, level)
				index1DRight =getIndex1DRight(i, j, level)
				! Get (i-1,j), (i+1,j), (i, j-1) and (i, j+1) values
				if (index1DUp == 0) then
					up = lhs1D(index1D)
				else
					up = lhs1D(index1DUp)
				end if
				if (index1DDown == 0) then
					down = lhs1D(index1D)
				else
					down = lhs1D(index1DDown)
				end if
				if (index1DLeft == 0) then
					left = lhs1D(index1D)
				else
					left = lhs1D(index1DLeft)
				end if
				if (index1DRight == 0) then
					right = lhs1D(index1D)
				else
					right = lhs1D(index1DRight)
				end if
				! Calculate residue
				res1D(index1D) = rhs1D(index1D) - ((up + down + left + right) - 4 * lhs1D(index1D)) / h**2
				!~ if (i == j) then
					!~ print *, i, j, index1D, h, res1D(index1D)	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!~ end if
			end do
		end do		
	
	end subroutine residue
	
		
	subroutine restrict(res1D, rhs1D, n1D, level)	

		integer, intent(in)	:: n1D
		real(8),	intent(in)	:: res1D(1 :  n1D)
		real(8)				:: rhs1D(1 : n1D)
		integer, intent(in)	:: level
		integer				:: nTemp
		integer				:: i, j, index1DCoarse, index1DFine

		nTemp = n / 2 ** (level + 1)
		Print *, "nTemp = ", nTemp
		
		do i = 1 , nTemp
			Do j = 1 , nTemp
				index1DCoarse = getIndex1D(i, j, level + 1)				
				!~ if (i == j) then
					!~ Print *, "index1DCoarse = ", index1DCoarse	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!~ end if
				! Lower left
				index1DFine = getIndex1D(2 * i - 1, 2 * j - 1, level)
				rhs1D(index1DCoarse) = res1D(index1DFine) / 4						
				!~ if (i == j) then
					!~ Print *, "index1DFine = ", index1DFine, "res1D = ", res1D(index1DFine) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!~ end if
				! Upper left
				index1DFine = getIndex1D(2 * i - 1, 2 * j, level)
				rhs1D(index1DCoarse) = rhs1D(index1DCoarse) + res1D(index1DFine) / 4					
				!~ if (i == j) then
					!~ Print *, "index1DFine = ", index1DFine, "res1D = ", res1D(index1DFine) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!~ end if
				! Lower right
				index1DFine = getIndex1D(2 * i, 2 * j - 1, level)
				rhs1D(index1DCoarse) = rhs1D(index1DCoarse) + res1D(index1DFine) / 4							
				!~ if (i == j) then
					!~ Print *, "index1DFine = ", index1DFine, "res1D = ", res1D(index1DFine) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!~ end if		
				! Upper right
				index1DFine = getIndex1D(2 * i, 2 * j, level)
				rhs1D(index1DCoarse) = rhs1D(index1DCoarse) + res1D(index1DFine) / 4					
				!~ if (i == j) then
					!~ Print *, "index1DFine = ", index1DFine, "res1D = ", res1D(index1DFine) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!~ end if
				if (i == j) then
					print *, i, j, index1DCoarse, rhs1D(index1DCoarse)	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				end if

			end do
		end do

	end subroutine restrict
	
	
	subroutine prolong(lhs1D, n1D, level)	

		integer, intent(in)	:: n1D
		real(8)				:: lhs1D(1 : n1D)
		integer, intent(in)	:: level
		integer				:: nTemp
		integer				:: i, j, index1DCoarse, index1DFine

		nTemp = n / 2 ** (level)
		
		do i = 1 , nTemp
			Do j = 1 , nTemp
				index1DCoarse = getIndex1D(i, j, level)
				! Upper left
				index1DFine = getIndex1D(2 * i - 1, 2 * j - 1, level - 1)
				lhs1D(index1DFine) = lhs1D(index1DFine) + lhs1D(index1DCoarse)
				! Upper right
				index1DFine = getIndex1D(2 * i - 1, 2 * j, level - 1)
				lhs1D(index1DFine) = lhs1D(index1DFine) + lhs1D(index1DCoarse)
				! Lower left
				index1DFine = getIndex1D(2 * i, 2 * j - 1, level - 1)
				lhs1D(index1DFine) = lhs1D(index1DFine) + lhs1D(index1DCoarse)		
				! Lower right
				index1DFine = getIndex1D(2 * i, 2 * j, level - 1)
				lhs1D(index1DFine) = lhs1D(index1DFine) + lhs1D(index1DCoarse)
			end do
		end do

	end subroutine prolong	
	
	
	subroutine vCycle(p, rhs, nRelax)
	
		integer,	intent(in)				:: nRelax
		real(8),	intent(in)				:: rhs(1 : n, 1 : n)
		real(8)							:: p(1 : n, 1 : n), pRes(1 : n, 1 : n)
		real(8)							:: nToReal
		integer 							:: nLevels
		integer							:: n1D
		integer							:: i, j, nTemp, thisIndex1D
		real(8), Dimension(:), Allocatable 	:: lhs1D, rhs1D, res1D
		
		nToReal = n * 1.0
		nLevels = int(log(nToReal) / log(2.0) + 0.1) - 2		
		
		print *, "No. of multi-grid levels: ", nLevels	
		
		n1D = 0
		
		! Allocate 1d arrays for all levels
		i = 0
		Do While (i <= nLevels)
			nTemp = n / (2 ** i)		
			n1D = n1D + (nTemp) ** 2	
			i = i + 1		
			print "(a10, i8, 10x, a10, i8)",  "nTemp = ", nTemp, "n1D = ", n1D		
		END Do
		
		Allocate(lhs1D(1 : n1D))
		Allocate(rhs1D(1 : n1D))
		Allocate(res1D(1 : n1D))		

		lhs1D = 0.
		rhs1D = 0.
		res1D = 0.
		
		! Copy from 2D array to 1D array
		do i = 1, n
			Do j = 1, n
				thisIndex1D = getIndex1D (i, j, 0)
				!print *, "1D index: ", thisIndex1D	!!!!!!!!!!!!!!!!!!!!!!!!	
				lhs1D(thisIndex1D) = p(i, j)
				!print *, "lhs: ", lhs1D(thisIndex1D)	!!!!!!!!!!!!!!!!!!!!!!!!
				rhs1D(thisIndex1D) = rhs(i, j)
				!print *, "rhs: ", rhs1D(thisIndex1D)	!!!!!!!!!!!!!!!!!!!!!!!!
			end do
		end do
		
		! Restrict down
		!Do i = 0, nLevels - 1
			i = 0	!!!!!!!!!!!!!!!!!!!!!!!
			call relaxGS(lhs1D, rhs1D, n1D, i, nRelax)
			call residue(lhs1D, rhs1d, res1d, n1D, i)
			call restrict(res1D, rhs1D, n1D, i)	
		!END Do
		
		!~ ! Lowest layer
		!~ call relaxGS(lhs1D, rhs1D, n1D, nLevels, nRelax)
		
		!~ ! Prolong up
		!~ Do i = nLevels, 1, -1
			!~ call prolong(lhs1D, n1D, i)
			!~ CALL relaxGS(lhs1D, rhs1D, n1D, i - 1, nRelax)	
		!~ END Do		
		
		! Copy from 1D array back to 2D array	
		do i = 1, n
			Do j = 1, n
				thisIndex1D = getIndex1D (i, j, 0)
				p(i, j) = lhs1D(thisIndex1D)
			end do
		end do

	end subroutine vCycle
	
	
	subroutine normaliseP(p)
	
		real(8)	:: p(1 : n, 1 : n)
		integer	:: i, j		
		real(8)	:: average = 0.
		
		do i = 1, n
			do j = 1, n
				average = average + p(i, j)
			end do
		end do
		
		average = average / n**2
		
		do i = 1, n
			do j = 1, n
				p(i, j) = p(i, j) - average
			end do
		end do		
	
	end subroutine normaliseP
	
	
	subroutine multiGridV(p, rhs)
			
		integer,	parameter	:: maxMGIt = 1
		integer,	parameter	:: nRelax = 5
		real(8), 	parameter	:: tol = 1.0e-5
		real(8),	intent(in)	:: rhs(1 : n, 1 : n)
		real(8)				:: p(1 : n, 1 : n)
		integer				:: i, j, nIt = 1 
		real(8)				:: resid = 1.0
		real(8)				:: pOld(1 : n, 1 : n)
		
		do while (nIt <= maxMGIt .and. resid > tol)
				nIt = nIt + 1
				pOld = p
				!print *, "Center pressure is: ", pOld(n / 2, n / 2)	!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call vCycle(p, rhs, nRelax)
				call normaliseP(p)
				resid = maxval(abs(pOld - p))
		end do
		
		print "(a17, 1x, i4, 5x, a9, 1x, f16.14)", "No. of iteraiton:", nIt, "Residual:", resid
	
	end subroutine multiGridV
	
	
end module moduleMultiGrid