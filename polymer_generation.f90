module polymer_generation
	implicit none
	private
	
	public AddBead, init_random_seed

contains

	subroutine AddBead(Polymer, PolWeight, counter, Ntheta, N, Ndim, beta, AverageDistance, RadiusGyration, PERM, Ntests, test)
		integer, intent(in) :: counter, N, Ndim, Ntheta, PERM, Ntests, test
		integer :: THETA_NUMBER
		real(8), intent(inout) :: PolWeight
		real(8) :: NewWeight
		real(8), intent(inout), dimension(N,Ndim) :: Polymer
		real(8), dimension(Ntheta**(Ndim-1)) :: weight_vector
		real(8) :: weight, theta, phi, beta, RandomTheta, RandomPhi, R, LowLim, UpLim
		real(8), parameter :: pi = 4*atan(1d0), highalpha = 2.2, lowalpha = 1.1
		real(8), intent(inout), dimension(N,3) :: AverageDistance
		real(8), intent(inout), dimension(N,3) :: RadiusGyration




		
		
		weight_vector = 0d0 		
		weight = 0d0 
		theta = 0d0
		phi = 0d0


		call CalcWeight(Polymer, N, Ndim , Ntheta, counter, beta, weight_vector, weight, RandomTheta, RandomPhi) 			! Calculate the weights for adding a bead at different angles
		THETA_NUMBER = PickThetaNumber(Ndim, Ntheta, weight_vector, weight)								! Pick a random number based on the calculated weights
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!	Assign angles for the new bead		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (Ndim == 2) then
			theta =  (THETA_NUMBER * 2 * pi) / Ntheta + RandomTheta
		else if (Ndim == 3) then
			theta =  (mod(THETA_NUMBER, Ntheta) * 2 * pi) / Ntheta + RandomTheta	
			phi = 	((THETA_NUMBER - mod(THETA_NUMBER, Ntheta))/Ntheta  * pi) / Ntheta + RandomPhi
		end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!	Assign angles for the new bead		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		call UpdatePolymer(Polymer, N, Ndim, counter, theta, phi)

		if (Ndim == 2) then 
			PolWeight = PolWeight * weight/(0.75*Ntheta)
		else if (Ndim == 3 ) then
			PolWeight = PolWeight*weight/(4*Ntheta)
		end if

		
		AverageDistance(counter,1) = AverageDistance(counter,1) + (Polymer(counter,1)**2 + Polymer(counter,2)**2) * PolWeight			! Record weighted end-to-end distance
		AverageDistance(counter,2) = AverageDistance(counter,2) + PolWeight									! Record the sum of the weights for normalization
		AverageDistance(counter,3) = AverageDistance(counter,3) + 1										! Record number of polymers with this length
		if (Ndim == 3) then
			AverageDistance(counter,1) = AverageDistance(counter,1) + Polymer(counter,3)**2 * PolWeight					! Record weighted end-to-end distance
		end if
		RadiusGyration(counter,1) = RadiusGyration(counter,1) + CalcRadiusGyrationSquared(Polymer,counter, N, Ndim) * PolWeight
		RadiusGyration(counter,2) = AverageDistance(counter,2)
		RadiusGyration(counter,3) = AverageDistance(counter,3)






if (PERM == 1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!	Pruning & enriching the polymer		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		

		if (counter == 3) then
			call AddBead(Polymer, PolWeight, (counter + 1), Ntheta, N, Ndim, beta, AverageDistance, RadiusGyration, PERM, Ntests, test)
		else if (counter < N) then
			if (PolWeight > 0 ) then
				if ( AverageDistance(counter,3) < 1000) then
					call AddBead(Polymer, PolWeight, (counter + 1), Ntheta, N, Ndim, beta, AverageDistance, RadiusGyration, PERM, Ntests, test)
				else
					UpLim = highalpha * AverageDistance(counter,2)*AverageDistance(3,3)/(AverageDistance(counter,3)* AverageDistance(3,2))
					LowLim = lowalpha * AverageDistance(counter,2)*AverageDistance(3,3)/(AverageDistance(counter,3)* AverageDistance(3,2))
					if (PolWeight > UpLim) then
						call random_number(R)
						if ( R < 50d0/counter**0.7 * test/Ntests ) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!The random number is used to control the population and prevent unlimited doubling or fast attrition. The last factor is to ensure that polymers aren't strongly correlated and doubling takes place across all simulations, the factor before that is so that most doubling happens early on in the chain.
!Good values found are:
!For T = 1.2, Ntheta = 9, 2D, N= 250,  70d0/counter**0.8 * test/Ntests
!For T = 1.5, Ntheta = 9, 2D, N= 250,  40d0/counter**0.8 * test/Ntests
!For T = 2.5, Ntheta = 5, 3D, N= 250,  50d0/counter**0.7 * test/Ntests
!For T = +2, Ntheta = 10, 2D, N= 250,  no number is necessary, unlimited doubling is ok

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

							NewWeight = 0.5 * PolWeight
							call AddBead(Polymer, NewWeight, (counter + 1), Ntheta, N, Ndim, beta, AverageDistance, RadiusGyration, PERM, Ntests, test)
							call AddBead(Polymer, NewWeight, (counter + 1), Ntheta, N, Ndim, beta, AverageDistance, RadiusGyration, PERM, Ntests, test)
						else 
							call AddBead(Polymer, PolWeight, (counter + 1), Ntheta, N, Ndim, beta, AverageDistance, RadiusGyration, PERM, Ntests, test)
						end if
					else if (PolWeight < LowLim) then
						CALL RANDOM_NUMBER(R)
						if (R < 0.5) then
							NewWeight = 2*PolWeight
							call AddBead(Polymer, NewWeight, (counter + 1), Ntheta, N, Ndim, beta, AverageDistance, RadiusGyration, PERM, Ntests, test)
						end if
					else	
						call AddBead(Polymer, PolWeight, (counter + 1), Ntheta, N, Ndim, beta, AverageDistance, RadiusGyration, PERM, Ntests, test)
					end if
				end if
			end if
		end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!	Pruning & enriching the polymer		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








else if (PERM == 0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!	No Pruning & enriching the Polymer	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (counter < N) then
			if (PolWeight > 0 ) then
				call AddBead(Polymer, PolWeight, (counter + 1), Ntheta, N, Ndim, beta, AverageDistance, RadiusGyration, PERM, Ntests, test)
			end if
		end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!	No Pruning & enriching the Polymer	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if
	end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!	Calculating the energy for adding a new segment to a polymer	!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function CalcEnergy(Theta, Phi, Polymer, counter, N, Ndim) result(energy)
		implicit none
		real(8), intent(in) :: Theta, Phi
		integer, intent(in) :: counter, N, Ndim
		real(8), intent(in), dimension(N,Ndim) :: Polymer
		integer :: i
		real(8), parameter :: e= 0.25, s = 0.64
		real(8), dimension(Ndim) :: NewSegment
		real(8) :: energy, distance_squared

		energy = 0d0
		distance_squared = 0d0

		if (Ndim == 2) then
			NewSegment(1) = cos(Theta) + Polymer(counter-1,1)
			NewSegment(2) = sin(Theta) + Polymer(counter-1,2)
		else if (Ndim == 3) then
			NewSegment(1) = cos(Theta) * sin(Phi) + Polymer(counter-1,1)
			NewSegment(2) = sin(Theta) * sin(Phi) + Polymer(counter-1,2)
			NewSegment(3) = cos(Phi) + Polymer(counter-1,3)
		end if 

		
		do i=1,(counter-1)
				if (Ndim == 2) then
					distance_squared = (NewSegment(1) - Polymer(i,1))**2 + (NewSegment(2) - Polymer(i,2))**2
				else if (Ndim == 3) then
					distance_squared = (NewSegment(1) - Polymer(i,1))**2 + (NewSegment(2) - Polymer(i,2))**2 + (NewSegment(3) - Polymer(i,3))**2
				end if 
				energy = energy + 4*e*((s/distance_squared)**6-(s/distance_squared)**3)
		end do

	end function CalcEnergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!	Calculating the energy for adding a new segment to a polymer	!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculating the Boltzmann factors for adding a new segment in different directions of a polymer !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine CalcWeight(Polymer, N, Ndim, Ntheta, counter, beta, weight_vector, sum_weight_vector, RandomTheta, RandomPhi)
		real(8), parameter :: pi = 4*atan(1d0)
		real(8) :: xs, ys, Theta, Phi
		integer :: i, j
		integer, intent(in) :: N, Ndim, Ntheta, counter
		real(8), intent(in) :: beta
		real(8), intent(in), dimension(N,2) :: Polymer
		real(8), dimension(Ntheta**(Ndim -1)) :: energy_vector
		real(8), intent(inout), dimension(Ntheta**(Ndim-1)) :: weight_vector
		real(8), intent(inout) :: sum_weight_vector	
		real(8), intent(out) :: RandomTheta, RandomPhi

		CALL RANDOM_NUMBER(xs)
		CALL RANDOM_NUMBER(ys)
		sum_weight_vector = 0d0
		weight_vector = 0d0
		Theta = 0d0
		Phi = 0d0

		RandomTheta = xs * 2 * pi
		RandomPhi = ys * 2 * pi

		if (Ndim == 2) then
			do i = 1,Ntheta
				Theta = RandomTheta + i * 2 * pi/Ntheta
				energy_vector(i) = CalcEnergy(Theta, Phi, Polymer, counter, N, Ndim)
				weight_vector(i) = exp(-1 * energy_vector(i) * beta)
				sum_weight_vector = sum_weight_vector + weight_vector(i)
				!print *, energy_vector(i), weight_vector(i), sum_weight_vector, counter
			end do
		else if (Ndim == 3) then 
			do j = 1,Ntheta
				Phi = RandomPhi + j * pi/Ntheta
				do i = 1, Ntheta	
					Theta = RandomTheta + i * 2 * pi/Ntheta
					energy_vector(i+(j-1)*Ntheta) = CalcEnergy(Theta, Phi, Polymer, counter, N, Ndim)
					weight_vector(i+(j-1)*Ntheta) = exp(-1 * energy_vector(i+(j-1)*Ntheta) * beta)
					sum_weight_vector = sum_weight_vector + weight_vector(i+(j-1)*Ntheta)
				end do
			end do
		end if 
		
	end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculating the Boltzmann factors for adding a new segment in different directions of a polymer !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	function PickThetaNumber(Ndim, Ntheta, weight_vector, total_weight) result(THETA_NUMBER)
		implicit none 
		integer, intent(in) :: Ntheta, Ndim
		real(8), intent(in), dimension(Ntheta**(Ndim-1)) :: weight_vector
		real(8), intent(in) :: total_weight
		real(8) :: weight, xs
		integer :: THETA_NUMBER, i
		
		THETA_NUMBER = 1
		
		CALL RANDOM_NUMBER(xs)
		weight = 0d0
		do i=1,Ntheta**(Ndim-1)
			weight = weight + weight_vector(i) / total_weight
			if (weight < xs) then
				THETA_NUMBER = THETA_NUMBER + 1
			end if
		end do

	end function PickThetaNumber



	subroutine UpdatePolymer(Polymer, N, Ndim, counter, theta, phi)
		integer, intent(in) :: N, Ndim, counter
		real(8), intent(inout), dimension(N,Ndim) :: Polymer
		real(8), intent(in) :: theta, phi

		if (Ndim == 2) then
			Polymer(counter, 1) = cos(theta) + Polymer(counter-1, 1)
			Polymer(counter, 2) = sin(theta) + Polymer(counter-1, 2)
		else if (Ndim == 3) then 
			Polymer(counter, 1) = cos(theta) * sin(phi) + Polymer(counter-1,1)
			Polymer(counter, 2) = sin(theta) * sin(phi) + Polymer(counter-1,2)
			Polymer(counter, 3) = cos(phi) + Polymer(counter-1,3)
		end if
		
	end subroutine

	function CalcRadiusGyrationSquared(Polymer,counter, N, Ndim) result(radius)
		integer, intent(in) :: counter, N, Ndim
		integer :: i, j
		real(8), dimension(Ndim) :: DifferenceSquared, center_mass
		real(8) :: radius
		real(8), intent(in), dimension(N,Ndim) :: Polymer
		
		DifferenceSquared = 0d0
		radius = 0d0
		center_mass = CalcCenterMass(Polymer, counter, N,Ndim)

		do i=1,counter
			do j=1,Ndim
				DifferenceSquared(j) = DifferenceSquared(j) + (Polymer(i,j)- center_mass(j))**2
			end do
		end do

		do j=1,Ndim
			radius = radius + DifferenceSquared(j)/counter 
		end do
	end function



	function CalcCenterMass (Polymer,counter, N,Ndim) result(CenterMass)
		integer, intent(in) :: N, Ndim, counter
		integer :: i, j
		real(8), dimension(Ndim) :: CenterMass
		real(8), intent(in), dimension(N,Ndim) :: Polymer

		CenterMass = 0d0
		
		do i=1,counter
			do j=1,Ndim
				CenterMass(j) = CenterMass(j) + Polymer(i,j)
			end do
		end do

		do j=1,Ndim
			CenterMass(j) = CenterMass(j)/counter
		end do

	end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!			Generating random numbers			!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!			Copied from ICCP coding-notes			!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine init_random_seed()
		implicit none
		integer, allocatable :: seed(:)
		integer :: i, n, un, istat, dt(8), pid, t(2), s
		integer(8) :: count, tms

		call random_seed(size = n)
		allocate(seed(n))
		open(newunit=un, file="/dev/urandom", access="stream",&
		form="unformatted", action="read", status="old", &
		iostat=istat)
		if (istat == 0) then
			read(un) seed
			close(un)
		else
			call system_clock(count)
			if (count /= 0) then
				t = transfer(count, t)
			else
				call date_and_time(values=dt)
				tms = (dt(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
				+ dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
				+ dt(3) * 24 * 60 * 60 * 60 * 1000 &
				+ dt(5) * 60 * 60 * 1000 &
				+ dt(6) * 60 * 1000 + dt(7) * 1000 &
				+ dt(8)
				t = transfer(tms, t)
			end if
			s = ieor(t(1), t(2))
			pid = getpid() + 1099279 ! Add a prime
			s = ieor(s, pid)
			if (n >= 3) then
				seed(1) = t(1) + 36269
				seed(2) = t(2) + 72551
				seed(3) = pid
				if (n > 3) then
					seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
				end if
			else
				seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
			end if
		end if
		call random_seed(put=seed)
	end subroutine init_random_seed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!			Generating random numbers			!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!			Copied from ICCP coding-notes			!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module 
