module polymer_dynamics
	implicit none
	private
	
	public init_random_seed, AddBead, CalcEnergy, UpdatePolymer, CalcWeight, PickThetaNumber

contains

	function CalcEnergy(Theta, Phi, Polymer, counter, N, Ndim) result(energy)
		implicit none
		real(8), intent(in) :: Theta, Phi
		integer, intent(in) :: counter, N, Ndim
		real(8), intent(in), dimension(N,Ndim) :: Polymer
		integer :: i
		real(8), parameter :: e= 125d-3, s = 64d-2
		real(8), dimension(Ndim) :: NewSegment
		real(8) :: energy, distance_squared

		if (Ndim == 2) then
			NewSegment(1) = cos(Theta) + Polymer(counter-1,1)
			NewSegment(2) = sin(Theta) + Polymer(counter-1,2)
		else if (Ndim == 3) then
			NewSegment(1) = cos(Theta) * sin(Phi) + Polymer(counter-1,1)
			NewSegment(2) = sin(Theta) * sin(Phi) + Polymer(counter-1,2)
			NewSegment(3) = cos(Phi) + Polymer(counter-1,3)
		end if 
		energy = 0d0
		
		distance_squared = 0d0
		
		do i=1,(counter-1)
				if (Ndim == 2) then
					distance_squared = (NewSegment(1) - Polymer(i,1))**2 + (NewSegment(2) - Polymer(i,2))**2
				else if (Ndim == 3) then
					distance_squared = (NewSegment(1) - Polymer(i,1))**2 + (NewSegment(2) - Polymer(i,2))**2 + (NewSegment(3) - Polymer(i,3))**2
				end if 
				energy = energy + 4*e*((s/distance_squared)**6-(s/distance_squared)**3)
		end do

	end function CalcEnergy
	

	subroutine AddBead(Polymer, PolWeight, counter, Ntheta, N, Ndim, beta)
		integer, intent(in) :: counter, N, Ndim, Ntheta
		integer :: THETA_NUMBER
		real(8), intent(inout) :: PolWeight
		real(8), intent(inout), dimension(N,2) :: Polymer
		real(8), dimension(Ntheta**(Ndim-1)) :: weight_vector
		real(8) :: weight, theta, phi, beta, RandomTheta, RandomPhi
		real(8), parameter :: pi = 4*atan(1d0)

		
		
		weight_vector = 0d0 ! Check if necessary		
		weight = 0d0 ! Check if necessary

		call CalcWeight(Polymer, N, Ndim , Ntheta, counter, beta, weight_vector, weight, RandomTheta, RandomPhi)
		!print *, weight_vector
		THETA_NUMBER = PickThetaNumber(Ndim, Ntheta, weight_vector, weight)
		if (Ndim == 2) then
			theta =  (THETA_NUMBER * 2 * pi) / Ntheta + RandomTheta
		else if (Ndim == 3) then
			theta =  (mod(THETA_NUMBER, Ntheta) * 2 * pi) / Ntheta + RandomTheta	
			phi = 	((THETA_NUMBER - mod(THETA_NUMBER, Ntheta))/Ntheta  * pi) / Ntheta + RandomPhi
		end if

		weight = weight_vector(THETA_NUMBER)

		call UpdatePolymer(Polymer, N, Ndim, counter, theta, phi)

		PolWeight = PolWeight * weight

		if (counter < N) then
			call AddBead(Polymer, PolWeight, (counter + 1), Ntheta, N, Ndim, beta)
		end if

	end subroutine

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
		print *, "Random number", xs
		sum_weight_vector = 0d0
		weight_vector = 0d0

		RandomTheta = xs * 2 * pi
		RandomPhi = ys * 2 * pi

		if (Ndim == 2) then
			do i = 1,Ntheta
				Theta = RandomTheta + i * 2 * pi/Ntheta
				energy_vector(i) = CalcEnergy(Theta, Phi, Polymer, counter, N, Ndim)
				weight_vector(i) = exp(-1 * energy_vector(i) * beta)
				sum_weight_vector = sum_weight_vector + weight_vector(i)
				!print *, "some"
			end do
		else if (Ndim == 3) then 
			do j = 1,Ntheta
				Phi = RandomPhi + j * pi/Ntheta
				do i = 1, Ntheta	
					Theta = RandomTheta + i * 2 * pi/Ntheta
					energy_vector(i+(j-1)*Ntheta) = CalcEnergy(Theta, Phi, Polymer, counter, N, Ndim)
					weight_vector(i+(j-1)*Ntheta) = exp(-1 * energy_vector(i+(j-1)*Ntheta) * beta)
					sum_weight_vector = sum_weight_vector + weight_vector(i+(j-1)*Ntheta)
					!print *, "Total weight is:", sum_weight_vector, "weight is ", weight_vector(i+(j-1)*Ntheta)
				end do
			end do
		end if 
		
		!print *, "Dimension = ", Ndim, " Total weight = ", sum_weight_vector
		
	end subroutine


	function PickThetaNumber(Ndim, Ntheta, weight_vector, total_weight) result(THETA_NUMBER)
		implicit none 
		integer, intent(in) :: Ntheta, Ndim
		real(8), intent(in), dimension(Ntheta**(Ndim-1)) :: weight_vector
		real(8), intent(in) :: total_weight
		real(8) :: weight, xs
		integer :: THETA_NUMBER, i
		
		THETA_NUMBER = 1
		
		CALL RANDOM_NUMBER(xs)
		print *, "Random number for picking an angle", xs
		weight = 0d0
		do i=1,Ntheta**(Ndim-1)
			!print *, "Weight vector for ", i , " is : " , weight_vector(i)
			weight = weight + weight_vector(i) / total_weight
			if (weight < xs) then
				THETA_NUMBER = THETA_NUMBER + 1
				!print *, "Theta number = ", THETA_NUMBER
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



	! copied from ICCP coding-notes
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

end module 
