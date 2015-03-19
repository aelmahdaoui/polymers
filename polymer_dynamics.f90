module polymer_dynamics
	implicit none
	private
	
	public init_random_seed, AddBead, CalcEnergy, UpdatePolymer, CalcWeight, PickThetaNumber

contains

	function CalcEnergy(Theta, Polymer, counter, N) result(energy)
		implicit none
		real(8), intent(in) :: Theta
		integer, intent(in) :: counter, N
		real(8), intent(in), dimension(N,2) :: Polymer
		integer :: i
		real(8), parameter :: e= 125d-3, s = 64d-2
		real(8), dimension(2) :: NewSegment
		real(8) :: energy, distance_squared

		NewSegment(1) = cos(Theta) + Polymer(counter-1,1)
		NewSegment(2) = sin(Theta) + Polymer(counter-1,2)
		energy = 0d0
		
		do i=1,(counter-1)
				distance_squared = (NewSegment(1) - Polymer(i,1))**2 + (NewSegment(2) - Polymer(i,2))**2
				energy = energy + 4*e*((s/distance_squared)**6-(s/distance_squared)**3)
		end do

	end function CalcEnergy
	

	subroutine AddBead(Polymer, PolWeight, counter, Ntheta, N, beta)
		integer, intent(in) :: counter, N, Ntheta
		integer :: THETA_NUMBER
		real(8), intent(inout) :: PolWeight
		real(8), intent(inout), dimension(N,2) :: Polymer
		real(8), dimension(Ntheta) :: weight_vector
		real(8) :: weight, theta, beta, RandomTheta
		real(8), parameter :: pi = 4*atan(1d0)

		
		
		weight_vector = 0d0 ! Check if necessary		
		weight = 0d0 ! Check if necessary

		call CalcWeight(Polymer, N, Ntheta, counter, beta, weight_vector, weight, RandomTheta)
		THETA_NUMBER = PickThetaNumber(Ntheta, weight_vector, weight)
		theta =  (THETA_NUMBER * 2 * pi) / Ntheta + RandomTheta
		weight = weight_vector(THETA_NUMBER)
		call UpdatePolymer(Polymer, N, counter, theta)

		PolWeight = PolWeight * weight

		if (counter < N) then
			call AddBead(Polymer, PolWeight, (counter + 1), Ntheta, N, beta)
		end if

	end subroutine

	subroutine CalcWeight(Polymer, N, Ntheta, counter, beta, weight_vector, sum_weight_vector, RandomTheta)
		real(8), parameter :: pi = 4*atan(1d0)
		real(8) :: xs, Theta
		integer :: i
		integer, intent(in) :: N, Ntheta, counter
		real(8), intent(in) :: beta
		real(8), intent(in), dimension(N,2) :: Polymer
		real(8), dimension(Ntheta) :: energy_vector
		real(8), intent(inout), dimension(Ntheta) :: weight_vector
		real(8), intent(inout) :: sum_weight_vector	
		real(8), intent(out) :: RandomTheta

		CALL RANDOM_NUMBER(xs)
		print *, "Random number", xs
		sum_weight_vector = 0d0

		RandomTheta = xs * 2 * pi
		do i = 1,Ntheta
			Theta = RandomTheta + i * 2 * pi/Ntheta
			energy_vector(i) = CalcEnergy(Theta, Polymer, counter, N)
			weight_vector(i) = exp(-1 * energy_vector(i) * beta)
			sum_weight_vector = sum_weight_vector + weight_vector(i)
		end do
		
	end subroutine


	function PickThetaNumber(Ntheta, weight_vector, total_weight) result(THETA_NUMBER)
		implicit none 
		integer, intent(in) :: Ntheta
		real(8), intent(in), dimension(Ntheta) :: weight_vector
		real(8), intent(in) :: total_weight
		real(8) :: weight, xs
		integer :: THETA_NUMBER, i
		
		THETA_NUMBER = 1
		
		CALL RANDOM_NUMBER(xs)
		print *, "Random number for picking an angle", xs
		weight = 0d0
		do i=1,Ntheta
			weight = weight + weight_vector(i) / total_weight
			if (weight < xs) then
				THETA_NUMBER = THETA_NUMBER + 1
			end if
		end do

	end function PickThetaNumber

	subroutine UpdatePolymer(Polymer, N, counter, theta)
		integer, intent(in) :: N, counter
		real(8), intent(inout), dimension(N,2) :: Polymer
		real(8), intent(in) :: theta


		Polymer(counter, 1) = cos(theta) + Polymer(counter-1, 1)
		Polymer(counter, 2) = sin(theta) + Polymer(counter-1, 2)
		
		
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
