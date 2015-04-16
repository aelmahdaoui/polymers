program polymers
	use polymers_init
	use polymer_dynamics
	use polymers_results	

	implicit none

	integer :: N, Ntheta, j, counter, Ndim, Ntests
	real(8) :: Temperature, PolWeight, beta
	real(8), dimension(:, :), allocatable :: Polymer, AverageDistance, AverageRadiusGyration, &
	RadiusGyration, TestVector, Distribution_vector, Weight_vector


	real ::  beg_cpu_time, end_cpu_time
	integer :: clck_counts_beg, clck_counts_end, clck_rate
	call system_clock ( clck_counts_beg, clck_rate )
	call cpu_time (beg_cpu_time)





	call ReadInputFile(N, Ndim, Ntheta, Temperature, Ntests)


	allocate ( Polymer(N, Ndim) )
	allocate ( AverageDistance(N,3) )
	allocate ( AverageRadiusGyration(N,3) )
	allocate ( RadiusGyration(N,3) )
	allocate ( TestVector(2, (N-2)) )
	allocate ( Distribution_vector(Ntests,(N-2)) )
	allocate ( Weight_vector(Ntests,(N-2)) )
		
	
	AverageDistance = 0d0
	AverageRadiusGyration = 0d0
	beta = 1d0/Temperature
	Distribution_vector = 0d0
	Weight_vector = 0d0

	call init_random_seed



do j=1,Ntests
	Polymer = 0d0
	Polymer(2,1) = 1d0
	RadiusGyration = 0d0
	counter = 3
	PolWeight = 1d0
	TestVector = 0d0

	
	if (mod(j,25) == 0) then
		print *, "Cycle ", j
	end if

	call AddBead(Polymer, PolWeight, counter, Ntheta, N, Ndim,  beta, AverageDistance, RadiusGyration, TestVector)

	AverageRadiusGyration(:,:) = AverageRadiusGyration(:,:) + RadiusGyration(:,:)
	!Distribution_vector(j,:) = Distribution_vector(j,:) + TestVector(1,:)
	!Weight_vector(j,:) = Weight_vector(j,:) + TestVector(2,:)
end do

	call WriteAverageDistanceToFile(AverageDistance,N)
	call WriteAverageRadiusGyrationToFile(AverageRadiusGyration,N)
	!call WriteTestVector(Distribution_vector, Ntests, (N-2))
	!call WriteTestVector2(Weight_vector, Ntests, (N-2))
	!print *, Ntests, N-2
	!call WritePolymer(Polymer,N)



	





	call system_clock ( clck_counts_end, clck_rate )
	call cpu_time (end_cpu_time)
	print *, "Clock time is: " , end_cpu_time - beg_cpu_time	
	print *, "Run time is: ", (clck_counts_end - clck_counts_beg) / real (clck_rate)

end program
