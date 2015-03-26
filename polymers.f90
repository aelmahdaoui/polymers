		

program polymers
	use polymers_init
	use polymer_dynamics	

	implicit none

	integer :: N, Ntheta, i, counter, Ndim
	real(8) :: Temperature, PolWeight, beta
	real(8), dimension(:, :), allocatable :: Polymer


	real ::  beg_cpu_time, end_cpu_time
	integer :: clck_counts_beg, clck_counts_end, clck_rate
	call system_clock ( clck_counts_beg, clck_rate )
	call cpu_time (beg_cpu_time)



	call ReadInputFile(N, Ndim, Ntheta, Temperature)
	print *, "N=", N, "Dimension =", Ndim,  "Ntheta=", Ntheta, "T=", Temperature

	allocate ( Polymer(N, Ndim) )

	Polymer = 0d0
	Polymer(2,1) = 1d0


	counter = 3
	PolWeight = 1d0
	beta = 1d0/Temperature

	call AddBead(Polymer, PolWeight, counter, Ntheta, N, Ndim,  beta)

	do i =1,N

		print *, Polymer(i,:)
	end do

	open (unit=2,file="beads_positions.dat",action="write")
	do i=1,N

		if (Ndim == 2) then
			write (2,"(4F18.10)")  Polymer(i,1), Polymer(i,2)
		else if (Ndim == 3) then
			write (2,"(4F18.10)")  Polymer(i,1), Polymer(i,2), Polymer(i,3)
		end if

	end do 
	
	


	call system_clock ( clck_counts_end, clck_rate )
	call cpu_time (end_cpu_time)
	print *, "Clock time is: " , end_cpu_time - beg_cpu_time	
	print *, "Run time is: ", (clck_counts_end - clck_counts_beg) / real (clck_rate)




end program
