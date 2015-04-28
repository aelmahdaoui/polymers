module polymers_results
	implicit none
	private
	
	public  WriteAverageDistanceToFile, WriteAverageRadiusGyrationToFile, WritePolymer

contains

	subroutine WriteAverageDistanceToFile(AverageDistance,N)
	integer, intent(in) :: N
	real(8), intent(in), dimension(N, 3) :: AverageDistance
	integer :: i

	open (unit=3,file="distances.dat",action="write")

	do i = 1,N
		if (i < 2) then
			write (3,"(I6, 4F35.5)")  i, AverageDistance(i,1), AverageDistance(i,2), AverageDistance(i,2), 0d0
		else if (i == 2) then
			write (3,"(I6, 4F35.5)")  i, 1d0 , AverageDistance(i,2), AverageDistance(i,2), 0d0	
		else if ( i > 2) then
			write (3,"(I6, 4F35.5)")  i, AverageDistance(i,1)/AverageDistance(i,2), AverageDistance(i,3), AverageDistance(i,2) &
			, log(AverageDistance(i,1)/AverageDistance(i,2)) 
		end if
	end do
	
	close(3)

	open (unit=3,file="logdistances.dat",action="write")

	do i = 1,(N-1)
		if (i < 2) then
			write (3,"(4F35.5)")  log10(i*1d0), log10(1d0) 
		else if ( i > 1) then
			write (3,"(4F35.5)")  log10(i+1d0), log10(AverageDistance((i+1),1))- log10(AverageDistance((i+1),2))
		end if
	end do
	
	close(3)

	end subroutine


	subroutine WriteAverageRadiusGyrationToFile(AverageRadiusGyration, N)
	integer, intent(in) :: N
	real(8), intent(in), dimension(N,3) :: AverageRadiusGyration
	integer :: i

	open (unit=4,file="radiusofgyration.dat",action="write")

	do i = 1,N
		if (i < 3) then
			write (4,"(I6, 4F35.5)")  i, 0d0	
		else
			write (4,"(I6, 4F35.5)")  i, AverageRadiusGyration(i,1)/AverageRadiusGyration(i,2), &
			 AverageRadiusGyration(i,2), AverageRadiusGyration(i,3) 
		end if
	end do
	
	close(4)


	end subroutine



	subroutine WritePolymer(vector, n)
		integer, intent(in) :: n
		integer :: i
		real(8), intent(in), dimension(n,2) :: vector


		open (unit=5,file="polymer.data",action="write")

		do i = 1,n
				write(5,*) vector(i,1), vector(i,2)
		end do
	
		close(5)

	end subroutine




end module 
