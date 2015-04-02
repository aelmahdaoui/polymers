module polymers_results
	implicit none
	private
	
	public  WriteAverageDistanceToFile, WriteAverageRadiusGyrationToFile, WriteTestVector

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

	end subroutine


	subroutine WriteAverageRadiusGyrationToFile(AverageRadiusGyration, N)
	integer, intent(in) :: N
	real(8), intent(in), dimension(N,2) :: AverageRadiusGyration
	integer :: i

	open (unit=4,file="radiusofgyration.dat",action="write")

	do i = 1,(N-2)
		if (i < 3) then
			write (4,"(I6, 4F35.5)")  i, 0d0	
		else
			write (4,"(I6, 4F35.5)")  i, AverageRadiusGyration(i,1)/AverageRadiusGyration(i,2)
		end if
	end do
	
	close(4)


	end subroutine

	subroutine WriteTestVector(vector, n, m)
		integer, intent(in) :: n,m
		integer :: i, j
		real(8), intent(in), dimension(n,m) :: vector


		open (unit=5,file="weightdistribution.data",action="write")

		do i = 1,n
				write(5,"(I2)",advance="no") i
				write (5,*)  ( vector(i,j), j=1,m )
		end do
	
		close(5)

	end subroutine




end module 
