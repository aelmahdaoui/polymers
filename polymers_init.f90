module polymers_init
	implicit none
	private
	
	public  ReadInputFile!, InitFirst2Beads

contains

	subroutine ReadInputFile(N, Ntheta, T)
	integer, intent(out) :: N, Ntheta
	real(8), intent(out) :: T	

	OPEN(1, FILE='PolymerInput.txt', STATUS='OLD', ACTION='READ')            ! opens the first file
  	READ(1,*), N
	read(1,*), Ntheta
	read(1,*), T 

	close(1)

	end subroutine




end module 
