module polymers_init
	implicit none
	private
	
	public  ReadInputFile!

contains

	subroutine ReadInputFile(N, Ndim, Ntheta, T)
	integer, intent(out) :: N, Ndim, Ntheta
	real(8), intent(out) :: T	

	open(1, file='PolymerInput.txt', status='old', action='read')            ! opens the first file
  	read(1,*), N
	read(1,*), Ntheta
	read(1,*), T 
	read(1,*), Ndim

	close(1)

	end subroutine




end module 
