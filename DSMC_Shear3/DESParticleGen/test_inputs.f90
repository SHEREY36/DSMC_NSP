
	program test_inputs
	
	implicit none
	
	character(len=255) :: ifilename, pfilename
	
	double precision :: xlength, ylength, zlength
	integer :: nx_dsmc, ny_dsmc, nz_dsmc
	integer :: pip
	double precision :: RealoSim
	double precision :: dia, rho
	double precision :: dt, tend
	logical :: periodic_x, periodic_y, periodic_z
	double precision, dimension(3,3) :: u_bc
	
	integer :: k
	
	character(len=255) :: filename_input
	double precision :: solfrac
	
	filename_input	= 'init_vars_input.dat'
	
	OPEN(UNIT=100, FILE=filename_input)
	
	read(100,*) xlength, ylength, zlength
	read(100,*) nx_dsmc, ny_dsmc, nz_dsmc
	read(100,*) pip
	read(100,*) RealoSim
	read(100,*) dia
	read(100,*) rho
	read(100,*) dt, tend
	read(100,*) periodic_x, periodic_y, periodic_z
	u_bc = 0.d0
	read(100,*) u_bc(1,:)
	read(100,*) u_bc(2,:)
	read(100,*) u_bc(3,:)
	write(*,*) xlength, ylength
	write(*,*) tend
	write(*,*) periodic_z
	write(*,*) u_bc
	return


	end program test_inputs
