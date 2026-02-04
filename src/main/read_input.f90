  subroutine read_input	
    use particle
    use run_param
    use dsmc
    use geometry
    use output
	
    implicit none
	
    integer :: i
    character(len=255) :: filename_input
    integer :: lunit = 10
    double precision :: solfrac
	
    filename_input	= 'system_input.dat'
	
    OPEN(unit=lunit, FILE=filename_input, status='old')

    read(10,*) xlength, ylength, zlength
    read(10,*) nx_dsmc, ny_dsmc, nz_dsmc
    read(10,*) nx_meas, ny_meas, nz_meas
    read(10,*) pip
    read(10,*) F_N ! simulated to real part ratio
    read(10,*) dia, AR
    read(10,*) rho
    read(10,*) alpha_pp
    read(10,*) dt, tend
    read(10,*) periodic_x, periodic_y, periodic_z
    u_bc = 0.d0
    read(10,*) (u_bc(1,i),i=1,3) ! x-normal
    read(10,*) (u_bc(2,i),i=1,3) ! y-normal
    read(10,*) (u_bc(3,i),i=1,3) ! z-normal
    read(10,*) rwalk

    FIRST_CYCLE = .True.
    close(10)
    return
  end subroutine read_input