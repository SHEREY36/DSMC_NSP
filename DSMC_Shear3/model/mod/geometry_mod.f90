  module geometry
    implicit none
	
    ! General Geometry
    double precision :: xlength, ylength, zlength
    double precision, dimension(3) :: length
	
    ! Boundary Flags
    ! Periodic
    logical :: periodic_x, periodic_y, periodic_z
    logical, dimension(3) :: periodic
    ! Uniform Shear flow
    logical :: shear_x, shear_y, shear_z
    logical, dimension(3) :: shear_bc
		
    ! Wall velocities (not defined for periodic)
    double precision, dimension(3,3) :: u_bc
    contains
	
  end module geometry