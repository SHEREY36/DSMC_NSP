!
!******************************************************************
!
!  Particle Generator
!  Generates Monodisperse Spherocylinder DSMC particles in random positions 
!  with a uniform speed dictated by an inital kinetic energy for 
!  the system 
!
!  Author: Muhammed Gbolasere
!******************************************************************

   Program DES_Particle_Generator
   implicit none

   ! Seed values
   ! Length of seed, Reference - time
   integer :: sz, time(8)
   ! Seed for random number
   integer, dimension(:), allocatable :: seed

   ! Values to be read (some are outputted)
   ! Solids Fraction
   double precision :: EPS
   ! Particle Radii
   double precision :: Diameter, radius, sigma_c, bmax
   ! Spherocylinder Aspect Ratio
   double precision :: AR, lcycl, mI, omI
   ! Particle Densities
   double precision :: DENSITY
   ! Particle Volume
   double precision :: pvol, n, mass, omass
   ! System Dimensions
   double precision :: xl, yl, zl
   ! Granular Temperatures
   double precision :: KTm, KTr, lam
   ! Time of simulation/time step
   double precision :: tend, dt
   ! Boundary Conditions
   logical :: periodic(3)
   double precision, dimension(3,3) :: u_bounds
   ! Measuring Grid
   integer, dimension(3) :: NDIM_DSMC, NDIM_MEAS
	! Restitution Coefficient
   double precision :: alpha
   ! Random walk flag
   logical :: rwalk

   ! Particle-Related Properties
   ! Number of Particles
   integer :: np
   ! Particle Velocity, STD of V_P
   double precision :: vs, pV
   double precision :: oSTDs, vmaxs
   double precision :: vmax, STD, oSTD
   double precision :: ktmhot
   double precision :: vmaxhot, STDhot, oSTDhot
   ! Temp Particle Velocities
   real*8 :: u1, u2, u3
   ! Temp Particle Positions
   real*8 :: xp, yp, zp
   ! Bulk Velocity
   double precision, dimension(3) :: BULK
   double precision :: beta
   DOUBLE PRECISION, PARAMETER :: PI = 4.D0*ATAN(1.D0)

   ! RDF
   double precision, parameter :: sol_frac_pack = 0.643D0
   double precision, parameter :: C1 = 2.5D0, C2 = 4.5904D0
   double precision, parameter :: C3 = 4.515439D0, C_POW = 0.67802D0
   double precision :: EPS_PACK=0.643d0, EPS_Ratio
   double precision :: numer, denom, chi_ma

   ! Values to be outputted
   ! Particle Position/Velocity Array
   real*8, dimension(:,:), allocatable :: r, u, omega
   ! Adjusted Particle Params for systems
   logical, dimension(:), allocatable :: skip
   
   ! Tolerances
   DOUBLE PRECISION :: SMALL_NUM

   ! Dummy Vars
   integer :: i, k
   double precision :: RR
  
   open(unit=10, file="param.in", status='old')
   open(unit=20, file="particle_input.dat", status='replace')
   open(unit=30, file="system_input.dat", status='replace')

   u_bounds = 0.d0
   read (10,*) EPS
   read (10,*) Diameter
   read (10,*) AR
   small_num = diameter/4.d0
   read (10,*) Density
   read (10,*) xl
   read (10,*) yl
   read (10,*) zl
   read (10,*) (NDIM_DSMC(k),k=1,3)
   read (10,*) (NDIM_MEAS(k),k=1,3)
   read (10,*) kTm
   read (10,*) kTr
   read (10,*) alpha
   read (10,*) np
   read (10,*) tend
   read (10,*) (periodic(k),k=1,3)
   read (10,*) u_bounds(1,2), u_bounds(1,3)
   read (10,*) u_bounds(2,1), u_bounds(2,3)
   read (10,*) u_bounds(3,1), u_bounds(3,2)
   read (10,*) rwalk
   
   ! Time end
   ! Periodic
   ! u_xy, v_xy
   ! v_yz, w_yz
   ! u_xz, w_xz

3000	format( d10.4 )
3001	format( I10 )
   
   ! Random seed for new particle files
   call random_seed( size = sz )
   allocate(seed(sz))
   call date_and_time(values=time)
   seed = (time(8)*2 + time(7)*3 + time(6)*5)*(/(I, I=1, sz)/)
   call random_seed(PUT=seed)

   ! Allocate Arrays for Output
   allocate( u(np,3) ) 
   allocate( r(np,3) )
   allocate(omega(np,3))
   
   ! Intialize Arrays
   u(:,:) = 0.d0 
   r(:,:) = 0.d0
   omega(:,:) = 0.d0

   ! Spherocylindrical properties
   lcycl = (AR-1.0d0)*diameter;
   radius = diameter/2.d0
   mI = (1.d0/48.d0*density*diameter**2.d0*lcycl**3.d0 + 3.d0/64.d0*density*diameter**4.d0*lcycl &
			+ 1.d0/60.d0*density*diameter**5.d0 + 1.d0/24.d0*density*diameter**3.d0*lcycl**2.d0)*pi; ! same as Iy

   pvol = ((diameter**3.d0)/6.d0 + lcycl*(radius**2.d0))*pi;
   mass = density*pvol


   sigma_c = (0.32d0*(AR**2.d0) + 0.694d0*AR - 0.014)*pi*(diameter**2)
   bmax = sqrt(sigma_c / pi)

   small_num = bmax/4.d0

   ! DSMC particles
   beta = EPS*XL*YL*ZL/PVOL	! Number of real particles
   beta = dble(beta)/dble(np) ! ! F_N

   ! Particle Positions
   CALL RANDOM_NUMBER(r)
   ! Buffer
   r = r*(1.d0 - 2.d0*small_num) + small_num
   r(:,1) = r(:,1)*xl
   r(:,2) = r(:,2)*yl
   r(:,3) = r(:,3)*zl

   ! Particle Velocities
   do i = 1,Np
      VMAX = 6.0D0
      k = 1
      do while(k.LE.3)
	 CALL RANDOM_NUMBER(RR)
         vs = (2.D0*RR-1.D0)*vmax
	 pV = exp(-(vs**2.D0)*0.5D0)
	 CALL RANDOM_NUMBER(RR)
	 IF(pV.GT.RR) THEN
	    u(i,k) = vs*SQRT(kTm / mass)
	    k = k + 1
	 end if
      end do
   end do

   ! Maxwell-Boltzmann Sampling for Angular Velocities
   do i = 1, Np
      VMAX = 6.0D0
      k = 1
      omega(i,1) = 0.d0  ! Zero out twist (x-axis) component
      do while(k.LE.2)
          call random_number(RR); 
          vs = (2.D0*RR-1.D0)*vmax
          pV = exp(-(vs**2.D0)*0.5)
          call random_number(RR)
          if(pV.GT.RR) then
             omega(i,k+1) = vs*SQRT(kTr / mI) 
             k = k + 1
          end if
      end do
   end do

   ! Shift bulk velocity to zero
   BULK(:) = 0.d0
   do i = 1,Np
   	BULK(:) = u(i,:) + BULK(:)
   end do
   BULK = BULK/dble(Np)
   
   do i = 1,Np
   	u(i,:) = u(i,:) - BULK
	write(20,12) r(i,1), r(i,2), r(i,3), u(i,1), u(i,2), u(i,3), omega(i,1), omega(i,2), omega(i,3)
   end do


   !n = EPS / pvol
   !lam = 1.d0 / (sqrt(2.d0) * n * sigma_c)
   !dt = lam/(sqrt(2.d0*kTm)*10.d0)
   dt = 1.0d-3
   
   ! System_input
   write(30,'(3(d12.6,2x))') xl, yl, zl
   write(30,'(3(I10,2x))') NDIM_DSMC
   write(30,'(3(I10,2x))') NDIM_MEAS
   write(30,'(I10)') np
   write(30,'(d12.6)') beta
   write(30,'(2(d12.6,2x))') diameter, AR
   write(30,'(d12.6)') density
   write(30,'(d12.6)') alpha
   write(30,'(2(d12.6,2x))') dt, tend
   write(30,*) periodic
   write(30,'(3(d10.4,2x))') u_bounds(1,:)
   write(30,'(3(d10.4,2x))') u_bounds(2,:)
   write(30,'(3(d10.4,2x))') u_bounds(3,:)
   write(30,*) rwalk

 12     FORMAT (9(d10.4,2x))
   end
