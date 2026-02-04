  subroutine initialize
    use run_param, only: first_cycle

    implicit none

    integer :: sz, time(8), I
    integer, dimension(:), allocatable :: seed
   
    ! Pick random seed
    call random_seed( size = sz )
    allocate(seed(sz))
    call date_and_time(values=time)
    seed = (time(8)*2 + time(7)*3 + time(6)*5)*(/(I, I=1, sz)/)
    call random_seed( PUT = seed )

    IF(FIRST_CYCLE) THEN
   	write(*,*) 'Setting up Grid'
	call init_geometry()
	write(*,*) 'Allocating Particle Arrays'
	call init_part()
	write(*,*) 'Initalizing DSMC Parameters'
	call init_dsmc()
	FIRST_CYCLE = .FALSE.
    ELSE
	call init_part()
    END IF
  end subroutine initialize

  ! Geometry
  !-----------------------------------------------------------------	
	
  SUBROUTINE INIT_GEOMETRY
    use geometry

    implicit none
    integer :: i,j

    PERIODIC = (/PERIODIC_X, PERIODIC_Y, PERIODIC_Z/)

    do i = 1,3
       j = mod(i+1,3); 
       if(j.eq.0) j = 3
       IF(u_bc(i,i).NE.0.d0.OR.u_bc(i,j).NE.0.d0) then
	  shear_bc(i) = .true.
       ELSE
	  shear_bc(i) = .false.
       END IF
    end do

    shear_x = shear_bc(1); 
    shear_y = shear_bc(2); 
    shear_z = shear_bc(3)
    LENGTH = (/XLENGTH, YLENGTH, ZLENGTH/)

  end subroutine INIT_GEOMETRY	
	
  ! Particles
  !-----------------------------------------------------------------	
	
  SUBROUTINE INIT_PART
    use particle
    use constants
    use dsmc
    use geometry
    use collision_tables_mod
    use gmm_cond_mod

    implicit none
	
    integer :: I, K, J
    integer :: lunit = 20
    double precision, dimension(3) :: pos
    double precision :: K_GRANULAR, ROT_GRANULAR
    character(255) :: lFILENAME
    character(len=256) :: gmm_file
	
    if(.not.allocated(des_pos_new)) allocate(DES_POS_NEW(PIP,3))
    if(.not.allocated(des_vel_new)) allocate(DES_VEL_NEW(PIP,3))
    if(.not.allocated(des_omega_new)) allocate(DES_OMEGA_NEW(PIP,3))
    if(.not.allocated(Er)) allocate(Er(PIP))

    des_pos_new(:,:) = 0.d0
    des_vel_new(:,:) = 0.d0
    des_omega_new(:,:) = 0.d0

    lcycl = dia*(AR - 1)
    pvol = 1.d0/6.d0*pi*dia**3.d0 + pi*lcycl*(dia/2.d0)**2.d0;
    pmass = pvol*rho
    sigma_c = (0.32d0*(AR**2.d0) + 0.694d0*AR - 0.014)*pi*(dia**2)
    bmax = sqrt(sigma_c / pi)
    ! Moment of Intertia for smooth spherocylinders
    mI(1) = pi/32.d0*rho*dia**4.D0*lcycl + pi/60.D0*rho*dia**5.D0
    mI(2) = pi/48.D0*rho*dia**2.D0*lcycl**3.D0 + 3.D0*pi/64.D0*rho*dia**4.D0*LCYCL+&
      pi/60.D0*rho*dia**5.D0 + pi/24.D0*rho*dia**3.D0*lcycl**2.D0
    mI(3) = mI(2); 

    COR_PP = (ALPHA_PP+1.D0)*0.5D0

    lFILENAME	= "particle_input.dat"
    OPEN(UNIT=lUNIT, FILE=lFILENAME, STATUS="OLD")
    do i = 1,PIP
	read (lunit,*) (des_pos_new(i,k),k=1,3), (des_vel_new(i,k),k=1,3), (des_omega_new(i,k),k=1,3)
    end do

    DO I = 1,PIP 
          Er(I) = 0.5d0 * mI(2) * (DES_OMEGA_NEW(I,2)**2.d0 + DES_OMEGA_NEW(I,3)**2.d0)
    END DO

    DO I = 1,PIP
       K_GRANULAR = K_GRANULAR + &
          DOT_PRODUCT(DES_VEL_NEW(I,:),DES_VEL_NEW(I,:))*PMASS
         DO J = 2,3
            ROT_GRANULAR = ROT_GRANULAR +mI(J)*DES_OMEGA_NEW(I,J)*DES_OMEGA_NEW(I,J)
         END DO
    END DO
    K_GRANULAR = K_GRANULAR/(3.D0*DBLE(PIP))
    ROT_GRANULAR = ROT_GRANULAR/(2.D0*DBLE(PIP))

    print *, "Translational temperature: ", K_GRANULAR
    print *, "Rotational temperature: ", ROT_GRANULAR

    des_pos_new(:,1) = mod(des_pos_new(:,1)+xlength,xlength)
    des_pos_new(:,2) = mod(des_pos_new(:,2)+ylength,ylength)
    des_pos_new(:,3) = mod(des_pos_new(:,3)+zlength,zlength)
    close(lUNIT)

    gmm_file = "/home/mgbolase/desktop/DSMC_Shear3/model/mod/gmm_cond_AR20.bin" 
    print *, "Initializing GMM..."
    call init_gmm_cond(trim(gmm_file))
    call print_gmm_debug()

    call init_pchi(AR, alpha_pp)

    call get_Zr_ref(AR, Zr_ref)

    call get_P1hit(alpha_pp, AR, P1hit)
    
    call get_delta_eps_max(alpha_pp, AR, delta_eps_max)

    return
  end subroutine INIT_PART
	
  ! DSMC
  !-----------------------------------------------------------------
	
  SUBROUTINE INIT_DSMC()
    use geometry, only: xlength, ylength, zlength
    use dsmc
    use particle
    use constants
    use run_param

    implicit none

    integer :: I, J, K, IJK

    !C_SELECT = 0.5d0 * F_N * 4.D0
    C_SELECT = 0.5d0 * F_N

    ! Collision Grid Discretization
    DSMC_DX = XLENGTH/DBLE(NX_DSMC)
    DSMC_DY = YLENGTH/DBLE(NY_DSMC)
    DSMC_DZ = ZLENGTH/DBLE(NZ_DSMC)
    VOL_DSMC = DSMC_DX*DSMC_DY*DSMC_DZ
    OVOL_DSMC = 1.d0/VOL_DSMC
    N_DSMC = Nx_DSMC*Ny_DSMC*Nz_DSMC
	
    allocate(CELL_LIST(N_DSMC))
    allocate(CCELL(N_DSMC))
    allocate(I_OF_CCELL(N_DSMC))
    allocate(J_OF_CCELL(N_DSMC))
    allocate(K_OF_CCELL(N_DSMC))
    do I=1,N_DSMC
       CCELL(I)%Np = 0
       CCELL(I)%NTrial = 0.D0
       CCELL(I)%WMAX = 0.D0
       CCELL(I)%sol_frac = 0.D0
       CCELL(I)%theta = 0.D0
       CELL_LIST(I) = I
    enddo
	
    allocate(PLARRAY(PIP))
    allocate(PART_CELLINDX(PIP))

    ! Compute coefficients for computing IJK from I, J, and K
    IJK_C1=Nx_DSMC; 
    IJK_C2=Nx_DSMC*Ny_DSMC

    ! Do some geometry on cells
    DO I=1, Nx_DSMC
	DO J=1, Ny_DSMC
	   DO K=1, Nz_DSMC
		IJK = DSMC_IJK(I,J,K)
		I_OF_CCELL(IJK)=I
		J_OF_CCELL(IJK)=J
		K_OF_CCELL(IJK)=K
	   ENDDO
	ENDDO
    ENDDO
	
    ! For smoothing solids fraction
    ! Weights
    CX_EPS = DIA/DSMC_DX; 
    CY_EPS = DIA/DSMC_DY; 
    CZ_EPS = DIA/DSMC_DZ

    ! Bounding sphere radius -> number of cells
    NX_EPS = CEILING(CX_EPS); 
    NY_EPS = CEILING(CY_EPS); 
    NZ_EPS = CEILING(CZ_EPS)

    ! Bound weights
    CX_EPS = MIN(CX_EPS,1.d0); 
    CY_EPS = MIN(CY_EPS,1.d0); 
    CZ_EPS = MIN(CZ_EPS,1.d0)

    ! Scaled EPS
    PART_EPS = PVOL*F_N/(DSMC_DX*DSMC_DY*DSMC_DX)
    CALL DSMC_INDEX
    CALL DSMC_INIT_WMAX

    RETURN
  END SUBROUTINE INIT_DSMC
	
  ! DSMC
  !-----------------------------------------------------------------
  SUBROUTINE DSMC_INIT_WMAX
    use dsmc
    use particle
    use constants
	
    implicit none

    integer :: I, J, K
    integer :: NPart
    double precision :: vrel(3), cr
    double precision :: wmax, wmaxtemp, T_temp
	
    WMAX = 0.d0
    !DO I = 2,PIP
    !	vrel = des_vel_new(I,:)*2.D0
    !	cr = sqrt(dot_product(VREL,VREL))
    !	wmaxtemp = cr * sigma_c
    !	IF(WMAXTEMP.GT.WMAX) WMAX = WMAXTEMP
    !END DO
    DO I = 1, PIP-1
       DO J = I+1, PIP
          vrel = DES_VEL_NEW(I,:) - DES_VEL_NEW(J,:)
          cr = sqrt(dot_product(vrel, vrel))
          wmaxtemp = cr * sigma_c
          IF(WMAXTEMP.GT.WMAX) WMAX = WMAXTEMP
       END DO
    END DO
    CCELL(:)%WMAX = WMAX * 1.5d0

    WMAX = 0.d0
    
    DO I = 1,N_DSMC
	WMAX = MAX(WMAX,DBLE(CCELL(I)%Np))
    END DO
    CCELL(:)%WMAX = CCELL(:)%WMAX*WMAX

    return
  end subroutine DSMC_INIT_WMAX
