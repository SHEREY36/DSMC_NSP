  subroutine measure_init()	
    use particle
    use constants
    use dsmc
    use geometry
    use output
    use run_param
	
    implicit none
	
    ! Initial Collision Frequency
    INTEGER :: I,J
    INTEGER :: SZ
    DOUBLE PRECISION :: VMEAN
    DOUBLE PRECISION :: VSQ(3)
    DOUBLE PRECISION :: CHI, max_sol_frac
	
    ! Output Files
    !---------------------------------------------------------
    open(unit=1001, status="replace", file="nu.txt")
    open(unit=1002, status="replace", file="tg.txt")
    open(unit=1003, status="replace", file="VoV0.txt")
    open(unit=1005, status='replace', file='Pijk.txt')
    open(unit=1006, status='replace', file='Pijc.txt')
    open(unit=1007, status='replace', file='ICS.txt')
    open(unit=1008, status='replace', file='KIso.txt')
    open(unit=1009, status='replace', file='VACF.txt')
    open(unit=1011, status='replace', file='vprof.txt')
    open(unit=1012, status='replace', file='nprof.txt')
    open(unit=1013, status='replace', file='tgprof.txt')

    ! Counters for recording frequency
    !---------------------------------------------------------	
    NCOLL = 0; 
    DT_RECORD = 0.1D0; 
    N_RECORD = 0.0D0
    PREV_TIME = 0.d0
   
    ! Discretize domain for measuring
    !---------------------------------------------------------	
    ! Measuring Grid
    N_MEAS = NX_MEAS*NY_MEAS*NZ_MEAS
    DX_MEAS = XLENGTH/DBLE(NX_MEAS)
    DY_MEAS = YLENGTH/DBLE(NY_MEAS)
    DZ_MEAS = ZLENGTH/DBLE(NZ_MEAS)

    IF(NX_MEAS.GT.1) THEN
	SZ = NX_MEAS; 
	DCELL = DX_MEAS;
	N_PROF_MEAS = NX_MEAS; 
	DIM_MEAS = 1

	IF(U_BC(1,2).NE.0.d0) THEN
		DIM_VEL = 2
	ELSEIF(U_BC(1,3).NE.0.d0) THEN
		DIM_VEL = 3
	ELSE
		DIM_VEL = 1
	END IF
		VOL_MEAS = DX_MEAS*YLENGTH*ZLENGTH

    ELSEIF(NY_MEAS.GT.1) THEN
	SZ = NY_MEAS; 
	DCELL = DY_MEAS;
	N_PROF_MEAS = NY_MEAS; 
	DIM_MEAS = 2

	IF(U_BC(2,1).NE.0.d0) THEN
		DIM_VEL = 1
	ELSEIF(U_BC(2,3).NE.0.d0) THEN
		DIM_VEL = 3
	ELSE
		DIM_VEL = 2
	END IF
		VOL_MEAS = DY_MEAS*XLENGTH*ZLENGTH

    ELSEIF(NZ_MEAS.GT.1) THEN
	SZ = NZ_MEAS; 
	DCELL = DZ_MEAS;
	N_PROF_MEAS = NZ_MEAS; 
	DIM_MEAS = 3

	IF(U_BC(3,1).NE.0.D0) THEN
		DIM_VEL = 1
	ELSEIF(U_BC(3,2).NE.0.D0) THEN
		DIM_VEL = 2
	ELSE
		DIM_VEL = 3
	END IF
		VOL_MEAS = DZ_MEAS*XLENGTH*YLENGTH

    ELSE
	SZ = 1; 
	DCELL = 1;
    END IF
	
    OVOL_MEAS = 1.D0/VOL_MEAS
	
    ! Set up output arrays
    !---------------------------------------------------------		
    pij_k = 0.d0; 
    pij_c = 0.d0

    IF(.NOT.ALLOCATED(vprof))	ALLOCATE(vprof(sz))
    IF(.NOT.ALLOCATED(nprof))	ALLOCATE(nprof(sz))
    IF(.NOT.ALLOCATED(Tgprof))	ALLOCATE(Tgprof(sz))
    IF(.NOT.ALLOCATED(bulk))	ALLOCATE(bulk(NX_MEAS,NY_MEAS,NZ_MEAS,3))
    IF(.NOT.ALLOCATED(NpCell))	ALLOCATE(NpCell(NX_MEAS,NY_MEAS,NZ_MEAS))
    IF(.NOT.ALLOCATED(V0))	ALLOCATE(V0(PIP,3))
	
    ! Reduced Velocity Distribution
    VMIN = 0; 
    VMAX = 10; 
    BINWIDTH = 0.1D0
    BINS = FLOOR((VMAX - VMIN)/BINWIDTH) + 1
    IF(.NOT.ALLOCATED(FREQ))	ALLOCATE(FREQ(BINS))

    ! Granular Temperature margins
    mTgErr = -(1.d0-ALPHA_PP**2.d0)/3.d0*0.01d0
		
    ! Flag for when critical length found
    LC_FOUND = .FALSE.
	
    ! VACF
    DO I = 1,PIP
	V0(I,:) = DES_VEL_NEW(I,:)
    END DO

    return
  end subroutine measure_init
	
