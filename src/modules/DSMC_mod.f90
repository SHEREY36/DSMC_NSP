  MODULE DSMC
	
    IMPLICIT NONE
    TYPE :: DSMC_CELL
	! Number of particles in cell
	integer :: Np
	! Number of particles
	! Particle lookup array start index
	integer :: PLA_START
	! Number of Selections
	double precision :: NTrial
	! MVCMAX in each DSMC cell: for all phases
	double precision :: WMAX
	! Solid Fraction
	double precision :: sol_frac
        ! Cell-averaged Temperature ratio
	double precision :: theta
	! Potential
	double precision :: uijk
    END TYPE DSMC_CELL
	
    ! Collision Cell Parameters
    TYPE(DSMC_CELL), DIMENSION(:), ALLOCATABLE ::  CCELL
    INTEGER, DIMENSION(:), ALLOCATABLE :: I_OF_CCELL, J_OF_CCELL, K_OF_CCELL
    INTEGER, DIMENSION(:), ALLOCATABLE :: CELL_LIST
	
    ! Particle Parameters
    DOUBLE PRECISION :: F_N, c_select
    INTEGER, DIMENSION(:), ALLOCATABLE :: PLARRAY
    INTEGER, DIMENSION(:), ALLOCATABLE :: PART_CELLINDX
		
    ! Collision Grid Discretization
    INTEGER :: NX_DSMC, NY_DSMC, NZ_DSMC, N_DSMC
    DOUBLE PRECISION :: DSMC_DX, DSMC_DY, DSMC_DZ
    DOUBLE PRECISION :: VOL_DSMC, OVOL_DSMC
    DOUBLE PRECISION :: IJK_C1, IJK_C2
	
    ! Smoothing Parameters
    DOUBLE PRECISION :: CX_EPS, CY_EPS, CZ_EPS
    INTEGER :: NX_EPS, NY_EPS, NZ_EPS
    DOUBLE PRECISION :: PART_EPS
	
    ! Random Walk Flag
    LOGICAL :: RWALK
	
    CONTAINS
    !-----------------------------------------------------------------	
    ! IJK value of Cell
    integer function DSMC_IJK(I,J,K)

      implicit none
      integer :: I,J,K
      DSMC_IJK = I + IJK_C1*(J-1) + IJK_C2*(K-1)
    end function DSMC_IJK

    !-----------------------------------------------------------------	
    double precision function chi_ma(sol_frac)

      implicit none
      double precision, intent(in) :: sol_frac
      double precision :: sol_frac_ratio
      double precision :: numer, denom
      double precision, parameter :: max_chi = 1000.d0 ! arbitrary limit
      double precision, parameter :: sol_frac_pack = 0.643D0
      double precision, parameter :: C1 = 2.5D0, C2 = 4.5904D0
      double precision, parameter :: C3 = 4.515439D0, C_POW = 0.67802D0

      IF(sol_frac.GT.sol_frac_pack)THEN
	 chi_MA = max_chi
	 RETURN
      ENDIF
      
      numer = 1.d0 + C1*sol_frac + C2*sol_frac**2.d0 + C3*sol_frac**3.d0
      numer = NUMER * 4.0D0*sol_frac
      sol_frac_ratio = sol_frac/sol_frac_pack
      DENOM = (1.d0 - sol_frac_ratio**3.0D0)**C_POW
      chi_ma = 1.d0 + NUMER/DENOM
      chi_ma = MIN(chi_ma,max_chi)

      RETURN
    end function chi_ma
	
  END MODULE DSMC