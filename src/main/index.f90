  SUBROUTINE DSMC_INDEX
    use dsmc
    use particle
    use constants

    implicit none

    integer :: I, J, IJK
    integer :: K, TEMP
    integer :: N
    integer :: PREV_START, PREV_END
    integer :: P, pid
    double precision :: POS(3), DCELL(3)
    double precision :: Eksum, Ersum
    double precision :: T_trans, T_rot
    double precision :: RR

    ! Zero out values
    DO I = 1,N_DSMC
       CCELL(I)%Np=0
       PART_CELLINDX(I) = 0	
    END DO

    ! Find cell each particle belongs to
    DCELL = (/DSMC_DX, DSMC_DY, DSMC_DZ/)
    DO I = 1, PIP
       POS = DES_POS_NEW(I,:)/DCELL
       IF(ANY(DES_POS_NEW(I,:).LT.0)) then
	  write(*,*) DES_POS_NEW(I,:), DCELL, I
	  stop
       END IF
       IJK = DSMC_IJK(FLOOR(POS(1))+1, FLOOR(POS(2))+1, FLOOR(POS(3))+1)
       CCELL(IJK)%Np = CCELL(IJK)%Np + 1; 
       PART_CELLINDX(I) = IJK
    ENDDO

    CALL CALC_SOLID_FRACTION()

    ! Create particle reference array
    PREV_START = 1
    DO I = 1, N_DSMC
	CCELL(I)%PLA_START = PREV_START
	PREV_START = PREV_START + CCELL(I)%Np
	CCELL(I)%Np = 0
    END DO
		
    DO I = 1, PIP
       IJK = PART_CELLINDX(I)
       J = CCELL(IJK)%PLA_START + CCELL(IJK)%Np
       PLARRAY(J) = I; 
       CCELL(IJK)%Np = CCELL(IJK)%Np + 1
    END DO

    ! Compute cell-averaged translational / rotational temperature ratio

    DO IJK = 1, N_DSMC
       IF (CCELL(IJK)%Np <= 0) THEN
          CCELL(IJK)%theta = 0.d0
          CYCLE
       END IF

       Eksum = 0.d0
       Ersum = 0.d0

       DO J = 1, CCELL(IJK)%Np
          pid = PLARRAY(CCELL(IJK)%PLA_START + J - 1)
          Eksum = Eksum + 0.5d0 * pmass * &
              DOT_PRODUCT(DES_VEL_NEW(pid,:), DES_VEL_NEW(pid,:))
          Ersum = Ersum + Er(pid)
       END DO

       ! Translational temperature
       T_trans = (2.d0 / (3.d0 * DBLE(CCELL(IJK)%Np))) * Eksum

       ! Rotational temperature
       T_rot = Ersum / DBLE(CCELL(IJK)%Np)

       ! Temperature ratio used in collision model
       IF (T_trans > SMALL_NUM) THEN
           CCELL(IJK)%theta = T_rot / T_trans
       ELSE
           CCELL(IJK)%theta = 0.d0
       END IF
    END DO
	
    ! Randomize order cells are visited
    DO I = 1,N_DSMC
       CALL RANDOM_NUMBER(RR)
       TEMP = CELL_LIST(I)
       IJK = FLOOR(RR*N_DSMC)+1
       CELL_LIST(I) = CELL_LIST(IJK)
       CELL_LIST(IJK) = TEMP
    END DO
	
    ! Randomize particle order within each cell
    DO IJK = 1,N_DSMC
	 I = CCELL(IJK)%PLA_START
	 N = CCELL(IJK)%Np
	 DO J = 1,N
	    TEMP = PLARRAY(I+J-1)
	    CALL RANDOM_NUMBER(RR)
	    K = CEILING(RR*N)
	    PLARRAY(I+J-1) = PLARRAY(I+K-1)
	    PLARRAY(I+K-1) = TEMP
         END DO
    END DO

    RETURN
  END SUBROUTINE DSMC_INDEX
	
  !---------------------------------------------------------------
  ! Approximate local solid fraction by smoothing all solid fractions
  ! of cells within the bounding sphere (radius = particle diameter)
  ! Approximate bounding sphere as a cube (some errors in corners)
  !---------------------------------------------------------------

  SUBROUTINE CALC_SOLID_FRACTION()
    use dsmc
    use particle, only: dia, pvol
    use geometry
    use constants

    IMPLICIT NONE

    INTEGER :: C, L
    INTEGER :: I, J, K, IJK
    INTEGER :: ICELL, JCELL, KCELL
    DOUBLE PRECISION :: X,Y,Z
    INTEGER :: NSample = 25
    DOUBLE PRECISION :: TOTAL_NP, NpMean
    DOUBLE PRECISION :: theta, phi, EIJ(3)

    DO C = 1,N_DSMC
       I = I_OF_CCELL(C); 
       J = J_OF_CCELL(C); 
       K = K_OF_CCELL(C)

       TOTAL_NP = 0
       DO L = 1,NSample
10        CALL RANDOM_NUMBER(theta); 
          theta = 2.d0 * pi * theta

	  CALL RANDOM_NUMBER(phi); 
          phi = acos(1.d0 - 2.d0*phi)
	  eij = (/sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)/)
	  eij = eij/sqrt(dot_product(eij,eij)); 
          eij = eij*dia

      	  CALL RANDOM_NUMBER(X); 
          X = (DBLE(I)-X)*DSMC_DX+EIJ(1)
      	  CALL RANDOM_NUMBER(Y); 
          Y = (DBLE(J)-Y)*DSMC_DY+EIJ(2)
      	  CALL RANDOM_NUMBER(Z); 
          Z = (DBLE(K)-Z)*DSMC_DZ+EIJ(3)

      	  ICELL = FLOOR(X/DSMC_DX) + 1
      	  JCELL = FLOOR(Y/DSMC_DY) + 1
      	  KCELL = FLOOR(Z/DSMC_DZ) + 1

      	  IF(.NOT.PERIODIC_X.AND.((ICELL.GT.NX_DSMC).OR.(ICELL.LT.1))) goto 10
      	  IF(.NOT.PERIODIC_Y.AND.((JCELL.GT.NY_DSMC).OR.(JCELL.LT.1))) goto 10
      	  IF(.NOT.PERIODIC_Z.AND.((KCELL.GT.NZ_DSMC).OR.(KCELL.LT.1))) goto 10

   	  ICELL = MOD(ICELL+NX_DSMC,NX_DSMC); IF(ICELL.EQ.0) ICELL = NX_DSMC
	  JCELL = MOD(JCELL+NY_DSMC,NY_DSMC); IF(JCELL.EQ.0) JCELL = NY_DSMC
	  KCELL = MOD(KCELL+NZ_DSMC,NZ_DSMC); IF(KCELL.EQ.0) KCELL = NZ_DSMC

      	  IJK = DSMC_IJK(ICELL,JCELL,KCELL)
      	  TOTAL_NP = TOTAL_NP + DBLE(CCELL(IJK)%Np)
       END DO
       NpMean = TOTAL_NP/DBLE(NSample)
       CCELL(C)%sol_frac = NpMean*PVOL*F_N*OVOL_DSMC
    END DO

    RETURN
  END SUBROUTINE CALC_SOLID_FRACTION
