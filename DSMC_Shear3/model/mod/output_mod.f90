  module output
	
    implicit none
	
    INTEGER :: NX_MEAS, NY_MEAS, NZ_MEAS, N_MEAS
    DOUBLE PRECISION :: DX_MEAS, DY_MEAS, DZ_MEAS
    !INTEGER :: ISTEP
		
    ! Timing
    DOUBLE PRECISION :: N_RECORD
    DOUBLE PRECISION :: DT_RECORD, PREV_TIME, FREQ_RECORD = 1.0D0
    DOUBLE PRECISION :: PREV_GRANULAR
    ! Collision Frequency
    INTEGER :: NCOLL
    DOUBLE PRECISION :: TAU, Tg_total
    DOUBLE PRECISION :: INIT_GRANULAR, GRANULAR, Ersum, Rot_Gran
    ! Pressure Tensor
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: BULK
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::NpCell
    DOUBLE PRECISION, DIMENSION(3,3) :: PIJ_K, PIJ_C
    ! 1d Profile
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vprof, Tgprof, nprof
    DOUBLE PRECISION :: DCELL
    DOUBLE PRECISION :: VOL_MEAS, OVOL_MEAS
    INTEGER :: DIM_MEAS, N_PROF_MEAS, DIM_VEL
    ! Instabiltiy Detection
    DOUBLE PRECISION :: EkEt0, NSQ0
    DOUBLE PRECISION :: mTgErr, Err0
    ! Distribution Function
    DOUBLE PRECISION :: VMIN, VMAX, BINWIDTH
    INTEGER :: BINS
    INTEGER, DIMENSION(:), ALLOCATABLE :: FREQ
    ! VACF
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: V0
	
    contains

    !-----------------------------------------------------------------		
    subroutine paramProfile
      use geometry
      use particle
      use dsmc, only: F_N

      implicit none
      INTEGER :: I,J
      DOUBLE PRECISION :: VOL_MEAS
	
      IF(SIZE(vprof).EQ.1) RETURN
      vprof = 0.d0 
      nprof = 0.d0 
      tgprof = 0.d0
      DO I = 1,PIP
	 J = FLOOR(DES_POS_NEW(I,DIM_MEAS)/DCELL) + 1
	 vprof(J) = vprof(J) + des_vel_new(I,DIM_VEL)
	 tgprof(J) = tgprof(J) + DOT_PRODUCT(des_vel_new(I,:),DES_VEL_NEW(I,:))
	 nprof(J) = nprof(J) + 1
      END DO
	
      DO I = 1,N_PROF_MEAS
	 vprof(I) = vprof(I)/dble(nprof(I))
	 tgprof(I) = Tgprof(I)/(3.d0*dble(nprof(I)))
      END DO
      write(1011,*) vprof
      write(1012,*) nprof*F_N*PVOL*OVOL_MEAS
      write(1013,*) tgprof
      flush(1011); flush(1012); flush(1013) ! Not necessary
      return
    end subroutine paramProfile

    !-----------------------------------------------------------------
    ! Pressure Tensor	
    subroutine Pij
      use particle
      use run_param
      use dsmc
	
      implicit none
      integer :: I, J, K, P
      double precision :: vrel(3)
      double precision :: GranularoInit, GRANULAR_TH
      double precision :: dTg, TgErr
      pij_k = 0.d0
      Ersum = 0.d0
      ! Compute Stress Tensor
      DO I = 1,PIP
	 pij_k(1,1) = DES_VEL_NEW(I,1)*DES_VEL_NEW(I,1) + pij_k(1,1)
	 pij_k(1,2) = DES_VEL_NEW(I,1)*DES_VEL_NEW(I,2) + pij_k(1,2)
	 pij_k(1,3) = DES_VEL_NEW(I,1)*DES_VEL_NEW(I,3) + pij_k(1,3)
	 pij_k(2,2) = DES_VEL_NEW(I,2)*DES_VEL_NEW(I,2) + pij_k(2,2)
	 pij_k(2,3) = DES_VEL_NEW(I,2)*DES_VEL_NEW(I,3) + pij_k(2,3)
	 pij_k(3,3) = DES_VEL_NEW(I,3)*DES_VEL_NEW(I,3) + pij_k(3,3)
         Ersum = Ersum + Er(I)
      END DO
      Rot_Gran = Ersum / dble(PIP)
      write(1005,'(6(d10.4,2x))') pij_k(1,1), pij_k(2,2), pij_k(3,3),&
	    pij_k(1,2), pij_k(1,3), pij_k(2,3)
      IF(PREV_TIME.NE.0.D0) then
	 pij_c = pij_c / (t-PREV_TIME) * bmax
	 write(1006,'(6(d10.4,2x))') pij_c(1,1), pij_c(1,2), pij_c(1,3),&
		pij_c(2,2), pij_c(2,3), pij_c(3,3)
      END IF
      pij_c = 0.d0
	
      ! Compute Bulk Velocities
      bulk = 0.d0
      NpCell = 0.d0
      DO P = 1,PIP
	 I = FLOOR(DES_POS_NEW(P,1)/DX_MEAS) + 1
	 J = FLOOR(DES_POS_NEW(P,2)/DY_MEAS) + 1
	 K = FLOOR(DES_POS_NEW(P,3)/DZ_MEAS) + 1
	 bulk(I,J,K,:) = bulk(I,J,K,:) + DES_VEL_NEW(P,:)
	 NpCell(I,J,K) = NpCell(I,J,K) + 1
      ENDDO
      
      DO I = 1,NX_MEAS
	 DO J = 1,NY_MEAS
	    DO K = 1,NZ_MEAS
	       bulk(I,J,K,:) = bulk(I,J,K,:)/NpCell(I,J,K)
	    END DO
	 END DO
      END DO
	
      ! Compute Granular Temperature
      !GRANULAR = 0.d0
      !DO P = 1,PIP
      !    I = FLOOR(DES_POS_NEW(P,1)/DX_MEAS) + 1
      !    J = FLOOR(DES_POS_NEW(P,2)/DY_MEAS) + 1
      !    K = FLOOR(DES_POS_NEW(P,3)/DZ_MEAS) + 1	
      !    VREL = DES_VEL_NEW(P,:) - BULK(I,J,K,:)
      !    GRANULAR = GRANULAR + DOT_PRODUCT(VREL,VREL)
      !END DO
      ! Tried to assess with 1D, did not seem like a good metric for quasi 1D
      !GRANULAR = GRANULAR/(3.D0*DBLE(PIP))

      GRANULAR = (pij_k(1,1)+pij_k(2,2)+pij_k(3,3))/(3.D0*DBLE(PIP))
      Tg_total = (3.d0*GRANULAR + 2.d0*Rot_Gran)/5.d0

      IF(PREV_TIME.EQ.0.D0) INIT_GRANULAR = GRANULAR
      GranularoInit = Granular/INIT_GRANULAR
      GRANULAR_TH = exp(-(1.d0-ALPHA_PP**2.d0)/3.d0*tau)
      dTg = log(GranularoInit) - log(Granular_TH); 
      dTg = ABS(dTg)
      IF(PREV_TIME.EQ.0.D0) Err0 = MAX(dTg,0.01D0) 
      dTg = dTg - Err0
      TgErr = mTgErr*tau; 
      TgErr = ABS(TgErr) + Err0
	
      write(1002,'(8(E10.4,2x))') t, TAU, Rot_Gran, GRANULAR, Tg_total, &
	    GranularoInit, dTg, TgErr
      IF(dTg.GT.TgErr.AND.tau.gt.20) LC_FOUND = .TRUE.

      flush(1002); flush(1005); flush(1006);

      return
    end subroutine Pij

    !-----------------------------------------------------------------
    ! Reduced Velocity Distribution
    subroutine fv
      use particle
      use run_param
	
      implicit none
      integer :: I, J
      DOUBLE PRECISION :: VScaled, SQGranular

      ! 0.5d0 m/s spacing
      ! min/max bin = [-10 10]
      SQGranular = SQRT(GRANULAR); 
      FREQ = 0
      DO I = 1,PIP
	 VScaled = ABS(DES_VEL_NEW(I,1))/SQGranular
	 J = FLOOR(VSCALED/BINWIDTH) + 1
	 J = MAX(J,1); 
         J = MIN(J,BINS)
	 FREQ(J) = FREQ(J) + 1
      END DO
      write(1003,'(101(d10.4,2x))') DBLE(FREQ)/PIP
      flush(1003)
      return
    end subroutine fv

    !-----------------------------------------------------------------
    ! Instability Detection
    subroutine EkEt_Cluster
      use particle
      use run_param
      use dsmc
	
      implicit none
      integer :: I, J, K
      integer :: p
      DOUBLE PRECISION, DIMENSION(3) :: VREL
      DOUBLE PRECISION :: Et, EkEt, NSQ
      DOUBLE PRECISION :: dEkEt, dNSQ

      Et = 0.d0; NSQ = 0.d0
      DO I = 1,NX_MEAS
	 DO J = 1,NY_MEAS
	    DO K = 1,NZ_MEAS
	       Et = Et + DOT_PRODUCT(bulk(I,J,K,:),bulk(I,J,K,:))
	       NSQ = NSQ + NpCell(I,J,K)**2.D0
	    END DO
	 END DO
      END DO
      
      EkEt = granular/Et
      IF(PREV_TIME.EQ.0.D0) THEN
	 EkEt0 = EkEt
	 NSQ0 = NSQ
      END IF
      ! IF either the bulk energy to temperature ratio
      ! or number fluctuations exceed 10% of the inital
      ! value, instability present
      dEkEt = ABS(EkEt-EkEt0)/EkEt0; dNSQ = (NSQ-NSQ0)/NSQ0
      write(1007,'(2(d10.4,2x))') EkEt/EkEt0, NSQ/NSQ0
      flush(1007)

      return
    end subroutine EkEt_Cluster

    !-----------------------------------------------------------------
    ! Isothermal Compressibility
    subroutine IsoCompressiblity
      use particle
      use run_param
      use dsmc
	
      implicit none
      integer :: I, J, K
      DOUBLE PRECISION :: NSQ, NMean, K_ISO

      NSQ = 0.d0
      DO I = 1,NX_MEAS
         DO J = 1,NY_MEAS
	    DO K = 1,NZ_MEAS
	       NSQ = NSQ + NpCell(I,J,K)**2.D0
	    END DO
	 END DO
      END DO

      NMean = DBLE(PIP)/DBLE(N_MEAS)
      K_ISO = (NSQ/DBLE(N_MEAS) - NMean**2.D0)/NMean
      write(1008,'(101(d10.4,2x))') K_ISO
      flush(1008)

      return
    end subroutine IsoCompressiblity

    !-----------------------------------------------------------------
    ! Velocity AutoCorrelation Function (VACF)
    subroutine VACF
      use particle
      use run_param
      use dsmc
	
      implicit none
      integer :: I, J, K
      DOUBLE PRECISION :: V0Vt

      V0Vt = 0.D0
      DO I = 1,PIP
	 V0Vt = DOT_PRODUCT(V0(I,:),DES_VEL_NEW(I,:)) + V0Vt
      END DO
      write(1009,'(101(d10.4,2x))') V0Vt/PIP
      flush(1009)
      return
    end subroutine VACF

  end module output	
