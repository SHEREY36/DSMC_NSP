  subroutine integrate_eom	
   use particle
   use run_param
   use geometry
   use dsmc
   use constants
   !USE output, ONLY: dx_meas, dy_meas, dz_meas, NX_MEAS, NY_MEAS, NZ_MEAS, bulk, NpCell
	
   implicit none
   integer :: I, J, K, L, M 
   integer :: P, PSTART, CC
   integer :: IJK, IJKNEW
   double precision :: r(3), r_dim
   double precision :: solfrac, chi
   double precision :: solfrac_old, solfrac_new
   double precision :: u_old, u_new, d_pot
   double precision :: RR, probCross
   ! Collisional loss
   double precision :: dE_total, GranTemp
   ! Thermostats
   !double precision :: zeta, prob_replace, vboost, stddev
   !DOUBLE PRECISION :: vrand(3), u1, u2

   !-----------------------------------------------------------------
   ! Traditional Streaming
   IF(.NOT.RWALK) THEN
      DES_POS_NEW(:,1) = DES_POS_NEW(:,1) + DES_VEL_NEW(:,1)*dt
      DES_POS_NEW(:,2) = DES_POS_NEW(:,2) + DES_VEL_NEW(:,2)*dt
      DES_POS_NEW(:,3) = DES_POS_NEW(:,3) + DES_VEL_NEW(:,3)*dt
	
      ! Currently do not have NSW
      ! Only moving boundaries in Lees-Edwards for shear flows
      ! X - Dim
      IF(PERIODIC_X) THEN
	WHERE(DES_POS_NEW(:,1).GT.XLENGTH) 
	   DES_POS_NEW(:,1) = MOD(DES_POS_NEW(:,1) + XLENGTH,XLENGTH)
	   DES_VEL_NEW(:,2) = DES_VEL_NEW(:,2) - 2.d0*u_bc(1,2)
	   DES_VEL_NEW(:,3) = DES_VEL_NEW(:,3) - 2.d0*u_bc(1,3)
	END WHERE

	WHERE(DES_POS_NEW(:,1).LT.0.D0)
	   DES_POS_NEW(:,1) = MOD(DES_POS_NEW(:,1) + XLENGTH,XLENGTH)
	   DES_VEL_NEW(:,2) = DES_VEL_NEW(:,2) + 2.d0*u_bc(1,2)
           DES_VEL_NEW(:,3) = DES_VEL_NEW(:,3) + 2.d0*u_bc(1,3)	
	END WHERE
	
      ELSE ! Reflect
	DO I = 1,PIP
	   r_dim = DES_POS_NEW(I,1)
	   IF(R_DIM.LT.DIA) THEN
	      r_dim = dia + (dia-r_dim)
	      des_vel_new(I,1) = -des_vel_new(I,1)
	   ELSE IF(R_DIM.LT.(XLENGTH-dia)) THEN
	      r_dim = xlength - dia - mod(r_dim,xlength-dia)
	      des_vel_new(I,1) = -des_vel_new(I,1)			
	   END IF
	END DO
      END IF

      ! Y - Dim
      IF(PERIODIC_Y) THEN
	 WHERE(DES_POS_NEW(:,2).GT.YLENGTH)
	    DES_POS_NEW(:,2) = MOD(DES_POS_NEW(:,2) + YLENGTH,YLENGTH)
	    DES_VEL_NEW(:,1) = DES_VEL_NEW(:,1) - u_bc(2,1)
	    DES_VEL_NEW(:,3) = DES_VEL_NEW(:,3) - u_bc(2,3)
         END WHERE

	 WHERE(DES_POS_NEW(:,2).LT.0.D0)
	    DES_POS_NEW(:,2) = MOD(DES_POS_NEW(:,2) + YLENGTH,YLENGTH)
	    DES_VEL_NEW(:,1) = DES_VEL_NEW(:,1) + u_bc(2,1)
	    DES_VEL_NEW(:,3) = DES_VEL_NEW(:,3) + u_bc(2,3)	
	 END WHERE
      ELSE
	 DO I = 1,PIP
	    r_dim = DES_POS_NEW(I,2)
	    IF(R_DIM.LT.DIA) THEN
	       r_dim = dia + (dia-r_dim)
	       des_vel_new(I,2) = -des_vel_new(I,2)
	    ELSE IF(R_DIM.LT.(YLENGTH-dia)) THEN
	       r_dim = ylength - dia - mod(r_dim,ylength-dia)
	       des_vel_new(I,2) = -des_vel_new(I,2)			
	    END IF
	 END DO
      END IF
	
      ! Z - Dim
      IF(PERIODIC_Z) THEN
	 WHERE(DES_POS_NEW(:,3).GT.ZLENGTH)
	    DES_POS_NEW(:,3) = MOD(DES_POS_NEW(:,3) + ZLENGTH,ZLENGTH)
	    DES_VEL_NEW(:,1) = DES_VEL_NEW(:,1) - u_bc(3,1)
	    DES_VEL_NEW(:,2) = DES_VEL_NEW(:,2) - u_bc(3,2)
	 END WHERE

	 WHERE(DES_POS_NEW(:,3).LT.0.D0)
	    DES_POS_NEW(:,3) = MOD(DES_POS_NEW(:,3) + ZLENGTH,ZLENGTH)
	    DES_VEL_NEW(:,1) = DES_VEL_NEW(:,1) + u_bc(3,1)
	    DES_VEL_NEW(:,2) = DES_VEL_NEW(:,2) + u_bc(3,2)	
	 END WHERE
      ELSE
	 DO I = 1,PIP
	    r_dim = DES_POS_NEW(I,3)
	    IF(R_DIM.LT.DIA) THEN
	       r_dim = dia + (dia-r_dim)
	       des_vel_new(I,3) = -des_vel_new(I,3)
	    ELSE IF(R_DIM.LT.(ZLENGTH-dia)) THEN
	       r_dim = zlength - dia - mod(r_dim,zlength-dia)
	       des_vel_new(I,3) = -des_vel_new(I,3)			
	    END IF
	 END DO
      END IF
	
      DO I = 1,PIP
	 IF(DES_POS_NEW(I,1).GT.XLENGTH.OR.DES_POS_NEW(I,2).GT.YLENGTH&
		.OR.DES_POS_NEW(I,3).GT.ZLENGTH.OR.ANY(DES_POS_NEW(I,:).LT.0.d0)) THEN
		write(*,*) des_pos_new(I,:)
		write(*,*) des_Vel_new(I,:)
		stop
	 end if
      END DO

   !-----------------------------------------------------------------
   ! Random Walk
   ELSEIF(RWALK) THEN
      DO IJK = 1,N_DSMC
	 CCELL(IJK)%sol_frac = CCELL(IJK)%Np*PART_EPS
      END DO

      ! Did not reduce number of simulators needed to converge in collision rate
      ! for high densities
      DO L = 1, N_DSMC
  	 IJK = CELL_LIST(L)
  	 I = I_OF_CCELL(IJK); 
	 J = J_OF_CCELL(IJK); 
	 K = K_OF_CCELL(IJK)
    	 PSTART = CCELL(IJK)%PLA_START - 1
    	 DO M = 1,CCELL(IJK)%Np
     	    P = PLARRAY(PSTART+M)
     	    r = des_pos_new(P,:) + des_vel_new(P,:)*dt

     	    r(1) = mod(r(1)+xlength,xlength)
     	    r(2) = mod(r(2)+ylength,ylength)
     	    r(3) = mod(r(3)+zlength,zlength)

     	    IJKNEW = DSMC_IJK(FLOOR(r(1)/dsmc_dx)+1, floor(r(2)/dsmc_dy)+1, floor(r(3)/dsmc_dz)+1)

            ! If particle enters new cell
            IF(IJKNEW.NE.IJK) THEN
		solfrac_old = CCELL(IJK)%sol_frac!-PART_EPS
		solfrac_new = CCELL(IJKNEW)%sol_frac! + PART_EPS

		u_old = 4.D0*solfrac_old*chi_ma(solfrac_old)
		u_new = 4.D0*solfrac_new*chi_ma(solfrac_new)

		d_pot = u_new-u_old
		probCross = min(1.D0,exp(-d_pot))

		CALL RANDOM_NUMBER(RR)
		IF(probCross.GT.RR) THEN
		   DES_POS_NEW(P,:) = r
		ELSE
	           CALL RANDOM_NUMBER(r)
		   r(1) = (I - r(1))*DSMC_DX
	           r(2) = (J - r(2))*DSMC_DY
	           r(3) = (K - r(3))*DSMC_DZ
		   DES_POS_NEW(P,:) = r
		END IF
	    ELSE
	     	CALL RANDOM_NUMBER(r)
		r(1) = (I - r(1))*DSMC_DX
	     	r(2) = (J - r(2))*DSMC_DY
	     	r(3) = (K - r(3))*DSMC_DZ

		DES_POS_NEW(P,:) = r	
	    END IF
    	 END DO
      END DO
   END IF
   RETURN
 end subroutine integrate_eom
