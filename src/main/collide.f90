  SUBROUTINE DSMC_COLLIDE
    use dsmc
    use run_param
    use particle
    use constants, only: pi
    use output, only: ncoll, pij_c
    use geometry
    use collision_tables_mod
    use gmm_cond_mod

    implicit none

    integer :: I,J,K ! Loop indices
    double precision :: RR
    integer :: IJK,C ! Cell ref
    integer :: p1, p2 ! Particle ref
    ! Selections
    integer :: NSelect
    double precision :: solfrac, chi, Np
    ! Collision partner sampling
    integer :: NPart, NeiPart
    double precision :: u(3)
    ! Collision acceptance
    double precision :: cr, vrel(3)
    double precision :: wij, wmax, wmaxtemp 
    double precision :: eij(3)
    ! Scattering
    double precision :: cor, theta, phi, eps
    double precision :: cos2theta, costheta, sintheta, tantheta
    double precision :: cr_loss_fac
    ! Collisional Pressure
    double precision :: v1(3), v2(3), v1com(3), v2com(3), vcom(3)
    double precision :: Etrans_i, Erot_i, Etrans_f, Erot_f, Etotal_i, Etotal_f
    double precision :: r
    double precision :: chi_angle, Zr, Pr
    double precision :: Rn1, Rn2
    logical :: relax_p1, relax_p2
    double precision :: eps_tr_in, eps_r_in, eps_tr_f, eps_r_f
    double precision :: epsilon_tr_i, epsilon_rot_1_i, epsilon_rot_2_i
    double precision :: epsilon_tr_f, epsilon_rot_1_f, epsilon_rot_2_f 
    double precision :: gamma, cr_new
    double precision, parameter :: a_fit = 1.211d0
    double precision, parameter :: b_fit = 3.672d0

    DO I=1,N_DSMC
        IJK = CELL_LIST(I)

	solfrac = CCELL(IJK)%sol_frac; 
	CHI=CHI_MA(solfrac)

	NPart = CCELL(IJK)%Np
	wmax = CCELL(IJK)%WMAX

        !r = CCELL(IJK)%theta
        !r = nint(r * 10.d0) /10.d0
        r = 1.d0

	CCELL(IJK)%NTrial = CCELL(IJK)%NTrial + c_select * DBLE(NPart)*DBLE(NPART-1) * wmax * OVOL_DSMC * dt !CHI 

	NSelect = FLOOR(CCELL(IJK)%NTrial)
	IF(NSelect.LT.1) CYCLE
	CCELL(IJK)%NTrial = CCELL(IJK)%NTrial - DBLE(NSelect)

	WMAX = CCELL(IJK)%wmax; 
	WMAXTEMP = 0.d0

	do J=1, NSelect
	   ! Particle from Cell
	   CALL RANDOM_NUMBER(RR)
	   P1 = FLOOR(RR*DBLE(NPart))
	   P1 = CCELL(IJK)%PLA_START+MIN(P1,NPart-1)
	   P1 = PLARRAY(P1)
			
	   ! Collision partner
	   call sampleCell(C,IJK,P1,EIJ,u)
	   !C = IJK ! For Bird DSMC (*)

	   DO WHILE(.TRUE.)
		CALL RANDOM_NUMBER(RR)
		P2 = FLOOR(RR*CCELL(C)%Np)+1
		P2 = CCELL(C)%PLA_START+MIN(P2,CCELL(C)%Np-1)
		P2 = PLARRAY(P2)
		IF(P1.NE.P2) exit
	   END DO
			
	   ! Accept/reject based on collisional cross section
	   vrel = DES_VEL_NEW(P1,:) - (DES_VEL_NEW(P2,:) + u)
	   !cr = SQRT(DOT_PRODUCT(vrel,vrel)) (*)
	   cr = DOT_PRODUCT(EIJ,VREL)
	   wij = cr* sigma_c *DBLE(CCELL(C)%Np)
	   IF(wij.GE.wmaxtemp) wmaxtemp = wij
	   CALL RANDOM_NUMBER(RR)
	   IF(wij.LT.(RR*wmax)) CYCLE
	      NCOLL = NCOLL + 2
	      V1 = DES_VEL_NEW(P1,:)
              V2 = DES_VEL_NEW(P2,:)
              
              VCOM = (V1 + V2)*0.5D0
              V1COM = V1 - VCOM
              V2COM = V2 - VCOM

              Etrans_i = 0.5d0 * pmass * (DOT_PRODUCT(V1COM,V1COM) + DOT_PRODUCT(V2COM,V2COM))

              Erot_i = Er(P1) + Er(P2)
              Etotal_i = Etrans_i + Erot_i

              epsilon_tr_i = Etrans_i/Etotal_i
              epsilon_rot_1_i = Er(P1)/Erot_i
              epsilon_rot_2_i = Er(P2)/Erot_i
   
              !Zr = Zr_ref*(0.39d0*(r**2.d0) + 0.09d0*r)
              Zr = 7.29d0
              Pr = 1.d0/Zr
              !Pr = min(Pr, 0.5d0)

              relax_p1 = .false.
              relax_p2 = .false.

              ! particle-selection prohibiting double relaxation
              call random_number(Rn1)
              if (Rn1 < Pr) then
		  eps_tr_in = epsilon_tr_i
                  eps_r_in = epsilon_rot_1_i
                  call sample_conditionals(r, eps_tr_in, eps_r_in, eps_tr_f, eps_r_f)
                  epsilon_tr_f = eps_tr_f
                  epsilon_rot_1_f = eps_r_f
                  epsilon_rot_2_f = 1.d0 - eps_r_f
              else
                  epsilon_tr_f = epsilon_tr_i
                  epsilon_rot_1_f = epsilon_rot_1_i
                  epsilon_rot_2_f = epsilon_rot_2_i
              end if

              !if (Rn1 < Pr) then
              !    relax_p1 = .true.
              !else
              !    call random_number(Rn2)
              !    if (Rn2 < Pr) then
              !        relax_p2 = .true.
              !    end if
              !end if
 
              ! apply relaxation to exactly one particle (or none)
              !if (relax_p1) then
              !    ! inputs are the *current* fractions for particle i (P1)
              !    eps_tr_in = epsilon_tr_i
              !    eps_r_in = epsilon_rot_1_i
              !    call sample_conditionals(r, eps_tr_in, eps_r_in, eps_tr_f, eps_r_f)
              !    epsilon_tr_f = eps_tr_f
              !    epsilon_rot_1_f = eps_r_f
              !    epsilon_rot_2_f = 1.d0 - eps_r_f

              !else if (relax_p2) then
                  ! inputs are the *current* fractions for particle j (P2)
              !    eps_tr_in = epsilon_tr_i
              !    eps_r_in = epsilon_rot_2_i   
              !    call sample_conditionals(r, eps_tr_in, eps_r_in, eps_tr_f, eps_r_f)
              !    epsilon_tr_f = eps_tr_f
              !    epsilon_rot_2_f = eps_r_f
              !    epsilon_rot_1_f = 1.d0 - eps_r_f

              !else
                  ! no relaxation occurs for this pair
              !    epsilon_tr_f = epsilon_tr_i
              !    epsilon_rot_1_f = epsilon_rot_1_i
              !    epsilon_rot_2_f = epsilon_rot_2_i
              !end if

              gamma = sample_dissp(a_fit, b_fit)
              gamma = gamma * delta_eps_max * P1hit

              Etotal_f = Etotal_i * (1.d0 - gamma)

              Etrans_f = Etotal_f * epsilon_tr_f
              Erot_f = Etotal_f - Etrans_f

              Er(P1) = epsilon_rot_1_f * Erot_f
              Er(P2) = epsilon_rot_2_f * Erot_f

              cr_new = SQRT(Etrans_f/pmass)

              ! Sample the post-collision angle
              chi_angle = sample_chi()
              chi_angle = chi_angle * pi

              call random_number(RR)
              eps = 2.d0*pi*RR

              call update_velocities(P1, P2, chi_angle, eps, cr_new)

	      ! Update Collisional Pressure
	      ! Only track first particle contribution and omit 1/2 from Pc
	      V1 = V1 - DES_VEL_NEW(P1,:)
	      Pij_c(1,1) = v1(1)*eij(1) + Pij_c(1,1)
	      Pij_c(1,2) = v1(1)*eij(2) + Pij_c(1,2)
	      Pij_c(1,3) = v1(1)*eij(3) + Pij_c(1,3)
	      Pij_c(2,2) = v1(2)*eij(2) + Pij_c(2,2)
	      Pij_c(2,3) = v1(2)*eij(3) + Pij_c(2,3)
	      Pij_c(3,3) = v1(3)*eij(3) + Pij_c(3,3)

	END DO ! End Selections loop
	IF(CCELL(IJK)%wmax<wmaxtemp) ccell(IJK)%wmax=wmaxtemp*1.2d0

    END DO ! End Cell Loop
    RETURN
  END SUBROUTINE DSMC_COLLIDE

!``````````````````````````````````````````````````````````````````````!
! Subroutine: sampleCell                                               !
!                                                                      !
! Purpose: Find cell to sample from                                    !
!``````````````````````````````````````````````````````````````````````!
  SUBROUTINE sampleCell(nCell, lCell, P, EIJ, u)

    use dsmc
    use geometry
    use particle
    use constants
    use dsmc

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: lCell, P
    INTEGER, INTENT(OUT) :: nCell
    DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: EIJ
	
    INTEGER ::  I, J, K
    DOUBLE PRECISION :: RR
    DOUBLE PRECISION :: THETA, PHI
    DOUBLE PRECISION, DIMENSION(3) :: RIJ
    DOUBLE PRECISION, DIMENSION(3) :: U

10  CALL RANDOM_NUMBER(theta); 
    theta = 2.d0 * pi * theta
    CALL RANDOM_NUMBER(phi); 
    phi = acos(1.d0 - 2.d0*phi)
    eij = (/sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)/)
    eij = eij/sqrt(dot_product(eij,eij))
    rij = DES_POS_NEW(P,:) + eij*bmax

    u = 0.d0
    IF(rij(1).GT.XLENGTH) THEN
	IF(.NOT.PERIODIC(1)) goto 10
	rij(1) = rij(1) - XLENGTH		
	u(2) = -u_bc(1,2); ! Top -> bottom: subtract
	u(3) = -u_bc(1,3)
    ELSE IF(rij(1).LT.0.D0) THEN
	IF(.NOT.PERIODIC(1)) goto 10
	rij(1) = rij(1) + XLENGTH		
	u(2) = u_bc(1,2); ! Bottom -> top: add
	u(3) = u_bc(1,3)
    END IF

    IF(rij(2).GT.YLENGTH.OR.rij(2).LT.0.D0) THEN
	IF(.NOT.PERIODIC(2)) goto 10
	u(1) = u_bc(2,1); u(3) = u_bc(2,3)
	u(1) = sign(u(1),rij(2)); u(3) = sign(u(3),rij(2))
	rij(2) = MOD(rij(2)+YLENGTH,YLENGTH)
    END IF

    IF(rij(3).GT.ZLENGTH.OR.rij(3).LT.0.D0) THEN
	IF(.NOT.PERIODIC(3)) goto 10
	u(1) = u_bc(3,1); u(2) = u_bc(3,2)
	u(1) = sign(u(1),rij(3)); u(2) = sign(u(2),rij(3))
	rij(3) = MOD(rij(3)+ZLENGTH,ZLENGTH)
    END IF

    rij(1) = rij(1)/DSMC_DX
    rij(2) = rij(2)/DSMC_DY
    rij(3) = rij(3)/DSMC_DZ

    I = FLOOR(rij(1))+1; 
    J = FLOOR(rij(2))+1; 
    K = FLOOR(rij(3))+1

    nCell = DSMC_IJK(I,J,K)
    IF(CCELL(nCell)%Np.LT.2&
	.OR.I.GT.NX_DSMC.OR.J.GT.NY_DSMC&
	.OR.K.GT.NZ_DSMC.OR.I.LT.1.OR.J.LT.1.OR.K.LT.1) THEN
	write(*,*) 'Bad cell sample'
	goto 10
    END IF

    RETURN
  END SUBROUTINE sampleCell


!``````````````````````````````````````````````````````````````````````!
! Subroutine: UPDATE_VELOCITIES                                        !
!                                                                      !
! Purpose: Creates particles lists with neighboring cells              !
!``````````````````````````````````````````````````````````````````````!
subroutine UPDATE_VELOCITIES(L, K, chi, eps, crmag)
  use constants, only: SMALL_NUM
  use particle
  implicit none

  ! I/O
  integer, intent(in) :: L, K
  double precision, intent(in) :: chi, eps
  double precision, intent(in) :: crmag   ! post-collision relative speed

  ! Local variables
  double precision :: coschi, sinchi, coseps, sineps
  double precision :: vrwr, norm_crel
  double precision :: ur, vr, wr

  double precision, dimension(3) :: velL, velK
  double precision, dimension(3) :: vcom
  double precision, dimension(3) :: crA
  double precision, dimension(3) :: crel

  ! Trig
  coschi = cos(chi)
  sinchi = sin(chi)
  coseps = cos(eps)
  sineps = sin(eps)

  ! Load velocities
  velL = des_vel_new(L,:)
  velK = des_vel_new(K,:)

  ! Center-of-mass velocity
  vcom = 0.5d0*(velL + velK)

  ! Relative velocity of particle L
  crA = velL - vcom

  ur = crA(1)
  vr = crA(2)
  wr = crA(3)

  vrwr = sqrt(vr*vr + wr*wr)

  ! Build scattered relative direction
  if (vrwr >= SMALL_NUM) then
     crel(1) = coschi*ur + sinchi*sineps*vrwr
     crel(2) = coschi*vr + sinchi*(crmag*wr*coseps - ur*vr*sineps)/vrwr
     crel(3) = coschi*wr - sinchi*(crmag*vr*coseps + ur*wr*sineps)/vrwr
  else
     crel(1) = coschi*vr + sinchi*(crmag*coseps - ur*sineps)
     crel(2) = coschi*wr - sinchi*(crmag*coseps + ur*sineps)
     crel(3) = 0.0d0
  end if

  ! Normalize crel
  norm_crel = sqrt(dot_product(crel, crel))
  if (norm_crel < SMALL_NUM) then
     write(*,*) "ERROR: crel normalization failed"
     stop
  end if
  crel = crel / norm_crel

  ! Update velocities (post-collision)
  des_vel_new(L,:) = vcom + crel * crmag
  des_vel_new(K,:) = vcom - crel * crmag

end subroutine UPDATE_VELOCITIES
