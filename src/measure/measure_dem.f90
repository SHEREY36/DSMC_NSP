  subroutine measure_dem()	
    use run_param
    use particle
    use output
    use geometry
    use DSMC
	
    implicit none
    integer :: i, j
    double precision :: solfrac, lam
	
    ! Count via collisions
    tau = dble(NCOLL)/PIP
    IF(tau.LT.N_RECORD) RETURN

    !IF (mod(istep, 100) == 0) THEN
    call Pij()
    !call IsoCompressiblity()
    call paramProfile()
    call VACF
    !call EkEt_Cluster()
    !call fv()
    !END IF
    !istep = istep + 1

    PREV_TIME = t
    N_RECORD = (FLOOR(tau/FREQ_RECORD)+1)*FREQ_RECORD
	
    IF(t.GT.tend) SIM_CONTINUE = .FALSE.
	
    ! Adjust time step to hasten at lower temperatures
    ! Conservative timestep (use max solids fraction)
    !SOLFRAC = DBLE(PIP)*PVOL*F_N/(XLENGTH*YLENGTH*ZLENGTH)
    !lam = DIA/(SQRT(2.D0)*6.D0*chi_ma(solfrac)*solfrac)
    !dt = lam/(sqrt(2.d0*granular)*10.d0)

    return	
  end subroutine measure_dem
