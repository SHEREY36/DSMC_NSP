  program main
    use dsmc
    use run_param
    use particle, only: pip
    use output

    implicit none
				
    write(*,*) 'Reading Inputs'
    call read_input()

    write(*,*) 'Initializing'
    call initialize()

    write(*,*) 'Recording Initial Conditions'	
    call measure_init()
    t = 0.D0
    sim_continue = .TRUE.

    write(*,*) 'Performing Collisions'
    !istep = 0
    DO WHILE(SIM_CONTINUE)
	call dsmc_index()
	call dsmc_collide()
	t = t + dt
	call integrate_eom()
	call measure_dem()
        write(*,*) 't =', t, 's' 
	write(*,*) 'Collisions per particle: ', DBLE(NCOLL)/PIP	
    END DO
    write(*,*) 'Recording Final State'
    call measure_final()

    return
  end program main