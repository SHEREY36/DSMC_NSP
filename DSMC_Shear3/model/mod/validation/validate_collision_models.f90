program validate_collision_models
    use gmm_cond_mod
    use collision_tables_mod
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: NCURVE = 1000
    integer, parameter :: NSAMP = 40000
    
    ! File handling variables
    character(len=256) :: gmm_file
    character(len=256) :: input_file
    character(len=256) :: out_file
    integer :: u_in, u_out, ios, n_in, i, j
    logical :: have_inputs
    
    ! Data variables
    real(dp), allocatable :: etr_in(:), er1_in(:)
    real(dp) :: r, e_tr, e_r1, eps_tr_f, eps_r1_f, u

    real(dp) :: AR, alpha, x, pval
    real(dp) :: P1hit, delta_eps_max

    integer :: ia, iAR
    integer :: u_chi, u_pchi, u_p1, u_de

    AR    = 2.0d0
    alpha = 0.95d0

    open(newunit=u_chi,  file="chi_samples.txt",        status="replace")
    open(newunit=u_pchi, file="pchi_curve.txt",         status="replace")
    open(newunit=u_p1,   file="p1hit_table.txt",        status="replace")
    open(newunit=u_de,   file="delta_eps_max_table.txt",status="replace")

    ! 1. Setup filenames
    gmm_file = "gmm_cond_AR20.bin"
    input_file = "inputs_r05.txt"
    out_file = "gmm_samples.txt"

    ! 2. Initialize GMM and Check
    print *, "Initializing GMM..."
    call init_gmm_cond(trim(gmm_file))
    call print_gmm_debug()  

    ! 3. Check Input File
    inquire(file=trim(input_file), exist=have_inputs)
    if (.not. have_inputs) then
        print *, "ERROR: Input file not found: ", trim(input_file)
        stop
    end if

    ! 4. Count Lines
    print *, "Reading inputs from: ", trim(input_file)
    open(newunit=u_in, file=trim(input_file), status="old", action="read", iostat=ios)
    if (ios /= 0) stop "Error opening input file"

    n_in = 0
    do
        read(u_in, *, iostat=ios) e_tr, e_r1
        if (ios /= 0) exit
        n_in = n_in + 1
    end do
    print *, "Found ", n_in, " input lines."
    rewind(u_in)

    ! 5. Allocate and Load
    allocate(etr_in(n_in), er1_in(n_in))
    
    do i = 1, n_in
        read(u_in, *) etr_in(i), er1_in(i)
    end do
    close(u_in)

    ! 6. Run Validation Loop
    print *, "Running validation sampling..."
    open(newunit=u_out, file=trim(out_file), status="replace", action="write")
    
    r = 0.5_dp 

    do i = 1, 20000
        ! Randomly select a row from input data
        call random_number(u)
        j = 1 + int(u * real(n_in, dp))
        if (j < 1) j = 1
        if (j > n_in) j = n_in

        e_tr = etr_in(j)
        e_r1 = er1_in(j)

        ! Sample
        call sample_conditionals(r, e_tr, e_r1, eps_tr_f, eps_r1_f)

        ! Write: r, input_tr, input_r1, pred_tr, pred_r1
        write(u_out, '(5(E14.6, 1x))') r, e_tr, e_r1, eps_tr_f, eps_r1_f
    end do

    close(u_out)
    print *, "Done! Results saved to ", trim(out_file)

    call random_seed()
    ! 1) p(chi) curve + samples
    call init_pchi(AR, alpha)

    do i = 0, NCURVE
       x = real(i,dp)/real(NCURVE,dp)
       pval = p_chi(x)          ! NOTE: p_chi expects x in [0,1]
       write(u_pchi,'(F12.8,1X,ES16.8)') x, pval
    end do

    do i = 1, NSAMP
       x = sample_chi()         ! returns x in [0,1]
       write(u_chi,'(F12.8)') x
    end do

    ! 3) P1hit table
    do iAR = 1, 5
       AR = 1.0d0 + 0.5d0*(iAR-1)
       do ia = 0, 8
          alpha = 0.60d0 + 0.05d0*ia
          call get_P1hit(alpha, AR, P1hit)
          write(u_p1,'(2F6.2,1X,F12.6)') AR, alpha, P1hit
       end do
    end do

    ! 4) delta_eps_max table
    do iAR = 1, 5
       AR = 1.0d0 + 0.5d0*(iAR-1)
       do ia = 0, 8
          alpha = 0.60d0 + 0.05d0*ia
          call get_delta_eps_max(alpha, AR, delta_eps_max)
          write(u_de,'(2F6.2,1X,F12.6)') AR, alpha, delta_eps_max
       end do
    end do

    close(u_chi); close(u_pchi); close(u_p1); close(u_de)
    print *, "Validation run completed."
    print *, "Files written:"
    print *, "  chi_samples.txt (x in [0,1])"
    print *, "  pchi_curve.txt (x, p(x))"
    print *, "  p1hit_table.txt"
    print *, "  delta_eps_max_table.txt"

end program validate_collision_models