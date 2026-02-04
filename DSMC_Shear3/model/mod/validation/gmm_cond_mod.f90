module gmm_cond_mod
    implicit none
    private
    public :: init_gmm_cond, sample_conditionals, print_gmm_debug

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: Dx = 3, Dy = 2, D = 5

    ! ---- numerical safety ----
    real(dp), parameter :: eps_logit = 1.0d-8

    integer :: M_gmm = 0  ! RENAMED to avoid conflicts with other modules
    logical :: gmm_ready = .false.

    real(dp), allocatable :: w(:)              ! (M_gmm)
    real(dp), allocatable :: means(:,:)        ! (M_gmm,5)
    real(dp), allocatable :: inv_xx(:,:,:)     ! (M_gmm,3,3)
    real(dp), allocatable :: logdet_xx(:)      ! (M_gmm)
    real(dp), allocatable :: A(:,:,:)          ! (M_gmm,2,3)
    real(dp), allocatable :: mu_y(:,:)         ! (M_gmm,2)
    real(dp), allocatable :: L(:,:,:)          ! (M_gmm,2,2)
    real(dp), allocatable :: sc_mean(:)        ! (5)
    real(dp), allocatable :: sc_scale(:)       ! (5)

    contains

  subroutine init_gmm_cond(binfile)
    character(len=*), intent(in) :: binfile
    integer :: unit, ios
    integer :: header(2)

    unit = 77
    open(unit=unit, file=binfile, access="stream", form="unformatted", status="old", action="read", iostat=ios)
    if (ios /= 0) stop "ERROR: cannot open GMM binary file."

    read(unit) header
    M_gmm = header(1)
    ! header(2) == D expected

    allocate(w(M_gmm), means(M_gmm,D), inv_xx(M_gmm,Dx,Dx), logdet_xx(M_gmm))
    allocate(A(M_gmm,Dy,Dx), mu_y(M_gmm,Dy), L(M_gmm,Dy,Dy))
    allocate(sc_mean(D), sc_scale(D))

    read(unit) w
    read(unit) means
    read(unit) inv_xx
    read(unit) logdet_xx
    read(unit) A
    read(unit) mu_y
    read(unit) L
    read(unit) sc_mean
    read(unit) sc_scale

    close(unit)

    gmm_ready = .true.
  end subroutine init_gmm_cond

  subroutine sample_conditionals(r, e_tr, e_r1, eps_tr_f, eps_r1_f)
    ! Inputs are your physical/unprocessed values:
    real(dp), intent(in)  :: r, e_tr, e_r1
    real(dp), intent(out) :: eps_tr_f, eps_r1_f

    real(dp) :: x_un(3), x_proc(3), x_sc(3)
    real(dp) :: dxv(3)
    real(dp) :: logp, maxlog, sumexp, u
    real(dp), allocatable :: loga(:)
    integer :: msel, m, i, j

    real(dp) :: mu_x(3), yc(2), z(2)
    real(dp) :: invS(3,3)

    if (.not. gmm_ready) stop "ERROR: GMM not initialized. Call init_gmm_cond first."

    allocate(loga(M_gmm))

    ! --- build x (unprocessed)
    x_un = (/ r, e_tr, e_r1 /)

    call preprocess_x(x_un, x_proc)
    call scale_x(x_proc, x_sc)

    ! --- compute log
    maxlog = -huge(1.0_dp)
    do m = 1, M_gmm
      mu_x = means(m,1:Dx)
      dxv = x_sc - mu_x
      
      ! Manually copy strided slice to contiguous local array
      do j=1,3
         do i=1,3
            invS(i,j) = inv_xx(m,i,j)
         end do
      end do
      
      logp = log(w(m)) + log_gauss3(dxv, invS, logdet_xx(m))
      loga(m) = logp
      if (logp > maxlog) maxlog = logp
    end do

    sumexp = 0.0_dp
    do m = 1, M_gmm
      sumexp = sumexp + exp(loga(m) - maxlog)
    end do

    call random_number(u)
    if (sumexp <= 1.0e-300_dp) then
      ! fallback uniform
      msel = 1 + int(u * real(M_gmm, dp))
      if (msel < 1) msel = 1
      if (msel > M_gmm) msel = M_gmm
    else
      ! sample categorical
      call sample_categorical_loga(loga, maxlog, sumexp, u, msel)
    end if

    ! Safety check: ensure msel is within valid bounds
    if (msel < 1 .or. msel > M_gmm) then
       write(*,*) "ERROR: invalid msel=", msel, " M_gmm=", M_gmm
       write(*,*) "u=", u, " sumexp=", sumexp, " maxlog=", maxlog
       msel = max(1, min(msel, M_gmm))
    end if

    ! --- conditional mean: mu_y + A*(x - mu_x)
    mu_x = means(msel,1:Dx)
    dxv = x_sc - mu_x
    yc(1) = mu_y(msel,1) + dot_product(A(msel,1,:), dxv)
    yc(2) = mu_y(msel,2) + dot_product(A(msel,2,:), dxv)

    ! --- sample y = yc + L*z, z~N(0,I)
    call randn2(z(1), z(2))
    ! L is lower-triangular
    yc(1) = yc(1) + L(msel,1,1)*z(1)
    yc(2) = yc(2) + L(msel,2,1)*z(1) + L(msel,2,2)*z(2)

    call postprocess_y(yc, eps_tr_f, eps_r1_f)

    deallocate(loga)

  end subroutine sample_conditionals

  pure subroutine scale_x(x_proc, x_sc)
    real(dp), intent(in)  :: x_proc(3)
    real(dp), intent(out) :: x_sc(3)
    integer :: i
    do i=1,3
      x_sc(i) = (x_proc(i) - sc_mean(i)) / sc_scale(i)
    end do
  end subroutine scale_x

  pure real(dp) function log_gauss3(dxv, invS, logdetS) result(val)
    real(dp), intent(in) :: dxv(3)
    real(dp), intent(in) :: invS(3,3)
    real(dp), intent(in) :: logdetS
    real(dp) :: q
    q = dxv(1)*(invS(1,1)*dxv(1) + invS(1,2)*dxv(2) + invS(1,3)*dxv(3)) &
      + dxv(2)*(invS(2,1)*dxv(1) + invS(2,2)*dxv(2) + invS(2,3)*dxv(3)) &
      + dxv(3)*(invS(3,1)*dxv(1) + invS(3,2)*dxv(2) + invS(3,3)*dxv(3))
    val = -0.5_dp*(3.0_dp*log(2.0_dp*acos(-1.0_dp)) + logdetS + q)
  end function log_gauss3

  subroutine sample_categorical_loga(loga, maxlog, sumexp, u, msel)
    real(dp), intent(in) :: loga(:), maxlog, sumexp, u
    integer, intent(out) :: msel
    real(dp) :: c, target
    integer :: m, Mloc
    
    Mloc = size(loga)
    target = u * sumexp
    c = 0.0_dp
    
    do m=1,Mloc
      c = c + exp(loga(m) - maxlog)
      if (c >= target) then
        msel = m
        return
      end if
    end do

    msel = Mloc
    
  end subroutine sample_categorical_loga

  subroutine randn2(z1, z2)
    ! Box-Muller
    real(dp), intent(out) :: z1, z2
    real(dp) :: u1, u2, r, theta
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, 1.0e-300_dp)
    r = sqrt(-2.0_dp*log(u1))
    theta = 2.0_dp*acos(-1.0_dp)*u2
    z1 = r*cos(theta)
    z2 = r*sin(theta)
  end subroutine randn2

  subroutine preprocess_x(x_un, x_proc)
    use iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in)  :: x_un(3)
    real(real64), intent(out) :: x_proc(3)
    real(real64) :: z

    ! ---- r : log ----
    x_proc(1) = log(max(x_un(1), eps_logit))

    ! ---- e_tr : logit ----
    z = min(max(x_un(2), eps_logit), 1.0d0 - eps_logit)
    x_proc(2) = log(z / (1.0d0 - z))

    ! ---- e_r1 : logit ----
    z = min(max(x_un(3), eps_logit), 1.0d0 - eps_logit)
    x_proc(3) = log(z / (1.0d0 - z))
  end subroutine preprocess_x

  subroutine postprocess_y(y_scaled, eps_tr_f, eps_r1_f)
    use iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in)  :: y_scaled(2)
    real(real64), intent(out) :: eps_tr_f, eps_r1_f
    real(real64) :: y_proc(2)
    real(real64) :: tmp

    ! ---- inverse StandardScaler ----
    y_proc(1) = y_scaled(1) * sc_scale(4) + sc_mean(4)
    y_proc(2) = y_scaled(2) * sc_scale(5) + sc_mean(5)

    ! ---- inverse preprocess (sigmoid) ----
    if (y_proc(1) > 100.0_dp) then
       eps_tr_f = 1.0_dp
    else if (y_proc(1) < -100.0_dp) then
       eps_tr_f = 0.0_dp
    else
       eps_tr_f = 1.0d0 / (1.0d0 + exp(-y_proc(1)))
    end if

    if (y_proc(2) > 100.0_dp) then
       eps_r1_f = 1.0_dp
    else if (y_proc(2) < -100.0_dp) then
       eps_r1_f = 0.0_dp
    else
       eps_r1_f = 1.0d0 / (1.0d0 + exp(-y_proc(2)))
    end if

    ! ---- hard safety clamp ----
    eps_tr_f = min(max(eps_tr_f, eps_logit), 1.0d0 - eps_logit)
    eps_r1_f = min(max(eps_r1_f, eps_logit), 1.0d0 - eps_logit)

  end subroutine postprocess_y

  subroutine print_gmm_debug()
    integer :: i
    
    if (.not. gmm_ready) then
        print *, "Debug Error: GMM not initialized!"
        return
    end if

    print *, "================ GMM DEBUG INFO ================"
    print *, "M (Components): ", M_gmm
    print *, "Dimensions: ", Dx, Dy
    print *, "--- Component 1 ---"
    print *, "Weight: ", w(1)
    print *, "Mean (first 3): ", means(1, 1:3)
    print *, "L Matrix (Row 1): ", L(1, 1, 1), L(1, 1, 2)
    print *, "L Matrix (Row 2): ", L(1, 2, 1), L(1, 2, 2)
    print *, "--- Component ", M_gmm, " ---"
    print *, "Weight: ", w(M_gmm)
    print *, "Mean (first 3): ", means(M_gmm, 1:3)
    print *, "================================================"

  end subroutine print_gmm_debug

end module gmm_cond_mod