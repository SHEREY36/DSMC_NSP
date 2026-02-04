module collision_tables_mod
    
    implicit none
    private
    public :: get_delta_eps_max, get_P1hit, init_pchi, sample_chi, p_chi

    integer, parameter :: dp = kind(1.0d0)

    ! Use UNIQUE names to avoid conflicts with other modules
    integer, parameter :: M_ELASTIC = 3   ! m = 0..2  → 3 values
    integer, parameter :: N_ELASTIC = 3   ! n = 0..2  → 3 values
    integer, parameter :: K_ELASTIC = 5   ! k = 0..4  → 5 values

    real(dp), parameter :: beta = 0.5d0

    ! Elastic coefficients - column-major order
    real(dp), parameter :: a_elastic(M_ELASTIC, K_ELASTIC) = reshape((/ &
      1.6873231093532046d-01, -2.3422025542515978d-01,  4.6990593038409800d-02, &
     -2.2335201400824332d+00,  7.7824615377092705d+00, -9.1796465773107383d-01, &
      2.5972229944444980d+01, -2.5437436684757500d+01,  2.3118578511070478d+00, &
     -3.8369445192346824d+01,  2.4542914992197733d+01, -1.3275712470782079d+00, &
      1.4414950493163628d+01, -6.5683799197148618d+00, -1.3361202626719099d-01  &
    /), (/M_ELASTIC, K_ELASTIC/))

    ! Delta-p coefficients - column-major order
    real(dp), parameter :: a_dp(M_ELASTIC, N_ELASTIC, K_ELASTIC) = reshape((/ &
      6.2414983673001086d-02, -5.0945613314919556d-02,  8.0399264059662010d-03, &
     -1.2166411325362669d+00,  6.3399165761527665d-01, -4.1559277887395225d-02, &
      2.6412861315526110d+00, -2.3392938619784158d+00,  1.2668998958113506d-02, &
     -7.3368333108724637d-01,  7.7906830826760654d-01, -1.4628184477457487d-01, &
      2.1824966567939089d+01, -5.5213896773089326d+00,  1.0276731605332390d+00, &
     -1.6414556363006179d+00,  5.5774538557737642d+01, -9.2319168509144927d+00, &
      2.4322928469077385d+00, -2.5784901055570839d+00,  4.8861229976985776d-01, &
     -2.4758309845418431d+01, -1.7172350210975793d+01,  6.2780114321076874d-01, &
     -6.9964611761633819d+01, -2.0607989894238793d+02,  4.0752937305957317d+01, &
     -3.0492092762828156d+00,  2.9615089265599890d+00, -5.3935769938669287d-01, &
     -4.9691898664908692d+01,  7.1411515099552773d+01, -7.5408087306280756d+00, &
      1.5677895634098536d+02,  2.5453928675442623d+02, -5.3978751623663598d+01, &
      1.2815327955067102d+00, -1.0982420739614940d+00,  1.8545748352060740d-01, &
      5.4879389762259606d+01, -5.0011678127184418d+01,  6.0201491869751136d+00, &
     -8.8316725032941164d+01, -1.0252047801167170d+02,  2.2570339021666808d+01  &
    /), (/M_ELASTIC, N_ELASTIC, K_ELASTIC/))

    real(dp) :: AR_c, alpha_c, phi_c, pchi_max
    logical  :: chi_initialized = .false.

  contains

    pure real(dp) function p_chi(chi)
      real(dp), intent(in) :: chi
      integer :: m_loop, n_loop, k_loop
      real(dp) :: s, chi_k

      s = 0.0_dp

      ! Elastic contribution
      do m_loop = 1, M_ELASTIC
        do k_loop = 1, K_ELASTIC
          chi_k = chi**(k_loop-1)
          s = s + a_elastic(m_loop, k_loop) * AR_c**(m_loop-1) * chi_k
        end do
      end do

      ! Delta-p contribution
      do m_loop = 1, M_ELASTIC
        do n_loop = 1, N_ELASTIC
          do k_loop = 1, K_ELASTIC
            chi_k = chi**(k_loop-1)
            s = s + a_dp(m_loop, n_loop, k_loop) * AR_c**(m_loop-1) * phi_c**(n_loop-1) * chi_k
          end do
        end do
      end do

      p_chi = max(s, 0.0_dp)
    end function p_chi

    subroutine init_pchi(AR, alpha)
      real(dp), intent(in) :: AR, alpha
      integer, parameter :: NCHI = 1000
      integer :: i_loop
      real(dp) :: chi, pval

      AR_c = AR
      alpha_c = alpha
      phi_c = 1.0d0 - alpha**beta

      pchi_max = 0.0d0
      do i_loop = 1, NCHI+1
         chi = real(i_loop-1,dp)/real(NCHI,dp)
         pval = p_chi(chi)
         if (pval > pchi_max) pchi_max = pval
      end do

      pchi_max = 1.05d0*pchi_max
      chi_initialized = .true.
    end subroutine init_pchi

    real(dp) function sample_chi()
      real(dp) :: chi_try, u

      if (.not. chi_initialized) then
         write(*,*) "ERROR: chi model not initialized"
         stop
      end if

      do
        call random_number(chi_try)
        call random_number(u)
        if (u*pchi_max <= p_chi(chi_try)) then
           sample_chi = chi_try
           return
        end if
      end do
    end function sample_chi

    subroutine get_delta_eps_max(alpha, AR, delta_eps_max)
      implicit none
      real(8), intent(in)  :: alpha, AR
      real(8), intent(out) :: delta_eps_max

      real(8), parameter :: alpha_vals(9) = (/ &
          0.60d0, 0.65d0, 0.70d0, 0.75d0, 0.80d0, 0.85d0, 0.90d0, 0.95d0, 1.00d0 /)
      real(8), parameter :: AR_vals(5) = (/ &
          1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0 /)

      real(8), parameter :: delta_tab(9,5) = reshape((/ &
      0.62843237d0, 0.56188775d0, 0.49842757d0, 0.42885752d0, 0.35272744d0, &
      0.27604604d0, 0.18801014d0, 0.09553369d0, 0.00002937d0, &
      0.59730676d0, 0.53706389d0, 0.47991723d0, 0.40989182d0, 0.34176128d0, &
      0.25641905d0, 0.17630623d0, 0.09206433d0, 0.00107540d0, &
      0.58311891d0, 0.51520812d0, 0.46263179d0, 0.40238170d0, 0.33524739d0, &
      0.24599273d0, 0.17893030d0, 0.08893593d0, 0.00157158d0, &
      0.56700614d0, 0.52261793d0, 0.45015398d0, 0.38835223d0, 0.33890285d0, &
      0.25566058d0, 0.17456592d0, 0.08878956d0, 0.00153673d0, &
      0.58887764d0, 0.52604992d0, 0.47435499d0, 0.38989718d0, 0.33055193d0, &
      0.25715936d0, 0.17037510d0, 0.08761990d0, 0.00189684d0  &
       /), (/9,5/))

      integer :: iAR, ia, i_loop

      iAR = -1
      do ia = 1, size(AR_vals)
          if (abs(AR - AR_vals(ia)) < 1.0d-12) then
              iAR = ia
              exit
          end if
      end do
      if (iAR < 0) stop "ERROR: AR out of scope"

      ia = -1
      do i_loop = 1, size(alpha_vals)
           if (abs(alpha - alpha_vals(i_loop)) < 1.0d-12) then
               ia = i_loop
               exit
           end if
      end do
      if (ia < 0) stop "ERROR: alpha out of scope"

      delta_eps_max = delta_tab(ia, iAR)
    end subroutine get_delta_eps_max

    subroutine get_P1hit(alpha, AR, P1hit)
      implicit none
      real(8), intent(in)  :: alpha, AR
      real(8), intent(out) :: P1hit

      real(8), parameter :: alpha_vals(9) = (/ &
            0.60d0, 0.65d0, 0.70d0, 0.75d0, 0.80d0, 0.85d0, 0.90d0, 0.95d0, 1.00d0 /)
      real(8), parameter :: AR_vals(5) = (/ &
            1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0 /)

      real(8), parameter :: P1hit_tab(9,5) = reshape((/ &
      1.000000d0, 1.000000d0, 1.000000d0, 1.000000d0, 1.000000d0, &
      1.000000d0, 1.000000d0, 1.000000d0, 1.000000d0, &
      0.845749d0, 0.849408d0, 0.869413d0, 0.878642d0, 0.886644d0, &
      0.892894d0, 0.897983d0, 0.904753d0, 0.912958d0, &
      0.727508d0, 0.749936d0, 0.762986d0, 0.773066d0, 0.783317d0, &
      0.789631d0, 0.798194d0, 0.808411d0, 0.817125d0, &
      0.665665d0, 0.697077d0, 0.713907d0, 0.725973d0, 0.732115d0, &
      0.737128d0, 0.748857d0, 0.759017d0, 0.764082d0, &
      0.650241d0, 0.674177d0, 0.682693d0, 0.706726d0, 0.707292d0, &
      0.714036d0, 0.723566d0, 0.738152d0, 0.745671d0  &
      /), (/9,5/))

      integer :: iAR, ia, i_loop

      iAR = -1
      do ia = 1, size(AR_vals)
          if (abs(AR - AR_vals(ia)) < 1.0d-12) then
              iAR = ia
              exit
          end if
      end do
      if (iAR < 0) stop "ERROR: AR out of scope"

      ia = -1
      do i_loop = 1, size(alpha_vals)
         if (abs(alpha - alpha_vals(i_loop)) < 1.0d-12) then
            ia = i_loop
            exit
         end if
      end do
      if (ia < 0) stop "ERROR: alpha out of scope"

      P1hit = P1hit_tab(ia, iAR)
    end subroutine get_P1hit

end module collision_tables_mod