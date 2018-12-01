!!$module precision
!!$  !_________________________________________________________________
!!$  !
!!$  ! Flags controlling the processes included in the calculations
!!$  !_________________________________________________________________
!!$  !
!!$  implicit none
!!$  integer,      parameter :: dp=selected_real_kind(p=15) ! real kind
!!$  integer,      parameter :: long=kind(1)                ! integer kind
!!$  real(kind=dp),parameter :: minus_infinity = -huge(1._dp) ! smallest < 0. number
!!$  real(kind=dp),parameter :: plus_infinity  =  huge(1._dp) ! greatest > 0. number
!!$end module precision

module module_parameters_flags
  !_________________________________________________________________
  !
  ! Flags controlling the processes included in the calculations
  !_________________________________________________________________
  !
  implicit none
  logical :: verbose     ! Display extensive information on terminal output
  logical :: do_thermal  ! Compute thermal losses from line emission (subroutine thermal_loss)
  logical :: do_shock    ! Compute evolution of MHD variables
  logical :: do_h2       ! Compute evolution of H2 levels
end module module_parameters_flags
!
!
module module_parameters_output
  !_________________________________________________________________
  !
  !
  !_________________________________________________________________
  !
  implicit none
  !
  ! Standard output : UNIX -> 6, MAC -> 9
  integer(kind=4), parameter :: screen = 6
  !
  ! Outputs control
  logical           :: speci_out_ad,speci_out_cd,speci_out_fd,speci_out_nr
  character(len=1), parameter :: comment="#"
  !
  ! Logical units
  integer(4) :: file_phys
  integer(4) :: file_spec_ad
  integer(4) :: file_spec_cd
  integer(4) :: file_spec_fd
  integer(4) :: file_spec_nr
  integer(4) :: file_h2_lev
  integer(4) :: file_h2_lin
  character(len=7) :: h2_lev_out
  character(len=10) :: h2_lin_out
  integer(4) :: file_cooling,file_energetics,file_intensity,file_populations,file_fe_pops,file_fe_lines
end module module_parameters_output
!
!
module module_gamma
  use module_precision
  !_________________________________________________________________
  !
  ! contains gamma (ratio of specific heats) and useful related
  ! variables.
  !_________________________________________________________________
  !
  implicit none
!  include "precision.f90"
  !
  ! gamma=5/3 for 3 degrees of freedom
  real(kind=dp), private, parameter :: freedom_deg=3._dp ! degrees of freedom
  real(kind=dp), parameter :: gamma=(freedom_deg+2._dp)/freedom_deg
  real(kind=dp), parameter :: gamma1=1._dp/(gamma-1._dp)
  real(kind=dp), parameter :: gamma2=gamma/(gamma-1._dp)
  real(kind=dp), parameter :: gamma3=0.5_dp*(gamma+1._dp)/(gamma-1._dp)
  !
  ! variables from "precision.f90" are private
!  private :: dp, long, minus_infinity, plus_infinity
end module module_gamma


module module_debug
  use module_precision
  !_________________________________________________________________
  !
  ! Control debugging
  !_________________________________________________________________
  implicit none
!  include "precision.f90"
  !
  logical            :: debug                   ! controls debugging (from input file)
  integer, parameter :: debug_ul_diffun_phys=81 ! logical unit for diffun outputs file
  integer, parameter :: debug_ul_diffun_chem=82 ! logical unit for diffun outputs file
  integer, parameter :: debug_ul_diffun_h2=83 ! logical unit for diffun outputs file
  integer :: call2diffun
  real(kind=dp) :: tr_1, tr_2, tr_3, tr_4, tr_5, tr_6
  real(kind=dp) :: tr_7, tr_8, tr_9, tr_10, tr_11, tr_12
  real(kind=dp) :: tr_13, tr_14, tr_15, tr_16, tr_17, tr_18, tr_19
  real(kind=dp) :: tr_20, tr_21, tr_22, tr_23, tr_24, tr_25, tr_26
  real(kind=dp) :: tr_27, tr_28, tr_29, tr_30
  !
end module module_debug


module module_grains
  use module_precision
  !_________________________________________________________________
  !
  ! Variables related to the dust grains
  !_________________________________________________________________
  !
  implicit none
  !
  real(kind=dp) :: grmin            ! [   cm] min. radius of the MRN distribution
  real(kind=dp) :: grmax            ! [   cm] max. radius
  real(kind=dp) :: rho_grain        ! [g/cm3] volumic mass
  real(kind=dp) :: rho_mantle       ! [g/cm3] volumic mass
  real(kind=dp) :: ratio_grain_gas=0.01_dp               ! default mass ratio (grains/gas)
  real(kind=dp) :: r_grain          ! [   cm] sqrt(<a^2>)
  real(kind=dp) :: rgrain1          ! [  cm2] <a>
  real(kind=dp) :: rgrain2          ! [  cm2] <a^2>
  real(kind=dp) :: rgrain3          ! [  cm3] <a^3>
  real(kind=dp) :: r_gr_scal23      ! [     ] <a^2> = r_gr_scale * <a^3>^2/3
  real(kind=dp) :: r_gr_scal12      ! [     ] = <a>/<a^2>
  real(kind=dp) :: dens_grain_init  ! [ cm-3] initial grain number density; calculated in initialize
  real(kind=dp) :: dens_grain       ! [ cm-3] grain number density calculated in diffun
  real(kind=dp) :: md_grain         ! [g/cm3] mass of grains per cm3 of gas
  real(kind=dp) :: vd_grain         ! [cm3/cm3] volume of grains per cm3 of gas
  real(kind=dp) :: rho_grain_charge ! [g/cm3] mass density of charged grains
  real(kind=dp) :: tgrain           ! [    K] grain temperature read in read_parameters
  real(kind=dp) :: teff_grain       ! [    K] effective temperature for sputtering reactions
  real(kind=dp) :: nlayers_grain    ! [     ] number of layers in grain mantle
  real(kind=dp), parameter :: nsites_grain = 1.0d6  ! number of sites/grain
  !
end module module_grains


module module_phys_var
  use module_precision
  !___________________________________________________________________
  !
  ! This module contains all the physical variable
  ! * yarray: all variables in the following order
  !   1/ MHD variables (nv_mhd)
  !   2/ denisty of each species (nspec)
  !   3/ density of h2 levels (nlevels_h2)
  ! * v_lvariab = log(yarray)
  !___________________________________________________________________
  !
  implicit none
  !
  integer(kind=long) :: nv_mhd  ! number of mhd variables (initialized in mhd)
  integer(kind=long) :: dimtot  ! total size of the Y array
  real(kind=dp),dimension(:), allocatable :: yarray
  real(kind=dp),dimension(:), allocatable :: yn
  !
  ! neutrals
  real(kind=dp) :: timen=0.d0  ! [  year] flow time calculated in mhd
  real(kind=dp) :: rhon        ! [g.cm-3] mass density
  real(kind=dp) :: densityn    ! [  cm-3] numeric density
  real(kind=dp) :: mun         ! [     g] mean mass
  real(kind=dp) :: tn          ! [     K] temperature read in read_parameters
  real(kind=dp) :: nh_init     ! [  cm-3] initial proton density, used in read_species
  real(kind=dp) :: vn          ! [  cm/s] velocity
  real(kind=dp) :: dvn         ! [??????] neutral velocity gradient
  real(kind=dp) :: xsphere     ! [      ] fractional radius of collapsing sphere
  real(kind=dp) :: xll         ! [    cm] characteristic viscous length
  real(kind=dp) :: grad_v      ! [??????] neutral velocity gradient too (with viscosity...)
  !
  ! positive ions
  real(kind=dp) :: timei=0.0_dp ! flow time (year) calculated in mhd
  real(kind=dp) :: rhoi        ! mass density (g.cm-3)
  real(kind=dp) :: densityi    ! numeric density (cm-3)
  real(kind=dp) :: mui         ! mean mass (g)
  real(kind=dp) :: ti          ! temperature (k)
  real(kind=dp) :: vi          ! velocity (cm/s)
  real(kind=dp) :: dvi         ! ions velocity gradient
  !
  ! negative ions
  real(kind=dp) :: rhoneg      ! mass density (g.cm-3)
  real(kind=dp) :: densityneg  ! numeric density (cm-3)
  real(kind=dp) :: muneg       ! mean mass (g)
  !
  ! electrons
  REAL(KIND=DP) :: Te            ! [     K] temperature (for e- and negative ions)
  !
  ! Timescales: dynamical, ambipolar diffusion (i_n, charged grain_neutral), collapse
  real(kind=dp) :: tau_ff        ! [     s] free-fall time
  real(kind=dp) :: tff_mult      ! [      ] Free fall timescale saling factor (input)
  real(kind=dp) :: tau_collapse  ! [     s] collapse timescale
  real(kind=dp) :: tau_dyn       ! [    yr] dynamical timescale
  real(kind=dp) :: tau_i_n       ! [    yr] timescale ambipolar diffusion i_n
  real(kind=dp) :: tau_grain_n   ! [    yr] timescale ambipolar diffusion charged grain_neutral
  !
  ! Environment
  real(kind=dp) :: rad  ! multipicative factor for the flux radiation
  real(kind=dp) :: av   ! exctinction (magnitudes)
  real(kind=dp) :: zeta ! cosmic ray ionization rate (s-1)
  !
  ! Dynamics parameters
  real(kind=dp) :: nh                ! [  cm-3] density of protons nh=n(h)+2n(h2)+n(h+)
  real(kind=dp) :: deltav            ! [  cm/s] vi-vn
  real(kind=dp) :: abs_deltav        ! [  cm/s] abs(vi-vn)
  real(kind=dp) :: vgrad             ! [km/s/cm] abs. value of the gradient of vn
  real(kind=dp) :: deltavmin=1d3     ! [  cm/s] if abs_deltav < deltavmin : stop
  real(kind=dp) :: vmagnet           ! [  cm/s] magnetosonic velocity
  real(kind=dp) :: vsound            ! [  cm/s] sound speed
  character(len=1)   :: shock_type   ! [      ] shock type ('C', 'J', 'S')
  integer(kind=long) :: nfluids      ! [      ] number of fluids (neutrals, ions, e-)
  real(kind=dp)      :: vs_km        ! [  km/s] shock velocity
  real(kind=dp)      :: vs_cm        ! [  cm/s] shock velocity
  real(kind=dp)      :: zinit        ! [    cm] for collapse: max radius
  real(kind=dp)      :: timej        ! [    yr] shock age
  logical            :: viscosity    ! [      ] .true. if we need artificial viscosity
  integer            :: force_i_c=1  ! [      ] flag=1: then force densityi = sum(ions)
  integer            :: cool_kn=0    ! [      ] flag=1: Kaufman & Neufeld cooling (h20 & co)
  !
  ! H2 levels and related variables
  real(kind=dp)      :: op_h2       ! h2-ortho/para, initialized in read_parameters
  integer(kind=long) :: nh2_lev     ! number of levels of h2 included
  integer(kind=long) :: nh2_lines_out  ! max number of h2 lines in output file
  character(len=4)   :: h_h2_flag   ! collision data (DRF/Flower, MM/Martin&Mandy, BOTH: DRF if possible and MM if not
  integer            :: iforh2 = 1  ! flag : h2 formation on grains 
                                    !   0: 1/3 of 4.4781 ev in internal energy (17249 k) (Allen, 1999)
                                    !   1: proportional to boltzman distribution at 17249 k
                                    !   2: dissociation limit : v = 14, j = 0,1 (4.4781 ev)
                                    !   3: v = 6, j = 0,1
                                    !   4: fraction = relative populations at t, initialised as h2_lev%density
                                    !                 and changed during integration
  integer            :: ikinh2 = 1  ! flag : h2 formation energy released as kinetic energy
                                    !   1: 0.5 * (4.4781 - internal)
                                    !   2: inf(1.4927 ev, 4.4781 - internal)
  real(kind=dp)     :: sum_h2
  !
  ! Magnetic field (gauss)
  real(kind=dp):: bfield
  real(kind=dp):: bbeta
  !
  ! Integration parameters
  integer(long) :: counter = 0      ! [      ] counts the number of calls to DVODE
  integer(8)    :: nstep_w          ! [      ] number of steps between 2 outputs
  integer(8)    :: nstep_max        ! [      ] max. number of integration steps
  real(kind=dp) :: distance=0.0_dp  ! [    cm] distance from beginning of the shock (cm)
  real(kind=dp) :: dist_step=0.0_dp ! [    cm] step between 2 consecutive calls to drive
  real(kind=dp) :: maxtime          ! [    yr] max. duration of the integration
  real(kind=dp) :: maxdist          ! [    cm] max. distance of cell particle / length of the shock
  real(kind=dp) :: masscloud        ! [    cm] max. distance of cell particle / length of the shock
  real(kind=dp) :: atol0            ! [      ] initial accuracy (argument of dvode)
  real(kind=dp) :: rtol             ! [      ] initial accuracy (argument of dvode)
  real(kind=dp) :: t0_v   = 0.0d0   ! [      ] initial value of the integration variable
  real(kind=dp) :: h0_v   = 1.d05   ! [      ] initial step legth
  real(kind=dp) :: tout   = 1.d01   ! [      ] value of t at wich output should occur
  real(kind=dp) :: tout_step        ! [      ] value of t at wich output should occur
  !
  ! Some usefule indices allowing to extract variables from the 2 vectors
  integer(kind=long) :: iv_xsphere, iv_vn, iv_vi
  integer(kind=long) :: iv_rhon, iv_rhoi, iv_rhoneg
  integer(kind=long) :: iv_tn, iv_ti, iv_te
  integer(kind=long) :: iv_densityn, iv_densityi
  integer(kind=long) :: iv_gv, iv_dens_grain
  integer(kind=long) :: bv_speci,ev_speci
  integer(kind=long) :: bv_neu,ev_neu
  integer(kind=long) :: bv_ion,ev_ion
  integer(kind=long) :: bv_neg,ev_neg
  integer(kind=long) :: bv_gra,ev_gra
  integer(kind=long) :: bv_cor,ev_cor
  integer(kind=long) :: bv_h2_lev,ev_h2_lev
  !
  ! Variables from "precision.f90" are private
!  private :: dp, long, minus_infinity, plus_infinity
  !
contains
  !
  subroutine read_parameters
    use module_grains, only    : tgrain,grmin,grmax,rho_grain,rho_mantle
    use module_constants, only : yearsec,msun
    use module_tools, only     : get_file_number
    use module_parameters_flags
    use module_parameters_output
    use module_debug
    !_________________________________________________________________
    !
    ! Read parameters from input file
    !_________________________________________________________________
    !
    implicit none
    character(len=*), parameter :: name_file_input='input/collapse.in'
    integer                     :: file_input
    character(len=32)           :: charact
    !
    !--- opening file ---
    file_input=get_file_number()
    open(file_input,file=name_file_input,status='old',&
         access='sequential',form='formatted',action='read')
    !
    !
    read(file_input,'(/)') !******************* Basic parameters
    read(file_input,*) tn
    read(file_input,*) nh_init
    read(file_input,*) zeta
    read(file_input,*) rad
    read(file_input,*) av
    read(file_input,*) tff_mult
    read(file_input,*) op_h2
    read(file_input,*) maxtime
    read(file_input,*) masscloud
    read(file_input,*) vs_km
    read(file_input,'(/)') !******************* Grain parameters
    read(file_input,*) grmin
    read(file_input,*) grmax
    read(file_input,*) rho_grain
    read(file_input,*) rho_mantle
    read(file_input,*) tgrain
    if (grmin.eq.grmax) then
       grmax = grmin*1.005d0
       grmin = grmax/(1.005d0)**2
    endif
    if (tgrain.le.0) then
       tgrain = tn ! Grain temperature set to gas temperature if Tgrain.le.0
    endif
    read(file_input,'(/)') !******************* Numerical parameters
    read(file_input,*) atol0
    read(file_input,*) rtol
    read(file_input,*) tout_step
    read(file_input,*) nstep_w
    read(file_input,*) nstep_max
    read(file_input,*) force_i_c
    read(file_input,'(/)') !******************* Output parameters
    read(file_input,'(l)') speci_out_ad
    read(file_input,'(l)') speci_out_cd
    read(file_input,'(l)') speci_out_fd
    read(file_input,'(l)') speci_out_nr
    read(file_input,'(/)') !******************* Flags
    read(file_input,'(l)') verbose
    read(file_input,'(l)') debug
    read(file_input,'(l)') do_thermal
    read(file_input,'(l)') do_shock
    read(file_input,'(l)') do_h2
    read(file_input,'(/)') !******************* Thermal balance
    read(file_input,*)     cool_kn
    read(file_input,'(/)') !******************* H2 parameters
    read(file_input,*) nh2_lev
    read(file_input,*) nh2_lines_out
    read(file_input,'(a4)') h_h2_flag
    read(file_input,*) iforh2
    read(file_input,*) ikinh2
    read(file_input,'(a7)') h2_lev_out
    read(file_input,'(a10)') h2_lin_out
    read(file_input,'(/)') !******************* Dynamics parameters
    read(file_input,*) shock_type
    read(file_input,*) nfluids
    read(file_input,*) bbeta
    read(file_input,*) deltavmin
    read(file_input,*) xll
    !******************* END
    read(file_input,'(/,a)') charact
    if (charact.ne."END") then
       write(*,*) "E-read_parameters: Incorrect input file",charact
       stop
    endif
    close(file_input)
    !
    ! Sanity checks
    if (do_thermal) then
       write(*,*) "E-read_parameters: do_thermal: option not allowed; call to DGESV was disabled in line_excit.f90"
       stop
    endif
    if (iforh2 < 0 .or. iforh2 > 4) then
      write(*,*) "E-read_parameters: iforh2 =", iforh2, "  is not allowed"
      stop
    endif
    if (ikinh2 /= 1 .and. ikinh2 /= 2) then
      write(*,*) "E-read_parameters: ikinh2 =", ikinh2, "  is not allowed"
      stop
    endif
    if (H_H2_flag /= "DRF" .AND. H_H2_flag /= "MM" .AND. H_H2_flag /= "BOTH") then
      write(*,*) "E-read_parameters: Wrong flag for H-H2 collisions: ",H_H2_flag
      stop
    endif
    !
    ! Initialize dynamics parameters
    masscloud = masscloud*msun ! conversion to gram
!!$    if (vs_km.le.0) then
!!$       if (maxdist.le.0) then
!!$          write(*,*) "E-read_parameters, invalid velocity or distance input. Can not be both zero."
!!$          stop
!!$       else
!!$          vs_cm = maxdist/maxtime/yearsec
!!$          vs_km = vs_cm/1d5
!!$       endif
!!$    else
!!$       vs_cm   = 1.d5*vs_km
!!$       maxdist = maxtime*vs_cm*yearsec
!!$    endif
    bfield  = bbeta*sqrt(nh_init)*1.0d-6
!!$    grad_v  = 1.0e-2_dp*vs_cm/xll
    if (shock_type == "J") then
      viscosity = .true.
    endif
    !
    call check_parameters
    if (verbose) write(*,*) "I-read_parameters: done"
  end subroutine read_parameters

  subroutine check_parameters
    use module_grains, only    : tgrain,grmin,grmax,rho_grain
    use module_constants, only : yearsec
    use module_tools, only     : get_file_number
    use module_parameters_flags
    use module_parameters_output
    !
    logical :: error
    error = .false.
    if (iforh2.eq.4.and..not.do_h2) then
       write(*,*) "Inconsistent parameters iforh2 and do_h2"
       error = .true.
    endif
    if (error) stop
  end subroutine check_parameters

end module module_phys_var
