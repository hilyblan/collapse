module module_energetics
  use module_precision
  implicit none
  !___________________________________________________________________
  !
  ! The module 'MODULE_ENERGETICS contains variables and one
  ! subroutines related to the fluxes of energy, mass and momentum.
  ! All the variables are calculated here.
  !
  ! 'energetic_fluxes' is needed only for shock models.
  !___________________________________________________________________
  !
  real(kind=dp) :: mass_flux            ! mass flux (g/s/cm2)
  real(kind=dp) :: mass_flux_init       ! initial mass flux (g/s/cm2)
  real(kind=dp) :: mass_cons            ! relative mass flux conservation
  real(kind=dp) :: momentum_flux_kin    ! kinetic momentum flux (erg/cm3)
  real(kind=dp) :: momentum_flux_the    ! thermal momentum flux (erg/cm3)
  real(kind=dp) :: momentum_flux_mag    ! magnetic momentum flux (erg/cm3)
  real(kind=dp) :: momentum_flux_vis    ! viscous momentum flux (erg/cm3)
  real(kind=dp) :: momentum_flux        ! total momentum flux (erg/cm3)
  real(kind=dp) :: momentum_flux_init   ! initial total momentum flux (erg/cm3)
  real(kind=dp) :: momentum_cons        ! relative momentum flux conservation
  real(kind=dp) :: energy_flux_kin      ! kinetic energy_flux (erg/s/cm2)
  real(kind=dp) :: energy_flux_the      ! thermal energy_flux (erg/s/cm2)
  real(kind=dp) :: energy_flux_int      ! h2 internal energy_flux (erg/s/cm2)
  real(kind=dp) :: energy_flux_mag      ! magnetic energy_flux (erg/s/cm2)
  real(kind=dp) :: energy_flux_vis      ! viscous energy_flux (erg/s/cm2)
  real(kind=dp) :: energy_flux          ! total energy_flux (erg/s/cm2)
  real(kind=dp) :: energy_flux_init     ! initial total energy_flux (erg/s/cm2)
  real(kind=dp) :: energy_cons          ! relative energy flux conservation
  real(kind=dp) :: heating_old = 0.0_dp ! heating (erg/s/cm3) before last call to drive
  real(kind=dp) :: energy_gain = 0.0_dp ! energy gain (erg/s/cm2) of the gas due to heating
!!$  !
!!$  ! variables from "precision.f90" are private
!!$  private :: dp, long, minus_infinity, plus_infinity

contains

  subroutine energetic_fluxes
    use module_phys_var
    use module_evolution, only : bn, bi, bneg, h2_int_energy
    use module_constants, only : kb, pi, mp, everg
    use module_chemical_species, only : speci, ind_h, ind_he,nions,nneg
    use module_debug
    use module_parameters_flags
    !___________________________________________________________________
    !
    ! Calculate the fluxes of mass(g/s/cm2), momentum (erg/cm3) and
    ! energy (erg/s/cm2). When this subroutine is called for the first
    ! time, it save the initial fluxes for comparison purpose.  It
    ! also check the conservation of those fluxes.
    !___________________________________________________________________
    !
    implicit none
    integer(kind=long), save :: ncall = 0 ! counts the number of call to this routine
    !
    !
    ncall = ncall + 1
    if (.not.do_shock) return

    !cccccccccccccccccccccccccccccccccccccccc
    ! rajouter force visqueuse pour choc J
    !cccccccccccccccccccccccccccccccccccccccc
    !--- mass (g/s/cm2) ---
    !--- momentum (erg/cm3) ---
    !--- energy flux (erg/s/cm2) ---
    mass_flux         = rhon*vn
    momentum_flux_kin = rhon*vn**2
    momentum_flux_the = rhon*kb*tn/mun
    momentum_flux_mag = (bfield * vs_cm/vi)**2/(8._dp*pi)
    energy_flux_kin = 0.5_dp*rhon*vn**3
    energy_flux_the = 2.5_dp*rhon*vn*kb*tn/mun
    energy_flux_int = h2_int_energy*vn*kb
    energy_flux_mag = bfield**2/(4._dp*pi)*vs_cm**2/vi
    if (nions.gt.0) then
       mass_flux = mass_flux+rhoi*vi
       momentum_flux_kin = momentum_flux_kin+rhoi*vi**2
       momentum_flux_the = momentum_flux_the+rhoi*kb*ti/mui
       energy_flux_kin = energy_flux_kin+0.5_dp*rhoi*vi**3
       energy_flux_the = energy_flux_the+2.5_dp*rhoi*vi*kb*ti/mui
    endif
    if (nneg.gt.0) then
       mass_flux = mass_flux+rhoneg*vi
       momentum_flux_kin = momentum_flux_kin+rhoneg*vi**2
       momentum_flux_the = momentum_flux_the+rhoneg*kb*te/muneg
       energy_flux_kin = energy_flux_kin+0.5_dp*rhoneg*vi**3
       energy_flux_the = energy_flux_the+2.5_dp*rhoneg*vi*kb*te/muneg
    endif
    if ((shock_type == "J") .AND. (viscosity)) then
       Momentum_flux_vis = &
            (RhoN + RhoI + RhoNEG) * (XLL * grad_V)**2.0_DP
    else
       Momentum_flux_vis = 0.0_dp
    endif
    Momentum_flux = &
         Momentum_flux_kin + &
         Momentum_flux_the + &
         Momentum_flux_mag + &
         Momentum_flux_vis
    if ((shock_type == "J") .AND. (viscosity)) then
       Energy_flux_vis = &
            (RhoN + RhoI + RhoNEG) * Vn * (XLL * grad_V)**2.0_DP
    else
       Energy_flux_vis = 0.0_dp
    endif
    Energy_flux = &
         Energy_flux_kin + &
         Energy_flux_the + &
         Energy_flux_int + &
         Energy_flux_mag + &
         Energy_flux_vis
    if (debug) then
       write(*,*) "D-energetic_fluxes: Momentum_flux_kin",Momentum_flux_kin
       write(*,*) "D-energetic_fluxes: Momentum_flux_the",Momentum_flux_the
       write(*,*) "D-energetic_fluxes: Momentum_flux_mag",Momentum_flux_mag
       write(*,*) "D-energetic_fluxes: Momentum_flux_vis",Momentum_flux_vis
       write(*,*) "D-energetic_fluxes: RhoN   , kB , Tn , muN  ", RhoN   , kB , Tn , muN  
       write(*,*) "D-energetic_fluxes: RhoI   , kB , Ti , muI  ", RhoI   , kB , Ti , muI  
       write(*,*) "D-energetic_fluxes: RhoNEG , kB , Te , muNEG", RhoNEG , kB , Te , muNEG
    endif
    !
    ! if it's the first time we call this routine, then save the fluxes
    if (ncall <= 2) then
       mass_flux_init     = mass_flux
       momentum_flux_init = momentum_flux
       energy_flux_init   = energy_flux
    end if
    !
    ! check conservation
    Mass_cons = (Mass_flux_init - Mass_flux) / Mass_flux_init
    Momentum_cons = (Momentum_flux_init - Momentum_flux) / Momentum_flux_init
    Energy_gain = Energy_gain + &
         0.5_DP * dist_step * (Heating_old + (Bn + Bi + Bneg))
    Heating_old = Bn + Bi + Bneg
    Energy_cons = (Energy_flux_init - Energy_flux + Energy_gain) &
         / Energy_flux_init
  end subroutine energetic_fluxes
end module module_energetics
