module module_initialize
  use module_parameters_flags
  !___________________________________________________________________
  !
  ! Initialize
  ! * the shock parameters
  ! * the chemical species
  ! * the chemical reactions
  ! * the H2 molecule (levels, collision rates, lines)
  ! * the SiO molecule (collision rates, Einstein coefficients)
  !___________________________________________________________________
  !
  implicit none
  !
contains
  !
  subroutine initialize
    use module_chem_react
    use module_phys_var
    use module_grains
    use module_constants
    use module_h2!, only:read_h2_levels,initialize_rovib_h2,read_h2_rates,read_h2_lines
    use module_sio, only : read_sio_rates, einstein_coeff_sio
    use module_debug
    use module_line_excit
    use module_read_fe_data
    !___________________________________________________________________
    !
    ! Read input parameters, chemical species and reactions,
    ! Initializes all physical variables
    !___________________________________________________________________
    !
    implicit none
    character(len=*), parameter :: species_dat='output/species.dat'
    integer :: file_dat
    integer(kind=long) :: i
    real(kind=dp) :: polari,fac,mcore,mmantle
    !
    ! Input parameteres
    call read_parameters
    !
    ! Chemical species
    ! Elemental abundances from Ander & Grevesse 1989 ---
    call initialize_elements
    ! read in file_input + add 5 species + initialize indexes ...
    call read_species
    ! check the set of species
    call check_species
    !
    ! Write chemical species in formatted output for bookkeeping
    ! purposes
    file_dat = get_file_number()
    open(file_dat,file=species_dat,status='unknown',access='sequential',form='formatted',action='write')
    do i=1,nspec+1
       write(file_dat,'(i4,a10,es12.3,f12.3)') speci(i)%index,speci(i)%name,speci(i)%density,speci(i)%enthalpy
    enddo
    !
    ! Other variables
    ! number of species in each fluid, indexes of beginning and end
    do i=1,nspec
       if (index(speci(i)%name,'+') >0) then
          nions=nions+1
          e_ion=i
       else if (index(speci(i)%name,'-') >0) then
          nneg=nneg+1
          e_neg=i
       else if (index(speci(i)%name,'**') >0) then
          noncores=noncores+1
          e_cor=i
       else if (index(speci(i)%name,'*') >0) then
          nongrains=nongrains+1
          e_gra=i
       else
          nneutrals=nneutrals+1
          e_neu=i
       endif
       ! verification : we must have speci(i)%index=i
       if (speci(i)%index /= i) then
          print *, speci(i)%name, "  i =", i, speci(i)%index
          stop "*** Error, problem in species indexes"
       end if
    end do
    !
    ! first and last indices for each type of species
    b_neu = e_neu-Nneutrals+1
    b_ion = e_ion-Nions+1
    b_neg = e_neg-Nneg+1
    b_gra = e_gra-Nongrains+1
    b_cor = e_cor-Noncores+1
    if (debug) then
       write(*,*) "first and last indices for each type of species: "
       write(*,'(a,(10i6))') &
            "b_neu,e_neu,b_ion,e_ion,b_neg,e_neg,b_gra,e_gra,b_cor,e_cor",&
            b_neu,e_neu,b_ion,e_ion,b_neg,e_neg,b_gra,e_gra,b_cor,e_cor
       write(*,'(a,(10i6))') "Nneutrals,Nions,Nneg,Nongrains,Noncores",&
            Nneutrals,Nions,Nneg,Nongrains,Noncores
    endif
    !
    ! use the initial value of nH, which is read in MODULE_PHYS_VAR
    nh = nh_init
    !
    ! grains : mean square radius (cm2), density (cm-3)
    ! Assume dn(a)/da = N0 a^-3.5
    ! rgrain2 = <a^2>
    ! rgrain3 = <a^3>
    ! rgrain1 = <a>
    fac = grmin**(-2.5)-grmax**(-2.5)
    rgrain1    = 5./3.*(grmin**(-1.5)-grmax**(-1.5))/fac
    rgrain2    = 5.*(grmin**(-0.5)-grmax**(-0.5))/fac
    rgrain3    = 5.*(grmax**0.5 -grmin**0.5)/fac
    r_gr_scal23= rgrain2 * rgrain3**(-2./3.0_dp)
    r_gr_scal12= rgrain1/sqrt(rgrain2)
    r_grain = sqrt(rgrain2)
    !
    ! dens_grain_init: mass density of solid material only --> dens_grain
    mcore    = dot_product(dble(speci(b_cor:e_cor)%density),dble(speci(b_cor:e_cor)%mass))
    mmantle  = dot_product(dble(speci(b_gra:e_gra)%density),dble(speci(b_gra:e_gra)%mass))
    md_grain = mcore+mmantle ! Total (core + mantle) mass of grain per unit volume of gas
    dens_grain_init          = mcore*3./(4.*pi*rho_grain*rgrain3)
    speci(ind_grain)%density = dens_grain_init
    dens_grain               = dens_grain_init
    ratio_grain_gas          = mcore/nh/mp
    !
    ! Update g0, g+, g- number densities (cm-3)
    ! Arbitrarily evenly distributed amongst the three species
    speci(ind_gg0)%density     = dens_grain/3.
    speci(ind_ggplus)%density  = dens_grain/3.
    speci(ind_ggminus)%density = dens_grain-speci(ind_ggplus)%density-speci(ind_gg0)%density
    if (debug) then
       write(*,*) "D-initialize: grain related quantities"
       write(*,*) "grmin, grmax rho_grain ratio_grain_gas r_gr_scal12 r_gr_scal23",&
            "MD_grain R_grain Dens_GRAIN rgrain2 rgrain3 r_grain:"
       write(*,*) grmin,grmax,rho_grain,ratio_grain_gas,&
            r_gr_scal12,r_gr_scal23,&
            md_grain,r_grain,dens_grain,rgrain2,rgrain3,r_grain
    endif
    !
    ! number density (cm-3) of each fluid =sum(density)
    fac = (sum(dble(speci(b_gra:e_gra)%density))+sum(dble(speci(b_cor:e_cor)%density)))
    densityn   = sum(dble(speci(b_neu:e_neu)%density)) + fac/3.
    densityi   = sum(dble(speci(b_ion:e_ion)%density)) + fac/3. ! ions >0
    densityneg = sum(dble(speci(b_neg:e_neg)%density)) + fac/3. ! ions <0
    speci(ind_e)%density = densityi - densityneg                 ! electrons (charge neutrality)
    !
    ! mass density (g.cm-3) of each fluid =sum(mass*density)
    ! include the contributions of the grains
    rhon   = dot_product(dble(speci(b_neu:e_neu)%density),dble(speci(b_neu:e_neu)%mass)) + md_grain/3.
    rhoi   = dot_product(dble(speci(b_ion:e_ion)%density),dble(speci(b_ion:e_ion)%mass)) + md_grain/3.
    rhoneg = dot_product(dble(speci(b_neg:e_neg)%density),dble(speci(b_neg:e_neg)%mass)) + md_grain/3.
    !
    ! mean mass (g) of each fluid mu=rho/density
    mun   = rhon/densityn
    mui   = rhoi/densityi
    muneg = rhoneg/densityneg
    ! !>>>>>> PHB
    ! if (nions.gt.0) epsion=1.d0
    ! if (nneg.gt.0) epsneg=1.d0
    ! mui   = epsion*rhoi/densityi
    ! muneg = epsneg*rhoneg/densityneg
    ! !<<<<<< PHB

    ! temperatures : all the same at the beginning of the shock
    Ti=Tn
    Te=Tn
    !
    ! Initial value of the fractional radius
    xsphere = 1._dp - 1.d-4
    !
    ! Free-fall time
    tau_ff = sqrt(3.d0*pi/(32.*gr_cons*(rhon+rhoi)))
    !
    ! Ion-neutral ambipolar diffusion timescale
    ! [equ. (5) of Walmsley et al. 2004, A&A, 418, 1035]
    ! PHB (11jul2014): introduce 'polari'
    polari = speci(ind_h)%density*alpha_h+(speci(ind_ph2)%density+speci(ind_oh2)%density)*alpha_h2
    polari = polari/(speci(ind_h)%density+speci(ind_ph2)%density+speci(ind_oh2)%density)
    tau_i_n = (2./pi/gr_cons) * MAX(2.41_dp*pi*qe*sqrt((mun+mui) &
         *polari/mun/mui), 1.0d-15*abs_deltav) &
         * rhoi / (mun+mui) / rhon
    !
    !  grain_n ambipolar diffusion timescale 
    tau_grain_n = (2.0/pi/gr_cons)*pi*rgrain2*sqrt(8.0_dp*kb*tn/(pi*mun)) & 
         * dens_grain / rhon &
         * (speci(ind_ggminus)%density + speci(ind_ggplus)%density)/ &
         (speci(ind_ggminus)%density + speci(ind_ggplus)%density + speci(ind_gg0)%density)

!!$    !  the collapse time is taken equal to the greatest of the local free fall time
!!$    !  and the two ambipolar diffusion time scales
!!$    tau_collapse = MAX(tf,tau_grain_n)
!!$    tau_collapse = MAX(tau_collapse,tau_i_n)
    !
    ! The collapse timescale is the free fall timescale multiplied by
    ! a scaling factor 'tff_mult' (input parameter)
    tau_collapse = tau_ff*tff_mult
    !
    ! In the case of collapse calculation, adjust the velocity to the
    ! user-defined sphere radius when supplied.
    ! maxtime = min(maxtime,tau_collapse/yearsec)
    if (.not.do_shock) then
       if (masscloud.gt.0) then
          maxdist = (3.d0*masscloud/(4.d0*pi*(rhon+rhoi)))**(1./3.) ! [cm]
          vs_cm = 1.0 !maxdist/maxtime/yearsec
          vs_km = vs_cm/1d5
       else
          vs_cm   = 1.d5*vs_km
          maxdist = maxtime*vs_cm*yearsec
       endif
       ! if (maxdist.gt.0) then
       !    vs_cm = maxdist/maxtime/yearsec
       !    vs_km = vs_cm/1d5
       ! else
       !    vs_cm   = 1.d5*vs_km
       !    maxdist = maxtime*vs_cm*yearsec
       ! endif
    endif
    !
    ! Dynamical properties
    ! * Velocity of neutrals and ions (cm/s)
    ! * Velocity gradient (km.s-1.cm-1)
    grad_v = 1.0e-2_dp*vs_cm/xll
    vgrad  = vs_km/tout
    vn=vs_cm
    if (nfluids == 1) then
       deltav=0.0_dp
    else
       deltav=deltavmin ! = vi-vn, to start the shock
    endif
    if (shock_type == 'S') then
       deltav = deltav * 10.0d3
    endif
    abs_deltav=abs(deltav)
    vi=vn-deltav
    !
    ! Read chemical network
    ! 1/ read reactions
    ! 2/ determine reaction types
    ! 3/ check reactions
    ! 4/ add missing reverse reactions for endothermic reactions
    ! 5/ Write the chemical network in machine readable format
    ! 6/ Write the chemical network for spin symetry conservation purposes
    call read_reactions
    call reaction_type
    call check_reactions
    call write_chemical_network
    call write_branching_network
    if (do_we_add_reactions) call add_reverse_reactions
    !
    ! Initialize spectroscopic data for H2, SiO, and Fe
    call h2_init
    if (do_thermal) then
       call read_sio_rates      ! collision rates for sio-h2
       call einstein_coeff_sio  ! Einstein coefficients for SiO
       call read_fe_data
       call line_excit_init
    endif
    !
    if (verbose) write(*,*) "I-initialize: done"
  end subroutine initialize
  

  subroutine initialize_indices
    use module_phys_var
    use module_chemical_species
    use module_debug
    !_________________________________________________________________
    !
    ! Initialize usefule indices
    !_________________________________________________________________
    !
    implicit none
    !
    nv_mhd = 13
    dimtot = nv_mhd + nspec + nh2_lev
    !
    ! Chemical species
    ! * index of neutrals, species on grain mantles, species on cores,
    ! * positive ions, negative ions
    bv_speci = nv_mhd + 1
    ev_speci = bv_speci + nspec - 1
    bv_neu = b_neu + bv_speci - 1
    ev_neu = e_neu + bv_speci - 1
    bv_gra = b_gra + bv_speci - 1
    ev_gra = e_gra + bv_speci - 1
    bv_cor = b_cor + bv_speci - 1
    ev_cor = e_cor + bv_speci - 1
    bv_ion = b_ion + bv_speci - 1
    ev_ion = e_ion + bv_speci - 1
    bv_neg = b_neg + bv_speci - 1
    ev_neg = e_neg + bv_speci - 1
    !
    ! H2 levels
    bv_h2_lev = nv_mhd + nspec + 1
    ev_h2_lev = bv_h2_lev + nh2_lev - 1
    !
    if (debug) then
       write(*,*) "D-initialize_indices:"
       write(*,*) "nv_mhd,nspec: ",nv_mhd,nspec
       print *, " bv_neu    =", bv_neu, ", ev_neu    =", ev_neu
       print *, " bv_gra    =", bv_gra, ", ev_gra    =", ev_gra
       print *, " bv_cor    =", bv_cor, ", ev_cor    =", ev_cor
       print *, " bv_ion    =", bv_ion, ", ev_ion    =", ev_ion
       print *, " bv_neg    =", bv_neg, ", ev_neg    =", ev_neg
    endif
  end subroutine initialize_indices
  !
  !
  subroutine initialize_arrays
    use module_phys_var
    use module_chemical_species
    use module_h2
    use module_grains
    use module_debug
    use module_constants
    !_________________________________________________________________
    !
    ! Initialize global arrays
    ! * memory allocation
    ! * initial values
    ! * physical variable indices
    !_________________________________________________________________
    !
    implicit none
    integer :: i
    !
    allocate(yarray(1:dimtot)) ! physical variables (y)
    allocate(yn    (0:nspec))  ! change in number density (cm-3.s-1), calculated in chemistry
    !
    yarray = 0.d0
    yn     = 0.d0
    !
    ! MHD variables
    i = 1  ; yarray(i) = xsphere     ; iv_xsphere = i     ! fractional radius of collapsing sphere
    i = 2  ; yarray(i) = vn          ; iv_vn = i          ! neutrals (cm/s)
    i = 3  ; yarray(i) = vi          ; iv_vi = i          ! ions (cm/s)
    i = 4  ; yarray(i) = rhon        ; iv_rhon = i        ! neutrals (g.cm-3)
    i = 5  ; yarray(i) = rhoi        ; iv_rhoi = i        ! ions (g.cm-3)
    i = 6  ; yarray(i) = tn          ; iv_tn = i          ! neutrals (k)
    i = 7  ; yarray(i) = ti          ; iv_ti = i          ! ions (k)
    i = 8  ; yarray(i) = te          ; iv_te = i          ! electrons (k)
    i = 9  ; yarray(i) = densityn    ; iv_densityn = i    ! neutrals (cm-3)
    i = 10 ; yarray(i) = densityi    ; iv_densityi = i    ! ions (cm-3)
    i = 11 ; yarray(i) = rhoneg      ; iv_rhoneg = i      ! negative ions (g.cm-3)
    i = 12 ; yarray(i) = grad_v      ; iv_gv     = i      ! - neutral velocity gradient
    i = 13 ; yarray(i) = dens_grain  ; iv_dens_grain = i  ! grain number density
    !
    ! Density of species (cm-3) Does not include special species
    yarray(bv_speci:ev_speci) = speci(1:nspec)%density
    !
    ! Density of H2 levels (cm-3)
    if (do_h2) yarray(bv_h2_lev:ev_h2_lev) = h2_lev(1:nh2_lev)%density
    !
    if (debug) then
       write(100,'(a15,i5)')"iv_xsphere   ", iv_xsphere     ! fractional radius of collapsing sphere
       write(100,'(a15,i5)')"iv_vi        ", iv_vi          ! ions (cm/s)
       write(100,'(a15,i5)')"iv_rhon      ", iv_rhon        ! neutrals (g.cm-3)
       write(100,'(a15,i5)')"iv_rhoi      ", iv_rhoi        ! ions (g.cm-3)
       write(100,'(a15,i5)')"iv_tn        ", iv_tn          ! neutrals (k)
       write(100,'(a15,i5)')"iv_ti        ", iv_ti          ! ions (k)
       write(100,'(a15,i5)')"iv_te        ", iv_te          ! electrons (k)
       write(100,'(a15,i5)')"iv_densityi  ", iv_densityi    ! ions (cm-3)
       write(100,'(a15,i5)')"iv_rhoneg    ", iv_rhoneg      ! negative ions (g.cm-3)
       write(100,'(a15,i5)')"iv_gv        ", iv_gv          ! - neutral velocity gradient
       write(100,'(a15,i5)')"iv_dens_grain", iv_dens_grain  ! grain number density
       close(100)
    endif
  end subroutine initialize_arrays

end module module_initialize
