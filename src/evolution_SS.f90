module module_evolution
  use module_precision
  use module_parameters_flags
  !___________________________________________________________________
  !
  ! * Source terms for the chemistry
  ! * Partial derivatives for DVODE integration package.
  !___________________________________________________________________
  !
  implicit none
  !
  ! source terms for each fluid (N -> neutrals, I -> ions, NEG ->
  ! negative ions) = rates of change in number, mass, momentum and
  ! energy densities
  real(kind=dp) :: cupn, cupi, cupneg        ! increase in number density (cm-3.s-1)
  real(kind=dp) :: cdownn, cdowni, cdownneg  ! decrease in number density (cm-3.s-1)
  real(kind=dp) :: ynn, yni, ynneg           ! change in number density (cm-3.s-1)
  real(kind=dp) :: sn, si, sneg, score       ! change in mass density (g.cm-3.s-1)
  real(kind=dp) :: an, ai, aneg              ! change in momentum density (g.cm-2.s-2)
  real(kind=dp) :: bn, bi, bneg              ! change in energy density (erg.cm-3.s-1)
  real(kind=dp) :: h2_int_energy             ! internal energy of h2 (K cm-3)
  !
  !
  private :: evolution_chemistry, evolution_physics
  !
contains
  !
  subroutine evolution_chemistry
    use module_chem_react
    use module_grains
    use h2_variables
    use module_h2, only: sel_ch_h2,sel_ne_h2,sel_io_h2,sel_el_h2,sel_rx_h2,h2_lev
    use module_phys_var
    use module_constants, only : pi, kb, me, everg, zero
    use module_debug
    !_________________________________________________________________
    !
    ! Calculates source terms related to chemistry. Each type of
    ! reaction gets its own reaction rate (with different meanings and
    ! unities for gamma, alpha and beta ).
    !
    ! For each specy :
    ! * YN : rate at which it is formed through chemical reactions [1], eqs. (14) and (20)
    ! * YS : rate of change in mass (g.cm-3.s-1)
    ! * YA : rate of change in momentum (g.cm-2.s-2)
    ! * YB : rate of change in energy (erg.cm-3.s-1)
    ! For each fluid (x= n -> neutrals, i -> ions, neg -> negative ions and e-)
    ! * CupX   : rate of increase in number density (cm-3.s-1)
    ! * CdownX : rate of decrease in number density (cm-3.s-1)
    ! * Nx     : rate of change in number density (cm-3.s-1)
    ! * Sx     : rate of change in mass density (g.cm-3.s-1)
    ! * Ax     : rate of change in momentum density (g.cm-2.s-2)
    ! * Bx     : rate of change in energy density (erg.cm-3.s-1)
    !
    ! References :
    !     [1] Flower et al., 1985, MNRAS 216, 775
    !     [2] Flower et al., 1986, MNRAS 218, 729
    !     [3] Pineau des Forets et al., 1986, MNRAS 220, 801
    !     [4] Viala et al., 1988, A&A 190, 215
    !_________________________________________________________________
    !
    implicit none
    !
    ! Local variables
    real(kind=dp), dimension(0:nspec_plus) :: up_n, up_mv, up_mv2, up_de
    real(kind=dp), dimension(0:nspec_plus) :: down_n, down_mv, down_mv2
    real(kind=dp), dimension(0:nspec) :: ys, ya, yb
    integer(kind=long) :: i, ir1, ir2, ip1, ip2, ip3, ip4
    real(kind=dp) :: massr1, massr2, densityr1, densityr2, v_r1, v_r2, mass_prod, nprod_m1
    real(kind=dp) :: vcm, vcm2
    real(kind=dp) :: creation_rate, destr_r1, destr_r2
    real(kind=dp) :: creation_vcm, creation_vcm2, creation_de
    real(kind=dp) :: destr1_vcm, destr2_vcm, destr1_vcm2, destr2_vcm2
    real(kind=dp) :: exp_factor, exp_factor1, exp_factor2, coeff, tr, ts, energy, teff
    real(kind=dp) :: react_rate, r_alph, r_beta, r_gamm
    real(kind=dp) :: frach2,nfree,nhtot
    real(kind=dp) :: ccf, cfact,fac,vdrift
    real(kind=dp) :: stick_h_gr = 1.0_dp
    integer :: lev
    !
    ! Sanity checks
    if (debug) then
       write(*,*) "D-evolution_chemistry: grain related quantities"
       write(*,*) "grmin, grmax rho_grain ratio_grain_gas r_gr_scal12",&
            "r_gr_scal23 MD_grain R_grain Dens_GRAIN rgrain2 rgrain3 r_grain:"
       write(*,*) grmin,grmax,rho_grain,ratio_grain_gas,&
            r_gr_scal12,r_gr_scal23,&
            md_grain,r_grain,dens_grain,rgrain2,rgrain3,r_grain
       write(*,'("h2_fo:",2a12,20a10)') "SP1","SP2","Tgrain","STICK_H","R^2",&
            "mass(R2)","C_Rate","Rate","Dens_gr","n(SP1)","n(SP2)","nH_0",&
            "r_gr_h2","n(H2)","n(pH2)","n(oH2)"
    endif
    !
    ! initialization
    creation_rate = 0.d0
    up_n     = zero
    up_mv    = zero
    up_mv2   = zero
    up_de    = zero
    down_n   = zero
    down_mv  = zero
    down_mv2 = zero
    if (do_h2) then
       sel_ch_h2(1:nh2_lev) = zero
       sel_ne_h2(1:nh2_lev) = zero
       sel_io_h2(1:nh2_lev) = zero
       sel_el_h2(1:nh2_lev) = zero
    endif
    !
    ! Density
    nfree = dens_h+dens_h2       ! number of free particles (cm-3)
    nhtot = dens_h+2.d0*dens_h2  ! proton density (cm-3)
    !
    ! Correction for the variation of the square of the grain radius
    ! relative to the mean square radius for the MRN distribution
    cfact = Rgrain2 / 4.088d-12
    vdrift = vi-vn
    !
    ! calculates creation and destruction terms
    ! for each reaction (i=index of the reaction)
    do i=1,nreact
       !
       ! Index, mass, abundance, velocity of each reactant and product
       ir1       = react(i)%r(1)       ; ir2       = react(i)%r(2)
       ip1       = react(i)%p(1)       ; ip2       = react(i)%p(2)
       ip3       = react(i)%p(3)       ; ip4       = react(i)%p(4)
       massr1    = speci(ir1)%mass     ; massr2    = speci(ir2)%mass
       densityr1 = speci(ir1)%density  ; densityr2 = speci(ir2)%density
       v_r1      = speci(ir1)%velocity ; v_r2      = speci(ir2)%velocity
       mass_prod = react(i)%mass_prod  ! sum of the masses of the products
       nprod_m1  = react(i)%nprod_m1   ! number of products minus one
       r_alph    = react(i)%alpha      ! 
       r_beta    = react(i)%beta       ! 
       r_gamm    = react(i)%gamma      ! 
       !
       ! Center of mass velocity (cm/s), and its square
       Vcm  = (massR1*V_R1 + massR2*V_R2) / (massR1+massR2)
       Vcm2 = Vcm*Vcm
       !
       ! now calculate the reaction rate (cm3.s-1) according
       ! to the type of the reaction (effective temperature,
       ! significance of gamma, alpha, beta, ...).
       select case (react(i)%type)
       CASE ('PHOTO')
          ! Photo reactions  (cm3.s-1)
          exp_factor = exp(-min(r_beta*av,maxexp))
          react_rate    = r_gamm*exp_factor*rad
          creation_rate = react_rate*densityr1*densityr2
       CASE ('CRIO_DIR','CRIO_IND')
          ! cosmic ray ionization or dissociation (cm3.s-1)
          ! We take into account :
          !     0) direct cosmic ray induced ionization/dissociation
          !     1) the CR induced secondary photons
          !     2) the UV photons produced by collisional excitation of H2 by electrons
          !        following shock heating.
          !     3) allow for the differential compression of the charged grains and the
          !        neutrals (Vi/Vn). This factor is applicable to cases (1) and
          !        (2) above, and not to case (0). However, the direct CR induced
          !        processes are completely negligible within the shock, wich is
          !        where Vi/Vn differs from 1.
          ! note that react(i)%beta is preset at a very large value (in REACTION_TYPE)
          ! for direct cosmic ray ionization : react(i)%beta = 1.0D8

          ! JLB + GP : 13 juin 2002  - correction de la fraction de H2 dans le gaz
          ! JLB     exp_factor = MIN(r_beta/Te,180._DP)
          !
          ! >> OLD CR
          exp_factor = r_beta / Te
          coeff      = sqrt(8.0_dp*kb*te/(pi*me)) * 1.0d-16 * (4.29d-1 + 6.14d-6*te)
          if (react(i)%type.eq.'CRIO_IND') then
             !
             ! DRF+GDF: Jan-2017
             ! proportionality factor, fracH2: 0.46 is the probability of cosmic ray ionization of H,
             ! 1.0 is the probability of cosmic ray ionization of H2
             ! The secondary photons are produced by radiative cascade, following collisional 
             ! excitation of Rydberg states of H2 by secondary electrons, generated by cosmic rays
             !
             ! Note: beta=140e3 K in SECPHO reactions corresponds to
             ! the energy of a 1000 ang photon capable of exciting H2
             ! into an electronic state;
             !
             ! Note: the correction coefficient 'coeff' is used to
             ! take into account electrons produced in strong J shocks
             ! which can also excite H2 into electronic states; this
             ! correction is, however, only significant at
             ! temperatures above ~1000 K.
             !
             frach2 = dens_h2/nhtot
             ! frach2 = (0.46*dens_h+dens_h2)/nfree*dens_h2/nfree
             ! if (te.gt.777) then
             !    exp_factor = 0.d0 !exp(-r_beta/te) !,maxexp))
             !    coeff = 0.d0
             !    ! exp_factor = exp(-min(r_beta/te,maxexp))
             !    ! coeff = sqrt(8.0_dp*kb*te/(pi*me))*1.0d-16*(4.29d-1 + 6.14d-6*te)
             !    ! coeff = coeff*exp_factor
             ! else
             !    coeff = 0.d0
             ! endif
          else if (react(i)%type.eq.'CRIO_DIR') then
             ! Direct cosmic ray induced ionization/dissociation
             frach2 = 1.0_dp
             coeff  = 0.d0
          endif
          !
          ! DRF+GDF: Jan-2017
          ! For the treatment of the grains (G, G-, G+), see F & PdesF 2003, MNRAS, 343, 390, Sect. 2.2.2
          ! Rate per unit volume of detachment of electrons from grains by the H_2 fluorescence photons (SECPHO):
          ! 0.15 zeta n_H f y (cm-3 s-1)
          ! where
          ! -- 0.15 is the fraction of the secondary electrons that collisionally excite the Rydberg states of H_2, 
          !    generating the secondary photons in the subsequent radiative cascade;
          ! -- zeta (s-1) is the rate of cosmic ray ionization of H_2;
          ! -- n_H = (n_H / n_g) n_g 
          ! -- f depends on the supposed fractional PAH abundance: f = 0.87 for n_PAH / n_H = 1e-6;
          !    f = 0.99 for n_PAH / n_H = 1e-7
          ! -- y = 0.03 for ionization of neutral grains, y = 0.2 for electron photo-detachment from negatively 
          !    charged grains (both values are uncertain).
          ! In the above paper, 
          ! n_H / n_g = 1 / 7.1e-11 = 1.4e10
          ! and hence the rates per unit volume of ionization/detachment of electrons from grains 
          ! are, when n_PAH / n_H = 0.0,
          ! zeta * n_g * 6.3e7 cm-3 s-1 for neutral grains
          ! and
          ! zeta * n_g- * 4.2e8 cm-3 s-1 for negatively charged grains
          ! The values of r_gamm must be scaled to the current value of n_H / n_g and also adapted 
          ! to the value of n_PAH / n_H
          ! >> NEW CR
          ! if (ir2.eq.ind_SECPHO) then
          !    if ((ir1.eq.ind_ggminus).or.(ir1.eq.ind_gg0).or.(react(i)%p(1) == ind_grain)) then
          !       ! See above
          !       ! Grain ionization and photodesorption
          !       r_gamm = r_gamm*7.1d-11*nhtot/dens_grain
          !    else
          !       ! Following the analysis of Gredel et al. 1989, ApJ, 347,289, the probabilities, p_M,
          !       ! of photo-ionization and photo-dissociation of atoms and molecules are inversely
          !       ! proportional to sigma_g, the grain extinction cross section per hydrogen nucleus;
          !       ! sigma_g = <n_g sigma> / n_H cm^2
          !       ! Gredel et al. adopt a value of sigma_g = 2e-21 cm^2, Thus, r_gamm must be scaled by
          !       ! a factor 2.e-21 / (<n_g sigma> / n_H)
          !       r_gamm = r_gamm*2.d-21/(dens_grain*pi*rgrain2/nhtot)
          !    endif
          !    if (react(i)%p(1) == ind_grain) then
          !       ! nlayers_grain = number of layers in the mantles
          !       ! IF nlayers_grain > 1 => the detachment probability is independent of nlayers_grain
          !       ! IF nlayers_grain < 1 => the detachment probability is proportional to nlayers_grain
          !       Nlayers_grain = Dens_ongrains/Nsites_grain/Dens_GRAIN
          !       r_gamm = r_gamm * min(nlayers_grain,1.0_dp) * dens_grain / dens_ongrains
          !    endif
          ! endif
          ! react_rate    = zeta*(tn/300._dp)**r_alph*r_gamm + dens_e*coeff/0.15*r_gamm*vi/vn
          ! react_rate    = react_rate*frach2
          ! creation_rate = react_rate*densityr1*densityr2
          react_rate    = zeta*(tn/300._dp)**r_alph*r_gamm+dens_e*coeff*exp(-exp_factor)/0.15_dp*r_gamm*vi/vn
          react_rate    = react_rate*frach2
          creation_rate = react_rate*densityr1*densityr2
       CASE ('CR_DE')
          ! cosmic ray induced desorption from grains (s-1)
          ! 1. values of the rates of desorption, relative to that for CO, are obtained from:
          !         exp[-(E_ads - 855)/T_g(max)], where E_ads is the adsorption energy (K)
          !  of the species (Aikawa et al. 1996, ApJ, 467, 684, Table 1; see also Bergin et al.
          !  2002, ApJ, 570, L101) and 855 K is the adsorption energy of CO on CO ice
          !  T_g(max) = 70 K (Hasegawa & Herbst 1983, MNRAS, 261, 83, equ. (15))
          !
          ! 2. Normalize to the current value of the CR ionization rate, in units of 1e-17 s-1
          react_rate = r_gamm * dens_grain / dens_ongrains * pi * rgrain2
          react_rate = react_rate * exp(-(r_beta-855.d0)/70d0)
          react_rate = react_rate * zeta / 1.d-17 
          creation_rate = react_rate * densityr1 * densityr2
       CASE ('H2_FO')
          ! H2 and HD formation on grains
          ! * Andersson & Wannier for H adsorption
          ! * Hollenbach & McKee (Eq. 3.8, 1979, ApJ Suppl, 41, 555)
          ! react_rate (cm3.s-1) = r_gamm * nH/Dens_H * ((Tn/300._DP)**r_alph)
          teff_grain = massr1*vdrift**2/(3.0_dp*kb) + speci(ir1)%temperature
          stick_h_gr = 1.0_dp / sqrt(1.0_dp + teff_grain / 30.0_dp)
          stick_h_gr = 1.0_dp / (1.d0 + 0.04_dp*(teff_grain+tgrain)**0.5_dp &
               + 2.0d-3*teff_grain + 8.0d-6*teff_grain*teff_grain)
          !  in order to select the appropriate species mass when
          !  computing the collision frequency, the reactants must be:
          !  R1 = H and R2 = D, in that order, when computing the HD
          !  formation rate
          react_rate = stick_h_gr * pi * rgrain2 * sqrt(8.0_dp*kb*teff_grain/massr2/pi)
          creation_rate = react_rate*dens_grain*densityr2
          if (ir1.eq.ind_h.and.ir2.eq.ind_h) then
             ! This is formation of H2
             for_gr_h2 = creation_rate
             r_gr_h2   = creation_rate/densityr2/nh_init ! R-factor: formation of h2 on grains
          endif
          ! r_gamm: branching ratio for ortho and para formed on grains
          creation_rate = r_gamm * creation_rate
          if (debug) &
               write(*,'("h2_fo:",2(a12),20es10.2)') speci(ir1)%name,speci(ir2)%name,&
               teff_grain,stick_h_gr,rgrain2,massr2,&
               creation_rate,react_rate,dens_grain,densityr1,densityr2,nh_init,&
               r_gr_h2,dens_h2,dens_ph2,dens_oh2
       CASE ('THREE')
          ! Three-body reactions on grains surface
          ! Reaction rate (cm3.s-1) (Andersson & Wannier for H adsorption)
          ! Effective temperature
          teff_grain = massr1*vdrift**2/(3.0_dp*kb)+speci(ir1)%temperature
          react_rate = r_gamm*pi*rgrain2*sqrt(8.0_dp*kb*teff_grain/massr1/pi)
          react_rate = react_rate/(teff_grain/r_beta+1.0_dp)**2.0_dp
          react_rate = react_rate*dens_grain/dens_ongrains
          creation_rate = react_rate*densityR1*densityR2
       CASE ('SPUTT')
          ! Sputtering of grain mantles
          ! Reaction rate (cm3.s-1)
          ! reaction  : X* + IR2 -> X + IR2 + GRAIN
          ! Nlayers_grain = number of layers in the mantles
          ! IF Nlayers_grain > 1 => the detachment probability is
          !                         gamma * n(X*)/n(ongrains)
          ! IF Nlayers_grain < 1 => the detachment probability is
          !                         gamma * n(X*)/(Nsites_grain*n(grains))
          ! Effective temperature of grain sputtering: teff_grain
          ts = massr2*(v_r2-v_r1)*(v_r2-v_r1)/(3.0_dp*kb)
          tr =  speci(ir2)%temperature
          teff_grain = ts + tr
          nlayers_grain = dens_ongrains/nsites_grain/dens_grain
          exp_factor1   = (4.0_dp*r_beta-1.5_dp*ts)/tr
          exp_factor2   = 4.0_dp*r_beta/teff_grain
          exp_factor    = min(max(exp_factor1,exp_factor2),maxexp)
          exp_factor    = exp(-exp_factor)
          react_rate = 4.0_dp * r_gamm * dens_grain/dens_ongrains * &
               pi*rgrain2 * sqrt(8.0_dp*kb*teff_grain/massr2/pi)
          react_rate = react_rate * (1.0_dp+2.0_dp*teff_grain/(4.0_dp*r_beta)) * &
               exp_factor
          react_rate = react_rate*min(nlayers_grain,1.0_dp)
          creation_rate = react_rate*densityR1*densityR2
       CASE ('EROSI')
          ! erosion of grain cores
          ! Reaction rate (cm3.s-1)
          ! The rate coefficient is proportional to the grain density
          ! and not to the core species density
          ! Absolute value of the ion-neutral velocity drift
          react_rate = zero
          energy = massr2 * abs_deltav*abs_deltav / 2.0_dp / everg ! energy (ev)
          if (energy > r_alph) react_rate = pi*rgrain2*abs_deltav*r_gamm*exp(-r_beta/(energy-r_alph))
          if (densityr1 > dens_grain) then
             creation_rate = react_rate*dens_grain*densityr2
          else
             creation_rate = react_rate*densityr1*densityr2
          endif
       CASE ('ADSOR')
          ! adsorption on to grains
          ! reaction rate (cm3.s-1) / Ref: Tielens & Hollenbach
          ! here, gamma = sticking coefficient
          ! react_rate = STICK * SIGMA(GRAIN) * V(MOYENNE)
          ! effective temperature (K) 
          teff_grain    = massr1*vdrift**2/(3.0_dp*kb) + speci(ir1)%temperature
          react_rate = r_gamm * pi * rgrain2 *sqrt(8.0_dp*kb*teff_grain/massr1/pi)
          if(r_beta /=zero) then
             react_rate = react_rate / (1.0_dp + 0.04_dp*(teff_grain+tgrain)**0.5_dp + &
                  2.0d-3 * teff_grain + 8.0d-6*teff_grain*teff_grain)
             write(*,*) "PHB //",i,ir1,ir2,ip1,ip2,r_beta,teff_grain,tgrain,react_rate
          end if
          creation_rate = react_rate * densityR1 * densityR2
       CASE ('DISSO')
          ! Collisional dissociation of H2
          ! effective temperature of the reaction                               
          ! see [1], eq. (43) and [3], eqs.(2.1)-(2.3)                          
          ! for reactions between members of the same fluid : Ts = 0, Teff = Tr 
          ! destruction is computed state by state                             
          ! reduce dissociation energy r_beta by the excitation energy         
          ! of the level
          ts   = massr1*massr2*(v_r1-v_r2)*(v_r1-v_r2)/(massr1+massr2)/3.0_dp/kb
          tr   = (massr1*speci(ir2)%temperature + massr2*speci(ir1)%temperature)/(massr1+massr2)
          teff = ts+tr
          do lev=1,nh2_lev
             ! exponential factor ([3], page 810)
             exp_factor1 = (r_beta - h2_lev(lev)%energy - 3.0_dp*ts) / tr
             exp_factor2 = (r_beta - h2_lev(lev)%energy) / teff
             exp_factor  = max(exp_factor1,exp_factor2)
             exp_factor  = min(max(exp_factor1,0.0_dp),maxexp)
             sel_rx_h2(lev) = exp(-exp_factor)
          end do
          ! reaction rate (cm3.s-1) 
          sel_rx_h2 = r_gamm * sel_rx_h2 * ((teff/300._dp)**r_alph)
          ! note: we avoid dividing by n(h2) just to be able to remultiply...
          react_rate = sum(dble(h2_lev(1:nh2_lev)%density * sel_rx_h2))
          creation_rate = react_rate * densityr2
          sel_ch_h2 = sel_ch_h2 - sel_rx_h2 * densityr2
          if (ir2 >= b_neu .and. ir2 <= e_neu) then
             sel_ne_h2 = sel_ne_h2 - sel_rx_h2 * densityr2
          else if (ir2 >= b_ion .and. ir2 <= e_ion) then
             sel_io_h2 = sel_io_h2 - sel_rx_h2 * densityr2
          else if (ir2 == ind_e) then
             sel_el_h2 = sel_el_h2 - sel_rx_h2 * densityr2
          endif
       CASE ('OTHER', 'REVER','NEUTGR','CHARGR')
          ! all other reactions
          ! reaction rate (cm3.s-1) 
          ! (including reverse reactions added in ADD_INVERSE_REACTIONS)
          ! effective temperature of the reaction                               
          ! see [1], eq. (43) and [3], eqs.(2.1)-(2.3)                          
          ! for reactions between members of the same fluid : Ts = 0, Teff = Tr 
          ts = massr1*massr2*(v_r1-v_r2)*(v_r1-v_r2)/(massr1+massr2)/3.0_dp/kb
          tr = (massr1*speci(ir2)%temperature + massr2*speci(ir1)%temperature) / &
               (massr1+massr2)
          teff = ts + tr
          ! exponential factor ([3], PAGE 810) 
          ! [3] Pineau des Forets et al., 1986, MNRAS 220, 801
          exp_factor1 = (r_beta-3.0_dp*ts) / tr
          exp_factor2 = r_beta / teff
          exp_factor  = min(max(exp_factor1,exp_factor2),maxexp)
          exp_factor = exp(-exp_factor)
          react_rate = r_gamm* (teff/300.0)**r_alph*exp_factor
          !
          !  include temperature dependent Coulomb enhancement factor
          !  in the cases of neutralization of charged grains, G+ and
          !  G-, by electrons and positive ions, respectively the
          !  enhancement factor is taken form Draine & Sutin 1987
          !  (equ. 3.4) correct the reaction rate coefficient by
          !  square of the ratio of the mean grain radius, as computed
          !  using the MRN grain size distribution, with the min and
          !  max values of the radius taken to be 0.01 and 0.3 micron
          !  (which yields a mean radius of 0.02 micron), and as
          !  determined within the program, using the specified values
          !  of r_min, r_max
          !
          if (react(i)%type.eq.'NEUTGR') then
             ! Product is a neutral grain
             fac = teff/835.*sqrt(cfact)
             ccf = (1.+1./fac)*(1.+sqrt(2./(2.+fac)))
             react_rate = react_rate*cfact*ccf
          else if (react(i)%type.eq.'CHARGR') then
             ! Product is a charged grain
             react_rate = react_rate*cfact
          endif
          creation_rate = react_rate*densityr1*densityr2
          ! if (ip1.eq.211.or.ir1.eq.211) then
          ! if (ip1.eq.65.or.ir1.eq.65) then
          !    write(*,'(i6,10es14.1e4)') i,ts,tr,teff,exp_factor1,exp_factor2,exp_factor,react_rate,densityr1,densityr2,creation_rate
          ! endif
          ! if (react(i)%type.eq.'NEUTGR'.or.react(i)%type.eq.'CHARGR') then
          !    write(*,'(i6,10es14.1e4)') i,teff,cfact,fac,ccf,react_rate,densityr1,densityr2,creation_rate
          ! endif
       end select
       !
       ! creation rate (cm-3.s-1), destruction rate (s-1) of each reactant 
       !  creation computed in select case above
       !  multiplication by density done already
       !      creation_rate = react_rate * densityR1 * densityR2
       !      destr_R1      = react_rate * densityR2
       !      destr_R2      = react_rate * densityR1
       destr_R1      = creation_rate
       destr_R2      = creation_rate
       !
       !
       ! source terms for destruction (multiplied by mass later) 
       destr1_Vcm   = destr_R1 * Vcm
       destr2_Vcm   = destr_R2 * Vcm
       destr1_Vcm2  = destr_R1 * Vcm2
       destr2_Vcm2  = destr_R2 * Vcm2
       !
       ! source terms for creation 
       creation_Vcm  = creation_rate * Vcm
       creation_Vcm2 = creation_rate * Vcm2
       creation_DE   = creation_rate * react(i)%DE
       !
       ! for the reactants, increase in the rate of destruction (s-1) of : 
       !      down_N(IR)   -> number                                       
       !      down_MV(IR)  -> momentum                                     
       !      down_MV2(IR) -> energy                                       
       down_n(ir1)   = down_n(ir1)   + destr_r1
       down_n(ir2)   = down_n(ir2)   + destr_r2
       down_mv(ir1)  = down_mv(ir1)  + destr1_vcm
       down_mv(ir2)  = down_mv(ir2)  + destr2_vcm
       down_mv2(ir1) = down_mv2(ir1) + destr1_vcm2
       down_mv2(ir2) = down_mv2(ir2) + destr2_vcm2
       !
       ! for the products, increase in the rate of creation (s-1) of :    
       !      up_N(IP)   -> number                                        
       !      up_MV(IP)  -> momentum                                      
       !      up_MV2(IP) -> kinetic energy                                
       !      up_DE(IP)  -> energy defect of the reaction                 
       !                    DE is distributed over the reaction products  
       !                    with a weighing factor : [1], eq. (31)        
       !                        =  DE/(n-1)*[1-mass(i)/MASS] if n > 1     
       !                        =  DE if n=1                              
       !                    where n=number of produts, i=index of product 
       !                    and MASS=sum of mass(i), i=1..n               
       !                    By this way, energy is conserved wathever     
       !                    the number of products in the reaction.       
       ! first product, always present !
       up_n(ip1)   = up_n(ip1)   + creation_rate
       up_mv(ip1)  = up_mv(ip1)  + creation_vcm
       up_mv2(ip1) = up_mv2(ip1) + creation_vcm2
       if (nprod_m1 == 0) then
          ! allow for reactions with only one product : H2_FO, ADSOR
          !  JLB - May 2001
          !  1/3 of H2 enthalpy formation goes To kinetic energy only
          !  Maybe should be corrected else-where ???
          !         up_DE(IP1) = up_DE(IP1) + creation_DE
          !  + (24 I 2002) account for lost of kinetic energy by H attaching to grain
          !    (proposed by David)
          if(react(i)%type=="H2_FO")then
             if(ikinh2==1)then
                up_de(ip1)=up_de(ip1)+0.5_dp*(creation_de-creation_rate*h2_int_e*kb/everg)&
                     -creation_rate*1.5_dp*kb*speci(ir2)%temperature/everg
             else
                up_de(ip1)=up_de(ip1)+min(creation_de/3.0_dp,creation_de-creation_rate*h2_int_e*kb/everg)&
                     -creation_rate*1.5_dp*kb*speci(ir2)%temperature/everg
             endif
          else
             up_de(ip1) = up_de(ip1) + creation_de
          endif
       else
          ! two or more products
          up_de(ip1) = up_de(ip1)+creation_de/nprod_m1*(1.0_dp-speci(ip1)%mass/mass_prod)
       endif
       !
       ! second product
       if (ip2 > 0) then
          up_n(ip2)   = up_n(ip2)   + creation_rate
          up_mv(ip2)  = up_mv(ip2)  + creation_vcm
          up_mv2(ip2) = up_mv2(ip2) + creation_vcm2
          ! in this case, the reaction has at least two products !
          up_de(ip2)  = up_de(ip2)  + creation_de/nprod_m1 * &
               (1.0_dp-speci(ip2)%mass/mass_prod)
       end if
       ! third product
       if (ip3 > 0) then
          up_n(ip3)   = up_n(ip3)   + creation_rate
          up_mv(ip3)  = up_mv(ip3)  + creation_vcm
          up_mv2(ip3) = up_mv2(ip3) + creation_vcm2
          ! in this case, the reaction has at least three products !
          up_de(ip3)  = up_de(ip3)  + creation_de / nprod_m1 * &
               (1.0_dp-speci(ip3)%mass/mass_prod)
       end if
       ! fourth product
       if (ip4 > 0) then
          up_n(ip4)   = up_n(ip4)   + creation_rate
          up_mv(ip4)  = up_mv(ip4)  + creation_vcm
          up_mv2(ip4) = up_mv2(ip4) + creation_vcm2
          ! in this case, the reaction has four products !
          up_de(ip4)  = up_de(ip4)  + creation_de / nprod_m1 * &
               (1.0_dp-speci(ip4)%mass/mass_prod)
       end if
    end do
    !
    !
    ! Calculate for every chemical species the change per unit volume
    ! and time in :
    !       YN(i) -> number (cm-3.s-1)
    !       YS(i) -> mass (g.cm-3.s-1)
    !       YA(i) -> momentum (g.cm-2.s-2)
    !       YB(i) -> energy (erg.cm-3.s-1)
    yn(1:nspec) = up_n(1:nspec) - down_n(1:nspec)
    ys(1:nspec) = yn(1:nspec) * speci(1:nspec)%mass
    ya(1:nspec) = (up_mv(1:nspec) - down_mv(1:nspec)) * speci(1:nspec)%mass
    yb(1:nspec) = (up_mv2(1:nspec) - down_mv2(1:nspec)) * speci(1:nspec)%mass
    !--------------------------------------------------------------------------
    ! Calculate the source terms (the changes per unit volume and time)
    ! of each fluid by summing over the source terms for the chemical species :
    !      Cup    -> increase in number (cm-3.s-1)
    !      Cdown  -> decrease in number (cm-3.s-1)
    !      YN     -> net change in number (cm-3.s-1) [1], eq. (16,17,20)
    !      S      -> mass (g.cm-3.s-1)               [1], eq. (18,19)
    !      A      -> momentum (g.cm-2.s-2)           [1], eq. (21)
    !      B      -> energy (erg.cm-3.s-1)           [1], eq. (25,26,31)
    !--------------------------------------------------------------------------
    !
    ! neutral species
    cupn   = sum(up_n(b_neu:e_neu))
    cdownn = sum(down_n(b_neu:e_neu))
    ynn    = sum(yn(b_neu:e_neu))
    sn     = sum(ys(b_neu:e_neu))
    an     = sum(ya(b_neu:e_neu))
    bn     = 0.5_dp * sum(yb(b_neu:e_neu)) + sum(up_de(b_neu:e_neu)) * everg
    !
    ! neutrals on mantles : add these terms to the neutrals
    cupn   = cupn   + sum(up_n(b_gra:e_gra))
    cdownn = cdownn + sum(down_n(b_gra:e_gra))
    ynn    = ynn    + sum(yn(b_gra:e_gra))
    sn     = sn     + sum(ys(b_gra:e_gra))
    an     = an     + sum(ya(b_gra:e_gra))
    bn     = bn     + 0.5_dp * sum(yb(b_gra:e_gra)) + sum(up_de(b_gra:e_gra)) * everg
    !
    ! neutrals eroded from grain cores : add these terms to the neutrals
    cupn   = cupn   + sum(up_n(b_cor:e_cor))
    cdownn = cdownn + sum(down_n(b_cor:e_cor))
    ynn    = ynn    + sum(yn(b_cor:e_cor))
    sn     = sn     + sum(ys(b_cor:e_cor))
    an     = an     + sum(ya(b_cor:e_cor))
    bn     = bn     + 0.5_dp * sum(yb(b_cor:e_cor)) + sum(up_de(b_cor:e_cor)) * everg
    !
    ! positive ions
    cupi   = sum(up_n(b_ion:e_ion))
    cdowni = sum(down_n(b_ion:e_ion))
    yni    = sum(yn(b_ion:e_ion))
    si     = sum(ys(b_ion:e_ion))
    ai     = sum(ya(b_ion:e_ion))
    bi     = 0.5_dp * sum(yb(b_ion:e_ion)) + sum(up_de(b_ion:e_ion)) * everg
    !
    ! negative ions
    cupneg   = sum(up_n(b_neg:e_neg))
    cdownneg = sum(down_n(b_neg:e_neg))
    ynneg    = sum(yn(b_neg:e_neg))
    sneg     = sum(ys(b_neg:e_neg))
    aneg     = sum(ya(b_neg:e_neg))
    bneg     = 0.5_dp * sum(yb(b_neg:e_neg)) + sum(up_de(b_neg:e_neg)) * everg
    !
    ! heating of the electrons by photo-ionization
    bneg     = bneg + up_de(ind_e) * everg
    !
    ! change in mass density of grains by erosion
    score    = sum(ys(b_cor:e_cor))

!!$    write(*,'(a,20(es12.1))') "PHB, Yn arrays: ",YN(1:Nspec)
!!$    write(*,'(a,20(es12.1))') "PHB, Ys arrays: ",Ys(1:Nspec)
!!$    write(*,'(a,20(es12.1))') "PHB, Ya arrays: ",Ya(1:Nspec)
!!$    write(*,'(a,20(es12.1))') "PHB, Yb arrays: ",Yb(1:Nspec)
!!$    write(*,'(a,20(es12.1))') "PHB, n  arrays: ",CupN,CdownN,YNn,sn,an,bn
!!$    write(*,'(a,20(es12.1))') "PHB, i  arrays: ",Cupi,Cdowni,YNi,si,ai,bi
!!$    write(*,'(a,20(es12.1))') "PHB, negarrays: ",Cupneg,Cdownneg,YNneg,sneg,aneg,bneg
!!$    write(*,*) dens_grain
  end subroutine evolution_chemistry



  subroutine evolution_physics
    use module_constants, only : pi,kb,me,everg,zero,alpha_h,alpha_h2,alpha_he,qe,gr_cons
    use module_gamma, only : gamma1
    use module_molecular_cooling
    use module_phys_var
    use module_chemical_species
    use module_h2
    use module_grains, only : tgrain, dens_grain, rgrain2
    use module_debug
    use module_line_excit
    use module_parameters_flags
    !___________________________________________________________________
    !
    ! First calculates the cooling of the different fluids, then
    ! computes the source terms needed in DIFFUN: momentum and
    ! energy. These terms are updated from there values calculated in
    ! CHEMISTRY. Note that source terms for number and mass density
    ! are not modified.
    !
    ! Source takes into account the cooling due to H2 lines, fine structure
    ! lines and other molecular cooling.
    ! It includes also other physical processes :
    !    - heat transfer between neutral and charged fluids;
    !    - photo-electric effect on grain;
    !    - thermalisation with grains;
    !    - momentum and energy transfer between neutral and charged fluids
    !      through elastic scattering;
    !    - energy transfer between charged fluids through
    !      Coulomb scattering.
    ! results :
    !     An,  Ai,  Aneg  -> source terms for momentum (g.cm-2.s-2)
    !     Bn,  Bi,  Bneg  -> source terms for energy (erg.cm-3.s-1)
    ! references :
    !     [1]  Flower et al., 1985, MNRAS 216, 775
    !     [2]  Flower et al., 1986, MNRAS 218, 729
    !     [3]  Pineau des Forets et al., 1986, MNRAS 220, 801
    !     [4]  Monteiro & Flower, 1987, MNRAS 228, 101
    !     [5]  Black, 1987, in Interstellar processes, eds. Hollenbach & !
    !          Thronson (Reidel, Dordrecht), P.731
    !     [6]  Viala et al., 1988, A&A 190, 215
    !     [7]  Spitzer, 1962, Physics of fully ionised gases (Wiley, New York)
    !     [8]  Spitzer & Scott, 1969, AP.J. 158, 161
    !     [9]  Flower & Watt, 1984, MNRAS 209, 25
    !     [10] Roueff 1990, A&A 234, 567
    !___________________________________________________________________
    !
    implicit none
    real(kind=dp) :: cooling_n, cooling_i, cooling_neg ! cooling rate (erg/cm3/s)
    real(kind=dp) :: teff, ve, alpha_n
    real(kind=dp) :: mun_plus_mui, mun_plus_muneg, rhon_times_rhoi, mun_neg
    real(kind=dp) :: sigma_in, sigma_en, sigma_negn, sigma_inelastic_en
    real(kind=dp) :: a_i_n, a_e_n, a_neg_n, a_grain_n
    real(kind=dp) :: b_i_n, b_neg_n, b_e_n, b_inelastic_e_n, b_grain_n
    real(kind=dp) :: b_heat_exchange_n, b_photo_grain_n, b_therm_grain_n, b_ioniz_rc_n
    real(kind=dp) :: lambda, b_i_e, b_i_neg, b_heat_exchange_i
    real(kind=dp) :: b_heat_exchange_neg
    real(kind=dp) :: sel_net_h2
    !
    !
    !----------------------------------------------------------
    ! in FINE_STRUCTURE_COOLING, cooling rates (erg/cm3/s) for
    ! fine structure lines of C+, Si+, C, Si, O, S+, and N+
    ! are calculated.
    ! results are :
    ! cooling_Cplus, cooling_Siplus, cooling_C, cooling_Si
    ! cooling_O, cooling_Splus, cooling_Nplus, cooling_Cplus_e,
    ! cooling_Siplus_e, cooling_Splus_e, cooling_Nplus_e
    !----------------------------------------------------------
    call LINE_THERMAL_BALANCE

    !---------------------------------------------------------
    ! in COMPUTE_MOLECULAR_COOLING, cooling rates (erg/cm3/s)
    ! 13CO, OH, NH3, CO and H2O are calculated.
    ! useful results for SOURCE is :
    !    molec_cool
    !---------------------------------------------------------
    !
    ! 29 mai 2001 - JLB : compute molec_cool directly in DIFFUN
    !
    !   CALL COMPUTE_MOLECULAR_COOLING

    !-----------------------------------------------------
    ! in COMPUTE_H2, cooling rate (erg/cm3/s) for
    ! ro-vibrational lines of H2 and source terms for
    ! H2 levels are calculated.
    ! results are :
    ! Cooling_H2, YN_rovib_H2, H2_Energy
    !-----------------------------------------------------
    CALL COMPUTE_H2

    !------------------------------------------------------------------
    ! Add different contributions to the cooling of the neutral fluid
    ! cooling_n (erg/cm3/s) :
    ! (1) fine structure lines (from FINE_STRUCTURE_COOLING)
    !     remark : we assume that the neutrals loss most of the energy
    !              in the collisions with C+ or Si+. This is why the
    !              terms cooling_Cplus and cooling_Siplus are here and
    !              not in the cooling of the positively charged fluid.
    ! (1) - bis : Superseded by LINE_THERMAL_BALANCE
    !     For each collision, energy loss is spread accurately to
    !     each fluid.
    ! (2) molecular (not H2) cooling (from COMPUTE_MOLECULAR_COOLING)
    ! (3) H2 rovibrational cooling (from COMPUTE_H2)
    ! (4) dissociation of H2 by collisions with H, He and H2
    !     we assume a rate 8 times smaller for H2 than for H
    !     and 10 times smaller for He
    !     remark :
    !        It is not taken into account in CHEMISTRY as DE=0._DP
    !        for endothermic reactions (see ENERGY_DEFECT_REACTION).
    !------------------------------------------------------------------
    cooling_n = Zero
    ! (1)
    cooling_n = cooling_n - Cool_n
    ! (2)
    !
    ! 29 mai 2001 - JLB : added directly in DIFFUN
    !
    !   cooling_n = cooling_n + &
    !        molec_cool
    ! (3)
    cooling_n = cooling_n + cooling_H2
    ! (4) to be supressed if endothermic reaction don't have DE = Zero
    ! Use real collisional dissociation by neutrals (56641.590 is energy of highest H2 level (in K))
    sel_net_h2 = sum(dble(h2_lev(1:nh2_lev)%density) * sel_ne_h2 &
         * (56641.590-h2_lev(1:nh2_lev)%energy))
    cooling_n = cooling_n-sel_net_h2*kb

    !   cooling_n = cooling_n + &
    !        Dens_H * Dens_H2 * 7.2D-22*EXP(-52000._DP/Tn) + &
    !        Dens_H2* Dens_H2 * 9.0D-23*EXP(-52000._DP/Tn)

    !------------------------------------------------------------------------
    ! Add the different contributions to the rate of radiative cooling of the
    ! negatively charged fluid cooling_neg (erg/cm3/s) through excitation by
    ! electron collision of :
    ! (1) O and C : 3P -> 1D; N 4S -> 2D (MENDOZA, 1983, PNE SYMPOSiUM)
    ! (2) C+, Si+, S+, and N+ (from FINE_STRUCTURE_COOLING)
    ! (1) and (2) - bis : This done in LINE_THERMAL_BALANCE now
    ! (3) rovibrational excitation of H2, rotational excitation of CO
    !-------------------------------------------------------------------------
    ! (1) and (2)
    cooling_neg = Zero
    cooling_neg = cooling_neg - Cool_e
    ! (3)
    Ve = SQRT(8.0_DP*kB*Te/(pi*me)) ! thermal velocity of e- (cm/s)
    cooling_neg = cooling_neg + &
         dens_e*dens_h2* everg * &
         ve * (4.0d-18*exp(-510._dp/te)*0.044_dp + 4.0d-17*exp(-6330._dp/te)*0.544_dp)
    cooling_neg = cooling_neg + &
         dens_e * dens_co * everg * ve * 1.0d-15*exp(-5.0_dp/te)*0.00043_dp

    !-------------------------------------------------------------------------
    ! Add the different contributions to the rate of radiative cooling of the
    ! positively charged fluid cooling_i (erg/cm3/s) :
    ! (1) ro-vibrational excitation of H2 by collisions with ions
    ! (2) dissociation of H2 by collisions with ions
    !     remark :
    !        It is not taken into account in CHEMISTRY as DE=Zero
    !        for endothermic reactions (see ENERGY_DEFECT_REACTION).
    ! (3) ionization of H2 by collisions with ions
    ! The energy loss by the positively charged fluid is proportional to
    !     mass_H2/(mass_H2 + muI); the corresponding loss by the neutral
    ! fluid is proportional to
    !     muI/(mass_H2 + muI)
    ! (4) ionization of H by collisions with ions
    ! The energy loss by the positively charged fluid is proportional to
    !     mass_H/(mass_H + muI); the corresponding loss by the neutral
    ! fluid is proportional to
    !     muI/(mass_H + muI)
    !-------------------------------------------------------------------------
    cooling_i = Zero
    ! (1)    
    Teff = & ! effective temperature, used in (1), (2) and (3)
         (mass_pH2*muI*ABS_DeltaV*ABS_DeltaV/3.0_DP/kB + (mass_pH2*Ti + muI*Tn)) / &
         (mass_pH2 + muI)
    cooling_i = cooling_i + &
         ABS_DeltaV * EVerg * dens_h2 * DensityI * &
         (4.0D-16*EXP(-510._DP/Teff)*0.044_DP + 4.0D-17*EXP(-6330._DP/Teff)*0.544_DP)
    ! (2) to be supressed if endothermic reaction don't have DE = Zero
    Sel_net_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density) * Sel_io_H2 &
         * (56641.590-H2_lev(1:NH2_lev)%energy))
    cooling_i = cooling_i - Sel_net_H2 * kB
    ! (3) to be supressed if endothermic reaction don't have DE = Zero
    cooling_i = cooling_i + &
         DensityI * dens_h2 * 1.10D-13*(Teff/300._DP)**0.5_DP* &
         EXP(-179160.0_DP/Teff)*15.43_DP * EVerg
    cooling_i = cooling_i * mass_ph2/(mass_pH2 + muI)
    cooling_n = cooling_n + cooling_i * muI/mass_pH2
    ! (4) to be supressed if endothermic reaction don't have DE = Zero
    Teff = & ! effective temperature used in (4)
         (mass_H*muI*ABS_DeltaV*ABS_DeltaV/3.0_DP/kB + (mass_H*Ti + muI*Tn)) / &
         (mass_H + muI)
    cooling_i = cooling_i + &
         DensityI * Dens_H * 1.30D-13*(Teff/300._DP)**0.5_DP* &
         EXP(-157890.0_DP/Teff)*13.60_DP * EVerg * mass_H/(mass_H + muI)
    cooling_n = cooling_n + &
         DensityI * Dens_H * 1.30D-13*(Teff/300._DP)**0.5_DP* &
         EXP(-157890.0_DP/Teff)*13.60_DP * EVerg * muI/(mass_H + muI)

    ! (0) add contribution from LINE_THERMAL_BALANCE
    ! Added after (1)-(4) so as not to add it also to cooling_n
    cooling_i = cooling_i - Cool_i

    !--- some useful variables ---
    muN_plus_muI    = muN  + muI           ! sum of neutral and positive ion mass (g)
    muN_plus_muNEG  = muN  + muNEG         ! sum of neutral and negative ion mass (g)
    RhoN_times_RhoI = RhoN * RhoI          ! product of neutral and positive ion mass densisities (g2/cm6)
    muN_NEG = muN * muNEG / muN_plus_muNEG ! reduced mass of a neutral/negative ion pair (g)

    !--- polarisability of the neutrals                    ---
    !--- = weighted average of polarisability of H and H2  ---
    alpha_n=(dens_h*alpha_h+dens_h2*alpha_h2) / (dens_h+dens_h2)

    !------------------------------------------------------------------
    ! calculate the cross-sections (cm-2)
    ! (1) Sigma_IN           -> elastic ions >0  - neutral scattering
    !                           ([1], eq.(23))
    ! (2) Sigma_eN           -> elastic electron - neutral scattering
    !                           ([1], eq.(34))
    ! (3) Sigma_NegN         -> elastic ions <0  - neutral scattering
    !                           ([1], eq.(23))
    ! (4) Sigma_inelastic_eN -> inelastic electron-neutral scattering
    !------------------------------------------------------------------
    ! (1)
    Sigma_IN = 2.41_DP * pi * qe * SQRT(muN_plus_muI*alpha_n/muN/muI)
    Sigma_IN = MAX(Sigma_IN, 1.0D-15*ABS_DeltaV)
    ! (2)
    Sigma_eN = 1.0D-15 * SQRT(8.0_DP*kB*Te/(pi*me))
    ! (3)
    Sigma_NegN = 2.41_DP * pi * qe * SQRT(alpha_n/muN_NEG)
    Sigma_NegN = MAX(Sigma_NegN, 1.0D-14*ABS_DeltaV)
    ! (4)
    Sigma_inelastic_eN = 1.0D-16 * (4.29D-1 + 6.14D-6*Te) * SQRT(8.0_DP*kB*Te/(pi*me))

    !--------------------------------------------------------------------
    ! Rates of production of momentum (g.cm-2.s-2) in the neutral fluid
    !   A_i_n     -> elastic ions >0  - neutral scattering
    !   A_e_n     -> elastic electron - neutral scattering
    !   A_neg_n   -> elastic ions <0  - neutral scattering
    !   A_grain_n -> elastic charged grains - neutral scattering
    !--------------------------------------------------------------------
    A_i_n     = Sigma_IN * DeltaV * RhoN_times_RhoI / muN_plus_muI
    !  i_n ambipolar diffusion timescale [equ. (5) of Walmsley et al. 2004, A&A, 418, 1035]
    tau_i_n   = (2/pi/gr_cons) * Sigma_IN * RhoI / muN_plus_muI / RhoN
    A_e_n     = DensityN * Dens_e * me * Sigma_eN * DeltaV
    A_neg_n   = DensityN * DensityNEG * muN_NEG * Sigma_NegN * DeltaV
    A_grain_n = RhoN * Dens_grain * pi * Rgrain2 * ABS_DeltaV * DeltaV &
         * (Dens_ggminus + Dens_ggplus)/ (Dens_ggminus + Dens_ggplus + Dens_gg0)

    !------------------------------------------------------------------
    ! Source terms of number and mass density are calculated in
    ! CHEMISTRY and not in SOURCE.
    ! These terms are (for neutrals, positive fluid, negative fluid) :
    ! YNn, YNi, YNneg -> change in number density (cm-3.s-1)
    ! Sn,  Si,  Sneg  -> change in mass density (g.cm-3.s-1)
    !------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Source terms of momentum (g.cm-2.s-2) of each fluid (update of
    ! the calculation in CHEMISTRY).
    !   An   -> neutral fluid
    !   Ai   -> positive fluid
    !   Aneg -> negative fluid
    !----------------------------------------------------------------
    An   = An   + A_i_n + A_e_n + A_neg_n + A_grain_n
    Ai   = Ai   - A_i_n
    Aneg = Aneg         - A_e_n - A_neg_n - A_grain_n

    !--------------------------------------------------------------------------
    ! Source term of energy Bn (erg.cm-3.s-1) in the neutral fluid
    ! obtained by adding up :
    ! +Bn                 -> exo/endothermicity of the reactions
    !                        (from CHEMISTRY)
    ! -cooling_n          -> energy loss through cooling
    ! +B_i_n              -> elastic diffusion of the neutrals by the
    !                        positive ions ([1], eq.(32))
    ! +B_neg_n            -> elastic diffusion of the neutrals by the
    !                        negative ions ([1], eq.(32))
    ! +B_e_n              -> elastic diffusion of the neutrals by the
    !                        electrons ([1], eq.(33))
    ! +B_inelastic_e_n    -> inelastic scattering of the electrons on H2;
    !                        mean excitation energy 12 eV = 140 000 K for the
    !                        singlet excited states, 10 eV = 116 300 K for the
    !                        triplet, of which 10 - 4.48 = 5.52 eV are recovered
    !                        by the neutrals through the dissociation products;
    !                        collisional excitation of the N=2 state of H.
    ! +B_grain_n          -> elastic diffusion of the neutrals by the
    !                        charged grains
    ! +B_heat_exchange_n  -> heat exchange
    !                        [1],eq.(27) with the 5/2 factors corrected to 3/2
    !                        (GAMMA -> GAMMA1)
    ! +B_photo_grain_n    -> heating through the photo-electric effect on grains
    !                        [5]
    ! +B_therm_grain_n    -> energy loss/gain through thermalisation with grains
    !                        TIELENS  & HOLLENBACH (1987, AP.J. 291, 722)
    ! +B_ioniz_RC_n       -> heating through cosmic ray ionization [8]
    !                        (not included in CHEMISTRY as DE = Zero for
    !                        reactions of the type 'CR_IO')
    !--------------------------------------------------------------------------

    B_i_n = (3.0_DP*kB*(Ti-Tn) + DeltaV*(muI*Vi+muN*Vn)) * &
         RhoN_times_RhoI/muN_plus_muI * Sigma_IN/muN_plus_muI
    B_neg_n = (3.0_DP*kB*(Te-Tn) + DeltaV*(muNEG*Vi+muN*Vn)) * &
         DensityN*DensityNEG * Sigma_NegN*muN_NEG/muN_plus_muNEG
    B_e_n = (4.0_DP*kB*(Te-Tn) + muN*Vn*DeltaV) * &
         DensityN*Dens_e * me/muN*Sigma_eN
    B_inelastic_e_n = - SUM(DBLE(H2_lev(1:NH2_lev)%density) * Sel_el_H2) * 5.52_DP*EVerg
    B_grain_n = A_grain_n*Vi
    B_heat_exchange_n = GAMMA1*kB*(CupN*(Ti+Te)/2._DP - CdownN*Tn)
    B_photo_grain_n = 4.0D-26*nH*RAD*EXP(-2.5_DP*AV)
    B_therm_grain_n = 3.5D-34*SQRT(Tn)*(Tgrain-Tn)*nH**2.0_DP
    B_ioniz_RC_n = Zeta * EVerg * DensityN * &
         MAX(5.7_DP,32.0_DP-7.0_DP*LOG10(DensityN/Dens_e))

    !--- sum of the different contributions to Bn ---
    Bn = Bn                  &
         - cooling_n         &
         + B_i_n             &
         + B_neg_n           &
         + B_e_n             &
         + B_inelastic_e_n   &
         + B_grain_n         &
         + B_photo_grain_n   &
         + B_therm_grain_n   &
         + B_ioniz_RC_n      !&
    !+ B_heat_exchange_n ! include only Bn and not B_heat_exchange_n (mars 2000)

    !--------------------------------------------------------------------------
    ! Source term of energy Bi (erg.cm-3.s-1) for the positive ions fluid
    ! obtained by adding up :
    !  +Bi                -> exo/endothermicity of the reactions
    !                        (from CHEMISTRY)
    !  -cooling_i         -> energy loss through cooling
    !  -B_i_e             -> diffusion of electrons by positive ions
    !                        [7], see also [1], eq.(35)
    !  -B_i_neg           -> diffusion of negative ions by positive ions
    !                        [7]
    !  -B_i_n             -> elastic diffusion of the neutrals by the
    !                        positive ions. See the calculation for Bn above
    !  +B_heat_exchange_i -> heat exchange
    !                        [1],eq.(28) with the 5/2 factors corrected to 3/2
    !                        (GAMMA -> GAMMA1)
    !--------------------------------------------------------------------------
    lambda = 1.5_DP/(qe**3.0_DP) * SQRT((kB*Te)**3.0_DP/(Dens_e*pi))
    B_i_e = 4.0_DP*qe**4._DP/muI/Te * SQRT(2.0_DP*me*pi/(kB*Te)) * &
         DensityI*Dens_e*LOG(lambda)*(Ti-Te)

    lambda = 1.5_DP/(qe**3.0_DP) * SQRT((kB*Te)**3.0_DP/(DensityI*pi))
    B_i_neg = qe**4.0_DP/(muI*muNEG) * (Ti/muI+Te/muNEG)**(-1.5_DP)* &
         (32.0_DP*pi/kB)**0.5_DP*DensityI*DensityNEG*LOG(lambda)*(Ti-Te)

    B_heat_exchange_i=GAMMA1*kB*(CupI*Tn-CdownI*Ti)

    !--- sum of the different contributions to Bi ---
    Bi = Bi                   &
         - cooling_i          &
         - B_i_e              &
         - B_i_neg            &
         - B_i_n              &
         + B_heat_exchange_i

    !--------------------------------------------------------------------------
    ! CALCULATE THE SOURCE TERM OF ENERGY Bneg (ERG.CM-3.S-1) FOR THE
    ! NEGATIVE FLUID BY ADDING UP
    !   B_i_e              -> diffusion of electrons by positive ions
    !                         [7], see also [1], eq.(35).
    !                         see calculation for Bi above.
    !   B_i_neg            -> diffusion of negative ions by positive ions
    !                         [7], see calculation for Bi above.
    !  -B_e_n              -> elastic diffusion of the neutrals by the
    !                         electrons ([1], eq.(33)).
    !                         see calculation for Bn above.
    !  -B_inelastic_e_n    -> inelastic scattering of the electrons on H2 and H.
    !                         For H, excitation of the N=2 state :
    !                         Aggarwal et al. 1991, J. PHYS. B, 24, 1385
    !                                ionization:
    !                         Hummer, 1963, MNRAS, 125, 461.
    !                         For H2, excitation of the repulsive B triplet state
    !                         and of the B and C singlet states
    !                                ionization:
    !                         Rapp & Englander-Golden, 1965, JCP, 43, 1464.
    !  -B_grain_n          -> elastic diffusion of the neutrals by the
    !                         charged grains, see calculation for Bn above.
    !  -B_neg_n            -> elastic diffusion of the neutrals by the
    !                         negative ions ([1], eq.(32)).
    !                         see calculation for Bn above.
    !  cooling_neg         -> energy loss through cooling
    !  B_heat_exchange_neg -> heat exchange
    !                        [1],eq.(28) with the 5/2 factors corrected to 3/2
    !                        (GAMMA -> GAMMA1)
    !--------------------------------------------------------------------------
    b_inelastic_e_n = sigma_inelastic_en*exp(-1.4d5/te)
    sel_net_h2 = sum(dble(h2_lev(1:nh2_lev)%density) * sel_el_h2 &
         * (116300.0_dp-h2_lev(1:nh2_lev)%energy))
    b_inelastic_e_n = dens_h2*dens_e * &
         b_inelastic_e_n * 12.0_dp * everg - &
         sel_net_h2 * kb
    b_inelastic_e_n = b_inelastic_e_n + &
         dens_h*dens_e*everg*3.8d-06*exp(-1.1842d5/te)*(te/300._dp)**(-0.5_dp)
    b_inelastic_e_n = b_inelastic_e_n + &
         dens_h*dens_e*everg*13.6_dp*9.2d-10*exp(-1.5789d5/te)*(te/300.0d0)**0.5_dp + &
         dens_h2*dens_e*everg*15.43_dp*1.40d-09*exp(-1.7916d5/te)*(te/300._dp)**0.5_dp
    b_heat_exchange_neg=0.5_dp*me*vi*vi*(yni-ynneg)-gamma1*kb*te*cdowni
    !
    ! Sum of the different contributions to Bneg
    Bneg = Bneg                 &
         - cooling_neg          &
         + B_i_e                &
         + B_i_neg              &
         - B_e_n                &
         - B_inelastic_e_n      &
         - B_grain_n            &
         - B_neg_n              &
         + B_heat_exchange_neg
  end subroutine evolution_physics



  subroutine diffun(n,z,y,dy)
    use module_phys_var
    use module_chemical_species
    use module_gamma
    use module_grains
    use module_molecular_cooling
    use module_h2
    use module_constants, only : kb,pi,zero,parsec,yearsec
    use module_debug
    use module_tools
    use module_parameters_flags
    !___________________________________________________________________________
    !
    ! Calculates the partial derivatives of the Y vector.
    !
    ! DY(i)/DZ= F(Y,Z), i=1,N
    !
    ! (see Flower et al., 1985, MNRAS 216, 775, equations (3)-(14))
    !
    ! Source terms (cm-3 s-1) from 'chemistry' and 'source'.    
    !___________________________________________________________________________
    !
    implicit none
    integer(kind=long),         intent(in)    :: n  ! dimension of y and dy
    real(kind=dp),              intent(in)    :: z  ! useless here : the dependence in z is implicit
    real(kind=dp),dimension(n), intent(in)    :: y  ! vector containing the variables
    real(kind=dp),dimension(n), intent(out)   :: dy ! vector containing the derivatives
    !
    ! Local variables
    integer            :: ii, i
    integer(kind=long) :: b_specy ! useful for densities of chemical species
    real(kind=dp)      :: bfield2, kt_mun, vi2, vn2, v_densityi ! useful for the calculation
    real(kind=dp)      :: aux1_dvn, aux2_dvn, dvn1, dvn2, gvn1, gvn2
    real(kind=dp)      :: xll2, bne, ane, sne, ynne, densityne, rhone
    real(kind=dp)      :: m_grain, v_grain, v_crit, v_turb, v_rel,dens_ggtot
    real(kind=dp)      :: mindens
    real(kind=dp)      :: fact,dilution
    real(kind=dp)      :: mcore,mmantle
    logical :: jshock,jshockvisc
    !
    !
    call2diffun = call2diffun+1
    !
    ! Sanity checks
    ! * NaN
    !
    call check_nan_dp(y,n,"diffun","yarray")
    yarray(1:n)=y(1:n)
    ! Set to the minimal density where needed
!!$    yarray = max(yarray,0._dp)
    ! Negative values
    if (any(yarray(bv_speci:ev_speci).lt.0.d0)) then
       write(*,*) "W-diffun: found negative abundances at"
       do i=bv_speci,ev_speci
          if (yarray(i).lt.0.d0) then
             write(*,'(i6,2es20.7e5)') i,yarray(i),dy(i)
             write(*,*) speci(i-bv_speci+1)%index
          endif
       enddo
       stop
    endif
    !
    ! Initialization
    jshock = do_shock.and.shock_type.eq."J"
    jshockvisc = jshock.and.viscosity
    !
    ! Initialize number densities for some species
    ! * H, oH2, pH2, H+, G0, G+, G-
    ! * number densities (cm-3) : neutrals, ions >0, ions <0, grains (core, adsorbed)
    ! * electrons
    call get_density(yarray,ind_h,dens_h)
    call get_density(yarray,ind_oh2,dens_oh2)
    call get_density(yarray,ind_ph2,dens_ph2)
    call get_density(yarray,ind_hplus,dens_hplus)
    call get_density(yarray,ind_gg0,dens_gg0)
    call get_density(yarray,ind_ggplus,dens_ggplus)
    call get_density(yarray,ind_ggminus,dens_ggminus)
    b_specy       = bv_speci-1
    dens_grain    = yarray(iv_dens_grain)
    dens_ggtot    = dens_gg0+dens_ggplus+dens_ggminus
    dens_ongrains = sum(yarray(bv_gra:ev_gra))
    dens_e        = speci(ind_e)%density
    dens_h2       = dens_oh2+dens_ph2
!!$    if (debug) write(*,'("D-evolution: grains properties",10es12.3)') &
!!$         dens_gg0,dens_ggplus,dens_ggminus,dens_grain,dens_ggtot,dens_ongrains
    ! Ensure a minimum amount of adsorbed species for DESORPTION reactions
    dens_ongrains = max(dens_ongrains,1d-11)
    dens_cor = sum(yarray(bv_cor:ev_cor))
    ! Check ion conservation
    if (force_i_c == 1) yarray(iv_densityi) = sum(yarray(bv_ion:ev_ion))
    densityn   = yarray(iv_densityn)
    densityi   = yarray(iv_densityi)
    densityneg = sum(yarray(bv_neg:ev_neg))
    speci(ind_e)%density     = densityi-densityneg
    !
    ! Proton density (cm-3)
    ! nH=n(H)+2(n(p-H2)+n(o-H2))+n(H+)
    nh = dens_h+2.d0*dens_h2+dens_hplus
    !
    ! Check H2 level populations
    if (do_h2) then
       sum_h2 = sum(yarray(bv_h2_lev:ev_h2_lev))
       if (abs(sum_h2-(yarray(b_specy+ind_ph2)+yarray(b_specy+ind_oh2)))/sum_h2.gt.0.01) &
            write(*,'("W-diffun: H2 level populations not conserved: ",10es12.3)') &
            sum_h2,yarray(b_specy+ind_ph2)+yarray(b_specy+ind_oh2),speci(ind_ph2)%density,speci(ind_oh2)%density
       !    sum_h2 = yarray(b_specy+ind_ph2)+yarray(b_specy+ind_oh2)
    endif
    !
    ! Saveguard against round-off errors in the derivation of Vi, Ti, Te
    ! for one fluid model, we have :
    !    Vn=Vi (velocities), Tn=Ti=Te (temperatures)
    ! for two fluids model, we have :
    !    Ti=Te (temperatures)
    ! This can be checked in printouts periodically performed in MHD
    select case (nfluids)
    case (1)
       yarray(iv_vi) = yarray(iv_vn) ! vi=vn
       yarray(iv_ti) = yarray(iv_tn) ! ti=tn
       yarray(iv_te) = yarray(iv_tn) ! te=tn
    case (2)
       yarray(iv_te) = yarray(iv_ti) ! te=ti
    case default ! do nothing for three-fluid
    end select
    !
    ! Extracts physical variables
    ! * fractional radius of collapsing cloud
    ! * temperatures (K) : neutrals, ions >0, electrons and ions <0 ---
    ! * velocities (cm/s) : neutrals, ions, drift velocity ---
    ! * mass densities (g/cm3) : neutrals, ions >0, ions <0 ---
    ! * average masses (g) : neutrals, ions >0, ions <0 ---
    xsphere = yarray(iv_xsphere)
    if (jshock) grad_V = yarray(iv_gv)
    tn = yarray(iv_tn)
    ti = yarray(iv_ti)
    te = yarray(iv_te)
    vn = yarray(iv_vn)
    vi = yarray(iv_vi)
    deltav     = vi-vn
    abs_deltav = abs(deltav)
    rhon = yarray(iv_rhon)
    rhoi = yarray(iv_rhoi)
    rhoneg = yarray(iv_rhoneg)
    mun = rhon/densityn
    !>>>>>>> PHB
    mui   = epsion*rhoi/densityi
    muneg = epsneg*rhoneg/densityneg
    !<<<<<<< PHB
    !
    !--- variables useful for calculation ---
    kt_mun     = kb*tn/mun
    vn2        = vn*vn
    vi2        = vi*vi
    v_densityi = vi*densityi
    bfield2    = (vs_cm/vi*bfield)**2/4./pi
    !
    !--- sound speed and ion magnetosonic speed (cm.s-1) ---
    vsound  = sqrt(gamma*kt_mun)
    vmagnet = gamma*kb*(ti+te)*(densityi+densityneg)/(densityi*mui+densityneg*muneg)+bfield2/rhoi
    vmagnet = sqrt(vmagnet)
    !
    ! Chemical species (used in CHEMISTRY)
    ! * density, velocity, temperature
    ! neutrals : V = Vn, T = Tn
    ! species on grain mantles : V = Vi, T = Tgrain (not used)
    ! species in grain cores : V = Vi, T = Tgrain (not used)
    ! positive ions : V = Vi, T = Ti
    ! negative ions : V = Vi, T = Te
    ! electrons, with charge neutrality, and V = Vi, T = Te
    ! photons, cosmic rays and secondary photons
    ! density=cst=1.0_DP, V=cst=Zero, T=cst=Zero (READ_SPECIES) => do nothing
    ! V=cst=Zero, T=cst=Zero (READ_SPECIES) => do nothing
    speci(1:nspec)%density         = yarray(bv_speci:ev_speci)
    speci(b_neu:e_neu)%velocity    = vn
    speci(b_neu:e_neu)%temperature = tn
    speci(b_gra:e_gra)%velocity    = vi
    speci(b_gra:e_gra)%temperature = tgrain
    speci(b_cor:e_cor)%velocity    = vi
    speci(b_cor:e_cor)%temperature = tgrain
    speci(b_ion:e_ion)%velocity    = vi
    speci(b_ion:e_ion)%temperature = ti
    speci(b_neg:e_neg)%velocity    = vi
    speci(b_neg:e_neg)%temperature = te
    speci(ind_e)%velocity          = vi
    speci(ind_e)%temperature       = te
    !
    ! Update grain mass density (g/cm3), average radius (cm), square radius (cm2)
    ! Warning: core and mantles are accounted for separately
    ! Recall: dens_grain = #grains per unit volume of gas
    mcore   = dot_product(yarray(bv_cor:ev_cor),dble(speci(b_cor:e_cor)%mass))
    mmantle = dot_product(yarray(bv_gra:ev_gra),dble(speci(b_gra:e_gra)%mass))
    md_grain = mcore+mmantle
    vd_grain = mcore/rho_grain + mmantle/rho_mantle
    rgrain3 = 3.0_dp*vd_grain/(4.*pi*dens_grain)
    rgrain2 = r_gr_scal23*rgrain3**(2./3.)
    r_grain = r_gr_scal12*sqrt(rgrain2)
    m_grain = md_grain/dens_grain
    !
    !
    ! thermal speed of neutral grain
    v_grain = sqrt(16.0_dp*kb*tn/pi/m_grain)
    ! grains, compressed with ions, and V = 0, T = 0
    speci(ind_grain)%density = dens_grain
    !  mass of a grain (same for neutral and charged grains)
    !            speci(ind_C60)%mass = MD_grain/Dens_GRAIN
    !            speci(ind_C60plus)%mass = MD_grain/Dens_GRAIN
    !            speci(ind_C60minus)%mass = MD_grain/Dens_GRAIN
    !  mass of charged grains
    rho_grain_charge = md_grain * (dens_ggminus+dens_ggplus)/dens_ggtot
    !  normalize the number densities of the grain species
    fact = dens_grain/dens_ggtot
    yarray(b_specy+ind_ggminus) = dens_ggminus * fact
    yarray(b_specy+ind_ggplus)  = dens_ggplus  * fact
    yarray(b_specy+ind_gg0)     = dens_gg0     * fact
    !
    !--- recompute magnetosonic speed (cm.s-1) allowing for charged grains ---
!!$    vmagnet = gamma * kb * (ti+te) / &
!!$         ((densityi*mui+densityneg*muneg) / (densityi+densityneg)) + &
!!$         (bfield*vs_cm/vi)**2 / (4._dp*pi*(rhoi+rhoneg+rho_grain_charge))
    vmagnet = gamma*kb*(ti+te)*(densityi+densityneg)/(densityi*mui+densityneg*muneg)&
         + bfield2/(rhoi+rhoneg+rho_grain_charge)
    vmagnet = sqrt(vmagnet)
    !
    ! density of each H2 level
    ! JLB + GF (25 VI 2002) - Update H2_lev%Form_gr if required
    if (do_h2) then
       h2_lev(1:nh2_lev)%density = yarray(bv_h2_lev:ev_h2_lev)
       if (iforh2.eq.4) then
          h2_lev(1:nh2_lev)%form_gr = h2_lev(1:nh2_lev)%density / (dens_ph2+dens_oh2)
          h2_int_e = sum(h2_lev(1:nh2_lev)%form_gr * h2_lev(1:nh2_lev)%energy)
       endif
    endif
    !
    ! Calculate the ortho:para H2 ratio, and the densities of ortho-H2
    ! and para-H2 (cm-3).  results : op_H2, Dens_orthoH2, Dens_paraH2
    call compute_op_h2
    !
    ! Source terms of ([1], equations(3)-(14)) are calculated as
    ! functions of the chemistry and the abundances of the species and
    ! the global properties of the fluids.  results in :
    !    YN
    !    CupN,   CupI,   CupNEG
    !    CdownN, CdownI, CdownNEG
    !    YNn,    YNi,    YNneg
    !    Sn,     Si,     Sneg
    !    An,     Ai,     Aneg
    !    Bn,     Bi,     Bneg
    call evolution_chemistry
    !
    !
    ! Source terms for MHD equations are computed depending on the
    ! results of subroutine CHEMISTRY and the function values :
    ! yarray(i).  SOURCE also computes H2 cooling and internal energy.
    ! results in :
    !   An, Ai, Aneg
    !   Bn, Bi, Bneg
    !   H2_energy, YN_rovib_H2
    !--------------------------------------------------------------------
    ! Compute MHD source terms only if required
    if (do_shock) call evolution_physics
    if (do_h2) then
       ! Add chemical contribution YN (from CHEMISTRY) to the source of
       ! internal energy of H2 H2_energy (from COMPUTE_H2 via SOURCE)
       ! Sel_tot_H2 is the rate of destruction of H2 by collisional
       ! dissociation
       sel_tot_h2 = sum(dble(h2_lev(1:nh2_lev)%density * sel_ch_h2))
       sel_tne_h2 = sum(dble(h2_lev(1:nh2_lev)%density * sel_ne_h2))
       sel_tio_h2 = sum(dble(h2_lev(1:nh2_lev)%density * sel_io_h2))
       sel_tel_h2 = sum(dble(h2_lev(1:nh2_lev)%density * sel_el_h2))
       h2_int_energy = sum(dble(h2_lev(1:nh2_lev)%density * h2_lev(1:nh2_lev)%energy))
       h2_energy = h2_energy+kb*(yn(ind_ph2)+yn(ind_oh2)-sel_tot_h2)/sum_h2*h2_int_energy &
            + kb*sum(dble(h2_lev(1:nh2_lev)%density*sel_ch_h2*h2_lev(1:nh2_lev)%energy))
    endif

    !
    !-----------------------------------------------------
    ! calculation of the derivatives of the MHD variables
    ! with respect to Z (distance into the cloud)
    !-----------------------------------------------------
    !
    ! 29 mai 2001 - JLB : molecular cooling not included in Bn any more
    !                     but computed directly here. Old expression was:
    !                     (with molec_cool in Bn, via cooling_n)
    ! one fluid model : Bne = Bn + Bi + Bneg
    select case (nfluids)
    case (1)
       ! one-fluid model : bne = bn + bi + bneg; ane = 0; sne = 0
       bne = bn + bi + bneg
       ane = 0.0_dp
       sne = 0.0_dp
       ynne = ynn + yni + yni
       densityne = densityn + densityi + densityi ! PHB Mistake here ?????????
       rhone = rhon + rhoi + rhoneg
    case (2)
       ! two-fluids model : bne = bn
       bne = bn
       ane = an
       sne = sn
       ynne = ynn
       densityne = densityn
       rhone = rhon
    case (3)
       ! three-fluids model : bne = bn
       bne = bn
       ane = an
       sne = sn
       ynne = ynn
       densityne = densityn
       rhone = rhon
    end select
    !
    ! Evolution of the fractional radius: dz = Vn * dt
    ! dy(X) := dX/dz = 1/vn * dX/dt 
    ! 
    ! The collapse time is taken equal to the largest of the local
    ! free fall time and the sum of the two local ambipolar diffusion
    ! time scales:
    ! tau_collapse = max(tf*x**1.5d0,(tau_grain_n+tau_i_n))
    !
    ! Evolution of the radius: xsphere = r/r0
    ! dy(iv_xsphere) = -0.5d0*pi*(1./xsphere-1.)**0.5d0/tau_collapse/vn
    dy(iv_xsphere) = 0.d0 
    !
    !
    !
    ! Vgrad should be larger than 1 km s-1 pc-1 = 1 / parsec = 1 / 3.0857d13
    if (do_thermal) then
       ! Molecular cooling
       ! * first call with old value of Vgrad
       call compute_molecular_cooling
       aux1_dvn = gamma3*sne*vn2-h2_energy-gamma2*ane*vn+bne
       aux2_dvn = (densityne*gamma2*kb*tn-rhone*gamma1*vn2)
       if (shock_type.eq."J") then
          vgrad = sqrt(grad_v*grad_v + 1.0d10 / (parsec*parsec)) * 1.0d-5
          call compute_molecular_cooling
          dvn2 = (aux1_dvn - molec_cool) / aux2_dvn
          dy(iv_vn) = dvn2
          if (nfluids.eq.1) aux2_dvn = aux2_dvn + gamma1 * bfield2 !  Add Magnetic term for J shocks
          if (viscosity) dy(iv_vn) = -grad_v
       else
          dvn1 = (aux1_dvn-molec_cool) / aux2_dvn
          vgrad = sqrt(dvn1*dvn1 + 1.0d10 / (parsec*parsec)) * 1.0d-5
          call compute_molecular_cooling
          dvn2 = (aux1_dvn - molec_cool) / aux2_dvn
          gvn1 = dvn1 - dvn2
          ii = 0
          do while ( abs(gvn1) > 1.0e-20 .AND. ii < 1000)
             vgrad = sqrt(dvn1*dvn1 + 1.0d10 / (parsec*parsec)) * 1.0d-5
             call compute_molecular_cooling
             gvn2 = dvn2 - (aux1_dvn - molec_cool) / aux2_dvn
             dvn1 = dvn2
             dvn2 = dvn2 - gvn2 * gvn1 / (gvn1 - gvn2)
             gvn1 = gvn2
             ii = ii + 1
          end do
          dy(iv_vn) = dvn2
          vgrad = sqrt(dvn1*dvn1 + 1.0d10 / (parsec*parsec)) * 1.0d-5
       endif
    else
       !
       ! constant density model
       dy(iv_vn) = 0.0_dp
    endif
    !
    ! Vi, according to the number of fluids
    select case (nfluids)
    case (1)
       ! one-fluid model : vi=vn
       dy(iv_vi) = dy(iv_vn)
    case default
       ! two- or three- fluids model : ions and electrons have the same velocity
       dy(iv_vi) = (-gamma3*sn*vi2 + gamma2*an*vi + bi + bneg) &
            / (gamma2*kb*densityi*(ti+te) + gamma1*bfield2 - gamma1*vi2*(rhoi+rhoneg+rho_grain_charge))
    end select
    !
    ! Mass densities
    dilution      = 3*dy(iv_xsphere)/xsphere ! for free-fall
    dy(iv_rhon)   =   sn/vn - dilution*rhon
    dy(iv_rhoi)   =   si/vi - dilution*rhoi
    dy(iv_rhoneg) = sneg/vi - dilution*rhoneg
    !
    ! Neutral temperature
    if (do_thermal) then
       if (jshockvisc) then
          xll2 = xll*xll
          dy(iv_tn) = bne - molec_cool - ane * vn + 0.5_dp * sne * vn2 &
               - gamma1 * kb * tn * ynne - h2_energy &
               - densityne*kb*tn * dy(iv_vn) &
               - rhone * xll2 * dy(iv_vn) * dy(iv_vn) * dy(iv_vn)
          dy(iv_tn) = dy(iv_tn) / (gamma1 * kb * vn * densityne)
       else
          dy(iv_tn) = bne-molec_cool-ane*vn+0.5_dp*sne*vn2      &
               -gamma1*kb*tn*ynne-h2_energy                     &
               -densityne*kb*tn*vn*dilution
          dy(iv_tn) = dy(iv_tn)/(gamma1*kb*vn*densityne)
       endif
    else
       ! constant temperature model
       dy(iv_tn) = 0.0_dp
    endif
    !
    !--- Ti and Te, according to the number of fluids ---
    select case (nfluids)
    case (1)
       ! one fluid model : tn=ti=te
       dy(iv_ti) = dy(iv_tn)
       dy(iv_te) = dy(iv_tn)
    case (2)
       ! two fluids model : ti=te
       !      dy(iv_ti)=-0.25_dp*sne*vi2-gamma1*yni*kb*ti+0.5_dp*ane*vi
       dy(iv_ti) = -0.25_dp*sn*vi2-gamma1*yni*kb*ti+0.5_dp*an*vi
       dy(iv_ti) = dy(iv_ti)+0.5_dp*(bi+bneg)-densityi*kb*ti*dy(iv_vi)
       dy(iv_ti) = dy(iv_ti)/(gamma1*kb*v_densityi)
       dy(iv_te) = dy(iv_ti)
    case (3)
       ! three fluids model
       dy(iv_ti) = (0.5_dp*si*vi2-gamma1*yni*kb*ti-ai*vi +bi &
            - densityi*kb*ti*dy(iv_vi))/(gamma1*kb*v_densityi)
       dy(iv_te) = (0.5_dp*sneg*vi2-gamma1*yni*kb*te-aneg*vi + bneg &
            - densityi*kb*te*dy(iv_vi))/(gamma1*kb*v_densityi)
    end select
    !
    ! DensityN, DensityI
    dy(iv_densityn) = ynn/vn-densityn*dilution
    dy(iv_densityi) = yni/vi-densityi*dilution
    !
    ! Neutral velocity gradient (only for J shocks)
    dy(iv_gv) = 0._dp
    if (jshockvisc) then
       if (nfluids == 1) then
          dy(iv_gv) = bne-molec_cool-0.5_dp*sne*vn2-xll2*dy(iv_vn)**2*sne &
               -gamma2*kb*tn*ynne-h2_energy-rhone*vn2*dy(iv_vn)                  &
               -gamma2*densityne*vn*kb*dy(iv_tn)+bfield2*dy(iv_vn)
          dy(iv_gv) = dy(iv_gv)/(2.0_dp*rhone*xll2*dy(iv_vn)*vn)
       else
          dy(iv_gv) = bne-molec_cool-0.5_dp*sne*vn2-xll2*dy(iv_vn)**2*sne &
               -gamma2*kb*tn*ynne-h2_energy-rhone*vn2*dy(iv_vn)                  &
               -gamma2*densityne*vn*kb*dy(iv_tn)
          dy(iv_gv) = dy(iv_gv)/(2.0_dp*rhone*xll2*dy(iv_vn)*vn)
       end if
       ! iv_gv is for -dvn/dz
       dy(iv_gv) = - dy(iv_gv)
    end if
    !
    ! Grain coagulation
    ! Steady-state calculation assumes no coagulation; v_crit=0
    v_crit = 0d0
    if (v_crit.eq.0.d0) then
       v_rel  = 0d0
       fact   = 0.d0
    else
       v_turb = vsound/4.d0
       fact = 2._dp*log(2.)*(v_crit/v_turb)**2
       v_rel  = (2.d0/pi/log(2.))**0.5d0*v_turb*(1.-(1.+fact)*exp(-fact))
       fact   = -4.*pi*rgrain2*yarray(iv_dens_grain)**2.*v_rel/vn
    endif
    dy(iv_dens_grain)= -yarray(iv_dens_grain)*dilution + fact
    !
    ! calculates the derivatives of the chemical species's densities
    ! with respect to Z (distance into the cloud)
    dy(bv_neu:ev_neu) = yn(b_neu:e_neu)/vn - yarray(bv_neu:ev_neu)*dilution ! neutrals (velocity = vn)
    dy(bv_gra:ev_gra) = yn(b_gra:e_gra)/vi - yarray(bv_gra:ev_gra)*dilution ! species on grains (velocity = vi)
    dy(bv_cor:ev_cor) = yn(b_cor:e_cor)/vi - yarray(bv_cor:ev_cor)*dilution ! species in grain cores (velocity = vi)    
    dy(bv_ion:ev_ion) = yn(b_ion:e_ion)/vi - yarray(bv_ion:ev_ion)*dilution ! ions >0 (velocity = vi)
    dy(bv_neg:ev_neg) = yn(b_neg:e_neg)/vi - yarray(bv_neg:ev_neg)*dilution ! ions <0 (velocity = vi)
    !
    !  redefine the gradients of the grain species to allow for grain coagulation
    fact = dy(iv_dens_grain)/yarray(iv_dens_grain)
    dy(b_specy+ind_gg0)     = yn(ind_gg0)/vn     + fact*yarray(b_specy+ind_gg0)
    dy(b_specy+ind_ggplus)  = yn(ind_ggplus)/vn  + fact*yarray(b_specy+ind_ggplus)
    dy(b_specy+ind_ggminus) = yn(ind_ggminus)/vn + fact*yarray(b_specy+ind_ggminus)
    !
    if (do_h2) then
       ! Densities of ro-vibrational levels of H2
       !
       ! The contribution of chemical reactions YN(ind_H2) is included in proportion
       ! to the fractional population of the level : yarray(J)/sum_H2
       ! except contribution from selective reactions which no ponderation
       dy(bv_h2_lev:ev_h2_lev) = (yn_rovib_h2(1:nh2_lev) &
            + (yn(ind_ph2)+yn(ind_oh2) - sel_tot_h2 - for_gr_h2) * yarray(bv_h2_lev:ev_h2_lev) / sum_h2 &
            + h2_lev(1:nh2_lev)%density * sel_ch_h2 &
            + for_gr_h2 * h2_lev(1:nh2_lev)%form_gr &
            - yarray(bv_h2_lev:ev_h2_lev)*vn*dilution) / vn
    endif
    !
    !
    ! Logarithmic derivatives
    ! divide DY by Y = yarray
    ! * MHD variables
    ! * Chemical species
    ! * H2 levels : tiny avoids overflow if the function value is very small
    ! * Evaluate the dynamical timescale rho_n/(d rho_n/dt)
    tau_dyn = 1._dp / dy(iv_rhon) / vn
    !
    ! Add molec_cool to Bn as before, so it may be used without
    ! modifications in other parts of the code
    Bn = Bn-molec_cool
!    write(*,*) "PHB: diffun 7",nv_mhd,bv_speci,ev_speci,&
!         minval(abs(dy),dy.ne.0.d0),maxval(abs(dy)),minval(abs(yarray),yarray.ne.0.d0),maxval(abs(yarray))
    !
!    write(debug_ul_diffun_chem,'(300(es15.3e3))') y
    if (debug) then
       write(debug_ul_diffun_phys,'(300(es10.1))') dy(1:nv_mhd) !(bv_h2_lev:ev_h2_lev)
       write(debug_ul_diffun_chem,'(300(es10.1))') dy(bv_speci:ev_speci) !(bv_h2_lev:ev_h2_lev)
       if (do_h2) write(debug_ul_diffun_h2,  '(300(es10.1))') dy(bv_h2_lev:ev_h2_lev) !(bv_h2_lev:ev_h2_lev)
    endif
  end subroutine diffun

end module module_evolution
