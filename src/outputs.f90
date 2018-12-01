module module_outputs
  use module_precision
  use module_parameters_flags
  implicit none
  character(len=*), parameter :: name_file_info='output/info_mhd.out'
  character(len=*), parameter :: format_headercol='(a1,i14,500(i15))'
  character(len=*), parameter :: format_header='(a15)'
  character(len=*), parameter :: format_out='(es12.3e3,1x)'
!!$!   CHARACTER(LEN=*), PARAMETER :: format_lng='(ES17.8E3,1X)'
!!$    CHARACTER(LEN=*), PARAMETER :: format_lng='(ES18.9E3,1X)'
  character(len=*), parameter :: format_lng='(es15.6e3)'
  character(len=*), parameter :: name_file_phys='output/mhd_phys.out'
  character(len=*), parameter :: name_ad='output/adens_speci.out'
  character(len=*), parameter :: name_cd='output/cdens_speci.out'
  character(len=*), parameter :: name_fd='output/fdens_speci.out'
  character(len=*), parameter :: name_nr='output/rates_speci.out'
  character(len=*), parameter :: name_file_h2_lev='output/h2_lev.out'
  character(len=*), parameter :: name_file_h2_lin='output/h2_line.out'
  character(len=*), parameter :: name_file_cooling='output/cooling.out'
  character(len=*), parameter :: name_file_energetics='output/energetics.out'
  character(len=*), parameter :: name_file_intensity='output/intensity.out'
  character(len=*), parameter :: name_file_populations='output/populations.out'
  character(len=*), parameter :: name_file_fe_pops='output/fe_pops.out'
  character(len=*), parameter :: name_file_fe_lines='output/fe_lines.out'
contains
  subroutine write_welcome
    !_________________________________________________________________
    !
    ! Screen message
    !_________________________________________________________________
    !
    if (verbose) then
       write(*,'("!",70("_"))')
       write(*,'("!")')
       write(*,'("!",30x,"Collapse")')
       write(*,'("!",70("_"))')
       write(*,'("!")')
    endif
  end subroutine write_welcome


  subroutine write_open_file(ul,name)
    use module_parameters_flags
    use module_tools, only : get_file_number
    !_________________________________________________________________
    !
    ! Open an output file
    !_________________________________________________________________
    !
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(out) :: ul
    ul = get_file_number()
    open(ul,file=name,status='replace',access='sequential',form='formatted',action='write')
  end subroutine write_open_file


  subroutine write_info
    use module_constants
    use module_chem_react
    use module_phys_var
    use module_h2
    use module_grains
    use module_energetics
    use module_parameters_flags
    use module_parameters_output
    !_________________________________________________________________
    !
    ! Write informations about shock parameters, H2 levels, grains,
    ! chemical species and reactions in the file file_info
    !_________________________________________________________________
    !
    implicit none
    ! Local variables
    integer(kind=long),dimension(8) :: date
    integer(kind=long) :: i
    integer                     :: ul
    !
    ! file opening
    call write_open_file(ul,name_file_info)
    call date_and_time(values=date)
    write(ul,'("! date : ",i0,"/",i0,"/",i4,"  time: ",i0,"h",i0)')&
         date(3),date(2),date(1),date(5),date(6)
    write(ul,'("!",70("_"))')
    write(ul,'(/,"! Physical conditions")')
    write(ul,'(t11,"nh         [   cm-3] :",es9.2)')nh
    write(ul,'(t11,"tn         [      K] :",f9.1)') Tn
    write(ul,'(t11,"rad        [       ] :",f9.2)') rad
    write(ul,'(t11,"av         [    mag] :",f9.2)') av
    write(ul,'(t11,"zeta       [    s-1] :",es9.2)')zeta
    write(ul,'(/,"! Gas")')
    write(ul,'(t11,"rho(neut)  [ g.cm-3] :", es9.2)') rhon
    write(ul,'(t11,"rho(ions)  [ g.cm-3] :", es9.2)') rhoi
    write(ul,'(t11,"mass/H     [    H-1] :", es9.2)') rhon/nh/mp
    write(ul,'(/,"! Grains (grain/gas computed from core+mantle abundances)")')
    write(ul,'(t11,"grmin      [     cm] :", es9.2)') grmin
    write(ul,'(t11,"grmax      [     cm] :", es9.2)') grmax
    write(ul,'(t11,"rho core   [ g.cm-3] :", es9.2)') rho_grain
    write(ul,'(t11,"rho mantle [ g.cm-3] :", es9.2)') rho_mantle
    write(ul,'(t11,"tgrain     [      K] :",f9.1)')   tgrain
    write(ul,'(t11,"r_grain    [     cm] :",es9.2)')  r_grain
    write(ul,'(t11,"md_grain   [ g.cm-3] :",es9.2)')  md_grain
    write(ul,'(t11,"dens_grain [   cm-3] :",es9.2)')  dens_grain
    write(ul,'(t11,"grain/gas mass ratio :", es9.2)') ratio_grain_gas
    write(ul,'(/,"! Timescales (yr)")')
    write(ul,'(t11,"Collapse         :",es9.2)')  tau_collapse/yearsec
    write(ul,'(t11,"Free fall        :",es9.2)')  tau_ff/yearsec
    write(ul,'(t11,"i/n amb. diff.   :",es9.2)')  tau_i_n/yearsec
    write(ul,'(t11,"g/n              :",es9.2)')  tau_grain_n/yearsec
    write(ul,'(/,"! Dynamical parameters")')
    write(ul,'(t11,"masscloud (msol) :",es9.2)') masscloud/msun
    write(ul,'(t11,"maxdist (cm)     :",es9.2)') maxdist
    write(ul,'(t11,"maxtime (yr)     :",es9.2)') maxtime
    write(ul,'(t11,"vs (km.s-1)      :",f9.3)')  vs_km
    write(ul,'(t11,"o/p h2           :",es9.2)') op_h2
    write(ul,'(t11,"b (micro gauss)  :",es9.2)') bfield*1.d6
    write(ul,'(t11,"do_shock         :",l9)')    do_shock
    write(ul,'(t11,"Shock type       :",a9)')    shock_type
    write(ul,'(t11,"nfluids          :",i9)')    nfluids
    write(ul,'(t11,"xll (cm)         :",es9.2)') xll
    write(ul,'(/,"! Chemical network summary")')
    call write_info_indices(ul)
    write(ul,'(/,"! Integration parameters")')
    write(ul,'(t11,"nstep_max        :",i9)')    nstep_max
    write(ul,'(t11,"nstep_w          :",i9)')    nstep_w
    write(ul,'(t11,"abs_err          :",es9.2)') atol0
    write(ul,'(t11,"rel_err          :",es9.2)') rtol
    write(ul,'(t11,"tout             :",es9.2)') tout
    write(ul,'(t11,"tout_step        :",es9.2)') tout_step
    write(ul,'(/,"! Outputs")')
    write(ul,'(t11,"species (AD,CD,FD,RATES) : ",l2,l2,l2,l2)') &
         speci_out_ad,speci_out_cd,speci_out_fd,speci_out_nr
    write(ul,'(/,"! Thermal balance")')
    write(ul,'(t11,"flag             :",l9)')    do_thermal
    write(ul,'(t11,"cool_kn          :",i9)')    cool_kn
    if (shock_type=='C') then
       write(ul,'(" (not used here)")')
    else
       write(ul,*)
    endif
    write(ul,'(/,"! H2 molecule ")')
    write(ul,'(t11,"flag             : ",l7)')     do_h2
    write(ul,'(t11,"iforh2           : ", i7)')iforh2
    if (do_h2) then
       write(ul,'(t11,"number of levels : ",i7)')   nh2_lev
       write(ul,'(t11,"E(v=",i2,",j=",i2,") (k) : ",f7.1,"  (last level)")') &
            h2_lev(nh2_lev)%v, &
            h2_lev(nh2_lev)%j, &
            h2_lev(nh2_lev)%energy
       write(ul,'(t11,"vmax             : ", i7)')vmax_h2
       write(ul,'(t11,"jmax             : ", i7)')jmax_h2
       write(ul,'(t11,"H-H2 collisions  : ", a4)')h_h2_flag
       write(ul,'(t11,"ikinh2           : ", i7)')ikinh2
       write(ul,'(t11,"h2_int_e         : ", es9.2," K")')h2_int_e 
       write(ul,'(t11,"h2_lev_out       : ",a)') h2_lev_out
       write(ul,'(t11,"h2_lin_out       : ",a)') h2_lin_out
    endif
    !
    ! --- energetics ---
    write(ul,'(/,"! Energetics")')
    write(ul,'(t11,"flag                    :",l9)')    do_shock
    write(ul,'(t11,"mass flux (g/s/cm2)     :",es9.2)') mass_flux_init
    write(ul,'(t11,"momentum flux (erg/cm3) :",es9.2)') momentum_flux_init
    write(ul,'(t11,"energy flux (erg/s/cm2) :",es9.2)') energy_flux_init
    !
    ! --- elemental abundances ---
    write(ul,'(/,"! Elemental abundances (gas + mantles + pah)")')
    do i=1,nelements
       write(ul,'(t11,a," : ",es8.2, "   (ref : ",es8.2,")")') &
            elements(i)%name, dble(elements(i)%ab_init), dble(elements(i)%ab_ref)
    end do
    !
    ! --- chemical species ---
    write(ul,*)
    write(ul,'("!",70("_"))')
    write(ul,'("! ",i0," chemical species (+",i0," added)")')nspec, nspec_plus-nspec
    write(ul,'("! ",i0," neutrals, ",i0," on grains, &
         &", i0, " positive ions, ", i0," negative ions")') &
         nneutrals, nongrains, nions, nneg
    write(ul,'("!NUM NAME      COMPOSITION   DENSITY      ENTHALPY (-99.999=unknown)")')
    write(ul,'("!                             [cm-3]       [kCal mol-1]")')
    write(ul,'("!",70("_"),/,"!")')
    do i=1, nspec
       call write_specy(ul,speci(i))
    enddo
    do i=nspec+1,nspec_plus
       write(ul,'(i3,2x,a7)') speci(i)%index,speci(i)%name
    end do

    ! --- chemical reactions ---
    write(ul,*)
    write(ul,'("!",70("_"))')
    write(ul,'("! ",i0," chemical reactions including:")') nreact
    write(ul,'("!   ",i0," photo, ",i0," cr_io, ", &
         &i0," cr_de, ",i0," h2_fo, ",i0," three, ",i0," sputt")') &
         &nphoto, ncr_io, ncr_de, nh2_fo, nthree, nsputt
    write(ul,'("!   ",i0," erosi, ",i0," adsor, ",i0,&
         &" disso, ",i0," other, and ",i0," rever")') nerosi, nadsor, ndisso, nother, nrever
    write(ul,'("!")')
    write(ul,&
         '("!R1     R2        P1      P2      P3      P4         &
         GAMMA       ALPHA       BETA        DE        TYPE  NUM!  REF")')
    write(ul,&
         '("!                                                    &
         [cm3 s-1]   []          [K]         [ev]")')
    write(ul,'("!",70("_"),/,"!")')
    do i=1,nreact
       call write_reaction(ul,i,react(i))
    end do
    !
    ! file closing
    close(ul)
    write(*,*) "I-write_info: done"
  end subroutine write_info


  subroutine write_info_indices(lunit)
    use module_parameters_output
    use module_chem_react
    use module_phys_var
    !---------------------------------------------------------------------------
    !
    ! Write basic informations on terminal output during integration.
    !
    !---------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: lunit
    !
    ! On terminal output:
    write(lunit,*)
    write(lunit,'(3x,(a),3i5)') "Total size of Y vector:         ",dimtot
    write(lunit,'(3x,(a),3i5)') "Number of physical variables:   ",nv_mhd,1,nv_mhd
    write(lunit,'(3x,(a),3i5)') "Number of chemical variables:   ",ev_speci-bv_speci+1,bv_speci,ev_speci
    write(lunit,'(3x,(a),3i5)') "Neutrals (total, indices):      ",e_neu-b_neu+1,bv_neu,ev_neu
    write(lunit,'(3x,(a),3i5)') "Grain mantles (total, indices): ",e_gra-b_gra+1,bv_gra,ev_gra
    write(lunit,'(3x,(a),3i5)') "Grain cores (total, indices):   ",e_cor-b_cor+1,bv_cor,ev_cor
    write(lunit,'(3x,(a),3i5)') "Positive ions (total, indices): ",e_ion-b_ion+1,bv_ion,ev_ion
    write(lunit,'(3x,(a),3i5)') "Negative ions (total, indices): ",e_neg-b_neg+1,bv_neg,ev_neg
    write(lunit,'(3x,(a),3i5)') "Number of H2 levels:            ",ev_h2_lev-bv_h2_lev+1,bv_h2_lev,ev_h2_lev
  end subroutine write_info_indices


  subroutine write_screen
    use module_parameters_output
    use module_phys_var
    use module_h2
    !---------------------------------------------------------------------------
    !
    ! Write basic informations on terminal output during integration.
    !
    !---------------------------------------------------------------------------
    implicit none
    !
    write(screen,'("step= ",i5, 2x)', advance='no') counter
    write(screen,'("z=", es9.1, 2x)', advance='no') distance
    write(screen,'("Tout=", es9.1, 2x)', advance='no') tout
    write(screen,'("timen=", es9.1, 2x)', advance='no') timen
    write(screen,'("Tn=", es9.1, 2x)', advance='no') tn
    write(screen,'("nh=", es9.1, 2x)', advance='no') nh
    write(screen,'("Dv=", es10.1, 2x)', advance='no') deltav
    write(screen,'("R_h2=",es10.1)', advance='yes') r_gr_h2
  end subroutine write_screen


  subroutine write_output(again)
    use module_tools, only : get_file_number
    use module_phys_var
    use module_constants
    use module_grains
    use h2_variables
    use module_h2, only: h2_lines,h2_lev
    use module_chemical_species, only : nspec,speci,ind_h,ind_hplus,ind_e,dens_h2
    use module_molecular_cooling
    use module_energetics
    use module_line_excit
    use module_parameters_flags
    use module_parameters_output
    !---------------------------------------------------------------------------
    !
    ! Open, write, and close output files
    !
    !---------------------------------------------------------------------------
    implicit none
    logical, intent(in) :: again ! FALSE=last calculation
    ! Local variables
    integer(kind=long) :: i
    character(len=15) :: str,str2
    character(len=2) :: str0
    character(len=64) :: fmt
    character(len=*), parameter :: rname='write_output'
    !
    ! On first call, open the different files and write the headers
    if (counter.eq.1) then
       call write_open_file(file_phys,name_file_phys)
       if (speci_out_ad) call write_open_file(file_spec_ad,name_ad)
       if (speci_out_cd) call write_open_file(file_spec_cd,name_cd)
       if (speci_out_fd) call write_open_file(file_spec_fd,name_fd)
       if (speci_out_nr) call write_open_file(file_spec_nr,name_nr)
       if (do_h2) then
          call write_open_file(file_h2_lev,name_file_h2_lev)
          call write_open_file(file_h2_lin,name_file_h2_lin)
       endif
       if (do_thermal) then
          call write_open_file(file_cooling,name_file_cooling)
          call write_open_file(file_energetics,name_file_energetics)
          call write_open_file(file_intensity,name_file_intensity)
          call write_open_file(file_populations,name_file_populations)
          call write_open_file(file_fe_pops,name_file_fe_pops)
          call write_open_file(file_fe_lines,name_file_fe_lines)
       endif
    endif

    IF (counter == 1) THEN
       !
       ! MHD variables
       write(file_phys,format_headercol) comment,(/(i,i=1,50)/)
       write(file_phys,'(a1)',advance='no') comment
       write(file_phys,'(a14)',advance='no')'distance'
       write(file_phys,format_header,advance='no')'timeN'
       write(file_phys,format_header,advance='no')'timeI'
       write(file_phys,format_header,advance='no')'t_dyn'
       write(file_phys,format_header,advance='no')'t_ff'
       write(file_phys,format_header,advance='no')'t_i_n'
       write(file_phys,format_header,advance='no')'t_grain_n'
       write(file_phys,format_header,advance='no')'xsphere'
       write(file_phys,format_header,advance='no')'nH'
       write(file_phys,format_header,advance='no')'op_H2'
       write(file_phys,format_header,advance='no')'Vn'
       write(file_phys,format_header,advance='no')'Vi'
       write(file_phys,format_header,advance='no')'Tn'
       write(file_phys,format_header,advance='no')'Ti'
       write(file_phys,format_header,advance='no')'Te'
       write(file_phys,format_header,advance='no')'Vsound'
       write(file_phys,format_header,advance='no')'Vmagnet'
       write(file_phys,format_header,advance='no')'rhoN'
       write(file_phys,format_header,advance='no')'rhoI'
       write(file_phys,format_header,advance='no')'rhoNEG'
       write(file_phys,format_header,advance='no')'DensityN'
       write(file_phys,format_header,advance='no')'DensityI'
       write(file_phys,format_header,advance='no')'DensityNEG'
       write(file_phys,format_header,advance='no')'Densitye'
       write(file_phys,format_header,advance='no')'muN'
       write(file_phys,format_header,advance='no')'muI'
       write(file_phys,format_header,advance='no')'muNEG'
       write(file_phys,format_header,advance='no')'-dVn/dz'
       write(file_phys,format_header,advance='no')'-dVi/dz'
       ! grains
       write(file_phys,format_header,advance='no')'Tgrain'
       write(file_phys,format_header,advance='no')'Teff_gr'
       write(file_phys,format_header,advance='no')'Nlayers_gr'
       write(file_phys,format_header,advance='no')'n(grain)'
       write(file_phys,format_header,advance='no')'MDens_gr'
       write(file_phys,format_header,advance='no')'MDens_gr_charge'
       write(file_phys,format_header,advance='no')'Rgrain'
       write(file_phys,format_header,advance='no')'-DissH2_t'
       write(file_phys,format_header,advance='no')'-DissH2_n'
       write(file_phys,format_header,advance='no')'-DissH2_i'
       write(file_phys,format_header,advance='no')'-DissH2_e'
       write(file_phys,format_header,advance='no')'F_gr_H2'
       write(file_phys,format_header,advance='no')'E_gr_H2_K'
       write(file_phys,format_header,advance='no')'x(H2)'
       write(file_phys,format_header,advance='no')'x(H)'
       write(file_phys,format_header,advance='no')'x(H+)'
       write(file_phys,format_header,advance='no')'sum(neut)'
       write(file_phys,format_header,advance='no')'sum(ions)'
       write(file_phys,format_header,advance='no')'sum(nega)'
       write(file_phys,format_header,advance='no')'grad_V'               !  Used only for J shocks
       write(file_phys,*)
       !
       ! Chemical species
       ! densities (cm-3) and/or column densities (cm-2) and/or fractional abundances
       write(fmt,'(a,i0,a,a)') '(',4+nspec,format_header,')'
       if (speci_out_ad) then
          write(file_spec_ad,'(a,1x,a)') comment,"Number density (cm-3)"
          write(file_spec_ad,format_headercol) comment,(/(i,i=1,nspec+5)/)
          write(file_spec_ad,'(a1,a14)',advance='no') comment,'distance'
          write(file_spec_ad,fmt) 'timeN','Tn','nH','e-',(/(adjustr(speci(i)%name),i=1,nspec)/)
       endif
       if (speci_out_cd) then
          write(file_spec_cd,'(a,1x,a)') comment,"Column density (cm-2)"
          write(file_spec_cd,format_headercol) comment,(/(i,i=1,nspec+5)/)
          write(file_spec_cd,'(a1,a14)',advance='no') comment,'distance'
          write(file_spec_cd,fmt) 'timeN','Tn','nH','e-',(/(adjustr(speci(i)%name),i=1,nspec)/)
       endif
       if (speci_out_fd) then
          write(file_spec_fd,'(a,1x,a)') comment,"Fractional density (n(X)/nH)"
          write(file_spec_fd,format_headercol) comment,(/(i,i=1,nspec+5)/)
          write(file_spec_fd,'(a1,a14)',advance='no') comment,'distance'
          write(file_spec_fd,fmt) 'timeN','Tn','nH','e-',(/(adjustr(speci(i)%name),i=1,nspec)/)
       endif
       if (speci_out_nr) then
          write(file_spec_nr,format_headercol) comment,(/(i,i=1,nspec+5)/)
          write(file_spec_nr,'(a1,a14)',advance='no') comment,'distance'
          write(file_spec_nr,fmt) 'timeN','Tn','nH','e-',(/(adjustr(speci(i)%name),i=1,nspec)/)
       endif
       !
       ! H2 levels and H2 lines
       ! populations
       if (do_h2) then
          write(file_h2_lev,'(a14)',advance='no')'distance'
          write(file_h2_lev,format_header,advance='no')'timeN'
          write(file_h2_lev,format_header,advance='no')'Tn'
          ! densities (cm-3)
          do i=1,nh2_lev
             write(str,'(i2)') h2_lev(i)%v
             write(str2,'(i2)') h2_lev(i)%j
             write(file_h2_lev,format_header,advance='no') '(' // trim(adjustl(str)) // &
                  ',' // trim(adjustl(str2)) // ')'
          end do
          write(file_h2_lev,*)
          ! intensities (erg/s/cm2/sr)
          write(file_h2_lin,'(a14)',advance='no')'distance'
          write(file_h2_lin,format_header,advance='no')'timen'
          write(file_h2_lin,format_header,advance='no')'tn'
          do i=1,nh2_lines_out
             write(file_h2_lin,format_header,advance='no')trim(h2_lines(i)%name)
          end do
          write(file_H2_lin,*)
       endif
       !
       ! Molecular + atomic cooling
       ! only emissivities (erg/s/cm3) ---
       if (do_thermal) then
          write(file_cooling,'(a14)',advance='no')'distance'
          write(file_cooling,format_header,advance='no')'Timen'
          write(file_cooling,format_header,advance='no')'Tn'
          write(file_cooling,format_header,advance='no')'E(H2)'
          write(file_cooling,format_header,advance='no')'E(13CO)'
          write(file_cooling,format_header,advance='no')'E(OH)'
          write(file_cooling,format_header,advance='no')'E(NH3)'
          write(file_cooling,format_header,advance='no')'E(r,CO)'
          write(file_cooling,format_header,advance='no')'E(v,CO)'
          write(file_cooling,format_header,advance='no')'E(CO)'
          write(file_cooling,format_header,advance='no')'E(r,o-H2O)'
          write(file_cooling,format_header,advance='no')'E(r,p-H2O)'
          write(file_cooling,format_header,advance='no')'E(v,H2O)'
          write(file_cooling,format_header,advance='no')'E(H2O)'
          write(file_cooling,format_header,advance='no')'E(C+)'
          write(file_cooling,format_header,advance='no')'E(Si+)'
          write(file_cooling,format_header,advance='no')'E(C)'
          write(file_cooling,format_header,advance='no')'E(Si)'
          write(file_cooling,format_header,advance='no')'E(O)'
          write(file_cooling,format_header,advance='no')'E(S+)'
          write(file_cooling,format_header,advance='no')'E(N+)'
          write(file_cooling,*)
          !
          ! Energetics
          write(file_energetics,'(a14)',advance='no')'distance'
          write(file_energetics,format_header,advance='no')'Timen'
          write(file_energetics,format_header,advance='no')'Tn'
          write(file_energetics,format_header,advance='no')'Mass'
          write(file_energetics,format_header,advance='no')'Mass_cons'
          write(file_energetics,format_header,advance='no')'Moment_kin'
          write(file_energetics,format_header,advance='no')'Moment_the'
          write(file_energetics,format_header,advance='no')'Moment_mag'
          write(file_energetics,format_header,advance='no')'Moment_vis'
          write(file_energetics,format_header,advance='no')'Moment'
          write(file_energetics,format_header,advance='no')'Mom_cons'
          write(file_energetics,format_header,advance='no')'Energy_kin'
          write(file_energetics,format_header,advance='no')'Energy_the'
          write(file_energetics,format_header,advance='no')'Energy_int'
          write(file_energetics,format_header,advance='no')'Energy_mag'
          write(file_energetics,format_header,advance='no')'Energy_vis'
          write(file_energetics,format_header,advance='no')'Energy'
          write(file_energetics,format_header,advance='no')'Energ_gain'
          write(file_energetics,format_header,advance='no')'Energ_cons'
          write(file_energetics,*)
          !
          ! Intensities---
          write(file_intensity,'(a14)',advance='no')'distance'
          write(file_intensity,format_header,advance='no')'Timen'
          write(file_intensity,format_header,advance='no')'Tn'
          write(file_intensity,format_header,advance='no') 'C+(158m)'                   ! 'C+(2-1)'
          write(file_intensity,format_header,advance='no') 'C+(2324.7A)'                ! 'C+(3-1)'
          write(file_intensity,format_header,advance='no') 'C+(2323.5A)'                ! 'C+(4-1)'
          write(file_intensity,format_header,advance='no') 'C+(2328.1A)'                ! 'C+(3-2)'
          write(file_intensity,format_header,advance='no') 'C+(2326.9A)'                ! 'C+(4-2)'
          write(file_intensity,format_header,advance='no') 'C+(2325.4A)'                ! 'C+(5-2)'
          write(file_intensity,format_header,advance='no') 'Si+(34.8m)'                 ! 'Si+(2-1)'
          write(file_intensity,format_header,advance='no') 'C(609.8m)'                  ! 'C(2-1)'
          write(file_intensity,format_header,advance='no') 'C(370.4m)'                  ! 'C(3-2)'
          write(file_intensity,format_header,advance='no') 'C(9850A)'                   ! 'C(4-3)'
          write(file_intensity,format_header,advance='no') 'C(9824A)'                   ! 'C(4-2)'
          write(file_intensity,format_header,advance='no') 'Si(129.7m)'                 ! 'Si(2-1)'
          write(file_intensity,format_header,advance='no') 'Si(68.5m)'                  ! 'Si(3-2)'
          write(file_intensity,format_header,advance='no') 'O(63.2m)'                   ! 'O(2-1)'
          write(file_intensity,format_header,advance='no') 'O(145.3m)'                  ! 'O(3-2)'
          write(file_intensity,format_header,advance='no') 'O(6300A)'                   ! 'O(4-1)'
          write(file_intensity,format_header,advance='no') 'O(6363A)'                   ! 'O(4-2)'
          write(file_intensity,format_header,advance='no') 'S+(6731A)'                  ! 'S+(2-1)'
          write(file_intensity,format_header,advance='no') 'S+(6716A)'                  ! 'S+(3-1)'
          write(file_intensity,format_header,advance='no') 'N+(205.3m)'                 ! 'N+(2-1)'
          write(file_intensity,format_header,advance='no') 'N+(121.8m)'                 ! 'N+(3-2)'
          write(file_intensity,format_header,advance='no') 'N+(6527A)'                  ! 'N+(4-1)'
          write(file_intensity,format_header,advance='no') 'N+(6548A)'                  ! 'N+(4-2)'
          write(file_intensity,format_header,advance='no') 'N+(6583A)'                  ! 'N+(4-3)'
          write(file_intensity,format_header,advance='no') 'N(5200A)'                   ! 'N(2-1)'
          write(file_intensity,format_header,advance='no') 'N(5197A)'                   ! 'N(3-1)'
          write(file_intensity,*)
          !
          ! Populations
          write(file_populations,'(a14)',advance='no')'distance'
          write(file_populations,format_header,advance='no')'Timen'
          write(file_populations,format_header,advance='no')'Tn'
          write(file_populations,format_header,advance='no') 'C(3P-J=0)'             ! C, lev = 1
          write(file_populations,format_header,advance='no') 'C(3P-J=1)'             ! C, lev = 2
          write(file_populations,format_header,advance='no') 'C(3P-J=2)'             ! C, lev = 3
          write(file_populations,format_header,advance='no') 'C(1D-J=2)'             ! C, lev = 4
          write(file_populations,format_header,advance='no') 'C(1S-J=0)'             ! C, lev = 5
          write(file_populations,format_header,advance='no') 'N(4S-J=3/2)'           ! N, lev = 1
          write(file_populations,format_header,advance='no') 'N(2D-J=5/2)'           ! N, lev = 2
          write(file_populations,format_header,advance='no') 'N(2D-J=3/2)'           ! N, lev = 3
          write(file_populations,format_header,advance='no') 'N(2P-J=1/2)'           ! N, lev = 4
          write(file_populations,format_header,advance='no') 'N(2P-J=3/2)'           ! N, lev = 5
          write(file_populations,format_header,advance='no') 'O(3P-J=2)'             ! O, lev = 1
          write(file_populations,format_header,advance='no') 'O(3P-J=1)'             ! O, lev = 2
          write(file_populations,format_header,advance='no') 'O(3P-J=0)'             ! O, lev = 3
          write(file_populations,format_header,advance='no') 'O(1D-J=2)'             ! O, lev = 4
          write(file_populations,format_header,advance='no') 'O(1S-J=0)'             ! O, lev = 5
          write(file_populations,format_header,advance='no') 'S(3P-J=2)'             ! S, lev = 1
          write(file_populations,format_header,advance='no') 'S(3P-J=1)'             ! S, lev = 2
          write(file_populations,format_header,advance='no') 'S(3P-J=0)'             ! S, lev = 3
          write(file_populations,format_header,advance='no') 'S(1D-J=2)'             ! S, lev = 4
          write(file_populations,format_header,advance='no') 'S(1S-J=0)'             ! S, lev = 5
          write(file_populations,format_header,advance='no') 'Si(3P-J=0)'            ! Si, lev = 1
          write(file_populations,format_header,advance='no') 'Si(3P-J=1)'            ! Si, lev = 2
          write(file_populations,format_header,advance='no') 'Si(3P-J=2)'            ! Si, lev = 3
          write(file_populations,format_header,advance='no') 'Si(1D-J=2)'            ! Si, lev = 4
          write(file_populations,format_header,advance='no') 'Si(1S-J=0)'            ! Si, lev = 5
          write(file_populations,format_header,advance='no') 'C+(2P-J=1/2)'          ! C+, lev = 1
          write(file_populations,format_header,advance='no') 'C+(2P-J=3/2)'          ! C+, lev = 2
          write(file_populations,format_header,advance='no') 'C+(4P-J=1/2)'          ! C+, lev = 3
          write(file_populations,format_header,advance='no') 'C+(4P-J=3/2)'          ! C+, lev = 4
          write(file_populations,format_header,advance='no') 'C+(4P-J=5/2)'          ! C+, lev = 5
          write(file_populations,format_header,advance='no') 'N+(3P-J=0)'            ! N+, lev = 1
          write(file_populations,format_header,advance='no') 'N+(3P-J=1)'            ! N+, lev = 2
          write(file_populations,format_header,advance='no') 'N+(3P-J=2)'            ! N+, lev = 3
          write(file_populations,format_header,advance='no') 'N+(1D-J=2)'            ! N+, lev = 4
          write(file_populations,format_header,advance='no') 'N+(1S-J=0)'            ! N+, lev = 5
          write(file_populations,format_header,advance='no') 'O+(4S-J=3/2)'          ! O+, lev = 1
          write(file_populations,format_header,advance='no') 'O+(2D-J=5/2)'          ! O+, lev = 2
          write(file_populations,format_header,advance='no') 'O+(2D-J=3/2)'          ! O+, lev = 3
          write(file_populations,format_header,advance='no') 'O+(2P-J=3/2)'          ! O+, lev = 4
          write(file_populations,format_header,advance='no') 'O+(2P-J=1/2)'          ! O+, lev = 5
          write(file_populations,format_header,advance='no') 'S+(4S-J=3/2)'          ! S+, lev = 1
          write(file_populations,format_header,advance='no') 'S+(2D-J=3/2)'          ! S+, lev = 2
          write(file_populations,format_header,advance='no') 'S+(2D-J=5/2)'          ! S+, lev = 3
          write(file_populations,format_header,advance='no') 'S+(2P-J=1/2)'          ! S+, lev = 4
          write(file_populations,format_header,advance='no') 'S+(2P-J=3/2)'          ! S+, lev = 5
          write(file_populations,format_header,advance='no') 'Si+(2P-J=1/2)'         ! Si+, lev = 1
          write(file_populations,format_header,advance='no') 'Si+(2P-J=3/2)'         ! Si+, lev = 2
          write(file_populations,format_header,advance='no') 'Si+(4P-J=1/2)'         ! Si+, lev = 3
          write(file_populations,format_header,advance='no') 'Si+(4P-J=3/2)'         ! Si+, lev = 4
          write(file_populations,format_header,advance='no') 'Si+(4P-J=5/2)'         ! Si+, lev = 5
          write(file_populations,*)
          !
          ! Fe+ level populations
          ! Relative to ground state
          write(file_fe_pops,'(a14)',advance='no')'distance'
          write(file_fe_pops, format_header, advance='no')'timeN'
          write(file_fe_pops, format_header, advance='no')'Tn'
          write(file_fe_pops, format_header, advance='no')'a6D-J9'
          write(file_fe_pops, format_header, advance='no')'a6D-J7'
          write(file_fe_pops, format_header, advance='no')'a6D-J5'
          write(file_fe_pops, format_header, advance='no')'a6D-J3'
          write(file_fe_pops, format_header, advance='no')'a6D-J1'
          write(file_fe_pops, format_header, advance='no')'a4F-J9'
          write(file_fe_pops, format_header, advance='no')'a4F-J7'
          write(file_fe_pops, format_header, advance='no')'a4F-J5'
          write(file_fe_pops, format_header, advance='no')'a4F-J3'
          write(file_fe_pops, format_header, advance='no')'a4D-J7'
          write(file_fe_pops, format_header, advance='no')'a4D-J5'
          write(file_fe_pops, format_header, advance='no')'a4D-J3'
          write(file_fe_pops, format_header, advance='no')'a4D-J1'
          write(file_fe_pops, format_header, advance='no')'a4P-J5'
          write(file_fe_pops, format_header, advance='no')'a4P-J3'
          write(file_fe_pops, format_header, advance='no')'a4P-J1'
          write(file_fe_pops, format_header, advance='no')'b4P-J5'
          write(file_fe_pops, format_header, advance='no')'b4P-J3'
          write(file_fe_pops, format_header, advance='no')'b4P-J1'
          write(file_fe_pops,*)
          !
          ! [FeII] lines
          ! Intensities (erg//s/cm2/sr)
          write(file_fe_lines,'(a14)',advance='no')'distance'
          write(file_fe_lines, format_header, advance='no')'timeN'
          write(file_fe_lines, format_header, advance='no')'Tn'
          write(file_fe_lines, format_header, advance='no')'1.644'
          write(file_fe_lines, format_header, advance='no')'1.275'
          write(file_fe_lines, format_header, advance='no')'1.810'
          write(file_fe_lines, format_header, advance='no')'1.677'
          write(file_fe_lines, format_header, advance='no')'1.321'
          write(file_fe_lines, format_header, advance='no')'1.249'
          write(file_fe_lines, format_header, advance='no')'1.295'
          write(file_fe_lines, format_header, advance='no')'1.328'
          write(file_fe_lines, format_header, advance='no')'1.543'
          write(file_fe_lines, format_header, advance='no')'1.800'
          write(file_fe_lines, format_header, advance='no')'1.248'
          write(file_fe_lines, format_header, advance='no')'1.279'
          write(file_fe_lines, format_header, advance='no')'1.298'
          write(file_fe_lines, format_header, advance='no')'1.600'
          write(file_fe_lines, format_header, advance='no')'1.271'
          write(file_fe_lines, format_header, advance='no')'1.664'
          write(file_fe_lines,*)
       endif
       !
       ! End of output file initialization
    end if
    !
    ! MHD variables + grains
    ! mhd variables
    write(file_phys,format_lng,advance='no')distance
    write(file_phys,format_lng,advance='no')timeN
    write(file_phys,format_lng,advance='no')timeI
    write(file_phys,format_lng,advance='no')tau_dyn/YEARsec
    write(file_phys,format_lng,advance='no')tau_ff/YEARsec
    !   write(file_phys,format_lng,advance='no')tf*xsphere**1.5d0/YEARsec
    write(file_phys,format_lng,advance='no')tau_i_n/YEARsec
    write(file_phys,format_lng,advance='no')tau_grain_n/YEARsec
    write(file_phys,format_lng,advance='no')xsphere
    write(file_phys,format_lng,advance='no')nH
    write(file_phys,format_lng,advance='no')op_H2
    write(file_phys,format_lng,advance='no')Vn
    write(file_phys,format_lng,advance='no')Vi
    write(file_phys,format_lng,advance='no')Tn
    write(file_phys,format_lng,advance='no')Ti
    write(file_phys,format_lng,advance='no')Te
    write(file_phys,format_lng,advance='no')Vsound
    write(file_phys,format_lng,advance='no')Vmagnet
    write(file_phys,format_lng,advance='no')rhoN
    write(file_phys,format_lng,advance='no')rhoI
    write(file_phys,format_lng,advance='no')rhoNEG
    write(file_phys,format_lng,advance='no')DensityN
    write(file_phys,format_lng,advance='no')DensityI
    write(file_phys,format_lng,advance='no')DensityNEG
    write(file_phys,format_lng,advance='no')(DensityI-DensityNEG)/nH
    write(file_phys,format_lng,advance='no')muN
    write(file_phys,format_lng,advance='no')muI
    write(file_phys,format_lng,advance='no')muNEG
    write(file_phys,format_lng,advance='no')-dVn
    write(file_phys,format_lng,advance='no')-dVi
    ! grains
    write(file_phys,format_lng,advance='no')Tgrain
    write(file_phys,format_lng,advance='no')Teff_grain
    write(file_phys,format_lng,advance='no')Nlayers_grain
    write(file_phys,format_lng,advance='no')Dens_grain
    write(file_phys,format_lng,advance='no')MD_grain
    write(file_phys,format_lng,advance='no')Rho_GRAIN_charge
    write(file_phys,format_lng,advance='no')r_grain
    write(file_phys,format_lng,advance='no')-Sel_tot_H2
    write(file_phys,format_lng,advance='no')-Sel_tne_H2
    write(file_phys,format_lng,advance='no')-Sel_tio_H2
    write(file_phys,format_lng,advance='no')-Sel_tel_H2
    write(file_phys,format_lng,advance='no')For_gr_H2
    write(file_phys,format_lng,advance='no')H2_int_E
    write(file_phys,format_lng,advance='no')dens_h2!speci(ind_H2)%Density/nH
    write(file_phys,format_lng,advance='no')speci(ind_H)%Density/nH
    write(file_phys,format_lng,advance='no')speci(ind_Hplus)%Density/nH
    write(file_phys,format_lng,advance='no')SUM(yarray(bv_neu:ev_neu))
    write(file_phys,format_lng,advance='no')SUM(yarray(bv_ion:ev_ion))
    write(file_phys,format_lng,advance='no')SUM(yarray(bv_neg:ev_neg))
    if (shock_type == "J") then
       write(file_phys,format_lng,advance='no')grad_V
    else
       write(file_phys,format_lng,advance='no')-dVn
    end if
    write(file_phys,*)
    !
    ! chemical species
    write(fmt,'(a,i4,a,a)') '(',5+nspec,format_lng,')'
    if (speci_out_ad) then
       write(file_spec_ad,fmt,advance='no') &
            distance,timen,tn,nh,speci(ind_e)%density,(/(speci(i)%density,i=1,nspec)/)
       write(file_spec_ad,*)
    endif
    if (speci_out_cd) then
       write(file_spec_cd,fmt,advance='no') &
            distance,timen,tn,nh,speci(ind_e)%col_dens,(/(speci(i)%col_dens,i=1,nspec)/)
       write(file_spec_cd,*)
    endif
    if (speci_out_fd) then
       write(file_spec_fd,fmt,advance='no') &
            distance,timen,tn,nh,speci(ind_e)%density/nh,(/(speci(i)%density/nh,i=1,nspec)/)
       write(file_spec_fd,*)
    endif
    if (speci_out_nr) then
       ! Net creation rates
       write(file_spec_nr,fmt,advance='no') &
            distance,timen,tn,nh,speci(ind_e)%netrate,(/(speci(i)%netrate,i=1,nspec)/)
       write(file_spec_nr,*)
    endif
    !
    ! H2 levels + H2 lines
    if (do_h2) then
       write(file_h2_lev,format_lng,advance='no') distance
       write(file_h2_lev,format_lng,advance='no') timen
       write(file_h2_lev,format_lng,advance='no') tn
       select case(trim(h2_lev_out))
       case('AD')
          ! densities (cm-3)
          do i=1,nh2_lev
             write(file_h2_lev,format_lng,advance='no')h2_lev(i)%density
          end do
       case('CD')
          ! column densities (cm-2)
          do i=1,nh2_lev
             write(file_h2_lev,format_lng,advance='no')h2_lev(i)%col_dens
          end do
       CASE('ln(N/g)')
          ! excitation diagram
          do i=1,nh2_lev
             write(file_h2_lev,format_lng,advance='no')log(h2_lev(i)%col_dens/h2_lev(i)%weight)
          end do
       case default
          write(*,*) rname//", wrong output specification"
          stop
       end select
       write(file_h2_lev,*)
       write(file_h2_lin,format_lng,advance='no') distance
       write(file_h2_lin,format_lng,advance='no') timen
       write(file_h2_lin,format_lng,advance='no') tn
       select case(trim(h2_lin_out))
       case('local')
          ! emissivities (erg/s/cm-3)
          do i=1,nh2_lines_out
             write(file_h2_lin,format_lng,advance='no')h2_lines(i)%emiss
          end do
       case('integrated')
          ! intensities (erg/s/cm2/sr)
          do i=1,nh2_lines_out
             write(file_h2_lin,format_lng,advance='no')h2_lines(i)%intensity
          end do
       case default
          write(*,*) rname//", wrong output specification"
          stop
       end select
       write(file_h2_lin,*)
    endif
    !
    ! molecular and atomic cooling rates (erg/s/cm3)
    if (do_thermal) then
       write(file_cooling,format_lng,advance='no')distance
       write(file_cooling,format_lng,advance='no')timen
       write(file_cooling,format_lng,advance='no')tn
       write(file_cooling,format_lng,advance='no')cooling_h2
       write(file_cooling,format_lng,advance='no')cooling_13co
       write(file_cooling,format_lng,advance='no')cooling_oh
       write(file_cooling,format_lng,advance='no')cooling_nh3
       write(file_cooling,format_lng,advance='no')cooling_rot_co
       write(file_cooling,format_lng,advance='no')cooling_vib_co
       write(file_cooling,format_lng,advance='no')cooling_co
       write(file_cooling,format_lng,advance='no')cooling_rot_o_h2o
       write(file_cooling,format_lng,advance='no')cooling_rot_p_h2o
       write(file_cooling,format_lng,advance='no')cooling_vib_h2o
       write(file_cooling,format_lng,advance='no')cooling_h2o
       write(file_cooling,format_lng,advance='no')sum(emicpl)
       write(file_cooling,format_lng,advance='no')sum(emisipl)
       write(file_cooling,format_lng,advance='no')sum(emicat)
       write(file_cooling,format_lng,advance='no')sum(emisiat)
       write(file_cooling,format_lng,advance='no')sum(emioat)
       write(file_cooling,format_lng,advance='no')sum(emispl)
       write(file_cooling,format_lng,advance='no')sum(eminpl)
       write(file_cooling,*)
       ! energetics
       write(file_energetics,format_lng,advance='no')distance
       write(file_energetics,format_lng,advance='no')timen
       write(file_energetics,format_lng,advance='no')tn
       write(file_energetics,format_lng,advance='no')mass_flux
       write(file_energetics,format_lng,advance='no')mass_cons
       write(file_energetics,format_lng,advance='no')momentum_flux_kin
       write(file_energetics,format_lng,advance='no')momentum_flux_the
       write(file_energetics,format_lng,advance='no')momentum_flux_mag
       write(file_energetics,format_lng,advance='no')momentum_flux_vis
       write(file_energetics,format_lng,advance='no')momentum_flux
       write(file_energetics,format_lng,advance='no')momentum_cons
       write(file_energetics,format_lng,advance='no')energy_flux_kin
       write(file_energetics,format_lng,advance='no')energy_flux_the
       write(file_energetics,format_lng,advance='no')energy_flux_int
       write(file_energetics,format_lng,advance='no')energy_flux_mag
       write(file_energetics,format_lng,advance='no')energy_flux_vis
       write(file_energetics,format_lng,advance='no')energy_flux
       write(file_energetics,format_lng,advance='no')energy_gain
       write(file_energetics,format_lng,advance='no')energy_cons
       write(file_energetics,*)
       ! Intensities
       write(file_intensity,format_lng,advance='no')distance
       write(file_intensity,format_lng,advance='no')timen
       write(file_intensity,format_lng,advance='no')tn
       write(file_intensity,format_lng,advance='no')intcpl(1)
       write(file_intensity,format_lng,advance='no')intcpl(5)
       write(file_intensity,format_lng,advance='no')intcpl(6)
       write(file_intensity,format_lng,advance='no')intcpl(2)
       write(file_intensity,format_lng,advance='no')intcpl(3)
       write(file_intensity,format_lng,advance='no')intcpl(4)
       write(file_intensity,format_lng,advance='no')intsipl(1)
       write(file_intensity,format_lng,advance='no')intcat(1)
       write(file_intensity,format_lng,advance='no')intcat(2)
       write(file_intensity,format_lng,advance='no')intcat(4)
       write(file_intensity,format_lng,advance='no')intcat(5)
       write(file_intensity,format_lng,advance='no')intsiat(1)
       write(file_intensity,format_lng,advance='no')intsiat(2)
       write(file_intensity,format_lng,advance='no')intoat(2)
       write(file_intensity,format_lng,advance='no')intoat(1)
       write(file_intensity,format_lng,advance='no')intoat(6)
       write(file_intensity,format_lng,advance='no')intoat(5)
       write(file_intensity,format_lng,advance='no')intspl(7)
       write(file_intensity,format_lng,advance='no')intspl(8)
       write(file_intensity,format_lng,advance='no')intnpl(1)
       write(file_intensity,format_lng,advance='no')intnpl(2)
       write(file_intensity,format_lng,advance='no')intnpl(6)
       write(file_intensity,format_lng,advance='no')intnpl(5)
       write(file_intensity,format_lng,advance='no')intnpl(4)
       write(file_intensity,format_lng,advance='no')intnat(6)
       write(file_intensity,format_lng,advance='no')intnat(7)
       write(file_intensity,*)
       ! Populations
       write(file_populations,format_lng,advance='no')distance
       write(file_populations,format_lng,advance='no')Timen
       write(file_populations,format_lng,advance='no')Tn
       write(file_populations,format_lng,advance='no')pop_cat(1)
       write(file_populations,format_lng,advance='no')pop_cat(2)
       write(file_populations,format_lng,advance='no')pop_cat(3)
       write(file_populations,format_lng,advance='no')pop_cat(4)
       write(file_populations,format_lng,advance='no')pop_cat(5)
       write(file_populations,format_lng,advance='no')pop_nat(1)
       write(file_populations,format_lng,advance='no')pop_nat(2)
       write(file_populations,format_lng,advance='no')pop_nat(3)
       write(file_populations,format_lng,advance='no')pop_nat(4)
       write(file_populations,format_lng,advance='no')pop_nat(5)
       write(file_populations,format_lng,advance='no')pop_oat(1)
       write(file_populations,format_lng,advance='no')pop_oat(2)
       write(file_populations,format_lng,advance='no')pop_oat(3)
       write(file_populations,format_lng,advance='no')pop_oat(4)
       write(file_populations,format_lng,advance='no')pop_oat(5)
       write(file_populations,format_lng,advance='no')pop_sat(1)
       write(file_populations,format_lng,advance='no')pop_sat(2)
       write(file_populations,format_lng,advance='no')pop_sat(3)
       write(file_populations,format_lng,advance='no')pop_sat(4)
       write(file_populations,format_lng,advance='no')pop_sat(5)
       write(file_populations,format_lng,advance='no')pop_siat(1)
       write(file_populations,format_lng,advance='no')pop_siat(2)
       write(file_populations,format_lng,advance='no')pop_siat(3)
       write(file_populations,format_lng,advance='no')pop_siat(4)
       write(file_populations,format_lng,advance='no')pop_siat(5)
       write(file_populations,format_lng,advance='no')pop_cpl(1)
       write(file_populations,format_lng,advance='no')pop_cpl(2)
       write(file_populations,format_lng,advance='no')pop_cpl(3)
       write(file_populations,format_lng,advance='no')pop_cpl(4)
       write(file_populations,format_lng,advance='no')pop_cpl(5)
       write(file_populations,format_lng,advance='no')pop_npl(1)
       write(file_populations,format_lng,advance='no')pop_npl(2)
       write(file_populations,format_lng,advance='no')pop_npl(3)
       write(file_populations,format_lng,advance='no')pop_npl(4)
       write(file_populations,format_lng,advance='no')pop_npl(5)
       write(file_populations,format_lng,advance='no')pop_opl(1)
       write(file_populations,format_lng,advance='no')pop_opl(2)
       write(file_populations,format_lng,advance='no')pop_opl(3)
       write(file_populations,format_lng,advance='no')pop_opl(4)
       write(file_populations,format_lng,advance='no')pop_opl(5)
       write(file_populations,format_lng,advance='no')pop_spl(1)
       write(file_populations,format_lng,advance='no')pop_spl(2)
       write(file_populations,format_lng,advance='no')pop_spl(3)
       write(file_populations,format_lng,advance='no')pop_spl(4)
       write(file_populations,format_lng,advance='no')pop_spl(5)
       write(file_populations,format_lng,advance='no')pop_sipl(1)
       write(file_populations,format_lng,advance='no')pop_sipl(2)
       write(file_populations,format_lng,advance='no')pop_sipl(3)
       write(file_populations,format_lng,advance='no')pop_sipl(4)
       write(file_populations,format_lng,advance='no')pop_sipl(5)
       write(file_populations,*)
       ! Fe+ level populations
       write(file_fe_pops,format_lng,advance='no')distance
       write(file_fe_pops,format_lng,advance='no')timen
       write(file_fe_pops,format_lng,advance='no')tn
       do i=1,nlvfepl
          write(file_fe_pops,format_lng,advance='no')pop_fepl(i)
       end do
       write (file_fe_pops,*)
       ! [FeII] lines
       write(file_fe_lines,format_lng,advance='no')distance
       write(file_fe_lines,format_lng,advance='no')timen
       write(file_fe_lines,format_lng,advance='no')tn
       write(file_fe_lines,format_lng,advance='no')intfepl(22)
       write(file_fe_lines,format_lng,advance='no')intfepl(19)
       write(file_fe_lines,format_lng,advance='no')intfepl(23)
       write(file_fe_lines,format_lng,advance='no')intfepl(30)
       write(file_fe_lines,format_lng,advance='no')intfepl(20)
       write(file_fe_lines,format_lng,advance='no')intfepl(26)
       write(file_fe_lines,format_lng,advance='no')intfepl(27)
       write(file_fe_lines,format_lng,advance='no')intfepl(28)
       write(file_fe_lines,format_lng,advance='no')intfepl(29)
       write(file_fe_lines,format_lng,advance='no')intfepl(31)
       write(file_fe_lines,format_lng,advance='no')intfepl(35)
       write(file_fe_lines,format_lng,advance='no')intfepl(36)
       write(file_fe_lines,format_lng,advance='no')intfepl(37)
       write(file_fe_lines,format_lng,advance='no')intfepl(38)
       write(file_fe_lines,format_lng,advance='no')intfepl(43)
       write(file_fe_lines,format_lng,advance='no')intfepl(44)
       write (file_fe_lines,*)
    endif
    !
    ! Last calculation step: close the files
    if (.not.(again)) then
       close(file_phys)
       if (speci_out_ad) close(file_spec_ad)
       if (speci_out_cd) close(file_spec_cd)
       if (speci_out_fd) close(file_spec_fd)
       if (speci_out_nr) close(file_spec_nr)
       if (do_h2) then
          close(file_h2_lev)
          close(file_h2_lin)
       endif
       if (do_thermal) then
          close(file_cooling)
          close(file_energetics)
          close(file_intensity)
          close(file_populations)
          close(file_fe_pops)
          close(file_fe_lines)
       endif
    end if
  end subroutine write_output


  subroutine write_excit
    use module_tools, only : get_file_number
    use module_phys_var, only : nh2_lev
    use module_h2
    !-----------------------------------------------------------------
    !
    ! Write H2 excitation diagram
    !
    !-----------------------------------------------------------------
    implicit none
    !
    integer(kind=long) :: i
    character(len=15) :: str,str2
    character(len=*), parameter :: format_head='(4x,"v   j   energy(k)   log(n/g)")'
    character(len=*), parameter :: format_i='(i5,1x)'
    character(len=*), parameter :: format_r='(es15.6e3,1x)'
    character(len=*), parameter :: name_file='output/excit.out'
    !
    !
    open(30,file=name_file,status='unknown',&
         access='sequential',form='formatted',action='write')
    write(30,format_head)
    do i=1,nh2_lev
       write(30,format_i,advance='no') h2_lev(i)%v
       write(30,format_i,advance='no') h2_lev(i)%j
       write(30,format_r,advance='no') h2_lev(i)%energy
       write(30,format_r,advance='no') log(h2_lev(i)%col_dens/h2_lev(i)%weight)
       write(30,*)
    end do
    close(30)
  end subroutine write_excit


  subroutine write_species
    use module_tools, only : get_file_number
    use module_phys_var
    use module_chemical_species
    !___________________________________________________________________
    !
    ! Write final abundances (same format as input file)
    !___________________________________________________________________
    !
    implicit none
    character(len=*), parameter :: name_file_in='input/species.in'
    character(len=*), parameter :: name_file_out='output/species.out'
    type(type_specy) :: sp
    integer(kind=long) :: file_in, file_out
    character(len=80) :: charact
    character(len=28) :: beg_l
    character(len=34) :: end_l
    integer :: i,ii
    real :: toto
    !
    ! Open output files
    file_in = get_file_number()
    open(file_in,file=name_file_in,status='old',access='sequential',form='formatted',action='read')
    file_out = get_file_number()
    open(file_out,file=name_file_out,status='replace',access='sequential',form='formatted',action='write')
    !
    ! Copy input species to species.out
    read(file_in,'(a80)') charact
    write(file_out,'(a,es8.1,a,f7.2,a)')   "! Species and abundances at t = ",timen," T=", tn, " K"
    write(file_out,'(a,es8.1,a,f7.2,a)')   "! Important: order = neutrals, mantles ($), cores (#), ions >0, ions <0"
    write(file_out,'(a,1x,12("(",a,")"))') "! Composition:","H","C","N","O","He","+","#","Mg","S","Si","-","Fe"
    write(file_out,'(a)') "! Name, composition, initial density(cm-3), formation enthalpy (kCal/mol)"
    do i=1,ncomment_specfile-1
       read(file_in,'(a80)') charact
    end do
    do i=1,nspec
       read(file_in,format_specy) &
            ii, &                          ! useless
            sp%name, &                     ! name
            sp%useless, &                  ! useless, replaced by CHEMICAL_FORMULA
            sp%density, &                  ! n(X) / nH
            sp%enthalpy                    ! kCal/mol
       write(file_out,format_specy) &
            i, &                           ! useless
            sp%name, &                     ! name
            sp%useless, &                  ! useless, replaced by CHEMICAL_FORMULA
            speci(i)%density/nh, &         ! n(X) / nH
            sp%enthalpy                    ! kCal/mol
!!$       read(file_in,'(a28,e10.3,a34)') beg_l, toto, end_l
!!$       write(file_out,'(a28,1p,d10.3,a34)') beg_l, speci(i)%density/nh, end_l
    end do
  end subroutine write_species


!!$  subroutine write_rates
!!$    use evolution, only : cupn,cupi,cupneg
!!$    !-----------------------------------------------------------------
!!$    !-----------------------------------------------------------------
!!$  end subroutine write_rates

END MODULE MODULE_OUTPUTS

