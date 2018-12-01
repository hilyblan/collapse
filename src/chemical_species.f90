module module_chemical_species
  use module_precision
  use module_io
  !___________________________________________________________________
  !
  ! Variables and routines related to chemical species
  !
  ! * the data type TYPE_SPECY
  ! * the vector containing the species
  ! * the numbers of species in each fluid
  ! * the index useful to find one specy in the vector
  ! * subroutine to write information on one specy
  ! * function to calculate the chemical formula of one specy
  ! * subroutine to calculate the elemental abundances
  !___________________________________________________________________
  !
  implicit none
  !
  ! integer(kind=4), parameter :: name_length = 7  ! (7: old, 8: new) length of specy/element name
  integer(kind=4), parameter :: nelements   = 13 ! number of 'basic' elements: 9 + 1(mg) + 1(gr)
  !
  ! data type for one element
  type type_element
     character(len=name_length) :: name    ! name of the element
     real(kind=dp)              :: mass    ! mass of the element (g)
     real(kind=dp)              :: ab      ! abundance
     real(kind=dp)              :: ab_init ! initial abundance
     real(kind=dp)              :: ab_ref  ! abundance from anders & grevesse 1989
  end type type_element
  !
  ! Elements are 'basic' components : H, C, N ..., initialized in INITIALIZE_ELEMENTS
  type (type_element),dimension(nelements) :: elements
  integer(kind=long) :: ind_elem_h ! index of 'h' in the vector 'elements'
  !
  ! data type for one specy
  type type_specy
     character(len=name_length) :: name                  ! name (ex : 'SiOH+')
     real(kind=dp)              :: density               ! density (cm-3)
     real(kind=dp)              :: netrate               ! formation-destruction (cm-3 s-1)
     real(kind=dp)              :: dens_old              ! density at last call to drive
     real(kind=dp)              :: col_dens              ! column density (cm-2)
     real(kind=dp)              :: mass                  ! mass (g)
     real(kind=dp)              :: enthalpy              ! enthalpy of formation (kcal/mol)
     real(kind=dp)              :: velocity              ! velocity (cm/s)
     real(kind=dp)              :: temperature           ! temperature (k)
     integer(kind=4), dimension(nelements) :: formula    ! chemical formula
     integer(kind=4), dimension(12) :: useless           ! (old system for chemical formula)
     integer(kind=4)            :: index                 ! index of the specy
  end type type_specy
  !
  ! Total number of species in the model and in each fluid
  integer(kind=4) :: nspec       ! number of chemical species
  integer(kind=4) :: nspec_plus  ! nspec + 5 (e-, photon, crp, grain, secpho)
  integer(kind=4) :: nneutrals=0 ! number of neutrals
  integer(kind=4) :: nions=0     ! number of positive ions
  integer(kind=4) :: nneg=0      ! number of negative ions
  integer(kind=4) :: nongrains=0 ! number of species in grain mantle (name with a '*')
  integer(kind=4) :: noncores=0  ! number of species in grain cores (name with a '**')
  real(kind=dp)   :: epsion=0.d0 ! a factor which is 1 if there are ions, 0 otherwise
  real(kind=dp)   :: epsneg=0.d0 ! a factor which is 1 if there are ions, 0 otherwise
  !
  ! vector containing the species, read in READ_SPECIES
  ! species commences at index zero, as tab_specy(none=0)%name is ''
  ! the species are between the index 1 and Nspec
  ! between Nspec+1 and Nspec_plus, we find added species
  ! (electrons, grains,...)
  type(type_specy), save, dimension(:), allocatable :: speci
  !
  ! index of some species
  ! remark : CP -> C+, HP -> H+, SiP -> Si+
  !          SP -> S+, NP -> N+
  integer(kind=4) :: ind_h=0,ind_ph2=0,ind_oh2=0,ind_he=0,ind_o=0,ind_oplus,ind_n=0,ind_c=0
  integer(kind=4) :: ind_s=0,ind_si=0,ind_sio=0,ind_h2o=0,ind_oh=0,ind_co=0,ind_nh3=0
  integer(kind=4) :: ind_cplus=0,ind_hplus=0,ind_siplus=0,ind_splus=0,ind_nplus=0,ind_feplus=0
  integer(kind=4) :: ind_gg0=0,ind_ggplus=0,ind_ggminus=0
  integer(kind=4) :: ind_e=0,ind_photon=0,ind_crp=0,ind_grain=0,ind_secpho=0
  integer(kind=4) :: ind_d=0
  !
  ! mass of some species
  ! remark : CP -> C+, HP -> H+, SiP -> Si+
  !          SP -> S+, NP -> N+
  real(kind=dp) :: mass_h=0._dp,mass_h2=0._dp,mass_ph2=0._dp,mass_oh2=0._dp,mass_he=0._dp,mass_o=0._dp,mass_oplus=0._dp
  real(kind=dp) :: mass_n=0._dp,mass_c=0._dp,mass_s=0._dp,mass_si=0._dp,mass_h2o=0._dp
  real(kind=dp) :: mass_oh=0._dp,mass_co=0._dp,mass_nh3=0._dp,mass_sio=0._dp
  real(kind=dp) :: mass_cplus=0._dp,mass_hplus=0._dp,mass_siplus=0._dp
  real(kind=dp) :: mass_gg0=0._dp,mass_ggplus=0._dp,mass_ggminus=0._dp
  real(kind=dp) :: mass_d=0._dp,mass_splus=0._dp,mass_nplus=0._dp,mass_feplus=0._dp
  !
  ! density of some species
  ! remark : CP -> C+, HP -> H+, SiP -> Si+
  !          SP -> S+, NP -> N+
  real(kind=dp) :: dens_h, dens_h2, dens_ph2, dens_oh2, dens_he, dens_o, dens_oplus
  real(kind=dp) :: dens_n, dens_c, dens_s, dens_si, dens_h2o
  real(kind=dp) :: dens_oh, dens_co, dens_nh3
  real(kind=dp) :: dens_cplus, dens_siplus, dens_hplus, dens_splus, dens_nplus, dens_feplus
  real(kind=dp) :: dens_gg0, dens_ggplus, dens_ggminus
  real(kind=dp) :: dens_ongrains, dens_cor, dens_e
  !
  ! index of beginning and end of each type of species (calculated in INITIALIZE)
  integer(kind=4) :: b_neu=0, e_neu=0 ! neutrals
  integer(kind=4) :: b_ion=0, e_ion=0 ! positive ions
  integer(kind=4) :: b_neg=0, e_neg=0 ! negative ions
  integer(kind=4) :: b_gra=0, e_gra=0 ! species on grain mantles
  integer(kind=4) :: b_cor=0, e_cor=0 ! species on grain mantles
  !
  ! read/write format for one specy
  character(len=*), parameter :: format_specy='(i3,2x,a7,2x,2i2,10i1,es10.3,f10.3)'
  integer(kind=4) :: ncomment_specfile
  !
!  ! variables from "precision.f90" are private
!  private :: dp, long, minus_infinity, plus_infinity
  !
contains
  !
  subroutine initialize_elements
    use module_constants, only : amu, zero
    !_________________________________________________________________________
    !
    ! Initialised the elements and their abundances
    ! from Anders & Grevesse 1989.
    ! Elements are 'basic' components of molecules : H, C, N ...
    !_________________________________________________________________________
    !
    implicit none
    integer(kind=4) :: i
    !
    ! initialization
    elements(:)%name    = ''
    elements(:)%mass    = Zero
    elements(:)%ab      = Zero
    elements(:)%ab_ref  = Zero
    elements(:)%ab_init = Zero
    !
    ! first with mass in amu
    i = 1 ; elements(i)%name = 'H' ; elements(i)%mass = 1.00797_DP ;elements(i)%ab_ref = 1._DP
    ind_elem_H = i ! index of H element in this vector
    i = 2 ; elements(i)%name = 'O' ; elements(i)%mass = 15.9994_DP ;elements(i)%ab_ref = 4.42D-4
    i = 3 ; elements(i)%name = 'C' ; elements(i)%mass = 12.0111_DP ;elements(i)%ab_ref = 3.55D-4
    i = 4 ; elements(i)%name = 'N' ; elements(i)%mass = 14.0067_DP ;elements(i)%ab_ref = 7.94D-5
    i = 5 ; elements(i)%name = 'He'; elements(i)%mass = 4.00260_DP ;elements(i)%ab_ref = 1.00D-1
    i = 6 ; elements(i)%name = 'Na'; elements(i)%mass = 22.9898_DP ;elements(i)%ab_ref = 2.06D-6
    i = 7 ; elements(i)%name = 'Mg'; elements(i)%mass = 24.3_DP    ;elements(i)%ab_ref = 3.70D-5
    i = 8 ; elements(i)%name = 'S' ; elements(i)%mass = 32.0640_DP ;elements(i)%ab_ref = 1.86D-5
    i = 9 ; elements(i)%name = 'Si'; elements(i)%mass = 28.0860_DP ;elements(i)%ab_ref = 3.37D-5
    i = 10; elements(i)%name = 'Fe'; elements(i)%mass = 55.8470_DP ;elements(i)%ab_ref = 3.23D-5
    i = 11; elements(i)%name = 'D' ; elements(i)%mass = 2.00000_DP ;elements(i)%ab_ref = 1.60D-5
    i = 12; elements(i)%name = 'Gr'; elements(i)%mass = 60*elements(3)%mass ;elements(i)%ab_ref = 1.00D-20
    i = 13; elements(i)%name = 'N.'; elements(i)%mass = 15.0_DP    ;elements(i)%ab_ref=elements(4)%ab_ref/440.
    !
    ! conversion of mass : AMU -> grammes
    elements(:)%mass = AMU * elements(:)%mass
    !
  end subroutine initialize_elements
  !
  !
  subroutine read_species
    use module_tools, only : get_file_number
    use module_constants, only : zero
    use module_phys_var, only : nh_init
    use module_parameters_flags
    !_________________________________________________________________________
    !
    ! Read chemical species from 'name_file_speci' in local input/ directory
    ! 
    ! * Compute an index for each species
    ! * Special species are: e-, photon, CRP, grain, SECPHO
    !_________________________________________________________________________
    !
    implicit none
    character(len=name_length) :: charact
    integer(kind=long) :: i,ii,error
    character(len=*), parameter :: name_file_speci='input/species.in'
    integer(4) :: ul
    real(kind=dp) :: density_limit
    !
    ! Determine number of commented lines
    ul = get_file_number()
    open(ul,file=name_file_speci,status='OLD',access='SEQUENTIAL',form='FORMATTED',action='READ')
    error = 0
    ncomment_specfile = 0
    do while (error.eq.0) ! stop at error or end of file
       read(ul,'(a)',iostat=error) charact
       if (charact(1:1).eq."!") then
          ncomment_specfile = ncomment_specfile+1
       else
          exit
       endif
    enddo
    !
    ! opening file
    close(ul)
    open(ul,file=name_file_speci,status='OLD',access='SEQUENTIAL',form='FORMATTED',action='READ')
    !
    ! counts the number of chemical species
    error = 0
    nspec = 0
    do i=1,ncomment_specfile
       read(ul,'(a1)') charact
    end do
    do while (error.eq.0) ! stop at error or end of file
       charact = ''
       read(ul,format_specy,iostat=error) ii,charact ! old
       select case (charact(1:1))
       case ('A':'Z')
          ! if the character is a capital letter, this is a specy
          nspec = nspec + 1
       case (' ')
          ! End of file
          exit
       case default
          ! if it's not a specy, nor end of file -> error
          if (error >=0) stop "E-read_species: error 1 in read_species"
       end select
    enddo
    !
    ! allocation and initialization of the vector species
    nspec_plus = nspec + 5          ! we'll add e-, photon, grains, crp, secpho
    allocate(speci(0:nspec_plus)) ! start a zero : undetermined specy
    speci(:)%name        = ''
    speci(:)%mass        = zero
    speci(:)%enthalpy    = zero
    speci(:)%density     = zero
    speci(:)%netrate     = zero
    speci(:)%dens_old    = zero
    speci(:)%col_dens    = zero
    speci(:)%velocity    = zero
    speci(:)%temperature = zero
    do i=0,nspec_plus
       speci(i)%formula = 0
    enddo
    !
    ! close and re-open file
    close(ul)
    open(ul,file=name_file_speci,status='OLD',access='SEQUENTIAL',form='FORMATTED',action='READ')
    !
    ! Skip the first commented lines
    do i=1,ncomment_specfile
       read(ul,'(a)') charact
    end do
    i = 0
    error = 0
    do while (i < nspec .and. error == 0)
       i = i + 1
       read(ul,format_specy,iostat=error) &
            ii, &                                ! useless
            speci(i)%name, &                     ! name
            speci(i)%useless, &                  ! useless, replaced by CHEMICAL_FORMULA
            speci(i)%density, &                  ! n(X) / nH
            speci(i)%enthalpy                    ! kCal/mol
       !
       ! Convert from abundances to densities [cm-3]
       speci(i)%density = speci(i)%density * nH_init
       !
       IF (error == 0) THEN
          speci(i)%index = i
          ! number of each element in the molecule
          speci(i)%formula = chemical_formula(speci(i)%name)
          ! mass is determined from the formula and the mass of each element (g)
          speci(i)%mass = DOT_PRODUCT(DBLE(elements(:)%mass),DBLE(speci(i)%formula))
          ! indexes and masses of most used chemical species
          SELECT CASE (TRIM(speci(i)%name))
          CASE ('H')
             ind_H  = i
             mass_H = speci(i)%mass
          CASE ('D')
             ind_D  = i
             mass_D = speci(i)%mass
          CASE ('H2')
             ind_ph2  = i
             mass_ph2 = speci(i)%mass
          CASE ('HH')
             ind_oh2  = i
             mass_oh2 = speci(i)%mass
          CASE ('He')
             ind_He  = i
             mass_He = speci(i)%mass
          CASE ('O')
             ind_O  = i
             mass_O = speci(i)%mass
          CASE ('O+')
             ind_Oplus  = i
             mass_Oplus = speci(i)%mass
          CASE ('N')
             ind_N  = i
             mass_N = speci(i)%mass
          CASE ('C')
             ind_C  = i
             mass_C = speci(i)%mass
          CASE ('S')
             ind_S  = i
             mass_S = speci(i)%mass
          CASE ('Si')
             ind_Si  = i
             mass_Si = speci(i)%mass
          CASE ('H2O')
             ind_H2O  = i
             mass_H2O = speci(i)%mass
          CASE ('OH')
             ind_OH  = i
             mass_OH = speci(i)%mass
          CASE ('CO')
             ind_CO  = i
             mass_CO = speci(i)%mass
          CASE ('NH3')
             ind_NH3  = i
             mass_NH3 = speci(i)%mass
          CASE ('SiO')
             ind_SiO  = i
             mass_SiO = speci(i)%mass
          CASE ('C+')
             ind_Cplus  = i
             mass_Cplus = speci(i)%mass
          CASE ('H+')
             ind_Hplus  = i
             mass_Hplus = speci(i)%mass
          CASE ('Si+')
             ind_Siplus  = i
             mass_Siplus = speci(i)%mass
          CASE ('S+')
             ind_Splus  = i
             mass_Splus = speci(i)%mass
          CASE ('N+')
             ind_Nplus  = i
             mass_Nplus = speci(i)%mass
          CASE ('Fe+')
             ind_Feplus  = i
             mass_Feplus = speci(i)%mass
          CASE ('Gr')
             ind_gg0  = i
             mass_gg0 = speci(i)%mass
          CASE ('Gr+')
             ind_ggplus  = i
             mass_ggplus = speci(i)%mass
          CASE ('Gr-')
             ind_ggminus  = i
             mass_ggminus = speci(i)%mass
          end select
       else
          if (error > 0) stop "E-read_species, error 2 in read_species"
       end if
    end do
    close(ul)

    !-------------------------------------------------------
    ! avoid too small numbers : lower limit to the densities
    !-------------------------------------------------------

    !--------------------------------------------------------------------------
    ! addition of e-, photon, CRP, GRAIN, SECPHO with : mass and enthalpy =0.0
    !--------------------------------------------------------------------------
    i = Nspec
    ! e-, initial density is set in INITIALIZE
    i = i + 1
    speci(i)%name  = 'ELECTR'
    speci(i)%index = i
    ind_e = i
    ! GRAIN, initial density is set in INITIALIZE
    i = i + 1
    speci(i)%name  = 'GRAIN'
    speci(i)%index = i
    ind_GRAIN = i
    ! photon, density = constant = 1._DP (cf CHEMISTRY)
    i = i + 1
    speci(i)%name  = 'PHOTON'
    speci(i)%index = i
    speci(i)%density = 1._DP
    ind_PHOTON = i
    ! CRP (cosmique ray proton), density = constant = 1._DP (cf CHEMISTRY)
    i = i + 1
    speci(i)%name  = 'CRP'
    speci(i)%index = i
    speci(i)%density = 1._DP
    ind_CRP = i
    ! SECPHO (secondary photon), initial = constant = 1._DP (cf CHEMISTRY)
    i = i + 1
    speci(i)%name  = 'SECPHO'
    speci(i)%index = i
    speci(i)%density = 1._DP
    ind_SECPHO = i
    !
    if (verbose) then
       write(*,'(1x,a,2i5,a)') "I-read_species: done"
       write(*,'(3x,"total: ",i5)') nspec_plus
       write(*,'(3x,"chemical species: ",i5)') nspec
       write(*,'(3x,"special species (e-, photon, grain, cr, secondary photons): ",i5)') nspec_plus-nspec
    endif
  end subroutine read_species



  subroutine check_species
    use module_debug
    use module_parameters_flags
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    tests if the set of species doesn't have identical species
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !---------------------------------------------------------------------------
    implicit none
    integer(kind=long) :: i, j
    !
    !
    if (debug) then
       do i=1, nspec
          write(*,*) speci(i)%name,speci(i)%mass,speci(i)%density
       enddo
    endif
    do i=1, nspec
       do j=i+1,nspec
          if (speci(i)%name==speci(j)%name) then
             write(*,*) "E-check_species: identical species at index ",i,j
             write(*,*) speci(i)%name,speci(i)%mass,speci(i)%density
             stop
          end if
       end do
    end do
    if (verbose) write(*,*) "I-check_species: done"
  end subroutine check_species

  SUBROUTINE WRITE_SPECY(num,one_specy)
    !--------------------------------------------------------------------------
    ! called by :
    !     WRITE_INFO
    ! purpose :
    !     write informations about one chemical specy
    ! subroutine/function needed :
    ! input variables :
    !     * num -> file number where to write the informations
    !     * one_specy -> (type TYPE_SPECY) chemical specy
    ! output variables :
    ! results :
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(type_specy), INTENT(in) :: one_specy
    INTEGER(KIND=LONG),INTENT(in) :: num
    INTEGER(KIND=LONG) :: i

    WRITE(num,format_specy) &
         one_specy%index, &
         one_specy%name, &
         one_specy%useless, &
         one_specy%density, &
         one_specy%enthalpy
  END SUBROUTINE WRITE_SPECY


  subroutine elemental_abundances
    use module_parameters_flags
    !_________________________________________________________________________
    !
    ! Compute elemental abundances in (H,C,N,...)
    !
    !_________________________________________________________________________
    !
    implicit none
    integer(kind=long) :: i
    !
    ! computation of abundances, using speci(i)%formula and speci(i)%density
    do i=1, nelements
       elements(i)%ab = sum(dble(speci(:)%formula(i)) * speci(:)%density)
    end do
    ! normalization to H abundance
    elements(:)%ab = elements(:)%ab / elements(ind_elem_H)%ab
    elements(:)%ab_init = elements(:)%ab
    !
    if (verbose) write(*,*) "I-abundances: done"
  end subroutine elemental_abundances


  function chemical_formula (name) result (formula)
    !_________________________________________________________________________
    !
    ! Computes for one molecule its composition in elements H, C, ...
    ! The result is a vector containing the number of each element in
    ! the molecule
    !_________________________________________________________________________
    !
    implicit none
    character(len=name_length), intent(in) :: name
    integer(kind=4), dimension(nelements) :: formula
    integer(kind=4) :: i,j,length
    integer(kind=4) :: digit1,digit2,number_element
    integer(kind=4) :: howmany,howmany_element,howmany_digit
    character(len=8), dimension(name_length) :: table_name
    character(len=2) :: charact
    !
    ! initialization
    table_name(:) = ''      ! describes what's in name(i)
    length = len_trim(name) ! length of name, without blanks
    formula = 0
    !
    ! new version with "real" names : He and not HE ...
    ! filling the vector table_name
    do i=1,length
       select case (name(i:i))
       case ('0':'9')! digit
          table_name(i)="digit"
       case ('a':'z')! 'He', 'Si' ...
          table_name(i)="2letters"
       case ('.',':')! Isotopologues
          table_name(i)="2letters"
       case ('+','-','*') ! ion or on grain
          table_name(i)="iongrain"
       case ('A':'Z') ! capital letters
          table_name(i)="element"
       case default ! other
          stop "*** warning, incorrect character in chemical_formula"
       end select
    end do
    !
    ! chemical formula
    i = 1
    do while (i <= length)
       digit1 = 1 ; digit2 = 0
       howmany_digit   = 1
       howmany_element = 1 ! 1 by default (if no digit)
       select case (trim(table_name(i)))
       case ("element")
          ! determine how many characters to take into acount -> howmany
          if (i < length) then
             if (table_name(i+1) == "2letters") howmany_element = 2
          end if
          howmany = howmany_element
          ! extracts from name
          charact = name(i:i+howmany-1)
          ! research of the element
          do j=1,nelements
             if (charact==trim(elements(j)%name)) exit
          end do
          ! if not found ...
          if ((j == nelements) .and. (charact /= trim(elements(j)%name))) then
             stop "*** warning, no chemical element for this molecule"
          end if
          number_element = j

       CASE ("digit")
          ! determine how many characters to take into acount -> howmany
          IF (i < length) THEN
             IF (table_name(i+1) == "digit") howmany_digit = 2
          END IF
          howmany = howmany_digit
          ! research of the number = digit1 or digit1*10+digit2
          digit1 = IACHAR(name(i:i)) - IACHAR('0')
          IF (howmany == 2) digit2 = IACHAR(name(i+1:i+1))-IACHAR('0')
       END SELECT
       !
       ! Determine the number of elements in the molecule
       select case (table_name(i))
       case('element')
          formula(number_element) = formula(number_element) + digit1
          if (howmany_digit == 2) then
            formula(number_element) = formula(number_element) + digit1 * 10 + digit2
          endif
       case('digit')
          formula(number_element) = digit1
          if (howmany_digit == 2) then
            formula(number_element) = digit1 * 10 + digit2
          endif
       case('iongrain') ! do nothing
       end select
       !
       ! look at the next characters in the name
       i = i + howmany
    end do
  end function chemical_formula

  !
  !
  subroutine get_density(array,index,dens)
    !_________________________________________________________________
    !
    ! Read the density for a given chemical species
    !_________________________________________________________________
    !
    use module_phys_var, only : bv_speci
    implicit none
    real(kind=dp), intent(in), dimension(:) :: array
    integer,       intent(in)  :: index
    real(kind=dp), intent(out) :: dens
    integer :: ii
    !
    ii = bv_speci-1+index
    if (index.eq.0) then
       dens = 0.d0
    else
       dens = array(ii)
    endif
  end subroutine get_density


END MODULE MODULE_CHEMICAL_SPECIES
