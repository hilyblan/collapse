module module_chem_react
  use module_precision
  use module_chemical_species
  use module_debug
  use module_io
  use module_constants
  !___________________________________________________________________
  !
  ! Read the chemical network
  ! * Decode reactions
  ! * Compute the mass of each species
  !___________________________________________________________________
  ! 
  implicit none
  !
  ! data type of one reaction
  type type_reaction
     integer(kind=4),dimension(2)    :: r                  ! index of 2 reactants
     integer(kind=4),dimension(4)    :: p                  ! index of 4 products
     real(kind=dp)                   :: gamma, alpha, beta ! arrhenius coefficients
     real(kind=dp)                   :: de                 ! exo(if >0)/endo(if <0)-thermicity
     real(kind=dp)                   :: mass_prod          ! sum of the product's masse (g)
     integer(kind=4)                 :: nprod_m1           ! number of products minus one
     character(len=32)               :: comment            ! field for comment
     character(len=12)               :: type               ! type of the reaction
     integer                         :: itype              ! integer-coded type of the reaction
  end type type_reaction
  !
  ! number of reactions in the model and in each type of reaction
  ! calculated in READ_REACTIONS
  ! Nreact is modified in ADD_INVERSE_REACTIONS
  ! as we add the missing endothermic reactions (if the boolean
  ! do_we_add_reactions is .TRUE.)
  integer(kind=4)           :: nreact          ! number of reactions
  logical, parameter        :: do_we_add_reactions=.true.   ! add the reverse reactions ?
  !
  ! number of reactions in each type of reaction (PHOTO, SPUTTER, ...)
  ! calculated in INITIALIZE, exept Nrever_new et Nrever : in ADD_INVERSE_REACTIONS
  integer(kind=4) :: nphoto=0, ncr_io=0, ncr_de=0, nh2_fo=0, nthree=0
  integer(kind=4) :: nsputt=0, nerosi=0, nadsor=0, nother=0, nrever=0
  integer(kind=4) :: ndisso=0

  ! index of beginning and end of each type of reaction
  ! calculated in INITIALIZE, exept b_rever et e_rever : in ADD_INVERSE_REACTIONS
  integer(kind=4) :: b_photo=0, e_photo=0
  integer(kind=4) :: b_cr_io=0, e_cr_io=0
  integer(kind=4) :: b_cr_de=0, e_cr_de=0
  integer(kind=4) :: b_h2_fo=0, e_h2_fo=0
  integer(kind=4) :: b_three=0, e_three=0
  integer(kind=4) :: b_sputt=0, e_sputt=0
  integer(kind=4) :: b_erosi=0, e_erosi=0
  integer(kind=4) :: b_adsor=0, e_adsor=0
  integer(kind=4) :: b_disso=0, e_disso=0
  integer(kind=4) :: b_other=0, e_other=0
  integer(kind=4) :: b_rever=0, e_rever=0
  !
  ! vector containing the chemical reactions
  ! stored betwenn the index 1 and Nreact
  type(type_reaction), dimension(nreact_max) :: react
  !
  ! useful for determining wich reactions are correct -> READ_REACTIONS
  character(len=name_length), parameter   :: charact_blank=''
  integer(kind=4), parameter :: rp_undefined = -1 ! if reactant or product undefined -> reaction is incorrect
  integer(kind=4), parameter :: p_none = 0        ! si product is not present in chemical species
  !
  ! read/write format for chemical reactions
  ! character(len=*), private, parameter  :: format_reaction      = '(2(a7,1x),"= ",4(a7,1x),4es12.3,i4,i6,"!",2x,a)'
  ! character(len=*), private, parameter  :: format_network_read  = '(2(a7,1x),2x  ,4(a7,1x),3es12.3,a)'
  ! character(len=*), private, parameter  :: format_network_write = '(2(a7,1x),"= ",4(a7,1x),3es12.3,i4,i6,"!",2x,a)'
  ! character(len=*), private, parameter  :: format_branching     = '(2(a7,1x),"= ",4(a7,1x),6(2i2,2x),3es12.3,i4,i6,"!",2x,a)'
  !
  !
contains
  !
  subroutine write_reaction(iu,num,reaction)
    !_________________________________________________________________
    !
    ! Write informations about one chemical reaction
    !_________________________________________________________________
    !
    implicit none
    integer(kind=long),  intent(in) :: iu
    integer,             intent(in) :: num
    type(type_reaction), intent(in) :: reaction
    !
    write(iu,format_reaction) &
         speci(reaction%r(1:2))%name, &
         speci(reaction%p(1:4))%name, &
         reaction%gamma,&
         reaction%alpha,&
         reaction%beta, &
         reaction%de,&
         reaction%itype,&
         num,&
         reaction%comment
  end subroutine write_reaction
  !
  !
  function reactions_identical(reaction1, reaction2) result(res)
    !___________________________________________________________________________
    !
    ! Test if 2 reactions are identical (same reactants and products)
    ! whatever the order in wich reactants and products are written
    !___________________________________________________________________________
    !
    implicit none
    type(type_reaction), intent(in) :: reaction1 , reaction2
    logical :: reactants_identical, produits_identical, res
    integer(kind=long) :: k,l
    logical,dimension(4) :: prod_id1, prod_id2
    !
    ! test if reactants are identical
    reactants_identical = &
         ((reaction1%r(1) == reaction2%r(1) .and. &
         reaction1%r(2) == reaction2%r(2)) .or. &
         (reaction1%r(1) == reaction2%r(2) .and. &
         reaction1%r(2) == reaction2%r(1)))

    ! if reactants identical, test if products are identical
    produits_identical = .false.
    if (reactants_identical) then
       ! test if each product of reaction 1 is in reaction 2
       prod_id1(:) = .false.
       do k=1,4
          do l=1,4
             if (reaction1%p(k) == reaction2%p(l)) prod_id1(k) = .true.
          end do
       enddo
       ! test if each product of reaction 2 is in reaction 1
       prod_id2(:) = .false.
       do k=1,4
          do l=1,4
             if (reaction2%p(k)==reaction1%p(l)) prod_id2(k)=.true.
          end do
       enddo
       ! produits_identical is .true. if all products are common to the 2 reactions
       produits_identical = all(prod_id1) .and. all(prod_id2)
    end if
    ! result is : (same reactants) and (sames products)
    res = reactants_identical .and. produits_identical
  end function reactions_identical
  !
  !
  !
  subroutine read_reactions
    use module_tools, only : get_file_number
    use module_constants, only : zero
    use module_parameters_flags
    !___________________________________________________________________________
    !
    ! Read chemical reactions in 'name_file_chemistry' in local input/ directory
    ! 
    ! Index each reactant and product
    !___________________________________________________________________________
    !
    implicit none
    character(len=5), parameter :: comment='!'
    character(len=5), parameter :: e_comment='!end'
    character(len=5), parameter :: e_file='END'
    character(len=32) :: ref=''
    character(len=name_length) :: r1, r2, p1, p2, p3, p4 ! names for reactants and products
    real(kind=dp) :: alpha, beta, gamma        ! reaction coefficients
    integer(kind=4) :: i, ir,error,start_comment
    character(len=*), parameter :: name_file_chemistry='input/chemistry.in'
    integer                     :: file_chemistry,chemical_network
    logical                     :: is_comment
    !
    ! initialization
    react(:)%type = ''
    react(:)%comment  = ''
    do i=1,2
       react(:)%r(i) = rp_undefined
    end do
    do i=1,4
       react(:)%p(i) = rp_undefined
    end do
    react(:)%gamma = zero
    react(:)%alpha = zero
    react(:)%beta = zero
    react(:)%de = zero
    react(:)%nprod_m1 = -1 ! don't forget the 'minus one' !
    react(:)%mass_prod = zero
    !
    ! Open the chemical network output file
    file_chemistry = get_file_number()
    open(file_chemistry,file=name_file_chemistry,status='old',&
         access='sequential',form='formatted',action='read')
    !
    nreact = 0
    r1  = ''
    ref = ''
    error = 0
    is_comment = .false.
    do while (r1.ne.e_file .and. nreact.lt.nreact_max .and. error.eq.0)
       r1 = ''; r2 = ''; p1 = ''; p2 = ''; p3 = ''; p4 = ''
       read(file_chemistry,format_network_read,iostat=error) r1,r2,p1,p2,p3,p4,gamma,alpha,beta,ref
       is_comment = (index(adjustl(r1),"!").eq.1)
       if (is_comment) then
          error = 0
          cycle
       endif
       ! if (debug) write(*,format_network_write) r1,r2,p1,p2,p3,p4,gamma,alpha,beta,-1,nreact,ref
       if (r1.eq.e_file.or.len_trim(r1)+len_trim(r2).eq.0) exit
       nreact = nreact + 1
       if (r1.ne.e_file.and.error.eq.0) then
          ! search for indexes of reactants and products
          do i=0, nspec_plus ! and not 1,nspec, as we search also none, e-, photon, crp, grain
             if (r1 == speci(i)%name) react(nreact)%r(1) = i
             if (r2 == speci(i)%name) react(nreact)%r(2) = i
             if (p1 == speci(i)%name) react(nreact)%p(1) = i
             if (p2 == speci(i)%name) react(nreact)%p(2) = i
             if (p3 == speci(i)%name) react(nreact)%p(3) = i
             if (p4 == speci(i)%name) react(nreact)%p(4) = i
          end do
          !
          ! Check for unidentified species
          if (any(react(Nreact)%R(:).eq.RP_undefined).or.any(react(Nreact)%P(:).eq.rp_undefined)) then
             write(*,'(a,i5,6a9,4es10.2)') "Error: unidentified species in reaction #",Nreact,r1,r2,p1,p2,p3,p4,gamma,alpha,beta
             write(*,'(a,i5,2i4)') "Error: unidentified species in reaction #",Nreact,react(Nreact)%R(:)
             write(*,'(a,i5,4i4)') "Error: unidentified species in reaction #",Nreact,react(Nreact)%P(:)
             stop
          endif
          !
          ! fill the fields of the reaction
          react(Nreact)%gamma = gamma
          react(Nreact)%alpha = alpha
          react(Nreact)%beta = beta
          !
          ! Read comment field
          start_comment = index(ref,"!")+1
          react(nreact)%comment = trim(adjustl(ref(start_comment:)))
          !
          !react(Nreact)%DE = DE      ! DE is not read, but calculated after
          ! calculates the number of products and the sum of their mass (used in CHEMISTRY)
          if (p1 /= charact_blank) then
             react(nreact)%nprod_m1 = react(nreact)%nprod_m1 + 1
             react(nreact)%mass_prod = react(nreact)%mass_prod + &
                  speci(react(nreact)%p(1))%mass
          endif
          if (p2 /= charact_blank) then
             react(nreact)%nprod_m1 = react(nreact)%nprod_m1 + 1
             react(nreact)%mass_prod = react(nreact)%mass_prod + &
                  speci(react(nreact)%p(2))%mass
          endif
          if (p3 /= charact_blank) then
             react(nreact)%nprod_m1 = react(nreact)%nprod_m1 + 1
             react(nreact)%mass_prod = react(nreact)%mass_prod + &
                  speci(react(nreact)%p(3))%mass
          endif
          if (p4 /= charact_blank) then
             react(nreact)%nprod_m1 = react(nreact)%nprod_m1 + 1
             react(nreact)%mass_prod = react(nreact)%mass_prod + &
                  speci(react(nreact)%p(4))%mass
          endif
       endif
    end do
    !
    if (verbose) then
       write(*,'(a,i5)') " I-read_reactions"
       write(*,'(3x,a,i5)') "#(reactions): ",nreact
    endif
    close(file_chemistry)
    !
  end subroutine read_reactions



  subroutine energy_defect_reaction(i)
    use module_constants, only : kcalev, zero
    !---------------------------------------------------------------------------
    ! called by :
    !    REACTION_TYPE
    ! purpose :
    !    Calculate the energy defect (DE) for the reaction i.
    !    If this reaction has one specy with and unkown enthalpie, set DE = Zero
    !    exept if this specy appears (with the same occurence) in BOTH the
    !    reactants and products. In this case, the uncertainty cancels and one
    !    can calculate DE = enthalpy(reactants)-enthalpy(products).
    !    REMARK :
    !       energy defect is set to zero for endothermic reactions
    !       (DE < Zero) => DE = Zero
    ! subroutine/function needed :
    ! input variables :
    !    i -> index of the chemical reaction
    ! ouput variables :
    ! results :
    !     react
    !---------------------------------------------------------------------------
    implicit none
    integer(kind=long), intent(in) :: i
    real(kind=dp), parameter       :: enthalpy_threshold=-99.9_dp
    logical                        :: enthalpy_unknown
    integer(kind=long), dimension(0:nspec_plus) :: occ_reactants
    integer(kind=long), dimension(0:nspec_plus) :: occ_products
    integer :: j
    !
    ! --- check if enthalpy of one specy is unknown         ---
    ! --- and if this specy appears with the same occurence ---
    ! --- in reactants and in products                      ---
    ! initialization
    occ_reactants = 0 ! counts how many time appears one specy in the reactants
    occ_products = 0  ! counts how many time appears one specy in the products
    do j=1,2
       occ_reactants(react(i)%r(j)) = occ_reactants(react(i)%r(j))+1
    end do
    do j=1,4
       occ_products(react(i)%p(j)) = occ_products(react(i)%p(j))+1
    end do
    enthalpy_unknown = .false.
    do j=1,2          ! check for enthalpy of reactants
       if (speci(react(i)%r(j))%enthalpy < enthalpy_threshold .and. &
            occ_reactants(react(i)%r(j)) /= occ_products(react(i)%r(j))) &
            enthalpy_unknown = .true.
    end do
    do j=1,4          ! check for enthalpy of reactants
       if (speci(react(i)%p(j))%enthalpy < enthalpy_threshold .and. &
            occ_reactants(react(i)%p(j)) /= occ_products(react(i)%p(j))) &
            enthalpy_unknown = .true.
    end do

    ! --- calculation of DE (eV, whereas speci%enthaply is in kCal/mol) ---
    IF (enthalpy_unknown) THEN
       ! DE = 0.0 in this case
       react(i)%DE = Zero
    ELSE
       ! DE = enthalpy (reactants) - enthalpy (products)
       react(i)%DE = kCaleV * ( &
            SUM(DBLE(speci(react(i)%R(:))%enthalpy)) - &
            SUM(DBLE(speci(react(i)%P(:))%enthalpy)))
    END IF

    ! set DE=Zero for endothermic reactions
    IF (react(i)%DE < Zero) react(i)%DE = Zero

  END SUBROUTINE ENERGY_DEFECT_REACTION


  subroutine reaction_type
    use module_constants, only : zero
    use module_phys_var, only:nh2_lev
    !---------------------------------------------------------------------------
    !
    !    Find the type of each reaction and re-order the entire set
    !    according to the different types. Compute the
    !    indexes (beginning, end, number) of each reaction type.
    !    Also calculate DE for each reaction, according to it's type.
    !    Order is :
    !       (1) : photo-reactions (type='PHOTO')
    !       (2) : direct cosmic ray ionization or dissociation (type='CRIO_DIR')
    !       (3) : indirect cosmic ray ionization or dissociation (type='CRIO_IND')
    !       (4) : cosmic ray induced desorption from grains (type='CR_DE')
    !       (5) : H2, HD and D2 formation (type='H2_FO')
    !       (6) : three body reactions on grain surfaces (type='THREE')
    !       (7) : sputtering of grain mantles (type='SPUTT')
    !       (8) : erosion of grain cores (type='EROSI')
    !       (9) : adsorption on grain surface (type='ADSOR')
    !      (10) : Collisional dissociation of H2 (type='DISSO') - see below
    !      (11) : all other reactions  (type='OTHER')
    !      (12) : Grain neutralization (type='NEUTGR')
    !      (13) : Grain charging       (type='CHARGR')
    ! (obsolete): reverse endothermic reactions (type='REVER')
    !             these reactions are defined and added in ADD_REVERSE_REACTIONS
    !
    !  22 juin 2001 - JLB - Add collisional dissociation of H2
    !       WARNING ! Required order in chemistry file : H2 + X -> X + H + H
    !                 identification is done on R1=H2, P2=H, P3=H, R2 = P1
    !---------------------------------------------------------------------------
    implicit none
    ! Local variables
    integer(kind=long) :: i,ir=0
    type (type_reaction), dimension(:), allocatable :: react_aux
    logical, dimension(:), allocatable :: done
    character(len=*), parameter :: rname='read_type'
    !
    ! Save the unsorted reaction set
    ! react_aux = old set (unsorted)
    allocate(react_aux(1:nreact), done(nreact))
    react_aux(1:nreact) = react(1:nreact)
    done = .false.
    !
    ! Initialize the vector 'react'
    react(:)%itype = -1000
    react(:)%type = ''
    react(:)%comment = ''
    do i=1,2
       react(:)%r(i) = rp_undefined
    end do
    do i=1,4
       react(:)%p(i) = rp_undefined
    end do
    react(:)%gamma = zero
    react(:)%alpha = zero
    react(:)%beta = zero
    react(:)%de = zero
    react(:)%nprod_m1 = -1 ! useless here ...
    react(:)%mass_prod = zero
    !
    ! Find the type of each reaction.
    ! Re-order the entire set and calculate the energy defect
    ! according to this type
    ! NOW : 'react' = sorted reactions set
    !
    !--- photo-reactions (type='PHOTO')---
    b_photo = 1 ! index of beginning
    do i=1, nreact
       if (any(react_aux(i)%r(1:2).eq.ind_photon)) then
          react_aux(i)%type = 'PHOTO'
          nphoto = nphoto + 1
          ir = ir+1
          react(ir) = react_aux(i)
          react(ir)%de = zero
          done(i) = .true.
       endif
    end do
    e_photo = b_photo+nphoto-1 ! index of end
    !
    ! Direct cosmic-ray ionization or dissociation (type='cr_io')
    b_cr_io = e_photo + 1 ! index of beginning
    do i=1,nreact
       if (done(i)) cycle
       if (any(react_aux(i)%p(:).eq.ind_grain)) cycle
       if (any(react_aux(i)%r(1:2).eq.ind_crp)) then
          react_aux(i)%type = 'CRIO_DIR'
          ncr_io = ncr_io + 1
          ir = ir+1
          react(ir) = react_aux(i)
          react(ir)%de = zero
          done(i) = .true.
       endif
    enddo
    !
    ! Indirect cosmic-ray ionization or dissociation (type='cr_io')
    do i=1,nreact
       if (done(i)) cycle
       ! if (any(react_aux(i)%p(:).eq.ind_grain)) cycle
       if (any(react_aux(i)%r(1:2).eq.ind_secpho).and.react_aux(i)%p(2).ne.ind_grain) then
          react_aux(i)%type = 'CRIO_IND'
          done(i) = .true.
       else
          done(i) = .false.
       endif
       if (done(i)) then
          ncr_io = ncr_io + 1
!!$          ir     = ncr_io+b_cr_io-1
          ir = ir+1
          react(ir) = react_aux(i)
          react(ir)%de = zero
       endif
    enddo
    e_cr_io = b_cr_io + Ncr_io - 1 ! index of end
    !
    ! Cosmic-ray induced desorption from grains (type='CR_DE') ---
    b_cr_de = e_cr_io + 1 ! index of beginning
    do i=1, nreact
       if (done(i)) cycle
       if (any(react_aux(i)%r(1:2).eq.ind_crp).and.react_aux(i)%p(2).eq.ind_grain) then
          react_aux(i)%type = 'CR_DE'
          ncr_de = ncr_de + 1
!!$          ir = ncr_de+ b_cr_de-1
          ir = ir+1
          react(ir) = react_aux(i)
          react(ir)%de = zero
          done(i) = .true.
       endif
    end do
    e_cr_de = b_cr_de + Ncr_de - 1 ! index of end
    !
    !--- H2, HD and D2 formation (type='H2_FO') ---
    b_h2_fo = e_cr_de + 1 ! index of beginning
    do i=1,nreact
       if (done(i)) cycle
       if (&
            all(react_aux(i)%r(1:2).eq.ind_h).or.&
            all(react_aux(i)%r(1:2).eq.ind_d).or.&
            (react_aux(i)%r(1).eq.ind_h.and.react_aux(i)%r(2).eq.ind_d)) then
          react_aux(i)%type = 'H2_FO'
          nh2_fo = nh2_fo + 1
          ir = ir+1
          react(ir) = react_aux(i)
          call energy_defect_reaction(ir) ! explicitely calculated
          done(i) = .true.
       endif
    enddo
    e_h2_fo = b_h2_fo + nh2_fo - 1 ! index of end
    !
    ! Three body reactions on grain surfaces (type='THREE') ---
    b_three = e_h2_fo + 1 ! index of beginning
    do i=1, nreact
       if (done(i)) cycle
       if (react_aux(i)%r(1).eq.ind_h.and.react_aux(i)%p(2).eq.ind_grain) then
          react_aux(i)%type = 'THREE'
          nthree = nthree + 1
          ir = ir+1
          react(ir) = react_aux(i)
          react(ir)%de = zero
          react(ir)%type = 'THREE'
!!$          react(nthree+ b_three-1) = react_aux(i)
!!$          react(nthree+ b_three-1)%de = zero
!!$          react(Nthree+ b_three-1)%type = 'THREE'
          done(i) = .true.
       endif
    end do
    e_three = b_three + Nthree - 1 ! index of end
    !
    ! Sputtering of grain mantles (type='SPUTT') ---
    b_sputt = e_three + 1 ! index of beginning
    do i=1, nreact
       if (trim(react_aux(i)%type) == '' .and. &
            react_aux(i)%p(3) == ind_grain .or. &
            react_aux(i)%p(4) == ind_grain) then
          nsputt = nsputt + 1
          react(nsputt+ b_sputt-1) = react_aux(i)
          react(nsputt+ b_sputt-1)%de = zero
          react_aux(i)%type = 'SPUTT'
          react(Nsputt+ b_sputt-1)%type = 'SPUTT'
          done(i) = .true.
       endif
    end do
    e_sputt = b_sputt + nsputt - 1 ! index of end
    !
    ! Erosion of grain cores (type='EROSI') ---
    b_erosi =  e_sputt + 1 ! index of beginning
    do i=1, nreact
       if (done(i)) cycle
       if (react_aux(i)%p(1).eq.ind_grain) then
          react_aux(i)%type = 'EROSI'
          nerosi = nerosi + 1
          ir = ir+1
          react(ir) = react_aux(i)
          react(ir)%de = zero
          react(ir)%type = 'EROSI'
          done(i) = .true.
       endif
    end do
    e_erosi = b_erosi + nerosi - 1 ! index of end
    !
    ! Adsorption on grain surface ---
    b_adsor = e_erosi + 1 ! index of beginning
    do i=1, nreact
       if (done(i)) cycle
       if (react_aux(i)%r(2) == ind_grain) then
          react_aux(i)%type = 'ADSOR'
          nadsor = nadsor + 1
          ir = ir+1
          react(ir) = react_aux(i)
          react(ir)%de = zero
          react(ir)%type = 'ADSOR'
          done(i) = .true.
       endif
    end do
    e_adsor = b_adsor + nadsor - 1 ! index of end
    !
    ! Collisional dissociation of H2
    b_disso = e_adsor + 1 ! index of beginning
    do i=1, nreact
       if (done(i)) cycle
       if (&
            react_aux(i)%r(1) == ind_ph2 .and. &
            react_aux(i)%p(2) == ind_h .and. &
            react_aux(i)%p(3) == ind_h .and. &
            react_aux(i)%p(1) == react_aux(i)%r(2) ) then
          ndisso = ndisso + 1
          react_aux(i)%type = 'DISSO'
          ir = ir+1
          react(ir) = react_aux(i)
          react(ir)%de = zero
          react(ir)%type = 'DISSO'
          if (nh2_lev.eq.0) then
             write(*,*) rname,"Collisional dissociation of H2 requires H2 levels."
             stop
          endif
          done = .true.
       endif
    end do
    e_disso = b_disso + ndisso - 1 ! index of end
    !
    ! Other-reactions ---
    b_other = e_disso + 1 ! index of beginning
    do i=1, nreact
       if (done(i)) cycle
       nother = nother + 1
       ir     = ir+1 !nother+b_other-1
       if (any(react_aux(i)%p(1:4).eq.ind_gg0))  then
          ! Grain neutralization
          react_aux(i)%type = 'NEUTGR'
       else if (&
            any(react_aux(i)%p(1:4).eq.ind_ggminus).or.&
            any(react_aux(i)%p(1:4).eq.ind_ggplus)) then
          ! Grain charging
          react_aux(i)%type = 'CHARGR'
       else
          react_aux(i)%type = 'OTHER'
       endif
       react(ir) = react_aux(i)
       call energy_defect_reaction(ir) ! explicitly calculated
    end do
    e_other = b_other + Nother - 1 ! index of end
    !
    ! integer-coded reaction types
    where(react%type=='PHOTO')    react%itype = 1
    where(react%type=='CRIO_DIR') react%itype = 2
    where(react%type=='CRIO_IND') react%itype = 3
    where(react%type=='CR_DE')    react%itype = 4
    where(react%type=='H2_FO')    react%itype = 5
    where(react%type=='THREE')    react%itype = 6
    where(react%type=='SPUTT')    react%itype = 7
    where(react%type=='EROSI')    react%itype = 8
    where(react%type=='ADSOR')    react%itype = 9
    where(react%type=='DISSO')    react%itype = 10
    where(react%type=='OTHER')    react%itype = 11
    where(react%type=='NEUTGR')   react%itype = 12
    where(react%type=='CHARGR')   react%itype = 13

    !--- Nrever, b_rever, e_rever -> see ADD_REVERSE_REACTIONS ---

    !--- deallocation of the temporary variable ---
    DEALLOCATE(react_aux,done)

  END SUBROUTINE REACTION_TYPE


  subroutine check_reactions
    use module_parameters_output, only : screen
    !---------------------------------------------------------------------------
    ! Tests if the set of reactions is correct
    !        * each reaction must have at least 2 reactants and 1 product
    !        * reactants and products must be in the set of chemical species
    !        * species have to be conserved (exept for reaction of the type
    !           THREE, EROSI or ADSOR)
    !        * charge has to be conserved
    !
    !---------------------------------------------------------------------------
    implicit none
    integer(kind=long)::i,j,ind
    integer(kind=long),dimension(2,nelements):: formula_reactants
    integer(kind=long),dimension(4,nelements):: formula_products
    integer(kind=long) :: nions_reactants, nneg_reactants, ne_reactants, charge_reactants
    integer(kind=long) :: nions_products, nneg_products, ne_products, charge_products
    character(len=name_length) :: name
    logical :: incorrect
    !
    do i=1,nreact
       ! --- checks if reaction is correct ---
       incorrect = &
            ! at least 2 reactants
       react(i)%R(1) == P_none .OR. &
            react(i)%R(2) == P_none .OR. &
            ! at least 1 product
       react(i)%Nprod_m1 < 0 .OR. &
            ! use only species in the current set
       react(i)%R(1) == RP_undefined .OR. &
            react(i)%R(2) == RP_undefined .OR. &
            react(i)%P(1) == RP_undefined .OR. &
            react(i)%P(2) == RP_undefined .OR. &
            react(i)%P(3) == RP_undefined .OR. &
            react(i)%P(4) == RP_undefined
       if (incorrect) then
          write(*,*) "E-check_reactions, incorrect reaction ",i
          call write_reaction(screen,i,react(i))
          stop
       end if

       ! --- species conservation ---
       formula_reactants(:,:) = 0
       formula_products(:,:)  = 0
       select case (react(i)%type)
       case ('THREE', 'EROSI', 'ADSOR')
          ! for these reaction, one can't chek the conservation of species
       case default
          do j=1,2 ! chemical formula of reactants
             name=speci(react(i)%r(j))%name
             select case (name)
             case ('PHOTON', 'CRP', 'ELECTR', 'GRAIN')
                ! do nothing
             case default
                formula_reactants(j,:) = speci(react(i)%r(j))%formula
             end select
          end do
          do j=1,4 ! chemical formula of products
             name=speci(react(i)%p(j))%name
             select case (name)
             case ('PHOTON', 'CRP', 'ELECTR', 'GRAIN')
                ! do nothing
             case default
                formula_products(j,:) = speci(react(i)%p(j))%formula
             end select
          end do
          ! checks the conservation of species between reactants and products
          do j=1,nelements
             if (sum(formula_reactants(:,j)) /= sum(formula_products(:,j))) then
                write(*,*) "E-check_reactions, element ",elements(j)%name," is not conserved in the reaction ",i
                call write_reaction(screen,i,react(i))
                stop
             end if
          enddo
       end select
       !
       ! --- conservation of the charge ---
       Nions_reactants=0 ; Nneg_reactants=0 ; Ne_reactants=0 ; charge_reactants=0
       Nions_products=0 ; Nneg_products=0 ; Ne_products=0 ; charge_products=0
       ! counts the number of ions and electrons in the reactants
       do j=1,2
          ind=react(i)%r(j)
          if (ind == ind_e) then ! electron
             ne_reactants = ne_reactants + 1
          else
             if (ind >= b_ion .and. ind <= e_ion) then ! ion
                nions_reactants = nions_reactants + 1
             else
                if (ind >= b_neg .and. ind <= e_neg) then ! negatif ion
                   nneg_reactants = nneg_reactants + 1
                end if
             end if
          end if
       end do
       ! counts the number of ions and electrons in the products
       do j=1,4
          ind = react(i)%p(j)
          if (ind == ind_e) then ! electron
             ne_products = ne_products + 1
          else
             if (ind >= b_ion .and. ind <= e_ion) then !ion
                nions_products = nions_products + 1
             else
                if (ind >= b_neg .and. ind <= e_neg) then ! negative ion
                   nneg_products = nneg_products + 1
                end if
             end if
          end if
       end do
       ! charge must be identical in reactants and in products
       charge_reactants = Nions_reactants-Ne_reactants-Nneg_reactants
       charge_products  = Nions_products-Ne_products-Nneg_products
       if (charge_reactants.ne.charge_products) then
          write(*,*) "E-check_reactions, charge isn't conserved in the reaction ",i,charge_reactants,charge_products
          write(*,*) nions_reactants,ne_reactants,nneg_reactants
          write(*,*) nions_products,ne_products,nneg_products
          call write_reaction(screen,i,react(i))
          stop
       endif

       ! --- look for the same reaction in the rest of the chemical set ---
       do j=i+1,nreact
          if (reactions_identical(react(i),react(j))) then
             write(*,*) "E-check_reactions, reactions ",i," and ",j," are identical"
             call write_reaction(screen,i,react(i))
             call write_reaction(screen,i,react(j))
             stop
          end if
       end do
    end do
  end subroutine check_reactions



  subroutine add_reverse_reactions
    use module_constants, only : zero
    use module_parameters_flags
    !---------------------------------------------------------------------------
    !
    ! Add reverse reactions when :
    !        * ion-neutral reaction (no barrier)
    !        * 2 reactants <-> 2 products
    !        * 0< DE <= DE_threshold
    !        * no radiative asspciation
    !        * no recombination
    !        * reverse reaction isn't already included
    ! The reverse reaction has the following :
    !        * same gamma, alpha as direct reaction
    !        * beta=DE(direct)*11600.
    !        * DE=0.0
    !
    !---------------------------------------------------------------------------
    implicit none
    real(kind=dp), parameter :: de_threshold=0._dp ! ev
    integer(kind=long) :: i,j
    logical :: eliminate
    type(type_reaction) :: reaction_aux
    !
    ! the reactions are added at the end of the set
    b_rever=Nreact+1

    ! looks for ion-neutral reaction only -> type 'OTHER'
    DO i=b_other,e_other
       ! eliminates the following reactions
       eliminate= &
            (react(i)%Nprod_m1 /= 1) .OR. & ! more or less than 2 products
            (react(i)%DE <= Zero) .OR. &       ! DE <= 0
            (react(i)%DE > DE_threshold) .OR. &    ! DE > DE_threshold
            ((react(i)%R(1)>=b_neu) .AND. (react(i)%R(1)<=e_neu).AND.&
            (react(i)%R(2)>=b_neu) .AND. (react(i)%R(2)<=e_neu)).OR. & ! neutral-neutral (possibility of a barrier)
            (react(i)%R(1)==ind_e .OR. react(i)%R(2)==ind_e) .OR. & ! recombination
            (react(i)%P(1)==ind_PHOTON .OR. react(i)%P(2)==ind_PHOTON) ! radiative association

       ! tests if the reverse reaction isn't already in the set
       j=1
       DO WHILE ((j<=e_other) .AND. (.NOT. eliminate))
          ! build the reverse reaction
          reaction_aux%R(1)=react(i)%P(1)
          reaction_aux%R(2)=react(i)%P(2)
          reaction_aux%P(1)=react(i)%R(1)
          reaction_aux%P(2)=react(i)%R(2)
          reaction_aux%P(3:4)=P_none
          ! test if it doesn't exist already
          eliminate=REACTIONS_IDENTICAL(react(j),reaction_aux)
          j=j+1
       END DO

       ! if the reaction is not eliminated, we add the corresponding reverse reaction
       IF (.NOT. eliminate) THEN
          ! update of Nreact
          Nreact=Nreact+1
          ! Nreact_max must be large enough
          IF (Nreact > Nreact_max) STOP "*** WARNING, Nreact_max is too small"
          ! update of Nrever (e_rever is calculated at the end of the subroutine)
          Nrever=Nrever+1

          ! definition of the reverse reaction
          react(Nreact)=reaction_aux ! reactants and products
          react(Nreact)%gamma=react(i)%gamma
          react(Nreact)%alpha=react(i)%alpha
          react(Nreact)%beta=react(i)%DE*11600._DP
          react(Nreact)%DE=Zero
          react(Nreact)%Nprod_m1=1 ! only 2 products
          react(Nreact)%mass_prod= &
               speci(react(Nreact)%P(1))%mass + &
               speci(react(Nreact)%P(2))%mass
          react(Nreact)%type='REVER'
          react(Nreact)%comment='REVER'
          react(Nreact)%itype=11
       END IF
    END DO
    ! calculate e_rever
    e_rever=b_rever+Nrever-1
    !
    if (verbose) write(*,'(3x,a,i5)') "#(Missing endothermic reactions): ",nrever
  end subroutine add_reverse_reactions

  subroutine write_chemical_network
    !---------------------------------------------------------------------------
    !
    ! Write the chemical network in a comprehensive and
    ! machine readable format
    !
    !---------------------------------------------------------------------------
    use module_tools, only : get_file_number
    implicit none
    integer :: chemical_network,i
    !
    !
    chemical_network = get_file_number()
    open(chemical_network,file='output/chemical_network.out',status='replace',&
         access='sequential',form='formatted',action='write')
    do i=1,nreact
       write(chemical_network,format_network_write) &
            speci(react(i)%r(1:2))%name, &
            speci(react(i)%p(1:4))%name, &
            react(i)%gamma,              &
            react(i)%alpha,              &
            react(i)%beta,               &
            react(i)%itype,              &
            i,&
            trim(react(i)%comment)
    enddo
  end subroutine write_chemical_network


  subroutine write_branching_network
    use module_tools, only : get_file_number
    !---------------------------------------------------------------------------
    !
    ! Write each chemical reaction, the number of H and D atoms of
    ! each reactant and product. Elements in species%elements are:
    !
    ! H O C N He Na Mg S Si Fe  D Gr N.
    ! 1 2 3 4  5  6  7 8  9 10 11 12 13
    !
    !---------------------------------------------------------------------------
    implicit none
    integer :: branching_network,i
    !
    !
    branching_network = get_file_number()
    open(branching_network,file='output/branching_network.out',status='replace',&
         access='sequential',form='formatted',action='write')
    do i=1,nreact
       write(branching_network,format_branching) &
            speci(react(i)%r(1:2))%name,  &
            speci(react(i)%p(1:4))%name,  &
            speci(react(i)%r(1))%formula(1),  &  ! # of H atoms
            speci(react(i)%r(1))%formula(11), &  ! # of D atoms
            speci(react(i)%r(2))%formula(1),  &  ! # of H atoms
            speci(react(i)%r(2))%formula(11), &  ! # of D atoms
            speci(react(i)%p(1))%formula(1),  &  ! # of H atoms
            speci(react(i)%p(1))%formula(11), &  ! # of D atoms
            speci(react(i)%p(2))%formula(1),  &  ! # of H atoms
            speci(react(i)%p(2))%formula(11), &  ! # of D atoms
            speci(react(i)%p(3))%formula(1),  &  ! # of H atoms
            speci(react(i)%p(3))%formula(11), &  ! # of D atoms
            speci(react(i)%p(4))%formula(1),  &  ! # of H atoms
            speci(react(i)%p(4))%formula(11), &  ! # of D atoms
!!$            speci(react(i)%r(1))%formula(1:13), &  ! complete formula
!!$            speci(react(i)%r(2))%formula(1:13), &  ! complete formula
!!$            speci(react(i)%p(1))%formula(1:13), &  ! complete formula
!!$            speci(react(i)%p(2))%formula(1:13), &  ! complete formula
!!$            speci(react(i)%p(3))%formula(1:13), &  ! complete formula
!!$            speci(react(i)%p(4))%formula(1:13), &  ! complete formula
            react(i)%gamma,               &
            react(i)%alpha,               &
            react(i)%beta,                &
            react(i)%itype,               &
            i,&
            react(i)%comment
    enddo
  end subroutine write_branching_network

END MODULE MODULE_CHEM_REACT
