module h2_variables
  use module_precision
  real(kind=dp), parameter :: h2_dissoc = 4.4781  ! h2 dissociation energy (in ev)
  real(kind=dp)            :: h2_int_e            ! h2 internal energy at formation on grains (in k)
  real(kind=dp)            :: h2_energy=0.0_dp    ! calculated in compute_h2
  real(kind=dp)            :: r_gr_h2             ! [  cm3 s-1] r-factor: formation of h2 on grains
  real(kind=dp)            :: for_gr_h2           ! [ cm-3 s-1] formation of h2 on grains = r nh n(h)
  real(kind=dp)            :: cooling_h2          ! [erg/cm3/s] cooling rate of the neutral fluid due to h2 lines 
  real(kind=dp)            :: sel_tot_h2       ! [cm-3 s-1] collisional dissociation of h2 (summed over all levels)
  real(kind=dp)            :: sel_tne_h2       ! [cm-3 s-1] collisional dissociation of h2 by neutrals (summed over all levels)
  real(kind=dp)            :: sel_tio_h2       ! [cm-3 s-1] collisional dissociation of h2 by ions (summed over all levels)
  real(kind=dp)            :: sel_tel_h2       ! [cm-3 s-1] collisional dissociation of h2 by electrons (summed over all levels)
end module h2_variables
module module_h2
  use module_precision
  use h2_variables
  use module_parameters_flags
  !------------------------------------------------------------------
  ! The module 'MODULE_H2' contains variables and
  ! subroutines related to the H2 molecule, except op_H2
  ! which is in MODULE_PHYS_VAR
  !     * levels, ortho/para, reaction coefficients
  !     * Aij, quadrupole lines
  ! It contains also variables and subroutines needed to
  ! compute H2 cooling and ortho/para ratio.
  !------------------------------------------------------------------
  implicit none
  !
  ! energy levels
  ! variables defined in READ_H2_LEVELS
  integer(kind=long)            :: vmin_h2, vmax_h2 ! min. and max. values for v (vibration)
  integer(kind=long)            :: jmin_h2, jmax_h2 ! min. and max. values for j (rotation)
  logical                       :: op_lte           ! .true. if o/p h2 has been initialized to lte value
  !
  ! data type of one level
  type type_h2_level
     integer(kind=long) :: v, j    ! numbers of vibration and rotation
     real(kind=dp)      :: weight  ! statistical weight : 2j+1 (para) or 3(2j+1) (ortho)
     real(kind=dp)      :: energy  ! energy (k)
     real(kind=dp)      :: density ! density of the level (cm-3)
     real(kind=dp)      :: dens_old ! density at last call to drive
     real(kind=dp)      :: col_dens ! column density of the level (cm-2)
     real(kind=dp)      :: form_gr  ! fraction of h2 formed on grains on level v j
  end type type_h2_level
  !
  ! initialized in READ_H2_LEVELS
  type(type_h2_level),dimension(:), allocatable         :: h2_lev           ! vector containing the h2 levels
  integer(kind=long), dimension(:,:), save, allocatable :: index_vj_h2      ! index (1..nh2_lev) of 1 level (v,j)
  integer(kind=long), private, parameter                :: num_undefined=-1 ! useful to initialize index_vj_h2
  real(kind=dp)                                         :: dens_parah2      ! density of para-h2 (cm-3)
  real(kind=dp)                                         :: dens_orthoh2     ! density of ortho-h2 (cm-3)
  !
  ! data type of one H2 line
  type type_h2_line
     character(len=9)   :: name      ! name of the line
     integer(kind=long) :: nup, nlow ! index of the upper and the lower levels
     real(kind=dp)      :: aij       ! aij (s-1)
     real(kind=dp)      :: deltae    ! energy of the line : eup-elow (k)
     real(kind=dp)      :: emiss     ! emiss of the line (erg/cm3/s)
     real(kind=dp)      :: emiss_old ! emissivity at the last call to drive
     real(kind=dp)      :: intensity ! intensity integrated along the shock (erg/cm2/s/sr)
  end type type_h2_line
  !
  !--- vector containing the H2 lines (dimension NH2_lines) ---
  integer(kind=long) :: nh2_lines      ! number of h2 lines
  type (type_h2_line), dimension(:), save, allocatable :: h2_lines ! lines
  !
  ! Collision rates
  ! read in READ_H2_RATES, used in EVOLUTION_H2
  ! their dimension used are (4,NH2_lev,NH2_lev), all values outside are rejected
  real(kind=dp), private, dimension(:,:,:), allocatable :: rate_h_h2     ! collisions h-h2
  logical, private, dimension(:,:),         allocatable :: mask_h_h2     ! collisions h-h2
  real(kind=dp), private, dimension(:,:,:), allocatable :: rate_he_h2    ! collisions he-h2
  real(kind=dp), private, dimension(:,:,:), allocatable :: rate_h2_h2    ! collisions h2-h2
  real(kind=dp), private, dimension(25,0:16,0:36)       :: r_raw_gr_h2   ! collisions gr-h2 (raw rates)
  real(kind=dp), private, dimension(25,0:16,0:36)       :: d2r_raw_gr_h2 ! collisions gr-h2 (raw rates)
  real(kind=dp), private, dimension(0:16,0:36)          :: vin_gr_h2     ! corresponding velocity
  real(kind=dp), private, dimension(0:16,0:36)          :: pgr0_gr_h2    ! corresponding rate
  !
  ! evolution terms for radiative transitions
  ! read in READ_H2_LINES
  ! Aij_H2(Nup,Nlow) : Aij (s-1) of the line Nup -> Nlow (N is the index of the level)
  real(kind=dp), dimension(:,:), allocatable :: aij_h2    ! [s-1] aij_h2(nup,nlow) for the transition nup -> nlow
  real(kind=dp),dimension(:), allocatable    :: sum_aij_h2! [s-1] sum of the aij over lines from level nup
  !
  !
  ! Data from Gerlich for ortho-para transfer by H+
  ! with David Flower modification for ortho/para
  ! alternation (3/1).
  real(kind=dp), private, parameter, dimension(29) :: gerlich = &
       (/2.024d-10, 9.575d-10, 2.729d-10, 7.910d-10 &
       , 2.246d-10, 6.223d-10, 1.828d-10, 4.975d-10, 1.484d-10 &
       , 3.000d-10, 1.000d-10, 3.000d-10, 1.000d-10, 3.000d-10 &
       , 1.000d-10, 3.000d-10, 1.000d-10, 3.000d-10, 1.000d-10 &
       , 3.000d-10, 1.000d-10, 3.000d-10, 1.000d-10, 3.000d-10 &
       , 1.000d-10, 3.000d-10, 1.000d-10, 3.000d-10, 1.000d-10 /)

  !--------------------------------------------------
  ! source terms for H2 level population (cm-3.s-1)
  ! and H2 internal energy (erg.cm-3.s-1)
  !--------------------------------------------------
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: YN_rovib_H2         ! change in H2 levels number density (cm-3.s-1)

  !--------------------------------------------------------------
  ! source terms of level selective chemical reactions involving H2
  ! 22 juin 2001 - JLB
  ! One type only yet: collisional dissociation of H2
  ! Allocated in : INITIALIZE_ROVIB_H2 (H2.f90)
  !---------------------------------------------------------------
  real(kind=dp), dimension(:), save, allocatable :: sel_ch_h2        ! total selective reaction rate
  real(kind=dp), dimension(:), save, allocatable :: sel_ne_h2        ! total selective reaction rate (neutrals only)
  real(kind=dp), dimension(:), save, allocatable :: sel_io_h2        ! total selective reaction rate (ions only)
  real(kind=dp), dimension(:), save, allocatable :: sel_el_h2        ! total selective reaction rate (electrons only)
  real(kind=dp), dimension(:), save, allocatable :: sel_rx_h2        ! partial selective reaction rate


contains

  subroutine h2_init
    use module_phys_var, only: nh2_lev
    !_________________________________________________________________
    !
    ! Initialize H2 related variables
    !_________________________________________________________________
    !
    implicit none
    ! initialization
    if (do_h2) then
       call read_h2_levels      ! levels
       call initialize_rovib_h2 ! population of the levels, according to the choice in ortho:para
       call read_h2_rates       ! read collision rates for H-H2, He-H2, H2-H2
       call read_h2_lines       ! read quadrupolar lines of H2
    else
       nh2_lev = 0
       return
    endif
  end subroutine h2_init

  subroutine h2_clean
    if (do_h2) then
       deallocate(index_vj_h2)
       deallocate(h2_lines)
       deallocate(h2_lev)
       deallocate(yn_rovib_h2)
       deallocate(sel_ch_h2)
       deallocate(sel_ne_h2)
       deallocate(sel_io_h2)
       deallocate(sel_el_h2)
       deallocate(sel_rx_h2)
       deallocate(rate_h_h2)
       deallocate(mask_h_h2)
       deallocate(rate_he_h2)
       deallocate(rate_h2_h2)
       deallocate(aij_h2)
       deallocate(sum_aij_h2)
    endif
  end subroutine h2_clean

  subroutine read_h2_levels
    use module_tools, only : get_file_number
    use module_constants, only : zero,data_dir
    use module_phys_var, only: nh2_lev, nh2_lines_out
    !_________________________________________________________________
    !
    ! Initialize the H2 levels: V, J, weight, energy.
    !_________________________________________________________________
    !
    implicit none
    character(len=1) :: charact
    integer(kind=long) :: i, ii, v, j
    character(len=*), parameter :: name_file_h2_lev=data_dir//'/h2_levels_evj.in'
    integer                     :: file_h2_lev
    !
    !
    allocate (h2_lev(nh2_lev))
    h2_lev(:)%v=0
    h2_lev(:)%j=0
    h2_lev(:)%weight=zero
    h2_lev(:)%energy=zero
    h2_lev(:)%density=zero
    h2_lev(:)%dens_old=zero
    h2_lev(:)%col_dens=zero
    h2_lev(:)%form_gr=zero
    !
    ! file opening
    file_h2_lev = get_file_number()
    open(file_h2_lev,file=name_file_h2_lev,status='old',&
         access='sequential',form='formatted',action='read')
    !
    ! comments
    do i=1,6
       read(file_h2_lev,'(a1)')charact
    end do
    !
    do i=1, nh2_lev
       read(file_h2_lev,*) &
            ii, &
            h2_lev(i)%v, &
            h2_lev(i)%j, &
            h2_lev(i)%weight, &
            h2_lev(i)%energy
       ! check if all lines have been correctly read, we must have i=ii
       if (i.ne.ii) then
          write(*,*) "E-read_h2_levels, error with H2 levels",i,ii
          stop
       endif
    end do

    ! min. and max. values of V and J for these levels
    vmin_h2=minval(h2_lev(:)%v)
    vmax_h2=maxval(h2_lev(:)%v)
    jmin_h2=minval(h2_lev(:)%j)
    jmax_h2=maxval(h2_lev(:)%j)
    !
    ! Allocate memory for the table index_VJ_H2
    allocate(index_vj_h2(vmin_h2:vmax_h2, jmin_h2:jmax_h2))
    index_vj_h2(:,:)=num_undefined
    do i=1,nh2_lev
       index_vj_h2(h2_lev(i)%v,h2_lev(i)%j)=i
    enddo

    ! file closure
    if (verbose) write(*,*) "I-read_h2_levels: done"
    close(file_h2_lev)
  end subroutine read_h2_levels


  subroutine initialize_rovib_h2
    use module_phys_var, only : tn, op_h2, nh2_lev, iforh2
    use module_chemical_species, only : speci,ind_h,ind_ph2,ind_oh2,dens_oh2,dens_ph2,dens_h2
    use module_constants, only : zero, kb, everg
    use module_debug
    !___________________________________________________________________
    !
    ! Initialize population of the ro-vibrational levels of H2
    !
    ! 2 possibilities :
    ! * op_H2=op_H2 read in READ_PARAMETERS
    ! * op_H2=op LTE at temperature Tn (if op_H2 = -1)
    !
    ! Note : op_H2 is (re)-calculated
    !___________________________________________________________________
    !
    implicit none
    character(len=7) :: op_choice
    integer(kind=long) :: i
    real(kind=dp) :: zortho, zpara, weight_ortho, weight_para
    real(kind=dp),parameter :: op_h2_lte=999._dp
!!$    real(kind=dp),parameter :: population_limit=1.d-20 ! lower limit for h2 population
    real(kind=dp),parameter :: population_limit=1.d-15 ! lower limit for h2 population
    real(kind=dp) :: int_enrg, t_ini, tform, tf_min, tf_max, res
    real(dp) :: densh2
    
    allocate (yn_rovib_h2(nh2_lev))
    allocate (sel_ch_h2(nh2_lev))
    allocate (sel_ne_h2(nh2_lev))
    allocate (sel_io_h2(nh2_lev))
    allocate (sel_el_h2(nh2_lev))
    allocate (sel_rx_h2(nh2_lev))
    yn_rovib_h2 = zero
    sel_ch_h2 = zero
    sel_ne_h2 = zero
    sel_io_h2 = zero
    sel_el_h2 = zero
    sel_rx_h2 = zero

    ! which o/p to take ?
    if (op_h2.eq.-1.0) then
       op_lte=.true.
    else
       op_lte=.false.
    end if
    !
    ! computes the populations, given the o/p choice
    if (op_lte) then
       ! LTE at temperature Tn
       h2_lev(:)%density=h2_lev(:)%weight * &
            exp(-(h2_lev(:)%energy-h2_lev(index_vj_h2(0,0))%energy)/tn)
    else
       ! o/p from the input file
       ! Partition functions and weights for ortho-H2 and para-H2
       ! oH2
       Zortho=Zero
       T_ini = Tn
       do i=index_vj_h2(0,1),nh2_lev,2 ! first ortho level : (v=0,j=1)
          zortho = zortho + h2_lev(i)%weight * &
               exp(-(h2_lev(i)%energy-h2_lev(index_vj_h2(0,0))%energy)/tn)
       end do
       weight_ortho=op_H2/(1._DP+op_H2)/Zortho
       !
       zpara=zero
       do i=index_vj_h2(0,0),nh2_lev,2 ! first para level : (v=0,j=0)
          zpara = zpara + h2_lev(i)%weight * &
               exp(-(h2_lev(i)%energy-h2_lev(index_vj_h2(0,0))%energy)/tn)
       end do
       weight_para=1._DP/(1._DP+op_H2)/Zpara
       !
       ! populations for ortho levels
       where (mod(h2_lev(:)%j,2) > 0)
          h2_lev(:)%density=weight_ortho * &
               h2_lev(:)%weight * &
               exp(-(h2_lev(:)%energy-h2_lev(index_vj_h2(0,0))%energy)/tn)
       end where
       !
       ! populations for para levels
       where (mod(h2_lev(:)%j,2) == 0)
          h2_lev(:)%density=weight_para * &
               h2_lev(:)%weight * &
               exp(-(h2_lev(:)%energy-h2_lev(index_vj_h2(0,0))%energy)/tn)
       end where
    end if
    !
    ! avoids too small numbers => set lower limit to the populations
    where (h2_lev(:)%density < population_limit)
       h2_lev(:)%density=population_limit
    end where

    if (debug) then
       write(*,'("D-rovib_h2, Initialization of H2 weights:",12es12.1)') h2_lev%weight
       write(*,'("D-rovib_h2, Initialization of H2 levels:",12es12.1)') H2_lev%density
    endif

!  Formation of H2 on grains :
!   0: 1/3 of 4.4781 eV in internal energy (=> 17249 K) (Allen, 1999)
!   1: Proportional to Boltzman distribution at 17249 K
!   2: Dissociation limit : v = 14, J = 0,1 (4.4781 eV)
!   3: v = 6, J = 0,1
!   4: fraction = relative populations at t, initialised as H2_lev%density
!                 and changed during integration

    IF (iforH2 == 0) then

      int_enrg = H2_dissoc * EVerg / (3.0_DP * kB)
      tform = int_enrg
      tf_min = 0.0_dp
      tf_max = 3.0_DP * int_enrg
      res = 1.0_dp

      do while (abs(res) > 1.0e-3_DP)
        res = SUM(H2_lev%weight * (int_enrg-H2_lev%Energy) &
            * EXP(-H2_lev%Energy/tform))
        if (res < 0.0_dp) then
          tf_max = tform
          tform = 0.5_DP * (tf_min + tf_max)
        else
          tf_min = tform
          tform = 0.5_DP * (tf_min + tf_max)
        endif
      end do

      int_enrg = tform

      H2_lev%Form_gr = H2_lev%weight * &
           EXP(-(H2_lev%Energy-H2_lev(index_VJ_H2(0,0))%Energy)/int_enrg)
      H2_lev%Form_gr = H2_lev%Form_gr / SUM(H2_lev%Form_gr)

    ELSE IF (iforH2 == 1) then

      int_enrg = H2_dissoc * EVerg / (3.0_DP * kB)
      H2_lev%Form_gr = H2_lev%weight * &
           EXP(-(H2_lev%Energy-H2_lev(index_VJ_H2(0,0))%Energy)/int_enrg)
      H2_lev%Form_gr = H2_lev%Form_gr / SUM(H2_lev%Form_gr)

    ELSE IF (iforH2 == 2) then

      H2_lev(index_VJ_H2(14,0))%Form_gr = 0.25_DP
      H2_lev(index_VJ_H2(14,1))%Form_gr = 0.75_DP

    ELSE IF (iforH2 == 3) then

      H2_lev(index_VJ_H2(6,0))%Form_gr = 0.25_DP
      H2_lev(index_VJ_H2(6,1))%Form_gr = 0.75_DP

    ELSE IF (iforH2 == 4) then

      H2_lev%Form_gr = H2_lev%density
      H2_lev%Form_gr = H2_lev%Form_gr / SUM(H2_lev%Form_gr)

    ENDIF

    H2_int_E = SUM(H2_lev%Form_gr * H2_lev%Energy)

    ! normalization at H2 density
    densh2 = speci(ind_ph2)%density+speci(ind_oh2)%density+speci(ind_h)%density/2._dp
    if (densh2.eq.0.) then
       h2_lev(:)%density=1.d-8
    else
       h2_lev(:)%density=h2_lev(:)%density*densh2/(sum(dble(h2_lev(:)%density)))
    endif
    ! op_H2 is (re)-calculated for checking
    call compute_op_h2
    if (op_lte) write(*,'("--- LTE chosen for ortho:para-H2, &
         &o/p = ",ES7.1," ---")')op_H2

  END SUBROUTINE INITIALIZE_ROVIB_H2



  SUBROUTINE READ_H2_RATES
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    read collision rates for H-H2, He-H2 and H2-H2.
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    rate_H_H2, rate_He_H2, rate_H2_H2
    !---------------------------------------------------------------------------
    use module_tools, only : get_file_number
    use module_constants, only : zero,data_dir
    use module_phys_var, only : nh2_lev, h_h2_flag
    use num_rec, only: spline
    implicit none
    !
    !
    integer(kind=long) :: nrate, i, levu, levl, vup, jup, vlow, jlow
    integer                     :: file_h_h2
    character(len=*), parameter :: name_file_he_h2=data_dir//'/coeff_he_h2.in'
    character(len=*), parameter :: name_file_h2_h2=data_dir//'/coeff_h2_h2.in'
    character(len=*), parameter :: name_file_gr_h2=data_dir//'/coeff_gr_h2.in'
    character(len=256)          :: name_file_h_h2 =data_dir//'/coeff_h_h2_drf.in' ! Default reaction coefficients
    integer                     :: file_he_h2
    integer                     :: file_h2_h2
    integer                     :: file_gr_h2
    integer :: ii, jj, kk
    integer :: iju, ivu
    REAL(KIND=DP) :: yp1, ypn
    REAL(KIND=DP) :: toto1, toto2, toto3, toto4
    REAL(KIND=DP), DIMENSION(1:25) :: vin
    data yp1, ypn / 1.0d30, 1.0d30 /
    data vin / 2.d0,4.d0,6.d0,8.d0,&
               10.d0,12.d0,14.d0,16.d0,18.d0,&
               20.d0,22.d0,24.d0,26.d0,28.d0,&
               30.d0,32.d0,34.d0,36.d0,38.d0,&
               40.d0,42.d0,44.d0,46.d0,48.d0,50.d0/


    ! initialization
    ALLOCATE (rate_H_H2(4,NH2_lev,NH2_lev))
    ALLOCATE (mask_H_H2(NH2_lev,NH2_lev))
    ALLOCATE (rate_He_H2(4,NH2_lev,NH2_lev))
    ALLOCATE (rate_H2_H2(4,NH2_lev,NH2_lev))
    rate_H_H2(1,:,:)=-50._DP
    rate_H_H2(2:4,:,:)=Zero
    mask_H_H2(:,:)= .FALSE.
    rate_He_H2(1,:,:)=-50._DP
    rate_He_H2(2:4,:,:)=Zero
    rate_H2_H2(1,:,:)=-50._DP
    rate_H2_H2(2:4,:,:)=Zero
    r_raw_GR_H2(:,:,:)=Zero
    d2r_raw_GR_H2(:,:,:)=Zero
    vin_GR_H2(:,:)=0

    !-----------------------
    ! rate H-H2
    !-----------------------

    ! We read either DRF file or MM file or BOTH
    ! If BOTH, then first MM then DRF so that Quantum rates override Semiclassical ones

    if (H_H2_flag == "MM" .OR. H_H2_flag == "BOTH") then
      name_file_H_H2 = data_dir//'/coeff_H_H2_MM.in'
      file_H_H2 = GET_FILE_NUMBER()
      OPEN(file_H_H2,file=name_file_H_H2,status='OLD',&
           access='SEQUENTIAL',form='FORMATTED',action='READ')
      DO i=1,4
         READ(file_H_H2,*)
      ENDDO

      READ(file_H_H2,*)Nrate
      DO i=1, Nrate
         READ(file_H_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
         IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
           rate_H_H2(1,LevU,LevL) = toto1
           rate_H_H2(2,LevU,LevL) = toto2
           rate_H_H2(3,LevU,LevL) = toto3
           rate_H_H2(4,LevU,LevL) = toto4
         ENDIF
      END DO
      CLOSE(file_H_H2)
    endif

    if (H_H2_flag == "DRF" .OR. H_H2_flag == "BOTH") then
      name_file_h_h2 = data_dir//'/coeff_h_h2_drf.in'
      file_h_h2 = get_file_number()
      open(file_h_h2,file=name_file_h_h2,status='old',&
           access='sequential',form='formatted',action='read')
      DO i=1,4
         READ(file_H_H2,*)
      ENDDO

      READ(file_H_H2,*)Nrate
      DO i=1, Nrate
         READ(file_H_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
         IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
           rate_H_H2(1,LevU,LevL) = toto1
           rate_H_H2(2,LevU,LevL) = toto2
           rate_H_H2(3,LevU,LevL) = toto3
           rate_H_H2(4,LevU,LevL) = toto4
           mask_H_H2(LevU,LevL) = .TRUE.
         ENDIF
      END DO
      CLOSE(file_H_H2)
    endif

    !-------------
    ! rate He-H2
    !-------------
    ! file opening
    file_He_H2 = GET_FILE_NUMBER()
    OPEN(file_He_H2,file=name_file_He_H2,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    ! commentss
    DO i=1,4
       READ(file_He_H2,*)
    ENDDO

    ! number of collision rates
    READ(file_He_H2,*)Nrate
    ! read des coefficients ligne par ligne
    DO i=1, Nrate
       READ(file_He_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
       IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
         rate_He_H2(1,LevU,LevL) = toto1
         rate_He_H2(2,LevU,LevL) = toto2
         rate_He_H2(3,LevU,LevL) = toto3
         rate_He_H2(4,LevU,LevL) = toto4
       ENDIF
    END DO
    ! file closure
    CLOSE(file_He_H2)

    !-------------
    ! rate H2-H2
    !-------------
    ! file opening
    file_H2_H2 = GET_FILE_NUMBER()
    OPEN(file_H2_H2,file=name_file_H2_H2,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    ! comments
    DO i=1,4
       READ(file_H2_H2,*)
    ENDDO

    ! number of collision rates
    READ(file_H2_H2,*)Nrate
    ! read des coefficients ligne par ligne
    DO i=1, Nrate
       READ(file_H2_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
       IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
         rate_H2_H2(1,LevU,LevL) = toto1
         rate_H2_H2(2,LevU,LevL) = toto2
         rate_H2_H2(3,LevU,LevL) = toto3
         rate_H2_H2(4,LevU,LevL) = toto4
       ENDIF
    END DO
    ! file closure
    CLOSE(file_H2_H2)

    ! Raw rates - file opening Para, then Ortho
    file_GR_H2 = GET_FILE_NUMBER()
    open(file_gr_h2,file=data_dir//"/ph2gr.in",status='old',&
         access='sequential',form='formatted',action='read')
    do ii=1,25
       read(file_gr_h2,*)
       do kk=0,36,2
         read(file_gr_h2,*) iju, (r_raw_gr_h2(ii,jj,kk),jj=0,16)
       end do
    end do
    ! file closure
    close(file_gr_h2)

    file_gr_h2 = get_file_number()
    open(file_gr_h2,file=data_dir//"/oh2gr.in",status='old',&
         access='sequential',form='formatted',action='read')
    do ii=1,25
       read(file_gr_h2,*)
       do kk=1,35,2
         read(file_gr_h2,*) iju, (r_raw_gr_h2(ii,jj,kk),jj=0,16)
       end do
    end do
    ! file closure
    close(file_gr_h2)
    !
    do ii=3,nh2_lev
      iju = h2_lev(ii)%j
      ivu = h2_lev(ii)%v
      do jj=1,25
        if (r_raw_gr_h2(jj,ivu,iju) > 0.0_dp) then
          vin_gr_h2(ivu,iju) = vin(jj)
          pgr0_gr_h2(ivu,iju) = r_raw_gr_h2(jj,ivu,iju)
          exit
        endif
      end do
      call spline (vin(jj:),r_raw_gr_h2(jj:,ivu,iju),yp1,ypn,d2r_raw_gr_h2(jj:,ivu,iju))
    end do
  end subroutine read_h2_rates


  SUBROUTINE READ_H2_LINES
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    reads H2 lines (levels, energies, Aij).
    ! subroutine/function needed :
    !    NAME_H2_LINE
    ! input variables :
    ! ouput variables :
    ! results :
    !    Aij_H2, sum_Aij_H2
    !---------------------------------------------------------------------------
    use module_tools, only : get_file_number
    use module_constants, only : zero,data_dir
    use module_phys_var, only : nh2_lev, nh2_lines_out
    implicit none
    !
    !
    real(kind=dp) :: aij_s, aij_q, aij_o
    integer(kind=long) :: i, vup, vlow, jup, jlow, nup, nlow, line_number
    character(len=1) :: line_type
    character(len=*), parameter :: name_file_aij_h2=data_dir//'/aij_h2.in'
    character(len=*), parameter :: format_aij_h2='(3i4,3d15.6)'
    integer                     :: file_aij_h2
    !
    ! initialization
    ALLOCATE (Aij_H2(NH2_lev,NH2_lev), sum_Aij_H2(NH2_lev))
    Aij_H2(:,:)=Zero
    sum_Aij_H2(:)=Zero

    ! file opening
    file_aij_h2 = get_file_number()
    open(file_aij_h2,file=name_file_aij_h2,status='old',&
         access='sequential',form='formatted',action='read')
    ! comments
    do i=1,5
       read(file_aij_h2,*)
    end do

    ! read line per line. There is max. 3 Aij per level (lines S, Q, O)
    ! the file is in increasing Vup order
    Vup=0
    DO WHILE (Vup <= Vmax_H2)
       READ(file_Aij_H2,format_Aij_H2) &
            Vup, Vlow, Jup, Aij_S, Aij_Q, Aij_O
       ! chock if the level is in the model
       IF ((Vup <= Vmax_H2) .AND. (Vlow <= Vmax_H2) .AND. (Jup <=Jmax_H2)) THEN
          Nup=index_VJ_H2(Vup,Jup) ! index of the upper level
          IF ((Nup > 0) .AND. (Nup <= NH2_lev)) THEN
             ! S line
             Jlow=Jup-2
             IF ((Jlow >=0) .AND. (Aij_S /= Zero)) THEN
                Nlow=index_VJ_H2(Vlow,Jlow)
                IF ((Nlow > 0) .AND. (Nlow <= NH2_lev)) Aij_H2(Nup,Nlow)=Aij_S
             ENDIF
             ! Q line
             Jlow=Jup
             IF (Aij_Q /= Zero) THEN
                Nlow=index_VJ_H2(Vlow,Jlow)
                IF ((Nlow > 0) .AND. (Nlow <= NH2_lev)) Aij_H2(Nup,Nlow)=Aij_Q
             ENDIF
             ! O line
             Jlow=Jup+2
             IF ((Jlow <=JMAX_H2) .AND. (Aij_O /= Zero)) THEN
                Nlow=index_VJ_H2(Vlow,Jlow)
                IF ((Nlow > 0) .AND. (Nlow <= NH2_lev)) Aij_H2(Nup,Nlow)=Aij_O
             ENDIF
          ENDIF
       ENDIF
    END DO

    ! file closure
    CLOSE(file_Aij_H2)

    ! computes the sum of the Aij of the lines starting from the level Nup
    ! starts at (V=0,J=2)
    DO Nup = index_VJ_H2(0,2), NH2_lev
       sum_Aij_H2(Nup)=SUM(Aij_H2(Nup,1:Nup-1))
    END DO
    ! counts the number of transitions in the model (Aij_H2 > 0)
    NH2_lines=COUNT(mask=(Aij_H2 > Zero))
    NH2_lines_out = min(NH2_lines,NH2_lines_out)
    ! allocation and initialization of H2_lines
    ALLOCATE(H2_lines(NH2_lines))
    H2_lines(:)%name=''
    H2_lines(:)%Nup=0
    H2_lines(:)%Nlow=0
    H2_lines(:)%emiss=Zero
    H2_lines(:)%emiss_old=Zero
    H2_lines(:)%intensity=Zero
    H2_lines(:)%DeltaE=Zero
    H2_lines(:)%Aij=Zero

    ! vector H2_lines is filled in order of increasing Eup
    line_number=0
    DO Nup=index_VJ_H2(0,2),NH2_lev
       DO Nlow=1, Nup-1
          IF (Aij_H2(Nup,Nlow) > Zero) THEN ! tests if line exists
             line_number=line_number+1
             Vup=H2_lev(Nup)%V
             Jup=H2_lev(Nup)%J
             Vlow=H2_lev(Nlow)%V
             Jlow=H2_lev(Nlow)%J
             ! finds the type of the line (S, Q, O)
             SELECT CASE(Jup-Jlow)
             CASE (2)
                line_type='S'
             CASE (0)
                line_type='Q'
             CASE (-2)
                line_type='O'
             CASE DEFAULT
             END SELECT
             ! fills H2_lines
             H2_lines(line_number)%name=NAME_H2_LINE(Vup=Vup,Vlow=Vlow,Jlow=Jlow,&
                  line_type=line_type)
             H2_lines(line_number)%Nup=Nup
             H2_lines(line_number)%Nlow=Nlow
             H2_lines(line_number)%Aij=Aij_H2(Nup,Nlow)
             H2_lines(line_number)%DeltaE=H2_lev(Nup)%energy-H2_lev(Nlow)%energy
          ENDIF
       END DO
    END DO

  END SUBROUTINE READ_H2_LINES


  subroutine compute_op_h2
    use module_phys_var, only : op_h2
    !---------------------------------------------------------------------------
    ! called by :
    !     INIT_ROVIB_H2
    !     DIFFUN
    ! purpose :
    !     Computes ortho:para H2 ratio from H2_lev%density. Densities of
    !     ortho-H2 and para-H2 are also calculated.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     op_H2, Dens_orthoH2, Dens_paraH2
    !---------------------------------------------------------------------------
    implicit none
    !
    if (.not.do_h2) return
    dens_parah2  = sum(h2_lev(:)%density,mask=(mod(h2_lev(:)%j,2)==0))
    dens_orthoh2 = sum(h2_lev(:)%density,mask=(mod(h2_lev(:)%j,2)>0))
    op_h2        = dens_orthoh2/dens_parah2
  end subroutine compute_op_h2


  FUNCTION NAME_H2_LINE(Vup,Vlow,Jlow,line_type) RESULT(name)
    !---------------------------------------------------------------------------
    ! called by :
    !    READ_H2_LINES
    ! purpose :
    !    computes the name of one line of the form 1-0S(1)
    ! subroutine/function needed :
    ! input variables :
    !    * Vup, Jup  -> (V,J) of the upper level  (two digits : < 100)
    !    * Jlow      -> J of th lower level       (two digits : < 100)
    !    * line_type -> 'S', 'O', 'Q'
    ! output variables :
    ! results :
    !    name -> name of the line, without blanks
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG), INTENT(in) :: Vup, Vlow, Jlow
    CHARACTER(LEN=1), INTENT(in)   :: line_type ! 'S', 'O' or 'Q'
    CHARACTER(LEN=2) :: str_Vup, str_Vlow,str_Jlow
    CHARACTER(LEN=9) :: name

    WRITE(str_Vup,'(I2)')Vup
    WRITE(str_Vlow,'(I2)')Vlow
    WRITE(str_Jlow,'(I2)')Jlow

    ! return the name without blanks
    name = TRIM(ADJUSTL(str_Vup)) // '-' // TRIM(ADJUSTL(str_Vlow))  // line_type // &
         '(' // TRIM(ADJUSTL(str_Jlow)) // ')'

  END FUNCTION NAME_H2_LINE



  subroutine evolution_h2
    use module_phys_var, only : tn, abs_deltav, nh2_lev, distance, h_h2_flag
    use module_grains, only : rgrain2, dens_grain
    use module_chemical_species, only : dens_h, dens_ph2, dens_oh2, dens_hplus, dens_he
    use module_constants, only: pi, zero
    use num_rec, only: splint
    use module_debug
    !---------------------------------------------------------------------------
    ! called by :
    !    COMPUTE_H2
    ! purpose :
    !    Calculate the source term (cm-3.s-1) for rovibrational levels of H2,
    !    using Aij_H2 (radiation, s-1) read in READ_H2_LINE and Cij_H2
    !    (collisions, s-1) calculated in this subroutine.
    !    The result, YN_rovib_H2, is used in COMPUTE_H2 and in DIFFUN.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !    YN_rovib_H2
    !---------------------------------------------------------------------------
    !
    implicit none
    real(kind=dp), dimension(nh2_lev,nh2_lev) :: cij_h2          ! probability of excitation by collision (s-1)
    real(kind=dp), dimension(nh2_lev)         :: sum_cij_h2      ! sum of cij_h2 (s-1) starting from a given level
    real(kind=dp), dimension(nh2_lev,nh2_lev) :: aij_plus_cij_h2 ! sum of the radiative and collisional terms ---
    integer(kind=long) :: nup, jup, vup                          ! upper level
    integer(kind=long) :: nlow, jlow, vlow                       ! lower level
    integer(kind=long) :: abs_deltaj, level1, level2, i
    real(kind=dp) :: deexcit_to_excit, excit_rate, deexcit_rate
    real(kind=dp) :: decal, tdd, tdd2, gam, gam1, gam2, delta_e
    real(kind=dp) :: coeff_h_nr, coeff_h_r, coeff_h, coeff_he, coeff_h2, coeff_hplus
    real(kind=dp) :: dv_kms, vth, prob_excit, excit_grain, vin_h2
    real(kind=dp), dimension(1:25) :: vin
    !
    data vin / 2.d0,4.d0,6.d0,8.d0,&
               10.d0,12.d0,14.d0,16.d0,18.d0,&
               20.d0,22.d0,24.d0,26.d0,28.d0,&
               30.d0,32.d0,34.d0,36.d0,38.d0,&
               40.d0,42.d0,44.d0,46.d0,48.d0,50.d0/
    !
    !--- initialization ---
    cij_h2          = zero ! collisions (2 dimensions)
    sum_cij_h2      = zero ! collisions (1 dimension)
    aij_plus_cij_h2 = zero ! collisions + radiation (2 dimensions)
    yn_rovib_h2     = zero ! source terms for h2-level populations (cm-3.s-1)
    !
    !
    ! A. Calculate the probability of excitation by collisions
    !    Cij_H2 (s-1)
    do nup = 2, nh2_lev ! loop on upper levels
       jup = h2_lev(nup)%j
       vup = h2_lev(nup)%v
       do nlow = 1, nup-1 ! loop on lower levels
          jlow = h2_lev(nlow)%j
          vlow = h2_lev(nlow)%v
          abs_deltaj = abs(jup - jlow)
          delta_e = h2_lev(nup)%energy - h2_lev(nlow)%energy
          !
          !--- deexcit_to_excit is the Boltzmann factor allowing to ---
          !--- get the excitation rate from the de-excitation rate. ---
          deexcit_to_excit = exp(-delta_e/tn) * h2_lev(nup)%weight/ h2_lev(nlow)%weight
          !
          ! 1) fit of H2 + H collision rates.
          !    Tdd is a reduced variable that can be "shifted" with a different value
          !    for reactive collisions (delta(J) odd) and for non-reactive ones.
          !    This allows for different asymptotic forms at low temperature.
          if (h_h2_flag == "DRF") then
            decal = 1.0_dp
          else
            if (mod(abs_deltaj,2) == 0) then
              decal = 1.0_dp
            else
              decal = 0.3_dp
            endif
          endif
          tdd = decal + tn * 1.0d-3
          tdd2 = tdd*tdd
          !
          ! 1.1) non-reactive collisions (always possible)
          !      Warning: reactive collision are mixed in Martin & Mandy,
          !               and are thus included in (badly-named) coeff_H_NR
          gam = rate_h_h2(1,nup,nlow)           &
               + rate_h_h2(2,nup,nlow) / tdd     &
               + rate_h_h2(3,nup,nlow) / tdd2    &
               + rate_h_h2(4,nup,nlow) * tdd
          coeff_h_nr = 10.0_dp**gam
          coeff_h_r = 0.0_dp
          !
          ! 1.2) reactive collisions are added (David Flower's rule) if H_H2_flag = DRF
          !      and for Delta J = 2 only if H_H2_flag =  BOTH (complement quantum rates)
          if (h_h2_flag == "DRF" .and. vup /= vlow) then
             if (mod(abs_deltaj,2) == 0) then
                coeff_h_r = 10.0_dp**gam * &
                     exp(-max(0.0_dp,(3900.0_dp-delta_e)/tn))
             else
                if (jlow > 0) then
                   level1 = index_vj_h2(vlow,jlow-1)
                else
                   level1 = index_vj_h2(vlow,jlow+1)
                endif
                level2 = index_vj_h2(vlow,jlow+1)
                if (level2 == num_undefined) then
                   level2 = index_vj_h2(vlow,jlow-1)
                endif

                gam1 = rate_H_H2(1,Nup,Level1)         &
                     + rate_H_H2(2,Nup,Level1) / Tdd   &
                     + rate_H_H2(3,Nup,Level1) / Tdd2  &
                     + rate_H_H2(4,Nup,Level1) * Tdd
                gam2 = rate_H_H2(1,Nup,Level2)         &
                     + rate_H_H2(2,Nup,Level2) / Tdd   &
                     + rate_H_H2(3,Nup,Level2) / Tdd2  &
                     + rate_H_H2(4,Nup,Level2) * Tdd
                coeff_h_r = (10.0_dp**gam1 + 10.0_dp**gam2) * 0.5_dp * &
                     exp(-max(0.0_dp,(3900.0_dp-delta_e)/tn))
                if (mod(jup,2) == 1) then
                   coeff_h_r = coeff_h_r / 3.0_dp
                endif
             endif
          endif
          if (h_h2_flag == "both" .and. vup /= vlow  .and. mod(abs_deltaj,2) == 0 .and. mask_h_h2(nup,nlow)) then
             coeff_h_r = 10.0_dp**gam * exp(-max(0.0_dp,(3900.0_dp-delta_e)/tn))
          endif
          !
          ! 1.3) ortho-para transfer with H from Schofield (only if MM are not included)
          !    + David Flower prescription (only for Delta v=0)
          if (h_h2_flag == "DRF" &
               .or. (h_h2_flag == "BOTH" .and. mod(abs_deltaj,2) == 0 .and. mask_h_h2(nup,nlow))) then
             if (vup == vlow .and. abs_deltaj < 3) then
                if (mod(jup,2) == 1 .and. mod(abs_deltaj,2) == 1) then
                   coeff_h_r = 8.0d-11 * exp(-3900.0_dp/tn) / 3.0_dp
                else
                   coeff_h_r = 8.0d-11 * exp(-3900.0_dp/tn)
                endif
             endif
          endif
          coeff_H = coeff_H_NR + coeff_H_R

          ! if((Nup.eq.26).and.((Nlow.le.21).or.(Nlow.eq.24))) then 
          !  print *, coeff_H_NR, coeff_H_R
          ! endif

          !--------------------------------------------------------------------------
          ! 2) fit of H2 + He collision rates.
          !    Tdd is a reduced variable that can be "shifted" with a different value
          !    for reactive collisions (delta(J) odd) and for non-reactive ones.
          !    This allows for different asymptotic forms at low temperature.
          !--------------------------------------------------------------------------
          gam = rate_He_H2(1,Nup,Nlow)          &
               + rate_He_H2(2,Nup,Nlow) / Tdd   &
               + rate_He_H2(3,Nup,Nlow) / Tdd2  &
               + rate_He_H2(4,Nup,Nlow) * Tdd
          coeff_He = 10.0_DP**gam

          !--------------------------------------------------------------------------
          ! 3) fit of H2 + para-H2 collision rates.
          !    Tdd is a reduced variable that can be "shifted" with a different value
          !    for reactive collisions (delta(J) odd) and for non-reactive ones.
          !    This allows for different asymptotic forms at low temperature.
          !--------------------------------------------------------------------------
          gam = rate_H2_H2(1,Nup,Nlow)          &
               + rate_H2_H2(2,Nup,Nlow) / Tdd  &
               + rate_H2_H2(3,Nup,Nlow) / Tdd2 &
               + rate_H2_H2(4,Nup,Nlow) * Tdd
          coeff_H2 = 10.0_DP**gam

          !------------------------------------------------
          ! 4) ortho-para transfer with ions (Gerlich data)
          !    only in v=0
          !------------------------------------------------
          IF (Vup == 0 .AND. ABS_deltaJ == 1) THEN
             coeff_Hplus = Gerlich(Jup)
          ELSE
             coeff_Hplus = 0.0_DP
          ENDIF

          !---------------------------
          ! 5) Calculate Cij_H2 (s-1)
          !---------------------------
          deexcit_rate= coeff_H2 * (dens_ph2+dens_oh2) &
               + coeff_He * Dens_He        &
               + coeff_Hplus * Dens_Hplus  &
               + coeff_H * Dens_H
          excit_rate  = deexcit_rate * deexcit_to_excit

          !----------------------------------
          ! 6) Calculate excitation by grains
          !----------------------------------
          ! Excitation probabilites are fitted by a spline,
          ! see f77 , for DeltaV > VIN (first impact
          ! velocity at which excitation probability becomes
          ! non zero), VTH = VIN-2km/s.  Exponential form at
          ! DeltaV < VIN, scaled by PRG0 to pass through
          ! point at VIN Raw excitation rates are in
          ! r_raw_GR_H2 Second derivatives, as computed by
          ! SPLINE are in d2r_raw_GR_H2
          dv_kms = ABS_DeltaV  * 1.0D-5
          Vin_H2 = vin_GR_H2(Vup,Jup)
          Vth = Vin_H2 - 2.0_DP

          excit_grain = 0.0_DP
          prob_excit = 0.0_DP
!         if ((Vlow==0) .and. ((Jlow==0) .or. (Jlow==1)) .and. (mod(Jup-Jlow,2)==0)) then
          if ((Vlow==0) .and. (Jlow <= 7) .and. (mod(Jup-Jlow,2)==0)) then
            if (dv_kms >= Vin_H2) then
               prob_excit = splint(vin,r_raw_GR_H2(:,Vup,Jup),d2r_raw_GR_H2(:,Vup,Jup),dv_kms)
               prob_excit = max (prob_excit,0.0d0)

            else if (dv_kms > 1.0d-20) then
!
! 8.02*dv_kms**2 = 0.1*Ts where 3/2kTs = 1/2mH2 dv**2
! extrapolation to velocity difference lower than the first point in the excitation grid
! at Vin = vin_GR_H2(Vup,Jup)

               prob_excit = EXP(delta_E/8.02*(1.0_DP/vin_GR_H2(Vup,Jup)**2 - 1.0_DP/dv_kms**2))
               prob_excit = prob_excit*pgr0_GR_H2(Vup,Jup)

            else

               prob_excit = 0.0_DP

            end if
            excit_grain = prob_excit * Dens_grain * pi * rgrain2 * ABS_DeltaV

          else
            excit_grain = 0.0_DP
          end if

          Cij_H2(Nlow,Nup) = Cij_H2(Nlow,Nup) + excit_rate + excit_grain
          Cij_H2(Nup,Nlow) = Cij_H2(Nup,Nlow) + deexcit_rate

       END DO ! end of loop on Nlow
    END DO ! end of loop on Nup
    ! stop

    !--- calculate sum_Cij_H2 ---
    DO Nup=1,NH2_lev
       sum_Cij_H2(Nup) = SUM(Cij_H2(Nup,1:NH2_lev))
    END DO

    !---------------------------------------------------------------
    ! B. Add (de)excitation probabilities from collisions (Cij_H2)
    !    and from radiation (Aij_H2)
    !---------------------------------------------------------------
    aij_plus_cij_h2 = aij_h2 + cij_h2 ! two dimensions
    do i=1,nh2_lev ! diagonal terms
       aij_plus_cij_h2(i,i) = aij_plus_cij_h2(i,i) - sum_aij_h2(i) - sum_cij_h2(i)
    end do

    !----------------------------------------------------------
    ! C. Calculate evolution terms for populations in cm-3.s-1
    !    = matrix multiplication of population density and
    !    (de)excitation probabilities.
    !----------------------------------------------------------
    yn_rovib_h2(1:nh2_lev) = matmul(h2_lev(1:nh2_lev)%density , aij_plus_cij_h2(1:nh2_lev,1:nh2_lev))

!!$    write(*,'(a,i6,/,(10es12.3))') ">>>>>>>>>> EVOLUTION_H2: ",nh2_lev,yn_rovib_h2(1:nh2_lev)
    yn_rovib_h2 = 0.d0
    ! PHB: why do we compute it if we set it to zero ?????
  end subroutine evolution_h2



  subroutine compute_h2
    use module_constants, only : kb
    use module_phys_var, only : nh2_lev
    use module_debug
    !---------------------------------------------------------------------------
    ! called by :
    !    SOURCE
    ! purpose :
    !    Calculate the source terms for H2 levels (cm-3.s-1), H2 internal energy
    !    (erg.cm-3.s-1), the H2 lines (erg.cm-3.s-1) and the corresponding
    !    cooling rate (erg.cm-3.s-1).
    ! subroutine/function needed :
    !    EVOLUTION_H2
    ! input variables :
    ! output variables :
    ! results :
    !    YN_rovib_H2, H2_Energy, H2_lines%emiss, cooling_H2
    !---------------------------------------------------------------------------
    implicit none
    !
    ! Compute H2 levels and related quantities only if required
    if (.not.do_h2) then
       yn_rovib_h2 = 0.d0
       h2_energy = 0.d0
       h2_lines%emiss = 0.d0
       cooling_h2 = 0.d0
       return
    endif

    !----------------------------------------------------------
    ! In EVOLUTION_H2, the source terms YN_rovib_H2 (cm-3.s-1)
    ! is calculated, with treatment of radiation and collision
    !----------------------------------------------------------
    CALL EVOLUTION_H2
    !
    !-----------------------------------------
    ! Calculate the emissivity (erg.cm-3.s-1)
    ! of each H2 quadrupolar line.
    ! result in H2_lines%emiss
    !-----------------------------------------
    ! first in K.cm-3.s-1
    H2_lines(1:NH2_lines)%emiss = &
         H2_lines(1:NH2_lines)%DeltaE * H2_lines(1:NH2_lines)%Aij * &
         H2_lev(H2_lines(1:NH2_lines)%Nup)%density
    ! change into erg.cm-3.s-1
    H2_lines%emiss=H2_lines%emiss*kB

    !-------------------------------------------------------------
    ! cooling rate cooling_H2 (erg.cm-3.s-1) = sum of emissivities
    !-------------------------------------------------------------
    cooling_H2=SUM(DBLE(H2_lines(1:NH2_lines)%emiss))

    !---------------------------------------------------------------
    ! source term for internal energy of H2 H2_energy (erg.cm-3.s-1)
    !---------------------------------------------------------------
    H2_energy = kB * SUM(YN_rovib_H2(1:NH2_lev) * H2_lev(1:NH2_lev)%energy)

  END SUBROUTINE COMPUTE_H2

END MODULE MODULE_H2
