!___________________________________________________________________
!
!
!                    T h e  C O L L A P S E   m o d e l
!
!
! history
!
!    * hiver/printemps 2001/2002 : DRF - GPdF - JLB
!        - Correction de nombreux bugs
!        - Chocs J
!
!    * mai/juin 2001 : Jacques Le Bourlot
!        - Remplacement de GEAR par DVODE (d'apres Pierre Hily-Blant)
!        - Revision generale (avec allegement)
!        - Correction de quelques bugs mineurs
!        - Suppression des variables 11 et 12
!          (Masse et rayon des grains
!
!    * nov 2000 : David Wilgenbus
!        - correction of grain compression
!        - H2O and CO cooling taken from Neufeld & Kaufman 1993
!
!    * oct 2000 : David Wilgenbus
!          translation into fortran 90 and revision of the mhdC code
!          from G. Pineau des Forets and D. Flower
!
!___________________________________________________________________
!
program collapse
  use module_grains,            only : md_grain, r_grain, dens_grain
  use module_constants,         only : yearsec,pi,zero,parsec,mp,kb,everg
  use module_initialize
  use module_chemical_species
  use module_h2!,                only : h2_lev,h2_lines,nh2_lines,index_vj_h2,op_lte
  use module_phys_var
  use dvode_f90_m
  use module_outputs
  use module_evolution
  use module_energetics,        only : energetic_fluxes
  use module_molecular_cooling, only : molecular_cooling_openfile,err_cool,f_err_cool,n_err_cool
  use module_debug
  use module_tools
  use module_line_excit
  use module_parameters_flags
  use module_parameters_output
  !
  !
  implicit none
  integer (kind=long) :: i
  logical             :: error=.false.    ! status flag
  logical             :: again = .true.   ! continue integration while it is .true.
  logical             :: writ_step        ! if .true. then write to output files
  logical             :: time_up = .false.! transition to j-shock when time_up = .true.
  character (len=80)  :: message_end      ! message written at the end of integration
  !
  ! Integration arguments
  real (kind=dp) :: distance_old = 0.0d0
  real (kind=dp) :: newstep,next_timen,gradn
  real (kind=dp) :: tn_old = 0.0d0
  real (kind=dp) :: vi_old = 0.0d0
  real (kind=dp) :: vn_old = 0.0d0
  real (kind=dp), dimension (:), allocatable :: dyarray,workarray
  real (kind=dp) :: dermax,dermin
  integer :: iflag, jj
  !
  ! DVODE arguments/options
  type (vode_opts) :: options
  integer :: nst,nfe,lenrw,leniw,nni,ncfn,netf,nfea,istate,itask,nje,nlu
  real(kind=dp), dimension(32) :: rstats
  integer,       dimension(32) :: istats
  real(kind=dp), dimension(:), allocatable :: atol
  !
  !
  !
  call write_welcome
  !
  ! Initializations parameters
  !    * shock parameters
  !    * chemical species (+ check)
  !    * chemical reactions (+ check)
  !    * molecule H2 (levels, collision rates, lines)
  !    * molecule SiO (collision rates, Einstein coeff)
  call initialize
  !
  ! Open file for molecular cooling
  call molecular_cooling_openfile
  !
  ! Initialize elemental abundances in (H, C, N, ...)
  call elemental_abundances
  !
  ! Intialize indices and global arrays.
  ! * allocation of the vectors containing physical variables
  ! * dimension = 1:dimtot
  ! * dealocation is made at the end of MHD
  call initialize_indices
  call initialize_arrays
  !
  ! Initialize local arrays
  allocate(workarray(1:dimtot)) ! workin array
  allocate(dyarray  (1:dimtot)) ! derivative of physical variables (y)
  dyarray   = 0.d0
  !
  ! Copy into a working array
  workarray = yarray
  !
  ! Calculate initial values of mass, momentum and energy fluxes
  call energetic_fluxes
  !
  ! Initialize the DVODE arguments
  allocate(atol(dimtot))
  atol           = atol0*nh !*1d-12  ! absolute error
  atol(1:nv_mhd) = 1d-10 !min(1d-12,abs(workarray(1:nv_mhd))*1d-6)
  if (do_h2) atol(bv_h2_lev:ev_h2_lev) = atol0
  itask   = 1
  istate  = 1
  tout    = 1d4
  options = set_opts(dense_j=.true.,abserr=atol0,relerr=rtol)
  !
  ! write informations about the current model (parameters, species,
  ! chemistry, H2, grains) in the file file_info
  if (verbose) then
     call write_info
     call write_info_indices(screen)
  endif
  !
  ! Start of the integration
  if (verbose) write(*,'(/," Start the integration")')
  if (debug) then
     write(*,*) "D-collapse: entering the integration loop:"
     write(*,'(a,/,(10es12.3))') "yarray",yarray
     write(*,'(a,/,(10es12.3))') "H2_lev(1:NH2_lev)%density",H2_lev(1:NH2_lev)%density
  endif
  !  
  do while (again)
     counter=counter+1 ! counts the number of integration steps
     !
     ! save some values before call to DRIVE
     distance_old = distance
     tn_old       = tn
     vn_old       = vn
     vi_old       = vi
     speci(1:nspec)%dens_old         = speci(1:nspec)%density
     h2_lev(1:nh2_lev)%dens_old      = h2_lev(1:nh2_lev)%density
     h2_lines(1:nh2_lines)%emiss_old = h2_lines(1:nh2_lines)%emiss
     !
     call2diffun = 0
     call dvode_f90(diffun,dimtot,workarray,t0_v,tout,itask,istate,options)
     call get_stats(rstats,istats) ! gather the integration statistics for this call
     if (istate<0) goto 10         ! check if an error occurred in vode_f90
     !
     if (debug) write(70,'(10es20.11e3)') tout, &
          tn, ti, te, &
          tr_1, tr_2, tr_3, tr_4, tr_5, &
          tr_6, tr_7, tr_8, tr_9, tr_10, &
          tr_11, tr_12, tr_13, tr_14, tr_15
     !
     ! Step size
     distance  = Tout
     dist_step = distance - distance_old
     call diffun(dimtot,1.d0,yarray,dyarray)
     gradn     = yarray(iv_densityn)/dyarray(iv_densityn)
     newstep   = min((tout_step-1.)*tout,abs(0.1*gradn))
     tout = tout+newstep
     ! call dvindy(tout,1,dyarray,iflag)
     speci(1:nspec)%netrate = dyarray(bv_speci:ev_speci)
     !
     !
     ! flow times for neutral and ions (year)
     timen = timen + 0.5_dp*dist_step*(1./vn_old+1./vn)/yearsec
     timei = timei + 0.5_dp*dist_step*(1./vi_old+1./vi)/yearsec
     next_timen = timen + 0.5_dp*newstep*(1./vn_old+1./vn)/yearsec
     !
     ! Velocity gradient (km.s-1.cm-1) ---
     !--- lower limit : 1 km s-1 / pc
     !   Vgrad = ABS(v_l_der(iv_Vn)) / 1.D5   ! Computed in DIFFUN (evolution.f90)
     vgrad = sqrt(dyarray(iv_vn)*dyarray(iv_vn) + 1.0d10 / (parsec*parsec)) * 1.0d-5
     !
     ! Column density (cm-2) of each chemical species and each H2
     ! energy level (note : trapezian rule)
     speci(1:nspec)%col_dens=speci(1:nspec)%col_dens+&
          0.5_dp*dist_step*(speci(1:nspec)%dens_old+speci(1:nspec)%density)
     h2_lev(1:nh2_lev)%col_dens=h2_lev(1:nh2_lev)%col_dens+&
          0.5_dp*dist_step*(h2_lev(1:nh2_lev)%dens_old+h2_lev(1:nh2_lev)%density)
     !
     ! Integrated intensity (erg/s/cm2/sr) of each H2 line and each
     ! fine structure line.
     H2_lines(1:NH2_lines)%intensity = H2_lines(1:NH2_lines)%intensity + &
          0.5_DP * dist_step / (4._DP*pi) * &
          (H2_lines(1:NH2_lines)%emiss_old + H2_lines(1:NH2_lines)%emiss)
     call line_integ
     !
     ! Mass, momentum and energy fluxes and check the conservation of
     ! those quantities
     call energetic_fluxes
     !
     ! STOP integration if :
     ! (1) Maximal number of steps has been reached
     ! (2) Maximal evolution time has been reached
     ! (3) Maximal distance has been reached
     ! (4) Numerical instabilities prevent conservations of ions
     if (counter.eq.nstep_max) then
        again = .false.
        message_end = "Normal end: maximal number of steps has been reached."
     end if
     if (yarray(iv_densityn).gt.1e10) then
        again = .false.
        message_end = "Normal end: maximal density (1e10) has been reached."
     end if
     if (next_timen.ge.maxtime) then
        if (verbose) write(*,*) timen, next_timen, maxtime
        again = .false.
        message_end = "Normal end: maximal time has been reached."
     end if
     if (distance+dist_step.ge.maxdist) then
        if (verbose) write(*,'(3es12.3)') distance,distance+dist_step,maxdist
        again = .false.
        message_end = "Normal end: maximal distance has been reached."
     end if
     if (abs((densityi-sum(yarray(bv_ion:ev_ion)))/densityi) > 1.0e-2_dp) then
        again = .false.
        message_end = "ions are not conserved"
        write(*,*) densityi,sum(yarray(bv_ion:ev_ion))
     endif
     !
     ! Outputs
     call write_output(again)
     if (verbose) call write_screen
  end do
  !
10 continue
  if (istate<0) then
     write(*,*) 'dvode_f90 returned istate = ', istate
  end if
  !     print the final integration statistics:
  nst = istats(11)
  nfe = istats(12)
  nje = istats(13)
  nlu = istats(19)
  lenrw = istats(17)
  leniw = istats(18)
  nni = istats(20)
  ncfn = istats(21)
  netf = istats(22)
  nfea = nfe
  nfea = nfe - nje
  if (verbose) write (6,90000) lenrw, leniw, nst, nfe, nfea, nje, nlu, nni, ncfn, netf
  !
  ! Write outputs
  ! * H2 excitation diagram
  ! * final abundances in input format
  call write_excit
  call write_species
  if (verbose) write(*,*) message_end
  IF (err_cool) WRITE(*,'("*** WARNING, check the file:", A)') n_err_cool
  !
  ! Terminate execution:
  if (istate==2) then
     if (verbose) write(*,*) 'Successful integration'
  else
     write(*,*) '! ***** The integration was not successful *****'
  end if
  !
  close(f_err_cool)
  deallocate(yarray,dyarray,workarray)
  deallocate(yn)
  deallocate(speci)
  call h2_clean
  !
  !     FORMATS:
90000 FORMAT (' Final statistics for this run:'/'   RWORK size      = ',I8, &
           '   IWORK size =',I8/'   Number of steps = ', &
           I8/'   Number of f-s   = ',I8/'   (excluding J-s) = ', &
           I8/'   Number of J-s   = ',I8/'   Number of LU-s  = ', &
           I8/'   Number of nonlinear iterations           = ', &
           I8/'   Number of nonlinear convergence failures = ', &
           I8/'   Number of error test failures            = ',I8/)
  stop
end program collapse

