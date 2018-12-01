module module_read_fe_data
  use module_precision
  use module_tools, only      : get_file_number
  use module_constants, only  : data_dir
  !___________________________________________________________________
  !
  ! Subroutine called in INITIALIZE to read in data
  ! Data used in LINE_EXCIT, subroutine COLFEPL
  !
  ! Module to read in data required for subroutine FE_LEVEL_POPULATIONS
  ! Reads in atomic data: 
  !                       effective collision strengths (gamfepl)
  !                       radiative decay rates         (aijfepl)
  !                       upper level                   (iupfepl)
  !                       lower level                   (jlofepl)
  !___________________________________________________________________
  !
  implicit none
!!$  include "precision.f90"
  !
  integer (kind=long),             public :: gam_n,aij_n
  integer (kind=long),             private :: i,j, up, down
  real (kind=dp),                  private :: gam, a
  real (kind=dp),dimension(19,19), public :: gamfepl
  real (kind=dp),dimension(110),   public :: aijfepl
  integer, dimension(110),         public :: iupfepl, jlofepl
  !
!!$  ! variables from "precision.f90" are private
!!$  private :: dp, long, minus_infinity, plus_infinity
  !
contains
  !
  subroutine read_fe_data
    use module_parameters_flags
    !-----------------------------------------------------------------
    ! Initialise
    !-----------------------------------------------------------------
    gamfepl  = 0e0_dp
    aijfepl  = 0e0_dp
    iupfepl  = 0e0_dp
    jlofepl  = 0e0_dp
    !
    ! Effective collisions strengths (Fe-electron collisions) for 1000 K
    gam_n=get_file_number()
    open (gam_n, file=data_dir//'/gamma_coll_efe.in', status='old')        
    do i=1,171
       read (gam_n,*) down, up, gam
       gamfepl(down, up) = gam
    end do
    close (gam_n)
    ! Radiative decay rates, upper and lower levels and transition energies
    aij_n=GET_FILE_NUMBER()           
    open (aij_n, file=data_dir//'/aij_fe.in', status='old') 
    do i=1,110
       read (aij_n,*) up, down, A
       aijfepl(i) = A
       iupfepl(i) = up
       jlofepl(i) = down
    end do
    close (aij_n)
    !
    if (verbose) write(*,*) "I-read_fe_data: done"
  end subroutine read_fe_data
end module module_read_fe_data
