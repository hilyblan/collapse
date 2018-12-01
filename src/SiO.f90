MODULE MODULE_SiO
  use module_precision
  !_________________________________________________________________
  !
  ! Read SiO data and initialize SiO-related variables
  !_________________________________________________________________
  !
  implicit none
!!$  include "precision.f90"
  !
  ! physical constants
  integer(kind=4), parameter :: jmin_sio=0, jmax_sio=20 ! min. and max. values of j
  real(kind=dp),   parameter :: brot_sio=21787.453_dp   ! rotational constant (mhz)
  real(kind=dp),   parameter :: moment_sio=3.0982_dp    ! dipolar moment (debye)
  !
  ! Einstein coefficients
  ! calculated in EINSTEIN_COEFF_SiO
  real(kind=dp), dimension(jmax_sio) :: a_sio  ! spontaneous emission (cgs)
  real(kind=dp), dimension(jmax_sio) :: bs_sio ! stimulated emission (cgs)
  real(kind=dp), dimension(jmax_sio) :: ba_sio ! absorption (cgs)
  !
  ! SiO-H2 collision rates
  ! variables read in READ_SiO_RATES
  ! * rates (cm3/s) : Rate_SiO(i,j,k) -> i=initial state, j=final state, k=temperature
  ! * temperatures (K) which define Rate_SiO(i,j,k)
  real(kind=dp),dimension(jmin_sio:jmax_sio,jmin_sio:jmax_sio,8) :: rate_sio
  real(kind=dp), dimension(8) :: temp_sio
  !
!!$  ! variables from "precision.f90" are private
!!$  private :: dp, long, minus_infinity, plus_infinity
  !
contains
  !
  subroutine read_sio_rates
    use module_tools, only : get_file_number
    use module_constants, only : zero,data_dir
    use module_parameters_flags
    !---------------------------------------------------------------------------
    !
    ! Initialize collision rates for SiO-H2.
    !
    !---------------------------------------------------------------------------
    implicit none
    integer(kind=long) :: i, ji, jf, error, nrates
    character(len=*), parameter :: name_file_sio=data_dir//'/coeff_sio.in'
    character(len=*), parameter :: format_sio_temp='(8x,8f9.1)'
    integer                     :: file_sio
    !
    ! initialization
    rate_sio=zero
    temp_sio=zero
    !
    ! file opening
    file_sio = get_file_number()
    open(file_sio,file=name_file_sio,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    ! comments
    do i=1,8
       read(file_sio,*)
    enddo
    ! temperatures
    read(file_sio,format_sio_temp)temp_sio
    ! comments
    do i=1,3
       read(file_sio,*)
    enddo
    ! collision rates
    
    nrates=(size(rate_sio,dim=1)-1)*size(rate_sio,dim=2) 
    error=0
    do i=1, nrates
       read(file_sio,*,iostat=error)ji,jf,rate_sio(ji,jf,:)
       if (error>0) stop "E-read_sio_rates: error in read_sio_rates"
    enddo
    !
    if (verbose) write(*,*) 'I-read_sio_rates: done'
  end subroutine read_sio_rates


  subroutine einstein_coeff_sio
    !---------------------------------------------------------------------------
    ! Computes Einstein coefficients for the SiO molecule
    ! method=simple rotator molecule
    ! ---------------------------------------------------------------------------
    use module_constants, only : zero
    implicit none
    integer(kind=4) :: j
    real(kind=dp)   :: jj, gl, gu

    ! initializations
    a_sio=zero ; bs_sio=zero ; ba_sio=zero

    ! computation
    do j=1,jmax_sio
       jj=dble(j) ! conversion to real
       gl        = 2._dp*(jj-1._dp)+1._dp
       gu        = 2._dp*jj+1._dp
       bs_sio(j) = 7.896d8 * moment_sio**2 * jj/gu ! stimulated emission
       ba_sio(j) = bs_sio(j) * gu/gl               ! absorption
       a_sio(j)  = 11.7944d-29 * brot_sio**3 *jj**3 * bs_sio(j) ! spontaneous emission
    end do
  end subroutine einstein_coeff_sio
  !
end module module_sio
