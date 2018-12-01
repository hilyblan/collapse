module module_precision
  implicit none
!  integer,       parameter :: dp=selected_real_kind(p=15)   ! real kind
  integer,       parameter :: dp = kind(1.0d0)
  integer,       parameter :: long=kind(1)                  ! integer kind
  real(kind=dp), parameter :: smallest = tiny(1._dp)        ! avoids overflow in diffun
  real(kind=dp), parameter :: maxexp = -log(smallest)       ! minimum exponential argument
end module module_precision

module module_constants
  use module_precision
  implicit none
  !___________________________________________________________________
  !
  ! The module 'MODULE_CONSTANTS' contains mathematical and physical
  ! constants, and unit conversion factors
  ! ___________________________________________________________________
  !
  !
  ! Mathematical constants
  real(kind=dp), parameter :: pi = 3.141592653589793238_dp
  real(kind=dp), parameter :: zero = 0.0_dp
  !
  ! Physical constants
  real(kind=dp), parameter :: mp = 1.67262d-24    ! [g            ] mass of the proton (g)
  real(kind=dp), parameter :: amu= 1.66054d-24    ! [g            ] atomic mass unit
  real(kind=dp), parameter :: me = 9.10939d-28    ! [g            ] mass of the electron (g)
  real(kind=dp), parameter :: qe = 4.80325d-10    ! [esu          ] charge of the electron (esu)
  real(kind=dp), parameter :: alpha_h  = 6.67d-25 ! [cm-3         ] polarisability of h (cm-3)
  real(kind=dp), parameter :: alpha_h2 = 7.70d-25 ! [cm-3         ] polarisability of h2 (cm-3)
  real(kind=dp), parameter :: alpha_he = 2.10d-25 ! [cm-3         ] polarisability of he (cm-3)
  real(kind=dp), parameter :: bohr = 5.29177d-9   ! [cm           ] bohr radius
  real(kind=dp), parameter :: kb = 1.380658d-16   ! [erg K-1      ] boltzmann's constant (erg.k-1)
  ! real(kind=dp), parameter :: r  = 8.314510d7     ! [erg.K-1.mol-1] molar gas constant (erg.k-1.mol-1)
  real(kind=dp), parameter :: gr_cons = 6.6704d-8 ! [erg cm g-2   ] gravitation constant (erg cm g-2)
  real(kind=dp), parameter :: msun    = 1.9885d33 ! [            g] Solar mass
  real(kind=dp), parameter :: everg   = 1.60218d-12 ! 1 ev       =  1.60218d-12 erg
  real(kind=dp), parameter :: yearsec = 3.15569d7   ! 1 year     =  3.15569d7 s
  real(kind=dp), parameter :: kcalev  = 4.3363d-2   ! 1 kcal/mol =  4.3363d-2 ev
  real(kind=dp), parameter :: aucgs   = 6.127d-9    ! conversion between au and cgs
  real(kind=dp), parameter :: parsec  = 3.0857d18   ! 1 pc = 3.0857d18 cm
  real(kind=dp), parameter :: arc_ster= 1/206265d0  ! conversion from arcsec to steradians 
  !
  ! Data directory (collisional cross sections, etc)
  character(len=*), parameter :: data_dir='/home/hilyblan/Codes/Chimies/collapse/data/'
  !
  integer(kind=4),parameter :: nreact_max=40000 ! max. number of reactions
end module module_constants

module module_io
  ! read/write format for chemical reactions
  ! integer(kind=4),  parameter :: name_length = 7  ! (7: old, 8: new) length of specy/element name
  ! character(len=*), parameter  :: format_reaction      = '(2(a7,1x),"= ",4(a7,1x),4es12.3,i4,i6,"!",2x,a)'
  ! character(len=*), parameter  :: format_network_read  = '(2(a7,1x),2x  ,4(a7,1x),3es12.3,a)'
  ! character(len=*), parameter  :: format_network_write = '(2(a7,1x),"= ",4(a7,1x),3es12.3,i4,i6,"!",2x,a)'
  ! character(len=*), parameter  :: format_branching     = '(2(a7,1x),"= ",4(a7,1x),6(2i2,2x),3es12.3,i4,i6,"!",2x,a)'
  integer(kind=4),  parameter  :: name_length = 8  ! (7: old, 8: new) length of specy/element name
  character(len=*), parameter  :: format_reaction      = '(2(a8,1x),"= ",4(a8,1x),4es12.3,i4,i6,"!",2x,a)'
  character(len=*), parameter  :: format_network_read  = '(2(a8,1x),2x  ,4(a8,1x),3es12.3,a)'
  character(len=*), parameter  :: format_network_write = '(2(a8,1x),"= ",4(a8,1x),3es12.3,i4,i6,"!",2x,a)'
  character(len=*), parameter  :: format_branching     = '(2(a8,1x),"= ",4(a8,1x),6(2i2,2x),3es12.3,i4,i6,"!",2x,a)'
end module module_io
