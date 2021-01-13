!=====================================================================
      Module bsr_se     
!=====================================================================
!     contains main variable and arrays for program bsr_se
!---------------------------------------------------------------------
      Implicit none

! ... files:

      Integer, parameter :: ma = 80
      Character(ma) :: AF

      Integer :: inp =  5; Character(ma) :: AF_inp  = 'bsr_se2.inp'
      Integer :: pri =  6; Character(ma) :: AF_log  = 'bsr_se2.log'
      Integer :: out =  7; Character(ma) :: AF_out  = 'projection2'
      
      Integer :: nut = 10; Character(ma) :: AF_tar  = 'target'
                           Character(ma) :: AF_ps   = 'target_ps'

      Integer :: nui = 11; Character(ma) :: AF_mat  = 'bsr_mat.nnn'
      Integer :: nur = 12; Character(ma) :: AF_rsol = 'rsol.nnn'
      Integer :: nub = 13; Character(ma) :: AF_bnd  = 'ubound.nnn'
      Integer :: nuo = 14; Character(ma) :: AF_rovl = 'rovl.nnn'
      Integer :: nuh = 15; Character(ma) :: AF_h    = 'h.nnn'

! ... main dimension parameters:

      Integer :: klsp, kch, khm, nopen

! ... pseudo-states parameters:

      Integer :: nps
      Integer, allocatable :: lps(:), sps(:), pps(:)
      Real(8), allocatable :: Eps(:)

! ... ASYPCK parameters:

      Integer :: debug = 0        ! debug level = 0,1,2
      Integer :: iauto = 1        ! automatic choice of method
      Integer :: mfgi  = 300      ! number of point in outer region
      Real(8) :: AC    = 0.01     ! accuracy, < 0.1, but < 0.0001
      Real(8) :: DELTA = 0.1      ! r1 - r2, < 0.2 

! ... scattering arrays:
 
      Real(8), allocatable :: RMAT(:,:), KMAT(:,:), ECH(:), &
                              F(:,:), G(:,:), FP(:,:), GP(:,:)
! ... overlaps arrays:

      Real(8), allocatable :: a(:,:),v(:),d(:)
      Real(8), allocatable :: dr(:),di(:) 
      Real(8), allocatable :: ddr(:,:),ddi(:,:),wt(:,:) 
      Integer, allocatable :: nop(:) 

! ... other parameters:

      Integer :: itarg = 1  ! debug level = 0,1,2
      Real(8) :: EL =-1.d0  ! input energy of ejected electron (in eV)
      Real(8) :: qqq        ! in respect to the itarg-ion state

      Real(8) :: EK = 0.d0  ! incident etargy in Ry.  
      Real(8) :: E0         ! incident etargy in eV.  

      Real(8) :: EL2= 0.d0  ! energy of ejected electron (in eV)
      Real(8) :: qq2        ! k^2 for ejected electron
                            ! in respect to the ground-ion state

      Real(8) :: EL1= 0.d0  ! energy of scattering electron (in eV)
      Real(8) :: qq1        ! k^2 for scattering electron
                            ! in respect to the ground-ion state

      Real(8) :: Eion_Ry    ! ionization energy (in Ry)
      Real(8) :: Eion_eV    ! ionization energy (in Ry)

      Real(8) :: E1         ! absolute energy (in a.u.) for ion ground state
      Real(8) :: EE         ! absolute energy (in a.u.) for e-A+ system
      Real(8) :: ET         ! absolute energy (in a.u.) for e-A  system

      Real(8) :: eps_e = 1.d-5

      Real(8) :: Ry = 13.6056d0

      End module bsr_se


