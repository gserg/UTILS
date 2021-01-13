!=======================================================================
!     UTILITY           o r d e r _ c                  
!
!               C O P Y R I G H T -- 2004
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!=======================================================================
!
!     Ordering the c-file  according configuration weights
!
!     Three input arguments:
!     1. initial.c
!     2. final.c 
!     3. eps_c - tolerance for coefficients
!
!     Default values:  cfg.out cfg.inp 0.d0
!     
!=======================================================================

      Implicit real(8) (A-H,O-Z)

      Character(80) :: AS
      Character(64), Allocatable :: AC(:)
      Character(60), Allocatable :: BC(:)

      Real(8), Allocatable :: WT(:)
      Integer, Allocatable :: IP(:)
 
      Integer :: nuc = 1   ! input.c
      Integer :: out = 2   ! result.c

      Character(40) :: AF = 'cfg.out'
      Character(40) :: BF = 'cfg.inp'

      Real(8) :: eps_c = 0.d0

!----------------------------------------------------------------------

      iarg = command_argument_count()
      if(iarg.ge.1)  Call get_command_argument(1,AF)
      if(iarg.ge.2)  Call get_command_argument(2,BF)
      if(iarg.ge.3)  then
       Call get_command_argument(3,AS);  read(AS,*) eps_c
      end if

      if(AF.eq.'?') then
       write(*,'(/a)') 'order_c   orders c-file according to configuration weights'   
       write(*,'(/a)') 'Call as:  order_c input.c output.c eps_c'
       write(*,'(/a)') 'Default:          cfg.out cfg.inp  0.d0'                  
       write(*,'(/a)') 'eps_c  -  tolerance for the weights'                  
       Stop ' '
      end if

      Open(nuc,file=AF,STATUS='OLD')
      ncfg=Idef_ncfg(nuc)
      Allocate(AC(ncfg),BC(ncfg),WT(ncfg),IP(ncfg))
      Open(out,file=BF)

! ... read configuration:
 
      rewind(nuc)
      read(nuc,'(a)') AS;   write(out,'(a)') AS
      read(nuc,'(a)') AS;   write(out,'(a)') AS
 
      i=0
    1 read(nuc,'(A)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
       i=i+1; AC(i)=AS(1:64)
       read(AS(65:75),'(f11.8)') WT(i)
       read(nuc,'(a)') BC(i)
      go to 1
    2 Continue
      if(i.ne.ncfg) Stop ' i <> ncfg'
 
! ... order the configurations:

      Do i=1,ncfg; IP(i)=i; End do
 
      Do i=1,ncfg-1
       m=i
       Do j=i+1,ncfg
        if(abs(WT(IP(j))).gt.abs(WT(IP(m)))) m=j
       End do
       if(m.ne.i) then
        k=IP(m); IP(m)=IP(i); IP(i)=k
       end if
      End do

! ... output:
 
      k=0
      Do m=1,ncfg
       i=IP(m); if(abs(WT(i)).lt.Eps_c) Exit
       write(out,'(a64,1x,f10.8)') AC(i),WT(i)
       AS=BC(i); ii=len_trim(AS)
       write(out,'(a)') AS(1:ii)
       k=k+1
      End do
      write(out,'(a1)') '*'
 
      Close(nuc); Close(out)
 
      End  ! interface ORDER_C


!======================================================================
      Integer Function Idef_ncfg(nu)
!======================================================================
!
!     gives the number of configuration in c-file (unit nu)
!
!----------------------------------------------------------------------

      IMPLICIT NONE
      
      INTEGER, INTENT(in) :: nu

      INTEGER :: ncfg
      Character(5) :: AS

      rewind(nu)
      ncfg=0
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(nu,'(a)') AS
      ncfg=ncfg+1
      go to 1
    2 rewind(nu)
      Idef_ncfg=ncfg

      End Function Idef_ncfg

