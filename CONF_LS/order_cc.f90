!=======================================================================
!     UTILITY           o r d e r _ c c                 
!
!               C O P Y R I G H T -- 2004
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!=======================================================================
!
!     Ordering the c-file  according to configuration weights
!
!     name.c  ->> name.c (reordered) + name.conf (with configurations only)
!                                                (no terms)
!
!     used to get the most important configurations
!
!=======================================================================

      Implicit real(8) (A-H,O-Z)

      Character(80) :: AS,BS
      Character(64), Allocatable, Dimension(:) :: AC
      Character(60), Allocatable, Dimension(:) :: BC

      Real(8), Allocatable, Dimension(:) :: WT,WTM
      Integer, Allocatable, Dimension(:) :: IP,IW
 
      Integer :: nua = 1   ! name.c
      Integer :: nub = 2   ! name.conf

      Character(40) :: AF = 'name.c'
      Character(40) :: BF = 'name.conf'

!----------------------------------------------------------------------

      iarg = command_argument_count()
      if(iarg.gt.0) Call get_command_argument(1,AF)
      if(iarg.lt.1.or.AF.eq.'?')  then
       write(*,'(/a)') 'order_cc  orders c-file according to configuration weights'   
       write(*,'(/a)') 'Call as:  order_cc  name.c '
       write(*,'(/a)') 'Action:   name.c  --> name.c (reordered) + name.conf (no terms)'                
       write(*,'(/a)') 'used to get the most important configurations'            
       Stop ' '
      end if

      i = LEN_TRIM(AF)
      if(AF(i-1:i).ne.'.c') AF(i+1:i+2)='.c'
      Open(nua,file=AF,STATUS='OLD')
      ncfg=Idef_ncfg(nua)
      Allocate(AC(ncfg),BC(ncfg),WT(ncfg),WTM(ncfg),IP(ncfg),IW(ncfg))

      i = INDEX(AF,'.',BACK=.TRUE.); BF=AF(1:i)//'conf'
      Open(nub,file=BF)

! ... read configuration:
 
      rewind(nua)
      read(nua,'(a)') AS;   write(nub,'(a)') AS
      read(nua,'(a)') BS;   write(nub,'(a)') BS
 
      i=0
    1 read(nua,'(A)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
       i=i+1; AC(i)=AS(1:64)
       read(AS(65:75),'(f11.8)') WT(i)
       read(nua,'(a)') BC(i)
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
 
      rewind(nua); read(nua,*); read(nua,*)

      Do m=1,ncfg; i=IP(m)
       write(nua,'(a64,1x,f10.8)') AC(i),WT(i)
       write(nua,'(a)') TRIM(BC(i))
      End do
      write(nua,'(a1)') '*'
 
      Close(nua)

! ... output configurations: 

      WT = WT * WT; IW = 1; wtm = WT

      Do i = 1,ncfg-1
       if(WT(i).eq.0.d0) Cycle
       Do j = i+1,ncfg
        if(WT(j).eq.0.d0) Cycle
        if(AC(j).eq.AC(i)) then
         if(WT(j).gt.WTM(i)) wtm(i) = WT(j) 
         WT(i) = WT(i) + WT(j); IW(i) = IW(i) + 1
         WT(j) = 0.0
        end if
       End do
      End do

      WT = sqrt(WT); wtm = sqrt(wtm)

      Call SORTr(ncfg,WT,IP)

      i = INDEX(AF,'.'); BF=AF(1:i)//'conf'
      Open(nub,file=BF)

      Do i = ncfg,1,-1;  ic=IP(i)
       if(WT(ic).eq.0.d0) Cycle
       write(nub,'(a64,3f8.4)') AC(ic), WT(ic),WTM(ic)
      End do
      write(nub,'(a)') '*'

      write(nub,'(/a/)') 'Above weights show:'
      write(nub,'(a)') '1. Total contribution of given configuration'
      write(nub,'(a)') '2. Maximum coefficient for given configuration'

      Close(nub)
 
      End  ! interface ORDER_CC




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


!====================================================================
      Subroutine SORTR(n,S,IPT)
!--------------------------------------------------------------------
!
!     gives sorting pointer IPT for real(8) array S
!
!--------------------------------------------------------------------

      IMPLICIT none

      Integer(4), intent(in) :: n
      Real(8), intent(in), dimension(n) :: S 
      Integer(4), intent(out), dimension(n) :: IPT 

      Integer(4) :: i,i1,j1,i2,j2

      Do i=1,n; IPT(i)=i;  End do

      Do i1=1,n-1; j1=IPT(i1)
       Do i2=i1+1,n;  j2=IPT(i2)
        if(S(j1).gt.S(j2)) then
         IPT(i2)=j1; j1=j2
        end if
       End do
       IPT(i1)=j1
      End do

      End Subroutine SORTR
