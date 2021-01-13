!======================================================================
!     UTILITY   L F I L E
!
!                   C O P Y R I G H T -- 2004
!
!     Written by:   Oleg Zatsarinny
!
!----------------------------------------------------------------------
!     find the average contributions for given configurations
!     (c-file) from the subset of solutions in given l-file
!
!      case.c + case.l -->  case.ca + case.cb
!
!     Example:
!
!     lfile 1.c 5   - analysing first 5 solutions in 1.l 
!
!                     and output in 1.ca all ACS in order of importance
!                     with c = sqrt[(c1^2 +c2^2 +...)/ nsol]
!
!                     and output in 1.cb only configurations (no terms!)
!                     with c = sqrt[(c1^2 +c2^2 +...)/nsol] 
!                     summed over all terms
!
!     This utility is used to find most important configurations
!     when we want to reduce the expansions to minimum
!
!----------------------------------------------------------------------

      Implicit double precision (A-H,O-Z)

      Logical :: EX

      Character(80) :: AS
      Character(40) :: AF,BF

      Character(64) :: CONFIG, LABEL
      Character(60) :: COUPLE
      Character(64), Allocatable :: CONF(:)
      Character(60), Allocatable :: COUP(:)
      Real(8), Allocatable :: C(:),CC(:)
      Integer, Allocatable :: IP(:)
      
      Integer :: nuc=1  !   c-file with CAS       
      Integer :: nul=2  !   l-files with weights
      Integer :: nua=3  !   ca-file - average contribution for CAS              
      Integer :: nub=3  !   cb-file - average contribution for configuration             

      Real(8) :: eps_c=0.d0
      
      Call inf_lfile

!----------------------------------------------------------------------
! ... input data:

      iarg = command_argument_count()
      Call get_command_argument(1,AF)        
      Call get_command_argument(2,BF); read(BF,*) nsol
      if(iarg.gt.2) then
       Call get_command_argument(3,BF); read(BF,*) eps_c
      end if
 
! ... input configurations:

      i=INDEX(AF,'.',BACK=.TRUE.); BF=AF(1:i)//'c'
      Open(nuc,file=BF,status='OLD')
      nc = Idef_ncfg(nuc)
      write(*,*) trim(AF),'  nc=',nc; if(nc.eq.0) Stop 'nc=0'
      Allocate(CONF(nc),COUP(nc),IP(nc),C(nc),CC(nc))
      C = 0.d0

      rewind(nuc); ic=0
      Do
       read(nuc,'(a)') AS
       if(AS(1:1).eq.'*') Exit
       if(AS(5:5).ne.'(') Cycle
       read(AS,'(a64,f11.8)') CONFIG,S
       read(nuc,'(a60)') COUPLE
       ic=ic+1; CONF(ic)=CONFIG; COUP(ic)=COUPLE
      End do
      if(ic.ne.nc) Stop ' ic <> nc '

! ... proceed solutions:

      C = 0.d0

      i=INDEX(AF,'.',BACK=.TRUE.); BF=AF(1:i)//'l'
      Open(nul,file=BF,status='OLD')
      i=Idef_st(nul); if(i.lt.nsol) nsol=i

      Do isol = 1,nsol
       CC = 0.d0
       Call Read_sol('l',nul,nc,CC,Label,E,jot)
       CC = abs(CC)
       Do i=1,nc; if(CC(i).gt.C(i)) C(i) = CC(i); End do
      End do
      
! ... write the configurations in ordered form:

      Call SORTR(nc,C,IP)

      i=INDEX(AF,'.'); BF=AF(1:i)//'ca'; Open(nua,file=BF)
      rewind(nuc)
      read(nuc,'(a)') AS; write(nua,'(a)') AS
      read(nuc,'(a)') AS; write(nua,'(a)') AS
      Do i = nc,1,-1;  ic=IP(i)
       if(CC(ic).lt.eps_c) Exit
       write(nua,'(a64,f11.8)') CONF(ic),C(ic)
       write(nua,'(a60)') COUP(ic)
      End do
      write(nua,'(a)') '*'

! ... max. configurations:

      i=INDEX(AF,'.'); BF=AF(1:i)//'cb'; Open(nub,file=BF)

      Do i = 2,nc
       Do j = 1,i-1
        if(C(j).eq.0.d0) Cycle
        if(CONF(j).eq.CONF(i)) then
         if(C(i).gt.C(j))  C(j)=C(i)
         C(i) = 0.0
        end if
       End do
      End do

      Call SORTR(nc,C,IP)

      rewind(nuc)
      read(nuc,'(a)') AS; write(nub,'(a)') AS
      read(nuc,'(a)') AS; write(nub,'(a)') AS
      Do i = nc,1,-1;  ic=IP(i)
       if(C(ic).eq.0.d0) Cycle
       write(nub,'(a64,f11.8)') CONF(ic),C(ic)
      End do
      write(nub,'(a)') '*'

      End !  Utility llfile



!======================================================================
      Subroutine inf_lfile
!======================================================================
!     provide screen information about utility  'lfile'
!----------------------------------------------------------------------
       
      Implicit real(8) (a-h,o-z)

      Character(80) :: A=' '

      iarg = command_argument_count()
      if(iarg.gt.0) Call get_command_argument(1,A)
      if(A.ne.'?')  Return

      write(*,'(a)') &
'                                                                  ',&
'lfile finds the average contributions for given configurations    ',&
'(c-file) from the subset of solutions in given l-file             ',&
'                                                                  ',&
' name.c + name.l -->  name.ca + name.cb                           ',&
'                                                                  ',&
'Example:                                                          ',&
'                                                                  ',&
'lile 1.l 5 [0.01] - analysing first 5 solutions in 1.l            ',&
'                                                                  ',&
'name.ca  -  all CAS in order of importance                        ',&
'            with c = sqrt[(c1^2 +c2^2 +...)/ nsol]                ',&
'            (CAS with weigts < 0.01 (eps_c) are ignored)          ',&
'                                                                  ',&
'name.cb  -  all configurations in order of importance             ',&
'            with c = sqrt[(c1^2 +c2^2 +...)/nsol]                 ',&
'            summed over all terms                                 ',&
'                                                                  ',&
'This utility is used to find most important CAS or configurations ',&
'when we want to reduce the expansions to minimum                  ',&
'                                                                  ',&
'                                                                  '
      Stop ' '

      End Subroutine inf_lfile


