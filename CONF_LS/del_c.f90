!======================================================================
!     utility       D E L _ C 
!
!                   C O P Y R I G H T -- 2011
!
!     Written by:   Oleg Zatsarinny
!======================================================================
!
!     Deleting the configurations in c-file
!
!     Four input arguments:
!
!     1. initial.c     (cfg.out)
!     2. final.c       (cfg.inp)
!     3. del_list      (one conf. index for a line, ending with * )
!     4. flag          indicate list beginning, optional
!
!     Call as:  del_c  1.c 2.c del_list [>>>]
!
!     Default:  del_c  cfg.out cfg.inp del_list
!
!     Additional:  log-file after ci_bnk may be used as del_list !!!
!
!=======================================================================

      Implicit real(8) (A-H,O-Z)

      Character(80) :: AS=' ',BS

      Integer :: nuc = 1;  Character(40) :: AF = 'cfg.out'   ! input.c
      Integer :: out = 2;  Character(40) :: BF = 'cfg.inp'   ! result.c
      Integer :: nud = 3;  Character(40) :: DF = 'del_list'  ! del_list

      Integer, Allocatable :: idel(:)

      Call inf_del_c
!----------------------------------------------------------------------

      iarg = IARGC()
      if(iarg.ge.1)  Call GETARG(1,AF)
      if(iarg.ge.2)  Call GETARG(2,BF)
      if(iarg.ge.3)  Call GETARG(3,DF)
      if(iarg.ge.4)  Call GETARG(4,AS)

      Open(nuc,file=AF,STATUS='OLD')
      Open(out,file=BF)
      Open(nud,file=DF,status='OLD')

      nc = Idef_ncfg(nuc)
      Allocate(idel(nc))
      idel = 0

      i = Ifind_position(nud,'Non-trivial total overlaps')
      if(i.gt.0) read(nud,*)
write(*,*) i

      if(LEN_TRIM(AS).gt.0) i = Ifind_position(nud,trim(DF))

   10 read(nud,'(a)',end=20) AS
      if(LEN_TRIM(AS).eq.0) go to 10
      if(AS(1:1).eq.'*') go to 20
write(*,*) trim(AS)
      read(AS,*) ii
write(*,*) ii
      if(ii.le.0.or.ii.gt.nc) go to 10
      idel(ii) = 1
      go to 10
   20 Continue

! ... read configuration:
 
      rewind(nuc)
      read(nuc,'(a)') AS;   write(out,'(a)') AS
      read(nuc,'(a)') AS;   write(out,'(a)') AS
 
      i=0
    1 read(nuc,'(A)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      i=i+1
      read(nuc,'(a)') BS
      if(idel(i).ne.0) go to 1
      write(out,'(a)') trim(AS)     
      write(out,'(a)') trim(BS)
      go to 1
    2 write(out,'(a)') '*'
 
      Close(nuc); Close(out)
 
      End  ! interface DEL_C


!======================================================================
      Integer Function Idef_ncfg(nu)
!======================================================================
!
!     gives the number of configuration in c-file (unit nu)
!
!----------------------------------------------------------------------

      IMPLICIT NONE
      
      INTEGER(4), INTENT(in) :: nu

      INTEGER(4) :: ncfg
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


!======================================================================
      Integer Function Ifind_position(nu,name)
!======================================================================
!     find position of line with "name" in the begining
!----------------------------------------------------------------------
      
      Implicit none

      Integer, Intent(in) :: nu
      Character(*), Intent(in) :: name

      Character(80) :: AS
      Integer :: i,j
 
      Ifind_position = 0

      i=LEN_TRIM(name); j=0
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      j=j+1
      if(AS(1:i).ne.name) go to 1
      Ifind_position = j
      Backspace(nu)
      Return
    2 rewind(nu)

      End Function Ifind_position


!======================================================================
      Subroutine inf_del_c
!======================================================================
!     provide screen information about DEL_C utility
!----------------------------------------------------------------------
       
      Implicit real(8) (a-h,o-z)

      Character(80) :: A

      iarg = IARGC()
      if(iarg.gt.0) Call GETARG(1,A)
      if(iarg.eq.0.or.A.ne.'?')  Return

      write(*,*) &
'                                                                 ',&
'del_c deletes the configurations in c-file                       ',&
'                                                                 ',&
'Four input arguments:                                            ',&
'                                                                 ',&
'1. initial.c                                                     ',&
'2. final.c                                                       ',&
'3. del_list      (one conf. index for a line, ending with * )    ',&
'4. flag          indicate beginning of the list, optional        ',&
'                                                                 ',&
'Call as:  del_c  1.c 2.c del_list [???]                          ',&
'                                                                 ',&
'Default:  del_c  cfg.out cfg.inp del_list                        ',&
'                                                                 ',&
'Additional feature:  log-file after ci_bnk may be used as del_list !',&
'                                                                  ',&
'                                                                 '
      Stop ' '

      End Subroutine inf_del_c

