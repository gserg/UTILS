!======================================================================
!  zf_cc - utility for the f-value calculations between set of c-files
!======================================================================
!
!     INPUT:    zf_cc.inp
!
!     OUTPUT:   zf_res from DBSR_DMAT program
!
!     SYSTEM CALL:  DBSR_MULT, DBSR_DMAT
!
!     Notes:  1. you would better delete all mult_bnk before running
!             2. results are appended to zf_res
!
!     zf_cc.inp containes:
!
!     atype = E1           or  E2, M1 and so on
!     nfiles = ...
!     list of c-files with energies and exp.coefficients with one J
!
!---------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)

      Integer :: inp=5; Character(20) :: AF_inp = 'zf_cc.inp'
      Integer :: nuc=1

      Character(1) :: blank = ' '
      Character(2) :: atype = 'E1'
      Character(4) :: cc = ' c c'
      Character(80) :: AS,AI,AJ
      Character(80), allocatable :: files(:) 
      Integer, allocatable :: jot(:), parity(:)

!----------------------------------------------------------------------
! ... check input data:

      Call Read_name(AS)

      if(AS.eq.'?') then
        write(*,'(a/)') 
        write(*,'(a/)') 'zf_cc - utility for the f-value calculations between set of c-files'
        write(*,'(a/)') 'INPUT:    zf_cc.inp, {name.c and name.bsw} '
        write(*,'(a/)') 'OUTPUT:   appended f-falues in zf_res file' 
        write(*,'(a/)') 'SYSTEM CALL:  DBSR_MULT, DBSR_DMAT'
        write(*,'(a/)') 'zf_cc.inp containes:'
        write(*,'(a )') 'atype = E1           or  E2, M1 and so on'
        write(*,'(a )') 'nfiles = ...'
        write(*,'(a/)') 'list of c-files with energies and exp.coefficients with one J'
        Stop ' '
       end if

      Call Check_file(AF_inp)
      open(inp,file=AF_inp)

      Call Read_apar(inp,'atype',atype)
      Read(atype,'(1x,i1)') kpol

      Call Read_ipar(inp,'nfiles',nfiles)
      if(nfiles.le.0) Stop 'nfiles = 0'

      Allocate(files(nfiles),jot(nfiles),parity(nfiles))

      Do i=1,nfiles
       read(inp,*) files(i)
       Call Check_file(files(i))
       open(nuc,file=files(i))
       Call Jdef_JP(nuc,jot(i),parity(i)) 
       close(nuc)
      End do

!----------------------------------------------------------------------

      Do i=1,nfiles; AI=files(i); ii=LEN_TRIM(AI) 
      Do j=i,nfiles; AJ=files(j); jj=LEN_TRIM(AJ)

       if(atype(1:1).eq.'E'.and.mod(kpol,2).eq.1.and. &
          parity(i).eq.parity(j)) Cycle
       if(atype(1:1).eq.'M'.and.mod(kpol,2).eq.0.and. &
          parity(i).eq.parity(j)) Cycle

       if(ITRA(jot(i),kpol+kpol,jot(j)).eq.0) Cycle
       
       AS = 'dbsr_mult '//AI(1:ii)//blank//AJ(1:jj)//blank//atype
       Call System(AS)

       AS = 'dbsr_dmat '//AI(1:ii)//blank//AJ(1:jj)//cc
       Call System(AS)

      End do
      End do 
       
      End ! program zf_cc

