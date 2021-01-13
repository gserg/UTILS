!=====================================================================
!     UTILITY  sec_cascb 
!=====================================================================
!
!     sec_cascb provides cross sections with cascade contribution
!
!     scattering information is given in file 'target'
!     
!     branching ratios are taken from file 'f_values'
!                                                                                               
!     direct cross sections are taken from 'tr_###_###'
!                                  or from 'cr_###_###'
!
!     Arguments: itr1,itr2, jtr1,jtr2 - range for transition indexes
!                                       to be considered
!                                       (itr -initial index, jtr - final) 
!
!                cr - 'tr' or 'cr'  - suppose tr_###_### or cr_###_### files
!                                     with direct cross sections 
!
!                steps   = 1 -> direct cascade 
!                          2 -> through additional intermediate level
!                          5 -> maximum intermediate levels
!                                    
!     Call as:   sec_cascb  itr1=.. itr2=.. jtr1=.. jtr2=.. cr=.. steps=..          
!
!     Results:   ca_###_###, where ### stand for transition indexes 
!
!     The folder should contain all needed tr- or cr-files, e.g., 
!     if we need cascade for transition from ground state, 
!     all  tr_001_###|cr_###_### should be presented
!---------------------------------------------------------------------
      Use target; Use channels

      Implicit real(8) (A-H,O-Z)

      Integer :: i11=1, i12=1, i21=2, i22=2, max_level=1

      Integer, allocatable :: ilv1(:),ilv2(:),itr(:)
      Real(8), allocatable :: STR(:), SS(:), E_eV(:), E_ry(:), sec(:), sec0(:)

      Logical EX

      Character(3)  :: AT
      Character(2)  :: AF, cr ='cr'

      Integer :: nut=1;  Character(20) :: AF_t  = 'target'
      Integer :: nub=1;  Character(20) :: AF_b  = 'f_values'
      Integer :: nui=2;  Character(20) :: AF_tr = 'cr_###_###'
      Integer :: out=3;  Character(20) :: AF_ca = 'ca_###_###'

      Real(8) :: Ry=13.60569172

      Real(8) :: Eps_E = 1.d-8

      Call Inf_sub
!----------------------------------------------------------------------
! ... target information:

      Inquire(file=AF_t,EXIST=EX)
      if(.not.EX) then
       Stop  ' No target file: run H_targb ! '
      else
       Open(nut,file=AF_t); Call R_target(nut); close(nut)
      end if

      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion
      E1 = etarg(1); etarg = (etarg-E1)*2.0
      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

      ntr=ntarg*(ntarg+1)/2;  Allocate(itr(10*ntr),SS(10*ntr))

!----------------------------------------------------------------------
! ... read branching ratios:

      Inquire(file=AF_b,EXIST=EX)
      if(.not.EX) then
       Stop  ' No branching ratios (f_values) ! '
      else
       Open(nub,file=AF_b,status='OLD')
      end if

      Call Read_ipar(nub,'nt',ntrans)
      write(*,*) ' # of branching ratios = ',ntrans

      Allocate(ilv1(ntrans), ilv2(ntrans), STR(ntrans))

      Do it = 1,ntrans
       read(nub,*) ilv1(it), ilv2(it), S,F,A,STR(it)
      End do
      Close(nub)

!----------------------------------------------------------------------
! ... read arguments from command line:

      i11=1; i12=ntarg-1;  i21=2; i22=ntarg

      Call Read_iarg('itr1',i11)
      Call Read_iarg('itr2',i12)

      Call Read_iarg('jtr1',i21)
      Call Read_iarg('jtr2',i22)

      max_level = 1;  cr = 'cr'
      Call Read_iarg('steps',max_level)
      Call Read_aarg('cr',cr)

!----------------------------------------------------------------------
! ... cycle over transitions:

      Do ktr1=i11,i12; Do ktr2=i21,i22

       if(ktr1.gt.ktr2) Cycle

! ... find cross section for main transition:

      write(AF_tr,'(a,a,i3.3,a,i3.3)') cr,'_',ktr1,'_',ktr2

      Inquire(file=AF_tr,EXIST=EX)
      if(.not.EX) then
       write(*,*) 'file ',AF_tr,' missed !'; Cycle
      else
       Open(nui,file=AF_tr,status='OLD')
      end if

! ... define number of energies, ne:

      ne=0
      rewind(nui); read(nui,*)
    1 read(nui,*,end=2) x;  ne = ne + 1; go to 1
    2 Continue

      Allocate(E_ev(ne),E_ry(ne),sec(ne),sec0(ne)) 

      rewind(nui); read(nui,*)
      Do ie = 1,ne
       read(nui,*) E_ev(ie),sec(ie),E_ry(ie)
      End do
      sec0 = sec

!----------------------------------------------------------------------
! ... find cascad contribution:

      SS(1)=1.d0;  itr(1)=ktr2;  ii=1; ilevel=1;  i1=1; i2=1

   10 Continue

      Do i = i1,i2; jtr1=itr(i)      ! cycle over new lower levels

      Do it=1,ntrans

       if(ilv1(it).ne.jtr1) Cycle; jtr2=ilv2(it);      if(jtr2.gt.31) Cycle     ! ???

       ii=ii+1;  itr(ii)=jtr2; SS(ii)=SS(i)*STR(it)

! ... add contribution: 

      write(AF_tr,'(a,a,i3.3,a,i3.3)') cr,'_',ktr1,'_',jtr2

      Inquire(file=AF_tr,EXIST=EX)
      if(.not.EX) then
       write(*,*) 'file ',AF_tr,' missed !'; Cycle
      else
       Open(nui,file=AF_tr,status='OLD')
      end if

      rewind(nui); read(nui,*)
    3 read(nui,*,end=4) e1,s,e2;  
      if(s.le.0.d0) go to 3

      met = 0
      Do ie = 1,ne
       if(abs(e2-E_Ry(ie)).gt.Eps_E) Cycle
       sec(ie) = sec(ie) + s * SS(ii)       
       met = 1
       Exit
      End do

      if(s.ne.0.d0.and.met.eq.0) then
       write(*,'(a,2f10.6,a)') trim(AF_tr),e1,e2,'  no match for this energy'
       Stop
      end if

      go to 3
    4 Continue

      End do  !  over it - transitions

      End do  !  over i -> states of this level

      ilevel=ilevel+1;  i1=i2+1;  i2=ii
 
      if(ilevel.le.max_level.and.i2.gt.i1) go to 10

      write(*,*) 'transition: ',ktr1,'  ->',ktr2,'  ilevel =',ilevel-1

!----------------------------------------------------------------------
! ... output data:

      write(AF_ca,'(a,i3.3,a,i3.3)') 'ca_',ktr1,'_',ktr2

      Open(out,file=AF_ca)

      write(out,'(6x,a,7x,a,9x,a,12x,a)') 'E(eV)',  'cascade', 'direct', 'E(Ry)'

      Do ie=1,ne
       write(out,'(f14.8,2e16.8,f12.6)') E_Ry(ie)*Ry,sec(ie),sec0(ie),E_Ry(ie)   
      End do

      Close(out)

      Deallocate(E_eV,E_Ry,sec,sec0)

!----------------------------------------------------------------------

      End do; End do  ! over transitions


      End  ! UTILITY  sec_cascb



!======================================================================
      Subroutine inf_sub
!======================================================================
!     provide screen information 
!----------------------------------------------------------------------
       
      Character :: A

      iarg = command_argument_count()
      if(iarg.eq.0) Return
      Call GET_COMMAND_ARGUMENT(1,A)        
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                              ',&
'     sec_cascb provides cross sections with cascade contribution              ',&
'                                                                              ',&
'     scattering information is given in file "target"                         ',&
'                                                                              ',&
'     branching ratios are taken from file "f_values"                          ',&
'                                                                              ',&                 
'     direct cross sections are taken from "tr_###_###"                        ',&
'                                  or from "cr_###_###"                        ',&
'                                                                              ',&
'     Arguments: itr1,itr2, jtr1,jtr2 - range for transition indexes           ',&
'                                       to be considered                       ',&
'                                       (itr -initial index, jtr - final)      ',&
'                                                                              ',&
'                cr - "tr" or "cr"  - suppose tr_###_### or cr_###_### files   ',&
'                                     with direct cross sections               ',&
'                                                                              ',&
'                steps   = 1 -> direct cascade                                 ',&
'                          2 -> through additional intermediate level          ',&
'                          5 -> maximum intermediate levels                    ',&
'                                                                              ',&
'     Call as:   sec_cascb  itr1=.. itr2=.. jtr1=.. jtr2=.. cr=.. steps=..     ',&     
'                                                                              ',&
'     Results:   ca_###_###, where ### stand for transition indexes            ',&
'                                                                              ',&
'     The folder should contain all needed tr- or cr-files, e.g.,              ',&
'     if we need cascade for transition from ground state,                     ',&
'     all  tr_001_###|cr_###_### should be presented                           ',&
'                                                                              ',&
'                                                                              '
      Stop ' '

      End Subroutine inf_sub

