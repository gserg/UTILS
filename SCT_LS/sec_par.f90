!======================================================================
!     UTILITY     S E C _ P A R   
!
!     provides partial cross sections for given transition i1 -> i2
!
!     target, zarm.kma -->   pa_i1_i2  
!
!     arguments: i1,i2  -  transition
!                ilsp1,ilsp2 - range of partial waves
!                i16 - 0: cross sections is atomic units (default) 
!                    - 1: in 10^-16
!                    - 2: in 10^-18 
!
!     remarks: if K-matrix file contains a dublicated results,
!              the program will use the last one
!
!======================================================================

      Use target;  Use channels

      Implicit real(8) (a-h,o-z)

      Integer, parameter   :: ke = 5000
      Integer, allocatable :: ipe(:)
      Real(8), allocatable :: e(:)

      Real(8) :: Ry = 13.6057

      Real(8), allocatable :: matr(:)
      Real(8), allocatable :: POM(:,:), fl(:,:)
      Integer, allocatable :: ifl(:)
      Character(4), allocatable :: Term(:)

      Logical EX

      Character(1) :: AL
      Character(3) :: AT

! ... files:

      Character(20) :: targ  = 'target';       Integer :: nut = 1
      Character(20) :: kname = 'zarm.kma';     Integer :: nuk = 2
      Character(20) :: pname = 'pa_00_00.dat'; Integer :: nup = 3

      Integer(4) :: nua = 99   ! scratch file


      Call inf_sec_par

!----------------------------------------------------------------------
! ... target information:

      Inquire(file=targ,EXIST=EX)
      if(.not.EX) then
       Stop  ' No target file: run H_targ ! '
      else
       Open(nut,file='target')
       Call R_target(nut)
       Call R_channels(nut)
       close(nut)
      end if
      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion
      E1 = etarg(1); etarg = (etarg-E1)*2.d0

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

!-----------------------------------------------------------------------
! ... define transition under consideration:

      iarg = IARGC()
      if(iarg.lt.4) &
       Stop 'Should be 5 arguments: itr1 itr2 ilsp1 ilsp2 [i16]'

      Call GETARG(1,AT); read(AT,*) itr1
      Call GETARG(2,AT); read(AT,*) itr2
      Call GETARG(3,AT); read(AT,*) ilsp1
      Call GETARG(4,AT); read(AT,*) ilsp2
      i16 = 0
      if(iarg.ge.5) then; Call GETARG(5,AT); read(AT,*) i16; end if

      if(itr1.le.0) itr1=1; if(itr1.gt.ntarg) itr1=ntarg
      if(itr2.le.itr1) itr2=itr1; if(itr2.gt.ntarg) itr2=ntarg

! ... name for output files:

      write(pname,'(a,i3.3,a,i3.3,a)') 'pa_',itr1,'_',itr2 
      open(nup,file=pname)

! ... allocations:

      Allocate(POM(mch,mch),term(nlsp))
      mdim = mch*(mch+1)/2; Allocate(matr(mdim))
      me = ke;  Allocate(e(me),fl(me,nlsp),ifl(nlsp))
      e = 0.d0; fl = 0.d0

!----------------------------------------------------------------------
! ... read K-matrices:

      Inquire (file=kname,EXIST=EX)
      if(.not.EX) Stop ' K-matrix file is absent '
      Open(nuk,file=kname)

      ne=0
    1 read(nuk,*,end=2,err=2) e1,nopen,ntr,ilsp
      read(nuk,*,err=1) matr(1:ntr)

      if(ilsp.lt.ilsp1.or.ilsp.gt.ilsp2) go to 1

      ie=0
      Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

! ... new energy:

      if(ie.eq.0) then; ne=ne+1; ie=ne; e(ie)=e1; end if

      if(ne.eq.me) then
       Open(nua,form='UNFORMATTED',status='SCRATCH')
       write(nua) (e(i),i=1,ne)
       write(nua) ((fl(i,j),i=1,ne),j=1,nlsp)
       Deallocate(fl,e); me=me+ke; Allocate(e(me),fl(me,nlsp)) 
       rewind(nua)
       read(nua) (e(i),i=1,ne)
       read(nua) ((fl(i,j),i=1,ne),j=1,nlsp)
       Close(nua)
      end if

! ... partial omega:

      g = (2.0*lpar(ilsp)+1) * iabs(ispar(ilsp)) / 2.0
      if(ispar(ilsp).eq.0) g = (lpar(ilsp)+1)/2.d0
      Call kma_om (mch,nopen,g,matr,POM)
      fl(ie,ilsp) = 0.d0
      Do i=1,nopen;  if(iptar(ilsp,i).ne.itr1) Cycle
      Do j=1,nopen;  if(iptar(ilsp,j).ne.itr2) Cycle
       fl(ie,ilsp) = fl(ie,ilsp) + POM(i,j)
      End do; End do

      go to 1
    2 Continue

!----------------------------------------------------------------------
! ... define terms for partial waves:

      Do i=1,nlsp
       if(ispar(i).eq.0) then
        i1=lpar(i)/10; i2=mod(lpar(i),10)
        write(term(i),'(a1,2i1,a1)') 'A', i1,i2
       else
        write(term(i),'(a1,i1,a1,a1)') 'A', &
                     iabs(ispar(i)),AL(lpar(i),2)
       end if
       if(ipar(i).eq. 1) term(i)(4:4)='e'
       if(ipar(i).eq.-1) term(i)(4:4)='o'
      End do

! ... energy order:

      Allocate(ipe(ne));  Call SortR(ne,E,ipe)

! ... statistical weight for initial state:

      g=iabs(IStarg(itr1))*(2.0*Ltarg(itr1)+1)
      if(IStarg(itr1).eq.0) g=Ltarg(itr1)+1

! ... choose non-zero partial waves:

      ifl=0; klsp = 0
      Do i = ilsp1,ilsp2
       if(SUM(fl(1:ne,i)).eq.0.d0) Cycle
       klsp = klsp + 1
       ifl(klsp) = i
      End do 

!----------------------------------------------------------------------
! ... output results:

      write(nup,'(6x,a2,10x,a2,10x,a3,5x,65(7x,a4,5x))') &
                    'eV','Ry','Sum',(term(ifl(i)),i=1,klsp)

      Do ie=1,ne; j=ipe(ie)
       S = SUM(fl(j,ilsp1:ilsp2))

        Do k=ilsp1,ilsp2
         if(fl(j,k).eq.0.d0) write(*,'(a,f10.6,a,i3,a)') 'e = ', e(j), &
          'partial wave =',k, '   data in zarm.kma are absent' 
        End do

       if(S.eq.0.d0) Cycle

       es=e(j)-etarg(itr1); es=es*zion

       SN = 3.1415926/g/es                    ! in ao^2
       if(i16.eq.1) SN = SN * 0.28003         ! in 10-16
       if(i16.eq.2) SN = SN * 28.003          ! in 10-18
       if(i16.eq.3) SN = SN * 2800.3          ! in 10-20 (pm)

       s = s * SN                

       fl(j,ilsp1:ilsp2) = fl(j,ilsp1:ilsp2) * SN

       i1=1; i2=klsp
       write(nup,'(2f12.6,64D16.8)') &
                   e(j)*Ry,e(j),S,(fl(j,ifl(i)),i=i1,i2)

      End do

      End  ! utility sec_par



!======================================================================
      Subroutine inf_sec_par
!======================================================================
!     provide screen information about sec_par utility
!----------------------------------------------------------------------
       
      Implicit real(8) (a-h,o-z)

      Character :: A

      Call get_command_argument(1,A)  
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                          ',&
'     SEC_PAR provides partial cross sections for given transition i1 -> i2',&
'                                                                          ',&
'                  zarm.kma + target  =>   par_i1_i2                       ',&
'                                                                          ',&
'     Call as:  sec_par  i1 i2 ilsp1 ilsp2 [i16]   (position arguments)    ',&
'                                                                          ',&
'     ilsp1,ilp2 - range of partial waves                                  ',&
'                                                                          ',&
'     i16 - control the output:                                            ',&
'           0 - sigma in a.u.  (default)                                   ',&
'           1 - sigma in 10-16 cm^2                                        ',&
'           2 - sigma in 10-18 cm^2                                        ',&
'                                                                          '
      Stop ' '

      End Subroutine inf_sec_par

