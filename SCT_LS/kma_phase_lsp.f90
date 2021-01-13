!======================================================================
!     UTILITY     K M A _ P H A S E _ L S P   
!
!     Eigenphases for given partial wave (based on the given K-matrix)
!
!     target, zarm.kma -->   phases.nnn
!
!     arguments: klsp - partial-wave index
!
!     remarks: if K-matrix file contains a few results for same
!              energy and ilsp, the program will use the last one
!======================================================================

      Use target; Use channels

      Implicit real(8) (a-h,o-z)

      Integer, parameter :: ke = 5000
      Integer, allocatable :: ipe(:)
      Real(8), allocatable :: e(:)

      Real(8) :: Ry = 13.6057

      Real(8), allocatable :: matr(:), ui(:)
      Real(8), allocatable :: kma(:,:), fl(:,:), uk(:,:)
      Integer, allocatable :: ifl(:,:)
      Character(4) :: Term
      Character(1), External :: AL

! ... files:

      Integer :: nut = 1;    Character(20) :: targ  = 'target'       
      Integer :: nuk = 2;    Character(20) :: kname = 'zarm.kma'     
      Integer :: nup = 3;    Character(20) :: pname = 'phases'       

      Integer :: nua = 99   ! scratch file

      Call get_command_argument(1,term)  
      if(term.eq.'?'.or.term.eq.'!') then

      write(*,'(a)') &
'                                                                      ',& 
'     kma_phase_LSP provides channel EIGENPHASES for given pratial waves ',& 
'                                                                      ',& 
'     zarm.kma + target  --> phahes                                    ',& 
'                                                                      ',& 
'     Call as:   kma_phase_lsp  [klsp=.. ek1=.. ek2=.. ]               ',& 
'                                                                      ',& 
'     klsp [1]      - index of partial wave                            ',& 
'                                                                      ',& 
'     ek1 [0]       - minimum electron energy in Ry                    ',& 
'     ek2 [1000]    - minimum electron energy in Ry                    ',& 
'                                                                      ',& 
'                                                                      ' 
      Stop ' '
      End if

!----------------------------------------------------------------------
! ... target information:

      Open(nut,file='target',status='OLD')
      Call R_target(nut)
      Call R_channels(nut)
      Close(nut)
      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion
      E1 = etarg(1); etarg = (etarg-E1)*2.d0

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

!-----------------------------------------------------------------------
! ... define partial waves under consideration:

      klsp = 1;   Call Read_iarg('klsp',klsp)

      ek1=0;      Call Read_rarg('ek1',ek1)
      ek2=1000;   Call Read_rarg('ek2',ek2)

! ... allocations:

      kch = nch(klsp)
      Allocate(kma(mch,mch), uk(mch,mch), ui(mch) )
      mdim = mch*(mch+1)/2; Allocate(matr(mdim) )
      me = ke;  Allocate(e(me),fl(me,0:kch),ifl(me,0:kch))
      e = 0.d0; fl = 0.d0; ifl = 0

!----------------------------------------------------------------------
! ... read K-matrices:

      Open(nuk,file=kname, status='OLD')

      ne=0
    1 read(nuk,*,end=2) e1,nopen,ntr,ilsp
      read(nuk,*) matr(1:ntr)

      if(e1.gt.ek2) go to 1
      if(e1.lt.ek1) go to 1

      if(ilsp.ne.klsp) go to 1

      ie=0; Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

! ... new energy:

      if(ie.eq.0) then; ne=ne+1; ie=ne; e(ie)=e1; end if

      if(ne.eq.me) then
       Open(nua,form='UNFORMATTED',status='SCRATCH')
       write(nua) (e(i),i=1,ne)
       write(nua) ((fl(i,j),i=1,ne),j=0,mch)
       Deallocate(fl,e); me=me+ke; Allocate(e(me),fl(me,0:mch)) 
       fl = 0.d0
       rewind(nua)
       read(nua) (e(i),i=1,ne)
       read(nua) ((fl(i,j),i=1,ne),j=0,mch)
       Close(nua)
      end if

! ... partial and total eigenphases:

      Do i=1,nopen        
       Do j=1,i
        ij=(i-1)*i/2+j
        KMA(i,j)=matr(ij)
        KMA(j,i)=matr(ij)
       End do
      End do

      Call EPHASE (nopen,mch,KMA,uk,ui,us)

      fl(ie,0) =  us   
      uk = uk * uk

      Do j=1,nopen
       S=0.d0; m=0
       Do i=1,nopen
        if(uk(i,j).gt.S) then; S=uk(i,j); m=i; end if
       End do
       fl(ie,m) = ui(j)
      End do
  
      go to 1
    2 Continue

!----------------------------------------------------------------------
! ... define the partial wave term:

      if(ispar(klsp).eq.0) then
       if(ipar(klsp).eq.+1) write(term,'(a1,i2.2,a1)') 'A', lpar(klsp),'e'
       if(ipar(klsp).eq.-1) write(term,'(a1,i2.2,a1)') 'A', lpar(klsp),'o'
      else
       if(ipar(klsp).eq.+1) write(term,'(a1,i1,3a1)') 'A',ispar(klsp),AL(lpar(klsp),2),'e'
       if(ipar(klsp).eq.+1) write(term,'(a1,i1,3a1)') 'A',ispar(klsp),AL(lpar(klsp),2),'o'
      end if

! ... energy order:

      Allocate(ipe(ne));  Call SortR(ne,E,ipe)

!----------------------------------------------------------------------
! ... adjust the phases:

      a = 0.5D0

      Do ich = 0,kch

      i = ipe(1)
    3 if(fl(i,ich).lt.0.d0) then
       fl(i,ich) = fl(i,ich) + 1.d0
       go to 3
      end if 

      Do ie=2,ne
       i =ipe(ie);   P =fl(i, ich)
       i1=ipe(ie-1); P1=fl(i1,ich)                             

       Do
        PP = P + a
        a1=abs(P-P1)
        a2=abs(PP-P1)
        if(a2.ge.a1) Exit
        P = PP
       End do

       Do
        PP = P - a
        a1=abs(P-P1)
        a2=abs(PP-P1)
        if(a2.ge.a1) Exit
        P = PP
       End do

       fl(i,ich) = P

      End do; End do

!----------------------------------------------------------------------
! ... output results:

 if(kch.gt.50) kch =50     !  ???

      i = LEN_TRIM(pname) + 1
      write(pname(i:i+3),'(a,i3.3)') '.',klsp

      Open(nup,file=pname)

      kopen = jopen(e(ipe(ne)),klsp)

if(kopen.gt.50) kopen=50

      write(nup,'(6x,a2,10x,a2,8x,50(2x,a4,4x))') &
                    'Ry','eV',term,(ELC(klsp,i),i=1,kopen)


!      fl = fl * 3.1415926  !  convert to radians

      Do ie=1,ne; j=ipe(ie);  fl(j,0)=sum(fl(j,1:kopen))

       write(nup,'(2f12.6,25(f10.3))') e(j),e(j)*Ry,(fl(j,i),i=0,kopen)

      End do

      Close(nup)


      End  ! utility kma_phase



!======================================================================
      Subroutine ephase (n,m,K,uk,ui,us)
!======================================================================
!     n - number of open channels
!     m - number of channels (max.dimension)
!     K - symmetrical matrix
!     uk - eigenvectors of K-matrix 
!     ui - eigenphases (on module of PI) ->  tan^-1 (eigenvalue) / PI
!     us - sum of eigenphases (on module of PI)
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in)  :: n,m
      Real(8), intent(in)  :: K(m,m)
      Real(8), intent(out) :: uk(m,m),ui(m),us
      Real(8) :: pi
      Integer :: i,info

! ... get the eigenvalues and eigenvectors of the K-matrix:

      uk=K;   Call LAP_DSYEV('V','L',n,m,uk,ui,info)

      if(info.ne.0) then
       ui=0; us=0; Return 
      end if

! ... eigenphases:

      pi = dacos(-1.d0)
      us = 0.d0
      Do i=1,n
       ui(i)  = datan(ui(i))/pi
       us = us + ui(i)
      End do

      End Subroutine ephase 

