!======================================================================
!     UTILITY     K M A _ P H A S E   
!
!     gives LSP eigenphases form given K-matrix
!
!     target, zarm.kma -->   phases
!
!     arguments: ilsp1,ilsp2 - range of partial waves
!
!     remarks: if K-matrix file contains a few results for same
!              energy and ilsp, the program will use the last one
!
!======================================================================
      Use target; Use channels                                                            

      Implicit real(8) (a-h,o-z)
      Integer :: ke = 5000
      Integer, allocatable :: ipe(:)
      Real(8), allocatable :: e(:)

      Real(8), allocatable :: matr(:), ui(:)
      Real(8), allocatable :: kma(:,:), fl(:,:), uk(:,:), ep(:), fl1(:), fd1(:)
      Integer, allocatable :: ifl(:)
      Character(5), allocatable :: Term(:)

      Logical EX

      Character(1) :: AL
      Character(3) :: AT

! ... files:

      Character(20) :: targ  = 'target';       Integer :: nut = 1
      Character(20) :: kname = 'zarm.kma';     Integer :: nuk = 2
      Character(20) :: pname = 'phases';       Integer :: nup = 11

      Integer :: nua = 99   ! scratch file

      Call get_command_argument(1,AT)  
      if(AT.eq.'?'.or.AT.eq.'!') then

      write(*,'(a)') &
'                                                                      ',& 
'     kma_phase provides SUM-of-EIGENPHASES for selected pratial waves ',& 
'                                                                      ',& 
'     zarm.kma + target  --> phahes.nnn_mmm  or  phases.nnn            ',& 
'                                                                      ',& 
'     Call as:   kma_phase   [ilsp1=.. ilsp2=.. E_min=.. E_max=.. ]    ',& 
'                                                                      ',& 
'     ilsp1 [1]     - index of first partial wave                      ',& 
'     ilsp2 [ntarg] - index of last partial wave                       ',& 
'                                                                      ',& 
'     E_min [0]     - minimum electron energy in Ry                    ',& 
'     E_max [1000]  - minimum electron energy in Ry                    ',& 
'                                                                      ',& 
'     if ilsp1=ilsp2 - we have additional output phase.nnn (nnn=ilsp1) ',& 
'                      with derivatives and gamma                      ',& 
'                                                                      ' 
      Stop ' '
      End if

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
! ... define partial waves under consideration:

      ilsp1=1;     Call Read_iarg('ilsp1',ilsp1)
      ilsp2=ntarg; Call Read_iarg('ilsp2',ilsp2)

      E_max=1000;   Call Read_rarg('E_max',E_max)
      E_min=0.d0;   Call Read_rarg('E_min',E_min)

! ... allocations:

      Allocate(kma(mch,mch), uk(mch,mch), ui(mch), term(nlsp) )
      mdim = mch*(mch+1)/2; Allocate(matr(mdim))
      me = ke;  Allocate(e(me),fl(me,nlsp),ifl(nlsp))
      e = 0.d0; fl = 0.d0

!----------------------------------------------------------------------
! ... read K-matrices:

      Inquire (file=kname,EXIST=EX)
      if(.not.EX) Stop ' K-matrix file is absent '
      Open(nuk,file=kname)

      ne=0
    1 read(nuk,*,end=2) e1,nopen,ntr,ilsp
      read(nuk,*) matr(1:ntr)

      if(e1.gt.E_max) go to 1
      if(e1.lt.E_min) go to 1


      if(ilsp.lt.ilsp1.or.ilsp.gt.ilsp2) go to 1

      ie=0; Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

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

! ... partial and total eigenphases:

      Do i=1,nopen        
       Do j=1,i
        ij=(i-1)*i/2+j
        KMA(i,j)=matr(ij)
        KMA(j,i)=matr(ij)
       End do
      End do

      Call EPHASE (nopen,mch,KMA,uk,ui,us)

      fl(ie,ilsp) =  us   

      go to 1
    2 Continue

!----------------------------------------------------------------------
! ... define terms for partial waves:

      Do i=1,nlsp
       if(ispar(i).eq.0) then
        i1=lpar(i)/10; i2=mod(lpar(i),10)
        write(term(i),'(a1,2i1,a1)') 'A', i1,i2,'_'
       else
        write(term(i),'(a1,a1,i1,a1,a1)') 'k',AL(lch(i,1),1), &
                     iabs(ispar(i)),AL(lpar(i),2)
       end if
       if(ipar(i).eq.+1) term(i)(5:5)='e'
       if(ipar(i).eq.-1) term(i)(5:5)='o'
      End do

! ... energy order:

      Allocate(ipe(ne));  Call SortR(ne,E,ipe)

! ... choose non-zero partial waves:

      ifl=0; klsp = 0
      Do i = ilsp1,ilsp2
       if(SUM(fl(:,i)).eq.0.d0) Cycle
       klsp = klsp + 1
       ifl(klsp) = i
      End do 

!----------------------------------------------------------------------
! ... adjust the phases:               
                                       
      a = 0.5D0                        
                                       
      Do ilsp = ilsp1,ilsp2            
                                       
      Do ie=2,ne                       
       i=ipe(ie); P=fl(i,ilsp)         
       i1=ipe(ie-1); P1=fl(i1,ilsp)                             
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
                                       
       fl(i,ilsp) = P                  
                                       
      End do; End do                   
                                       
!----------------------------------------------------------------------
! ... output results:


      if(ilsp1.lt.ilsp2) then

       write(pname,'(a,a,i3.3,a,i3.3)') trim(pname),'.',ilsp1,'_',ilsp2
       Open(nup,file=pname)
       write(nup,'(6x,a2,10x,a2,3x,125(5x,a5,4x))') &
                     'Ry','eV',(term(ifl(i)),i=1,klsp)
     
       Do ie=1,ne; j=ipe(ie)
     
        S = SUM(fl(j,ilsp1:ilsp2)); if(S.eq.0.d0) Cycle
     
        i1=1; i2=klsp
     
        write(nup,'(2f12.6,125E14.6)') e(j),e(j)*Ry,(fl(j,ifl(i)),i=i1,i2)
     
       End do
     
       Close(nup)

      else

       Allocate(fl1(ne),fd1(ne),ep(ne))
     
       ke=0
       Do ie=1,ne; j=ipe(ie)
     
        S = SUM(fl(j,ilsp1:ilsp2)); if(S.eq.0.d0) Cycle
     
        i1=1; i2=klsp
     
        ke = ke + 1;  fl1(ke)=fl(j,ifl(i1)); ep(ke)=e(j)
     
       End do
     
! ... derivatives:

       Call DERIVATIVE1_NONUNIFORM(ke,ep,fl1,fd1)
     
       write(pname,'(a,a,i3.3)') trim(pname),'.',ilsp1
     
       Open(nup,file=pname)
     
       write(nup,'(6x,a2,10x,a2,3x,3(5x,a5,4x),5x,a)') &
                 'Ry','eV',term(ilsp1),'deriv','gamma','phase in PI, deriv over ek in Ry, gamma in mev'
     
       write(nup,'(2f12.6,3E14.6)') (ep(j),ep(j)*Ry,fl1(j),fd1(j),3.14*2/fd1(j)*Ry*1000,j=1,ke)
     
       Close(np) 
     
       end if

      End  ! utility kma_phase



!===================================================================
  SUBROUTINE DERIVATIVE1_NONUNIFORM (N,X,F,F1)
!===================================================================
! 1st order derivatives with the three-point non-uniform formula. 
! Linear extrapolations are made at the boundaries.  
! X: input x; F: input f(x);  F1: f'
!===================================================================

  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I
  REAL(8), INTENT (IN)  :: X(N),F(N)
  REAL(8), INTENT (OUT) :: F1(N)
  REAL(8) :: h,hh, h1, hh1

! f' from three-point non-uniform formulas

  h1 = X(2)-X(1);  hh1=h1*h1
  DO i = 2, N-1
    h = X(i+1)-X(i); hh = h*h
    F1(i) = (hh1*F(i+1)+(hh-hh1)*F(i)-hh*F(i-1))/(h*h1*(h+h1))
    h1 = h; hh1 = hh
  END DO

! Linear extrapolation for the boundary points

  F1(1) = 2.d0*F1(2)-F1(3)
  F1(N) = 2.d0*F1(N-1)-F1(N-2)

  END SUBROUTINE DERIVATIVE1_NONUNIFORM
