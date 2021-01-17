!============================================================================
!      Differential cross sectins and scattering amplitudes (LS coupling) 
!      for scattering on neutral atoms with closed subshells (as He).
!      Considered all transitions from the ground state for one input energy.
!============================================================================
!                      zarm.tmb  +  target   -->  difsec_He
!============================================================================
!
!      Input:        T-matrixes  (zarm.tmb), target information (target)
!
!      Output:       difsec_He: amplitudes, diff. cross sections, 10^16 cm^2,
!                               angle-integrated cr.sec. for each final state      
!                    
!      Call:         dif_sec_He   ek=...  res=...  tm=...
!                                 g0=...  hg=...   ng=...   diff=...
!
!      This program, in contrast to the previous version, read all T-matrix 
!      in memeory. It may require a big memory !
!============================================================================
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)

      Integer, parameter :: ma = 80
      Character(ma) :: AF
      Integer :: nut=1; Character(ma) :: AF_targ = 'target_ps'
                        Character(ma) :: AF_tmat = 'zarm.tmb'
      Integer :: nur=2; Character(ma) :: AF_res  = 'difsec_He'
      Integer :: nus=3  ! scratch file
      Integer :: nud=4  ! DCS if any
      Integer :: pri=9; Character(ma) :: AF_log  = 'dif_sec_He.log'

      Integer, allocatable :: L1(:), L2(:), LJ(:), Lcase(:), Scase(:), &
                              kop(:),nop(:), kj(:)
      Real(8), allocatable :: TR(:),TI(:), matr(:,:), mati(:,:), fr(:,:),fi(:,:)

      Integer :: ng = 181
      Real(8), allocatable :: G(:), DCS(:), DD(:)

      Integer, external :: Def_ij

      Call Inf_sub

      PI  = abs(ACOS(-1.d0))
      open(pri,file=AF_log)

! ... set up the degree grid:

      g0 = 0.0d0;  Call Read_rarg('g0',g0)
      hg = 1.d00;  Call Read_rarg('hg',hg)
      ng = 181;    Call Read_iarg('ng',ng)
      Allocate(G(ng), DCS(ng), DD(ng))

      Do i=1,ng; G(i)=g0 + (i-1)*hg; End do

! ... target:

      Call Check_file(AF_targ)
      Open(nut,file=AF_targ)
      Call R_target(nut)
      Call R_channels(nut)
      np = ntarg; Call Read_ipar(nut,'np',np)
      ni = ntarg; Call Read_ipar(nut,'ni',ni)
      Close(nut)
      e1=etarg(1); Do i=1,ntarg; etarg(i)=(etarg(i)-e1)*2; End do
      ion = nz - nelc; if(ion.gt.0) Stop 'ion > 0'

      istyle=1; if(np.eq.ntarg.and.ni.eq.ntarg) istyle=0
                 
! ... T-matrix file:

      Call Read_aarg('tm',AF_tmat)
      Call Check_file(AF_tmat)
      Open(nut,file=AF_tmat)

      if(istyle.eq.0) then

       mtr = mch*(mch+1)/2

      else

      rewind(nut);  mtr=0                  
      Do 
       read(nut,*,end=20) ee,nopen,kopen,ilsp,i1,i2,nj
       ntr = kopen*(kopen+1)/2
       read(nut,*) (t1,t2,i=1,ntr)
       if(ntr.gt.mtr) mtr=ntr
       ktr = (nopen-kopen)*nj
       if(ktr.gt.0) then 
        read(nut,*) (t1,t2,i=ntr+1,ntr+ktr)
        if(ntr+ktr.gt.mtr) mtr=ntr+ktr
       end if
      End do
   20 Continue

      end if

      Allocate(matr(nlsp,mtr),mati(nlsp,mtr))
      Allocate(tr(mtr),ti(mtr))
      Allocate(nop(nlsp),kop(nlsp),kj(nlsp)); nop=0; kop=0; kj=0

! ... initial energy in Ry:

      ek=0.d0;  Call Read_rarg('ek',ek) 
      if(ek.le.0.d0) Stop 'define ek=...'

      idiff = 0;  Call Read_iarg('diff',idiff)
               
! ... initial state is supposed to be the ground states:

      itr1=1; IL1=ltarg(itr1);  IS1=istarg(itr1)

! ... output file:

      Call Read_aarg('res',AF_res)
      open(nur,file=AF_res)
      write(nur,'(a)') 'Differential cross sections for energy:'
      write(nur,'(a,f10.5)') 'EK =',EK
      write(nur,'(a,i5,a)') 'ng =',ng,'  -  number of degrees'
      write(nur,'(f8.2)') G

! ... read T-matrix:

      rewind(nut);  matr=0.d0; mati=0.d0                  
   1  if(istyle.eq.0) then
       read(nut,*,end=2) ee,nopen,ntr,ilsp   
       if(ntr.gt.mtr) Stop 'ntr > mtr'
       read(nut,*) (tr(i),ti(i),i=1,ntr)
      else
       read(nut,*,end=2) ee,nopen,kopen,ilsp,i1,i2,nj
       if(kopen.gt.nopen)     Stop ' kopen > nopen'
       ntr = kopen*(kopen+1)/2
       read(nut,*) (tr(i),ti(i),i=1,ntr)
       if(ntr.gt.mtr) Stop 'ntr > mtr'
       ktr = (nopen-kopen)*nj
       if(ntr+ktr.gt.mtr) Stop 'ntr +ktr > mtr'
       if(ktr.gt.0) read(nut,*) (tr(i),ti(i),i=ntr+1,ntr+ktr)
      end if

      if(ee.ne.ek) go to 1

      if(ilsp.gt.nlsp)       Stop ' ilsp  > nlsp'
      if(ilsp.lt.1)          Stop ' ilsp < 0'
      if(nopen.gt.nch(ilsp)) Stop ' nopen > nch(ilsp)'

      nop(ilsp) = nopen
      kop(ilsp) = kopen
      kj(ilsp) = nj 

      matr(ilsp,1:ntr+ktr)=tr(1:ntr+ktr)
      mati(ilsp,1:ntr+ktr)=ti(1:ntr+ktr)

      go to 1
    2 Continue
      Deallocate(tr,ti)

!----------------------------------------------------------------------
! ... main loop other states:

      Open(nus,FORM='UNFORMATTED',status='SCRATCH')

      Do itr2=1,ntarg; IL2=ltarg(itr2); IS2=istarg(itr2); IP2=iptarg(itr2)
       if(etarg(itr2).gt.ek) Exit

       write(*,*) 'itr2=',itr2           
       nstates = itr2

       ncase=0; rewind(nus)                  

       Do ilsp = 1,nlsp
        kopen = kop(ilsp)
        nopen = nop(ilsp)
        nj = kj(ilsp) 
                
       Do ich1=1,nopen
        if(iptar(ilsp,ich1).ne.itr1) Cycle
        Do ich2=1,nopen
         if(iptar(ilsp,ich2).ne.itr2) Cycle
      
         if(itr1.ne.itr2) then
          if(ich2.lt.ich1) Stop ' ich2 < ich1 '
          itr=ich2*(ich2-1)/2+ich1
         else
          i1=min(ich1,ich2); i2 = max(ich1,ich2); itr=i1*(i1-1)/2+i2
         end if

         if(istyle.eq.1.and.ich2.gt.kopen) &
          itr = kopen*(kopen+1)/2 + nj*(ich2-kopen-1) + ich1

         ncase = ncase + 1
         write(nus) ispar(ilsp),lpar(ilsp),lch(ilsp,ich1),lch(ilsp,ich2), &
                    matr(ilsp,itr),mati(ilsp,itr)
        End do
       End do

       End do  ! ilsp

       if(ncase.eq.0) Stop 'number of T-elements = 0'

       Allocate(TR(ncase),TI(ncase), l1(ncase), l2(ncase), &
                Lcase(ncase), Scase(ncase))
       rewind(nus)
       Do i = 1,ncase
        read(nus) Scase(i),Lcase(i),l1(i),l2(i),TR(i),TI(i)
       End do

! ... m-dependent amplitudes:

       if(allocated(fr)) Deallocate(FR); Allocate(FR(-IL2:IL2,ng)); FR=0.d0
       if(allocated(fi)) Deallocate(FI); Allocate(FI(-IL2:IL2,ng)); FI=0.d0

       Do i=1,ncase

        S= DFLOAT(2*l1(i)+1); S=sqrt(S)

        Do m = -IL2,IL2   ! over small m
         SM = S * CLEBCH(IL2,-m,l2(i),m,l1(i),0)
         Do ig=1,ng; x = DCOS(G(ig)/180*pi)
          SG = SM * ALEGFM (L2(i),m,X,1)
          ar = TR(I)*SG; ai = TI(I)*SG; 
          Call Fano_tma (l1(i),l2(i),ar,ai,br,bi)     !    i ** (l1-l2)        
          FR(-m,ig) = FR(-m,ig)+br; FI(-m,ig) = FI(-m,ig) + bi; 
         End do
        End do

       End do

! ... diff. cross-sections:

       SN = 1.d0/(2*ek)  * 0.529177**2     !  in 10-16

       Do ig = 1,ng
        S = 0.d0
        Do m = -IL2,IL2 
         S = S + FR(m,ig)*FR(m,ig) +  FI(m,ig)*FI(m,ig) 
        End do
        DCS(ig) = S * SN
        DD(ig) = DCS(ig) * DSIN(G(ig)/180*pi)
       End do 
       
       ST = ZIS(hg,ng,dd) * pi/180 * 2*pi

       EP = E1 + etarg(itr2)/2
       write(nur,'(a,i5,E20.8,3i5,E15.5)') 'is,E,ilsp,sec =', &
                                           itr2,EP,IL2,IS2,IP2,ST
       Do ig=1,ng
        write(nur,'(f10.2,E15.5,10(2x,2E15.5))') G(ig),DCS(ig), &
                                  (FR(m,ig),FI(m,ig),m=-IL2,IL2)
       End do

       Deallocate(TR, TI, l1, l2, Lcase, Scase )

       if(idiff.eq.0) Cycle

       if(idiff.gt.0.and.idiff.ne.itr2) Cycle
       write(AF,'(a,i3.3)') 'dcs_',itr2
       open(nud,file=AF)
       write(nud,'(a,a,f10.6,a,E15.5)') '  angle   DCS    (in 10^-16 cm^2) ', &
            '  k^2 = ', ek, ' sigma =', ST
       Do ig = 1,ng
        write(nud,'(f9.2,E16.5)') G(ig),DCS(ig)
       End do

      End do   ! over states (itr2)
!---------------------------------------------------------------
      write(nur,'(/a,i5/)') 'nstates =',nstates
      lmax = maxval(ltarg)
      write(nur,'(/a,i5/)') 'lmax =',lmax
      ismax = maxval(istarg)
      write(nur,'(/a,i5/)') 'smax =',ismax

      End  ! program dif_secT_He


!======================================================================
      Subroutine inf_sub
!======================================================================
!     provide screen information about add_farm utility
!----------------------------------------------------------------------
      Character(80) :: A

      Call get_command_argument(1,A)  
      if(A.ne.'?') Return

      write(*,'(a)') &
'                                                                  ',&
'This utility calulates DCS and scattering amplitudes (LS coupling)',&
'for scattering on neutral atoms with closed subshells (as He).    ',&
'                                                                  ',&
'Input:     T-matrixes (zarm.tmb), target information (target)     ',&
'                                                                  ',&
'Call as:   dif_sec_He   ek=..  res=..  tm=.. g0=.. hg=.. ng=..    ',&
'                                                                  ',&
'ek         - initial electron energy in Ry                        ',&
'                                                                  ',&
'g0 [0]     - initial degree                                       ',&
'hg [1]     - step over defrees                                    ',&
'ng [181]   - number of dergee points                              ',&
'                                                                  ',&
'tm [zarm.tmb] - T-matrix file                                     ',&
'res [difsec_He] - file for final DCS and scattering amplitudes    ',&
'                                                                  '
    Stop 
                         
    End Subroutine inf_sub
