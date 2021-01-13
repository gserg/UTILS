!=====================================================================
      Module tdcs_amplitude     
!=====================================================================
!     contains main arrays for calculation of tdcs amplitudes
!---------------------------------------------------------------------
      Use target_ion     
      Use channels_ion

      Implicit none

      Integer :: ng = 180                          ! number of angles
      Real(8) :: hg = 1                            ! angle step

      Integer :: Lmax = 0                          ! max. L for (e-A+) system
      Integer :: Smax = 0                          ! max. S for (e-A+) system
      Integer :: llmax = 0                         ! max. small l for ejected electrons
      Integer :: LT_max = 0                        ! max. L for ion states

      ! this above information can be get from target and channels modules

      Real(8), allocatable :: fi(:,:,:,:,:)        ! projected amplitudes (direct)
      Real(8), allocatable :: fj(:,:,:,:,:)        ! projected amplitudes (exchange)

      ! fi, fj ->  (1:2,-lmax:lmax,mch,nlsp,0:ng)

      Real(8), allocatable :: FSM1(:,:,:,:)      ! direct  ionization amplitudes 
      Real(8), allocatable :: FSM2(:,:,:,:)      ! exchage ionization amplitudes 

      ! FSM  ->  (1:2,-LTmax:LTmax,smax,ntarg)     ! for given angles: teta1,fi1,teta2,fi2        
      
      Real(8), allocatable :: DCS(:)               ! dcs(1:ntarg) - dcs for given angles 

      Real(8), allocatable :: CP1(:,:), CP2(:,:)   ! Coulomb phases:

      ! CP   ->  (0:llmax,ntarg)                 

      Real(8), allocatable :: plmn(:,:,:)          ! normalized Legendre functions

      ! plmn ->  (0:ng,0:llmax,-llmax:llmax)  

      Real(8) :: pi

      End Module tdcs_amplitude


!=====================================================================
!     UTILITY  TDCS
!=====================================================================
!
!     triple differential cross section for ionization of He or other
!     closed-shell atoms with projection method
!
!     Input files:
!
!     target_ion  - description of (e-A+) scattering
!     difsec_He   - collection of scattering amplitudes for (e+A) 
!                   for each pseudosctate;
!                   here comes energy of incident electron - EK0
!     projection2 - projection coefficients of each psedostate to 
!                   to the continuum (e-A+) at EL2 = E_ejected 
!                   and at  at EL1 = E_scattering
!                   (projection files define the EL2 and EL1 energies)
!
!     Output files:
!
!     ddcs_E0_E2      
!
!---------------------------------------------------------------------
      Use tdcs_amplitude

      Implicit real(8) (A-H,O-Z)

      Character(200) :: AS,AU

      Real(8), allocatable :: g(:)                 ! teta grid:  -180 : 180
      Real(8), allocatable :: f(:,:,:,:)           ! amplidutes for p.s. 

      Real(8), allocatable :: tdcs(:,:)            ! (ig2,itarg) or (ig2,ifi)

      Real(8), allocatable :: Eps(:)               ! energies of p.s.
      Integer, allocatable :: Lps(:),Sps(:),Pps(:) ! terms of p.s.

      ! ... projections:

      Real(8), allocatable :: wt1(:,:), di1(:,:), dr1(:,:)
      Real(8), allocatable :: wt2(:,:), di2(:,:), dr2(:,:)

      ! ... files involved:

      Integer :: nut=1;  Character(40) :: AF_i   = 'target_ion'
                         Character(40) :: AF_p   = 'projection2'
                         Character(40) :: AF_dcs = 'difsec_He'

      Integer :: out=3;  Character(40) :: AF_td  = 'tdcs.out'

      Real(8) :: Ry = 13.6056
      Real(8) :: k0,k1,k2         ! electrons momentums

      Call inf_tdcs
       
      Call CPU_time(t1)

      pi  = abs(ACOS(-1.d0))

!----------------------------------------------------------------------
! ... define input parameters:

      ig1  = 10;      Call Read_iarg('ig'   ,ig1  )
      ifi  =  0;      Call Read_iarg('fi'   ,ifi  )
      mode =  2;      Call Read_iarg('mode' ,mode )
      itarg=  1;      Call Read_iarg('itarg',itarg)
      igeom=  1;      Call Read_iarg('geom' ,igeom) 

      idte =  1;      Call Read_iarg('idte' ,idte )      
      idfi =  1;      Call Read_iarg('idfi' ,idfi )      
      idel = 90;      Call Read_iarg('delta',idel ) 
      i16  =  0;      Call Read_iarg('i16'  ,i16  )

      icorr=  1;      Call Read_iarg('corr' ,icorr) 
      iexch=  1;      Call Read_iarg('exch' ,iexch) 
      icong=  1;      Call Read_iarg('cong' ,icong) 

      if(ig1.lt.-180.or.ig1.gt.180) then
       write(*,*) 'we suppose  -180 < ig1 < 180 '
       Stop 'Call as:   tdcs_all ig=20 i16=20 '
      end if
      iga = iabs(ig1)
      write(*,'(a,i5)') 'teta_1 = ',ig1
      write(*,'(a,i5)') 'ifi    = ',ifi
      write(*,'(a,i5)') 'icorr = ',icorr
      write(*,'(a,i5)') 'mode  = ',mode
      write(*,'(a,i5)') 'iexch = ',iexch
      write(*,'(a,i5)') 'icong = ',icong
      write(*,'(a,i5)') 'geom  = ',igeom
      write(*,'(a,i5)') 'delta = ',idel

!----------------------------------------------------------------------
! ... ion target information:

      Call Read_aarg('ion',AF_i)
      Call Check_file(AF_i)
      Open(nut,file=AF_i)
      Call R_target_ion(nut)         !   target_ion - target    ???
      Call R_channels_ion(nut)
      Close(nut)

      Z = nz;  ZC = Z - nelc       

!----------------------------------------------------------------------
! ... read amplidudes for all pseudostates at given angle: 
! ... we supposed universal angle grid: 0,1,2,3,...,180

      Call Read_aarg('dcs',AF_dcs)
      Call Check_file(AF_dcs)
      Open(nut,file=AF_dcs)

      Call Read_rpar(nut,'EK',EK)      
      write(*,'(a,f10.5)') 'EK =',EK
      k0 = sqrt(EK)

      Call Read_ipar(nut,'ng',ng)      
      ng = ng-1
      write(*,'(a,i10)') 'ng = ', ng
      if(ng.ne.180) Stop 'ng <> 180 in dif_sec file'
      Allocate(g(0:ng))
      read(nut,*) g(0:ng)

      Call Read_ipar(nut,'nstates',nps)      
      Call Read_ipar(nut,'lmax',lmax)      
      Call Read_ipar(nut,'smax',smax)      

      Allocate(Eps(nps),Lps(nps),Sps(nps),Pps(nps))
      Allocate(F(2,-lmax:lmax,nps,0:ng))

      F=0.d0
      rewind(nut)
    1 read(nut,'(a)',end=2) AS
      if(AS(1:2).ne.'is') go to 1
      i=INDEX(AS,'=')+1 
      read(AS(i:),*) is,Eps(is),Lps(is),Sps(is),Pps(is); l=Lps(is)
      Do i=0,ng; read(nut,*) s,d,((f(k,m,is,i),k=1,2),m=-l,l); End do
      go to 1
    2 Continue
      Close(nut)

!----------------------------------------------------------------------
! ... index of ionizing pseudostates:

      iion = 0
      Do it=1,nps
       if(Eps(it).gt.Etarg(1)) Exit;  iion = it
      End do
      write(*,'(a,i3)') 'bound p.s. = ', iion
      Call Read_iarg('iion',iion) 
      write(*,'(a,i3)') 'bound p.s. = ', iion

      E_ion = (Etarg(1)-Eps(1))*2*Ry
      write(*,'(a,f10.2)') 'E_ion = ', E_ion

!----------------------------------------------------------------------
! ... projection coefficient (direct ionization):

      Call Read_aarg('p',AF_p)
      Call Check_file(AF_p)
      Open(nut,file=AF_p)

      Call Read_rpar(nut,'EK',EKK)      
      if(abs(EK-EKK).gt.0.001) Stop 'different EK in dif_sec and projection files'

      Call Read_ipar(nut,'nps',npp) 
      if(npp.ne.nps) Stop 'different nps in dif_sec and projection files' 

      Call Read_ipar(nut,'ion',itarg) 
      write(*,'(a,i5)') 'itarg = ',itarg

      Call Read_rpar(nut,'EL2',EL1) 
      write(*,'(a,2f10.2)') 'EL1 =',EL1, EK*Ry-EL1-E_ion
      Allocate(dr1(mch  ,nps)); dr1=0.d0
      Allocate(di1(mch  ,nps)); di1=0.d0
      Do it=1,nps
       read(nut,*) is,nop
       read(nut,*) dr1(1:nop,it)
       read(nut,*) di1(1:nop,it)
      End do
      di1 = -di1         ! complex conjugate,     ??? 

! ... Coulomb phases: cp(l+1) = s(l+1)-s(0) = cp(l) + tn^-1(q/(l+1)):

      llmax = maxval(lch)

      Allocate(CP1(0:llmax,ntarg)); CP1 = 0.d0
      Do it=1,ntarg
       E2 = EL1/Ry/2 + (etarg(1) - etarg(it)) ! total energy of ionized electron
       if(E2.lt.0.d0) Cycle
       k2 = sqrt(2*E2) 
       q = -1.d0/k2; CP1(0,1)=0.d0
       Do l=1,llmax; s=l; CP1(l,it)=CP1(l-1,it)+DATAN2(q,s); End do
      End do

! ... define projected amplitudes summed over pseudostates for the same symmetry (L1,S1)

      Allocate(FI(2,-lmax:lmax,mch,nlsp,0:ng)); FI = 0.d0

      ET = Eps(1) + EK/2  
      E1 = ET - EL1/Ry/2 - etarg(1)
      k1 = sqrt(2*E1)

      Do ip = iion+1,nps; ilsp = iipar_ion (lps(ip),Sps(ip),Pps(ip))

       s1 = EK/2 - (Eps(ip)-Eps(1)); s1 = sqrt(2*s1)
       s  = 1.d0
       if(icorr.gt.0) s  = sqrt(k1/s1)                 !  ???

       Do ich=1,nch(ilsp)
        or = dr1(ich,ip) * S
        oi = di1(ich,ip) * S

        Do m=-lps(ip),lps(ip)
         FI(1,m,ich,ilsp,:)=FI(1,m,ich,ilsp,:) + F(1,m,ip,:)*or - F(2,m,ip,:)*oi 
         FI(2,m,ich,ilsp,:)=FI(2,m,ich,ilsp,:) + F(1,m,ip,:)*oi + F(2,m,ip,:)*or 
        End do

       End do ! ich
      End do  ! ip  -> p.s.
      Deallocate (dr1,di1)


!----------------------------------------------------------------------
! ... projection coefficients (exchange ionization):

      Call Read_rpar(nut,'EL1',EL2) 
      write(*,'(a,2f10.2)') 'EL2 =',EL2,EK*Ry-EL2-E_ion
      Allocate(dr2(mch  ,nps)); dr2=0.d0
      Allocate(di2(mch  ,nps)); di2=0.d0
      Do it=1,nps
       read(nut,*) is,nop
       read(nut,*) dr2(1:nop,it)
       read(nut,*) di2(1:nop,it)
      End do
      di2 = -di2         ! complex conjugate,     ??? 

! ... Coulomb phases: cp(l+1) = s(l+1)-s(0) = cp(l) + tn^-1(q/(l+1)):

      llmax = maxval(lch)
      Allocate(CP2(0:llmax,ntarg)); CP2 = 0.d0
      Do it=1,ntarg
       E2 = EL2/Ry/2 + etarg(1) - etarg(it)
       if(E2.lt.0.d0) Cycle
       k2 = sqrt(2*E2) 
       q = -1.d0/k2
       Do l=1,llmax; s=l; CP2(l,it)=CP2(l-1,it)+DATAN2(q,s); End do
      End do

! ... define projected amplitudes summed for the same symmetry for given theta:

      Allocate(FJ(2,-lmax:lmax,mch,nlsp,0:ng)); FJ = 0.d0

      ET = Eps(1) + EK/2  
      E1 = ET - EL2/Ry/2 - etarg(1)
      k1 = sqrt(2*E1)

      Do ip = iion+1,nps; ilsp = iipar_ion (lps(ip),Sps(ip),Pps(ip))

       s1 = EK/2 - (Eps(ip)-Eps(1)); s1 = sqrt(2*s1)
       s  = 1.d0
       if(icorr.gt.0) s  = sqrt(k1/s1)                 !  ???

       Do ich=1,nch(ilsp)
        or = dr2(ich,ip) * S
        oi = di2(ich,ip) * S
        Do m=-lps(ip),lps(ip)
         FJ(1,m,ich,ilsp,:)=FJ(1,m,ich,ilsp,:) + F(1,m,ip,:)*or - F(2,m,ip,:)*oi 
         FJ(2,m,ich,ilsp,:)=FJ(2,m,ich,ilsp,:) + F(1,m,ip,:)*oi + F(2,m,ip,:)*or 
        End do
       End do ! ich

      End do  ! ip  -> p.s.

      Deallocate (dr2,di2)

      write(*,'(a,f10.3)') 'electrons energy:',EL1+EL2

!----------------------------------------------------------------------
! ... normalized Legendre functions:
!----------------------------------------------------------------------

      llmax = maxval(lch)
      Allocate(plmn(0:ng,0:llmax,-llmax:llmax)); plmp=0.d0

      Do ig=0,ng; y=pi*ig/180; x=COS(y) 
       Do l=0,llmax; Do m=-llmax,llmax
        plmn(ig,l,m) = ALEGFM (l,m,X,1)
       End do; End do
      End do

!----------------------------------------------------------------------
! ... calculations:
!----------------------------------------------------------------------

      Deallocate(f)

      Allocate(DCS (ntarg));  DCS=0.d0

      LT_max=maxval(ltarg)
      Allocate(FSM1(2,-LT_max:LT_max,smax,ntarg)); FSM1 = 0.d0
      Allocate(FSM2(2,-LT_max:LT_max,smax,ntarg)); FSM2 = 0.d0

      Select case(igeom)

      Case(1)           !  xz or scattering plane

       Allocate(tdcs(-ng:ng,ntarg))
       Do ig2=-ng,ng
        FSM1 = 0.d0; if(iexch.ne.2)  Call TDCS1(ig1,ig2,ifi)
        FSM2 = 0.d0; if(iexch.gt.0)  Call TDCS2(ig1,ig2,ifi)
        Call Sum_ion_ampl(mode)
        tdcs(ig2,:)=DCS(:)
       End do

      Case(2)           !  yz or perpendicular plane

       Allocate(tdcs(-ng:ng,ntarg))
       ifi = 90
       Do ig2=-ng,ng
        FSM1 = 0.d0; if(iexch.ne.2)  Call TDCS1(ig1,ig2,ifi)
        FSM2 = 0.d0; if(iexch.gt.0)  Call TDCS2(ig1,ig2,ifi)
        Call Sum_ion_ampl(mode)
        tdcs(ig2,:)=DCS(:)
       End do

      Case(3)           !  xy or full-perpendicular plane

       Allocate(tdcs(0:2*ng,ntarg))
       ig2=90
       Do ifi=0,2*ng
        FSM1 = 0.d0; if(iexch.ne.2)  Call TDCS1(ig1,ig2,ifi)
        FSM2 = 0.d0; if(iexch.gt.0)  Call TDCS2(ig1,ig2,ifi)
        Call Sum_ion_ampl(mode)
        tdcs(ifi,:)=DCS(:)
       End do

      Case(4)           !  perpendicular plane

       Allocate(tdcs(0:2*ng,ntarg))
       ig1=90; ig2=90
       Do ifi=0,2*ng
        FSM1 = 0.d0; if(iexch.ne.2)  Call TDCS1(ig1,ig2,ifi)
        FSM2 = 0.d0; if(iexch.gt.0)  Call TDCS2(ig1,ig2,ifi)
        Call Sum_ion_ampl(mode)
        tdcs(ifi,:)=DCS(:)
       End do

      Case(5)           !  symmetric theta1 and theta2

       Allocate(tdcs(-ng:ng,ntarg))
       Do ig2=0,ng;  ig1=-ig2; ifi=0
        FSM1 = 0.d0; if(iexch.ne.2)  Call TDCS1(ig1,ig2,ifi)
        FSM2 = 0.d0; if(iexch.gt.0)  Call TDCS2(ig1,ig2,ifi)
        Call Sum_ion_ampl(mode)
        tdcs(ig2,:)=DCS(:)
       End do

      Case(6)           !  fixed theta1-theta2

       Allocate(tdcs(-ng:ng,ntarg))
       Do ig2=-ng,ng;  ig1 = ig2-idel
                       if(ig1.lt.-180) ig1=ig1+360
                       if(ig1.gt. 180) ig1=ig1-360
        FSM1 = 0.d0; if(iexch.ne.2)  Call TDCS1(ig1,ig2,ifi)
        FSM2 = 0.d0; if(iexch.gt.0)  Call TDCS2(ig1,ig2,ifi)
        Call Sum_ion_ampl(mode)
        tdcs(ig2,:)=DCS(:)
       End do          

      Case(10)          !  full 3D

      Allocate(tdcs(-ng:ng,0:360))

      Do ifi = 0,360,idfi
       Do ig2= 0,ng,idte
        FSM1 = 0.d0; if(iexch.ne.2)  Call TDCS1(ig1,ig2,ifi)
        FSM2 = 0.d0; if(iexch.gt.0)  Call TDCS2(ig1,ig2,ifi)
        Call Sum_ion_ampl(mode)
        tdcs(ig2,ifi)=DCS(itarg)
       End do
      End do

      End Select

!----------------------------------------------------------------------
! ... normalization:
!----------------------------------------------------------------------

      tdcs = tdcs /  (2*Pi*EK)     !  2*pi from spherical harmonics

! ... units:

      write(AU,'(a,i2,a)') '   in  a.u.'       
      SN = 1.d0
      if(i16.gt.0) then
       SN = SN * 0.529177**2  * 10**(i16-16) / (2*Ry)
       write(AU,'(a,i2,a)') '  units:  10^-',i16,'cm^2 sr^-2 eV^-1'       
       tdcs = tdcs * SN
      end if      

!----------------------------------------------------------------------
! ... output:
!----------------------------------------------------------------------

      Call Read_aarg('out',AF_td)
      if(igeom.ne.10) open(out,file=AF_td)

      iout = 1;   Call Read_iarg('iout',iout) 

      Select case(igeom*100+iout)

       Case(101)

        write(out,'(1x,a,4x,a,5x,a)') 'theta2','xz-plane',trim(AU)

        Do ig=0,180;  write(out,'(i5,E15.5)') ig,tdcs(ig,itarg); End do
        Do ig=-179,0; write(out,'(i5,E15.5)') ig+360,tdcs(ig,itarg); End do

       Case(102)

        write(out,'(1x,a,4x,a,5x,a)') 'theta2','xz-plane',trim(AU)

        Do ig=-180,180;  write(out,'(i5,E15.5)') ig,tdcs(ig,itarg); End do

       Case(201)

        write(out,'(1x,a,4x,a,5x,a)') 'theta2','yz-plane',trim(AU)

        Do ig=0,180;  write(out,'(i5,E15.5)') ig,tdcs(ig,itarg); End do
        Do ig=-179,0; write(out,'(i5,E15.5)') ig+360,tdcs(ig,itarg); End do

       Case(202)

        write(out,'(1x,a,4x,a,5x,a)') 'theta2','yz-plane',trim(AU)

        Do ig=-180,180;  write(out,'(i5,E15.5)') ig,tdcs(ig,itarg); End do

       Case(301)

        write(out,'(3x,a,5x,a,5x,a)') 'fi2','xy-plane',trim(AU)

        Do ig=0,2*ng;  write(out,'(i5,E15.5)') ig,tdcs(ig,itarg); End do

       Case(302)

        write(out,'(3x,a,5x,a,5x,a)') 'fi2','xy-plane',trim(AU)

        Do ig=0,179;  write(out,'(i5,E15.5)') ig,tdcs(2*ng-ig,itarg); End do
        Do ig=0,ng;   write(out,'(i5,E15.5)') ig,tdcs(ig,itarg); End do

       Case(401)

        write(out,'(3x,a,5x,a,5x,a)') 'fi2','xy-plane',trim(AU)

        Do ig=0,2*ng;  write(out,'(i5,E15.5)') ig,tdcs(ig,itarg); End do

       Case(501)

        Do ig=0,ng
         write(out,'(i5,E15.5)') ig,tdcs(ig,itarg)     
        End do

       Case(601)
        Do ig=-ng,ng
         write(out,'(i5,E15.5)') ig+ng,tdcs(ig,itarg)     
        End do

!------------------------------------------------------------
      Case(1001)

! ... all angle data (spherical coordinates):

      ia = LEN_TRIM(AF_td)
      AF_td = AF_td(1:ia)//'.theta_fi'
      open(out,file=AF_td)

      Do jg=0,ng,idte; Do ifi = 0,360,idfi 
       write(out,'(2i5,E15.5)') jg, ifi,tdcs(jg,ifi)
      End do; End do

! ... all angle data (Carthesin coordinates):

      AF_td = AF_td(1:ia)//'.xyz'
      open(out,file=AF_td)

      Do jg=0,ng,idte; do ifi = 0,360,idfi
       T=pi*jg/180; R=pi*ifi/180; 
       if(jg.eq.  0.and.ifi.ne.0) Cycle
       if(jg.eq.180.and.ifi.ne.0) Cycle
       S = tdcs(jg,ifi)
       x = S*sin(T)*cos(R)
       y = S*sin(T)*sin(R)
       z = S*cos(T)
       write(out,'(3E15.5)') x,y,z
      End do; End do

! ... sections in three perpendicular planes:

      AF_td = AF_td(1:ia)//'.sec'
      open(out,file=AF_td)
      write(out,'(a,7x,a,13x,a,13x,a,10x,a)') 'teta','xz','yz','xy',AU
      Do jg=0,ng,idte 
       write(out,'(i3,3E15.5)') jg, tdcs(jg, 0), &
                                    tdcs(jg,90), &
                                    tdcs(90,jg)
      End do

      Do jg=ng+idte,ng+ng,idte;  ig = 360-jg;  
       write(out,'(i3,3E15.5)')  jg,tdcs(ig,180), &
                                    tdcs(ig,270), &
                                    tdcs(90,jg )
      End do



!------------------------------------------------------------
       Case default
        Stop 'unknown case'
      End Select

      Call CPU_time(t2)
      write(*,'(/a,f6.2,a)')  ' time = ',(t2-t1)/60,' min'

     End  ! UTILITY  TDCS_all



!======================================================================
      Subroutine TDCS1(ig1,ig2,fi2)
!======================================================================

      Use tdcs_amplitude
      Implicit real(8) (A-H,O-Z)
      Integer, intent(in) :: ig1,ig2,fi2

      kg1 = iabs(ig1); kg2=iabs(ig2)

      sfi = pi*fi2/180; if(ig2.lt.0) sfi=sfi+pi; if(ig1.lt.0) sfi=sfi-pi

! ... cycle over momentums of ejected electron:

      Do ilsp=1,nlsp; LL=Lpar(ilsp); LS=ispar(ilsp);    ! over L1 for given S1
                                                        
      Do ML=-LL,LL                                      ! over  M1

       Do ich = 1,nch(ilsp); it=iptar(ilsp,ich)         ! over l2
                             l =lch  (ilsp,ich)
         s = 0.d0
         s = s - pi/2*l                                 ! (-i)^l               (sign ???) 
         s = s + CP1(l,it)                              ! Coulomb phase        (sign ???)

       Do m=-l,l                                        ! over m2

         ss = s +  m*sfi                                ! exp(m2*fi)           (sign ???)

         br = DCOS(ss)                                  ! we have complex congugated in comparison 
         bi = DSIN(ss)                                  ! to the notes  !!!

        LT = ltarg(it)                         
        Do MT = -LT,LT                                 !  ion-target Mf
          yc =  CLEBCH(LT,MT,l,m,LL,ML); if(yc.eq.0.d0) Cycle
          y = yc * plmn(kg2,l,m)                                              
          ar = FI(1,ML,ich,ilsp,kg1)*y; ai = FI(2,ML,ich,ilsp,kg1)*y 
          FSM1(1,MT,LS,it) = FSM1(1,MT,LS,it) + ar*br - ai*bi
          FSM1(2,MT,LS,it) = FSM1(2,MT,LS,it) + ar*bi + ai*br
        End do  ! MT

       End do; End do  ! l,m
      End do; End do   ! ilsp, ML
      

      End  Subroutine TDCS1 


!======================================================================
      Subroutine TDCS2(ig1,ig2,fi2)
!======================================================================

      Use tdcs_amplitude

      Implicit real(8) (A-H,O-Z)

      Integer, intent(in) :: ig1,ig2,fi2

      kg1 = iabs(ig1); kg2=iabs(ig2)

      sfi = pi*fi2/180; if(ig2.lt.0) sfi=sfi+pi

! ... cycle over momentums of ejected electron:

      Do ilsp=1,nlsp; LL=Lpar(ilsp); LS=ispar(ilsp)      ! over L1 for given S1
      Do ML=-LL,LL                                       ! over  M1   

       sm = -ML*sfi                                      ! exp(m*fi) for scattering electron 
                                                         ! M1=-m1
      Do ich = 1,nch(ilsp); it = iptar(ilsp,ich)         ! over l2
                             l = lch  (ilsp,ich)
       s = sm
       s = s - pi/2*l                                    ! (-i)^l               (sign ???)
       s = s + CP2(l,it)                                 ! Coulomb phase        (sign ???)
       br = DCOS(S);  bi = DSIN(S)

      Do m=-l,l                                          ! over m2

       LT = ltarg(it)
       Do MT = -LT,LT                                    !  ion-target Mf

        y =  CLEBCH(LT,MT,l,m,LL,ML)
        if(y.eq.0.d0) Cycle
        y = y*plmn(kg1,l,m); if(ig1.lt.0) y=y*(-1)**m

        ar = FJ(1,ML,ich,ilsp,kg2)*y; ai = FJ(2,ML,ich,ilsp,kg2)*y 
        FSM2(1,MT,LS,it) = FSM2(1,MT,LS,it) + ar*br - ai*bi
        FSM2(2,MT,LS,it) = FSM2(2,MT,LS,it) + ar*bi + ai*br

       End do  ! MT

       End do; End do  ! lm
       End do; End do  ! ilsp

     End  Subroutine TDCS2


!======================================================================
      Subroutine Sum_ion_ampl(mode)
!======================================================================

      Use tdcs_amplitude

      Implicit real(8) (A-H,O-Z)

      Select case(mode)
!----------------------------------------------------------------------
      Case(0)        ! incoherent adding

      Do it=1,ntarg
       DCS(it)= SUM(FSM1(:,:,:,it)*FSM1(:,:,:,it)) +   &
                SUM(FSM2(:,:,:,it)*FSM2(:,:,:,it))
      End do
       
!----------------------------------------------------------------------
      Case(1)        ! coherent adding

      FSM1 = FSM1 + FSM2
      Do it=1,ntarg
       DCS(it)= SUM(FSM1(:,:,:,it)*FSM1(:,:,:,it)) 
      End do

!----------------------------------------------------------------------
      Case(2)        ! spin-dependent adding

      c1 = 1.d0/2.d0;  c2 = sqrt(3.d0)/2.d0

      FSM1(:,:,1,:) = FSM1(:,:,1,:) - c1* FSM2(:,:,1,:) &
                                    - c2* FSM2(:,:,3,:)
      FSM1(:,:,3,:) = FSM1(:,:,3,:) - c2* FSM2(:,:,1,:) &
                                    + c1* FSM2(:,:,3,:)
      Do it=1,ntarg
       DCS(it)= SUM(FSM1(:,:,:,it)*FSM1(:,:,:,it)) 
      End do

!----------------------------------------------------------------------

      Case default; Stop 'unknown mode in sum_ion_ampl'

      End Select

      End Subroutine Sum_ion_ampl


!======================================================================
      Subroutine inf_tdcs
!======================================================================
!     provide screen information about add_farm utility
!----------------------------------------------------------------------
      Implicit real(8) (a-h,o-z)

      Character(80) :: A

      Call get_command_argument(1,A)  
      if(A.ne.'?') Return

      write(*,*) &
'                                                                  ',&
'tdcs calculates triple differential cross section for ionization  ',&
'of He or other closed-shell atoms, using the projection method    ',&
'                                                                  ',&
'Input parameters and files (provided as par=... ):                ',&
'                                                                  ',&
'ion    [target_ion]  - description of (e-A+) scattering           ',&
'dcs    [difsec_He ]  - collection of scattering amplitudes for    ',&
'                       (e+A), for each pseudosctate at given E0   ',&
'p1     [projection2] - projection coefficients for each psedostate',&
'                       to the continuum (e-A+) at E = E_ejected   ',&
'p2     [projection1] - projection coefficients for each psedostate',&
'                       to the continuum (e-A+) at E = E_scattering',&
'ig1    [10]          - theta_1 (in degrees)                       ',&
'i16    [0]           - units for TDCS, 0 -> a.u.                  ',&
'mode   [2]           - mode for adding of direct and exchange     ',&
'                       0 - incoherent (debug option)              ',&
'                       1 - coherent   (debug option)              ',&
'                       2 - spin-dependent                         ',&
'df     [10]          - step in fi angle  (in output .all)         ',&
'dt     [10]          - step in teta angle (in output .all)        ',&
'itagr  [1]           - ion target state                           ',&
'                                                                  ',&
'geom [1]             - geometry mode                              ',&
'iout [1]             - output mode                                ',&
'                                                                  ',&
'                                                                  '
    Stop ' '
                         
    End Subroutine inf_tdcs
                                         