!=====================================================================
!     UTILITY  DDCS
!=====================================================================
!
!     double differential cross section for ionization of He or other
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

      Use target_ion             ! we can work with usual target ???
      Use channels_ion

      Implicit real(8) (A-H,O-Z)

      Character(200) :: AS,AU

      Real(8), Allocatable :: g(:)                 ! teta grid:  -180 : 180
      Real(8), Allocatable :: f(:,:,:,:)           ! amplidutes for p.s. 

      Real(8), Allocatable :: fi(:,:,:,:,:)        ! projected amplitudes (direct)
      Real(8), Allocatable :: fj(:,:,:,:,:)        ! projected amplitudes (exchange)

      Real(8), Allocatable :: FSM1(:,:,:,:,:)      ! direct  ionization amplitudes 
      Real(8), Allocatable :: FSM2(:,:,:,:,:)      ! exchage ionization amplitudes 

      Real(8), Allocatable :: dcs (:,:)            ! tdcs for given fi
      Real(8), Allocatable :: tdcs(:,:,:)          ! complete tdcs

      Real(8), Allocatable :: sec(:), dd(:,:)      ! auxiliary

      Real(8), Allocatable :: ddcs1(:),ddcs2(:)    ! double dcs

      Real(8), Allocatable :: Eps(:)               ! energies of p.s.
      Integer, Allocatable :: Lps(:),Sps(:),Pps(:) ! terms of p.s.

      ! ... projections:

      Real(8), Allocatable :: wt1(:,:), di1(:,:), dr1(:,:)
      Real(8), Allocatable :: wt2(:,:), di2(:,:), dr2(:,:)

      ! ... Coulomb phases:

      Real(8), Allocatable :: CP1(:,:), CP2(:,:)

      ! ... normalized Legendre functions:

      Real(8), Allocatable :: plmn(:,:,:)

      ! ... files involved:

      Integer :: nut=1;  Character(40) :: AF_i   = 'target_ion'
                         Character(40) :: AF_p   = 'projection2'
                         Character(40) :: AF_dcs = 'difsec_He'

      Integer :: out=3;  Character(40) :: AF_dd  = 'ddcs.out'

      Real(8) :: Ry = 13.6056
      Real(8) :: hg = 1
      Real(8) :: k0,k1,k2         ! electrons momentums
      Integer :: lmax, smax       ! max. L and 2S+1 for p.s. 

      Call inf_ddcs
       
      Call CPU_time(t1)

      pi  = abs(ACOS(-1.d0))
!----------------------------------------------------------------------
! ... define input parameters:

      i16  = 18;      Call Read_iarg('i16'  ,i16)
      mode =  2;      Call Read_iarg('mode' ,mode)
      itarg=  1;      Call Read_iarg('itarg',itarg)
      icorr=  1;      Call Read_iarg('corr' ,icorr) 
      iexch=  1;      Call Read_iarg('exch' ,iexch) 

      write(*,'(a,i5)') 'mode  = ',mode
      write(*,'(a,i5)') 'iexch = ',iexch
      write(*,'(a,i5)') 'itarg = ',itarg
      write(*,'(a,i5)') 'i16   = ',i16
      write(*,'(a,i5)') 'icorr = ',icorr

!----------------------------------------------------------------------
! ... ion target information:

      Call Read_aarg('ion',AF_i)
      Call Check_file(AF_i)
      Open(nut,file=AF_i)
      Call R_target_ion(nut)
      Call R_channels_ion(nut)
      Close(nut)

       Z = nz;  AWT = 0.d0;  ZC = Z - nelc       
!      Call Conv_au (Z,AWT,au_cm,au_eV,0)
!      Ry = au_eV/2.d0                           !  ???
!      write(*,'(a,f10.5)') 'Ry = ', Ry

!----------------------------------------------------------------------
! ... read amplidudes for all pseudostates at given angle: 
! ... we supposed universal angle grid: 0,1,2,3,...180

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
      if(abs(EK-EKK).gt.0.001) Stop 'different EK in difsec and projection files'

      Call Read_ipar(nut,'nps',npp) 
      if(npp.ne.nps) Stop 'different nps in difsec and projection files' 

      Call Read_rpar(nut,'EL2',EL1) 
      write(*,'(a,2f10.2)') 'EL1 =',EL1, EK*Ry-EL1-E_ion
      Allocate(dr1(mch  ,nps)); dr1=0.d0
      Allocate(di1(mch  ,nps)); di1=0.d0
      Allocate(wt1(ntarg,nps)); wt1=0.d0
      Do it=1,nps
       read(nut,*) is,nop
       read(nut,*) dr1(1:nop,it)
       read(nut,*) di1(1:nop,it)
!       read(nut,*) wt1(:,it)
      End do
      di1 = -di1         ! complex conjugate,     ??? 

! ... Coulomb phases: cp(l+1) = s(l+1)-s(0) = cp(l) + tn^-1(q/(l+1)):

      llmax = maxval(lch)
      Allocate(CP1(0:llmax,ntarg)); CP1 = 0.d0
      Do it=1,ntarg
       E2 = EL1/Ry/2 + (etarg(1) - etarg(it)) ! total energy of ionized electron
       if(E2.lt.0.d0) Cycle
       k2 = sqrt(2*E2) 
       q = -1.d0/k2
       Do l=1,llmax; s=l; CP1(l,it)=CP1(l-1,it)+DATAN2(q,s); End do
      End do

! ... define projected amplitudes summed for the same symmetry for given teta:

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
      Deallocate (dr1,di1,wt1)


!----------------------------------------------------------------------
! ... projection coefficients (exchange ionization):

      Call Read_rpar(nut,'EL1',EL2) 
      write(*,'(a,2f10.2)') 'EL2 =',EL2,EK*Ry-EL2-E_ion
      Allocate(dr2(mch  ,nps)); dr2=0.d0
      Allocate(di2(mch  ,nps)); di2=0.d0
      Allocate(wt2(ntarg,nps)); wt2=0.d0
      Do it=1,nps
       read(nut,*) is,nop
       read(nut,*) dr2(1:nop,it)
       read(nut,*) di2(1:nop,it)
!       read(nut,*) wt2(:,it)
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

! ... define projected amplitudes summed for the same symmetry for given teta:

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

      Deallocate (dr2,di2,wt2)

      write(*,'(a,f10.3)') 'electrons energy:',EL1+EL2

!----------------------------------------------------------------------
! ... normalized Legendre functions:

      llmax = maxval(lch)
      Allocate(plmn(0:ng,0:llmax,-llmax:llmax))

      Do ig=0,ng; x=DCOS(DFLOAT(ig)/180*pi) 
       Do l=0,llmax; Do m=-llmax,llmax
        plmn(ig,l,m) = ALEGFM (l,m,X,1)
       End do; End do
      End do

!----------------------------------------------------------------------
! ... calculations:
!----------------------------------------------------------------------

      Deallocate(f)

      Allocate(DCS (0:ng,ntarg))
      Allocate(tdcs(0:ng,0:2*ng,ntarg))

      LT=maxval(ltarg)
      Allocate(FSM1(2,-LT:LT,smax,ntarg,0:ng))
      Allocate(FSM2(2,-LT:LT,smax,ntarg,0:ng))

      Allocate(sec(0:2*ng), dd(0:ng,0:ng))

      Do ig1 = 0,ng

      Do ifi = 0,2*ng;  sfi = pi*ifi/180.d0
       FSM1 = 0.d0; if(iexch.ne.2)  Call TDCS1
       FSM2 = 0.d0; if(iexch.gt.0)  Call TDCS2
       Call Sum_ion_ampl(mode)
       Do it=1,ntarg;  tdcs(:,ifi,it)=DCS(:,it);  End do
      End do

      Do ig2 = 0,ng
       sec(:) = tdcs(ig2,:,itarg) 
       dd(ig1,ig2) = ZIS(hg,2*ng+1,sec) * (pi/ng)  
      End do

      End do ! over ig1

      Deallocate(FI, FJ, FSM1, FSM2, tdcs)

!----------------------------------------------------------------------
! ... normalization:

      dd = dd / EK  / (2*pi)  

! ... units:

      write(AU,'(a,i2,a)') 'in  a.u.'       
      SN = 1.d0
      if(i16.gt.0) then
       SN = SN * 0.529177**2  * 10**(i16-16) / (2*Ry)
       write(AU,'(a,i2,a)') '  units:  10^-',i16,'cm^2 sr^-2 eV^-1'       
       dd = dd * SN
      end if      

!----------------------------------------------------------------------
! ... integration over one teta:

      Allocate (ddcs1(0:ng),ddcs2(0:ng)) 

      Do ig1=0,ng
       Do ig2=0,ng; sec(ig2) = dd(ig1,ig2)*DSIN(DFLOAT(ig2)/180*pi); End do
       ddcs1(ig1) = ZIS(hg,ng+1,sec) * (pi/ng)  
      End do

      Do ig2=0,ng
       Do ig1=0,ng; sec(ig1) = dd(ig1,ig2)*DSIN(DFLOAT(ig1)/180*pi); End do
       ddcs2(ig2) = ZIS(hg,ng+1,sec) * (pi/ng)  
      End do

! ... single differential cross sections at EL1 and EL2

      Do ig1=0,ng; sec(ig1) = ddcs1(ig1)*DSIN(DFLOAT(ig1)/180*pi); End do
      s1 = ZIS(hg,ng+1,sec) * (pi/ng)  * 2*pi ! from integration over fi 
     
      Do ig2=0,ng; sec(ig2) = ddcs2(ig2)*DSIN(DFLOAT(ig2)/180*pi); End do
      s2 = ZIS(hg,ng+1,sec) * (pi/ng)  * 2*pi

!----------------------------------------------------------------------
! ... output:

      IE = EK*Ry + 0.2
      JE = EL1 + 0.2

      write(AF_dd,'(a,i3.3,a,i3.3,a,i1,a)') &
                   'ddcs_',IE,'eV_',JE,'eV_t',itarg,'.sec'
      Call Read_aarg('out',AF_dd) 


      open(out,file=AF_dd)

      write(out,'(a,10x,a,12x,a,12x,a,a)') &
                    ' teta','dcs1','dcs2','dcs',AU
      Do ig=0,ng 
       write(out,'(i5,2E15.5)') ig,ddcs1(ig),ddcs2(ig)
      End do

      write(out,'(//a5,3F10.3)') '#EL', EL1,EL2,EL1+EL2
      write(out,'(//a5,3E15.5)') '#sdcs', s1,s2
    
      Call CPU_time(t2)
      write(*,'(/a,f6.2,a)')  ' time = ',(t2-t1)/60,' min'


!======================================================================
CONTAINS
!======================================================================


!======================================================================
      Subroutine TDCS1
!======================================================================

      Implicit real(8) (A-H,O-Z)

! ... cycle over angles of ejected electron:

      Do ilsp=1,nlsp; LL=Lpar(ilsp); LS=ispar(ilsp);
      Do ML=-LL,LL             !  p.s. M
        sign = 1.d0;  if(ig1.lt.0) sign = (-1)**ML        ! ig1 ???

       Do ich = 1,nch(ilsp); it=iptar(ilsp,ich)
                             l =lch  (ilsp,ich)
         s = 0.d0
         s = s - pi/2*l         !  (-i)^l
         s = s + CP1(l,it)      !  Coulomb phase

       Do m=-l,l                !  is it m for electron ???

         ss = s +  m*sfi        !  exp(m*fi) 
         br = DCOS(ss)
         bi = DSIN(ss)

         LT = ltarg(it)                         
        Do MT = -LT,LT         !  ion-target M

          yc =  CLEBCH(LT,MT,l,m,LL,ML); if(yc.eq.0.d0) Cycle

        Do ig=0,ng 
          y = yc * sign * plmn(ig,l,m)          
          ar = FI(1,ML,ich,ilsp,ig1)*y; ai = FI(2,ML,ich,ilsp,ig1)*y 
          FSM1(1,MT,LS,it,ig) = FSM1(1,MT,LS,it,ig) + ar*br - ai*bi
          FSM1(2,MT,LS,it,ig) = FSM1(2,MT,LS,it,ig) + ar*bi + ai*br
        End do   ! ig


        End do  ! MT

       End do; End do  ! l,m
      End do; End do  ! ilsp, ML
      

      End  Subroutine TDCS1 


!======================================================================
      Subroutine TDCS2
!======================================================================

      Implicit real(8) (A-H,O-Z)

! ... we have one angle for ejected electron:

      Do ilsp=1,nlsp; LL=Lpar(ilsp); LS=ispar(ilsp)
      Do ML=-LL,LL    ! it is M of p.s., small m is opposite  

       sm = -ML*sfi   !  exp(m*fi) for scattering electron 

      Do ich = 1,nch(ilsp); it = iptar(ilsp,ich)
                             l = lch  (ilsp,ich)
       s = sm
       s = s - pi/2*l  
       s = s + CP2(l,it) 
       br = DCOS(S);  bi = DSIN(S)

      Do m=-l,l

       LT = ltarg(it)
       Do MT = -LT,LT   !  ion-target M

        y =  CLEBCH(LT,MT,l,m,LL,ML)
        if(y.eq.0.d0) Cycle
        y = y*plmn(ig1,l,m)

        Do i=0,ng
         ar = FJ(1,ML,ich,ilsp,i)*y; ai = FJ(2,ML,ich,ilsp,i)*y 
         FSM2(1,MT,LS,it,i) = FSM2(1,MT,LS,it,i) + ar*br - ai*bi
         FSM2(2,MT,LS,it,i) = FSM2(2,MT,LS,it,i) + ar*bi + ai*br
        End do

       End do  ! MT

       End do; End do  ! lm
       End do; End do  ! ilsp

     End  Subroutine TDCS2


!======================================================================
      Subroutine Sum_ion_ampl(mode)
!======================================================================

      Implicit real(8) (A-H,O-Z)

      Select case(mode)
!----------------------------------------------------------------------
      Case(0)        ! incoherent adding

      Do it=1,ntarg; Do i=0,ng
       DCS(i,it)= SUM(FSM1(:,:,:,it,i)*FSM1(:,:,:,it,i)) +   &
                  SUM(FSM2(:,:,:,it,i)*FSM2(:,:,:,it,i))
      End do; End do
       
!----------------------------------------------------------------------
      Case(1)        ! coherent adding

      FSM1 = FSM1 + FSM2
      Do it=1,ntarg; Do i=0,ng
       DCS(i,it)= SUM(FSM1(:,:,:,it,i)*FSM1(:,:,:,it,i)) 
      End do; End do

!----------------------------------------------------------------------
      Case(2)        ! spin-dependent adding

      c1 = 1.d0/2.d0;  c2 = sqrt(3.d0)/2.d0

      FSM1(:,:,1,:,:) = FSM1(:,:,1,:,:) - c1* FSM2(:,:,1,:,:) &
                                        - c2* FSM2(:,:,3,:,:)
      FSM1(:,:,3,:,:) = FSM1(:,:,3,:,:) - c2* FSM2(:,:,1,:,:) &
                                        + c1* FSM2(:,:,3,:,:)
      Do it=1,ntarg; Do i=0,ng
       DCS(i,it)= SUM(FSM1(:,:,:,it,i)*FSM1(:,:,:,it,i)) 
      End do; End do

!----------------------------------------------------------------------

      Case default; Stop 'unknown mode in sum_ion_ampl'

      End Select

      End Subroutine Sum_ion_ampl


     End  ! UTILITY  TDCS_all


!======================================================================
      Subroutine inf_ddcs
!======================================================================
!     provide screen information about add_farm utility
!----------------------------------------------------------------------
       
      Implicit real(8) (a-h,o-z)

      Character(80) :: A

      iarg = IARGC()
      if(iarg.gt.0) Call GETARG(1,A)      
      if(iarg.eq.0) A='?'
      if(A.ne.'?') Return

      write(*,*) &
'                                                                  ',&
'tdcs calculates triple differential cross section for ionization  ',&
'of He or other closed-shell atoms, using the projection method    ',&
'                                                                  ',&
'Input parameters and files (provided as  par=... ):               ',&
'                                                                  ',&
'ion    [target_ion]  - description of (e-A+) scattering           ',&
'dcs    [difsec_He ]  - collection of scattering amplitudes for    ',&
'                       (e+A), for each pseudosctate at given E0   ',&
'p1     [projection2] - projection coefficients for each psedostate',&
'                       to the continuum (e-A+) at E = E_ejected   ',&
'p2     [projection1] - projection coefficients for each psedostate',&
'                       to the continuum (e-A+) at E = E_scattering',&
'ig1    [10]          - teta_1 (in degrees)                        ',&
'i16    [0]           - units for TDCS, 0 -> a.u.                  ',&
'mode   [2]           - mode for adding of direct and exchange     ',&
'                       0 - incoherent                             ',&
'                       1 - coherent                               ',&
'                       2 - spin-dependent                         ',&
'df     [10]          - step in fi angle  (in output .all)         ',&
'dt     [10]          - step in teta angle (in output .all)        ',&
'itagr  [1]           - ion target state                           ',&
'icorr  [1]           - flag for correlation k1/k1_p_s             ',&
'                                                                  ',&
'Output files:                                                     ',&
'                                                                  ',&
'tdcs_E0_E2_teta_targ.all  (total output: fi, teta, tdcs)          ',&
'tdcs_E0_E2_teta_targ.xyz  (total output as surface in x,y,z       ',&  
'tdcs_E0_E2_teta_targ.sec  (cross sections in zx, zy, xy planes)   ',&
'                                                                  '
    Stop 
                         
    End Subroutine inf_ddcs
                                         