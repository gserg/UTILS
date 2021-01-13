!======================================================================
!    Differential cross section (case of neutral atoms in JK scheme)
!======================================================================
!
!    Input:       zarm.tma     - T-matrix elements for all LSp 
!                 zarm.tmb     - new format (more compact)
!                 target_jj    - description of channels
!
!    Output:      mt_nn_nn     - angle-integrated direct and momentum 
!                                trnsfer cross-sections 
!                 dcs_nn_nn    - DCS cross_sections
!
!    Call as:     sec_dif_JK itr1=...  itr2=... 
!                            [i16=... ifano=...]
!                            [ek1=...  ek2=... ]
!                            [gmin=... gmax=... gstep=...]
!======================================================================

      Use tmat_list
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)

      Integer :: nut= 1;  Character(80) :: AF_tar='target'
      Integer :: num= 2;  Character(80) :: AF_tm ='zarm.tma'
                          Character(80) :: AF_tmb='zarm.tmb'
      Integer :: nus=10;  Character(80) :: AF_td ='tmat.done'
                          Character(80) :: AF_inp='tmat.done_inp'
                          Character(80) :: AF_out='tmat.done_out'
      Integer :: nur=11;  Character(80) :: AF_sec='mt_nn_nn'
      Integer :: nud=12;  Character(80) :: AF_dcs='dcs_nn_nn'
      Integer :: pri= 6;  Character(80) :: AF_log='sec_dif_JK.log'

      Real(8), Allocatable :: ev(:),TRmat(:),TImat(:), G(:)
      Real(8), allocatable :: DCS(:), DS(:), DM(:), CP(:,:)
      Integer, Allocatable :: ipe(:)
      Real(8), allocatable :: FR(:,:,:,:,:), FI(:,:,:,:,:)
    
      Character(200) :: AS 
      Real(8) :: Ry = 13.6056,  a0 = 0.529177
      PI  = abs(ACOS(-1.d0))

      Call Inf_sub  

!----------------------------------------------------------------------
! ... input parameters:

      ifano=0;      Call Read_iarg('ifano',ifano)
      i16=16;       Call Read_iarg('i16'  ,i16  )
      idcs = 1;     Call Read_iarg('dcs'  ,idcs )
      itdone = 0;   Call Read_iarg('tdone',itdone)
      JJ_top=0;     Call Read_iarg('JJ_top',JJ_top)

      Glow  = 0;    Call Read_rarg('Glow' ,Glow )
      Ghigh = 180;  Call Read_rarg('Ghigh',Ghigh)
      Gstep = 1;    Call Read_rarg('Gstep',Gstep)
      ng = NINT( (Ghigh-Glow)/Gstep ) + 1
      hg = Gstep /180.d0 * pi
      Allocate(G(ng),DCS(ng),DS(ng),DM(ng))
      Do i=1,ng; G(i) = Glow + (i-1)*Gstep; End do

      open(pri,file=AF_log)
      write(pri,'(a,i5)') 'ng    = ',ng
      write(pri,'(a,i5)') 'ifano = ',ifano
      write(pri,'(a,i5)') 'i16   = ',i16
      write(pri,'(a,i5)') 'dcs   = ',idcs
      write(pri,'(a,i5)') 'tdone = ',itdone

      if(JJ_top.gt.0) Call T_extend(JJ_top,AF_inp,AF_out)

!----------------------------------------------------------------------
! ... target:

      Call Check_file(AF_tar)
      Open(nut,file=AF_tar)
      Call R_target(nut)
      Call R_channels(nut)
      np = ntarg; Call Read_ipar(nut,'np',np)
      ni = ntarg; Call Read_ipar(nut,'ni',ni)
      Close(nut)

      e1=etarg(1); Do i=1,ntarg; etarg(i)=(etarg(i)-e1)*2; End do
      ion = nz - nelc; zion=ion*ion; if(ion.eq.0) zion = 1.d0
      if(ion.ne.0) Stop 'ion <> 0, not a neutral atom'

!----------------------------------------------------------------------
! ... check transition:

      itr1=1;   Call Read_iarg('itr1' ,itr1 )
      itr2=1;   Call Read_iarg('itr2' ,itr2 )

      if(itr1.lt.1)     itr1=1
      if(itr1.gt.ntarg) itr1=ntarg
      if(itr1.gt.np)    Stop 'itr1 > np, # physical states'
      if(itr2.lt.itr1)  itr2=itr1
      if(itr2.gt.ntarg) itr2=ntarg

      JT1=jtarg (itr1); JT2=jtarg (itr2)
      ek1 = 0.d0;  Call Read_rarg('ek1' ,ek1)
      ek2 = 0.d0;  Call Read_rarg('ek2' ,ek2)
      if(ek2.lt.ek1) ek2=ek1

      ekk = 0.d0;  Call Read_rarg('ekk' ,ekk)
      if(ekk.ne.0.d0) ek1=ekk
      if(ekk.ne.0.d0) ek2=ekk

      E1=etarg(itr1);  E2=etarg(itr2)

      if(itdone.ne.0) go to 10

!----------------------------------------------------------------------
! ... read T-matrix:

      if(icheck_file(AF_tmb).eq.1) then
       istyle=1; Open(num,file=AF_tmb)
      elseif(icheck_file(AF_tm).eq.1) then
       istyle=0; Open(num,file=AF_tm)
      else
       Stop 'cannot open T-matrx file'
      end if       

      mdim=mch*(mch+1)/2;  Allocate (trmat(mdim),timat(mdim))

      rewind(num)
   1  Continue

      if(istyle.eq.0) then
       read(num,*,end=2) e1,nopen,ntr,ilsp
       read(num,*) (TRmat(i),TImat(i),i=1,ntr)
      else
       read(num,*,end=2) e1,nopen,kopen,ilsp,i1,i2,nj
       if(kopen.gt.nopen)     Stop ' kopen > nopen'
       ntr = kopen*(kopen+1)/2
       read(num,*) (trmat(i),timat(i),i=1,ntr)
       ktr = (nopen-kopen)*nj
       if(ktr.gt.0) read(num,*) (trmat(i),timat(i),i=ntr+1,ntr+ktr)
      end if

      if(e1.le.etarg(itr2)) go to 1
      if(ek1.gt.0.d0.and.e1.lt.ek1) go to 1
      if(ek2.gt.0.d0.and.e1.gt.ek2) go to 1

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


        Call Add_tmat_list(e1,TRMAT(itr),TIMAT(itr),Jpar(ilsp)-1,Ipar(ilsp), &
                            lch(ilsp,ich1),lch(ilsp,ich2), jkch(ilsp,ich1)-1,jkch(ilsp,ich2)-1 )
       End do
      End do

      go to 1
    2 write(*,'(/a,i5)') 'number of T-elements:',ndata
      if(ndata.eq.0) Stop 'nothing to do '
    
      Call Sort_tmat_list(0)

      ne=1;  Do i=2,ndata; if(ek(i-1).ne.ek(i)) ne=ne+1;  End do
      write(*,'( a,i5)') 'number of energies:  ',ne

!----------------------------------------------------------------------
! ... create tmat.done:      

      open(nut,file=AF_td)

      ip1 = 0; if(iptarg(itr1).eq.-1) ip1=1
      ip2 = 0; if(iptarg(itr2).eq.-1) ip2=1

      write(nut,'(6i6,2i8,2D15.7)')  &
       itr1,itr2,JT1-1,JT2-1,ip1,ip2,ndata,ne,Etarg(itr1),Etarg(itr2)

      Do i=1,ndata
        write(nut,'(D15.7,6I5,2D15.7)')  &
                   ek(i),jjpar(i),ip(i),l1(i),jk1(i),l2(i),jk2(i),TR(i),TI(i)
      End do

      Close(nut)

!============================================================================
! ... insertion for the T-Matrix symmetries:

      if(ekk.eq.0.d0) go to 5     

      Call Def_patten

      Call Sort_tmat_list(1)

! ...  record tmat.done_inp:

      open(nut,file='tmat.done_inp')
      write(nut,'(6i6,2i8,2D15.7)')  itr1,itr2,JT1-1,JT2-1,ip1,ip2,ndata,ne,etarg(itr1),etarg(itr2)

      Do i=1,ndata
       a=1.d0; aa=1.d0
       if(i.gt.1) then
         a=TR(i)/TR(i-1); if(a.lt.0.001.or.a.gt.0.999) a=0.d0  
         aa=TI(i)/TI(i-1); if(aa.lt.0.001.or.aa.gt.0.999) aa=0.d0 
       end if
       write(nut,'(D15.7,6I5,2D15.7,2f10.5)')  &
         ek(i),jjpar(i),ip(i),l1(i),jk1(i),l2(i),jk2(i),TR(i),TI(i),a,aa
      End do

      write(nut,'(/a,i3)') 'kcase =',kcase
      Do i=1,ndata
       if(jjpar(i).ne.jj_max) Cycle
       a=1.d0; aa=1.d0
       if(i.gt.1) then
         a=TR(i)/TR(i-1); if(a.lt.0.001.or.a.gt.0.999) a=0.d0  
         aa=TI(i)/TI(i-1); if(aa.lt.0.001.or.aa.gt.0.999) aa=0.d0 
       end if
       write(nut,'(D15.7,6I5,2D15.7,2f10.5,5i5)')  &
         ek(i),jjpar(i),ip(i),l1(i),jk1(i),l2(i),jk2(i),TR(i),TI(i),a,aa,jo(:,ip(i)) 
      End do

      Close(nut)

      Stop 'tmat.done_inp is written'

    5 Continue

!-------------------------------------------------------------------------------
! ... read tmat.done

   10 Continue
      Call alloc_tmat_list(-1)

      open(nut,file=AF_td)
      read(nut,*) itr1,itr2,JT1,JT2,ip1,ip2,ncase,ne,E1,E2
      if(ncase.le.0) Stop 'ncase = 0'
      if(ne.le.0) Stop 'ne = 0'

      Allocate(ev(ne),ipe(ncase)); nn=0

      Do i=1,ncase
        read(nut,'(D15.7,6I5,2D15.7)')  &
                   e1, jj, ii, i1,k1,i2,k2,ar,ai
        Call Add_tmat_list(e1,ar,ai,jj,ii,i1,i2,k1,k2)
        ie = 0
        Do j=1,nn; if(e1.eq.ev(j)) then; ie=j; Exit; end if; End do  
        if(ie.eq.0) then; nn=nn+1;  ev(nn)=e1; ie=nn; end if
        ipe(i)=ie 
      End do
      if(nn.ne.ne) Stop 'ne in tmat.done is incorrect'
      if(ncase.ne.ndata) Stop 'ncase in tmat.done is incorrect'

      llmax = max(MAXVAL(l1(1:ncase)),MAXVAL(l2(1:ncase)))

      write(pri,'(/a,i5)') 'llmax   = ',llmax
      
      write(pri,'( a,i5)') 'number of energies:  ',ne

      write(pri,'(a,i5)') 'number of T-elements:',ndata

!------------------------------------------------------------------------
! ... Coulomb phases: cp(l+1) = s(l+1)-s(0) = cp(l) + tn^-1(q/(l+1)):

      if(ion.gt.0) then

       Allocate(CP(0:llmax,ne));  CP = 0.d0
       Do ie = 1,ne
        q = ion/sqrt(ev(ie)-E1)
        Do l=1,llmax
         s=l; CP(l,ie)=CP(l-1,ie)+DATAN2(q,s)
        End do
       End do

      end if

!------------------------------------------------------------------------------
! ... files

      AS = 'in  a_o^2'
      if(i16.gt.0) write(AS,'(a,i2,a)') 'in 10^',i16,' cm^2'
      Snorm = 1.d0; if(i16.ne.0) Snorm = a0**2  *  10**(i16-16)

      if(ion.eq.0.or.itr1.ne.itr2) then
       write(AF_sec,'(a,i2.2,a,i2.2)') 'mm_',itr1,'_',itr2
       open(nur,file=AF_sec)
       write(nur,'(3x,5(a,10x))')  'E_eV', 'sec ', 'mt  ', 'E_Ry', trim(AS)         
      end if

      write(AF_dcs,'(a,i2.2,a,i2.2)') 'dcs_',itr1,'_',itr2
      open(nud,file=AF_dcs)
      write(nud,'(3x,a,7x,a,7x,a,5x,a)')  'Angle',  'E_Ry', 'DCS', trim(AS)

!-------------------------------------------------------------------------------
! ... amplitude:  (we are working in 2J representation
   
      Allocate(FR(2,2,JT1+1,JT2+1,ng),FI(2,2,JT1+1,JT2+1,ng))

      Do ie = 1,ne; E = ev(ie)-etarg(itr1); xk = sqrt(E);  q = ion/xk
       FR = 0.d0
       FI = 0.d0

       Do i = 1,ncase; if(ipe(i).ne.ie) Cycle
        JP = jjpar(i); ll1=l1(i)+l1(i); ll2=l2(i)+l2(i); j1=jk1(i); j2=jk2(i) 
        Call Fano_tma (l1(i),l2(i),TR(i),TI(i),br,bi)

        ! introduce factor -i ->  may be important for Coulomb case (?)             
        ar = bi;  ai = - br
 
        if(ion.gt.0) then      
         sigma = CP(l1(i),ie) + CP(l2(i),ie)
         sr = cos(sigma); si = sin(sigma)
         br = ar*sr - ai*si; bi = ar*si + ai*sr               
         ar = br; ai = bi
        end if
 
        S=sqrt((l1(i)+l1(i)+1.d0)/2)                           !!! 1/2
 
        Do imt1 = 1,JT1+1;  mt1 = -JT1 + (imt1-1)*2;  mj1 = mt1   
         S1 = S  * CLEBSH2(JT1,mt1,ll1,0,j1,mj1)
         if(S1.eq.0) Cycle;   

 
        Do imu1 = 1,2; mu1 = -1; if(imu1.eq.2) mu1=1; M   = mj1 + mu1
         S2 = S1 * CLEBSH2(j1,mj1,1,mu1,JP,M)            
         if(S2.eq.0) Cycle 
       
        Do imu2 = 1,2; mu2 = -1; if(imu2.eq.2) mu2=1; mj2 = M - mu2 
         S3 = S2 * CLEBSH2(j2,mj2,1,mu2,JP,M)
         if(S3.eq.0) Cycle 
 
        Do imt2 = 1,JT2+1;  mt2 = -JT2 + (imt2-1)*2;  ml2 = mj2 - mt2
         S4 = S3 *CLEBSH2(JT2,mt2,ll2,ml2,j2,mj2)
         if(S4.eq.0) Cycle 
        Do ig=1,ng; x = COS(G(ig)*pi/180)
         P =  S4 * ALEGFM (l2(i),ml2/2,x,1)  
         FR(imu1,imu2,imt1,imt2,ig) = FR(imu1,imu2,imt1,imt2,ig) + P*ar 
         FI(imu1,imu2,imt1,imt2,ig) = FI(imu1,imu2,imt1,imt2,ig) + P*ai 
        End do
 
        End do; End do; End do; End do
 
       End do  !  ncase

!----------------------------------------------------------------------------------

       if(ion.gt.0.and.itr1.eq.itr2) then

        Do ig = 1,ng
         T = G(ig)*PI/180;   s = sin(T/2); s = s * s;  if(s.eq.0.d0) Cycle
         sn = q/(2*xk*s)
         sigma = q*log(s)      
         sr = cos(sigma)*sn; si = sin(sigma)*sn
         Do imu = 1,2; Do imt=1,JT1+1              
          FR(imu,imu,imt,imt,ig) = FR(imu,imu,imt,imt,ig) + sr 
          FI(imu,imu,imt,imt,ig) = FI(imu,imu,imt,imt,ig) + si 
         End do; End do
        End do

       end if
 
!----------------------------------------------------------------------------------
! ... diff. cross-sections:

       SN = 1/(2*E*(JT1+1))  * snorm

       Do ig = 1,ng
        S = SUM(FR(:,:,:,:,ig)*FR(:,:,:,:,ig)) + SUM(FI(:,:,:,:,ig)*FI(:,:,:,:,ig))
        DCS(ig) = S * SN 
        DS(ig) = DCS(ig) * SIN(G(ig)/180*pi)
        DM(ig) = DS(ig) * (1.d0-COS(G(ig)/180*pi))
        if(ion.eq.0.or.G(ig).gt.0.d0) &
        write(nud,'(f8.1,f12.6,E15.5)') G(ig),ev(ie),DCS(ig)
       End do 
       
!  ... angle-integrated cross-sections:
  
       if(ion.eq.0.or.itr1.ne.itr2) then

        sig = ZIS(hg,ng,DS) *  2*pi    
        sigM = ZIS(hg,ng,DM) *  2*pi    
        write(nur,'(f10.5,2e15.5,f10.6)')  E*Ry,sig,sigM,ek(ie)
       end if
  
      End do    ! over energies, ie

      End  ! program sec_dif_JK


!======================================================================
      Subroutine inf_sub
!======================================================================
!     provides screen information about utility
!----------------------------------------------------------------------

      Character :: A

      iarg = command_argument_count()
      if(iarg.eq.0) Return
      Call GET_COMMAND_ARGUMENT(1,A)        
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                            ',&
' sec_dif_JK provides differential and angle-integrated      ',&
' (ordinary and momentum transfer) cross sections            ',&
' for given transition ii -> ff                              ',&
'                                                            ',&        
'        zarm.tma + target  =>   dcs_ii_ff, mt_ii_ff         ',&
' or     tmat.done          =>   dcs_ii_ff, mt_ii_ff         ',&
'                                                            ',&        
' Call as:                                                   ',&
'          sec_dif_JK_ampl itr1=ii itr2=jj i16=.. ifano=0|1  ',&
'                          ek1=...  ek2=... or ekk=..        ',&
'                          gmin=... gmax=... gstep=...       ',&
'                          dcs=0|1  tdone=0|1                ',&
'                                                            ',&        
' ii - index of initial state  (default - 1)                 ',&
' jj - index of final state    (default - 1)                 ',&
' ek1 - min. electron energy                                 ',&
' ek2 - max. electron energy                                 ',&
' ekk - if /= 0, ekq=ek2=ekk and only tmat.done_inp created  ',&
' i16 - 0 - sigma in a.u.                                    ',&
'      >0 - sigma in 10-i16 cm^2 (default i16=16)            ',&
' ifano = 1  - Fano phase convention                         ',&
'       = 0  - Condon-Shortly phase convention, default      ',&
' dcs   = 1  - dif. cross sections in dcs_nn_nn              ',&
'       = 0  - only ICS in mt_nn_nn                          ',&
' tdone = 1  - tmat.done is supposed to exist                ',&
'       = 0  - tmat.done is created                          ',&
' JJ_top  - if =1 - tmat.done_inp ->> tmat.done_out          ',&
'           where JJ values extended to JJ_top value         '
      Stop ' '

      End Subroutine inf_sub

