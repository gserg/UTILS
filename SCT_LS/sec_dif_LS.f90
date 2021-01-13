!======================================================================
!      Differential cross section (case of neutral atoms in LS scheme)
!======================================================================
!      Input:        zarm.tma     - T-matrix elements for all LSp 
!                    zarm.tmb     - new format (more compact)
!                    target       - description of channels
!
!      Output:       sec_nn_nn    - angle-integrated direct and momentum 
!                                   trnsfer cross-sections 
!                    dcs_nn_nn    - DCS cross_sections
!
!      Call as:      sec_dif_LS  itr1=...  itr2=... 
!                               [i16=... ifano=...]
!                               [ek1=...  ek2=... ]
!                               [Glow=... Ghigh=... Gstep=...]
!======================================================================
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)

      Integer :: nut= 1;  Character(80) :: AF_tar='target'
      Integer :: num= 2;  Character(80) :: AF_tm ='zarm.tma'
                          Character(80) :: AF_tmb='zarm.tmb'
      Integer :: nur=11;  Character(80) :: AF_sec='sec_nn_nn'
      Integer :: nud=12;  Character(80) :: AF_dcs='dcs_nn_nn'
      Integer :: nus= 9;  ! scratch file  

      Real(8), allocatable :: ek(:),ev(:),TR(:),TI(:),MR(:),MI(:)
      Real(8), allocatable :: Ak(:,:),B(:),g(:)
      Integer, allocatable :: ipe(:),l1(:),l2(:),klsp(:)
      Integer, allocatable :: ml1(:),ml2(:),mlj(:),ms(:)
      
      Character(200) :: AS 
      Real(8) :: Ry = 13.6057,  a0 = 0.529177
      Integer :: ke = 5000

      Call Inf_sec_dif_LS  
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

      llmax = MAXVAL(lch)
      km=llmax+llmax+1
      kl=llmax+1

      istyle=1; if(np.eq.ntarg.and.ni.eq.ntarg) istyle=0

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

      idif=1;   Call Read_iarg('idif',idif)
     
!----------------------------------------------------------------------
! ... check transition:

      itr1=1;   Call Read_iarg('itr1' ,itr1 )
      itr2=1;   Call Read_iarg('itr2' ,itr2 )
      ifano=0;  Call Read_iarg('ifano',ifano)
      i16=0;    Call Read_iarg('i16'  ,i16  )

      if(itr1.lt.1)     itr1=1
      if(itr1.gt.np)    Stop 'itr1 > nphysical states'
      if(itr2.lt.itr1)  itr2=itr1
      if(itr2.gt.ntarg) itr2=ntarg

      if(ion.ne.0.and.itr1.eq.itr2) Stop 'ion <> 0, and itr1=itr2: not a neutral atom'

      IL1=ltarg (itr1); IL2=ltarg (itr2)
      IS1=istarg(itr1); IS2=istarg(itr2)

      ekk = 0.d0;  Call Read_rarg('ek' ,ekk)
      if(ekk.le.0.d0) then
       ek1 = 0.d0;  Call Read_rarg('ek1' ,ek1)
       ek2 = 0.d0;  Call Read_rarg('ek2' ,ek2)
      else
       ek1=ekk; ek2=ekk
      end if
      if(ek2.lt.ek1) ek2=ek1

      Glow  = 0;    Call Read_rarg('Glow' ,Glow)
      Ghigh = 180;  Call Read_rarg('Ghigh',Ghigh)
      Gstep = 1;    Call Read_rarg('Gstep',Gstep)

      ng = NINT( (Ghigh-Glow)/Gstep ) + 1
      Allocate(G(ng))
      Do i=1,ng; G(i) = Glow + (i-1)*Gstep; End do

!----------------------------------------------------------------------
! ... find energies:

      if(istyle.eq.0) then
       Call Check_file(AF_tm); Open(num,file=AF_tm)
      else
       Call Check_file(AF_tmb); Open(num,file=AF_tmb)
      end if

      me = ke; Allocate(ek(me)); ek = 0.d0
      ne=0
   10 Continue
      if(istyle.eq.0) then
       read(num,*,end=20) e1,nopen,ntr,ilsp
       read(num,*) (ar,ai,i=1,ntr)
      else
       read(num,*,end=20) e1,nopen,kopen,ilsp,i1,i2,nj
       read(num,*) ((S1,S2,j=1,i),i=1,kopen)
       if(nopen.gt.kopen.and.nj.gt.0) &
        read(num,*) ((S1,S2,j=1,nj),i=kopen+1,nopen)
      end if

       if(e1.le.etarg(itr2)) go to 10
       if(ek1.gt.0.d0.and.e1.lt.ek1) go to 10
       if(ek2.gt.0.d0.and.e1.gt.ek2) go to 10

      k = 0;  Do ie=1,ne; if(ek(ie).ne.e1) Cycle; k=1; Exit; End do
       
      if(k.eq.0) then; ne=ne+1; ek(ne)=e1; end if

      if(ne.eq.me) then
       Allocate(ev(me)); ev=ek; Deallocate(ek)
       me=ne+ke; Allocate(ek(me)); ek(1:ne)=ev(1:ne); Deallocate(ev)
      end if

      go to 10
   20 write(*,*) 'number of energies:   ',ne
      if(ne.eq.0) Stop 'Stop: ne = 0'

      Call Rsort(ne,ek)
!----------------------------------------------------------------------
! ... read T-matrix:

      mdim=mch*(mch+1)/2;  Allocate (TR(mdim),TI(mdim))

      ncase=0                   !  different T-matrix elements
      rewind(num)
   1  if(istyle.eq.0) then
       read(num,*,end=2) e1,nopen,ntr,ilsp
       read(num,*) (TR(i),TI(i),i=1,ntr)
      else
       read(num,*,end=2) e1,nopen,kopen,ilsp,i1,i2,nj
       if(kopen.gt.nopen)     Stop ' kopen > nopen'
       ntr = kopen*(kopen+1)/2
       read(num,*) (tr(i),ti(i),i=1,ntr)
       ktr = (nopen-kopen)*nj
       if(ktr.gt.0) read(num,*) (tr(i),ti(i),i=ntr+1,ntr+ktr)
      end if

      if(ilsp.gt.nlsp)       Stop ' ilsp  > nlsp'
      if(ilsp.lt.1)          Stop ' ilsp < 0'
      if(nopen.gt.nch(ilsp)) Stop ' nopen > nch(ilsp)'

      ie=0
      Do i=1,ne; if(ek(i).ne.e1) Cycle; ie=i; Exit; End do
      if(ie.eq.0) go to 1

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
        write(nus) ie,ilsp,lch(ilsp,ich1),lch(ilsp,ich2),TR(itr),TI(itr)
       End do
      End do

      go to 1
    2 write(*,*) 'number of T-elements: ',ncase
      if(ncase.eq.0) Stop 'number of T-elements = 0'
    
      Deallocate(TR,TI)
      Allocate(TR(ncase),TI(ncase),l1(ncase),l2(ncase),klsp(ncase), &
               ipe(ncase))
      rewind(nus)
      Do i = 1,ncase
	read(nus) ipe(i),klsp(i),l1(i),l2(i),TR(i),TI(i)
      End do

!----------------------------------------------------------------------
! ... define  M-elements:

      Allocate(MR(km),MI(km))
      mcase=0                 
      rewind(nus)
    3 m=0
      Do i = 1,ncase
       if(klsp(i).gt.0) m=i; if(m.gt.0) exit
      End do
      if(m.eq.0) go to 4

      m1=l1(m); m2=l2(m); ilsp=klsp(m); IS=ispar(ilsp); ie=ipe(m)

      S = (m1+m1+1)*(m2+m2+1)*IS
      S = sqrt(S)*(-1)**(m1+m2)

      mj1=max(iabs(m1-m2),iabs(IL1-IL2))
      mj2=min(iabs(m1+m2),iabs(IL1+IL2))

      MR = 0.d0;  MI = 0.d0
      Do i=1,ncase
        if(ipe(i).ne.ie) Cycle
        if(klsp(i).eq.0) Cycle
        if(l1(i).ne.m1.or.l2(i).ne.m2) Cycle
        ilsp=klsp(i)
        if(ISPAR(ilsp).ne.IS) Cycle
        IL=LPAR(ilsp)
        SL = S * (IL+IL+1) * (-1)**IL
        Do mj = mj1,mj2
          SS = Z_6jj(IL1,IL2,mj,m2,m1,IL)
          if(SS.eq.0.d0) Cycle
          SS = SS*SL
          k=mj+1; MR(k)=MR(k)+SS*TR(i); MI(k)=MI(k)+SS*TI(i)
        End do
        klsp(i)=0
      End do

      Do mj = mj1,mj2
        k=mj+1; if(abs(MR(k))+abs(MI(k)).eq.0.d0) Cycle
        write(nus) m1,m2,mj,IS,MR(k),MI(k),ie
        mcase = mcase + 1
      End do

      go to 3
    4 Continue
      Deallocate(TR,TI,l1,l2,klsp,ipe,MR,MI)
      write(*,*) 'number of M-elements: ',mcase
      Allocate(MR(mcase),MI(mcase),ml1(mcase),ml2(mcase),mlj(mcase), &
               ms(mcase),ipe(mcase))
      rewind(nus)
      Do i = 1,mcase
       read(nus)  ml1(i),ml2(i),mlj(i),ms(i),MR(i),MI(i),ipe(i)
      End do

!----------------------------------------------------------------------
! ... define Ak coefficients:
      
      Allocate(Ak(0:km,ne)); Ak = 0.d0

      Do i = 1,mcase;  Do j = i,mcase
	  
        if(ipe(i).ne.ipe(j)) Cycle;  ie = ipe(i)

        if(ms(i) .ne.ms(j) ) Cycle     ! total spin
        if(mlj(i).ne.mlj(j)) Cycle     ! momentum transfer
 
        kmin=max(iabs(ml1(i)-ml1(j)),iabs(ml2(i)-ml2(j)))
        kmax=min(iabs(ml1(i)+ml1(j)),iabs(ml2(i)+ml2(j)))
        Do k = kmin,kmax,2
		 
        S =  Z_3j0(ml1(i),ml1(j),k) * Z_3j0(ml2(i),ml2(j),k)
        S =  S * Z_6jj(ml1(i),ml2(i),mlj(i),ml2(j),ml1(j),k)
        if(S.eq.0.d0) Cycle
        S =  S * (mlj(i)+mlj(i)+1) * (k+k+1) * (-1)**(mlj(i)+k)
        SR =  MR(i)*MR(j) + MI(i)*MI(j)
 
        if(ion.eq.0.or.itr1.ne.itr2) then
          SM =  SR
        else
          Stop ' ion <> 0 '
        end if
 
        S = S*SM;  if(i.ne.j) S = S * 2.d0
 
        kz = ml1(i)-ml2(i)+ml2(j)-ml1(j)
        if(mod(kz,2).ne.0) Stop ' Summa over l is odd !'
        kz = kz/2;  if(ifano.eq.0) S=S*(-1)**kz
        
        Ak(k,ie) = Ak(k,ie) + S
        End do

      End do; End do   ! over M-elements

!----------------------------------------------------------------------
! ... diff. cross-sections:

      write(AF_dcs,'(a,i2.2,a,i2.2)') 'dcs_',itr1,'_',itr2
      Call Read_aarg('out',AF_dcs)
      iout = 1
      Call Read_iarg('style',iout)

      open(nud,file=AF_dcs)

      write(AF_sec,'(a,i2.2,a,i2.2)') 'sec_',itr1,'_',itr2
      open(nur,file=AF_sec)

      AS = 'in  a_o^2'
      if(i16.gt.0) write(AS,'(a,i2,a)') 'in 10^',i16,' cm^2'
      if(iout.eq.1) &
      write(nud,'(4(3x,a,4x))')  'Angle', 'DCS',  'E_Ry', trim(AS)
      if(iout.eq.2) &
      write(nud,'(4(3x,a,4x))')  'E_eV', 'DCS',  'Angle', trim(AS)

      write(nur,'(5(6x,a,3x))')  'E_eV', 'sec ', 'mt  ', 'E_Ry', trim(AS)         

      Allocate(B(0:km))

      PI  = abs(ACOS(-1.d0))
      Snorm = 1.d0; if(i16.ne.0) Snorm = 0.529177**2  *  10 ** (i16-16)

      Do ie = 1,ne

       S = 8.d0*(ek(ie)-etarg(itr1))*(IL1+IL1+1)*IS1     
       Ak(:,ie) = Ak(:,ie) / S
       A0 = Ak(0,ie)
       A1 = Ak(1,ie)
       sig  = 4*PI * Snorm * A0
       sigM = 4*PI * Snorm * (A0 - A1/3.d0)

       Do ig=1,ng
        if(idif.eq.0) Exit
        T = G(ig) / 180. * PI
        Call LEGPOL(T,km+1,B)   
        SM = SUM(Ak(0:km,ie)*B(0:km)) * Snorm
        if(iout.eq.1) &
        write(nud,'(f8.2,E14.5,f15.6)') G,SM,ek(ie)
        if(iout.eq.2) &
        write(nud,'(f15.6,E14.5,f8.2)') ek(ie)*Ry,SM,G
       End do   ! over degrees

       write(nur,'(f10.5,2e15.5,f10.6)')  ek(ie)*Ry,sig,sigM,ek(ie)

      End do   ! over energies

      close(nus,status='DELETE')

      End  ! program sec_dif_LS


!======================================================================
      Subroutine inf_sec_dif_LS
!======================================================================
!     provides screen information about sec_dif_LS utility
!----------------------------------------------------------------------
      Character :: A

      iarg = command_argument_count()
      if(iarg.eq.0) Return
      Call GET_COMMAND_ARGUMENT(1,A)        
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                            ',&
' sec_dif_LS provides differential and angle-integrated      ',&
' (ordinary and momentum transfer) cross sections            ',&
' for given transition ii -> ff                              ',&
'                                                            ',&
'        zarm.tma + target  =>   dcs_ii_ff, mt_ii_ff         ',&
'                                                            ',&
' Call as:                                                   ',&
'                                                            ',&
'        sec_dif_LS itr1=ii itr2=jj i16=-1|0|. ifano=0|1     ',&
'                   ek1=...  ek2=...                         ',&
'                   gmin=... gmax=... gstep=...              ',&
'                                                            ',&
' ii - index of initial state  (default - 1)                 ',&
' jj - index of final state    (default - 1)                 ',&
'                                                            ',&        
' i16 - controls the output units:                           ',&
'       0 - sigma in a.u.  (default)                         ',&
'      >0 - sigma in 10^-i16 cm^2                            ',&
'                                                            ',&        
' ifano = 1  - Fano phase convention                         ',&
'       = 0  - Condon-Shortly phase convention, default      ',&
'                                                            '

      Stop ' '

      End Subroutine inf_sec_dif_LS

