!======================================================================
!     UTILITY      S E C _ T O P _ C B E  
!
!     zarm.omb_par --> zarm.omb_top
!
!======================================================================
!     genearte top-up omegas for all included transitions 
!     (based on the information in 'zarm.omb_par')
!
!     Call as:  sec_top  [par= top= ek1= ek2= eps_tail= eps_fail= eps_x=] 
!---------------------------------------------------------------------

      Use target; Use channels

      Implicit real(8) (a-h,o-z)

      Real(8), allocatable ::  fl(:), om(:), e(:), fail(:),coef(:)
      Real(8), allocatable ::  om_sum(:,:), om_top(:,:)
      Real(8), allocatable ::  fom(:,:,:), OS(:,:)

      Integer, allocatable ::  iop(:),jop(:), met(:), ic(:), jc(:)
      Integer, allocatable ::  IL2(:),IS2(:), IP2(:), nch2(:), ilsp2(:)

      Integer :: ke = 10240  !  initial number of energies 

      Real(8) :: eps_tail = 0.001
      Real(8) :: eps_x    = 0.025

      Integer :: np=0, ni=0
!----------------------------------------------------------------------
! ... files:

      Character(80) :: AF, label=' '
      Integer :: nup=11;  Character(80) :: targ  = 'target'
      Integer :: nut=12;  Character(80) :: par   = 'zarm.omb_par'
      Integer :: nuo=14;  Character(80) :: top   = 'zarm.omb_top'
      Integer :: nuq=15;  Character(80) :: oms   = 'zarm.omb_sum'
      Integer :: pri=16;  Character(80) :: out   = 'sec_top_omb.log'
      Integer :: nuf=17;  Character(80) :: inf   = 'sec_top_omb.fail'
      Integer :: nuc=18;  Character(80) :: ccc   = 'sec_top_coef_fail'

      Integer :: nus=19;  Character(80) :: AF_OS = 's_values'

      Integer :: nua=99;  ! scratch file

!----------------------------------------------------------------------
      Call Inf_sub

      Open(pri,file=out)

      Call CPU_time(t1)

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(targ)
      Open(nup,file=targ)
      Call R_target(nup)
      Call R_channels(nup)
      Close(nup)
      e1=etarg(1); Do i=1,ntarg; etarg(i)=(etarg(i)-e1)*2; End do
      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion

      Z = nz
      AWT = 0.d0;  Call Read_rarg('AWT',AWT)
      Call Conv_au (Z,AWT,au_cm,au_eV,pri)
      Ry = au_eV/2.d0

!-----------------------------------------------------------------------
! ... define other parameters:

      ek1 = 0.d0; Call Read_rarg('ek1',ek1)
      ek2 = 0.d0; Call Read_rarg('ek2',ek2)

      Call Read_rarg('tail',eps_tail)
      Call Read_rarg('x',eps_x)

      jtr1=0; Call Read_iarg('jtr1',jtr1)
      jtr2=0; Call Read_iarg('jtr2',jtr2)

      Call Read_aarg('label',label)

!-----------------------------------------------------------------------
! ... get s-values:

      open(nus,file=AF_OS) 
      read(nus,*) nsvalues     
      Allocate(OS(ntarg,ntarg)); OS =0.d0
      
      Do n=1,nsvalues
       read(nus,*) i,j,OS(i,j)
       OS(j,i) = OS(i,j)
      End do

!----------------------------------------------------------------------
! ... find energies:

      if(len_trim(label).gt.0) par = trim(par)//'_'//trim(label)
      Call Read_aarg('par',par)
      Call Check_file(par)
      Open(nut,file=par)

      me = ke; Allocate(e(me))

      ne=0; mom=0
    1 read(nut,*,end=2) e1,ilsp,nom,k1,k2,i1,i2
      read(nut,*) (S,i=1,nom)

      if(nom.gt.mom) mom=nom

      if(ek1.gt.0.d0.and.e1.lt.ek1) go to 1
      if(ek2.gt.0.d0.and.e1.gt.ek2) go to 1

      if(np.eq.0) np = i1       
      if(ni.eq.0) ni = i2       
      if(np.ne.i1) Stop 'diferent np'
      if(ni.ne.i2) Stop 'diferent ni'

      ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

      if(ie.eq.0) then; ne=ne+1; e(ne)=e1;  end if

      if(ne.eq.me) then
       open(nua,form='UNFORMATTED',status='SCRATCH')
       rewind(nua);   write(nua) (e(i),i=1,ne)
       Deallocate(e); me=me+ke; Allocate(e(me))
       rewind(nua);   read(nua) (e(i),i=1,ne)
      end if       

      go to 1
    2 write(*,'(a,i5,a,i10)') 'ne =',ne,'   mom = ',mom
      if(ne.eq.0) Stop 'nothing to do !'

      Call RSort(ne,e)

      write(*,'(a,2f15.8)') 'e(1),e(ne)', e(1), e(ne)
      write(*,'(a,2i10)') 'np,ni =',np,ni
!----------------------------------------------------------------------
! ... allocations:

      maxl=maxval(lpar)
      write(*,'(/a,3i10/)') 'maxl = ', maxl 

      maxll = maxval(lch(nlsp,1:nch(ilsp)))    ! ???


      S = mom * (nlsp+3) * ne * 8.0 / (1024*1024) 
      write(*,'(/a,f10.2,a/)') 'Memory required:  ', S, ' Mb' 
      if(S.gt.60000.d0) Stop ' > 60 GB'
      
      Allocate(fom(mom,nlsp,ne), fl(0:maxl), iop(ne), jop(ne), om(mom) )

      Allocate(om_top(mom,ne), om_sum(mom,ne) )

      fom = 0.d0;  iop = 0;  jop =0; om_top = 0.d0; om_sum = 0.d0; 

      mmm = np*(np+1)/2;  if(ion.ne.0) mmm=np*(np-1)/2

      if(np.lt.ntarg) mmm = mmm + (ntarg-np)*ni

      Allocate( fail(mmm), coef(mmm), met(mmm), ic(mmm), jc(mmm) )
      fail = 0.d0; coef = 0.d0; met = -2; 
      
      Do itr1 = 1,np
      Do itr2 = itr1,ntarg
       itr = Index_TR(ion,itr1,itr2,np,ni) 
       if(itr.eq.0) Cycle         
       ic(itr) = itr1
       jc(itr) = itr2
      End do; End do
      
!-----------------------------------------------------------------------------------
! ... check continuation: 

      i =  Icheck_file(ccc)
      if(i.gt.0) then
       open(nuc,file=ccc)
       read(nuc,*) mm;  if(mmm.ne.mm) Stop 'mmm <> mm'
       Do j=1,mmm
        read(nuc,*) ic(j),jc(j), fail(j), coef(j), met(j)
       End do
       close(nuc)
      end if

!----------------------------------------------------------------------
! ... read partial OM:

      rewind(nut)
    3 read(nut,*,end=4) e1,ilsp,nom,k1,k2,i1,i2
      read(nut,*) (om(i),i=1,nom)

      ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do
      if(ie.eq.0) go to 3

      if(iop(ie).eq.0) iop(ie)=k1
      if(jop(ie).eq.0) jop(ie)=k2

      if(iop(ie).ne.k1) then
       write(*,'(i5,f14.6,10i5)') ilsp, e1, ie, k1, iop(ie), IOPEN(ntarg,e1,etarg)
       Stop 'different iop'
      end if 

      if(jop(ie).ne.k2) Stop 'different jop'

      Do i=1,nom; fom(i,ilsp,ie) = om(i); End do

      go to 3
    4 Close(nut)

!-----------------------------------------------------------------------
! ... cycle over energies:

      if(len_trim(label).gt.0) top = trim(top)//'_'//trim(label)
      Call Read_aarg('top',top)
      Open(nuo,file=top)

      if(len_trim(label).gt.0) oms = trim(oms)//'_'//trim(label)
      Call Read_aarg('oms',oms)
      Open(nuq,file=oms)
     
      Do  ie=1,ne
       ek = e(ie)
       om = 0.d0
       if(iop(ie).le.0) Cycle
!-----------------------------------------------------------------------
! ... cycle over transitions:

      ntr1 = iop(ie);  ntr = ntr1*(ntr1+1)/2
      if(ion.ne.0)     ntr = ntr1*(ntr1-1)/2
      ntr2 = jop(ie);  nom = ntr + (ntr2-ntr1)*ni

      Do itr1 = 1,ntr1
      Do itr2 = itr1,ntr2

       if(jtr1.ne.0.and.jtr1.ne.itr1) Cycle
       if(jtr2.ne.0.and.jtr2.ne.itr2) Cycle

       itr = Index_TR(ion,itr1,itr2,np,ni) 
       if(itr.eq.0) Cycle         

       write(pri,'(15(''-''))')
       write(pri,'(a,2i5,f15.8)') 'transition: ',itr1,itr2,ek
       write(pri,'(15(''-''))')

!-----------------------------------------------------------------------
! ... find met:

       if(IStarg(itr1).ne.IStarg(itr2)) then
        met(itr) = -1                                      ! exchange
       elseif(ISTARG(itr1).ne.0.and.   &
              iptarg(itr1).ne.iptarg(itr2).and. &
              ltarg(itr1)+ltarg(itr2).ne.0.and. &
              iabs(ltarg(itr1)-ltarg(itr2)).le.1) then
        met(itr) =  0                                      ! dipole, LS
       elseif(ISTARG(itr1).eq.0.and.   &
              iptarg(itr1).ne.iptarg(itr2).and. &
              ltarg(itr1)+ltarg(itr2).ne.0.and. &
              iabs(ltarg(itr1)-ltarg(itr2)).le.2) then
        met(itr) =  0                                      ! dipole, JK
       else
        met(itr) =  1                                      ! non-dipole
       end if

       y = (EK-Etarg(itr1))/(Etarg(itr2)- Etarg(itr1))
 
       fl = 0.d0
       Do ilsp = 1,nlsp
        l = lpar(ilsp);  fl(l) = fl(l) + fom(itr,ilsp,ie)
       End do

       S=SUM(fl)
       om_sum(itr,ie)=S
       om_top(itr,ie)=S
       if(S.eq.0.d0) Cycle

       f1=0.d0; f2=0.d0; f3=0.d0
       Do il=0,maxl
        if(fl(il).eq.0.d0) Cycle
        f1=f2; f2=f3; f3=fl(il)        
        if(f2.eq.0.d0) f2=f3
        write(pri,'(i2,a1,D12.4,f10.3)') il,'.',f3,f2/f3-1.d0
       End do
       write(pri,'(15(''-''))')
       write(pri,'(a3,D12.4)') 'SUM',S
       write(pri,'(15(''-''))')

       if(IStarg(itr1).ne.IStarg(itr2)) then
        write(pri,'(a3,D12.4,a)') 'TOP',S, ' -  exchange, no top-up'

       elseif(fail(itr).gt.0.d0.and.ek.ge.fail(itr)) then
        S = coef(itr);  if(met(itr).eq.0) S = S * log(y)
        om_top(itr,ie)=S
        write(pri,'(a3,D12.4,a)') 'TOP',S, ' -  extrapolated'

       elseif(S*eps_tail.gt.f1+f2+f3) then
        write(pri,'(a3,D12.4,a)') 'TOP',S, ' -  small correction, no top-up'

       else

! ...  non-dipole:

       if(met(itr).eq.1) then
        x1=f1/f2-1.d0; x2=f2/f3-1.d0; xx=(x1+x2)/2
        if(x1.le.0.d0.or.x2.le.0.d0) fail(itr)=ek        
        xe=(EK-Etarg(itr1))/(EK-Etarg(itr2)) - 1.d0
        if(xx.lt.eps_x) fail(itr)=ek

        x = xx
        S = S + f3/x;   om_top(itr,ie)=S

        if(fail(itr).eq.0.d0)  coef(itr)= S

        write(pri,'(a3,D12.4,3F10.3,a)') 'TOP',S,xx,xe,x, &
         '  average, energy-driven, chosen '

        if(fail(itr).gt.0.d0.and.ek.ge.fail(itr)) then
         S = coef(itr)
         om_top(itr,ie)=S
         write(pri,'(a3,D12.4,a)') 'TOP',S, ' -  extrapolated'
        end if
       end if

! ... dipole:

       if(met(itr).eq.0) then

         ekk1 = EK-Etarg(itr1)
         ekk2 = EK-Etarg(itr2)
   
         lamda = maxl + min(Ltarg(itr1),Ltarg(itr2))   ! - ???
         if(IStarg(1).eq.0) lamda = maxll-2
            
         eps = 1.d-5
         F1 = Fdip0(ekk1,lamda,ekk2,lamda-1,eps,ifail)
         F2 = Fdip0(ekk1,lamda-1,ekk2,lamda,eps,ifail)
         if(ifail.ne.0) write(*,'(a,i5,2f10.5)') 'fail: ',lamda,ekk1,ekk2

         om1 = 16.d0/3 * lamda * OS(itr1,itr2) * F1**2
         om2 = 16.d0/3 * lamda * OS(itr1,itr2) * F2**2
   
         C = (1+lamda**2*ekk2)*om1 - (1+lamda**2*ekk1)*om2
         C = C /lamda**2 / (ekk1-ekk2) 

! ... it is neutral case:  ???

         write(pri,'(a3,D12.4,a)') 'TOP CBE',C, ' -  CBE'
         CC = C

         ii = 0;
         if(ii.eq.1) then   !  it is not working yet ???
         Lmin = lamda - max(Ltarg(itr1),Ltarg(itr2))
         LT1 =  Ltarg(itr1);  LT2 =  Ltarg(itr2); 
         Do ilsp  = 1,nlsp
          if(Lpar(ilsp).lt.lmin) Cycle;  ILT = lpar(ilsp); IST = ispar(ilsp)
          Do ich1 = 1,nch(ilsp); if(iptar(ilsp,ich1).ne.itr1) Cycle
          l1=lch(ilsp,ich1); if(l1.le.lamda) Cycle
          Do ich2 = 1,nch(ilsp); if(iptar(ilsp,ich2).ne.itr2) Cycle
          l2=lch(ilsp,ich2); if(l2.le.lamda) Cycle
          if(iabs(l1-l2).ne.1) Cycle
          F = Fdip0(ek1,l1,ek2,l2,eps,ifail); if(ifail.ne.0) write(*,*) 'ifail',itr1,itr2
          sll = 16.d0/3 * max(l1,l2) * OS(itr1,itr2) * F**2
          sfl =  Z_6jj(LT1,l1,ILT,l2,LT2,1)**2 * 3 * (2*ILT+1) 
          if(sfl.eq.0.d0) Cycle 
          sll = sll * sfl * IST/2/IStarg(itr1)
          C = C - sll
          End do; End do
         End do
         end if

         om_top(itr,ie)=om_sum(itr,ie) + CC
         write(pri,'(a3,2D12.4)') 'TOP CBE',CC,om_top(itr,ie)

       end if


       end if


       write(pri,'(15(''-''))')


      End do         !  over  itr2
      End do         !  over  itr1

      End do        ! over ie  

!-----------------------------------------------------------------------
!... fail information:

      open(nuc,file=ccc)
      rewind(nuc)
      write(nuc,'(i10,a)') mmm, ' - number of transitions'

      i = 0
      Do j= 1,mmm

        if(fail(j).eq.0.d0) &
        write(nuc,'(2i8,f16.8,E16.8,2i5)') ic(j),jc(j), fail(j), coef(j),met(j)
        if(fail(j).ne.0.d0) i = i + 1
        if(fail(j).ne.0.d0) &
        write(nuc,'(2i8,f16.8,E16.8,2i5)') ic(j),jc(j), fail(j), coef(j),met(j),i
        
      End do
      write(nuc,*) 'failed transitions: ',i
      Close(nuc)
 
      write(*,*) 'failed transitions: ',i
      
!----------------------------------------------------------------------
! ... output new 'topped' om:

      Do ie=1,ne 
       if(iop(ie).lt.0) Cycle

       ntr1 = iop(ie);  ntr = ntr1*(ntr1+1)/2
       if(ion.ne.0)     ntr = ntr1*(ntr1-1)/2
       ntr2 = jop(ie);  nom = ntr + (ntr2-ntr1)*ni

       write(nuo,'(F10.6,5i8,a)')  e(ie),nom,iop(ie),jop(ie),np,ni, &
                              '   e(ie),nom,iopen,jopen,np,ni'
       write(nuo,'(5D16.8)') (om_top(i,ie),i=1,nom)

       write(nuq,'(F10.6,5i8,a)')  e(ie),nom,iop(ie),jop(ie),np,ni, &
                              '   e(ie),nom,iopen,jopen,np,ni'
       write(nuq,'(5D16.8)') (om_sum(i,ie),i=1,nom)

      End do        ! over ie  

      write(*,*) 'failed energies: ',count(iop.lt.0) 
 
      Call CPU_time(t2)
      write(*,'(a,f10.1,a)') 'time = ',(t2-t1)/60,' min'
      write(pri,'(a,f10.1,a)') 'time = ',(t2-t1)/60,' min'
                                                                 
      End   !  program sec_top


!======================================================================
      Subroutine inf_sub
!======================================================================
!     provide screen information about sec_top_TM utility
!----------------------------------------------------------------------
       
      Character :: A=' '

      Call get_command_argument(1,A)  
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                  ',&
'SEC_TOP tops up the zarm.omb_par file based on the                ',&
'geometric series approximation. For dipole transitions it requiers',&
'the "s_values" file.                                              ',&
'                                                                  ',&
'Arguments: par  - file with partial omega"s    [zarm.omb_par]     ',&
'           top  - file with topped up omega"s  [zarm.omb_top]     ',&
'           ek1, ek2 - energy interval in Ry                       ',&
'           eps_tail - tolerence for top-up contribution  [0.01]   ',&
'           eps_x    - tolerence for geom. series parameter [0.001]'
      Stop ' '

      End Subroutine inf_sub


