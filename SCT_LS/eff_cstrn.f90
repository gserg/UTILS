!======================================================================
!     UTILITY     E F F _ C S T R   - effective collision strengths
!
!     target, zarm.om, temperature  -->  eff_str.dat, eff_str.tab, [cs_ii_ff]
!
! Call as   eff_cstr  itr1=.. itr2=.. jtr1=.. jtr2=..  [cs=cs]  
!
! i11 i12 i21 i22  -  range of transitions
! cs - create separate  files  cs_ii_ff for transition ii-ff    
!---------------------------------------------------------------------

      Use target

      Implicit real(8) (a-h,o-z)

      Real(8), Allocatable ::  e(:), om(:), T(:), TR(:), ecs(:), matr(:)
      Real(8), Allocatable ::  eom(:), ome(:)
      Integer, Allocatable ::  iom(:), ipt(:)
      Real(8) :: sdif(-1:5)

      Real(8), parameter :: Ry = 13.6057, K_eV = 11604.5 

      Character(3) :: AT
      Character(8) :: aa, unit = 'K'
      Character(11), allocatable :: a(:)
      Character(80) :: AS

! ... files:

      Integer :: nut=1;   Character(20) :: targ   = 'target'
                          Character(20) :: temp   = 'temperature'
      Integer :: nuo=3;   Character(20) :: omega  = 'zarm.om'
                          Character(20) :: omegb  = 'zarm.omb'
                          Character(20) :: omegt  = 'zarm.omb_top'

      Integer :: nuq=8;   Character(20) :: dat    = 'eff_str.dat'
      Integer :: nuc=9;   Character(20) :: cs     = ' '
      Integer :: pri=66;  Character(20) :: eff    = 'eff_cstrn.log'
                                                  
      Call Inf_sub

      open(pri,file=eff)

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(targ)
      Open(nut,file=targ)
      Call R_target(nut)
      np=ntarg; Call Read_ipar(nut,'np',np)
      ni=ntarg; Call Read_ipar(nut,'ni',ni)
      Close(nut)
      ion = nz-nelc; zion=1.d0; if(ion.ne.0) zion=ion*ion
      E1 = etarg(1); etarg = (etarg-E1)*2.0

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)

      Call Read_aarg('cs',cs) 
      ics = Len_TRIM(cs)

!----------------------------------------------------------------------
! ... temperature information:

      Call Check_file(temp)
      Open(nut,file=temp)
      Call Read_ipar(nut,'nt',nt)  
      Allocate(T(nt),TR(nt),ecs(nt),a(nt))
      Do i=1,nt; read(nut,*) T(i); End do

      Call Read_apar(nut,'unit',unit)  
      Select case(unit)
       case('K');                 TR = T / (K_eV * Ry)
       case('log');               TR = 10**T / (K_eV * Ry)
       case('Ry');                TR = T            
       case('eV');                TR = T / Ry   
      End Select

      write(*,*) 'nt =',nt,'  - number of temperatures'

!-----------------------------------------------------------------------
! ... define transition under consideration:

      i11 = 1; Call Read_iarg('itr1',i11)  
      i12 = 0; Call Read_iarg('itr2',i12);  if(i12.eq.0) i12=ntarg
      i21 = 2; Call Read_iarg('jtr1',i21)  
      i22 = 0; Call Read_iarg('jtr2',i22);  if(i22.eq.0) i22=ntarg

      if(i12.gt.np) i12=np
      if(i22.gt.np) i22=np      

!-----------------------------------------------------------------------
! ... read omega:

      mdim = ntarg*(ntarg+1)/2;  Allocate(matr(mdim))
      if(Icheck_file(omega).eq.1) then
       open(nut,file=omega)
      elseif(Icheck_file(omegb).eq.1) then
       open(nut,file=omegb)
      elseif(Icheck_file(omegt).eq.1) then
       open(nut,file=omegt)
      else
       Stop 'no OMEGA-related file'
      end if

      ip = 0; me = 0
    1 read(nut,*,end=2) ek,ns 
      read(nut,*) (matr(i),i=1,ns)
      me = me +1
      ip = ip + ns
      go to 1
    2 Continue

      Allocate(ome(ip), eom(me), iom(0:me), ipt(me))
      write(*,*) 'ip = ',ip, '  memory = ', (ip*8.d0 + me*8.d0 +me*4.d0)/1024.0/1024,  '  Mb'

      ip = 0; ne=0; iom=0
      rewind(nut)
    3 read(nut,*,end=4) ek,ns 
      read(nut,*) (matr(i),i=1,ns)
      ne = ne + 1
      jp = ip + 1
      ip = ip + ns
      eom(ne) = ek
      iom(ne) = ip
      ome(jp:ip) = matr(1:ns)
      go to 3
    4 Close(nut)

      Call SortR(me,eom,IPT)

      write(*,*) 'me =',me,'  - number of energies'

!-----------------------------------------------------------------------
! ... open output files:

       Open(nuq,file=dat)
       write(nuq,*)
       write(nuq,'( 7x, 13f11.0 )')  T 
       write(nuq,*)

!-----------------------------------------------------------------------
! ... cycle over transitions:

      Allocate(e(0:me),om(0:me))

      Do itr1 = i11,i12
      Do itr2 = i21,i22

       ecs = 0.d0
       if(itr1.ge.itr2) Cycle

! ...  extract omega:

       itr = itr2*(itr2-1)/2 + itr1 
       if(ion.ne.0) itr = itr - itr2 + 1

       om = 0.d0;  e(0) = etarg(itr2);  om(0)=0;  ne=0
       Do j=1,me; i=IPT(j)
        if(iom(i)-iom(i-1).lt.itr) Cycle
        ne = ne + 1
        e(ne)  = eom(i)
        om(ne) = ome(iom(i-1)+itr)
       End do

       write(*,*) 'ne =',ne,'  - number of energies'

       e(0:ne) = e(0:ne) - Etarg(itr2);  de = Etarg(itr2)-Etarg(itr1)

       if(SUM(om(1:ne)).eq.0.d0) go to 20

!----------------------------------------------------------------------
! ... define met:

       if(IStarg(itr1).ne.IStarg(itr2)) then
        met = 2
       elseif(IStarg(itr1).eq.0.and.  &
              iptarg(itr1).ne.iptarg(itr2).and. &
              ltarg(itr1)+ltarg(itr2).ne.0.and. &
              iabs(ltarg(itr1)-ltarg(itr2)).le.2) then
        met = -1
       elseif(IStarg(itr1).ne.0..and.  &
              iptarg(itr1).ne.iptarg(itr2).and. &
              ltarg(itr1)+ltarg(itr2).ne.0.and. &
              iabs(ltarg(itr1)-ltarg(itr2)).le.1) then
        met = -1
       else
        met =  0
       end if

      write(6,'(a,i3,a,i3,a,i3)')  '   transition ',jtr1,' ---> ',jtr2, &
               '   met = ', met
     
!----------------------------------------------------------------------
! ... define effective collision strength for each temperature:

      ecs = 0.d0

      Do it = 1,nt

       tt=TR(it);  s = 0.d0;  e1=E(0);  s1 = OM(1) 

       Do ie=1,ne

        e2=E(ie); h = e2-e1; s2 = OM(ie) 

        u1 = e1/tt; u2 = e2/tt; du = u2-u1; ed = exp(-du)

!        S = S + (s1 + s2*ed) * exp(-u1) * h / 2 

        S = S + (s1-s2*ed-(s1-s2)*(1-ed)/du)*exp(-u1) * tt

!        s = s + h*(s2+s1)/2.d0 * exp(-(e1+h/2)/TT)
 

       e1 = e2; s1 = s2

!        e2=E(ie); h = e2-e1; s2 = OM(ie); if(s2.lt.0.d0) s2=s1 
!        s = s + h*(s2+s1)/2.d0 * exp(-(e1+h/2)/TT)
!        e1 = e2; s1 = s2

       End do

       e1 = E(ne); om1 = OM(ne)

       h = 0.1;  if(h.gt.0.01*TT) h=0.01*TT
       h = 0.01*tt
       kk=0
       Do 
        e2 = e1 + h
        if(met.eq.-1) then
          om2= om1/LOG(e1)*LOG(e2)
        else
          om2= om1*e1**met/e2**met
        end if

        u1 = e1/tt; u2 = e2/tt; du = u2-u1; ed = exp(-du)
        s1 = om1; s2 =om2
        ds =  (s1-s2*ed-(s1-s2)*(1-ed)/du)*exp(-u1) * tt

!        ds = h*(om1+om2)/2.d0 * exp( -(e1+h/2)/TT )

        s = s + ds

        e1=e2; om1=om2

        if(ds.lt.1e-3*S) Exit
        kk=kk+1
        if(kk.gt.10000000) then;  write(6,*) itr1,itr2,kk; Exit; end if
       End do

       ecs(it) = S /TT

      End do


   20 Continue

! ... output:


      write(nuq,'(i3,a1,i3,1P50E11.2)') &
                   itr1,' ',itr2,(ecs(i),i=1,nt)


      if(ics.gt.0) then
       write(AS,'(a,i3.3,i3.3)') trim(cs),itr1,itr2
       Open(nuc,file=AS)
       Do i=1,nt
        tx =  TR(i) * K_eV * Ry
        write(nuc,'(f16.8,f12.0,1Pd12.3)') log10(tx),tx,ecs(i)
       End do 
      end if

!----------------------------------------------------------------------

      End do; End do  ! over transitions


      End   !  program eff_str


!======================================================================
      Subroutine inf_sub
!======================================================================
!     provide screen information about add_farm utility
!----------------------------------------------------------------------
       
      Character :: A

      iarg = command_argument_count()
      if(iarg.eq.0) Return
      Call GET_COMMAND_ARGUMENT(1,A)        
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                       ',&
'     utility  EFF_CSTRN  generates set of effective collision strengths',&
'     based on the collision strengths given in one of the files        ',&
'                                                                       ',&
'       zarm.om,  zarm.omb, zarm,omb_top   (plus file "target")         ',&
'                                                                       ',&
'     for set of temperatures given in file "temperatures"              ',&
'                                                                       ',&
'     Call as:  eff_cstrn  [itr1=.. itr2=.. jtr1=.. jtr2=.. cs=..]      ',&
'                                                                       ',&
'     i11 i12 i21 i22  -  range of transitions                          ',&
'     cs - if given also created separate files cs_ii_ff                ',&
'     for transition ii-ff which are ready for plotting                 ',&
'                                                                       ',&
'     output files  - eff_str.dat, eff_str.tab, [cs_ii_ff]              ',&
'                                                                       '
      Stop ' '

      End Subroutine inf_sub



