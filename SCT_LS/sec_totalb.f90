!=====================================================================
!     UTILITY  sec_totalb
!=====================================================================
!
!     total collision strengths or cross sections for the given initial
!     target state
!
!     Arguments:      is   -  initial state            [1]
!                     iion -  number of bound states   [0]
!                     i16  -  output units -> 10^-i16  [16]
!
!     Input files:
!
!      target     -  description of scattering model
!      zarm.omb   -  collection for collision strengths
!                    (or zarm.omb_top, zarm.om, om=..)
!
!     Output files:
!
!      tt_nn_nn.dat  -   nn stands for value of 'is'      
!
!     Call as:  sec_totalb  is=.. iion=.. i16=.. [om=..   np=.. ni=..]
!                
!---------------------------------------------------------------------

      Use target; Use channels

      Implicit real(8) (A-H,O-Z)

      Real(8), Allocatable :: e(:), y(:), om_tot(:), om_ion(:), om_jon(:), &
                                          om_els(:), om_exc(:), om_dex(:), &
                                          om_nl(:), om_set(:,:)
      Integer, Allocatable :: ipt(:), iset(:)
      Character(200) :: AS     
                              
! ... files:                  
                              
      Integer :: nut=1;  Character(80) :: AF_t  = 'target'
      Integer :: nuo=2;  Character(80) :: AF_om = 'zarm.om'
                         Character(80) :: AF_omb= 'zarm.omb'
                         Character(80) :: AF_top= 'zarm.omb_top'
      Integer :: out=3;  Character(80) :: AF_tt = 'tt_nn_nn.dat'
      Integer :: nui=8;  Character(80) :: AF_ion= 'ion_nn_nn'
      Integer :: nus=9;  Character(80) :: AF_set= 'target_set'

      Call inf_totalb
!----------------------------------------------------------------------
! ... target information:

      Call Check_file(AF_t)
      Open(nut,file=AF_t)
      Call R_target(nut)
      Call R_channels(nut)
      np = ntarg;  Call Read_ipar(nut,'np',np)
      ni = ntarg;  Call Read_ipar(nut,'ni',ni)
      Close(nut)

      ion = nz-nelc; if(ion.ne.0) Stop ' ion <> 0 '

      E1 = etarg(1); etarg = (etarg-E1)*2.0

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

      Eion = 0.d0;  Call Read_rarg('Eion',Eion) 
!----------------------------------------------------------------------
! ... define transition under consideration:

      is = 1;    Call Read_iarg('is',is);   if(is.gt.ni) Stop 'is > ni'

      i16 = 16;  Call Read_iarg('i16',i16)
                 if(i16.eq.1) i16=16
                 if(i16.eq.2) i16=18

      iion = 1;  Call Read_iarg('iion',iion)

      itr1 = 0;    Call Read_iarg('itr1',itr1)
      itr2 = 0;    Call Read_iarg('itr2',itr2)

! ... define input file:

      if(Icheck_file(AF_top).eq.1) then
       Open(nuo,file=AF_top)
      elseif(Icheck_file(AF_omb).eq.1) then   
       Open(nuo,file=AF_omb)
      elseif(Icheck_file(AF_om).eq.1) then   
       Open(nuo,file=AF_om)
       np=ntarg; ni=ntarg
      else
       AS=' '; Call Read_aarg('om',AS)
       if(Icheck_file(AS).eq.0) Stop 'No omega file given'
       Open(nuo,file=AS)
       np=ntarg; Call Read_iarg('np',np)
       ni=ntarg; Call Read_iarg('ni',ni)
      end if

! ... statistical weight for initial state:

      g=iabs(IStarg(is))*(2.0*Ltarg(is)+1)
      if(IStarg(is).eq.0) g=Ltarg(is)+1

!----------------------------------------------------------------------
! ... define energies:

      me=0
      rewind(nuo)
    1 read(nuo,*,end=2) x,ns
      read(nuo,*) (S,ix=1,ns)
      if(x.le.etarg(is)) go to 1
      me = me + 1
      go to 1
    2 Continue
      write(*,*) ' me = ',me

! ... allocations:

      ntr = ntarg*(ntarg+1)/2
      Allocate( y(ntr),E(me),ipt(me) , om_tot(me), om_ion(me), om_jon(me), &    
                                       om_els(me), om_exc(me), om_dex(me), &
                                       om_nl(me)   )                           
! ... subsets of states:

      Allocate(iset(ntarg)); iset=0
      if(Icheck_file(AF_set).eq.1) then
       open(nus,file=AF_set)
       Do i=1,ntarg; read(nus,*) iset(i); End do
       kset = 1
      end if
      nset = maxval(iset)
      if(nset.gt.0) then
       Allocate(om_set(me,nset))
       om_set=0.d0
      end if


!----------------------------------------------------------------------
! ... read data:
      
      ne=0
      rewind(nuo)
    5 read(nuo,*,end=15) x,ns
      read(nuo,*) (y(i),i=1,ns)
      if(x.le.etarg(is)) go to 5

      ie=0;  Do i=1,ne; if(x.ne.e(i)) Cycle; ie=i; Exit; End do
      if(ie.eq.0) then
       ne=ne+1;  if(ne.gt.me) Stop 'ne > me';  ie=ne
      end if

      S_tot=0.d0; S_jon=0.d0; S_ion=0.d0; S_els=0.d0; S_dex=0.d0; S_exc=0.d0; S_nl=0.d0

      Do i=1,ntarg
       if(is.le.i) then
        itr=(i-1)*i/2+is
        if(i.gt.np)  itr = (np+1)*np/2 + (i-np-1)*ni + is
       else
        itr=(is-1)*is/2+i
       end if

       if(itr.gt.ns) Exit

       S_tot = S_tot + y(itr)                      ! total

       if(i.gt.iion) then
         S_ion = S_ion + y(itr)                    ! ion
       else
         S_jon = S_jon + y(itr)                    ! total-ion
       end if 

       if(i.eq.is) then
         S_els = S_els + y(itr)                    ! els
       elseif(i.lt.is) then                       
         S_dex = S_dex + y(itr)                    ! dexc
       elseif(i.le.iion) then
         S_exc = S_exc + y(itr)                    ! exc

         if(i.ge.itr1.and.i.le.itr2) S_nl=S_nl+y(itr)

       end if

       if(nset.gt.0) then
        ii = iset(i)
        if(ii.gt.0.and.i.gt.is) om_set(ie,ii) = om_set(ie,ii) + y(itr) 
       end if

      End do

      e(ie)=x
      om_tot(ie)=S_tot
      om_ion(ie)=S_ion
      om_jon(ie)=S_jon
      om_els(ie)=S_els
      om_exc(ie)=S_exc
      om_dex(ie)=S_dex
      om_nl(ie) =S_nl

      go to 5
   15 Continue

      write(*,*) ' ne = ',ne

! ... ordering the data:

      Call Sortr(ne,e,ipt)

! ... output data:

      write(AF_tt,'(a,i2.2,a,i2.2,a)') 'tt_',is,'_',is,'.dat'
      Open(out,file=AF_tt)

      write(AS,'(8(6x,a3,6x),a,i2,a,i4)') &
      ' eV','els','exc','ion','dex','tot','jon',' Ry','in  10^-',i16,' cm^2;    iion =',iion

      write(out,'(a)') trim(AS)

! ...

      write(AF_ion,'(a,i2.2)') 'ion_',is
      Open(nui,file=AF_ion)

      write(AS,'(3(6x,a3,6x),a,i2,a)') &
      ' eV','ion','Ry','in  10^-',i16,' cm^2'

      write(nui,'(a)') trim(AS)
      if(Eion.gt.0) &
       write(nui,'(f14.8,e15.6,f14.8)')  (Eion-etarg(is))*Ry, 0.d0,Eion 

      S = 3.1415926/g  

      if(i16.gt.0) S = S * 0.28003 * 10**(i16-16)

      Do je=1,ne; ie=IPT(je)
       es=e(ie)-etarg(is)

       S_tot = om_tot(ie) * S / es
       S_ion = om_ion(ie) * S / es
       S_jon = om_jon(ie) * S / es
       S_els = om_els(ie) * S / es
       S_exc = om_exc(ie) * S / es
       S_dex = om_dex(ie) * S / es

       write(out,'(f14.8,6e15.6,f14.8)') &
         es*Ry,S_els,S_exc,S_ion,S_dex,S_tot,S_jon,e(ie)

       if(S_ion.gt.0.d0) &
       write(nui,'(f14.8,e15.6,f14.8)')  es*Ry,S_ion,e(ie) 

      End do

      Close(out)


! ... output nl:

      if(itr1.ne.0) then

      write(AF_tt,'(a,i2.2,a,i2.2,a,i2.2,a)') 'nl_',is,'_',itr1,'_',itr2,'.dat'
      Open(out,file=AF_tt)

      write(AS,'(a,i2,a)') &
      ' eV  sigma  Ry   omega   in  10^-',i16,' cm^2'

      write(out,'(a)') trim(AS)
!      es = etarg(itr1)-etarg(is)
!      write(out,'(f14.8,e15.6,f14.8,e15.6)') es*Ry,0.e0,es,0.d0

      S = 3.1415926/g  

      if(i16.gt.0) S = S * 0.28003 * 10**(i16-16)

      Do je=1,ne; ie=IPT(je)

       es=e(ie)-etarg(is)

       if(e(ie).le.etarg(itr1)) Cycle

       S_nl = om_nl(ie) * S / es

       write(out,'(f14.8,e15.6,f14.8,e15.6)') &
         es*Ry,S_nl,e(ie),om_nl(ie)

      End do

      Close(out)

      end if  !  nl output

! ... output sets:

      if(nset.gt.0) then

      write(AF_tt,'(a,i2.2,a)') 'sets_',is,'.dat'
      Open(out,file=AF_tt)

      write(AS,'(a,i2,a)') &
      ' eV  Ry   sigma  in  10^-',i16,' cm^2'

      write(out,'(a)') trim(AS)

      S = 3.1415926/g  

      if(i16.gt.0) S = S * 0.28003 * 10**(i16-16)

      Do je=1,ne; ie=IPT(je)

       es=e(ie)-etarg(is)
       if(sum(om_set(ie,:)).eq.0.d0) Cycle
       om_set(ie,:) = om_set(ie,:) * S / es

       write(out,'(e15.6,f14.8,10e15.6)')  es*Ry, e(ie), om_set(ie,:), sum(om_set(ie,:))

      End do

      end if ! sets output


      End  ! UTILITY  sec_totalb                                       


!======================================================================
      Subroutine inf_totalb
!======================================================================
!     provide screen information about sec_totalb utility
!----------------------------------------------------------------------
      Character :: A=' '

      Call get_command_argument(1,A)  
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                       ',&
'     total collision strengths or cross sections for the given initial ',&
'     target state                                                      ',&
'                                                                       ',&
'     Arguments:      is   -  initial state            [1]              ',&
'                     iion -  number of bound states   [0]              ',&
'                     i16  -  output units -> 10^-i16  [16]             ',&
'                                                                       ',&
'     Input files:                                                      ',&
'                                                                       ',&
'      target     -  description of scattering model                    ',&
'      zarm.omb   -  collection for collision strengths                 ',&
'                    (or zarm.omb_top, zarm.om, om=..)                  ',&
'                                                                       ',&
'     Output files:                                                     ',&
'                                                                       ',&
'      tt_nn_nn.dat  -   nn stands for value of "is"                    ',&
'                                                                       ',&  
'      Call as:  sec_totalb  is=.. iion=.. i16=.. om=..                 ',&
'                                                                       '        
      Stop ' '

      End Subroutine inf_totalb
