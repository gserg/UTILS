!=====================================================================
!     UTILITY  sec_ionb
!=====================================================================
!     provides ionization+excitation cross sections
!     (based on the psedo-state de-composition)
!
!     Arguments: is   -  initial state   (default: is=1)
!
!     Input files:
!
!      target_ps     - description of e - A scattering
!      zarm.omb_top  - collection for collision strenths
!      bound_ovl     - pseudostates  projections
!
!     Output files:
!
!     Example of call:  sec_ionb  is=... i16=...  om=...
!
!---------------------------------------------------------------------
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)

      Real(8), allocatable :: e(:),om(:),y(:) 
      Real(8), allocatable :: st(:,:), om_ion(:,:), pom(:,:) 
      Real(8), allocatable :: wtarg(:,:) 

      Character(20)  :: AF
      Character(400) :: AS

      Real(8) :: Ry = 13.6057

      Integer :: nut=1;  Character(20) :: AF_t = 'target_ps'
                         Character(20) :: AF_o = 'bound_ovl'
      Integer :: nuo=2;  Character(20) :: AF_om = 'zarm.omb_top'
      Integer :: nus=4;  Character(20) :: AF_st = 'sec_ion_##'

!----------------------------------------------------------------------
      Call get_command_argument(1,AF)  
      if(AF.eq.'?') then
      write(*,'(a)') &
'                                                                ',&
'     SEC_IONB provides ionization+excitation cross section      ',&
'     based on the psedo-state de-composition                    ',&
'                                                                ',&
'       target_ps, bound_ovl, zarm.omb_par  -->  sec_ion_nn      ',&
'                                                                ',&
'     Call as:    sec_ionb  [is=... i16=...  om=...]             ',&
'                                                                ',&
'     is  -  index of initial state                              ',&
'     i16 -  control the units for output cross sections         ',&
'     om  -  input omega file                                    ',&
'                                                                '
      Stop ' '                                                                    
      end if
 
!----------------------------------------------------------------------
! ... target information:

      Call Check_file(AF_t)
      Open(nut,file=AF_t)
      Call R_target(nut)
      Call R_channels(nut)
      np = ntarg;  Call Read_ipar(nut,'np',np)
      ni = ntarg;  Call Read_ipar(nut,'ni',ni)
      write(*,*) 'np =', np
      write(*,*) 'ni =', ni
      Close(nut)

      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion
      if(ion.ne.0) Stop ' ion <> 0 '
      E1 = etarg(1); etarg = (etarg-E1)*2.0

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,9)
      Ry = au_eV/2.d0

! ... projection coefficients:

      Call Check_file(AF_o)
      Open(nut,file=AF_o)
      Call Read_ipar(nut,'ntarg_ion',ntarg_ion)
      Call Read_ipar(nut,'nps',nps)
      if(nps.ne.ntarg) Stop 'nps <> ntarg'
      Allocate(wtarg(ntarg_ion,ntarg))
      Do i = 1,nps
       read(nut,*) j,IL,IS,IP,EP, wtarg(:,i)
      End do 
      Close(nut)

! ... define transition under consideration:

      is = 1;  Call Read_iarg('is',is)
      if(is.gt.ntarg) is=ntarg
      if(is.lt.1) is=1
      i16=0;  Call Read_iarg('i16',i16)

! ... statistical weight for the initial state:

      g=iabs(IStarg(is))*(2.0*Ltarg(is)+1)
      if(IStarg(is).eq.0) g=Ltarg(is)+1

!----------------------------------------------------------------------
! ... define energies:

      Call Read_aarg('om',AF_om)
      Open(nuo,file=AF_om,status='OLD')

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
      Allocate( y(ntr),E(me),OM(me),st(ntarg,me),om_ion(ntarg_ion,me))
      st = 0.d0; y = 0.d0; om = 0.d0
!----------------------------------------------------------------------
! ... read data:

      ne=0
      rewind(nuo)
    5 read(nuo,*,end=15) x,ns

      read(nuo,*) (y(i),i=1,ns)
      if(x.le.etarg(is)) go to 5

      ie=0
      Do i=1,ne; if(x.ne.e(i)) Cycle; ie=i; Exit; End do
      if(ie.eq.0) then
       ne=ne+1; if(ne.gt.me) Stop 'ne > me';  ie=ne
      end if

      S=0.d0

      Do i=1,ntarg

       if(is.le.i) then
        itr=(i-1)*i/2+is
        if(i.gt.np)  itr = (np+1)*np/2 + (i-np-1)*ni + is
       else
        itr=(is-1)*is/2+i
       end if

       if(itr.gt.ns) Exit

       S = S + y(itr)
       st(i,ie)=y(itr)

      End do

       e(ie)=x;  om(ie)=S
      go to 5

   15 write(*,*) 'ne =', ne
!--------------------------------------------------------------------------
! ... ordering the energies:

      Do i=1,ne-1; Do j=i+1,ne
       if(e(i).gt.e(j)) then
        S=e(i); e(i)=e(j); e(j)=S
        S=om(i); om(i)=om(j); om(j)=S
        Do it = 1,ntarg
         S=st(it,i); st(it,i)=st(it,j); st(it,j)=S
        End do
       end if
      End do; End do
!--------------------------------------------------------------------------
! ... transform to cross-sections:

      Do ie=1,ne
       es=e(ie)-etarg(is)
       ss = 3.1415926/g/es                             ! in ao^2

       if(i16.ne.0) ss = SS * 0.28003 * 10**(i16-16)

       om(ie) = om(ie) * ss
       st(:,ie) = st(:,ie) * ss
      End do
!-------------------------------------------------------------------------
! ... decomposition:

      Allocate(pom(ntarg_ion,ntarg))

      write(AF_st,'(a,i2.2)') 'sec_ion_',is;  Open(nus,file=AF_st)

      ia = 1      
      write(AS(ia:),'(a)') '    eV';   ia = ia+14
      write(AS(ia:),'(a)') '    Ry';   ia = ia+14
      Do i=1,ntarg_ion
       write(AS(ia:),'(a,i1)') '  t',i;  ia=ia+16
      End do
      write(AS(ia:),'(a)') '  S_ion';  ia = ia+16
      write(AS(ia:),'(a)') '  S_els';  ia = ia+16
      write(AS(ia:),'(a)') '  S_tot';  ia = ia+16

      if(i16.eq.0) write(AS(ia:),'(a)') '   in au '
      if(i16.gt.0) write(AS(ia:),'(a,i2.2,a)') '   in 10^-',i16,'  cm^2'

      write(nus,'(a)') trim(AS)

      om_ion = 0.d0
      Do ie=1,ne
       ss = 0.d0
       Do it = 1,ntarg
        om_ion(:,ie) =  om_ion(:,ie) + ST(it,ie) * wtarg(:,it)
        pom(:,it) =  ST(it,ie) * wtarg(:,it)
       End do  ! it

       S_ion = sum(om_ion(:,ie))
       S_els = ST(is,ie)
       S_tot = om(ie)
       if(s_ion.gt.0) &
       write(nus,'(f14.8,f12.6,80e16.8)') (e(ie)-etarg(is))*Ry, e(ie),om_ion(:,ie),S_ion,S_els,S_tot

      End do  ! ie

      Close(nus)

      End  ! UTILITY  sec_ionb


