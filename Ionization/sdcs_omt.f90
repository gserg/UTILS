!=====================================================================
!     UTILITY  sdcs_omt
!=====================================================================
!     single differential cross sections for electron-impact ionization
!     based on the excitation collision strengths for pseudostates and
!     their projections to the real ionic continuum
!
!     incident and ejected electron energies defined from file "projections" 
!
!     INPUT FILES:
!
!     target_ps   - description of e - A  scattering 
!     target_ion  - description of e - A+ scattering
!     zarm.omb    - collision strenths for e - A problem
!     prodections - ion-states projections for given initial energy (EK)
!                   and set of ejected energies (nq)
!                   (obtained with program BSR_SEn)
!           
!     OUTPUT FILES:  sdcs_omt_nn (nn indicates final ionic state)      
!
!     File names can be re-definded through the arguments as:
!        AF_t=...  AF_i=... om=... p=... out=...
!
!     PARAMETERS:
!
!     is [1]     -  initial atomic state
!     iion       -  possible redifinition of number of bound pseudostates
!     eps_ek     -  tollerance for initial electron energy
!     i16  [16]  -  define units for resulting sdcs (i16=0 mens a.u) 
!     ion_state [0]  - indicate specific final ion state;
!                      0 means all states 
!---------------------------------------------------------------------

      Use target
      Use channels
      Use target_ion,   only: ntarg_ion => ntarg, etarg_ion => etarg
      Use channels_ion, only: nlsp_ion  => nlsp,  mch_ion => mch, &
                              iptar_ion => iptar

      Implicit real(8) (A-H,O-Z)

      Character(3)  :: AT
      Character(20) :: AF
      Character(200) :: AS, BS

      Real(8), allocatable :: om(:),y(:),sd(:),sdd(:), S_ion(:)
      Real(8), allocatable :: Eq(:), qq(:), dd(:,:,:)

      Integer :: nut=1;  Character(20) :: AF_t  = 'target_ps'
                         Character(20) :: AF_i  = 'target_ion'
                         Character(20) :: AF_p  = 'projections'
      Integer :: nuo=2;  Character(20) :: AF_om = 'zarm.omb'
      Integer :: out=3;  Character(20) :: AF_sd = 'sdcs_om.out'

      Real(8) :: Ry = 13.6056
      Real(8) :: eps_ek = 1.d-6

      Integer :: i16 = 16, is = 1, klagr = 3, ion_state=0

      Call inf_sdcs_omt

!----------------------------------------------------------------------
! ... target information:

      Call Read_aarg('AF_t',AF_t)
      Call Check_file(AF_t)
      Open(nut,file=AF_t)
      Call R_target(nut)
      Call R_channels(nut)

      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion
      if(ion.ne.0) Stop ' ion for target <> 0 '

      np = ntarg; Call Read_ipar(nut,'np',np)
      ni = ntarg; Call Read_ipar(nut,'ni',ni)

      Close(nut)

!----------------------------------------------------------------------
! ... ion information:

      Call Read_aarg('AF_i',AF_i)
      Call Check_file(AF_i)
      Open(nut,file=AF_i)
      Call R_target_ion(nut)
      Call R_channels_ion(nut)
      Close(nut)

     
!----------------------------------------------------------------------
! ... projection coefficients:

      Call Read_aarg('p',AF_p)
      Call Check_file(AF_p)
      Open(nut,file=AF_p)

      Call Read_rpar(nut,'EK',EK) 
      Call Read_ipar(nut,'nps',nps) 
      Call Read_ipar(nut,'nion',nion) 
      Call Read_ipar(nut,'nq',nq) 
 
      Allocate(dd(0:2*nq,nps,ntarg_ion)); dd=0.d0
      Allocate(qq(0:2*nq)); qq=0.d0
      Allocate(eq(0:2*nq)); eq=0.d0
      Allocate(s_ion(nion)); sion=0.d0

!----------------------------------------------------------------------
! ... define input parameters:

      Call Read_rarg('eps_ek',eps_ek)
      Call Read_iarg('i16',i16)
      Call Read_iarg('is',is)
      Call Read_iarg('ion_state',ion_state)
      if(is.gt.ni) then 
        write(*,*) 'is =', is, ' - initial atomic states'
        write(*,*) 'ni =', ni, ' - max. atomic states in zarm.om file'
        Stop 'Stop: is > ni'
      end if

!----------------------------------------------------------------------
! ... index of ionization psedostates:

      iion = 0
      Do it=1,ntarg
       if(Etarg(it).gt.Etarg_ion(1)) Exit
       iion = it
      End do

      Call Read_iarg('iion',iion)
      write(*,*) 'iion = ',iion,'  -  number  of bound pseudostates'

!----------------------------------------------------------------------
! ... read omegas for given energy and initial state (is=1):

      ntr = ntarg*(ntarg+1)/2
      Allocate( y(ntr), om(ntarg))
      om = 0.d0

      Call Read_aarg('om',AF_om)
      open(nuo,file=AF_om,status='OLD')
      rewind(nuo)
    5 read(nuo,*,end=15) x,ns   
      read(nuo,*) (y(i),i=1,ns) 

      if(abs(x-EK).gt.eps_ek) go to 5

      Do js = is+1,ntarg
       itr = Index_TR(ion,is,js,np,ni)
       if(itr.eq.0) Cycle
       if(itr.gt.ns) Exit
       om(js)=y(itr)
      End do

   15 Continue

      S_om = SUM(om(iion+1:))
      if(S_om.eq.0.d0)  Stop 'no data in ZARM.OMB for given energy'

! ... transfer omega to cross-sections

      g=iabs(IStarg(is))*(2.0*Ltarg(is)+1)
      if(IStarg(is).eq.0) g=Ltarg(is)+1
      es=EK - (etarg(is)-etarg(1))*2
      S=g/es * 3.1415926                            ! in ao^2

      if(i16.gt.0) then 
       S = S * 0.28003                              ! in 10^-16
       ii = i16-16
       if(ii.gt.0) S = S * 10**ii
      end if
      S_om = S_om * S;   SN=S

      AS = ' in a_0^3'
      if(i16.ne.0) write(AS,'(a,i2,a)') ' in 10^-',i16,' cm^2'

!----------------------------------------------------------------------
! ... define sdcs:

      Allocate(sd(0:2*nq),sdd(0:2*nq))

      E1 = etarg_ion(1)
      etarg_ion = (etarg_ion-E1)*Ry*2

! ... cycle over ion states:
      
      Do ion=1,nion

       read(nut,'(a)') AF;  write(*,*) AF

       Do jq = 1, 2*nq
        read(nut,*) iq,qq(iq),eq(iq)
        Do is=1,nps
         read(nut,*) js,nopen
         read(nut,*) (S,i=1,nopen)
         read(nut,*) (S,i=1,nopen)
         read(nut,*) (dd(iq,js,it),it=1,ntarg_ion)
        End do
       End do

       if(ion_state.gt.0.and.ion_state.ne.ion) Cycle

! ... cycle over energies of ionized electron:

      Do iq=1,2*nq
       S = 0.d0
       Do js=iion+1,nps                       
        S = S + om(js)*dd(iq,js,ion)
       End do
       Sd(iq) = S
      End do

      eq(0) =  eq(1) - (eq(2)-eq(1)) 

      sd(0) =  XLAGR(klagr,nq,eq(1),sd(1),eq(0))

      Do iq=0,2*nq; sdd(iq)=sd(iq)+sd(2*nq-iq); End do
    
      dq = (qq(2)-qq(1))   

      S_ion(ion) =  ZIS(dq,2*nq+1,SDD(0)) / 2.d0 * SN


      SD  = SD * SN  / (Ry)    !  dS/dE with E in eV
      SDD = SDD * SN  / (Ry) 

      write(AF_sd,'(a,i2.2)') 'sdcs_omt_',ion
      open(out,file=AF_sd)
      write(out,'(a,5x,a)') &
        '       E_eV         SDCS_tot       SDCS1_dir      SDCS2_exch', trim(AS)
      Do iq = 0,2*nq
       E1 = eq(iq)-etarg_ion(ion); if(iq.eq.0) E1 = 0.d0
       write(out,'(4E15.5)') E1, sdd(iq),SD(iq),SD(2*nq-iq)
      End do

      End do   !  ion states

!---------------------------------------------------------------------------------
! ... debug information:

      write(*,'(/a/)') 'total ion target sec.'
      write(*,'(20E15.5)') EK*Ry,(s_ion(it),it=1,nion)

      write(*,'(/a/)') 'accuracy estimation:  SUM(S_ion),S_om, SUM(S_ion)/S_om'
      write(*,'(2E15.5,f12.3)') SUM(S_ion),S_om, SUM(S_ion)/S_om

     End  ! UTILITY  SDCS_omt


!======================================================================
      Subroutine inf_sdcs_omt
!======================================================================
!     provide screen information about sdcs_omt utility
!----------------------------------------------------------------------
      Implicit real(8) (a-h,o-z)

      Character :: A

      Call get_command_argument(1,A)  
      if(A.ne.'?') Return

      write(*,'(a)') &
'                                                                            ',&
'     sdcs_omt:                                                              ',&
'                                                                            ',&
'     single differential cross section for electron-impact ionization       ',&
'     based on the excitation collision strengths for pseudostates and       ',&
'     their projections to the real ionic continuum                          ',&
'                                                                            ',&
'     incident and ejected electron energies defined from file "projections" ',&
'                                                                            ',&
'     INPUT FILES:                                                           ',&
'                                                                            ',&
'     target_ps   - description of e - A  scattering                         ',&
'     target_ion  - description of e - A+ scattering                         ',&
'     zarm.omb    - collision strenths for e - A problem                     ',&
'     prodections - ion-states projections for given initial energy (EK)     ',&
'                   and set of ejected energies (nq)                         ',&
'                   (obtained with program BSR_SEn)                          ',&
'                                                                            ',&
'     OUTPUT FILES:  sdcs_omt_nn (nn indicates final ionic state)            ',&
'                                                                            ',&
'     File names can be re-definded through the arguments as:                ',&
'        AF_t=...  AF_i=... om=... p=... out=...                             ',&
'                                                                            ',&
'     PARAMETERS:                                                            ',&
'                                                                            ',&
'     is [1]     -  initial atomic state                                     ',&
'     iion       -  possible redifinition of number of bound pseudostates    ',&
'     eps_ek     -  tollerance for initial electron energy                   ',&
'     i16  [16]  -  define units for resulting sdcs (i16=0 mens a.u)         ',&
'     ion_state [0]  - indicate specific final ion state;                    ',&
'                      0 means all states                                    ',&
'                                                                            '
    Stop ' '
                         
    End Subroutine inf_sdcs_omt


