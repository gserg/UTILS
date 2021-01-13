!=======================================================================
      Subroutine Read_data
!=======================================================================

      USE bsr_se
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)
 
      Integer :: i
      Integer, external :: Icheck_file

!----------------------------------------------------------------------
! ... read arguments:
    
      Open(pri,file=AF_log)
 
      if(icheck_file(AF_inp).ne.0) then
       Open(inp,file=AF_inp)

       Call Read_ipar(inp,'iauto',iauto)
       Call Read_ipar(inp,'mfgi' ,mfgi )
       Call Read_rpar(inp,'ac'   ,AC   )
       Call Read_rpar(inp,'delta',delta)
       Call Read_ipar(inp,'debug',debug)

       Call Read_rpar(inp,'EK',EK)
       Call Read_ipar(inp,'nq',nq)
       Call Read_apar(inp,'ps',AF_ps)
       Call Read_apar(inp,'out',AF_out)

      end if

      Call Read_iarg('iauto',iauto)
      Call Read_iarg('mfgi' ,mfgi )
      Call Read_rarg('ac'   ,AC   )
      Call Read_rarg('delta',delta)
      Call Read_iarg('debug',debug)

      Call Read_rarg('EK',EK)
      Call Read_iarg('nq',nq)
      Call Read_aarg('ps',AF_ps)
      Call Read_aarg('out',AF_out)

      if(EK.le.0.d0) &
       Stop 'provide EK=... for incident electron (in eV)'
      if(nq.le.0) &
       Stop 'provide nq=... - number of energies for ejected electron'

!----------------------------------------------------------------------
! ... pseudo-state data:

      Call Check_file(AF_ps)
      open(nut,file=AF_ps)
      Call R_target(nut)
      Close(nut)

      ET = etarg(1) + EK/2
      Do i=1,ntarg
       if(etarg(i).gt.ET) Exit; nps=i                   
      End do       
      Allocate(Eps(nps),lps(nps),sps(nps),pps(nps))
      Eps = etarg(1:nps)
      lps=ltarg(1:nps); sps=istarg(1:nps); pps = iptarg(1:nps)    

!----------------------------------------------------------------------
! ... e + ion scattering parameters: 

      Open(nut,file=AF_tar,status='OLD')
      Call Allocate_target(0)
      Call R_target (nut)
      Call R_channels(nut)
      Close(nut)
      Allocate(Eion(ntarg));  Eion=Etarg
      E1 = Etarg(1)
      Do i = 1,ntarg; Etarg(i) = 2.d0*(Etarg(i)-E1); End do

!----------------------------------------------------------------------
! ... set up energies:

      Z = nz;  AWT = 0.d0

!      Call Conv_au (Z,AWT,au_cm,au_eV,0)   ! ???
!      Ry = au_eV/2.d0

      E0 = EK*Ry
      write(pri,'(/a,f10.5,a)') 'EK   = ',EK,' - incident energy in Ry' 
      write(pri,'(/a,f10.5,a)') 'E0   = ',E0,' - incident energy in eV' 

      Eion_Ry = (E1 - Eps(1))*2
      Eion_eV = (E1 - Eps(1))*2*Ry
      if(Eion_Ry.lt.0) Stop 'E_ionization < 0'
      write(pri,'(/a,f10.5,a)') 'Eion = ',Eion_eV,' - ionization energy in eV' 

      Emax = E0-Eion_eV
      write(pri,'(/a,f10.5,a)') 'Emax = ',Emax,' - max. ejected energy in eV' 
      if(Emax.lt.0) Stop 'E_ejected < 0'

      write(pri,'(/a,i10,a)')   'nps  = ',nps,' - number of pseudostates' 
      write(pri,'(/a,i10,a)')   'nq   = ',nq,' - number of energies (E2)' 

! ... allocations for projections (at one energy):

      Allocate(ddr(mch,nps),ddi(mch,nps),nop(nps),wt(ntarg,nps))
      
! ... find open target states: 

      nopen_ion = 0  
      Do i=1,ntarg
       Eion(i) = (Eion(i) - Eps(1))*2
       if(EK-Eion(i).gt.0.d0) nopen_ion=i
      End do

      End Subroutine Read_data

