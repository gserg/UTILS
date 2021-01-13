!=======================================================================
      Subroutine Read_data
!=======================================================================
      Use bsr_se
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

       Call Read_rpar(inp,'EL',EL)
       Call Read_apar(inp,'AF_ps',AF_ps)

      end if

      Call Read_iarg('iauto',iauto)
      Call Read_iarg('mfgi' ,mfgi )
      Call Read_rarg('ac'   ,AC   )
      Call Read_rarg('delta',delta)
      Call Read_iarg('debug',debug)

      Call Read_rarg('EL',EL)
      Call Read_aarg('AF_ps',AF_ps)

! ... pseudo-state data:

      Call Check_file(AF_ps)
      open(nut,file=AF_ps)
      Call R_target(nut)
      Close(nut)

      nps=ntarg
      Allocate(Eps(nps),lps(nps),sps(nps),pps(nps))
      Eps = etarg; lps=ltarg; sps=istarg; pps = iptarg    

! ... e + ion scattering parameters: 

      Open(nut,file=AF_tar,status='OLD')
      Call Allocate_target(0)
      Call R_target (nut)
      Call R_channels(nut)
      Close(nut)
      E1 = Etarg(1)
      Do i = 1,ntarg; Etarg(i) = 2.d0*(Etarg(i)-E1); End do

! ... set up energies:

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

      EE = E1 + EL/Ry/2;  qq = EL/Ry

      Allocate(ddr(mch,nps),ddi(mch,nps),wt(ntarg,nps),nop(nps))
      ddr = 0.d0; ddi = 0.d0; wt = 0.d0; nop = 0
      
      End Subroutine Read_data
