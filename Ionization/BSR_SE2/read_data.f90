!=======================================================================
      Subroutine Read_data
!=======================================================================

      USE bsr_se
      Use target
      Use channels

      Implicit none
 
      Integer :: i, ntarg_ps
      Integer, external :: Icheck_file
      Real(8) :: z, AWT

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
       Call Read_rpar(inp,'EL',EL)
       Call Read_apar(inp,'AF_ps',AF_ps)

      end if

      Call Read_iarg('iauto',iauto)
      Call Read_iarg('mfgi' ,mfgi )
      Call Read_rarg('ac'   ,AC   )
      Call Read_rarg('delta',delta)
      Call Read_iarg('debug',debug)

      Call Read_rarg('EK',EK)
      Call Read_rarg('EL',EL)
      Call Read_aarg('AF_ps',AF_ps)
      Call Read_iarg('itarg',itarg)
      Call Read_aarg('out',AF_out)

! ... pseudo-state data:

      Call Check_file(AF_ps)
      open(nut,file=AF_ps)
      Call R_target(nut)
      Close(nut)

      nps = ntarg; ntarg_ps=ntarg
      Allocate(Eps(nps),lps(nps),sps(nps),pps(nps))
      Eps = etarg(1:nps)
      lps=ltarg(1:nps); sps=istarg(1:nps); pps = iptarg(1:nps)    

! ... e + ion scattering parameters: 

      Open(nut,file=AF_tar,status='OLD')
      Call Allocate_target(0)
      Call R_target (nut)
      Call R_channels(nut)
      Close(nut)

! ... set up energies:

      Z = nz;  AWT = 0.d0
!      Call Conv_au (Z,AWT,au_cm,au_eV,0)   ! ???
!      Ry = au_eV/2.d0

      Eion_Ry = (Etarg(itarg) - Eps(1))*2
      Eion_eV = (Etarg(itarg) - Eps(1))*2*Ry

      if(EL.lt.0.d0) then
       ET = etarg(itarg) - EL /Ry
       EK = (ET - eps(1)) * 2
       EL = -EL
      else
       ET = eps(1) + EK/2
      end if

      Do i=1,ntarg_ps
       if(eps(i).gt.ET) Exit; nps=i                   
      End do     

      E0 = EK*Ry

      qqq = EL/Ry
      if(EL.eq.0.d0) then    ! equal energy case
       qqq = (EK - Eion_Ry) / 2
       EL = qqq * Ry
      end if

      qq2 = qqq + (Etarg(itarg)-Etarg(1))*2 
      EL2 = qq2*Ry

      qq1 = EK - Eion_Ry - qqq + (Etarg(itarg)-Etarg(1))*2
      EL1 = qq1*Ry

      write(pri,'(/a,f10.5,a)') 'EK   = ',EK, ' - incident energy in Ry' 
      write(pri,'(/a,f10.5,a)') 'E0   = ',E0, ' - incident energy in eV' 

      write(pri,'(/a,f10.5,a)') 'EL   = ',EL, ' - ejected energy in eV' 
      write(pri,'(/a,f10.5,a)') 'Eion = ',Eion_eV,' - ionization energy in eV' 

      write(pri,'(/a,i10,a)')   'nps  = ',nps,' - number of pseudostates' 

      if(qqq+Eion_Ry.ge.EK) then
       write(pri,'(/a)') 'Ejected energy is too big'
       Stop 'Ejected energy is too big'
      end if

! ... allocations for final projections:

      Allocate(ddr(mch,nps),ddi(mch,nps),wt(ntarg,nps),nop(nps))
      
      End Subroutine Read_data

