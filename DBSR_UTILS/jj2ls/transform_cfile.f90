!======================================================================
      Subroutine transform_cfile
!======================================================================
      Use jj2ls
      Use conf_jj,      ncfg_jj  => ncfg
      Use conf_LS,      only: ncfg_LS  => ncfg

      Implicit none 

      Integer :: ic
      Integer, allocatable :: ipt(:)
      Real(8) :: CM, CN1, CN2

      Allocate(C1(ncfg_jj),C2(ncfg_LS))
      C1 = WC(1:ncfg_jj)
      CN1 = SUM(C1*C1); if(CN1.eq.0.d0) Return

      Do ic = 1,ncfg_LS                                               
       C2(ic) = SUM(C1(:)*C_trans(:,ic))
      End do

      CN2 = SUM(C2*C2)
      write(pri,'(/a,a,f12.5,a,f12.5)') &
       'solution in c-file:','  JJ_norm =',CN1,'  LS_norm =',CN2

      Allocate(ipt(ncfg_ls))
      Call SORTA(ncfg_LS,C2,IPT) 
      Call write_LS_conf_c(nuc,C2,ipt)

      Deallocate(C1,C2,ipt)
      AF = trim(name)//'.c'
      write(pri,'(/a,a,a/)') 'file ',trim(AF),' with coefficients is recorded'
	  
      End Subroutine transform_cfile 
