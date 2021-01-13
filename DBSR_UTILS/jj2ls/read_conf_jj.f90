!======================================================================
      Subroutine read_conf_jj(muc,kshift,job,check)
!======================================================================
!     read and add configurations to the list "conf_jj"
!     job  =  'add'     -  just add
!          =  'detect'  -  return -ic if exist and 
!          =   others   -  add if not exist 
!     check = 'check'   -  check the configurations for number of 
!                          electrons and parity 
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Character(*), intent(in) :: job, check 
      Integer, intent(in) :: muc,kshift
      Integer, external :: Jadd_cfg_jj
      Integer :: nuc,i,ic

      nuc=iabs(muc); if(muc.gt.0) rewind(nuc)
      if(check.eq.'check') then; ne=0; parity=0; end if
    1 read(nuc,'(a)',end=2) CONFIG
      if(CONFIG(1:3).eq.'*') go to 2
      if(CONFIG(6:6).ne.'(') go to 1
      read(nuc,'(a)') SHELLJ
      read(nuc,'(5x,a)') INTRAJ
      Call Decode_cj
      in = in + kshift
      if(check.eq.'check') Call Test_cj
      ic = Jadd_cfg_jj(job)

      if(ic.lt.0) Stop 'Read_conf_jj: repeated states?'
      WC(ic)=0.d0
      i=INDEX(CONFIG,')',BACK=.TRUE.)+1
      if(LEN_TRIM(CONFIG(i:)).ne.0) Read(CONFIG(i:),*) WC(ic)
      go to 1
    2 Continue

      End Subroutine read_conf_jj


!======================================================================
      Subroutine read_config_jj(nuc)
!======================================================================
!     Read only configurations from GRASP c-file (unit 'nuc')
!     (no weights)
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer, intent(in) :: nuc
      Integer, external :: Jadd_cfg_jj
      Integer :: i

      SHELLJ =' '
      INTRAJ =' '
      Jshell=0; Vshell=0; Jintra=0

      rewind(nuc); ne=0; parity=0
    1 read(nuc,'(a)',end=2) CONFIG
      if(CONFIG(1:1).eq.'*') go to 2
      if(CONFIG(6:6).ne.'(') go to 1
      Call Decode_confj
      i = Jadd_cfg_jj('find')
      go to 1
    2 Continue

      End Subroutine read_config_jj


!======================================================================
      Integer Function Jadd_cfg_jj(job)
!======================================================================
!     add new or detect existing CAS in conf_jj list;
!     returns position of given CAS in the list
!     job  =  'add'     -  just add
!          =  'detect'  -  return -ic if exist
!          =   others   -  add if not exist 
!
!     it is special version to jj2ls program
!
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Character(*), intent(in) :: job 
      Integer :: i,ic,ip
      Integer, external :: Ifind_jjorb, Iadd_symc, Iadd_symt

      Jadd_cfg_jj = 0
      if(no.le.0) Return

      if(mcfg.eq.0) Call Alloc_cfg(icfg)
      Jtotal = Jintra(no)

      Do i=2,no
       if(nn(i-1).eq.nn(i).and.ln(i-1).eq.ln(i)) iq(i) = -iq(i)         ! !!!
      End do

      iconf = Iadd_symc(Jtotal,no,iq,kn)

      iterm = Iadd_symt(iconf,no,Jshell,Vshell,Jintra)

      Do i = 1,no; np(i)=Ifind_jjorb(nn(i),kn(i),in(i),2); End do

! ... check if we already have such state:

      if(job.ne.'add') then
       Do ic = 1,ncfg
        if(IS_term(ic).ne.iterm) Cycle
        ip = ip_state(ic); Jadd_cfg_jj = ic; if(job.eq.'detect') Jadd_cfg_jj=-ic 
        Do i = 1,no; ip=ip+1
         if(np(i).ne.IP_orb(ip)) then; Jadd_cfg_jj=0; Exit; end if
        End do
        if(Jadd_cfg_jj.ne.0) Return
       End do
      end if

      ncfg=ncfg+1
      if(ncfg.eq.mcfg.or.lcfg+no.gt.kcfg) Call Alloc_cfg(mcfg+icfg)
      IS_term(ncfg)=iterm
      ip_state(ncfg)=lcfg
      Do i=1,no; lcfg=lcfg+1; IP_orb(lcfg)=np(i); End do
      Jadd_cfg_jj = ncfg

      End Function Jadd_cfg_jj


