!=====================================================================
!      target_jj --> target 
!---------------------------------------------------------------------
      Use target    
      Use channels
                                                                         
      Use target_jj, ntarg_jj => ntarg, nelc_jj => nelc, nz_jj => nz,  &
        nct_jj => nct, nwt_jj => nwt, AFT_jj => AFT, BFT_jj => BFT,    &
        etarg_jj => etarg, jtarg_jj => jtarg, ptarg_jj => ptarg,       &  
        nctarg_jj => nctarg, nwtarg_jj => nwtarg, ictarg_jj => ictarg

      Use channels_jj, nlsp_jj => nlsp, mch_jj => mch, jpar_jj => jpar, &
        ipar_jj => ipar, nch_jj => nch, iptar_jj => iptar, kch_jj => kch,& 
        lch_jj => lch, jch_jj => jch, ipch_jj => ipch, ipconf_jj => ipconf, &
        ELC_jj => ELC, AFP_jj => AFP, BFP_jj => BFP, ncp_jj => ncp,     &
        nwp_jj => nwp

      Implicit none

      Integer :: inp=1;   Character(80) :: t_jj    = 'target_jj'
      Integer :: out=2;   Character(80) :: t_LS    = 'target'

      Integer :: i,j,n
      Character(4), external :: ELF4
      Character :: AF

      Call get_command_argument(1,AF)  
      if(AF.eq.'?') then    !  help section 
       write(*,*)
       write(*,*) 'targ_jj converts target_jj(DBSR) in target(BSR) formats'
       write(*,*)
       write(*,*) 'Call as:   targ_jj'
       write(*,*)
       Stop ' '
      end if

      Call Check_file(t_jj)
      Open(inp,file=t_jj) 
      Call Read_target_jj(inp)
      Call Read_channels_jj(inp)
      Close(inp)

      Call allocate_target(ntarg_jj)

      nelc  = nelc_jj
      nz    = nz_jj
      nct   = nct_jj
      nwt   = nwt_jj
      AFT   = AFT_jj
      BFT   = BFT_jj
      etarg = etarg_jj
      ltarg = jtarg_jj
      istarg = 0
      iptarg = ptarg_jj
      COUPLING  = 'JJ'
      nctarg = nctarg_jj
      nwtarg = nwtarg_jj

      Open(out,file=t_LS) 
      Call Write_target_LS(out)

      nlsp = nlsp_jj
      mch  = mch_jj
      Call Allocate_channels

      ispar = 0
      lpar  = jpar_jj
      ipar  = ipar_jj
      AFP   = AFP_jj     
      BFP   = BFP_jj
      ncp   = ncp_jj
      nwp   = nwp_jj

      nch   = nch_jj
      lch   = lch_jj
      iptar = iptar_jj
      ipconf = ipconf_jj
      jkch   = jch_jj
      
      max_nc = 0
      max_wf = 0
      Do i=1,nlsp;  n=nch(i)
       if(ipconf(i,n).gt.max_nc) max_nc=ipconf(i,n)
       Do j=1,n
        ELC(i,j) = ELF4(107,lch(i,j),0)      
       End do
      End do

      Call Write_channels_LS (out,1,max_nc,max_wf)

      Call R_target(out)


      End ! utilty targ_jj

