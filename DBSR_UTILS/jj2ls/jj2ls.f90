!======================================================================
!     utility     j j 2 l s
!
!                 C O P Y R I G H T -- 2009
!
!     Written by:   Oleg Zatsarinny
!
!======================================================================
!
!     transform name.c and/or name.j in jj-coupling to 
!     name_LS.c and name_LS.j in LS-coupling
!
!     Call as:  jj2ls  name     
!
!----------------------------------------------------------------------
      Use jj2ls

      Use symc_list,    only: JJ_nsymc => nsymc
      Use symc_list_LS, only: LS_nsymc => nsymc

      Use symt_list,    only: JJ_nsymt => nsymt
      Use symt_list_LS, only: LS_nsymt => nsymt

      Use conf_jj,      only: ncfg_jj  => ncfg
      Use conf_LS,      only: ncfg_LS  => ncfg, core

      Implicit real(8) (A-H,O-Z)
      Character :: AS*80

!----------------------------------------------------------------------
! ... files

      Call GET_COMMAND_ARGUMENT(1,name)

      if(name.eq.'?'.or.command_argument_count().lt.1)  &
      Stop 'Call as: jj2ls name, or name.c, or.name.j'

      i=LEN_TRIM(name)
      if(name(i-1:i).eq.'.c'.or.name(i-1:i).eq.'.j') name(i-1:i)='  '

      AF = trim(name)//'.c'
      Call Check_file(AF)
      Open(nuc,file=AF)
      
      AF = trim(name)//'_LS.log'
      Open(pri,file=AF)

      AF = trim(name)//'.scr'
      Open(nua,file=AF,form='UNFORMATTED')

      Call Read_iarg('debug',debug)

      write(pri,'(a,a)') 'name of case:   ',trim(name)

!----------------------------------------------------------------------
! ... read jj-configurations:

      Call Read_core_jj(nuc)

      Call Read_conf_jj(nuc,0,'add','check');  JJ_nterms = JJ_nsymt

      AF = trim(name)//'.c'
      write(pri,'(/a,a/)') trim(AF),' contains:'

      write(pri,'(a,i5)') 'number of atomic states   = ',ncfg_jj
      write(pri,'(a,i5)') 'number of configurations  = ',JJ_nsymc
      write(pri,'(a,i5)') 'number of ang. symmetries = ',JJ_nsymt

! ... define all possible LS-terms:	    
	  
      Call get_LS_terms;    LS_nterms = LS_nsymt

      write(pri,'(/a/)') 'corresponding LS values:' 

      write(pri,'(a,i5)') 'number of configurations  = ',LS_nsymc
      write(pri,'(a,i5)') 'number of ang. symmetries = ',LS_nsymt

! ... define the recoupling coefficiens between jj and LS terms: 

      Allocate(C_term(JJ_nterms,LS_nterms));  C_term =  zero

      Call get_tr_coefs

      ncoef = 0 
      Do i=1,jj_nterms; Do j=1,ls_nterms
       if(C_term(i,j).ne.zero) ncoef=ncoef+1
      end do; end do

      write(pri,'(/a/)') 'get_tr_coef provide:'
      write(pri,'(a,i8)') 'total number of transition coefficients = ', &
                           jj_nterms*ls_nterms
      S = ncoef; S = S/(jj_nterms*ls_nterms)*100
      write(pri,'(a,i8,f8.1,a)') 'non-zero transition coefficients = ',ncoef, S,' %'

! ... define the LS configurations:

      Call get_LS_conf

      write(pri,'(/a/)') 'get_LS_conf provide:'
      write(pri,'(a,i5)') 'number of LS atomic states   = ',ncfg_LS

! ... record LS configurations:

      AF = trim(name)//'_LS.c'
      open(nuc,file=AF)
      Call get_LS_core(core)
      Call write_LS_conf(nuc,core)
      Call R_term(nuc)

      write(pri,'(/a,a,a/)') 'file ',trim(AF),' is recorded'

! ... read transformation matrix from scratch file:

      Deallocate(C_term)
      Allocate (C_trans(ncfg_jj,ncfg_LS)); C_trans = zero
      rewind(nua)
    1 read(nua,end=2) i,j,C; C_trans(i,j)=C; go to 1
    2 Close(nua,status='DELETE') 
      
      write(pri,'(/a/)') 'file for transformation coefficients is created'
      write(pri,*) 

! ... transform solutions in j-file or c-file:

      AF = trim(name)//'.j'

      if(Icheck_file(AF).eq.1) then
       Open(nuj,file=AF)
       Call transform_jfile
       write(pri,'(/a,a,a/)') 'file ',trim(AF),' is recorded'
      else
       Call transform_cfile
      end if

      End  !  utility-program   " j j 2 l s "

