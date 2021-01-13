!======================================================================
      Subroutine get_tr_coefs
!======================================================================
!     find transformation coeffocients for each angular symmetry in row
!----------------------------------------------------------------------
      Use jj2ls, only: n,JJ_term,nmom,ncup,JJ_coupling,LS_coupling
      Use conf_jj
      Use symt_list, only: nsymt

      Implicit none

      Do iterm=1,nsymt

       Call Get_symt(iterm,iconf,no,Jshell,Vshell,Jintra)
       Jtotal = Jintra(no)
       Call Get_symc(iconf,Jtotal,no,nn,kn,ln,jn,iq,in)

       JJ_term = iterm

       Call Check_jj_shells

       Call Make_couplings
       if(n.gt.1) Call RECUP(nmom,ncup,JJ_coupling,LS_coupling)

       Call Recup_jj_shells

      End do
 
      End  Subroutine get_tr_coefs

