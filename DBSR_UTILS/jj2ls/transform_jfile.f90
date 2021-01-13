!======================================================================
      Subroutine transform_jfile
!======================================================================
      Use jj2ls
      Use conf_jj,      ncfg_jj  => ncfg
      Use conf_LS,      only: ncfg_LS  => ncfg

      Implicit none 

      Integer :: i,j, jj,k,m, ic,ic1,ic2, is,js, nc, nsol,jsol
      Integer, allocatable :: jsol_index(:), jsol_found(:)
      Real(8) :: CM, CN1, CN2, E, gvJ,gvLS
      Integer, external :: Ifind_position, Ipointer

      Call Read_ipar(nuj,'ncfg',nc) 
      if(nc.ne.ncfg_jj) then 
       write(*,*) 'nc_j, ncfg_jj =', nc, ncfg_jj
       Stop ' ncfg_jj in j-file inconsistent'
      end if
      Call Read_ipar(nuj,'nsol',nsol); if(nsol.eq.0) Return 

      i=Ifind_position(nuj,'Solutions');  read(nuj,*)  

      AF = trim(name)//'_LS.j'
      open(nul,file=AF)

      Call Jdef_ne(nuc)
      write (nul,'(8X,A,F5.1,A,I3,A,I7)' ) &
        '  Z = ',one ,'  NEL = ',  ne, '   NCFG = ',ncfg_LS

      Allocate(jsol_index(nsol),jsol_found(nsol))

      jsol_index = 0

      Call Read_iarr('jsol',nsol,jsol_index)       

      jsol=0; if(jsol_index(1).ne.0) jsol=1
      Allocate(C1(ncfg_jj),C2(ncfg_LS))

      Call Def_jblocks

      Do j = 1,njbl; Jtotal = JJc(j)
       i=Ifind_position(nuj,'Solutions');  read(nuj,*)

      m = 0
      Do is = 1,nsol
       read(nuj,*) i
       read(nuj,*) E,jj,ic1,ic2
       read(nuj,*) C1(ic1:ic2)

       if(jj.ne.Jtotal) Cycle      
       if(jsol.ne.0) then
        if(Ipointer(nsol,jsol_index,i).eq.0) Cycle
       end if
       m = m + 1;  jsol_found(m)=i
      End do

      write (nul,'(//A8,I4,2X,A8,I4)') '  2*J = ',Jtotal,'NUMBER =',m

      i=Ifind_position(nuj,'Solutions');  read(nuj,*)

      js = 0
      Do is=1,nsol
       C1 = zero
       read(nuj,*) i
       read(nuj,*) E,jj,ic1,ic2
       read(nuj,*) C1(ic1:ic2)
       if(Ipointer(nsol,jsol_found,i).eq.0) Cycle
       js = js + 1
       CN1 = SUM(C1*C1)
       k=0; CM = zero
       Do ic = 1,ncfg_LS
        C2(ic) = SUM(C1(:)*C_trans(:,ic))
        if(abs(C2(ic)).gt.CM) then; k=ic; CM=abs(C2(ic)); end if         
       End do
       Call Get_cfg_LS(k)
       Call Label_c (AS,1,0)
       Call g_factor(Jtotal+1,k,ncfg_LS,C2,gvJ,gvLS) 
       write(nul,'(3x,a,f15.10,4x,a,f15.10,2x,a,f15.10)') &
        'Ssms=',0.d0,'g_J=',gvj,'g_JLS=',gvls
       write(nul,'(i6,f16.8,3x,a)') js,E,trim(AS)
       write(nul,'(7f11.8)') C2
       CN2 = SUM(C2*C2)

       write(pri,'(/a,i5,2x,a,f12.5,10x,a,f12.5/)') &
        'solution:',is,'  JJ_norm =',CN1,'  LS_norm =',CN2
     
       Call Compare(pri,ncfg_jj,C1,ncfg_LS,C2)

      End do ! over is

      End do ! over j-blocks

      write(nul,'(a)') '***'
	  
      End Subroutine transform_jfile 
