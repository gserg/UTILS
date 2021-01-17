!======================================================================
!     PROGRAM       B S R _ S E n
!
!               C O P Y R I G H T -- 2012
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!======================================================================
!
!     Calculation of overlaps between pseudostates and  continuum wave
!     functions obtained in the R-matrix approach
!
!     INPUT FILES:
!
!     target_ps   - information about psedostates ( e + atom )
!
!     target      - information about ionized system ( e + ion )
!
!     bsr_mat.nnn - overlap matrix for partial wave 'nnn'
!                   (after program bsr_mat)
!
!     h.nnn       - standard H.DAT file for partial wave 'nnn' 
!     rsol.nnn    - R-matrix solutions for partial wave 'nnn'
!                   (after program bsr_hd with itype=1)
!
!     ubound.nnn  - bound-like solutions
!                   (after program bsr_hd with itype=-1)
!
!     OUTPUT FILES:
!
!     projections - resulting overlaps for given scattering energy   
!                   and a set of ejected energies (1:2*nq)
!                   for each final target state
!=====================================================================
!
!     ARGUMRNTS:
!
!     EK - energy of scattering electron in Ry 
!     nq - number of ejected energies
!
!     ASYMPCK parameters (optional):
!
!     debug - level of debug printing (0,1,2,3)
!     AC, iauto, mfgi - the parameter for asymptotic package 'ASYMPT'
!                       by Creese, CPC, 1982, ...
!                       Interface module for call ASYMPT --> ZAFACE 
!     AC   - accuracy ( < 0.01)
!     iauto - 2 
!     mfgi  - number of points in external region
!
!     CALL AS:  bsr_sen EK=7.35 nq=20
!
!=====================================================================
      Use bsr_se
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)

      Call CPU_time(t1)

      Call bsr_sen_inf

! ... read arguments:
    
      Call Read_data
 
! ... output basic parameter: 

      open(out,file=AF_out)
      write(out,'(a,f10.5,a,f12.2,a)') 'EK =',EK,' Ry  ',E0,' eV' 
      write(out,'(a,i5,a)') 'nps  =',nps,' - number of psedostates' 
      write(out,'(a,i5,a)') 'nion =',nopen_ion, &
                            ' - number of open ion states'  
      write(out,'(a,i5,a)') 'nq   =',nq,' - number of energy poins'  

      Do it = 1,nopen_ion

      write(out,'(a,i5,a)') 'targ =',it,' - final ion state'  

! ... set up the energy grid for ejected electron:

      dq = (EK - Eion(it))/nq/2                

      Do iq=1,2*nq  !  we escape zero energy

      qq = etarg(it) + dq*iq;  EL = qq*Ry

      write(out,'(i3,E15.5,f10.5,a)') iq,qq,EL, '    -> iq,qq,EL '

      Call Sub1(qq)

      Do is = 1,nps
       if(nop(is).eq.0) nop(is)=1
       write(out,'(2i5,a )') is,nop(is),'  p.s., nopen:  dr,di' 
       write(out,'(5D16.8)') ddr(1:nop(is),is)
       write(out,'(5D16.8)') ddi(1:nop(is),is)
       write(out,'(5D16.8)') wt(1:ntarg,is)
      End do

      End do  !  iq

      End do  !  it

      Call CPU_time(t2)
      write(pri,'(/a,f6.2,a)')  'time = ',(t2-t1)/60,' min'

      END  ! program bsr_SE

 
!======================================================================
      Subroutine Sub1(qq)
!======================================================================
! ... calculations for one ejected energy qq
!----------------------------------------------------------------------
      Use bsr_se
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)

      ddr = 0.d0; ddi = 0.d0; wt = 0.d0; nop = 0

      EE = E1 + qq/2  

! ... loop over ion partial waves:

      Do klsp=1,nlsp

       i=Len_trim(AF_mat ); write(AF_mat (i-2:i),'(i3.3)') klsp
       i=Len_trim(AF_rsol); write(AF_rsol(i-2:i),'(i3.3)') klsp
       i=Len_trim(AF_rovl); write(AF_rovl(i-2:i),'(i3.3)') klsp
       i=Len_trim(AF_bnd ); write(AF_bnd (i-2:i),'(i3.3)') klsp
       i=Len_trim(AF_h   ); write(AF_h   (i-2:i),'(i3.3)') klsp

! ...  check the R-matrix overlaps:

       if(Icheck_file(AF_rovl).eq.0) Call rovl_out(AF_mat,AF_rsol,AF_rovl)

       Open(nuo,file=AF_rovl,form='UNFORMATTED')
       read(nuo) khm,nhm,kch
       if(allocated(a)) Deallocate(a,d,v,dr,di)
       Allocate(a(nhm,khm),v(nhm),d(khm),dr(kch),di(kch))
       Do i=1,khm; read(nuo) a(1:nhm,i); End do

! ...  scattering calculations: 

       Call Sub_sct(qq)
       if(nopen.eq.0)  Cycle

! ...  open ubound.nnn (bound.nnn ?):

       Open(nub,file=AF_bnd,form='UNFORMATTED')
       read(nub) ns_b,kch_b,kcp_b,nhm_b,nbound
       if(nhm_b.ne.nhm) Stop 'BSR_SE: different nhm in bound.nnn'

! ... cycle over pseudostates:

       Do ib=1,nbound
        read(nub) jb, labl; read(nub) E; read(nub) v
        if(E.gt.Eps(nps)+eps_e) Exit

! ... check if the bound state belong to pseudostates:

       is = 0
       Do i=1,nps
        if(lps(i).ne.lpar (klsp))  Cycle 
        if(sps(i).ne.ispar(klsp))  Cycle
        if(pps(i).ne.ipar (klsp))  Cycle
        if(abs(Eps(i)-E).gt.eps_e) Cycle
        is = i; Exit
       End do
       if(is.eq.0) Cycle 

! ... warning:

       if(abs(EE-E).lt.0.001) write(*,'(a,i5,a,F10.6)') &
        'Pseudostate',is,'  is too close:',abs(EE-E) 

! ... ovelap vector with RM_solutions:

       d = MATMUL(v,a)

! ... final calculations of overlaps:

      Call RM_OVLF (qq,nopen,d,dr,di,F,G,RMAT,KMAT) 

      nop(is) = nopen
      ddr(1:nopen,is) = dr(1:nopen)
      ddi(1:nopen,is) = di(1:nopen)

      wt(:,is) = 0.d0
      Do i=1,nopen
       it = iptar(klsp,i)
       wt(it,is) = wt(it,is) + dr(i)**2 + di(i)**2
      End do

      End do    ! over pseudo-states

      End do    ! over partial waves

! ... check if all pseudostates are covering:

      Do i=1,nps
       if(nop(i).ne.0) Cycle
       write(*,*) 'Psedostate ',i,'  is not covered',Eps(i) 
      End do

      End Subroutine Sub1


!======================================================================
      Subroutine bsr_sen_inf
!======================================================================
!     provide screen information about bsr_sen utility
!----------------------------------------------------------------------
      Character(80) :: A

      Call get_command_argument(1,A)  
      if(A.ne.'?') Return

      write(*,'(a)') &
'                                                                                  ',&
'     BSR_SE2 provide the overlaps between the pseudostates and continuum          ',&
'     functions for the given incident energy (output file projections),           ',&
'     including both direct and exchange cases, for set of ejected energies        ',&
'                                                                                  ',&
'     Call as:    bsr_sen  EK=..  EL=.. itarg=..  out=..                           ',&
'                                                                                  ',&
'     EK  -  energy of incident electron in Ry                                     ',&
'     nq  -  number of ejected energies equaly devided all posssible interval      ',&
'     itarg - index of ionic state                                                 ',&
'     out   - name of output file [projection2]                                    ',&
'                                                                                  ',&
'     Calculations should be done in the folder of BSR calculations  for the       ',&
'     ejected electron - many files are using, see remarks in bsr_sen.f90          ',&
'                                                                                  '
      Stop ' '                                                                    
                                                                                   
      End Subroutine bsr_sen_inf                                                 


