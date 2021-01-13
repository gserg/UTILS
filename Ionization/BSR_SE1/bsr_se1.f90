!======================================================================
!     PROGRAM       B S R _ S E 1
!
!               C O P Y R I G H T -- 2011
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!======================================================================
!
!     Calculation of overlaps between pseudostates and  ontinuum wave
!     functions obtained in the R-matrix approach for one given energy
!     of the ejjected electron 
!
!     INPUT FILES:
!
!     target      - information about scattering system
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
!     projection1 - resulting overlaps for given incident energy   
!
!=====================================================================
!
!     OTHER PARAMETERS
!
!     EL - energy of ejected electron in eV 
!
!     ASYMPCK parameters:
!     debug - level of debug printing (0,1,2,3)
!     AC, iauto, mfgi - the parameter for asymptotic package 'ASYMPT'
!                       by Creese, CPC, 1982, ...
!                       Interface module for call ASYMPT --> ZAFACE 
!     AC   - accuracy ( < 0.01)
!     iauto - 2 
!     mfgi  - number of points in external region
!
!=====================================================================
      Use bsr_se
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)

      Call bsr_se1_inf

      Call CPU_time(t1)

! ... read arguments:
    
      Call Read_data
 
! ... loop over ion partial waves:

      Do klsp=1,nlsp

       i=Len_trim(AF_mat );  write(AF_mat (i-2:i),'(i3.3)') klsp
       i=Len_trim(AF_rsol);  write(AF_rsol(i-2:i),'(i3.3)') klsp
       i=Len_trim(AF_rovl);  write(AF_rovl(i-2:i),'(i3.3)') klsp
       i=Len_trim(AF_bnd );  write(AF_bnd (i-2:i),'(i3.3)') klsp
       i=Len_trim(AF_h   );  write(AF_h   (i-2:i),'(i3.3)') klsp

! ...  check the R-matrix overlaps:

       if(Icheck_file(AF_rovl).eq.0) Call rovl_out(AF_mat,AF_rsol,AF_rovl)

       Open(nuo,file=AF_rovl,form='UNFORMATTED')
       read(nuo) khm,nhm,kch
       if(allocated(a)) Deallocate(a,d,v,dr,di)
       Allocate(a(nhm,khm),v(nhm),d(khm),dr(kch),di(kch))
       Do i=1,khm; read(nuo) a(1:nhm,i); End do

! ...  scattering calculations: 

       Call Sub_sct(qq)

! ...  open ubound.nnn (bound.nnn option ?):

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
       if( lps(i).ne.lpar (klsp))  Cycle 
        if(sps(i).ne.ispar(klsp))  Cycle
        if(pps(i).ne.ipar (klsp))  Cycle
        if(abs(Eps(i)-E).gt.eps_e) Cycle
        is = i; Exit
       End do
       if(is.eq.0) Cycle 

! ... warning:

       if(abs(EE-E).lt.0.001) then
        write(*,'(a,i5,a)') 'Pseudostate',is,'  is too close' 
        write(*,'(2E16.8)') E,EE 
       end if

! ... ovelap vector with RM_solutions:

       d = MATMUL(v,a)

! ... final calculations of overlaps:

      Call RM_OVLF (qq,nopen,d,dr,di,F,G,RMAT,KMAT) 

      nop(is) = nopen
      ddr(1:nopen,is) = dr(1:nopen)
      ddi(1:nopen,is) = di(1:nopen)

      Do i=1,nopen
       it = iptar(klsp,i)
       wt(it,is) = wt(it,is) + dr(i)**2 + di(i)**2
      End do

      End do    ! over pseudo-states

      End do    ! over partial waves

! ... check if all pseudostates are covering:

      Do i=1,nps
       if(nop(i).ne.0) Cycle
       write(*,*) 'Psedostate ',i,'  is not covered' 
      End do

! ... output:

      open(out,file=AF_out)
      write(out,'(a,f12.5,a,f12.5)') 'EL =',EL,'     Ry =',Ry 
      write(out,'(a,i5)') 'nps =',nps 
      Do is = 1,nps
       write(out,'(2i5)') is,nop(is) 
       write(out,'(5D16.8)') ddr(1:nop(is),is)
       write(out,'(5D16.8)') ddi(1:nop(is),is)
       write(out,'(5D16.8)') wt(:,is)
      End do

      Call  CPU_time(t2)
      write(pri,'(/a,f6.2,a)')  'time = ',(t2-t1)/60,' min'
 
      End  ! program bsr_SE
 

 
!======================================================================
      Subroutine bsr_se1_inf
!======================================================================
!     provide screen information about bound_tab utility
!----------------------------------------------------------------------
      Character(80) :: A

      Call get_command_argument(1,A)  
      if(A.ne.'?') Return

      write(*,'(a)') &
'                                                                                  ',&
'     BSR_SE1 provide the overlaps between the pseudostates and continuum          ',&
'     functions for the given energy (output file projection1)                     ',&
'                                                                                  ',&
'     Call as:    BSR_SE1 EL=...;   EL - electron energy in eV                     ',&
'                                                                                  ',&
'     Calculations should be done in the folder of BSR calculations  for the       ',&
'     ejected electron                                                             ',&
'                                                                                  '
      Stop ' '                                                                    
                                                                                   
      End Subroutine bsr_se1_inf                                                 

