!======================================================================
!     PROGRAM       B S R _ S E 2
!
!               C O P Y R I G H T -- 2011
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
!     projection2  - resulting overlaps for given ejected-electron energy   
!                    and scattering-electron energy
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
 
      USE bsr_se
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)
 
      Real(8) :: t1,t2

      Call CPU_time(t1)

! ... read arguments:
    
      Call Read_data
 
      open(out,file=AF_out)
      write(out,'(a,f10.5,f12.5)') 'EK =',EK,E0 
      write(out,'(a,i5)') 'nps =',nps 
      write(out,'(a,i5)') 'ion =',itarg 

! ... ejected electron:

      Call Sub1(qq2)

      write(out,'(a,f10.5)') 'EL2 =',EL2 
      Do is = 1,nps
       if(nop(is).eq.0) nop(is)=1
       write(out,'(2i5)') is,nop(is) 
       write(out,'(5D16.8)') ddr(1:nop(is),is)
       write(out,'(5D16.8)') ddi(1:nop(is),is)
      End do

! ... scattering electron:

      Call Sub1(qq1)

      write(out,'(a,f10.5)') 'EL1 =',EL1 
      Do is = 1,nps
       if(nop(is).eq.0) nop(is)=1
       write(out,'(2i5)') is,nop(is) 
       write(out,'(5D16.8)') ddr(1:nop(is),is)
       write(out,'(5D16.8)') ddi(1:nop(is),is)
      End do

      Call CPU_time(t2)
      write(pri,'(/a,f6.2,a)')  'time = ',(t2-t1)/60,' min'

      END  ! program bsr_SE

 
!======================================================================
      Subroutine Sub1(qq)
!======================================================================
!     calculations for electron energy qq 
!     (in Ry, relative to ground state)
!----------------------------------------------------------------------
      USE bsr_se
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)

      ddr = 0.d0; ddi = 0.d0; wt = 0.d0; nop = 0

      EE = E1 + qq/2  ! total energy

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

      nop(is) = nopen

      Call RM_OVLF (qq,nopen,d,dr,di,F,G,RMAT,KMAT) 

      ddr(1:nopen,is) = dr(1:nopen)
      ddi(1:nopen,is) = di(1:nopen)

      Do i=1,nopen
       it = iptar(klsp,i)
       wt(it,is) = wt(it,is) + dr(i)**2 + di(i)**2
      End do

      End do    ! over pseudo-states

      End do    ! over partial waves

! ... check if all pseudostates are covering:

      ii = 0
      Do i=1,nps
       if(nop(i).eq.0) ii=ii+1
      End do
      write(*,'(i5,a)') ii,' psedostates are not covered' 

      End Subroutine Sub1
