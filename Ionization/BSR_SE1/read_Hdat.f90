!======================================================================
      Module h_nnn
!======================================================================
      Implicit none

      Integer :: nelc          ! number of target electrons
      Integer :: nz            ! nuclear charge
      Integer :: lm            ! max. small l
      Integer :: km            ! max. multipole
      Integer :: ion           ! ionization index 

      Integer :: ntarg         ! # of target states
      Integer, allocatable :: Ltarg(:),IStarg(:)
      Real(8), allocatable :: Etarg(:)

      Real(8) :: RA,RB         ! R-matrix radius

      Integer :: IL2,IS2,IP2   ! total term
 
      Integer :: nch           ! # of channels
      Integer :: nhm           ! # of solutions

      Integer, allocatable :: NCONAT(:)
      Integer, allocatable :: LCH(:)
      Real(8), allocatable :: CF(:,:,:)  
      Real(8), allocatable :: VALUE(:)  
      Real(8), allocatable :: WMAT(:,:)  

      End Module h_nnn


!======================================================================
      Subroutine Read_h_nnn(AF_h)
!----------------------------------------------------------------------
!     read h.nnn file for partial wave nnn = klsp
!----------------------------------------------------------------------      

      Use h_nnn

      Character(*) :: AF_h
      Integer :: i, nuh
      Real(8) :: C,E1
      
      Call Check_file(AF_h)
      Call Find_free_unit(nuh) 
      Open(nuh,file=AF_h,form='UNFORMATTED')

      read(nuh) nelc,nz,lm,km,ntarg,RA,RB
      ion = nz - nelc
 
      if(allocated(Etarg)) Deallocate(Etarg,Ltarg,IStarg)
      Allocate(Etarg(ntarg),Ltarg(ntarg),IStarg(ntarg))

      read (nuh) Etarg
      E1  = Etarg(1)
      DO i = 1,ntarg; Etarg(i)=2.d0*(Etarg(i)-E1); END DO

      read (nuh) Ltarg
      read (nuh) IStarg;  IStarg = iabs(IStarg)  
 
! ... skip BUTTLE CORRECTION    ( don't use here! )
 
      READ (nuh) ((C,i=1,3),j=1,LM)
 
! ... look for the required partial wave

      READ (nuh) IL2,IS2,IP2,NCH,NHM
      IP2 = (-1)**IP2
                       
      if(allocated(NCONAT)) Deallocate(NCONAT,LCH,CF,VALUE,WMAT)
      Allocate(NCONAT(ntarg),LCH(nch),CF(nch,nch,km),VALUE(nhm), &
               WMAT(nch,nhm))

      READ (nuh) (NCONAT(i),i=1,ntarg)
      READ (nuh) (LCH(i),i=1,nch)
      READ (nuh) (((CF(i,j,k),i=1,nch),j=1,nch),k=1,km)

      READ (nuh) (VALUE(i),i=nhm,1,-1)
      DO i = 1,NHM; VALUE(i) = 2.0D0*(VALUE(i)-E1); END DO

      READ (nuh) ((WMAT(i,j),i=1,nch),j=nhm,1,-1)

      Close(nuh)
  
      End Subroutine Read_h_nnn

