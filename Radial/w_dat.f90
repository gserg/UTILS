!========================================================================
!      name.w   -->   name.dat  - file for plotting (e.q., for ORIGIN)
!
!      INPUT ARGUMENTS:  name.w  
!
!========================================================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

! ... radial grid parameters:

      Real(8), parameter :: RHO = -4.D0   ! initial value of logarifmic grid
      Real(8), parameter :: H = 1./16.D0  ! step of logarifmic grid
      Integer, parameter :: NR = 220      ! max. number of radial points
      Real(8) :: Z  = 0.d0                ! charge of nuclear
      Real(8) :: R(nr), R2(nr)

! ... functions:

      Real(8), Allocatable :: P(:,:)
      Character(6) :: atom, term
      Character(3) :: EL3
      Character(3), Allocatable :: ELN(:)

      Character(80) AF
      Integer :: in=1, out=2

      iarg = IARGC();  if(iarg.gt.0) Call GETARG(1,AF)

      if(iarg.le.0.or.AF.eq.'?') then
        write(*,'(/a)') 'w_dat converts file  name.w  to name.dat'
        write(*,'(/a)') 'w - unformatted  MCHF (CFF) format'
        write(*,'(/a)') 'dat - column representation of orbitals, ready for plotting'
        write(*,'(/a)') 'Call as:  w_dat  name.w  '
        write(*,'(/a)') 'Results:  name.dat'
        Stop ' '
      end if        

      OPEN(in,file=AF,form='UNFORMATTED',status='OLD')
    
! ... define number of orbitals:

      nrf=0;  mm=0
   1  read(in,end=2) atom,term,el3,m,Z
      if(m.eq.0) go to 2
      nrf=nrf+1
      if(m.gt.mm) mm=m
      go to 1
   2  Continue
      if(nrf.eq.0) Stop 'nrf=0 -> nothing to export ! '
      if(Z.le.0.d0) Stop ' Z <= 0 '

! ... define Z-dependent radial grid:

      DO I=1,NR
       R(I)=EXP(RHO+(I-1)*H)/Z
       R2(I)=SQRT(R(I))
      END DO

! ... read and renormalize radial functions:

      Allocate(P(nr,nrf), eln(nrf));  P=0.d0
      rewind(in)
      Do i=1,nrf
       read(in,end=2) atom,term,eln(i),m,Z,E,EK,AZ, (P(j,i),j=1,m)
       P(:,i) = P(:,i) * R2(:) 
      End do

! ... output in "column" format:

      i = INDEX(AF,'.',BACK=.TRUE.); if(i.eq.0) i=LEN_TRIM(AF)  
      AF = AF(1:i)//'dat'

      OPEN(out,file=AF)
      write(out,'(i3,a)') nrf,'  - number of radial functions'

      write(out,'(6x,a1,4x,50(5x,a3,5x))') 'R',(eln(i),i=1,nrf)

      Do j=1,mm
       write (out,'(50D13.5)') R(j),(P(J,I),I=1,nrf)
      End do

      END  ! utility w_dat

