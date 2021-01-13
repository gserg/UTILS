!=====================================================================
!
!     utility     s l a t e r _ w
!
!                 C O P Y R I G H T -- 2005
!
!     Written by:   Oleg Zatsarinny
!
!======================================================================
!     Transfer slater orbitals to MCHF format:   name.slw  --> name.w
!
!     Slater orbitals are supposed to be described as in R_matrix code:
!
!        SUM (J=1:NCO) [ C(J) * r ^ IRAD(J) * exp(-ZE(J)) ]
!
!======================================================================
!
!     File name.slw is supposed has the structure:
!
!     Z, atom, term
!     nwf
!                       list:  (1:nwf)
!     n,l,k
!     NCO
!     IRAD(1:NCO)
!     ZE(1:NCO)
!     C(1:NCO)
!
!======================================================================

      Use radial

      IMPLICIT REAL(8) (A-H,O-Z)

      Real(8), parameter :: Emax=100.d0
      Real(8), parameter :: EPS_p=1.d-16

      Integer(4), Allocatable, Dimension(:) :: IRAD
      Real(8), Allocatable, Dimension(:) ::  ZE,C

      Character(40) :: AF
      Character(3), External :: ELF3
      Character(4), External :: ELF4
      Character(3) :: EL3

      Integer(4) :: in = 8
      Integer(4) :: iout = 9

! ... input name:

      iarg = IARGC()

      if(iarg.gt.0) Call GETARG (1,AF)
      if(iarg.le.0.or.AF.eq.'?') then
        write(*,'(/a)') 'slater_w  converts file  name.slw  to name.w'
        write(*,'(/a)') 'w - unformatted  MCHF (CFF) format'
        write(*,'(/a)') 'slw - slater radial wave functions (CIV3 format)'
        write(*,'(/a)') '      (see code for description of slw-format)'
        write(*,'(/a)') 'Call as:  slater_w  name.slw '
        write(*,'(/a)') 'Results:  name.w'
        Stop ' '
      end if        


      Call Check_file(AF)
      open(in,file=AF,status='OLD')

      read(in,*) Z,atom,term  
      read(in,*) nn

      Call INITR
      Call Alloc_radial(nn)

      Do iw = 1,nn

        read(in,*) nro(iw),lro(iw),kro(iw)
        nwf = iw
        read(in,*) NCO
        if(Allocated(IRAD)) Deallocate(IRAD,ZE,C)
        Allocate(IRAD(nco),ZE(nco),C(nco))
        read(in,*) IRAD(1:NCO)
        read(in,*) ZE(1:NCO)
        read(in,*) C(1:NCO)

        ero(iw) = ELF4(nro(iw),lro(iw),kro(iw))

        Do j = 1,NCO
         SS = (2.d0*ZE(j))**(IRAD(j)+0.5d0)
         Do i=2,2*IRAD(j)
          S=dble(i); S=sqrt(S); SS=SS/S
         End do
         C(j)=C(j)*SS
        End do

        Do i = 1,nr
         S = R(i);  PR = 0.d0
         Do j = 1,NCO
          SS = ZE(j)*S; if(SS.gt.Emax) Cycle
          PR = PR + C(j) * S**IRAD(j) * exp(-SS)
         End do
         PR=PR/R2(i); if(abs(PR).lt.EPS_p) PR = 0.d0; P(i,iw)=PR
        End do

        Do i = NR,1,-1
         if(abs(P(i,iw)).eq.0.d0) Cycle;  mro(iw)=i; Exit
        End do

        S = 0.d0
        Do j = 1,NCO
         if(IRAD(j).ne.lro(iw)+1) Cycle;  S = S + C(j)
        End do
        AZ(iw) = S

        m = iw
        mexp(m) = lro(m) + 1
        aexp(m) = - D1/mexp(m)
        ZR = Z*R(1)
        bexp(m) = P(1,m)/(AZ(m)*R2(1)*R(1)**lro(m))
        bexp(m) = (bexp(m) - D1 + ZR/(lro(m)+1) )/ZR**2

        PN = QUADR(iw,iw,0)
         write(*,'(1x,a,a4,a,a4,a,F10.6)') &
          '<',ero(iw),'|', ero(iw),'> = ',PN

        Do jw = 1,iw-1
         if(lro(iw).ne.lro(jw)) Cycle
         S = QUADR(iw,jw,0)
         write(*,'(1x,a,a4,a,a4,a,F10.6)') &
          '<',ero(iw),'|', ero(jw),'> = ',S
        End do

       End do  !  over input orbitals

! ... output results:

       i=INDEX(AF,'.',BACK=.TRUE.); AF=AF(1:i)//'w'
       Open(iout,file=AF,form='UNFORMATTED')

       Do i = 1,nwf
        EI = -HL(i,i)/2.d0
        ZETA = 0.d0;  if(lro(i).gt.0) ZETA = Z / QUADR(i,i,-3)
        el3 = ELF3(nro(i),lro(i),kro(i))
        write(iout) Atom,Term,el3,mro(i),Z,EI,ZETA,AZ(i),P(1:mro(i),i)
       End do


       End     !  utility  slater_w

