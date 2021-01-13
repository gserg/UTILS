!=======================================================================
!     merging the set of w-files with choice of orbitals and 
!     with optional changing the spectroscopic notation
! 
!          1.w + 2.w + 3.w + ... --> res.w 
!
!=======================================================================

      Use radial, E => Zk, EK => Yk

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      Character EL3*3, ans*3, EL4*4
      Character(3), External :: ELF3
      Character(4), External :: ELF4
      Character(40) :: AF
      Integer(4) :: nu=1

      Z = 1.d0; Call INITR
      Call Alloc_radial(irf)
      i=1

      iarg = IARGC()
      if(iarg.gt.0) Call GETARG(1,AF)

      if(AF.eq.'?') then
        write(*,'(/a)') 'w123 merges the set of w-files with choice of orbitals'
        write(*,'(/a)') '     1.w + 2.w + 3.w + ... --> res.w'
        write(*,'(/a)') 'program has interactive input/output' 
        write(*,'(/a)') 'Call as:  w123  or  w123  1.w '
        write(*,'(/a)') 'Results:  res.w'
        Stop ' '
      end if        

      if(iarg.gt.0)  go to 11
   10 write(*,*) 'Enter file-name for w-file  or  end: '
      read(*,'(a)') AF
   11 ii=LEN_trim(AF)
      if(AF(1:3).eq.'end') go to 3
      if(ii.eq.0) go to 10
      Call OpenF(nu,AF,'UNFORMATTED','OLD')

    1 READ (nu,end=2) ATOM,TERM,EL3,mro(i),Z,E(I),EK(I),AZ(I), &
                      (P(J,I),J=1,mro(i))
      write(*,'(a8,a8,a6)') ATOM,TERM,EL3
      write(*,'(a)') ' new EL ?  (d|EL,a3): '
      read(*,'(a)') ans 
      if(ans(1:1).eq.'d') go to 1
      ii=LEN_TRIM(ans)
      if(ii.gt.1) then
       Call EL3_nlk(ans,n1,l1,k1); ero(i)=ELF4(n1,l1,k1)
      else
       Call EL3_nlk(EL3,n1,l1,k1); ero(i)=ELF4(n1,l1,k1)
      end if

      Do j=1,i-1
       if(ero(j).eq.ero(i)) then; i=i-1; Exit; end if
      End do
      i=i+1

      if(i.gt.mrf) Call Alloc_radial(mrf+irf)
      go to 1

    2 Close(nu)
      go to 10

    3 nrf=i-1
      write(*,'(a,i3,a)')  'Now there is ', nrf, '  w.f.:'
      write(*,'(15(1x,a4))') (ero(i),i=1,nrf)
      write(*,*) 'Enter file-name for result w-file : '
      read(*,'(a)') AF
      Call OpenF(nu,AF,'UNFORMATTED','UNKNOWN')

      DO i=1,nrf
       Call EL4_nlk(ero(i),n1,l1,k1); el3=ELF3(n1,l1,k1)
       Write(nu) ATOM,TERM,EL3,mro(i),Z,E(I),EK(I),AZ(I), &
                (P(J,I),J=1,mro(i))
      End do
      Close(nu)

      END  ! utility w123


!-----------------------------------------------------------------------
      Subroutine OpenF (nu,AF,AFOR,AST)
!-----------------------------------------------------------------------


      Character(*) AF,AFOR,AST

    1 OPEN(nu,file=AF,form=AFOR,status=AST,err=2)
      Return

    2 write(*,*) 'There is no file  ',TRIM(AF)
      write(*,*) 'Try agin: enter file-name  or  end :'
      read(*,*) AF
      if(AF(1:3).eq.'end') Stop
      go to 1

      End Subroutine OpenF

