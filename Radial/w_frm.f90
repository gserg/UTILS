!---------------------------------------------------------------------
!
!      name.w  -->  name.frm
!
!      change the unformatted representation to formatted
!
!---------------------------------------------------------------------

      IMPLICIT REAL*8(A-H,O-Z)

      CHARACTER AT*6,TT*6,EL*3,NAME*20,AF*24
      DIMENSION PT(220)

      iarg = IARGC()

      if(iarg.gt.0) Call GETARG(1,name)
      if(iarg.le.0.or.name.eq.'?') then
        write(*,'(/a)') 'w_frm  converts file  name.w  to name.frm'
        write(*,'(/a)') 'w - unformatted, frm - formated  MCHF (CFF) formats'
        write(*,'(/a)') 'Call as:  w_frm  name.w '
        write(*,'(/a)') 'Results:  name.frm'
        Stop ' '
      end if        


      ii = INDEX(name,'.',BACK=.TRUE.)

      in=1
      AF = name(1:ii)//'w'
      Open(in,file=AF,form='UNFORMATTED',status='OLD')

      iout=2
      AF = name(1:ii)//'frm'
      Open(iout,file=AF)

    1 read(in,end=2) AT,TT,EL,M,ZT,ETI,EKI,AZI,(PT(J),J=1,M) 
      if(m.eq.0) go to 2
      WRITE(iout,'(A6,A6,3X,A3,I6,F6.2,3(1PE18.10))')  &
                   AT,TT,EL,M,ZT,ETI,EKI,AZI
      WRITE(iout,'(4(1PE18.10))') (PT(J),J=1,M)
      go to 1
    2 Continue

      End
