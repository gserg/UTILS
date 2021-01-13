!---------------------------------------------------------------------
!     name.frm  -->  name.w
!---------------------------------------------------------------------

      IMPLICIT REAL*8(A-H,O-Z)

      CHARACTER AT*6,TT*6,EL*3,NAME*20,AF*24
      DIMENSION PT(220)

      iarg = IARGC()

      if(iarg.gt.0) Call GETARG(1,name)
      if(iarg.le.0.or.name.eq.'?') then
        write(*,'(/a)') 'frm_w  converts file  name.frm  to name.w'
        write(*,'(/a)') 'w - unformatted, frm - formated  MCHF (CFF) formats'
        write(*,'(/a)') 'Call as:  frm_w  name.frm '
        write(*,'(/a)') 'Results:  name.w'
        Stop ' '
      end if        
            
      ii = INDEX(name,'.',BACK=.TRUE.)
      if(name(ii+1:ii+3).ne.'frm') Stop  ' Extension in input file shoud be ".frm" '

      in=1
      AF = name
      Open(in,file=AF,status='OLD')

      iout=2
      AF = name(1:ii)//'w'
      Open(iout,file=AF,form='UNFORMATTED')


    1 read(in,'(A6,A6,3X,A3,I6,F6.2,3(1PE18.10))',end=2)  &
                   AT,TT,EL,M,ZT,ETI,EKI,AZI
      read(in,'(4(1PE18.10))') (PT(J),J=1,M)
      write(iout) AT,TT,EL,M,ZT,ETI,EKI,AZI,(PT(J),J=1,M) 
      go to 1
    2 Continue

      End
