!=====================================================================
!     UTILITY  gettdonejk
!=====================================================================
!     get 'tmat.done' file from 'zarm.tma' 
!---------------------------------------------------------------------

      Use target; Use channels; Use tm_list

      Implicit real(8) (A-H,O-Z)

      Real(8), allocatable :: matr(:),mati(:) 

      Real(8) ::  eps_e = 1.d-7

      Call Inf_sub

!----------------------------------------------------------------------
! ... target information:

! ... read channel information:

      nut=1; Open(nut,file='target',status='OLD')
      Call R_target (nut)
      Call R_channels(nut)
      Close(nut)

      ion = nz-nelc;  zion=1.d0; if(ion.gt.1) zion=ion*ion

! ... define transition and energy under consideration:

      itr1=1; Call Read_iarg('itr1',itr1)
      itr2=1; Call Read_iarg('itr2',itr2)

      ek1=0.d0; Call Read_rarg('ek1',ek1)
      ek2=0.d0; Call Read_rarg('ek2',ek2)

      Call Read_rarg('eps_e',eps_e)

!----------------------------------------------------------------------
! ... define input file and parameters:

      in = 1; Open(in,file='zarm.tma',status='OLD')

      E1=(ETARG(itr1)-ETARG(1))*2; E2=(ETARG(itr2)-ETARG(1))*2;
      JT1 = jtarg(itr1)-1; JT2 = jtarg(itr2)-1
      ip1=0; ip2=0

! ... allocations:

      mtr = mch*(mch+1)/2;  Allocate(matr(mtr),mati(mtr))

!----------------------------------------------------------------------
! ... define T-matrix elements:

      jj_max = 0; ee_max = 0.d0

      rewind(in)
    1 read(in,*,end=2) ee,nopen,ntr,ii
      read(in,*) (matr(i),mati(i),i=1,ntr)
      write(*,*) ii
      if(ek1.gt.0.d0.and.ee.lt.ek1) go to 1
      if(ek2.gt.0.d0.and.ee.gt.ek2) go to 1

      jj = jpar(ii)-1; ipp = (-ipar(ii)+1)/2

      if(jj_max.lt.jj) jj_max = jj
      if(ee_max.lt.ee) ee_max = ee

      Do i = 1,nch(ii); if(iptar(ii,i).ne.itr2) Cycle
      Do j = 1,i;       if(iptar(ii,j).ne.itr1) Cycle

        ktr = (i-1)*i/2+j; if(ktr.gt.ntr) Stop 'ktr > ntr'

        lch1=lch(ii,j); lch2=lch(ii,i)  
        kk1= jkch(ii,j)-1; kk2=jkch(ii,i)-1

        Call Fano_tma(lch1,lch2,matr(ktr),mati(ktr),ar,ai)

        Call Add_tm_list(ee,ar,ai,jj,ipp,lch1,lch2,kk1,kk2)

        if(itr1.eq.itr2.and.i.ne.j) &
        Call Add_tm_list(ee,ar,ai,jj,ipp,lch2,lch1,kk2,kk1)

      End do; End do

      go to 1
    2 Close(in)

!----------------------------------------------------------------------
! ... output:

      Call Def_io(jj_max,ee_max)

      Call Sort_tm_list

!======================================================================

      if(Icheck_file('T_extra').ne.0) then         
       open(in,file='T_extra')
       Do k=1,kcase
        read(in,*) ee,j1,j2,i1,i2,i3,i4,t1,t2,jj,s1,s2

        Do jj=j1,j2,2
         lch1 = (jj + jo(2,k))/2
         lch2 = (jj + jo(3,k))/2
         kk1  = jj + jo(4,k)
         kk2  = jj + jo(5,k)
         t1=t1*s1; t2=t2*s2
         iip = mod(lch1,2)
         Call Add_tm_list(ee,t1,t2,jj,ipp,lch1,lch2,kk1,kk2)
        End do
       End do
      end if

!      Call Sort_tm_list

!======================================================================

      iout = 3; Open(iout,file='tmat.done')
      ne=1
      Do i=2,ndata; ee = EK(i)
       m=1
       Do j=1,i-1;if(abs(ee-EK(j)).lt.1.d-7) then; m=0; Exit; end if; End do
       ne = ne + m
      End do

      kk = ndata/ne
      write(*,*) 'ndata,ne*kk,ne,kk =',ndata,ne*kk,ne,kk 

      write(iout,'(8I5,2D15.7)') itr1,itr2,JT1,JT2,ip1,ip2,ndata,ne,E1,E2

      Do i=1,ndata

       a=1.d0; b=1.d0
       if(i.gt.1) then
        if(io(i-1).eq.io(i)) then
         a=TR(i)/TR(i-1)  
         b=TI(i)/TI(i-1)  
        end if
       end if
       write(iout,'(D15.7,6I5,2D15.7,i3,2f10.5)')  &
                    EK(i),jjpar(i),ip(i), &
                    l1(i),jk1(i),l2(i),jk2(i),TR(i),TI(i), &
                    io(i),a,b
      End do

      Close(iout)

      End  ! UTILITY  gettdonejk




!=========================================================================
      Subroutine Def_io(jj,ee)
!=========================================================================
! ... define parameters different symmetries in T-matrix elements
! ... for partial wave with total JJ and energy ee  
!-------------------------------------------------------------------------

      Use tm_list

      Implicit none
      Integer :: i,j,k, jj
      Real(8) :: ee

! ... define number of cases:

      k=0
      Do i=1,ndata
       if(jjpar(i).ne.jj)  Cycle
       if(EK(i).ne.ee) Cycle
       k=k+1
      End do
      kcase = k   

! ... define case pattens

      Allocate(jo(5,k))
      k=0
      Do i=1,ndata
       if(jjpar(i).ne.jj)  Cycle
       if(EK(i).ne.ee) Cycle
       k=k+1
       jo(1,k) = ip(i)
       jo(2,k) = 2*l1(i) - jj
       jo(3,k) = 2*l2(i) - jj
       jo(4,k) = jk1(i) - jj
       jo(5,k) = jk2(i) - jj
      End do

! ... define case index:

      io = 0
      Do i=1,ndata

       if(jjpar(i).lt.2) Cycle
       Do j=1,k
!        if(ip(i).ne.jo(1,j)) Cycle
        if(2*l1(i)-jjpar(i).ne.jo(2,j)) Cycle
        if(2*l2(i)-jjpar(i).ne.jo(3,j)) Cycle
        if(jk1(i)-jjpar(i).ne.jo(4,j)) Cycle
        if(jk2(i)-jjpar(i).ne.jo(5,j)) Cycle
        io(i)=j
        Exit
       End do
       if(io(i).eq.0) Stop 'Unknown case'

      End do

      End Subroutine Def_io


!======================================================================
      Subroutine inf_sub
!======================================================================
!     provide screen information about add_farm utility
!----------------------------------------------------------------------
       
      Character :: A

      iarg = command_argument_count()
      if(iarg.eq.0) Return
      Call GET_COMMAND_ARGUMENT(1,A)        
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                       ',&
'     gettdonejk utility generates "tmat.done" file in JK-coupling      ',&
'     based on the T-matrix elements given in "zarm.tma" file           ',&
'     obtained with BSR in the JK-coupling, to be able to use them      ',&
'     in the program MJK                                                ',&
'                                                                       ',&                                                                           '                                                                       ',&
'              zarm.tma + target   ->    tmat.done                      ',&
'                                                                       ',&
'     Call as: gettdonejk  itr1=.. itr2=.. ek1=.. ek2-.. eps_e=..       ',&
'                                                                       ',&
'     itr1, itr2  -  transition itr1 -> itr2                            ',&
'     ek1, ek2    -  range of electron energy                           ',&
'     eps_e       -  energy tolerance                                   ',&
'                                                                       '
      Stop ' '

      End Subroutine inf_sub



