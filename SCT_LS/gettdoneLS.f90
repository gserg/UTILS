!=====================================================================
!     UTILITY  gettdoneLS
!=====================================================================
!     get 'tmat.done' file from 'zarm.tma' or 'zarm.tmb' 
!
!     gettdoneLS itr=.. itr1=.. itr2=.. ek=.. ek2=.. eps_e=..
!---------------------------------------------------------------------
      Use target; Use channels; Use tm_LS

      Implicit real(8) (A-H,O-Z)

      Real(8), allocatable :: matr(:),mati(:) 

      Character(40) :: AF_tma = 'zarm.tma'
      Character(40) :: AF_tmb = 'zarm.tmb'
      Character(40) :: AF_tdone = 'tmat.done'

      Integer :: ST1,ST2,LT1,LT2,IP1,IP2

      Call Inf_sub

!----------------------------------------------------------------------
! ... define input parameters if any:

      itr =0; Call Read_iarg('itr' ,itr )
      if(itr.eq.0) then 
       itr1=1; Call Read_iarg('itr1',itr1)
       itr2=1; Call Read_iarg('itr2',itr2)
      else
       itr1=itr; itr2=itr
      end if
      if(itr2.lt.itr1) itr2=itr1

      ekk=0.d0; Call Read_rarg('ek' ,ekk)
      if(ekk.le.0.d0) then
       ek1=0.d0; Call Read_rarg('ek1',ek1)
       ek2=0.d0; Call Read_rarg('ek2',ek2)
      else
       ek1=ekk; ek2=ekk
      end if
      if(ek2.lt.ek1) ek2=ek1

      eps_e = 1.d-7; Call Read_rarg('eps_e',eps_e)

      ifano=0; Call Read_iarg('ifano',ifano)

!----------------------------------------------------------------------
! ... target and channels information:

      nut=1; Open(nut,file='target',status='OLD')
      Call R_target  (nut)
      Call R_channels(nut)
      np = ntarg; Call Read_ipar(nut,'np',np)
      ni = ntarg; Call Read_ipar(nut,'ni',ni)
      Close(nut)

      ion = nz-nelc;  zion=1.d0; if(ion.gt.1) zion=ion*ion

!----------------------------------------------------------------------
! ... transition parameters:

      E1=(ETARG(itr1)-ETARG(1))*2
      E2=(ETARG(itr2)-ETARG(1))*2
      LT1 = ltarg(itr1)
      LT2 = ltarg(itr2)
      ST1 = istarg(itr1)
      ST2 = istarg(itr2)
      ip1=(-iptarg(itr1)+1)/2
      ip2=(-iptarg(itr2)+1)/2                   

! ... allocations:

      mtr = mch*(mch+1)/2;  Allocate(matr(mtr),mati(mtr))

!----------------------------------------------------------------------
! ... define input file and parameters:

      in = 1 
      if(Icheck_file(AF_tmb).eq.1) then
       Open(in,file=AF_tmb,status='OLD')
       istyle = 1
      elseif(Icheck_file(AF_tma).eq.1) then
       Open(in,file=AF_tma,status='OLD')
       istyle = 0
      else
       Stop 'Can not find T-matrix file'  
      end if

!----------------------------------------------------------------------
! ... define T-matrix elements:

      L_max = 0; ee_max = 0.d0

      rewind(in)
   1  if(istyle.eq.0) then
       read(in,*,end=2) ee,nopen,ntr,ii
       read(in,*) (matr(i),mati(i),i=1,ntr)
      else
       read(in,*,end=2) ee,nopen,kopen,ii,i1,i2,nj
       if(kopen.gt.nopen)     Stop ' kopen > nopen'
       ntr = kopen*(kopen+1)/2
       read(in,*) (matr(i),mati(i),i=1,ntr)
       ktr = (nopen-kopen)*nj
       if(ktr.gt.0) read(in,*) (matr(i),mati(i),i=ntr+1,ntr+ktr)
      end if

      if(ek1.gt.0.d0.and.ee.lt.ek1) go to 1
      if(ek2.gt.0.d0.and.ee.gt.ek2) go to 1

      IL=lpar(ii); IS=ispar(ii)

      if(L_max.lt.IL) L_max = IL
      if(ee_max.lt.ee) ee_max = ee

      Do i = 1,nch(ii); if(iptar(ii,i).ne.itr2) Cycle
      Do j = 1,i;       if(iptar(ii,j).ne.itr1) Cycle


       if(itr1.ne.itr2) then
        ktr=ich2*(ich2-1)/2+ich1
       else
        i1=min(ich1,ich2); i2 = max(ich1,ich2); itr=i1*(i1-1)/2+i2
       end if

       ktr = (i-1)*i/2+j
       if(istyle.eq.1.and.i.gt.kopen) &
        ktr = kopen*(kopen+1)/2 + nj*(i-kopen-1) + j

        lch1=lch(ii,j)
        lch2=lch(ii,i)  

        if(ifano.gt.0) then 
         Call Fano_tma(lch1,lch2,matr(ktr),mati(ktr),ar,ai)
        else
         ar = matr(ktr)
         ai = mati(ktr)
        end if
        Call Add_tm_LS(ee,ar,ai,IL,IS,lch1,lch2,0)

        if(itr1.eq.itr2.and.i.ne.j) &
        Call Add_tm_LS(ee,ar,ai,IL,IS,lch2,lch1,0)        
      End do; End do

      go to 1
    2 Close(in)

!----------------------------------------------------------------------
! ... output:

      ip = 0
      Call Def_ip_LS(L_max,ee_max)


      if(Icheck_file('T_extra').ne.0) then         
       open(in,file='T_extra')
       Do k=1,kcase
        read(in,*) ee,IS,i0,j1,j2,t1,t2,jj,s1,s2
        Do jj=j1,j2
         lch1 = jj + jp(2,k)
         lch2 = jj + jp(3,k)
         IL = jj
         IS = jp(1,k)                
         t1=t1*s1; t2=t2*s2
         Call Add_tm_LS(ee,t1,t2,IL,IS,lch1,lch2,k)
        End do
       End do
      end if

      Call Sort_tm_LS

      iout = 3; Open(iout,file=AF_tdone)
      ne=1
      Do i=2,ndata; ee = EK(i)
       m=1
       Do j=1,i-1;if(abs(ee-EK(j)).lt.1.d-7) then; m=0; Exit; end if; End do
       ne = ne + m
      End do

      kk = ndata/ne
      write(*,*) 'ndata,ne*kk,ne,kk =',ndata,ne*kk,ne,kk 

      write(iout,'(8I5,2D15.7,a)') itr1,itr2,LT1,LT2,ST1,ST2,ndata,ne,E1,E2, &
                             '     itr1,itr2,LT1,LT2,ST1,ST2,ndata,ne,E1,E2'

      Do i=1,ndata

       a=1.d0; b=1.d0
       if(i.gt.1) then
        if(ip(i-1).eq.ip(i)) then
         a=TR(i)/TR(i-1)  
         b=TI(i)/TI(i-1)  
        end if
       end if
       write(iout,'(E15.7,4I5,2E15.7,i3,2f10.5)')  &
                    EK(i),STpar(i),LTpar(i), &
                    l1(i),l2(i),TR(i),TI(i), &
                    ip(i),a,b
      End do

      Close(iout)

      End  ! UTILITY  gettdoneLS

 

!=========================================================================
      Subroutine Def_ip_LS(ll,ee)
!=========================================================================
!     define the different pattens for the given total L and energy ee;
!     kcase - number of pattens; jp(kcase) - patten itself  
!     then define the ip-pointer for other records (ip arrays)
!-------------------------------------------------------------------------
      Use tm_LS

      Implicit none
      Integer :: i,j,k, ll
      Real(8) :: ee

! ... define number of cases:

      k=0
      Do i=1,ndata
       if(LTpar(i).ne.ll)  Cycle
       if(EK(i).ne.ee) Cycle
       k=k+1
      End do
      kcase = k   

! ... define case pattens

      Allocate(jp(3,k))
      k=0
      Do i=1,ndata
       if(LTpar(i).ne.ll)  Cycle
       if(EK(i).ne.ee) Cycle
       k=k+1
       jp(1,k) = STpar(i)
       jp(2,k) = l1(i) - ll
       jp(3,k) = l2(i) - ll
      End do

! ... define case index:

      ip = 0
      Do i=1,ndata
       if(LTpar(i).lt.5) Cycle

       Do j=1,k
        if(STpar(i).ne.jp(1,j)) Cycle
        if(l1(i)-LTpar(i).ne.jp(2,j)) Cycle
        if(l2(i)-LTpar(i).ne.jp(3,j)) Cycle
        ip(i)=j
        Exit
       End do

!       if(ip(i).eq.0) Stop 'Unknown case'

      End do

      End Subroutine Def_ip_LS



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
'     gettdoneLS utility generates "tmat.done" file in LS-coupling      ',&
'     based on the T-matrix elements given in "zarm.tma(b)" file        ',&
'     obtained with BSR in the LS-coupling                              ',&
'                                                                       ',&          
'          zarm.tma (zarm.tmb) + target   ->    tmat.done               ',&
'                                                                       ',&
'     Call as: gettdoneLS  itr1=.. itr2=.. ek1=.. ek2-.. eps_e=..       ',&
'                                                                       ',&
'     itr1, itr2  -  transition itr1 -> itr2                            ',&
'     ek1, ek2    -  range of electron energy                           ',&
'     eps_e       -  energy tolerance                                   ',&
'                                                                       '
      Stop ' '

      End Subroutine inf_sub




