!======================================================================
!     UTILITY           M E R G E                  
!
!               C O P Y R I G H T -- 2011
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!======================================================================
!
!  Merge c(w)-files in one file with consistent set indexes for 
!  one-electron rdial functions
! 
!  In order to define the same orbitals, the following criteria are 
!  used:  
!              | <p1|p2> - 1 |  <  eps 
!              | <p1|r|p1> - <p2|r|p2> | <  eps 
!              | <p1|1/r|p1> - <p2|1/r|p2> | <  eps 
!
!======================================================================
!
!  INPUT FILES:
!
!     name.c   -  set of c-files parameters of calculation
!     name.w   -  set of corresponding w-files 
!
!  OUTPUT FILES:
!
!     merge.c     -  resulting list of CSFs
!     merge.w     -  common list of radial faunctions
!     merge.log   -  running information
!     
!  INPUT DATA are defined as arguments in command line:
!
!     merge  name1 name2 ... merge=... eps=...  jmin=... jmax=...
!
!
!  Default values:  merge=merge; eps=1.d-7;  jmin=jmax=-1 (N/A)
!
!======================================================================
     
      Use conf_LS
      Use orb_LS
      Use radial, n => nro, l => lro, ks => kro, mx => mro, EL => ero
     
      Implicit real(8) (A-H,O-Z)

      Logical :: EX
      Character( 3) :: EL3
      Character( 4) :: EL4
      Character( 3), external :: ELF3
      Character( 4), external :: ELF4

! ... default value for parameter: 

      Real(8) :: eps = 1.d-6
      Integer :: jmin = -1
      Integer :: jmax = -1
      Integer :: icheck = 0

! ... files:
      
      Integer, parameter :: ma=80, mt=1000
      Character(ma) :: AF,AS, merge = 'merge',  files(mt)

      Integer :: inp= 5;  Character(ma) :: AFi    = 'merge.inp'
      Integer :: pri= 6;  Character(ma) :: AF_log = 'merge.log'
      Integer :: nuc=11;  Character(ma) :: AFC
      Integer :: nuw=12;  Character(ma) :: AFW
      Integer :: muc=13;  Character(ma) :: BFC    = 'merge.c'
      Integer :: muw=14;  Character(ma) :: BFW    = 'merge.w'

!----------------------------------------------------------------------

      Call inf_merge

      Open(pri,file=AF_log)

      nt=0
!----------------------------------------------------------------------
! ... read data from input file if any:

      Call Read_aarg('inp',AFi)
      if(Icheck_file(AFi).ne.0) then

       open(inp,file=AFi)
       Do 
        read(inp,'(a)',end=1) AF
        if(LEN_TRIM(AF).eq.0) Cycle
        if(INDEX(AF,'=').ne.0) Cycle
        if(AF(1:1).eq.'*') Exit
        nt=nt+1; files(nt)=AF
       End do
    1  Continue

       Call Read_apar(inp,'merge',merge)
       Call Read_rpar(inp,'eps',eps)
       Call Read_iarg('jmin',jmin)
       Call Read_iarg('jmax',jmax)
       Call Read_iarg('icheck',icheck)

      end if

!----------------------------------------------------------------------
! ... define input data from command line if any:


      Call Read_aarg('merge',merge)
      write(pri,'(a,a/)') 'merge =  ',merge

      Call Read_rarg('eps',eps)
      write(pri,'(a,d12.5/)') 'eps   = ',eps

      Call Read_iarg('jmin',jmin)
      write(pri,'(a,i3/)') 'jmin  = ',jmin

      Call Read_iarg('jmax',jmax)
      write(pri,'(a,i3/)') 'jmax  = ',jmax

      Call Read_iarg('icheck',icheck)
      write(pri,'(a,i3/)') 'icheck  = ',icheck

      iarg = command_argument_count()
      Do ia=1,iarg
       Call get_command_argument(ia,AF); if(INDEX(AF,'=').ne.0) Cycle
       nt=nt+1; files(nt)=AF
      End do

      if(nt.eq.0) then 
       write(pri,'(/a/)') 'There are no names in the command line '       
       Stop 'There are no names in the command line '       
      else
       write(pri,'(/a,i5/)') 'number of files =',nt       
       write(*,'(/a,i5/)') 'number of files =',nt       
      end if

      write(pri,'(/72(''-''))')

!----------------------------------------------------------------------
! ... files for final results:

      i=LEN_TRIM(merge); if(merge(i-1:i).eq.'.c') i=i-2 
      BFC = merge(1:i)//'.c'
      BFW = merge(1:i)//'.w'

      open(muc,file=BFC)
      open(muw,file=BFW,form='UNFORMATTED')

      Call alloc_orb_LS(-1)
      Call Alloc_radial(-1)

!----------------------------------------------------------------------
!                                                     proceed the files:
       Do it = 1,nt

        AF = files(it)

        i=LEN_TRIM(AF); if(AF(i-1:i).eq.'.c') i=i-2
        AFC=AF(1:i)//'.c'
        AFW=AF(1:i)//'.w'       

        write(pri,'(/i5,a,a)') it,'.  ',AF(1:i)
        Call Check_file(AFC)
        open(nuc,file=AFC)
        Call Check_file(AFW)
        open(nuw,file=AFW,form='UNFORMATTED')

        Call R_CLOSED(nuc)
        if(it.eq.1) then
         core = CLOSED
         write(muc,'(/a)') trim(core)
         read(nuw) Atom,Term,el3,MX(1),Z
         Call INITR; rewind(nuw)
        else
         if(CLOSED.ne.core) Stop 'check core !!!' 
        end if
        CALL  SUB1

       End do  ! over it
       write(muc,'(a)') '*'

!-----------------------------------------------------------------------
! ... create the final w file:

      Do i = 1,nrf 
       if(ks(i).gt.61) then
        el3 = 'xxx'
        write(atom,'(2i3)') n(i),l(i)
        write(term,'(i6)') ks(i)
       else
        el3 =  ELF3(N(i),L(i),ks(i))
       end if
       write(muw) Atom,Term,el3,MX(i),Z,EI,ZETA,AZ(i),P(1:MX(i),i)
      End do
      write(pri,'(72(''-''))')

! ... final list of target one-electron orbitals: 

      write(pri,'(/a,i5/)') 'nrf = ',nrf
      write(pri,'(14a4)') (EL(j),j=1,nrf)
      write(pri,'(72(''-''))')

!----------------------------------------------------------------------
! ... remove close configurations:

      if(icheck.ne.0) then

       ncfg=0;  lcfg=0; Call R_conf_LS(muc,0) 
       Do ic = 1,ncfg-1       
        Call Get_cfg_LS(ic)
        Call Save_cfg_LS(1)
       Do jc = ic+1,ncfg
        Call Get_cfg_LS(jc)
        if(iterm.ne.iterm1) Cycle 
        if(iconf.ne.iconf1) Cycle
        if(abs(WC(ic)).gt.abs(WC(jc))) WC(jc)=0.d0
        if(abs(WC(jc)).gt.abs(WC(ic))) WC(ic)=0.d0
       End do
       End do

       rewind(muc)
       write(muc,*) 
       write(muc,'(a)') trim(closed)
       Do ic=1,ncfg
        if(WC(ic).eq.0.d0) Cycle
        Call Pri_conf (muc,ic,WC(ic))
       End do

      end if

      CONTAINS
     

!======================================================================
      SUBROUTINE SUB1
!======================================================================
!
!     this routine analize one set of c- and w- files: AFC and AFW
!     and output information in files BFC and BFW
!
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)
      
!----------------------------------------------------------------------
! ... read configurations and define the list of orbitals:

      nwf=nclosd;  kset=0

      ncfg=0;  lcfg=0; Call R_conf_LS(nuc,0) 

      write(pri,'(/a,i5/)') ' nwf =',nwf
      write(pri,'(18a4)') ELF(1:nwf)
      write(pri,*)

      IEF = 0

! ... read radial functions:

      rewind(nuw)
    1 m=nrf+1; if(m.ge.mrf) CALL Alloc_radial(mrf+jrf)

      read(nuw,end=2) Atom,Term,el3,MX(m),ZZ,EI,ZETA,AZ(m), &
                      (P(j,m),j=1,MX(m))

      EL(m) = el3

      if(Z.ne.ZZ) then; write(pri,*) 'Z,ZZ =', Z,ZZ; Stop ' Check Z !'
      end if

      P(MX(m)+1:NO,m) = 0.d0

      Call EL4_nlk(EL(m),N(m),L(m),KS(m))
      EL4 = EL(m); nnn=N(m); lll=L(m); kkk=KS(m)
      mexp(m) = L(m) + 1
      aexp(m) = - D1/mexp(m)
      ZR = Z*R(1)
      bexp(m) = P(1,m)/(AZ(m)*R2(1)*R(1)**L(m))
      bexp(m) = (bexp(m) - D1 + ZR/(L(m)+1) )/ZR**2
 
      ii = Ifind_nlk(N(m),L(m),KS(m),0)
      if(ii.eq.0) then
       write(pri,'(a4,a)') EL(m),'  the excessive orbital'
       go to 1   ! go to next orbital
      end if

! ...  define OVERLAPSs with existing orbitals:

       Do i = 1,m                  
        OVERLAPS(i,m)=0.d0; if(L(i).eq.L(m)) OVERLAPS(i,m)=QUADR(i,m,0)
       End do

! ...  store the core orbitals in case of first w-file:

       if(m.le.nclosd) then; nrf=m; IEF(ii)=nrf; go to 1; end if

!---------------------------------------------------------------------
!                                  compare with the existing orbitals:
       
       Do i = 1,nrf;  if(abs(OVERLAPS(i,m)).lt.eps) Cycle

        ! ... check orthogonality to core: 

        if(i.le.nclosd.and.ii.gt.nclosd) then
         write(pri,'(a,a,a,f10.5)') '  orbital ', EL4, &
              ' does not orthogonal to core:', OVERLAPS(i,m)
         Stop 'problem with orthogonality to core'
        end if

        ! ... define is that approximately the same orbital:

        S  = abs(OVERLAPS(i,m)-1.d0)
        S1 = abs(QUADR(i,i, 1)-QUADR(m,m, 1))
        S2 = abs(QUADR(i,i, 2)-QUADR(m,m, 2))

        if(S.lt.eps.and.S1.lt.eps.and.S2.lt.eps) then           
         IEF(ii) = i
         write(pri,'(a,a,a,a,3f15.9)')  &
               EL4,' --> ',EL(i),'    same                  ',S
         go to 1
        end if
       
       End do

       ! ...  core orbitals should be the same:   

       if(it.gt.1.and.ii.le.nclosd) then
        write(pri,'(a,a,a)') 'file ',AFW,'  has another core orbital'
        Stop ' anoter core orbital? '
       end if

!---------------------------------------------------------------------
!                                    assign set index for new orbital: 
       ks(m)=-1
       knew = New_index(L(m),ksmax,nrf,L,KS)
       Do k = 1,knew-1
        S=0.d0           
        Do i = 1,nrf
         if(L(m).ne.L(i).or.k.ne.ks(i)) Cycle
         S=max(S,abs(OVERLAPS(i,m)))
        End do
        if(S.lt.eps) then; ks(m)=k; Exit; end if
       End do  

       if(ks(m).eq.-1) then  ! the orbital belongs to new set  

        ks(m) = knew


        EL(m)=ELF4(N(m),L(m),ks(m))
        write(pri,'(a,a,a,a,f15.9)') &
              EL4,' --> ',EL(m),'    new orbitals and new k',S

       else                   ! check the same label for diff.orbitals            


        EL(m)=ELF4(N(m),L(m),ks(m))
        Do i = 1,nrf
         if(EL(m).ne.EL(i)) Cycle
         write(pri,'(a)') 'the same n for orthogonal orbitals ? '
         write(pri,'(a5,a5,f10.5)')  EL(m),EL(i),OVERLAPS(i,m)
         Stop ' the same n for orthogonal orbitals ?'
        End do
        write(pri,'(a,a,a,a,f15.9)') &
              EL4,' --> ',EL(m),'    new orbitals but old k',S
       
       end if

       nrf=m; IEF(ii)=m

       go to 1    ! go to next orbital
   2  Continue

      Do i=nclosd+1,nwf
       if(IEF(i).eq.0) then
        write(pri,'(a,i3,a,a,a)') ' target ',it,'  orbital ',ELF(i), &
                                   ' not found in the w-file'
        Stop ' unknown orbitals ! '
       end if
      End do

!----------------------------------------------------------------------
!                                                     write new c-file:
      rewind(nuc)
      read(nuc,'(a)') AS
      read(nuc,'(a)') AS

    3 read(nuc,'(a)',end=4) AS; if(AS(5:5).ne.'(') go to 3
      read(AS,'(a64,f11.8)') CONFIG,C;  read(nuc,'(a60)') COUPLE
      
      Call Decode_c

      Do i=1,no
       ii=Ifind_nlk(nn(i),ln(i),kn(i),0); ii=IEF(ii)
       nn(i)=N(ii); ln(i)=L(ii); kn(i)=ks(ii)
      End do
      
      Call Incode_c

      ILT=LS(no,4)
      IST=LS(no,5)
      if(jmin.ne.-1.and.jmin.gt.ILT+IST-2) go to 3 
      if(jmax.ne.-1.and.jmax.lt.iabs(ILT-IST)) go to 3 

      write(muc,'(a64,f11.8)') CONFIG,C
      write(muc,'(a)') trim(COUPLE)
       
      go to 3
    4 Continue

! ... additional check

      Do i=1,nrf-1
       Do j=i+1,nrf

        if(L(i).ne.L(j)) Cycle
        S = abs(OVERLAPS(i,j))

        if(ks(i).eq.ks(j).and.S.gt.eps)  then
         write(pri,*) & 
              ' k_set1 = k_set2, but orbitals are not orthogonal ???'
         write(pri,'(2a5,f12.8)') EL(i),EL(j),S
         Stop ' k_set1 = k_set2, but orbitals are not orthogonal ???'
        end if

       End do
      End do

      write(pri,'(72(''-''))')
      
      End Subroutine Sub1
      

      End  ! program NO_PREP


!======================================================================
      Subroutine inf_merge
!======================================================================
!     provide screen information about MERGE utility
!----------------------------------------------------------------------
       
      Implicit real(8) (a-h,o-z)

      Character(80) :: A

      iarg = command_argument_count()
      if(iarg.gt.0) then
       Call get_command_argument(1,A);  if(A.ne.'?')  Return
      end if

      write(*,'(a)') &
'                                                                  ',&
'Merging c(w)-files in one file with consistent set indexes for    ',&
'one-electron radial functions                                     ',&
'                                                                  ',&
'In order to define the same orbitals, the following criteria are  ',&
'used:                                                             ',&
'            | <p1|p2> - 1 |  <  eps                               ',&
'            | <p1|r|p1> - <p2|r|p2> | <  eps                      ',&
'            | <p1|1/r|p1> - <p2|1/r|p2> | <  eps                  ',&
'                                                                  ',&
'INPUT DATA are defined as arguments in command line:              ',&
'                                                                  ',&
'   merge  name1 name2 ... merge=... eps=...  jmin=... jmax=...    ',&
'                                                                  ',&
'Default values:  merge=merge; eps=1.d-7;  jmin=jmax=-1 (N/A)      ',&
'                                                                  ',&
'Result: merge.c, merge.w (see merge.log for running information)  ',&
'                                                                  ',&
'jmin,jmax restrict output configurations for 2J-values            ',&
'                                                                  '
      Stop ' '

      End Subroutine inf_merge

