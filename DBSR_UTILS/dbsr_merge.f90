!======================================================================
!     PROGRAM       D B S R _ M E R G E              
!
!               C O P Y R I G H T -- 2008
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!======================================================================
!  Analizes the orthogonality conditions for one-electron orbitals 
!  and generates new c- and bsw-files  with consistent
!  indication of the set indexes 
! 
!  In order to define the same orbitals, the following criteria are 
!  used:  
!              | <p1|p2> - 1 |  <  eps_ovl
!              | <p1|r|p1> - <p2|r|p2> | <  eps_ovl
!              | <p1|1/r|p1> - <p2|1/r|p2> | <  eps_ovl
!
!  where the parameter eps_ovl can be read from bsr_par, 
!  or as argument: eps_ovl=...
!======================================================================
!
!  INPUT FILES:
!
!     c- and bsw  -  configuration and w.f. files for each case
!
!  OUTPUT FILES:
!
!     merge.c      -  modified list of target states
!     merge.bsw    -  all target bsw-functions
!
!======================================================================
      Use conf_jj; Use orb_jj
      Use DBS_grid; Use DBS_orbitals_pq; Use DBS_gauss
      
      Implicit real(8) (A-H,O-Z)

      Integer, parameter :: ma = 80
      Character(ma) :: AFC,AFW, BFC,BFW, line
      Character(ma) :: AF_inp='no', AF_out='merge', AF_log='dbsr_merge.log'
      Character(ma), allocatable :: files(:)
      Integer :: pri = 66;         ! log file
      Integer :: inp = 1;          ! list of c-files
      Integer :: nuc = 11;         ! input c-file
      Integer :: nuw = 12;         ! input w-file
      Integer :: muc = 13;         ! temporaly c-file

      Real(8), allocatable :: pbw(:),qbw(:),tw(:)

! ... default value for overlap parameter: 

      Real(8) :: eps_ovl  = 1.d-7
      Real(8) :: eps_core = 1.d-5
      Integer :: nfile = 0

      Call Read_name(line)
      if(line.eq.'?') then
       write(*,*)
       write(*,*) 'Merging a set of (name.c, name.bsw) files to one pair (merge.c, merge.bsw)'
       write(*,*) 'with consistent set-indexes if orbitals are not-orthogonal'
       write(*,*)
       write(*,*) 'Call as:  dbsr_merge AF1.c AF2.c ... nfile=.. merge=..'
       write(*,*)
       write(*,*) 'list of merging c-files are defined as first position arguments'
       write(*,*)
       write(*,*) 'list of merging files (AF1.c, AF.2.c, ...) '
       write(*,*) 'also can be given in input file defined through inp=...  '
       write(*,*) 'one filename per line'  
       write(*,*) 'In this case call as: dbsr_merge inp=...  merge=...'
       write(*,*)  
       write(*,*) 'nfile  - number of merging files'
       write(*,*) 'merge  - name for resulting c-file: merge.c and merge.bsw'
       write(*,*) 'jjmin  - minimum 2J value (optional)'
       write(*,*) 'jjmax  - maximum 2J value (optional)'
       write(*,*) 'eps_ovl- tolerance for one-electron overlaps (optional)'
       write(*,*) 'eps_core - tolerance for overlaps with core functions (optional)'
       Stop ' '
      end if

!----------------------------------------------------------------------
      Open(pri,file=AF_log)

      Call Read_aarg('inp',AF_inp)
      if(AF_inp.ne.'no') then

       Call Check_file(AF_inp)
       open(inp,file=AF_inp)
       nfile=0       
    11 read(inp,'(a)',end=12) line
       if(Index(line,'=').ne.0) go to 11
       if(LEN_TRIM(line).eq.0) go to 11
       nfile=nfile+1
       go to 11
    12 rewind(inp)
       if(nfile.eq.0) Stop 'nfile = 0 --> nothing to do'
       write(pri,'(/a,i6/)') 'nfile    =',nfile
       Allocate(files(nfile))
       i = 0
    13 read(inp,'(a)',end=14) line
       if(Index(line,'=').ne.0) go to 13
       if(LEN_TRIM(line).eq.0) go to 13
       i = i + 1; files(i)=line
       go to 13
    14 rewind(inp)

       jjmin=-1;    Call Read_ipar(inp,'jjmin', jjmin )
       jjmax=jjmin; Call Read_ipar(inp,'jjmax', jjmax )

       Call Read_rpar(inp,'eps_ovl' ,eps_ovl )
       Call Read_rpar(inp,'eps_core',eps_core)

       Call Read_apar(inp,'out',AF_out)

      else

       iarg = command_argument_count()
       nfile = 0
       Do i = 1,iarg
        Call GET_COMMAND_ARGUMENT(i,AFC)
        if(Index(AFC,'=').ne.0) Cycle
        if(LEN_TRIM(AFC).eq.0) Cycle
        nfile = nfile + 1
       End do
       if(nfile.eq.0) Stop 'nfile = 0 --> nothing to do'
       write(pri,'(/a,i6/)') 'nfile    =',nfile
       Allocate(files(nfile))
       j = 0
       Do i = 1,iarg
        Call GET_COMMAND_ARGUMENT(i,AFC)
        if(Index(AFC,'=').ne.0) Cycle
        if(LEN_TRIM(AFC).eq.0) Cycle
        j = j + 1; files(j)=AFC
       End do
      end if

      Call Read_rarg('eps_ovl', eps_ovl )
      Call Read_rarg('eps_core',eps_core)

      write(pri,'(a,1Pd12.2)') 'eps_ovl  = ',eps_ovl
      write(pri,'(a,1Pd12.2)') 'eps_core = ',eps_core

      jjmin=-1;    Call Read_iarg('jjmin', jjmin )
      jjmax=jjmin; Call Read_iarg('jjmax', jjmax )

      write(pri,'(/a,i5)') 'jjmin    = ',jjmin
      write(pri,'( a,i5)') 'jjmax    = ',jjmax

      Call Read_aarg('out',AF_out)

      write(pri,'(/a,a)') 'AF_out   = ',AF_out

! ... check style of names: (just name or name.c)

      Do i=1,nfile
       Call Check_BSR_name(files(i))
      End do
      Call Check_BSR_name(AF_out)

! ... initial allocations for one-electron orbitals

      nwf=0; Call alloc_orb_jj(iwf)

!----------------------------------------------------------------------
      Do it=1,nfile
 
       AFC = trim(files(it))//'.c'

       Call Check_file(AFC)
       Open(nuc,file=AFC)
              
       if(it.eq.1) then
        Call Read_core_jj(nuc)
       else
        read(nuc,'(/a)') CLOSED
        if(CLOSED.ne.core)  then
         write(pri,'(/a,a,a)') 'case ',AFC,' has different core'
         Stop ' ' 
        end if
       end if

       write(pri,'(/a,a/)') 'case ',AFC
      
       AFW = trim(files(it))//'.bsw'

       CALL  SUB1

      End do 


!----------------------------------------------------------------------
!                                                     write new c-file:

      AFC = trim(AF_out)//'.c'

      open(nuc,file=AFC)
      write(nuc,'(a)') 'Core subshells:'
      write(nuc,'(a)')  core(1:ncore*5)
      write(nuc,'(a)') 'Peel subshells:'
      write(nuc,'(20a5)') (ebs(i),i=ncore+1,nbf)
      write(nuc,'(a)') 'CSF(s):' 

      ncfg=0; lcfg=0; nwf=ncore; Call read_conf_jj(muc,0,'no','no') 

! ... re-define jj_min and jj_max:

      jj_min = 1000;  jj_max = 0
      Do ic=1,ncfg
       jj = jot_ic(ic)
       if(jj.lt.jj_min) jj_min=jj
       if(jj.gt.jj_max) jj_max=jj
      End do

      Do jj = jj_min,jj_max,2
       Do ic=1,ncfg
        if(jot_ic(ic).ne.jj) Cycle
        Call print_conf_jj (nuc,ic,wc(ic))
       End do
       if(jj.lt.jj_max) write(nuc,'(a)') ' *'
       if(jj.eq.jj_max) write(nuc,'(a)') '*'
      End do

      Close(nuc)
      Close(muc,status='DELETE')

!-----------------------------------------------------------------------
! ... create the final.bsw file:

      BFW = trim(AF_out)//'.bsw'

      Open(nuw,file=BFW,form='UNFORMATTED')
      write(nuw) grid_type,ns,ks,t(1:ns+ks),ksp,ksq
      Do i = 1,nbf
       write(nuw) ebs(i),mbs(i)
       write(nuw) pq(1:mbs(i),1,i)
       write(nuw) pq(1:mbs(i),2,i)
      End do
      Close(nuw)

CONTAINS

!======================================================================
      SUBROUTINE SUB1
!======================================================================
!     this routine analize one set of c- and w- files: AFC and AFW
!     and output information in files BFC and BFW
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Character(5) :: elw
      Character(5), External :: ELi

!----------------------------------------------------------------------
! ... read configurations and define the list of orbitals:

      ncfg=0; nwf=ncore; Call read_conf_jj(nuc,0,'add','check') 
      ipef=0

! ... read radial functions:

      Call Check_file(AFW);  Open(nuw,file=AFW,form='UNFORMATTED')

      if(it.eq.1) then
        Call Read_grid_bsw(nuw)
        Call alloc_DBS_gauss
        Allocate(pbw(ns),qbw(ns),tw(ns+ks))
        nbf=0; CALL alloc_DBS_orbitals_pq(ibf,ns)
      else
        read(nuw) igrid,nsw,ksw,tw(1:nsw+ksw),kp,kq
        if(igrid.ne.grid_type) Stop 'Another knot grid?'      
        if(ksw.ne.ks) Stop ' Read_bsw:  ksw <> ks'
        if(nsw.ne.ns) Stop ' Read_bsw:  nsw <> ns'
        if(ksp.ne.kp) Stop ' Read_bsw:  ksp <> kp'
        if(ksq.ne.kq) Stop ' Read_bsw:  ksq <> kq'
        k=1
        Do i=1,ns+ks
         if(abs(t(i)-tw(i)).lt.1.d-10) Cycle; k=0; Exit
        End do    
        if(k.eq.0) Stop 'Another knot grid t(1:ns+ks) ?'
      end if

    1 m=nbf+1; if(m.gt.mbf) CALL Alloc_dbs_orbitals_pq(mbf+ibf,ns)
      Call Idel_obs(m)

      read(nuw,end=2) elw,mbw

      pbw=0.d0; read(nuw) (pbw(i),i=1,mbw)
      qbw=0.d0; read(nuw) (qbw(i),i=1,mbw)
     
      Call EL_nljk(elw,n,kappa,l,j,k); ii=Ifind_jjorb(n,kappa,k,0)

      if(ii.gt.0) then                              
       mbs(m)=ns; nbs(m)=n; kbs(m)=kappa; ibs(m)=k; ebs(m)=elw
       pq(1:ns,1,m)=pbw(1:ns)
       pq(1:ns,2,m)=qbw(1:ns)
      else
       write(pri,'(a5,a)') elw,'  the excessive orbital'
       go to 1
      end if

! ... define overlaps with existing orbitals:

      Do i = 1,m
       if(kbs(i).ne.kbs(m)) Cycle 
       S = QUADR_pq(m,i,0)
       if(abs(S).lt.eps_ovl) Cycle
       Call Iadd_obs(m,i,S)
      End do
      SM1 = QUADR_pq(m,m,1)
      SM2 = QUADR_pq(m,m,2)

! ... store the core orbitals in case of first w-file:

      if(m.le.ncore) then; nbf=m; ipef(ii)=nbf; go to 1; end if

!---------------------------------------------------------------------
!                                  compare with the existing orbitals:
       
      Do i = 1,nbf;  if(abs(OBS(i,m)).lt.eps_ovl) Cycle

       ! ... check orthogonality to core: 

       if(i.le.ncore.and.ii.gt.ncore) then
        if(abs(OBS(i,m)).gt.eps_core) then
        write(pri,'(a,a,a,f10.5)') &
        '  orbital ', elw, ' does not orthogonal to core:', OBS(i,m)
        Stop 'problem with orthogonality to core'
        end if
       end if

       ! ... define is that approximately the same orbital:

       S  = abs(OBS(i,m)-OBS(m,m)) +  abs(OBS(i,m)-OBS(i,i))
       S1 = abs(QUADR_pq(i,i,1)-SM1)
       S2 = abs(QUADR_pq(i,i,2)-SM2)


       if(S.lt.eps_ovl.and.S1.lt.eps_ovl.and.S2.lt.eps_ovl) then           
        ipef(ii) = i
        write(pri,'(a,a,a,a)') elw,' --> ',ebs(i),'    the same'
        go to 1
       end if
       
      End do

      ! ...  core orbitals should be the same:   

      if(it.gt.1.and.ii.le.ncore) then
       write(pri,'(a,a,a)') 'file ',AFW,'  has another core orbital'
       Stop ' anoter core orbital? '
      end if

!---------------------------------------------------------------------
!                                    assign set index for new orbital: 
      ibs(m)=-1
      inew = New_index(kappa,ksmax,nbf,kbs,ibs)
      S = 0.d0

! ... check existing orthogonal subsets:

      Do i = 1,inew-1
       S=0.d0           
       Do j = 1,nbf
        if(kappa.ne.kbs(j).or.i.ne.ibs(j)) Cycle
        S=max(S,abs(OBS(j,m)))
       End do
       if(S.lt.eps_ovl) then; ibs(m)=i; Exit; end if
      End do  

      if(ibs(m).eq.-1) then  ! the orbital belongs to new set  

       ibs(m) = inew
       EBS(m)=ELi(nbs(m),kbs(m),ibs(m))
       write(pri,'(a,a,a,a,f15.9)') &
             elw,' --> ',EBS(m),'    new orbitals and new set index',S

      else              ! check the same label for diff.orbitals            

       EBS(m)=ELi(nbs(m),kbs(m),ibs(m))
       Do i = 1,nbf
        if(EBS(m).ne.EBS(i)) Cycle
         write(pri,'(a)') 'the same n for orthogonal orbitals ? '
         write(pri,'(a5,a5,f10.5)')  EBS(m),EBS(i),OBS(i,m)
         Stop ' the same n for orthogonal orbitals ?'
       End do
       write(pri,'(a,a,a,a,f15.9)') &
             elw,' --> ',EBS(m),'    new orbitals but old set index',S
       end if

       nbf=m; ipef(ii)=m

       go to 1    ! go to next orbital
   2  Continue

! ... check if we have all radial functions:

      Do i=ncore+1,nwf
       if(ipef(i).eq.0) then
        write(pri,'(a,a,a)') 'orbital ',ELF(i), &
                            ' not found in the w-file'
        Stop ' unknown orbitals ! '
       end if
      End do

!-----------------------------------------------------------------------

      rewind(nuc)
    3 read(nuc,'(a)',end=4) AS; if(AS(6:6).ne.'(') go to 3
      i = INDEX(AS,')',BACK=.TRUE.)
      if(i.le.72) then
       read(AS,'(a72,f12.8)') CONFIG,C
      else
       CONFIG = AS(1:i); read(AS(i+1:i+12),'(f12.8)') C
      end if
      read(nuc,'(a)') SHELLJ
      read(nuc,'(9x,a)') INTRAJ
      
      Call Decode_cj

      if(jjmin.ge.0.and.Jtotal.lt.jjmin) go to 3 
      if(jjmax.ge.0.and.Jtotal.gt.jjmax) go to 3 

      Do i=1,no
       ii=Ifind_jjorb(nn(i),kn(i),in(i),1); ii=ipef(ii)
       nn(i)=nbs(ii); kn(i)=kbs(ii); in(i)=ibs(ii); 
      End do
      
      Call Incode_cj

      ii = 72; if(no.gt.8) ii=9*no
      write(muc,'(a,f12.8)') CONFIG(1:ii),C
      write(muc,'(a)') SHELLJ
      write(muc,'(9x,a)') INTRAJ
       
      go to 3
    4 Continue

      write(pri,'(72(''-''))')
      
      End Subroutine Sub1
      
      End  ! program dbsr_merge



!======================================================================
      Subroutine Check_BSR_name(AF)
!======================================================================
!     rename old-version file-name as name.c just to 'name'
!----------------------------------------------------------------------
      Character(*) :: AF
      i = LEN_TRIM(AF)
      if(i.le.2) Return
      if(AF(i-1:i).eq.'.c') AF(i-1:i) = '  '
      End Subroutine Check_BSR_name

