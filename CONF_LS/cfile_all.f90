!======================================================================
!     utility       C F I L E _ A L L
!
!                   C O P Y R I G H T -- 2015
!
!     Written by:   Oleg Zatsarinny
!
!======================================================================
!
!     Extract the states from  l(j)-file in separate c-files
!
!     arguments: 1. name for l- or j-file
!                2. msol = number of solutions (optional)
!                3. jj = 2J value
!                3. eps = cut_off factor (optional)
!
!     When eps > 0, configurations in result c-file are ordered
!     according to their weights
!
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)

      Logical :: EX
      Character(80) :: AS, clousd
      Character(40) :: name, AN, AF ,BF

      Real(8), Allocatable :: WT(:)
      Integer, Allocatable :: IP(:)
      Character(64), Allocatable :: AC(:)
      Character(60), Allocatable :: AU(:)

      Integer :: nuc =1       !   name.c         
      Integer :: nuj =2       !   name.l or name.j
      Integer :: iout=3       !   result.c

      Real(8) :: eps_c = 1.d-7
      Integer :: jj = -1, msol = 0

!----------------------------------------------------------------------
!                                                           input data:
      iarg = command_argument_count() 
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)

      if(AF.eq.'?'.or.iarg.lt.1) then
        write(*,'(/a)') 'cfile_all extracts the solutions from l(j)-file in separate c-file' 
        write(*,'(/a)') 'Call as:    cfile name.l(j) msol=... eps= '
        write(*,'(/a)') 'name.l(j) - input file with solutions'
        write(*,'( a)') 'msol      - number of solutions needed '
        write(*,'( a)') 'jj        - 2J value '
        write(*,'( a)') 'eps       - tolerance for weights, default: 1e-7'
        write(*,'(/a)') 'When eps > 0, configurations in result c-file are ordered'
        write(*,'( a)') 'according to their weights' 
        Stop ' '
      else
        Call Read_name(AN)
        Call Read_iarg('msol',msol)
        Call Read_iarg('jj',jj)
        Call Read_rarg('eps',eps_c) 
      end if

!----------------------------------------------------------------------
! ... read list of configuration:

      i = LEN_TRIM(AN)
      name = AN(1:i-2)
      AF=AN;  AF(i:i)='c'
      Open(nuc,file=AF,STATUS='OLD')

      ncfg=Idef_ncfg(nuc)
      if(allocated(WT)) Deallocate(WT,AC,AU,IP)
      Allocate(WT(ncfg), AC(ncfg), AU(ncfg), IP(ncfg))

      rewind(nuc)
      read(nuc,'(a)') AS
      read(nuc,'(a)') Clousd
      i=0
    2 read(nuc,'(a)',end=3) AS
      if(AS(1:1).eq.'*') go to 3
      if(AS(5:5).ne.'(') go to 2
      i=i+1; AC(i)=AS(1:64)
      read(nuc,'(a)') AU(i)
      go to 2
    3 Close(nuc)
!----------------------------------------------------------------------
! ... read solution:

      Open(nuj,file=AN,STATUS='OLD')
      rewind(nuj)
    5 read(nuj,'(a)',end=9) AS
      if(AS(3:5).ne.'2*J') go to 5
      read(AS,'(7x,i5,10x,i4)') j,nj
      if(jj.ge.0.and.j.ne.jj) go to 5

      if(msol.gt.0.and.msol.lt.nj) nj = msol

      Do is=1,nj
       read(nuj,'(/i6,f16.8,3x,a)') n,E,AF
       read(nuj,'(7F11.8)') (WT(i),i=1,ncfg)
       Call SORTA(ncfg,WT,IP)  

! ... output the c-file:

      write(AF,'(a,a,i3.3,a,i3.3)') trim(name),'_',jj,'_',is

      BF = trim(AF)//'.c'
      open(iout,file=BF)
      if(jj.ge.0) then 
       write(iout,'(15x,f16.8,a,i3)') E,'   2*J = ',jj
      else
       write(iout,'(15x,f16.8)') E
      end if
      write(iout,'(a)') Clousd
 
      Do m=1,ncfg; i=m; if(eps_c.gt.0.d0) i=IP(m)
       if(abs(WT(i)).lt.Eps_c) Cycle
       write(iout,'(a64,f11.8)') AC(i),WT(i)
       AS=AU(i); ii=len_trim(AS)
       write(iout,'(a)') AS(1:ii)
      End do
      write(iout,'(a)') '*'
      Close(iout)

      i = LEN_TRIM(AN); AN(i:i)='w'
      if(Icheck_file(AN).eq.1) then
       BF = trim(AF)//'.w'
       write(AS,'(a,a,a,a)') 'cp ',trim(AN),' ',trim(BF)
       i = SYSTEM(AS)
      end if

      End do  !  over solutions
      Stop
    9 write(*,*) ' No solution for jj =', jj  
     
      End  ! program cfile_all
