!======================================================================
!     utility       D E L _ C O R E
!
!                   C O P Y R I G H T -- 2004
!
!     Written by:   Oleg Zatsarinny
!======================================================================
!
!     This utility moves the orbital from core to configurations
!
!     Call as   del_core  el=...  inp=... out=...
!
!----------------------------------------------------------------------
      Use conf_ls

      Implicit real(8) (A-H,O-Z)
 
      Character(80) :: AF, BF, AS
      Character(4) :: EL
      Character(4), External :: ELF4

      Integer :: nu1=1; Character(40) :: AF_inp      ! initial c-file
      Integer :: nu2=2; Character(40) :: AF_out      ! resulting c-file

      iarg = command_argument_count() 
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)
 
      if(AF.eq.'?'.or.AF.eq.'!') then
       write(*,'(/a)') 'This utility moves the orbital from core to configurations'
       write(*,'(/a)') 'Call as:  del_core  el=...  inp=... out=...'
       write(*,'(/a)') 'el  - orbital to move'
       write(*,'( a)') 'inp - input c-file'
       write(*,'(a/)') 'out - output c-file'
       Stop
      end if
!----------------------------------------------------------------------
!                                                           input data:

      EL = ' ';      Call Read_aarg('el',el)
      AF_inp = ' ';  Call Read_aarg('inp',AF_inp)   
      AF_out = ' ';  Call Read_aarg('out',AF_out)   

      if(len_trim(EL).eq.0.or.len_trim(AF_inp).eq.0.or.len_trim(AF_out).eq.0) then
        write(*,*) 'You should provide core orbital to add and input and output files as: ' 
        write(*,*) 'del_core el=... inp=... out=...'
        Stop ' '
      end if

      Call EL4_nlk(EL,n1,l1,k1); EL=ELF4(n1,l1,k1)

      Call Check_file(AF_inp)       

      Open(nu1,file=AF_inp,action='read')
      Open(nu2,file=AF_out,action='write')

!----------------------------------------------------------------------
!                                                            read core:
      read(nu1,'(a)') AS
      write(nu2,'(a)') AS

      read(nu1,'(a)') CLOSED
      i=LEN_TRIM(CLOSED); NCLOSD=i/4;  i1=0
      Do i=1,NCLOSD
       j=(i-1)*4; if(CLOSED(j+1:j+4).eq.EL) i1=i
       if(i1.gt.0) Exit
      End do
      if(i1.eq.0) then
       write(*,*) 'EL=', EL
       write(*,*) trim(CLOSED)
       write(*,*) 'failed to find specified orbital'
      end if

      Do i=i1,NCLOSD
       j=(i-1)*4; CLOSED(j+1:j+4)=CLOSED(j+5:j+8)
      End do
      write(nu2,'(a)') CLOSED

!----------------------------------------------------------------------
!                                                  read configurations:
   10 read(nu1,'(a)',end=20) AS
      if(AS(1:1).eq.'*') go to 20
      if(AS(5:5).ne.'(') then
       write(nu2,'(a)') AS; go to 10
      end if
      read(AS,'(a64,F11.8)') CONFIG, C
      read(nu1,'(a60)') COUPLE

      Call Decode_c

!     add deleted shell from core

      Do i=no,1,-1
       i1=i+1
       nn(i1)=nn(i); ln(i1)=ln(i); iq(i1)=iq(i); kn(i1)=kn(i)
       LS(i1,1:5)=LS(i,1:5)
      End do

      no=no+1
      if(no.gt.msh) Stop ' no > nsh'
      nn(1)=n1; ln(1)=l1; iq(1)=4*l1+2; kn(1)=k1
      LS(1,1)=0; LS(1,2)=1; LS(1,3)=1; LS(1,4)=1; LS(1,5)=1
      Call Incode_c

      write(nu2,'(a64,F11.8)') CONFIG, C
      write(nu2,'(a)') COUPLE

      go to 10
   20 Close(nu1)

      write(nu2,'(a)') '*'      
      Close(nu2)

      End ! program DEL_CORE
