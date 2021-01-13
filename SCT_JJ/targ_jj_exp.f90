!=====================================================================
!     UTILITY  convert_target_jj
!=====================================================================
!     target_jj, threshoulds -> target_jj.exp
!
!     Call as:  indicated are the optional (default) paremeters:
!
!     convert_target  [ targ1=target_jj.theory
!                       targ2=target_jj.exp
!                       thresholds=thereholds ]
!---------------------------------------------------------------------

      Use target_jj; Use channels_jj

      Implicit real(8) (A-H,O-Z)

      Real(8), allocatable :: E_exp(:)
      Integer, allocatable :: ip_exp(:), jp_exp(:), npch(:), nptar(:)
      Character(80) :: AF

      Integer :: nu1=1;  Character(80) :: AF1  = 'target_jj.theory'
      Integer :: nu2=2;  Character(80) :: AF2  = 'target_jj.exp'
      Integer :: nu3=3;  Character(80) :: AF3  = 'thresholds'

!----------------------------------------------------------------------
      iarg = COMMAND_ARGUMENT_COUNT()
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)

      if(AF.eq.'?') then
        write(*,'(/a)') 'change target file to fit new experimental state energies'
        write(*,'(/a)') 'given in file thresholds'
        write(*,'(/a)') 'default:   target_jj.theory + thrsholds => target_jj.exp'
        write(*,'(/a)') 'change default names as:' 
        write(*,'(/a)') 'convert_target_jj  [targ1=...  targ2=... thresholds=...]'  
        Stop 
      end if        

      Call Read_aarg('targ1',AF1)
      Call Read_aarg('targ2',AF2)
      Call Read_aarg('thresholds',AF3)

!----------------------------------------------------------------------
! ... original target and channel information:

      Open(nu1,file=AF1,status='OLD')
      Call Read_target_jj (nu1)
      Call Read_channels_jj(nu1)
      rewind(nu1); read(nu1,'(a)') title
      Close(nu1)

      Z = nz
      Call Conv_au (Z,0.d0,au_cm,au_eV,0)

!----------------------------------------------------------------------
! ... new thresholds:

      Open(nu3,file=AF3,status='OLD')

      Allocate(E_exp(ntarg))
      Do i=1,ntarg; read(nu3,*) E_exp(i); End do
      Allocate(ip_exp(ntarg),jp_exp(ntarg))
      Call SORTR(ntarg,E_exp,ip_exp)
      iiexp=0
      Do i=1,ntarg; if(ip_exp(i).ne.i) iiexp=1; End do
      etarg = E_exp

      if(iiexp.eq.0) write(*,*) 'target states order is not changed'
      if(iiexp.eq.1) write(*,*) 'target states order is changed !!!'

!----------------------------------------------------------------------
! ... new target states:

      Open(nu2,file=AF2);  nut=nu2

      write(nut,'(a)') TRIM(title) 
      write(nut,'(80(''-''))')
      write(nut,'(a,i4,5x,a)') &
                'nz    = ',nz,   ' !   nuclear charge' 
      write(nut,'(a,i4,5x,a)') &
                'nelc  = ',nelc, ' !   number of electrons'
      write(nut,'(80(''-''))')
      write(nut,'(a,i4,5x,a)') &
                'ntarg = ',ntarg,' !   number of target states'
      write(nut,'(80(''-''))')

      ia = 0; ib = 0
      Do i=1,ntarg
       ia = max(ia,LEN_TRIM(AFT(i)))
       ib = max(ib,LEN_TRIM(BFT(i)))
      End do
      ia = max(ia,30)
      ib = max(ib,10)

      Do j=1,ntarg; i=ip_exp(j); jp_exp(i) = j
       E_Ry = (etarg(i)-etarg(1))*2
       E_eV = (etarg(i)-etarg(1))*au_eV
       write(nut,'(a,2x,a,2x,2i4,f18.8,2i5,f12.6,f10.3)') AFT(i)(1:ia),BFT(i)(1:ib), &
        jtarg(i),ptarg(i),etarg(i),nctarg(i),nwtarg(i), E_Ry,E_ev
      End do

      write(nut,'(80(''-''))')
      write(nut,'(a,i5,T18,a)') 'nct    =',nct, &
       '!   total number of target configurations' 
      write(nut,'(a,i5,T18,a)') 'nwt    =',nwt, &
       '!   total number of target orbitals' 
      write(nut,'(80(''-''))')

!----------------------------------------------------------------------------
! ... partial waves:

      write(nut,'(a,i4,5x,a)') 'nlsp  = ',nlsp,' !   number of partial waves' 
      write(nut,'(80(''-''))')

      Do i = 1,nlsp
       if(ncp(i).eq.0) then
        write(nut,'(i3,a,2i4,6x,a40,2x,a10,2i5)') &
         i,'.',jpar(i),ipar(i)
       else
        write(nut,'(i3,a,2i4,6x,a40,2x,a10,2i5)') &
         i,'.',jpar(i),ipar(i),AFP(i),BFP(i),ncp(i),nwp(i)
       end if
      End do
      write(nut,'(80(''-''))')
 
!---------------------------------------------------------------------------------
! ... new channels:

      write(nut,'(a,i5)') 'channels:' 
      write(nut,'(80(''-''))')

      Do ilsp = 1,nlsp
       write(nut,'(i3,a,a,i6)')  ilsp,'.',' nch =',nch(ilsp)
       write(nut,'(80(''-''))')

       ! ... new channel order:

       allocate( npch(nch(ilsp)), nptar(nch(ilsp)) )
       k = 0
       Do i=1,ntarg; it=ip_exp(i)      
        Do ich=1,nch(ilsp); if(iptar(ilsp,ich).ne.it) Cycle
         k=k+1; npch(k) = ich; nptar(k) = i
        End do
       End do
       if(k.ne.nch(ilsp)) Stop 'Problems with new channel order'

       Do j = 1,nch(ilsp); i = npch(j)
        write(nut,'(i4,2x,a,2x,2i6,i8)') &
         i,ELC(ilsp,i),kch(ilsp,i),nptar(j),ipconf(ilsp,i)
       End do
       write(nut,'(80(''-''))')

       deallocate(npch,nptar)

      End do ! ilsp

      write(nut,'(a,i9)') 'max_ch =',mch 
      write(nut,'(80(''-''))')

      End ! program convert_target_jj
