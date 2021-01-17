!======================================================================
!     H.DAT  -->  target  (BSR-format)
!
!     for case when when different partial H.DAT(or h.nnn) are recorded 
!     in the different directories nnn
!======================================================================
      Use target, nconat => ictarg
      Use channels

      Implicit real(8) (A-H,O-Z)

      Allocatable  llch(:), jjch(:), jptar(:), eval(:), aval(:,:), &
                   etarg1(:), ltarg1(:), istarg1(:), iptarg1(:)

      Character(40) :: AF, BF 
      Character(1), external :: AL

      Real(8), parameter ::  Ry = 13.6057

      ii=0; Call Read_iarg('ii',ii)

      klsp=0; Call Read_iarg('klsp',klsp)

      Call get_command_argument(1,BF)  
      if(BF.eq.'?') then
        write(*,'(/a)') 'h_targn creates "target" file based on the "nnn/H.DAT" ', &
                         ' or nnn/h.nnn files plased in the set of sub-directories nnn'
        write(*,'(/a)') '   {nnn/H.DAT} or {nnn/h.nnn}  -->  target.h  (BSR-format)'
        write(*,'(/a)') 'Call as:  h_targn  klsp=...'
        write(*,'(/a/)') 'klsp - number of partial waves'
        Stop ' '
      end if

!----------------------------------------------------------------------
! ... files:

      isch=2; Open(isch,status='SCRATCH',form='UNFORMATTED')
      iout=3; Open(iout,file='target.hhh')

!----------------------------------------------------------------------    
! ... Cycle over different H-files:

      istart= 1

      Do ih = 1,klsp

       write(AF,'(i3.3,a,i3.3)') ih,'/H.DAT'
       write(BF,'(i3.3,a,i3.3)') ih,'/h.',ih
       i = Icheck_file(AF)
       j = Icheck_file(BF)
       if(Icheck_file(AF).gt.0) then
        Open(in,file=AF,form='UNFORMATTED',status='OLD')
       elseif(Icheck_file(BF).gt.0) then
        Open(in,file=BF,form='UNFORMATTED',status='OLD')
       else
        Cycle
       end if

!----------------------------------------------------------------------
! ... target states information:

      if(istart.eq.1) then

      read(in) NELC, NZ, LRANG2, LAMAX, ntarg, RA, BSTO

      m=ntarg; Call Allocate_target(m)
      Allocate(etarg1(ntarg), ltarg1(ntarg), istarg1(ntarg), iptarg1(ntarg))

      read(in) etarg
      read(in) ltarg
      read(in) istarg,iptarg
 
      nlsp=0; mch=0; istart=0

! ... Buttle corrections:

      read(in) ((C,k=1,3),l=1,LRANG2)   

      coupling = 'LS'; if(istarg(1).eq.0) coupling= 'JK'

      else

      rewind(in)
      read(in) NELC1, NZ1, LRANG3, LAMAX1, ntarg1, RA1, BSTO1

      if(NELC.ne.NELC1)    Stop 'NELC <> NELC1'
      if(NZ  .ne.NZ1  )    Stop 'NZ   <> NZ1'
      if(LAMAX.ne.LAMAX1)  Stop 'LAMAX  is different'
      if(ntarg.ne.ntarg1)  Stop 'ntarg  is different'
      if(RA1.ne.RA)        Stop 'RA is different'
      if(BSTO.ne.BSTO1)    Stop 'RA is different'

      read(in) etarg1
      read(in) ltarg1
      read(in) istarg1,iptarg1      

       Do N=1,ntarg
       if(etarg1(N).ne.etarg(N)) then 
        write(*,*) N, etarg1(N), etarg(N)
        write(*,*) ' ENAT <> ENAT1 for  ',AF; !Stop
       end if
       if(ltarg1(N).ne.ltarg(N)) then 
        write(*,*) N, ltarg1(N), ltarg(N)
        write(*,*) ' LAT <> LAT1 for  ',AF; !Stop
       end if
       if(istarg1(N).ne.istarg(N)) then 
        write(*,*) N, istarg1(N), istarg(N)
        write(*,*) ' ISAT <> ISAT1 for  ',AF; !Stop
       end if
       if(iptarg1(N).ne.iptarg(N)) then 
        write(*,*) N, iptarg1(N), iptarg(N)
        write(*,*) ' iptar <> iptar1 for  ',AF; ! Stop
       end if
       End do

       read(in) ((C,k=1,3),l=1,LRANG3)   

      end if

!----------------------------------------------------------------------
! ... read symmetries and record useful information in scratch file:

      read(in) LRGL, NSPN, NPTY, NCHAN, MNP2, MORE

      nlsp=nlsp+1; if(nchan.gt.mch) mch=nchan
 
      read(in) NCONAT

      Allocate(llch(NCHAN), jjch(NCHAN), jptar(NCHAN), eval(MNP2))

      read(in) llch,jjch

! ... asymptotic coefficients:

      read(in) C

! ... (N+1)-electron eigenvalues:

      read(in) eval 

! ... surface amplitudes:

      read(in) C

! ... find pointer 'channel --> target'

      i1=1
      Do it=1,ntarg
       if(NCONAT(it).eq.0) Cycle
       i2=i1+NCONAT(it)-1
       Do i=i1,i2; jptar(i)=it; End do
       i1=i2+1
      End do

! ... record for a while:

      write(isch) LRGL, NSPN, NPTY, NCHAN
      write(isch) llch,jjch
      write(isch) jptar
      write(isch) (eval(i),i=MNP2,MNP2-4,-1)

      Deallocate (llch,jjch,jptar,eval)

      Close(in)

      End do

!----------------------------------------------------------------------
! ... read the whole information from scratch file:

   2  Call Allocate_channels    
      Allocate(aval(nlsp,5))
      rewind(isch)
      Do i=1,nlsp
       read(isch) lpar(i), ispar(i), ipar(i), nch(i)
       read(isch) (lch(i,j),j=1,nch(i)),(jkch(i,j),j=1,nch(i))
       read(isch) (iptar(i,j),j=1,nch(i))
       read(isch) (aval(i,j),j=1,5)

       write(*,'(5i6,5F16.5)') i,lpar(i), ispar(i), ipar(i), nch(i), aval(i,1:5)

      End do

!----------------------------------------------------------------------
! ... find jk-values in JK-case:

      go to 10   ! original way to define jk-values
      jkch=0
      if(coupling.eq.'JK') then
       Do i = 1,nlsp 
        K_min=iabs(lpar(i)-1)+1
        K_max=iabs(lpar(i)+1)+1 
        Do j = 1,nch(i)
         ll = 2*lch(i,j)+1; JT = ltarg(iptar(i,j))+1
         Do KK = K_min,K_max,2
          if(ITRI(JT,ll,kk).eq.0) Cycle
          m = 0
          Do jj = 1,j-1
           if(lch(i,j).eq.lch(i,jj).and. &
              iptar(i,j).eq.iptar(i,jj).and. &
              jkch(i,jj).eq.kk) m=1
           if(m.eq.1) Exit
          End do
          if(m.eq.1) Cycle
          jkch(i,j) = kk; Exit
         End do
        End do
       End do
      end if
   10 Continue

!----------------------------------------------------------------------
! ... record 'target' file:

      write(iout,'(a)') '   e + ...   '
      write(iout,'(72(''-''))')

      write(iout,'(a,a2,a,a2,a)') 'coupling = ',coupling,'    !   ', &
                                   coupling,' - case'    
      write(iout,'(a,i3,a)') 'nz = ',nz,'         !   nuclear charge'    
      write(iout,'(a,i3,a)') 'nelc =',nelc, &
                             '        !   number of electrons'
      write(iout,'(72(''-''))')

      write(iout,'(a,i3,a)') 'ntarg =',ntarg, &
                             '       !   number of target states'
      write(iout,'(72(''-''))')
      Do i = 1,ntarg
       write(AF,'(a,i3.3)') 'targ_',i
       write(BF,'(a,i3.3)') 't',i
       E_Ry = 2*(Etarg(i)-Etarg(1))
       E_eV = E_Ry * Ry
       write(iout,'(a,2x,a,3i5,F16.8,2i3,F10.6,f11.5)') &
         trim(BF),trim(AF),ltarg(i),istarg(i),iptarg(i),Etarg(i),0,0,E_Ry,E_eV
      End do
      write(iout,'(72(''-''))')

      write(iout,'(a,i3,5x,a)') 'nlsp =',nlsp, &
                             '   !   number of partial waves'
      write(iout,'(72(''-''))')

      Do i = 1,nlsp
       AF = 'no'
       i3=1-2*ipar(i)
       write(iout,'(i3.3,3i5,5x,2a10,2i5)') &
                    i+ii,lpar(i),ispar(i),i3,AF,AF,0,0

      End do
      write(iout,'(72(''-''))')

      write(iout,'(a)') 'channels:'
      write(iout,'(72(''-''))')

      Do i = 1,nlsp
       write(iout,'(i3,a,i3.3,a,i6,a,2i10)')  &
                    i+ii,'.  ',i+ii,'  nch =',nch(i),'  nc =',0,0 
       Do j = 1,nch(i)
        write(iout,'(2x,a1,a1,3i5,i10,i5)') &
              'k',AL(lch(i,j),1),lch(i,j),iptar(i,j),j,0,jkch(i,j) 
       End do
       write(iout,'(72(''-''))')
      End do

      write(iout,'(a,i5,T20,a)')   'max_ch =',mch,   '!  max. number of channels'
      write(iout,'(a,f8.3,T20,a)') 'RA = ',RA,       '!  bouder radius'    
      write(iout,'(a,i2,T20,a)')   'lamax = ',lamax, '!  max. multipole index'    

      End  ! utility H_targb
