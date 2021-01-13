!======================================================================
!     H.DAT  -->  target  (bsr-format)
!======================================================================
      Use target, nconat => ictarg
      Use channels

      Implicit real(8) (A-H,O-Z)

      Character(40)  :: AF = 'H.DAT', BF
      Character( 1), external :: AL

      Allocatable :: llch(:), jjch(:), jptar(:), eval(:), aval(:,:)

      Real(8), parameter ::  Ry = 13.6056

!----------------------------------------------------------------------
! ... files:

      Call get_command_argument(1,BF)  

      if(BF.eq.'?') then
        write(*,'(/a)') 'h_targb creates "target" file based on the "H.DAT" file'
        write(*,'(/a)') '     H.DAT  -->  target.h  (BSR-format)'
        write(*,'(/a)') 'Call as:  h_targb  [h=...]'
        write(*,'(/a/)') 'default name:  h = H.DAT'
        Stop ' '
      end if

      Call Read_aarg('h',AF)
      Call Check_file(AF)

      in=1;   Open(in,file=AF,form='UNFORMATTED',status='OLD')
      isch=2; Open(isch,status='SCRATCH',form='UNFORMATTED')
      iout=3; Open(iout,file='target.h')

      ii=0; Call Read_iarg('ii',ii)  ! shift of partial-wave indexes, rear used oftion

!----------------------------------------------------------------------
! ... target states information:

      read(in) NELC, NZ, LRANG2, LAMAX, ntarg, RA, BSTO

      m=ntarg; Call Allocate_target(m)

      read(in) etarg
      read(in) ltarg
      read(in) istarg,iptarg

! ... Buttle corrections:

      read(in) ((C,k=1,3),l=1,LRANG2)   

      coupling = 'LS'; if(istarg(1).eq.0) coupling= 'JK'

!----------------------------------------------------------------------
! ... read symmetries and record useful information in scratch file:

      nlsp=0; mch=0

    1 read(in,end=2) LRGL, NSPN, NPTY, NCHAN, MNP2, MORE

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

! ... find iptarg:

!      Do i = 1,nchan; it = jptar(i)
!       if(mod(llch(i),2).eq.0) then
!        iptarg(it) = NPTY 
!       else
!        iptarg(it) = mod(NPTY+1,2) 
!       end if
!      End do

! ... record for a while:

      write(isch) LRGL, NSPN, NPTY, NCHAN
      write(isch) llch,jjch
      write(isch) jptar
      write(isch) (eval(i),i=MNP2,MNP2-4,-1)

      Deallocate (llch,jjch,jptar,eval)
      if(more.eq.1) go to 1

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
      End do

!----------------------------------------------------------------------
! ... find jk-values in JK-case:

      go to 10  ! old version

      jkch = 0
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

! ... convert iptarg to BSR format:

!      Do i = 1, ntarg
!       if(iptarg(i).eq.0) then
!        iptarg(i) = 1
!       else
!       iptarg(i) = -1
!       end if
!      End do

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

      write(iout,'(a,i3,5x,a)') 'nlsp =',nlsp,'   !   number of partial waves'
      write(iout,'(72(''-''))')
      Do i = 1,nlsp
       AF = 'no'
       ip=1-2*ipar(i)
       write(iout,'(i3.3,3i5,5x,2a10,2i5)') &
                    i+ii,lpar(i),ispar(i),ip  ! ,AF,AF,0,0
      End do
      write(iout,'(72(''-''))')

      write(iout,'(a)') 'channels:'
      write(iout,'(72(''-''))')

      Do i = 1,nlsp
       write(iout,'(i3,a,i3.3,a,i4,a,2i6)')  &
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
