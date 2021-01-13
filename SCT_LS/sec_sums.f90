!=====================================================================
!     UTILITY  sec_sums
!=====================================================================
!
!     provides sums of cross sections from given target state
!
!     Arguments: is   -  initial state   (default: is=1)
!
!     Input files:
!
!      target  - description of scattering model
!      zarm.om - collection for collision strenths
!
!     Output files:
!
!      sum_sec_##, where ## stands for 'is'.      
!
!     Example of call:  sec_sums is=.. index=.. i16=.. om=..  out=..
!
!---------------------------------------------------------------------

      Use target

      Implicit real(8) (A-H,O-Z)

      Real(8), allocatable :: e(:),om(:),y(:) 

      Logical EX

      Character(3)  :: AT
      Character(40) :: AF,BF
      Character(80) :: AS

      Real(8) :: Ry = 13.6057

      Integer :: nut=1;  Character(20) :: AF_t   = 'target'
      Integer :: nuo=2;  Character(20) :: AF_om  = 'zarm.om'

      Integer :: out=4;  Character(20) :: AF_out = 'sum_sec_01'

      Integer, allocatable :: index(:)

      Call Inf_sec_sums

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(AF_t)
      Open(nut,file=AF_t)
      Call R_target(nut)
      np=ntarg; Call Read_ipar(nut,'np',np)
      ni=ntarg; Call Read_ipar(nut,'ni',ni)
      Close(nut)

      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion
      E1 = etarg(1); etarg = (etarg-E1)*2.0

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

! ... define transition under consideration:

      is = 1;   Call Read_iarg('is',is)

      i16 = 1;  Call Read_iarg('i16',i16)

      Allocate(index(ntarg))

      Call Read_aarg('index',AS)

      Call Read_iarr_string(AS,ii,index)

! ... statistical weight for initial state:

      g=iabs(IStarg(is))*(2.0*Ltarg(is)+1)
      if(IStarg(is).eq.0) g=Ltarg(is)+1

!----------------------------------------------------------------------
! ... define energies:

      Call Read_aarg('om',AF_om)
      Open(nuo,file=AF_om,status='OLD')

      me=0
      rewind(nuo)
    1 read(nuo,*,end=2) x,ns
      read(nuo,*) (S,ix=1,ns)
      me = me + 1
      go to 1
    2 Continue
      write(*,*) ' me = ',me

! ... allocations:

      ntr = ntarg*(ntarg+1)/2
      Allocate( y(ntr),E(me),OM(me) );  om=0.d0

!----------------------------------------------------------------------
! ... read data:
      
      ne=0
      rewind(nuo)
    5 read(nuo,*,end=15) x,ns
      read(nuo,*) (y(i),i=1,ns) 
      if(x.le.etarg(is)) go to 5

      ie=0
      Do i=1,ne; if(x.ne.e(i)) Cycle; ie=i; Exit; End do
      if(ie.eq.0) then
       ne=ne+1; if(ne.gt.me) Stop 'ne > me';  ie=ne
      end if
      e(ie)=x; S = 0.d0

      i1=is
      Do j2=1,ii;  i2=index(j2)

       if(x.le.etarg(i2)) Cycle
       if(ion.ne.0.and.i1.eq.i2) Cycle
       if(i1.gt.i2) Cycle
       if(i2.gt.np.and.i1.gt.ni) Cycle

! ... position of transition in arrays:

       itr=Index_TR(ion,i1,i2,np,ni)

       if(itr.eq.0) Cycle

       es = e(ie)-etarg(is)

       S = S + y(itr)/g/es * 3.1415926

      End do ! over i2

      om(ie) = S

      go to 5

   15 Continue

! ... ordering the data:

      Do i=1,ne-1; Do j=i+1,ne
       if(e(i).gt.e(j)) then
        S=e(i); e(i)=e(j); e(j)=S
        S=om(i); om(i)=om(j); om(j)=S
       end if
      End do; End do

! ... output data:

      write(AF_out,'(a,i2.2)') 'sum_sec_',is
      Call Read_aarg('out',AF_out)
      open(out,file=AF_out)

      BF  =  ' in ao^2'                     
      if(i16.gt.0) write(BF,'(a,i2.2,a)') ' in 10^-',i16,' cm^2'
      write(out,'(7x,a,10x,a,12x,a,15x,a,a,100i3)') &
               'eV','sigma','Ry',trim(BF),'   indexes: ', &
               index(1:ii) 
    
      Do ie=1,ne
       if(om(ie).eq.0.d0) Cycle
       s = om(ie)
       if(i16.gt.0) s = s * 0.28003 * 10**(i16-16)     ! in 10^-i16
       write(out,'(f14.8,e16.8,f14.6)') e(ie)*Ry,s,e(ie)   
      End do             

      Close(out)

      End  ! UTILITY  sec_sums


!======================================================================
      Subroutine inf_sec_sums
!======================================================================
!     provide screen information about sec_top_TM utility
!----------------------------------------------------------------------
       
      Character :: A=' '

      Call get_command_argument(1,A)  
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                      ',&
'     SEC_SUMS provides sums of cross sections from given target state ',&
'                                                                      ',&
'     Input files:                                                     ',&
'                                                                      ',&
'      target  - description of scattering model                       ',&
'      zarm.om - collection for collision strenths                     ',&
'                                                                      ',&
'     Output files:                                                    ',&
'                                                                      ',&
'      sum_sec_##, where ## stands for ''is''                          ',&
'                                                                      ',&
'     Call as:  sec_sums is=.. index=.. i16=.. om=..  out=..           ',&
'                                                                      ',&
'      is [1]  - index of initial state                                ',&
'      index   - indexes of the final states, e.g. 2,3,5-7,9 ...       ',&
'      i16 [0] - units for cross sections                              ',&
'                0 - in a_0^2                                          ',&
'               >0 - in 10^-i16  cm^2                                  ',&
'      om [zarm.om]  - name of omega file, e.g. zarm.omb, zarm.omb_top ',&
'      out [sec_sum_##] - name of output file                          ',&
'                                                                      '
      Stop ' '

      End Subroutine inf_sec_sums
