!======================================================================
!     UTILITY      S E C _ P A R B  
!
!     zarm.omb_par --> pa_nnn_nnn
!
!======================================================================
!
!     generate partial cross sections for given transitions 
!     (based on the information in 'zarm.omb_par')
!
!     Call as:  sec_parb ilsp1= ilsp2=  itr1= itr2=  par=  i16=  
!
!---------------------------------------------------------------------

      Use target; Use channels

      Implicit real(8) (a-h,o-z)

      Real(8), allocatable ::  om(:), e(:), fom(:,:)
      Integer, allocatable ::  ipe(:), ind(:)
      Character(4), allocatable :: term(:)
      Character(10) :: aaa

      Integer :: ke = 10000  !  initial number of energies 
      Character(1), external :: AL

!----------------------------------------------------------------------
! ... files:

      Integer :: nup=11;  Character(20) :: targ  = 'target'
      Integer :: nut=12;  Character(20) :: par   = 'zarm.omb_par'
      Integer :: nuo=14;  Character(20) :: out   = 'pa_nnn_nnn'

      Integer :: nua=99;  ! scratch file

      Call Inf_par_omb
 
!----------------------------------------------------------------------
! ... target information:

      Call Check_file(targ)
      Open(nup,file=targ)
      Call R_target(nup)
      Call R_channels(nup)
      np=0; Call Read_ipar(nup,'np',np)
      ni=0; Call Read_ipar(nup,'ni',ni)
      Close(nup)
      e1=etarg(1); Do i=1,ntarg; etarg(i)=(etarg(i)-e1)*2; End do
      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion

      Z = nz
      AWT = 0.d0;  Call Read_rarg('AWT',AWT)
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

      itr1=1; Call Read_iarg('itr1',itr1)
      itr2=1; Call Read_iarg('itr2',itr2)

      ilsp1=1; Call Read_iarg('ilsp1',ilsp1)
      ilsp2=1; Call Read_iarg('ilsp2',ilsp2)

      i16=0; Call Read_iarg('i16',i16)

!----------------------------------------------------------------------
! ... find energies:

      Call Read_aarg('par',par)
      Call Check_file(par)
      Open(nut,file=par)

      me = ke; Allocate(e(me))

      ne=0; np=0; ni=0; mom=0
    1 read(nut,*,end=2) e1,ilsp,nom,iop,jop,i1,i2
      read(nut,*) (S,i=1,nom)
      if(jop.lt.itr2) go to 1

      if(np.eq.0) np = i1       
      if(ni.eq.0) ni = i2       
      if(np.ne.i1) Stop 'diferent np'
      if(ni.ne.i2) Stop 'diferent ni'

      if(nom.gt.mom) mom=nom

      ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

      if(ie.eq.0) then; ne=ne+1; e(ne)=e1;  end if

      if(ne.eq.me) then
       open(nua,form='UNFORMATTED',status='SCRATCH')
       rewind(nua);   write(nua) (e(i),i=1,ne)
       Deallocate(e); me=me+ke; Allocate(e(me))
       rewind(nua);   read(nua) (e(i),i=1,ne)
      end if       

      go to 1
    2 write(*,'(a,i5)') 'ne =',ne
      if(ne.eq.0) Stop 'nothing to do !'

      Allocate(fom(nlsp,ne), om(mom) )

      fom = 0.d0

!----------------------------------------------------------------------
! ... read partial OM:

      rewind(nut)
    3 read(nut,*,end=4) e1,ilsp,nom,iop,jop
      read(nut,*) (om(i),i=1,nom)
      if(jop.lt.itr2) go to 3

      ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do
      if(ie.eq.0) Stop 'ie = ? '

      ntr = iop*(iop+1)/2
      if(ion.ne.0)     ntr = iop*(iop-1)/2
      ktr = ntr + (jop-iop)*ni
      if(ktr.ne.nom) then 
       write(*,*) 'ilsp = ',ilsp
       write(*,*) 'nom  = ',nom
       write(*,*) 'iop,jop = ',iop,jop
       write(*,*) 'ntr,ktr = ',ntr,ktr
       Stop 'nom should be equal ktr - ?'
      end if 

      itr = (itr2-1)*itr2/2 + itr1
      if(ion.ne.0) itr = itr - itr2 + 1
      if(itr2.gt.iop) itr = ntr + (itr2-np-1)*ni + itr1
      if(itr.lt.0) Stop 'itr < 0'
      if(itr.gt.nom) Stop 'itr > nom'

      fom(ilsp,ie) = om(itr)
      go to 3
    4 Close(nut)

      Deallocate(om); Allocate(om(ne))
      Do ie=1,ne; om(ie)=SUM(fom(ilsp1:ilsp2,ie)); End do

      if(i16.ge.-1) then

       g=iabs(IStarg(itr1))*(2.0*Ltarg(itr1)+1)
       if(IStarg(itr1).eq.0) g=Ltarg(itr1)+1

       s =  3.1415926/g                !  in a_0^2
       if(i16.eq.-1) S = 1.0/g         !  in pi a_0^2
       if(i16.gt.0)  S = S * 0.28003 * 10.0**(i16-16)      

       Do ie=1,ne
        om(ie) = om(ie) * s / (e(ie)-etarg(itr1))
        fom(:,ie) = fom(:,ie) * s / (e(ie)-etarg(itr1))
       End do

      end if

      aaa = 'omega'
      if(i16.eq.-1) aaa = 'pi*a_0^2'
      if(i16.eq. 0) aaa = 'a_0^2'
      if(i16.gt. 0) write(aaa,'(a,i2)') '10^-',i16

!----------------------------------------------------------------------
! ... output:

! ... define terms for partial waves:

      Allocate(term(nlsp),ind(nlsp))
      klsp = 0
      Do i=1,nlsp
       if(ispar(i).eq.0) then
        i1=lpar(i)/10; i2=mod(lpar(i),10)
        write(term(i),'(a1,2i1,a1)') '_', i1,i2
       else
        write(term(i),'(a1,i1,a1,a1)') '_', &
                     iabs(ispar(i)),AL(lpar(i),2)
       end if

       if(ipar(i).eq. 1) term(i)(4:4)='e'
       if(ipar(i).eq.-1) term(i)(4:4)='o'

       S = SUM(fom(i,:))
       if(S.gt.0.d0.and.i.ge.ilsp1.and.i.le.ilsp2) then
        klsp=klsp+1; ind(klsp)=i
       end if

      End do

      if(klsp.eq.0) Stop 'klsp = 0, number of partial waves for output'
      Allocate(ipe(ne))
      Call SortR(ne,E,IPE)

      write(out,'(a,i3.3,a,i3.3,a)') 'pa_',itr1,'_',itr2 

      open(nuo,file=out)

      write(nuo,'(6x,a,13x,a,13x,a,8x,50(3x,a,i3.3,a,4x))')  &
                 'E','Ry','SUM', ('p',ind(i),term(ind(i)),i=1,klsp),aaa

      Do je=1,ne; ie=IPE(je); ee = e(ie)-etarg(itr1)
       write(nuo,'(258E15.5)') ee*Ry,ee, om(ie),(fom(ind(i),ie),i=1,klsp)
      End do

      Close(nuo)

      End   !  program sec_parb


!======================================================================
      Subroutine inf_par_omb
!======================================================================
!     provide screen information about sec_top_TM utility
!----------------------------------------------------------------------
       
      Character :: A=' '

      Call get_command_argument(1,A)  
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                          ',&
'     SEC_PARB generates partial cross sections for transition itr1->itr2, ',&
'     based on the information in "zarm.omb_par"                           ',&
'                                                                          ',&
'            zarm.omb_par + target  =>  par_i1_i2                          ',&
'                                                                          ',&
'     Call as:  sec_parb ilsp1=.. ilsp2=.. itr1=.. itr2=..  [par=.. i16=..]',&
'                                                                          ',&
'     ilsp1,ilsp2 - range of partial waves for output                      ',&
'                                                                          ',&
'     par - can be used to cahnge default name for par-file                ',&
'                                                                          ',&
'     i16 - control the output:                                            ',&
'           0 - sigma in a.u.  (default)                                   ',&
'          >0 - sigma in 10^-i16 cm^2                                      ',&
'                                                                          '

      Stop ' '

      End Subroutine inf_par_omb

