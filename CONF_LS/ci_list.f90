!====================================================================
!
!     Utility    C I _ L I S T 
!
!====================================================================
!
!     creates list of states in set of j-files (l-files) 
!
!     name.c + name.j + ci_list.inp  ->  name.list
!
!     additional parameters:
!
!     mco    -   max.number of configuration to show [1]
!     eps_c  -   tolerance for weights [0.2]'           
!     unit   -   au, Ry, eV, cm [au]'                       
!     shift  -   overall shift [0.d0]'                  
!     msol   -   max.number of solutions for output [0->all]'         
!
!     Call as:    ci_list  a.j
!                 ci_list  inp=b.inp
!                 ci_list  a.j b.j  c.j ... unit=eV  msol=10 
!
!     Example ci_list.inp:
!
!     name1.j(l)
!     name2.j(l)
!     .........
!     nameN.j(l)
!     *
!     mco = 2
!     eps_c = 0.1
!     unit = eV
!     msol = 10
!     shift = 1000   ! (shift=1 means shift=-E1)
!
!     All parameters are optional    
!
!--------------------------------------------------------------------

      Use conf_LS

      Implicit double precision (A-H,O-Z)

      Integer :: mco   = 1
      Integer :: iev   = 0
      Integer :: msol  = 0
      Real(8) :: eps_c = 0.2d0      
      Real(8) :: shift = 0.d0
      Character(40) :: mode='state'
      Real(8) :: cut_off = 0.01      
      Real(8) :: AWT = 0.d0

      Character(1) :: ctype='j'

      Real(8), Allocatable :: EE(:), C(:), gfact(:)
      Real(8), Allocatable :: CC(:,:)

      Integer, Allocatable :: IP(:), IPE(:), JJ(:), kco(:), isol(:)
      Integer, Allocatable :: IPC(:,:)

      Character(64),Allocatable :: Lab(:,:), files(:)
      Character( 6),Allocatable :: terms(:)

      Character(64)  :: Label1,Label2,AF
      Character(200) :: AS
      Character(2)   :: unit='au'
      Character(6)   :: term1

      Integer :: inp= 5; Character(64) :: AFi = 'ci_list.inp'
      Integer :: nuc=11; Character(64) :: AFc = ' '
      Integer :: nuj=12; Character(64) :: AFj = ' '
      Integer :: nur=13; Character(64) :: AFr = 'ci_list'
      Integer :: nub=14; Character(64) :: AFb = ' '
      Integer :: nua=21; ! scratch file

      Logical :: EX

      Call inf_ci_list

!----------------------------------------------------------------------
! ... input parameters:

      Call Read_aarg('inp',AFi)

      if(Icheck_file(AFi).ne.0) then

       Open(inp,file=AFi)

       Call Def_files_inp

       Call Read_rpar(inp,'eps_c',eps_c)
       Call Read_rpar(inp,'shift',shift)
       Call Read_ipar(inp,'mco'  ,mco  )
       Call Read_ipar(inp,'msol' ,msol )
       Call Read_apar(inp,'unit' ,unit )
       Call Read_apar(inp,'out'  ,AFr  )
       Call Read_rpar(inp,'awt'  ,AWT  )

       Call Read_apar(inp,'AFb',AFb)
       Call Read_apar(inp,'mode',mode)
       Call Read_rpar(inp,'cut_off',cut_off)

      end if

      Call Def_files_arg

      Call Read_rarg('eps_c',eps_c)
      Call Read_rarg('shift',shift)
      Call Read_iarg('mco'  ,mco  )
      Call Read_iarg('msol' ,msol )
      Call Read_aarg('unit' ,unit )
      Call Read_aarg('out'  ,AFr  )
      Call Read_rarg('awt'  ,AWT  )

      Call Read_aarg('AFb',AFb)
      Call Read_aarg('mode',mode)
      Call Read_rarg('cut_off',cut_off)
 
      if(ifile.eq.0) Stop 'no input file to process'

      Open(nua,status='SCRATCH',form='UNFORMATTED') 

!----------------------------------------------------------------------
! ... cycle over j-files:

      ilabel = 1
      nsol=0

      Do jfile=1,ifile
       AFj = files(jfile)

       open(nuj,file=AFj,status='OLD')
       i=LEN_TRIM(AFj); ctype=AFj(i:i); AFc=AFj(1:i-1)//'c'
       open(nuc,file=AFc,status='old')

! ... find conversion factor:
       
      if(nsol.eq.0) then
       rewind(nuj)
       read(nuj,'(a)') AF
       read(AF(14:19),*) Z
       Call Conv_au (Z,AWT,au_cm,au_eV,6)
       Select case(unit)
        case('au'); conv = 1.d0; iev=0
        case('Ry'); conv = 2.d0; iev=1
        case('eV'); conv = au_eV; iev=2
        case('cm'); conv = au_cm; iev=3
       End select
      end if

! ... find the labels for config.s:

      ncfg=Idef_ncfg(nuc);  if(ncfg.eq.0) Stop ' ncfg = 0'

      Call R_label(nuc,1,0)      

      Allocate(C(ncfg),IP(ncfg))

      Call R_term(nuc)

      Do ic=1,ncfg
       i=LEN_TRIM(LABEL(ic)); if(i.gt.ilabel) ilabel=i
      End do      

! ... read solutions:

      ns = Idef_st(nuj);  if(ns.eq.0) Cycle

      Rewind(nuj)
      Do is = 1,ns

       Call Read_sol(ctype,nuj,ncfg,C,Label2,es,jot,gf,js)
       Call Sorta(ncfg,C,IP)

! ... define leading weights:

       nco=1                                 
       if(ncfg.gt.1) then
       Do i=2,mco
        if(abs(C(IP(i))).lt.eps_c.and.i.gt.2) Exit; nco=i
       End do
       end if

! ... leading term from label:

       k = IP(1);  Label1=Label(k)
       j = LEN_TRIM(Label1)
       i = INDEX(Label1,'_',BACK = .TRUE.)
       term1 = Label1(i+1:j)

! ... save information:
 
       write(nua) es,jot-1,nco,term1,gf,js,AFJ
       Do i=1,nco
        j=IP(i); write(nua) C(j),LABEL(j)
       End do
       nsol=nsol+1         

      End do   

      Close(nuj)

      Deallocate(C,IP)

      End do  ! jfile

      write(*,*) ' nsol = ',nsol
      Deallocate(files)

!----------------------------------------------------------------------
! ... restore the solutions:             

      ns = nsol
      Allocate(EE(ns),IPE(ns),kco(ns),JJ(ns),CC(ns,mco),IPC(ns,mco))
      Allocate(gfact(ns),files(ns),isol(ns),terms(ns))
      Allocate(LAB(ns,mco))

      rewind(nua)
      Do i=1,ns
       read(nua) EE(i),JJ(i),nco,terms(i),gfact(i),isol(i),files(i)
       kco(i)=nco
       Do j=1,nco
        read(nua) CC(i,j),LAB(i,j)
       End do
      End do

! ... sort of results:

      Call SortR(ns,EE,IPE)

      if(shift.eq.1.d0) shift = EE(IPE(1)) 
      
      Do i=1,ns; EE(i)=(EE(i) - shift)*conv; End do


!----------------------------------------------------------------------
! ... print the results:

      if(ifile.eq.1) AFr=trim(AFj)//'_list'

      if(LEN_TRIM(AFb).gt.0) then
       Open(nub,file=AFb)       
      else
       nub=0 
      end if

      Open(nur,file=AFr)
      ms=ns; if(ms.gt.msol.and.msol.ne.0) ms=msol

      ilabel1 = 1
      Do i=1,ms; ilabel1=max(LEN_TRIM(LAB(i,1)),ilabel1); End do

      ilabel2 = 1
      Do i=1,ifile; ilabel2=max(LEN_TRIM(files(i)),ilabel2);  End do

!----------------------------------------------------------------------

      Do ii=1,ms;  write(AS,'(i3,a)') ii,'.'

       i=IPE(ii);  nco=kco(i);  LABEL1=LAB(i,1)

       k=6; write(AS(k:),'(a)') LABEL1(1:ilabel1)
       k=k+ilabel1
       jot = JJ(i)
       if(jot/2*2.eq.jot) then
         write(AS(k:),'(i4)') jot/2
       else
         write(AS(k:),'(i2,a2)') jot,'/2'
       end if  
       k = k + 6

       Select case(iev)
        case(0,1); write(AS(k:),'(f16.8)') EE(i); k=k+20
        case(2);   write(AS(k:),'(f10.3)') EE(i); k=k+12 
        case(3);   write(AS(k:),'(f10.1)') EE(i); k=k+12  
       End Select

       write(AS(k:),'(f8.3)') gfact(i); k=k+13  

       write(AS(k:),'(f5.1)') CC(i,1)**2*100; k=k+12  

       Do j=2,nco
        LABEL1 = LAB(i,j)
        write(AS(k:),'(f4.1,2x,a)') CC(i,j)**2*100,LABEL1(1:ilabel)
        k=k+5+ilabel+3
       End do
       write(nur,'(a,i10)') trim(AS), isol(i)

       if(nub.eq.0) Cycle

       k=1
       write(AS(k:),'(a,a)') 'cfile  ',trim(files(i))  
       k=k+ilabel2+10
       write(AS(k:),'(2i5)') isol(i),jot  
       k=k+13

       if(mode.eq.'state') then
        write(AS(k:),'(a,i3.3)') 'state_',ii  
        k=k+11
       else
        write(AS(k:),'(a)') LABEL1(1:ilabel1)
        k=k+ilabel1+2
       end if
       write(AS(k:),'(f8.5)') cut_off

       write(nub,'(a)') trim(AS)


      End do


CONTAINS

!=======================================================================      
      Subroutine Def_files_inp
!=======================================================================

      ifile = 0
      rewind(inp)
      Do 
       read(inp,'(a)',end=1) AF
       if(LEN_TRIM(AF).eq.0) Cycle
       if(INDEX(AF,'=').ne.0) Cycle
       if(AF(1:1).eq.'*') Exit
       ifile=ifile+1
      End do
    1 if(ifile.eq.0) Return

      Allocate(files(ifile))
      i = 0
      rewind(inp)
      Do 
       read(inp,'(a)',end=2) AF
       if(LEN_TRIM(AF).eq.0) Cycle
       if(INDEX(AF,'=').ne.0) Cycle
       if(AF(1:1).eq.'*') Exit
       i=i+1; files(i)=AF
      End do
    2 Continue

      End Subroutine Def_files_inp


!=======================================================================      
      Subroutine Def_files_arg
!=======================================================================

      ifile = 0
      Do i=1,iarg
       Call GETARG(i,AF)
       if(INDEX(AF,'=').ne.0) Cycle
       ifile=ifile+1
      End do

      if(ifile.eq.0) Return

      Allocate(files(ifile))
      j = 0
      rewind(inp)
      Do i=1,iarg
       Call GETARG(i,AF)
       if(INDEX(AF,'=').ne.0) Cycle
       j=j+1; files(j)=AF
      End do

      End Subroutine Def_files_arg



      End ! utility CI_LIST_NIST




!======================================================================
      Subroutine Read_sol(ctype,nu,nc,C,Label,E,jot,gv,nn)
!======================================================================
!
!     read one next solution from c-,l-,j- or b-files 
!
!     nu    -  file unit
!     ctype -  type of file (c,l,j,b)
!     nc    -  dimention of solution
!     C     -  solution
!     Label -  spectroscopic notation
!     E     -  energy 
!     jot   -  (2J+1) for j-type calculations
!
!     the first line in files 'l,j,b' is supposed to have already 
!     been read !
!----------------------------------------------------------------------


      Implicit none

      Character(1), Intent(in) :: ctype
      Character(64), Intent(out) :: Label
      Integer(4), Intent(in) :: nu,nc 
      Integer(4), Intent(inout) :: jot,nn 
      Real(8), Intent(out) :: E,gv
      Real(8),Dimension(nc),Intent(out) :: C
      Integer(4) :: i
 
      Character(80) :: AS

      Select Case(ctype)

       Case('b')

        read(nu,'(7x,a64)') Label
        read(nu,'(D15.8)') E  
        read(nu,'(5D15.8)') C

      Case('c')
	  
       Call R_expn(nu,nc,C)
       Call Idef_LabelC(nu,0,0,Label)        
       rewind(nu); read(nu,'(a)') AS; read(AS,'(15x,f16.8)') E
       jot=0
       i=INDEX(AS,'='); if(i.gt.0) read(AS(i+1:),*) jot; jot=jot+1

      Case('j') 

      1 read(nu,'(a80)',end=2) AS
        i=INDEX(AS,'g_J=')
        if(i.ne.0) then
         read(AS(i+4:),*) gv 
         read(nu,*) nn,E, Label
         read(nu,'(7F11.8)') C
        elseif(AS(5:5).eq.'J') then
         read(AS(8:),*) jot; jot=jot+1; go to 1
        else
         go to 1
        end if
        Return
      2 write(*,*) 'Read_sol: read the file end '
        Stop

       Case('l') 

      3 read(nu,'(a80)',end=4) AS
	if(AS(14:14).eq.'.') then
         read(AS,*) nn, E, Label
         read(nu,'(7F11.8)') C
        else
         go to 3
        end if
        Return
      4 write(*,*) 'Read_sol: read the end of file'
        Stop

      Case Default
	    
	Stop ' Read_sol: unknown case '

      End Select
 
      End Subroutine Read_sol


!======================================================================
      Subroutine inf_ci_list
!======================================================================
!     provide screen information about utility  'ci_list'
!----------------------------------------------------------------------
       
      Implicit real(8) (a-h,o-z)

      Character(80) :: A=' '

      iarg = command_argument_count() 
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,A)
 
      if(A.ne.'?')  Return

      write(*,'(a)') &
'                                                              ',&
'ci_list creates list of states in set of j-files (l-files)    ',&
'                                                              ',&
'name.c + name.j + [ci_list.inp]  ->  name.j_list              ',&
'                                                              ',&
'additional parameters:                                        ',&
'                                                              ',&
'mco    -   number of leading configurations to show [1]       ',&
'eps_c  -   tolerance for weights [0.2]                        ',&
'unit   -   au, Ry, eV, cm [au]                                ',&
'shift  -   overall shift [0.d0]                               ',&
'msol   -   max.number of solutions for output [0->all]        ',&
'                                                              ',&
'Call as:    ci_list  a.j                                      ',&
'            ci_list  inp=b.inp                                ',&
'            ci_list  a.j b.j  c.j ... unit=eV  msol=10        ',&
'                                                              ',&
'Example of input file [ci_list.inp]:                          ',&
'                                                              ',&
'name1.j(l)                                                    ',&
'name2.j(l)                                                    ',&
'.........                                                     ',&
'nameN.j(l)                                                    ',&
'                                                              ',&
'mco = 2                                                       ',&
'eps_c = 0.1                                                   ',&
'unit = eV                                                     ',&
'msol = 10                                                     ',&
'shift = 1000   ! (shift=1 means shift=-E1)                    ',&
'                                                              ',&
'All parameters are optional                                   ',&
'                                                               '
      Stop ' '

      End Subroutine inf_ci_list

