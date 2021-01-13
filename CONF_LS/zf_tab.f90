!======================================================================
!     zf_res   -->   zf_tab, zf_tau
!======================================================================
      IMPLICIT REAL(8) (A-H,O-Z)
      
      Integer :: ires=1;  Character(20) :: AF_res = 'zf_res'
      Integer :: itab=2;  Character(20) :: AF_tab = 'zf_tab'
      Integer :: itau=3;  Character(20) :: AF_tau = 'zf_tau'
      Integer :: ista=4;  Character(20) :: AF_sta = 'zf_states'
      Integer :: idat=9;  Character(20) :: AF_dat = 'zf_dat'

      Logical :: EX

      Character(85) :: AS
      Character(64) :: lab,lab1,lab2
      Character( 2) :: conv
      
      Character(64), allocatable :: LABEL(:)
      Character(2), allocatable ::  atype(:)
      Character(2) ::  btype

      REAL(8), Allocatable :: ET(:), E_exp(:), LAMDA(:)
      REAL(8), Allocatable :: FL(:),FV(:), SL(:),SV(:), AL(:),AV(:)
      REAL(8), Allocatable :: ATL(:),ATV(:)

      Integer, Allocatable :: IPT(:),IPT1(:),IPT2(:),jot(:)
     
      Real(8) :: eps_e = 1.D-8
      Real(8) :: eps_f = 1.D-26

      Real(8) :: time_au = 2.4189d-17
      Real(8) :: c_au    = 137.03599976
      Real(8) :: Ry_cm   = 109737.31568549
      Real(8) :: Ry_eV   = 13.60569172
      Real(8) :: au_eV   = 27.2113834
      Real(8) :: au_cm   = 219471.62

      Call inf_zf_tab
!----------------------------------------------------------------------      
! ... open zf_res:       

      Call Check_file(AF_res)
      Open(ires,file=AF_res,status='OLD')

      Call Read_rarg('eps_e',eps_e)
      Call Read_rarg('eps_f',eps_f)


! ... number of transitions:
       
       nt = 0; rewind(ires) 
       Do 
        read(ires,'(a)',end=1) AS
        if(AS(10:10).ne.'.') Cycle
        read(ires,'(a)',end=1) AS
        nt=nt+1        
       End do    
     1 write(*,*) ' nt = ',nt; if(nt.eq.0) Stop ' nt = 0 '

      Allocate(IPT(nt), IPT1(nt), IPT2(nt),  &
                FL(nt), FV(nt), SL(nt), SV(nt), AL(nt), AV(nt) )

      ms = 2*nt 
      Allocate(ET(ms), LABEL(ms), jot(ms), ATL(ms), ATV(ms),atype(ms), LAMDA(ms) ) 

! ... number of states:
       
       ns = 0; rewind(ires) 
       Do 
        read(ires,'(a)',end=2) AS
        if(AS(10:10).ne.'.') Cycle
        read(AS,'(i4,f14.8,2x,a)') j,E,lab
        m = 0
        Do i = 1,ns
         if(j.ne.jot(i).or.abs(E-ET(i)).gt.eps_e) Cycle
         m=i; Exit
        End do
        if(m.eq.0) then 
         ns=ns+1; ET(ns)=E; jot(ns)=j; Label(ns)=lab 
        end if
       End do    
     2 write(*,*) ' ns = ',ns; if(ns.eq.0) Stop ' ns = 0 '

! ... sort levels according to energy:

      Do i = 1,ns-1; Do j = i+1,ns
       if(ET(i).le.ET(j)) Cycle
       E=ET(i); ET(i)=ET(j); ET(j)=E
       lab=Label(i); Label(i)=Label(j); Label(j)=lab
       jj=jot(i); jot(i)=jot(j); jot(j)=jj
      End do; End do

! ... read all f-values:                              

       m = 0
       rewind(ires) 
       Do n=1,nt
     3  read(ires,'(a)',end=4) AS; if(AS(10:10).ne.'.') go to 3
        read(AS,'(i4,f14.8,2x,a)') j1,E1,lab1
        read(ires,'(i4,f14.8,2x,a)') j2,E2,lab2
        read(ires,'(a)') AS; S=0; if(AS(20:20).ne.'*') read(AS(20:),*) S
        read(ires,'(1x,a2,5x,D13.5,7x,D13.5,8x,D13.5)') btype,xsl,xfl,xal
        if(btype.eq.'E1') read(ires,'(8x,D13.5,7x,D13.5,8x,D13.5)') xsv,xfv,xav

        i1=0; i2=0
        Do i = 1,ns
         if(j1.eq.jot(i).and.abs(E1-ET(i)).le.eps_e) i1 = i
         if(j2.eq.jot(i).and.abs(E2-ET(i)).le.eps_e) i2 = i
        End do

        if(i1.eq.0) then
         write(*,*) j1,E1,trim(Lab1)
         Stop 'Cannot find this state' 
        end if

        if(i2.eq.0) then
         write(*,*) j2,E2,trim(Lab2)
         Stop 'Cannot find this state' 
        end if

        k = 0
        Do i = 1,m
         if(atype(i).ne.btype) Cycle
         if(i1.ne.IPT1(i)) Cycle
         if(i2.ne.IPT2(i)) Cycle
         k = 1; Exit
        End do

        if(k.eq.0) then
        m = m + 1
        IPT1(m) = i1; IPT2(m) = i2
        FL(m) = xfl; FV(m) = xfl; if(btype.eq.'E1') FV(m) = xfv  
        SL(m) = xsl; SV(m) = xsl; if(btype.eq.'E1') SV(m) = xsv  
        AL(m) = xal; AV(m) = xal; if(btype.eq.'E1') AV(m) = xav         
        atype(m) = btype;  LAMDA(m) = S
        end if

       End do

     4 write(*,*) ' nt = ',nt,m
       nt = m

       Close(ires)

!       Call sortI2(nt,IPT1,IPT2,IPT)

       Call sortI2(nt,IPT2,IPT1,IPT)

!----------------------------------------------------------------------
! ... read new labels if any:

      if(Icheck_file(AF_sta).ne.0) then
       open(ista,file=AF_sta)
    5 read(ista,*,err=5,end=6) lab,E           
      Do i=1,ns
       if(abs(E-ET(i)).gt.eps_e) Cycle
       LABEL(i) = lab
       Exit
      End do
      go to 5
    6 Continue
      end if

!----------------------------------------------------------------------
! ... output in zf_tab:

      Open(itab,file=AF_tab)

      ilab1=0; ilab2=0
      Do j=1,nt; i=IPT(j); i1=IPT1(i); i2=IPT2(i)
       ilab1=max(ilab1,LEN_TRIM(Label(i1)))
       ilab2=max(ilab2,LEN_TRIM(Label(i2)))
      End do

       lab1='initial'; lab2='final'
       write(itab,'(6x,a7,5x,a,a3,a,2a11,a6,2x,4a11)') 'indexes', &
            lab1(1:ilab1),' - ',lab2(1:ilab2), &
            '     SL    ', '     SV    ', '  err ',  &
            '     FL    ', '     FV    ','     AL    ', '     AV    '

      j1 = 0
      Do j=1,nt; i=IPT(j)

!       if(FL(i).lt.EPS_f.or.FV(i).lt.EPS_f) Cycle   !!!

       err = abs(SL(i)-SV(i))/(SL(i)+SV(i))*200   

       i1=IPT1(i); i2=IPT2(i)
       if(j1.ne.i2) then; write(itab,*); j1=i2; end if 

       lab1=LABEL(i1); lab2=LABEL(i2)

       write(itab,'(2i6,2x,a,2x,a,a3,a,2i4,2(1PE12.3),0Pf6.0,2x,1P4E12.3)') i1,i2, atype(i), &
            lab1(1:ilab1),' - ',lab2(1:ilab2), jot(i1), jot(i2), &
            SL(i),SV(i), abs(SL(i)-SV(i))/(SL(i)+SV(i))*200, &
            FL(i),FV(i), AL(i),AV(i)
 
      End do

      Close(itab)
!----------------------------------------------------------------------
! ... output in zf_dat:

      Open(idat,file=AF_dat)

      Do j=1,nt; i=IPT(j)

       if(FL(i).lt.EPS_f.or.FV(i).lt.EPS_f) Cycle   !!!

       if(Lamda(I).eq.0.d0) Cycle 

       i1=IPT1(i); i2=IPT2(i)
       
                                                                                                                                                
!       write(idat,'(2i5,2x,a2,f16.2,1P3E15.2)') i1,i2, atype(i), Lamda(i), SL(i), FL(i), AL(i)

       write(idat,'(2i5,2x,a2,f16.2,1P3E15.3)') i1,i2, atype(i), Lamda(i), FL(i), FV(i), AL(i)

 
      End do

      Close(itab)


!----------------------------------------------------------------------
! ...  output total Ak and life-times in zf_tau:


! ... define total decay rates:

       ATL = 0.d0; ATV = 0.d0
       Do n = 1,nt
        i=IPT2(n)
        ATL(i) = ATL(i) + AL(n)
        ATV(i) = ATV(i) + AV(n)
       End do      

       ilab = max(ilab1,ilab2)
       Call SortR(ns,ET,IPT)

       Open(itau,file=AF_tau)

       lab='State'
       write(itau,'(11a)') &
            lab(1:ilab),' 2J', &
            '       TL   ', '       TV   ',   &
            '       AL   ', '       AV   ', '      err  ',  &
            '    E_eV   ', '      E_cm   ','    E_au '

       write(itau,*)


       Do j=1,ns; i = IPT(j)

        E_eV = (ET(i)-ET(1))*au_eV
        E_cm = (ET(i)-ET(1))*au_cm

        TL = 0.d0; if(ATL(i).ne.0.d0) TL=1.d+0/ATL(i)
        TV = 0.d0; if(ATV(i).ne.0.d0) TV=1.d+0/ATV(i)
        DT=0.d0; if(TL+TV.gt.0.d0) DT = abs(ATL(i)-ATV(i))/(ATL(i)+ATV(i))*200

        lab = LABEL(i)
        write(itau,'(a,i3,2x,1P4E12.2,0Pf8.0,f10.3,f12.1,f20.8)') & 
         LAB(1:ilab), jot(i), TL,TV, ATL(i),ATV(i), DT,  E_eV, E_cm, ET(i)

       End do

      Close(itau)


      End  ! utility zf_tab            



!======================================================================
      Subroutine inf_zf_tab
!======================================================================
!     provide screen information about utility  'zf_tab'
!----------------------------------------------------------------------
       
      Implicit real(8) (a-h,o-z)

      Character(80) :: A=' '

      iarg = IARGC()
      if(iarg.gt.0) Call GETARG(1,A)
      if(A.ne.'?')  Return

      write(*,*) &
'                                                              ',&
'zf_tab creates tables for f-values and tau:                   ',&
'                                                              ',&
'   zf_res   -->   zf_tab, zf_tau                              ',&
'                                                              ',&
'Call as:   zf_tab  [eps_e=...  eps_f=...]                     ',&
'                                                              ',&
'eps_e  - tolerance for equal energies     [1.d-8]             ',&
'eps_f  - tolerance for output of f-values [1.d-6]             ',&
'                                                              ',&
'                                                              ',&
'                                                               '
      Stop ' '

      End Subroutine inf_zf_tab

