!======================================================================
      MODULE data_list
!======================================================================
!     Containes the allocatable data with three real and six integer
!     parameters
!----------------------------------------------------------------------

      Implicit none

      Integer(4) :: ndata = 0       ! number of determinants
      Integer(4) :: mdata = 0       ! current dimentsion of list
      Integer(4) :: idata = 2**10   ! supposed max. dimentsion

      Integer(4) :: nus   = 99      ! unit for reallocation
    
      Integer(4), Allocatable :: jjpar(:),ip(:),l1(:),l2(:),jk1(:),jk2(:)   
      Real(8), Allocatable :: EK(:),TR(:),TI(:)

      End MODULE data_list

!======================================================================
      Subroutine alloc_data_list(m)
!======================================================================

      Use data_list

      Implicit none
      Integer(4) :: m,i

      if(m.le.0) then
       if(allocated(EK)) Deallocate(EK,TR,TI,jjpar,ip,l1,l2,jk1,jk2)
       mdata = 0; ndata = 0
      elseif(.not.allocated(EK)) then
       Allocate(EK(m),TR(m),TI(m),jjpar(m),ip(m),l1(m),l2(m),jk1(m),jk2(m))
       mdata = m; ndata = 0
      elseif(m.le.mdata) then
       Return
      elseif(ndata.gt.0) then
       Open(nus,form='UNFORMATTED',status='SCRATCH')
       rewind(nus)
       Do i=1,ndata
        write(nus) EK(i),TR(i),TI(i),jjpar(i),ip(i),l1(i),l2(i),jk1(i),jk2(i)
       End do
       Deallocate(EK,TR,TI,jjpar,ip,l1,l2,jk1,jk2)
       Allocate(EK(m),TR(m),TI(m),jjpar(m),ip(m),l1(m),l2(m),jk1(m),jk2(m))
       rewind(nus)
       Do i=1,ndata
        read(nus) EK(i),TR(i),TI(i),jjpar(i),ip(i),l1(i),l2(i),jk1(i),jk2(i)
       End do
       Close(nus)
      end if

      End Subroutine alloc_data_list



!=======================================================================
      Subroutine Iadd_data(e,ar,ai,jj,ii,m1,m2,k1,k2)
!=======================================================================
!     add the data to the list 
!-----------------------------------------------------------------------

      Use data_list

      Implicit none
      Integer(4) :: jj,ii,m1,m2,k1,k2, i
      Real(8) :: e,ar,ai

      if(mdata.eq.0) Call Alloc_data_list(idata)

! ... check if the same data are already in the list:

      Do i=1,ndata
       if(jj.ne.jjpar(i)) Cycle
       if(ii.ne.ip(i)) Cycle
       if(m1.ne.l1(i)) Cycle
       if(m2.ne.l2(i)) Cycle
       if(k1.ne.jk1(i)) Cycle
       if(k2.ne.jk2(i)) Cycle
       if(e.ne.EK(i)) Cycle              !  ???
       TR(i)=TR(i)+ar
       TI(i)=TI(i)+ai
       Return
      End do

! ... add new integral:

      if(ndata.eq.mdata) Call Alloc_data_list(mdata+idata)

      ndata = ndata+1 
      jjpar(ndata) = jj
      ip(ndata) = ii
      l1(ndata) = m1
      l2(ndata) = m2
      jk1(ndata) = k1
      jk2(ndata) = k2
      EK(ndata) = e
      TR(ndata) = ar
      TI(ndata) = ai

      End Subroutine Iadd_data

!=====================================================================
!     UTILITY  get1tdonejj
!=====================================================================
!     get 'tmat.done' file for 1 energy from 'zarm.tma'
!     with jj-coupling
!---------------------------------------------------------------------
      Use target_jj; Use channels_jj; Use data_list

      Implicit real(8) (A-H,O-Z)
      Real(8), allocatable :: matr(:),mati(:) 
      Real(8) ::  eps_e = 1.d-7

      Call Inf_screen
!----------------------------------------------------------------------
! ... target and channels information:

      nut=1; Open(nut,file='target_jj',status='OLD')
      Call Read_target_jj (nut)
      Call Read_channels_jj(nut)
      Close(nut)

      ion = nz-nelc
      zion=1.d0; if(ion.gt.1) zion=ion*ion

! ... define transition and energy under consideration:

      itr1=1; Call Read_iarg('itr1' ,itr1 )
      itr2=1; Call Read_iarg('itr2' ,itr2 )
      Call Read_rarg('ek1'  ,ek1  )
      Call Read_rarg('eps_e',eps_e)

!----------------------------------------------------------------------
! ... define input file and parameters:

      in = 1; Open(in,file='zarm.tma',status='OLD')

      E1=(ETARG(itr1)-ETARG(1))*2; E2=(ETARG(itr2)-ETARG(1))*2;
      JT1 = jtarg(itr1); JT2 = jtarg(itr2)
      ip1=(-ptarg(itr1)+1)/2; ip2=(-ptarg(itr2)+1)/2                    !  ???

! ... allocations:

      mtr = mch*(mch+1)/2;  Allocate(matr(mtr),mati(mtr))

!----------------------------------------------------------------------
! ... define T-matrix elements:

      rewind(in)
    1 read(in,*,end=2) ee,nopen,ntr,ii
      read(in,*) (matr(i),mati(i),i=1,ntr)

      if(abs(ee-ek1).gt.eps_e) go to 1

      jj = jpar(ii); ipp = (-ipar(ii)+1)/2

      Do i = 1,nch(ii); if(iptar(ii,i).ne.itr2) Cycle
      Do j = 1,i;       if(iptar(ii,j).ne.itr1) Cycle

        ktr = (i-1)*i/2+j; if(ktr.gt.ntr) Stop 'ktr > ntr'

        lch1=lch(ii,j); lch2=lch(ii,i); ll1=lch1+lch1; ll2=lch2+lch2  
        jch1=jch(ii,j); jch2=jch(ii,i); jj1=jch1; jj2=jch2

        Call Fano_tma(lch1,lch2,matr(ktr),mati(ktr),br,bi)

!br=matr(ktr); bi=mati(ktr)

        k1_min = iabs(JT1-ll1); k1_max = JT1+ll1
        Do kk1 = k1_min,k1_max
         S1 = T_jj_jk(ll1,1,JT1,jj1,kk1,JJ)
         if(S1.eq.0.d0) Cycle

        k2_min = iabs(JT2-ll2); k2_max = JT2+ll2
        Do kk2 = k2_min,k2_max
         S2 = T_jj_jk(ll2,1,JT2,jj2,kk2,JJ)
         if(S2.eq.0.d0) Cycle

         S = S1*S2;   ar = br * S; ai = bi * S

         m=1
         if(itr1.eq.itr2) then
           if(lch1.eq.lch2.and.kk1.gt.kk2) m=0
           if(lch1.gt.lch2) m=0
         end if

!         if(m.ne.0) 
        
          Call Iadd_data(ee,ar,ai,jj,ipp,lch1,lch2,kk1,kk2)


       End do; End do


       if(i.eq.j.or.itr1.ne.itr2) Cycle


        lch1=lch(ii,i); lch2=lch(ii,j); ll1=lch1+lch1; ll2=lch2+lch2  
        jch1=jch(ii,i); jch2=jch(ii,j); jj1=jch1; jj2=jch2

!        Call Fano_tma(lch1,lch2,matr(ktr),mati(ktr),br,bi)

        k1_min = iabs(JT1-ll1); k1_max = JT1+ll1
        Do kk1 = k1_min,k1_max
         S1 = T_jj_jk(ll1,1,JT1,jj1,kk1,JJ)
         if(S1.eq.0.d0) Cycle

        k2_min = iabs(JT2-ll2); k2_max = JT2+ll2
        Do kk2 = k2_min,k2_max
         S2 = T_jj_jk(ll2,1,JT2,jj2,kk2,JJ)
         if(S2.eq.0.d0) Cycle

         S = S1*S2;   ar = br * S; ai = bi * S

         m=1
         if(itr1.eq.itr2) then
           if(lch1.eq.lch2.and.kk1.gt.kk2) m=0
           if(lch1.gt.lch2) m=0
         end if

!         if(m.ne.0) 
          Call Iadd_data(ee,ar,ai,jj,ipp,lch1,lch2,kk1,kk2)


        End do; End do


      End do; End do


      go to 1
    2 Close(in)

!----------------------------------------------------------------------
! ... output:

      iout = 3; Open(iout,file='tmat.done')
      ne = 1
      write(iout,'(8I5,2D15.7)') itr1,itr2,JT1,JT2,ip1,ip2,ndata,ne,E1,E2

      Do i=1,ndata
       write(iout,'(D15.7,6I5,2D15.7)') EK(i),jjpar(i),ip(i), &
                    l1(i),jk1(i),l2(i),jk2(i),TR(i),TI(i)
      End do

      Close(iout)

      End  ! UTILITY  get1tdonejj


!======================================================================
      Subroutine inf_screen
!======================================================================
!     provides screen information about the utility
!----------------------------------------------------------------------
      Character :: A

      Call get_command_argument(1,A)  
      if(A.ne.'?') Return

      write(*,'(a)') &
'                                                            ',&
' get1tdonejj gets "tmat.done" file for 1 energy from        ',&
' "zarm.tma" file in case of jj-coupling                     ',&
'                                                            ',&        
' Call as: get1tdonejj itr1=.. itr2=..  ek1=..  [eps_e=..]   ',&
'                                                            ',&        
' itr1  - index of initial state (default - 1)               ',&
' itr2  - index of final state  (default - 1)                ',&
' ek1   - electron energy in Ry                              ',&
' eps_e - tolerance for electron energy [1.d-7]              ',&
'                                                            '        
      Stop ' '

      End Subroutine inf_screen

