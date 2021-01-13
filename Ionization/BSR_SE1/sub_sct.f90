!======================================================================
      Subroutine Sub_sct(EK) 
!======================================================================
! ... scattering calculation for partial wave "nnn" and energy EK (k^2);
! ... results: R-matrix, K-matrix and boundary values (F,G,FP,GP)   
!----------------------------------------------------------------------
      Use bsr_se
      Use h_nnn

      Implicit none
      Integer :: i,info
      Real(8) :: EK,dmn,dmx,dav
      Real(8), external :: DET

! ... read h.nnn:

      Call Read_h_nnn(AF_h)

      ! kch, khm, klsp - in module bsr_sev
      ! nch, nhm       - in module h_nnn

      if(kch.ne.nch) Stop 'different nch in h.dat'
      if(khm.ne.nhm) Stop 'different nhm in h.dat'

! ... allocate arrays:

      if(allocated(RMAT)) Deallocate(RMAT,KMAT,ECH,F,G,FP,GP)

      Allocate(RMAT(nch,nch),KMAT(nch,nch),ECH(nch), &
               F(nch,nch),G(nch,nch),FP(nch,nch),GP(nch,nch))

! ... define open channels:
 
      Call ZOPEN (EK, ntarg, etarg, NCONAT, nch, nopen, ECH)

      if(nopen.le.0) Stop 'Stop in Sub_sct: number of open channels <= 0'

! ... define asymptotic:

      Call ZAFACE(debug,AC,ION,km,RA,delta,nch,nopen,LCH,ECH,CF, &
                  iauto,mfgi,info,F,G,FP,GP)                                  

      if(info.ne.0) Stop 'Stop in Sub_sct: info in ASYMPT /= 0'

! ... define R-matrix:
 
      Call ZRMAT (EK,RA,nch,nhm,RMAT,VALUE,WMAT)

      if(debug.gt.0) then
       write(pri,'(/a,f10.5/)') 'RMAT matrix:  EK = ',EK
       Do i=1,nch; write(pri,'(5E15.5)') RMAT(i,1:nch); End do
       write(pri,'(/a,E15.5/)') 'RM determinant:',DET(nch,RMAT)
       Call ZRMAT (EK,RA,nch,khm,RMAT,VALUE,WMAT)
      end if

! ... define K-matrix:
 
      Call ZKMAT(nch,nopen,RA,RB,RMAT,F,G,FP,GP,KMAT) 

      if(debug.gt.0) then
       write(pri,'(/a/)') 'K-matrix:'
       Do i=1,nch; write(pri,'(5E15.5)') KMAT(i,1:nch); End do
      end if

! ... evaluate the symmetry of the K-matrix:

      Call Sym_mat (nopen,nch,KMAT,dmn,dmx,dav)
 
      if(debug.gt.0) write(pri,'(/a,3f7.1/)') &
          ' K-matrix asymmetry, LOG(Min/Max/Av): ',dmn,dmx,dav

      END Subroutine Sub_sct 
 
