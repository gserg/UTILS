!=======================================================================
      SUBROUTINE RM_OVLF (EK,nopen,d,dr,di,F,G,RMAT,KMAT) 
!=======================================================================
!     input:   d - R-matrix overlaps and scattering paramters (R,K,F,G)
!     results: dr,di (1:nopen)  - overlap m.e. for each open channel
!
!     calculations are done for given partial wave 'nnn'
!-----------------------------------------------------------------------
      Use h_nnn
      Use bsr_se, only: debug,pri

      Implicit none

      Integer, intent(in) :: nopen 
      Real(8), intent(in) :: EK, d(nhm),                  &
                             RMAT(nch,nch),KMAT(nch,nch), &
                             F(nch,nch),G(nch,nch) 

      Real(8), intent(out) :: dr(nch),di(nch)

      Real(8), parameter :: zero = 0.d0, one = 1.d0 

      Real(8) :: P,CK,s1,s2, dmn,dmx,dav,  &
                 FF(nch,nch), AA(nch,nch), BB(nch,nch)

      Real(8), external :: DET 
      Integer :: i,k,info

!-----------------------------------------------------------------------
!     DETERMINE THE ASYMPTOTIC SOLUTIONS WITH INGOING WAVE 
!     BOUNDARY CONDITIONS
!
!     F = S + C*K;    F- = -iF / (1 - iK); 
!
!     F-  -->  -i / sqrt(PI*k)  [ sin  + cos * K ] / (1+K^2) * (1+iK)  
!-----------------------------------------------------------------------

! ... define F function

      FF(1:nch,1:nopen) = MATMUL(G(1:nch,1:nch),KMAT(1:nch,1:nopen))
      FF(1:nch,1:nopen) = F(1:nch,1:nopen) + FF(1:nch,1:nopen)

      if(debug.gt.0) then
       write(pri,'(/a/)') 'FF-matrix:'
       Do i=1,nch; write(pri,'(5E15.5)') FF(i,1:nopen); End do
      end if
                       
! ... F- --> F / sqrt(pi)

      P = ACOS(-one)
      P = one / sqrt(P)
      FF(1:nch,1:nopen) = P * FF(1:nch,1:nopen)

!-----------------------------------------------------------------------
!     matrix elements for given solution can be obtained from the
!     matrix elements between initial state and R-matrix basis 
!     states (array DK) by weighting them with expansion coefficients 
!
!         A(k) = (1/2a) (E(k) - E) SUM(i) [ w(i,k) (R^-1) F(i) ]
!
!     for given solution j.
!-----------------------------------------------------------------------

! ... AA -->  R^(-1) *  F-

      BB=RMAT; 
 
      Call INV(nch,nch,BB)

      Call Sym_mat (nopen,nch,BB,dmn,dmx,dav)
 
      if(debug.gt.0) then
       write(pri,'(/a/)') 'RMAT-1 matrix:'
       Do i=1,nch; write(pri,'(5E15.5)') BB(i,1:nch); End do
      end if

      AA(1:nch,1:nopen) = MATMUL(BB(1:nch,1:nch),FF(1:nch,1:nopen))            

! ... overlap m.e.

      di = zero
      Do K = 1,NHM
       CK = D(k) / (VALUE(K)-EK)
       di(1:nch) = di(1:nch) + CK*WMAT(1:nch,k)
      End do
      di = di / RA 
      dr(1:nopen) = MATMUL(di(1:nch),AA(1:nch,1:nopen))
      di = zero

! ... define  1 / (1+K^2),  where K - open-open part

      AA(1:nopen,1:nopen) = MATMUL(KMAT(1:nopen,1:nopen),KMAT(1:nopen,1:nopen))
      Do i = 1,nopen;  AA(i,i) = AA(i,i) + one;  End do
      Call INV(nopen,nch,AA)
      Call Sym_mat (nopen,nch,AA,dmn,dmx,dav)
 
      dr(1:nopen) =  MATMUL(AA(1:nopen,1:nopen),dr(1:nopen))
      di(1:nopen) = -MATMUL(KMAT(1:nopen,1:nopen),dr(1:nopen))

      End Subroutine RM_OVLF
 
