!======================================================================
      Subroutine rovl_out (AF_mat,AF_rsol,AF_rovl)
!======================================================================
!     R-matrix overlaps vectors:     a = S b,     where
!
!     S - overlap matrix, read from AF_mat file [bsr_mat.nnn] 
!
!     b - R-matrix solution, read from AF_rsol file [rsol.nnn] 
!
!     a - recoded in AF_rovl file (for each R-matrix solution)
!----------------------------------------------------------------------
      Implicit none

      Character(*) :: AF_mat,AF_rsol,AF_rovl
      Real(8), allocatable :: c(:,:),cc(:,:),v(:),w(:)
      Real(8) :: S                                                                                                                                     !
      Integer :: i,j,i1,i2,j1,j2, ic,jc, is,ich,ii, nui,nur,nuo,&
                 ns,kch,kcp,nhm,khm, ns1,kch1,kcp1,nhm1

!----------------------------------------------------------------------
! ... read overlap matrix in full:

      Call Check_file(AF_mat )
      Call Find_free_unit(nui)
      Open(nui,file=AF_mat ,form='UNFORMATTED')
 
      rewind(nui)
      read(nui) ns,kch,kcp;   nhm = ns*kch + kcp

      Allocate(cc(ns,ns),c(nhm,nhm),v(ns));  c=0.d0; ii=ns*kch-kch

      ! ... diagonal blocks:    
 
      Do ich=1,kch
       read(nui) cc
       i1=(ich-1)*ns+1; i2=ich*ns; c(i1:i2,i1:i2) = cc
      End do 

      ! ... other blocks:

      Do 
       read(nui) ic,jc;  if(ic.le.0) Exit

       if(ic.gt.kch.and.jc.gt.kch) then  !  pert-pert

        read(nui) S;  c(ic+ii,jc+ii) = S; c(jc+ii,ic+ii) = S 
        
       elseif(ic.gt.kch) then            !  ch-pert

        read(nui) v
        ic = ic + ii;  j1=(jc-1)*ns+1; j2=jc*ns
        C(ic,j1:j2)=v; C(j1:j2,ic) = v

       else                              !  ch-ch

        read(nui) cc
        i1=(ic-1)*ns+1; i2=ic*ns
        j1=(jc-1)*ns+1; j2=jc*ns
        C(i1:i2,j1:j2)=cc; C(j1:j2,i1:i2) = cc

       end if

      End do

      Close(nui); Deallocate(cc,v)

!----------------------------------------------------------------------
! ... check rsol file and multiply overlap matrix on solutions:
        
      Call Check_file(AF_rsol)
      Call Find_free_unit(nur) 
      Open(nur,file=AF_rsol,form='UNFORMATTED')

      Call Find_free_unit(nuo) 
      Open(nuo,file=AF_rovl,form='UNFORMATTED')

      rewind(nur)
      read(nur) nhm1,khm,kch1,kcp1,ns1

      if(ns1 .ne.ns ) Stop ' BSR_OVL: different ns  in rsol file'
      if(kch1.ne.kch) Stop ' BSR_OVL: different kch in rsol file'
      if(kcp1.ne.kcp) Stop ' BSR_OVL: different kcp in rsol file'

      read(nur) (S,i=1,khm)  ! skip energies ?
             
      Allocate(v(1:nhm),w(1:nhm))

      write(nuo) khm,nhm,kch,kcp,ns

      Do is=1,khm
       read(nur) v;  w=MATMUL(c,v); write(nuo) w
      End do
      
      Close(nur); Close(nuo) 

      Deallocate(c,v,w)

      End Subroutine rovl_out
