!======================================================================
!     Utility:   delete_data
!======================================================================
!
!     deletes  data for given energies from K-matrix, T-matrix or
!     Omega files, if any
!
!     Call:   delete_kma  del=  type=  inp=  out=   ek1=  ek2=  ekk=
!
!     del [delete_list] - contains list of deleting energies in Ry (as a column)
!
!     type [all]        - type of data: kma,tma,tmb,om,omb,par,all      
!                                                                       
!     inp - input file for given type if different from standard name   
!     out - output file for given type if different from standard name  
!                                                                       
!     ek1 - minimum electron energy to keep                             
!     ek2 - maximum electron energy to keep                             
!                                                                       
!     ekk - the only energy to keep                                     
!                                                                       
!     This utility is used to clean the data files, e.g., K-matrix file 
!     if we have  incomplete data for some energies                      
!                                                                       
!     The information about incomplete data may be obtained with        
!     utility kma_om                                                    
!----------------------------------------------------------------------

      Implicit real(8) (a-h,o-z)

      Character(80) :: AF, type, inp, out
      Real(8), allocatable :: matr(:), ed(:)

! ... file units and standard names:

      Integer :: nu1 = 11   ! input data 
      Integer :: nu2 = 12   ! output data
      Integer :: nud = 13   ! delete list

      Character(80) :: AF_del ='delete_list'

      Character(80) :: AF_kma1='zarm.kma'
      Character(80) :: AF_kma2='zarm.kma_'

      Character(80) :: AF_tma1='zarm.tma'
      Character(80) :: AF_tma2='zarm.tma_'

      Character(80) :: AF_om1 ='zarm.om'
      Character(80) :: AF_om2 ='zarm.om_'

      Character(80) :: AF_tmb1='zarm.tmb'
      Character(80) :: AF_tmb2='zarm.tmb_'

      Character(80) :: AF_omb1='zarm.omb'
      Character(80) :: AF_omb2='zarm.omb_'

      Character(80) :: AF_par1='zarm.omb_par'
      Character(80) :: AF_par2='zarm.omb_par_'

      Call inf_del 

!-----------------------------------------------------------------------
! ... get arguments:

      Call Read_aarg('del',AF_del)

      type = 'all'; Call Read_aarg('type',type)

      inp = ' ';  Call Read_aarg('inp',inp)
      out = ' ';  Call Read_aarg('out',out)

      ek1=0.d0; Call Read_rarg('ek1',ek1)
      ek2=0.d0; Call Read_rarg('ek2',ek2)
      if(ek2.lt.ek1) ek2 = ek1

      ekk=0.d0; Call Read_rarg('ekk',ekk)

! ... initial estimation:

      mdim = 10000; Allocate(matr(mdim))

!-----------------------------------------------------------------------
! ... read deleted energies:

      ne = 0
      if(Icheck_file(AF_del).gt.0) then

      Open(nud,file=AF_del)
    1 read(nud,*,end=2) E; ne=ne+1; go to 1
    2 write(*,*) 'ne =',ne

      Allocate(ED(ne))
      rewind(nud)
      Do i=1,ne; read(nud,*) ED(i); End do   

      end if

!-----------------------------------------------------------------------
! ... select data and clean:

      Select case(type) 

       Case('kma');  Call Clean_kma
       Case('tma');  Call Clean_tma
       Case('tmb');  Call Clean_tmb
       Case('om ');  Call Clean_om
       Case('omb');  Call Clean_omb
       Case('par');  Call Clean_par

       Case('all');  Call Clean_kma
                     Call Clean_tma
                     Call Clean_tma
                     Call Clean_om 
                     Call Clean_omb
                     Call Clean_par
      End Select


 CONTAINS


!======================================================================
      Subroutine Clean_kma
!======================================================================

      if(len_trim(inp).eq.0) inp = AF_kma1
      if(len_trim(out).eq.0) out = AF_kma2

      if(Icheck_file(inp).eq.0) Return

      Open(nu1,file=inp)
      Open(nu2,file=out)

   11 read(nu1,*,end=12) e,nopen,ntr,ilsp
      if(ntr.gt.mdim) then
       Deallocate(matr); mdim=ntr; Allocate(matr(mdim))
      end if
      read(nu1,*) (matr(i),i=1,ntr)

      if(ek1.gt.0.d0.and.e.lt.ek1) go to 11
      if(ek2.gt.0.d0.and.e.gt.ek2) go to 11

      if(ekk.ne.0.d0.and.e.ne.ekk) go to 11

      if(ne.gt.0) then
       m=0; Do i=1,ne; if(e.ne.ed(i)) Cycle; m=1; Exit; End do
       if(m.eq.1) go to 11
      end if

      write(nu2,'(F10.6,3i8)') e,nopen,ntr,ilsp
      write(nu2,'(5D16.8)') (matr(i),i=1,ntr)

      go to 11
   12 Close(nu1); Close(nu2)

      End Subroutine Clean_kma


!======================================================================
      Subroutine Clean_tma
!======================================================================

      if(len_trim(inp).eq.0) inp = AF_tma1
      if(len_trim(out).eq.0) out = AF_tma2

      if(Icheck_file(inp).eq.0) Return

      Open(nu1,file=inp)
      Open(nu2,file=out)

   11 read(nu1,*,end=12) e,nopen,ktr,ilsp
      ntr=2*ktr
      if(ntr.gt.mdim) then
       Deallocate(matr); mdim=ntr; Allocate(matr(mdim))
      end if
      read(nu1,*) (matr(i),i=1,ntr)

      if(ek1.gt.0.d0.and.e.lt.ek1) go to 11
      if(ek2.gt.0.d0.and.e.gt.ek2) go to 11

      if(ekk.ne.0.d0.and.e.ne.ekk) go to 11

      if(ne.gt.0) then
       m=0; Do i=1,ne; if(e.ne.ed(i)) Cycle; m=1; Exit; End do
       if(m.eq.1) go to 11
      end if

      write(nu2,'(F10.6,3i8)') e,nopen,ktr,ilsp
      write(nu2,'(5D16.8)') (matr(i),i=1,ntr)

      go to 11
   12 Close(nu1); Close(nu2)

      End Subroutine Clean_tma


!======================================================================
      Subroutine Clean_tmb
!======================================================================

      if(len_trim(inp).eq.0) inp = AF_tmb1
      if(len_trim(out).eq.0) out = AF_tmb2

      if(Icheck_file(inp).eq.0) Return

      Open(nu1,file=inp)
      Open(nu2,file=out)

   11 read(nu1,*,end=12) e,nopen,kopen,ilsp,i1,i2,nj
      ntr = kopen*(kopen+1)/2
      ktr = (nopen-kopen)*nj

      kdim = 2 * (ntr+ktr)
      if(kdim.gt.mdim) then
       Deallocate(matr); mdim=kdim; Allocate(matr(mdim))
      end if
      read(nu1,*) (matr(i),i=1,2*ntr)
      if(ktr.gt.0) read(nu1,*) (matr(i),i=2*ntr+1,kdim)

      if(ek1.gt.0.d0.and.e.lt.ek1) go to 11
      if(ek2.gt.0.d0.and.e.gt.ek2) go to 11

      if(ekk.ne.0.d0.and.e.ne.ekk) go to 11

      if(ne.gt.0) then
       m=0; Do i=1,ne; if(e.ne.ed(i)) Cycle; m=1; Exit; End do
       if(m.eq.1) go to 11
      end if

      write(nu2,'(F10.6,6i6,a)') e,nopen,kopen,ilsp,i1,i2,nj ,'   ee,nopen,kp,ilsp,np,ni,nj'
      write(nu2,'(6D16.8)')  (matr(i),i=1,2*ntr)
      if(ktr.gt.0) write(nu2,'(6D16.8)') (matr(i),i=2*ntr+1,kdim)

      go to 11
   12 Close(nu1); Close(nu2)

      End Subroutine Clean_tmb


!======================================================================
      Subroutine Clean_om
!======================================================================

      if(len_trim(inp).eq.0) inp = AF_om1
      if(len_trim(out).eq.0) out = AF_om2

      if(Icheck_file(inp).eq.0) Return

      Open(nu1,file=inp)
      Open(nu2,file=out)

   11 read(nu1,*,end=12) e,ntr
      if(ntr.gt.mdim) then
       Deallocate(matr); mdim=ntr; Allocate(matr(ntr))
      end if
      read(nu1,*) (matr(i),i=1,ntr)

      if(ek1.gt.0.d0.and.e.lt.ek1) go to 11
      if(ek2.gt.0.d0.and.e.gt.ek2) go to 11

      if(ekk.ne.0.d0.and.e.ne.ekk) go to 11

      if(ne.gt.0) then
       m=0; Do i=1,ne; if(e.ne.ed(i)) Cycle; m=1; Exit; End do
       if(m.eq.1) go to 11
      end if

      write(nu2,'(F10.6,i8)') e,ntr
      write(nu2,'(5D16.8)') (matr(i),i=1,ntr)

      go to 11
   12 Close(nu1); Close(nu2)

      End Subroutine Clean_om


!======================================================================
      Subroutine Clean_omb
!======================================================================

      if(len_trim(inp).eq.0) inp = AF_omb1
      if(len_trim(out).eq.0) out = AF_omb2

      if(Icheck_file(inp).eq.0) Return

      Open(nu1,file=inp)
      Open(nu2,file=out)

   11 read(nu1,*,end=12) e,ntr,k1,k2,i1,i2
      if(ntr.gt.mdim) then
       Deallocate(matr); mdim=ntr; Allocate(matr(ntr))
      end if
      read(nu1,*) (matr(i),i=1,ntr)

      if(ek1.gt.0.d0.and.e.lt.ek1) go to 11
      if(ek2.gt.0.d0.and.e.gt.ek2) go to 11

      if(ekk.ne.0.d0.and.e.ne.ekk) go to 11

      if(ne.gt.0) then
       m=0; Do i=1,ne; if(e.ne.ed(i)) Cycle; m=1; Exit; End do
       if(m.eq.1) go to 11
      end if

      write(nu2,'(F10.6,5i8)') e,ntr,k1,k2,i1,i2
      write(nu2,'(5D16.8)') (matr(i),i=1,ntr)

      go to 11
   12 Close(nu1); Close(nu2)

      End Subroutine Clean_omb


!======================================================================
      Subroutine Clean_par
!======================================================================

      if(len_trim(inp).eq.0) inp = AF_par1
      if(len_trim(out).eq.0) out = AF_par2

      if(Icheck_file(inp).eq.0) Return

      Open(nu1,file=inp)
      Open(nu2,file=out)

   11 read(nu1,*,end=12) e,ntr,k1,k2,i1,i2
      if(ntr.gt.mdim) then
       Deallocate(matr); mdim=ntr; Allocate(matr(ntr))
      end if
      read(nu1,*) (matr(i),i=1,ntr)

      if(ek1.gt.0.d0.and.e.lt.ek1) go to 11
      if(ek2.gt.0.d0.and.e.gt.ek2) go to 11

      if(ekk.ne.0.d0.and.e.ne.ekk) go to 11

      if(ne.gt.0) then
       m=0; Do i=1,ne; if(e.ne.ed(i)) Cycle; m=1; Exit; End do
       if(m.eq.1) go to 11
      end if

      write(nu2,'(F10.6,5i8)') e,ntr,k1,k2,i1,i2
      write(nu2,'(5D16.8)') (matr(i),i=1,ntr)

      go to 11
   12 Close(nu1); Close(nu2)

      End Subroutine Clean_par


      End  ! utility delete_kma


!======================================================================
      Subroutine inf_del
!======================================================================
!     provide screen information about add_farm utility
!----------------------------------------------------------------------
       
      Character(80) :: A

      iarg = command_argument_count()
      if(iarg.eq.0) Return
      Call GET_COMMAND_ARGUMENT(1,A)        
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                       ',&
'     deletes all data for given energies from K-matrix, T-matrix or    ',&
'     Omega files, if any                                               ',&
'                                                                       ',&
'     Call:   delete_data  [del=..  type=..  inp=..  out=..             ',&
'                           ek1=..  ek2=..  ekk=..]                     ',&
'                                                                       ',&
'     del [delete_list] - contains list of deleting energies in Ry      ',&
'                         (as a column)                                 ',&
'                                                                       ',&
'     type [all]        - type of data: kma,tma,tmb,om,omb,par,all      ',&
'                                                                       ',&
'     inp - input file for given type if different from standard name   ',&
'     out - output file for given type if different from standard name  ',&
'                                                                       ',&
'     ek1 - minimum electron energy to keep                             ',&
'     ek2 - maximum electron energy to keep                             ',&
'                                                                       ',&
'     ekk - the only energy to keep                                     ',&
'                                                                       ',&
'     This utility is used to clean the data files, e.g., K-matrix file ',&
'     if we have  incomplete data for some energies                     ',&  
'                                                                       ',&
'     The information about incomplete data may be obtained with        ',&
'     utility kma_om                                                    ',&
'                                                                       '

      Stop ' '

      End Subroutine inf_del

