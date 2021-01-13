!======================================================================
      Subroutine Compare(pri,n1,C1,n2,C2)
!======================================================================

      Implicit none 
      Integer :: pri,n1,n2,i,j,n
      Integer :: ipt1(n1), ipt2(n2)
      Real(8) :: C1(n1), C2(n2)
      Character(200) :: A1(5), A2(5), AS

      A1 = ' '
      Call SORTA(n1,C1,IPT1) 
      n = min(5,n1)
      Do j = 1,n; i = ipt1(j)
       Call Get_cfg_jj(i)
       Call Label_jj (200,AS,0)
       write(A1(j),'(F8.5,2x,a)') C1(i),trim(AS) 
      End do

      A2 = ' '
      Call SORTA(n2,C2,IPT2) 
      n = min(5,n2)
      Do j = 1,n; i = ipt2(j)
       Call Get_cfg_LS(i)
       Call Label_c (AS,1,0)
       write(A2(j),'(F8.5,2x,a)') C2(i),trim(AS) 
      End do

      Do i = 1,5
       write(pri,'(a,T55,a)') trim(A1(i)), trim(A2(i))
	 End do 

      End Subroutine Compare 

