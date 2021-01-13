      Implicit real(8) (A-H,O-Z)

      Integer, parameter :: ns = 340, nt = 26
      Real(8) :: upsilon(nt,ns,ns), A_value(ns,ns), E(ns), SJ(ns), T(nt)
      Integer :: LI(ns), SI(ns)
      Character(18) :: AC, BC, conf(ns)
      Character(1) :: AL
      Character(2) :: type

      Z = 26.0
      IZ = Z + 0.1
      IZ0 = 1
      IZ1 = 2
      eion = 130655.4

      AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)      
      ev_cm = au_cm/au_eV
      write(*,*) 'ev_cm =', ev_cm
      write(*,*) 'au_cm =', au_cm

!---------------------------------------------------------------------
      open(1,file='Table_E')

      Do j=1,ns

       read(1,'(i3,3x,a18,i1,a1,2x,i2,2x,f10.5)') i, AC, IS, AL,JJ,EE

       if(len_trim(AC).eq.0) AC=BC
       conf(i) = AC
       BC = AC

       E(i) = EE * ev_cm

       if(len_trim(AL).gt.0)  IL = LA(AL)
       if(len_trim(AL).eq.0)  IL = JL
       JL = IL

       SJ(i) = JJ/2.0
       LI(i) = IL

       if(IS.le.0) IS = JS
       JS = IS
       SI(i) = IS

      End do 
       
      close(1)
!----------------------------------------------------------------------
      open(1,file='target')

      Do i=1,ns
       read(1,*) E(i)
      End do
      close(1)
!----------------------------------------------------------------------


      open(2,file='FeI_adf04')

      write(2,'(1a3,i2,2i10,f10.1)') 'Fe+',IZ0,IZ,IZ1,eion
      Do i = 1,ns
       write(2,'(i5,1x,a18,a1,i1,a1,i1,a1,f3.1,a1,f10.1)') &
           i,conf(i),'(',SI(i),')',LI(i),'(',SJ(i),')',E(i)
      End do

      write(2,'(i5)') -1

!----------------------------------------------------------------------
      open(1,file='Table_f')

      A_values = 0.d0
      Do 
       read(1,*,end=10) i,j, type, S1,S2,S3, A
       A_value(i,j) = A_value(i,j) + A
      End do

   10 close(1)
!---------------------------------------------------------------------
      open(1,file='temperatures')

      Do i=1,nt
       read(1,*) T(i)
      End do

      close(1)
!---------------------------------------------------------------------
      open(1,file='Table_upsilon')

      Do i=1,ns-1
       Do j=i+1,ns
        read(1,*) ii,jj, upsilon(:,i,j)
       End do
      End do

      close(1)
!----------------------------------------------------------------------

      zeff = iz1
      write(2,'(f5.1,i5,8x,1P20e9.2)') zeff, 3, (T(i),i=5,nt-4,2) 
      Do j = 2,ns
       Do i = 1,j-1
        write(2,'(2i5,1P20e9.2)') j, i, A_value(i,j), (upsilon(k,i,j),k=5,nt-4,2)
       End do
      End do
      write(2,'(i5)') -1
      write(2,'(2i5)') -1,-1

      End ! program







