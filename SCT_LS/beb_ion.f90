!=================================================================================
!     Use BEB formula to estimate direct ionization from the given level
!
!     Call as:   beb_ion  AF    
!
!     where AF - file with input data (see beb_inp as example)
!     Results are added to the same file
!---------------------------------------------------------------------------------

      Real :: Ry=13.6056,  PI=3.1415926
      Character(80) :: AF
      Character(4), allocatable :: EL(:)
      Real, allocatable :: EE(:),B(:),U(:),N(:)

      Call get_command_argument(1,AF)
      nu = 1
      open(nu,file=AF,status='OLD')

      read(nu,*) 
      read(nu,*) ns
      read(nu,*) 
      Allocate(EL(ns), B(ns), U(ns), N(ns))
      Do i = 1,ns
       read(nu,*)  EL(i), B(I),U(i),N(i)
      End do
      read(nu,*) 
      read(nu,*)  E1, E2, de

!-----------------------------------------------------------------------------
      E = E1
      Do

       Stotal=0.0;  ER = E/Ry
       Do i = 1,ns
        if(ER.lt.B(i)) Cycle
        uu = U(i)/B(I)*2     !  ???
        t = ER/B(i)
        S =  4*PI * N(i) /B(i)**2 
        tl = log(t)
        sec = tl/2*(1-1./t/t) + 1 - 1/t - tl/(1+t)
        sec = sec / (t+uu+1)
        Stotal = Stotal + sec * S
       End do
       Stotal = Stotal * 0.28
       if(Stotal.ne.0.0) write(nu,'(f10.2,E15.5)')  E, Stotal

       E = E + de
       if(E.gt.E2) Exit
      End do

      END  ! program
                 
       
        
       




      







