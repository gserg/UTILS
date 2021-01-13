        program momtrans
        implicit double precision (a-h,o-z)

        Ry = 13.6054d0
        pi = acos(-1.0d0)

        write(*,*) 'input E0 and E1 (in eV) and theta1 (in deg)'
        read(*,*) e0,e1,theta1

        xk0 = sqrt(e0/13.605d0)
        xk1 = sqrt(e1/13.605d0)
        xk1x = xk1*cos(theta1*pi/180.0d0)
        xk1y = xk1*sin(theta1*pi/180.0d0)
        qx = xk0-xk1x
        qy = -xk1y
        thetaq = atan2(qy,qx)
        eq = (qx**2+qy**2)*13.605d0

        write(*,'(a,f6.2,a,f6.2)') 'momentum transfer angle (deg):', &
                  thetaq*180.0/pi,'   EQ =',eq

        end ! program momtrans

