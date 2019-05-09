MODULE LORENZ_63

           REAL(kind=8), parameter :: sigma = 10.0, r = 28.0
           REAL(kind=8), parameter :: beta=8.0/3.0

CONTAINS

!=============================================================
      SUBROUTINE LORENZ63_JACOBIAN(state,dstate,N,dstate_dt)
      ! Generates the Tangent Linear Model for the Euler Forward Scheme
      INTEGER,intent(in):: N
      REAL(kind=8), DIMENSION(N),intent(in) :: state,dstate
      REAL(kind=8), DIMENSION(N),intent(out) :: dstate_dt
      REAL(kind=8) x, y,z, dx, dy, dz
      
      x = state(1)
      y = state(2)
      z = state(3)
      
      dx = dstate(1)
      dy = dstate(2)
      dz = dstate(3)
      
      dstate_dt = 0D0
      dstate_dt(1) =  sigma * (dy - dx) 
      dstate_dt(2) =  (r-z)*dx-dy-x*dz   
      dstate_dt(3) = y*dx+x*dy-beta*dz
      
      END SUBROUTINE LORENZ63_JACOBIAN
!=============================================================

!=============================================================
      SUBROUTINE LORENZ63_JACOBIAN_TRANSPOSE(state,dstate,N,dstate_dt)
      ! Generates the Tangent Linear Model for the Euler Forward Scheme
      INTEGER,intent(in):: N
      REAL(kind=8), DIMENSION(N),intent(in) :: state,dstate
      REAL(kind=8), DIMENSION(N),intent(out) :: dstate_dt
      REAL(kind=8) x, y,z, dx, dy, dz
      
      x = state(1)
      y = state(2)
      z = state(3)
      
      dx = dstate(1)
      dy = dstate(2)
      dz = dstate(3)
      
      dstate_dt = 0D0
      dstate_dt(1) = -sigma*dx + (r-z)*dy + y*dz
      dstate_dt(2) =  sigma*dx - dy + x*dz   
      dstate_dt(3) = -x*dy - beta*dz
      
      END SUBROUTINE LORENZ63_JACOBIAN_TRANSPOSE

!=============================================================

      SUBROUTINE LORENZ63(state,N,dstate_dt)
            ! Lorenz 63 nonlinear function
           INTEGER, intent(in) :: N
           REAL(kind=8), DIMENSION(N),intent(in) :: state
           REAL(kind=8), DIMENSION(N),intent(out) :: dstate_dt
           REAL(kind=8) x, y,z
           
           x = state(1)
           y = state(2)
           z = state(3)
           
           dstate_dt = 0D0
           dstate_dt(1) = sigma * ( y - x )
           dstate_dt(2) = x * (r - z) - y
           dstate_dt(3) = x * y - beta *  z
            
      END SUBROUTINE LORENZ63
END MODULE