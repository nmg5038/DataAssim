%%file ASSIMILATION.f95

MODULE ASSIMILATION
USE LORENZ_63
CONTAINS


      SUBROUTINE RUN_FORWARD(X,Nsteps,Tstep,N,Xnew)
      ! Steps a number of states forward using the Lorenz 63 equations
      INTEGER N,Nsteps
      REAL(kind=8) :: Tstep
      REAL(kind=8), DIMENSION(N) :: X
      REAL(kind=8), DIMENSION(Nsteps,N) :: Xnew
      INTEGER :: i
      REAL(kind=8), DIMENSION(N) :: Xhold
      REAL(kind=8), DIMENSION(N) :: dXdt
!f2py intent(in) X,Nsteps, Tstep
!f2py integer intent(hide),depend(X) :: n=shape(x,0)
!f2py intent(out) Xnew
      Xhold = X
      Xnew(1,:) = X
      do i = 2, Nsteps
      
           ! Run the non-nonlinear model
           CALL LORENZ63(Xhold,N,dXdt)
           ! Euler Forward time-differencing scheme executed over Nsteps
           Xhold = Xhold + dXdt * tstep 
           Xnew(i,:) = Xhold
      end do 
      
      END SUBROUTINE RUN_FORWARD




      SUBROUTINE CALC(state,dstate,innov,o_freq,first,last,inv_R,&
    &inv_B,Tstep,Nsteps,N,Nobs,cost)
      ! Calculate the gradient and full cost function
      INTEGER N,Nsteps,Nobs,o_freq,last,first,first_ind
      REAL(kind=8) :: Tstep
      REAL(kind=8) :: t_co
      
      REAL(kind=8), DIMENSION(Nobs,N) :: innov
      REAL(kind=8), DIMENSION(N,N) :: inv_R,inv_B
      REAL(kind=8), DIMENSION(N) :: dstate, adj_inno
      REAL(kind=8), DIMENSION(N) :: J_obs, lambda
      REAL(kind=8), DIMENSION(N) :: term1
      REAL(kind=8), DIMENSION(N+1) :: cost
      REAL(kind=8), DIMENSION(N,N) :: H
      REAL(kind=8), DIMENSION(Nsteps,N) :: state, state_tlm
      INTEGER :: i
      
!f2py intent(in) state, Tstep, Nsteps, innov
!f2py intent(in) inv_R,inv_B,o_freq,last,first
!f2py integer intent(hide),depend(state) :: nsteps=shape(state,0)
!f2py integer intent(hide),depend(state) :: n=shape(state,1)
!f2py integer intent(hide),depend(innov) :: nobs=shape(innov,0)
!f2py intent(out) cost
       
       ! If first is = 0, then no observation at initial time
       ! 		   is = 1, then observation at initial time
       ! Run the Tangent Linear Model Forward over the window
       CALL RUN_FORWARD_TLM(state,dstate,Tstep,Nsteps,N,&
    & state_tlm)

       ! Make an identity matrix    
       H= 0D0
       do i = 1,N,1
          H(i,i) = 1D0
       enddo
       t_co = 0D0
       
       last = last + 1

        ! If first is = 0, then no observation at initial time
        ! 		   is = 1, then observation at initial time
       if (first .eq. 1) then
           first_ind = 1
           adj_inno =  innov(1+last/o_freq,:)-matmul(H,state_tlm(last,:))
       else
           first_ind = 2
           adj_inno =  innov(last/o_freq,:)-matmul(H,state_tlm(last,:))
       endif

       ! Generate the contribution to the gradient of the cost function
       !  from the innovation and observational uncertainty
       J_obs = matmul(H, matmul(inv_R,adj_inno))
       
       ! Generate the contribution to the cost function from the 
       !  innovation and observational uncertainty
       t_co = dot_product(adj_inno, matmul(inv_R,adj_inno))
       
       ! Define the first lambda (at the last observation
       lambda = J_obs
       
       ! Step backwards in time toward the start of the window
       do i = last-o_freq, first_ind, -o_freq
            
            if (first .eq. 1) then
            adj_inno =  innov(1+i/o_freq,:)-matmul(H,state_tlm(i,:))
            else     
            adj_inno =  innov(i/o_freq,:)-matmul(H,state_tlm(i,:))
            endif
            ! Generate the contribution to the gradient of the cost function
            J_obs = matmul(H, matmul(inv_R,adj_inno))
            
            ! Generate the contribution to the cost function at this observation
             t_co = t_co+dot_product(adj_inno, matmul(inv_R,adj_inno))
            
            ! Run the adjoint model backward between observation times
            CALL RUN_BACKWARD_ADJ(state(i:i+o_freq,:),&
     &lambda,Tstep,o_freq+1,N,term1)
            
            ! Generate a new lambda
            lambda = J_obs + term1
       end do
       
       ! Run the adjoint model backward between first obs time and first timestep
       !     IFF no observation at first time
       if (first .eq. 0) then
           CALL RUN_BACKWARD_ADJ(state(1:1+o_freq,:),&
     &lambda,Tstep,o_freq+1,N,term1)
           lambda = term1
       endif

       ! Generate the contribution to the gradient of the cost function from
       !    backward steps and initial uncertainty/state increment
       term1=matmul(inv_B,dstate)-lambda
       do i = 1,N,1
            cost(i) = term1(i)
       end do
       ! Finalize cost function value
       t_co = t_co+dot_product(dstate, matmul(inv_B,dstate))
       cost(N+1) = t_co / 2.0D0
      
      END SUBROUTINE CALC

      SUBROUTINE RUN_FORWARD_TLM(state,dstate,Tstep,Nsteps,N,state_tlm)
      ! Steps a number of states forward using the Lorenz 63 Tangent Linear Model
      INTEGER N,Nsteps
      REAL(kind=8) :: Tstep
      REAL(kind=8), DIMENSION(N) :: dstate, dstate_dt
      REAL(kind=8), DIMENSION(Nsteps,N) :: state, state_tlm
      INTEGER :: i
      
!f2py intent(in) state, Tstep, Nsteps
!f2py integer intent(hide),depend(state) :: nsteps=shape(state,0)
!f2py integer intent(hide),depend(state) :: n=shape(state,1)
!f2py intent(out) state_tlm

      state_tlm(1,:) = dstate
      do i = 1, Nsteps-1
           ! Call the tangent linear model 
           CALL LORENZ63_JACOBIAN(state(i,:),state_tlm(i,:),&
    & N,dstate_dt)
           ! Deploy a euler forward time difference to the tlm
           state_tlm(i+1,:) = state_tlm(i,:)+dstate_dt * tstep
      end do 
      
      END SUBROUTINE RUN_FORWARD_TLM


      SUBROUTINE RUN_BACKWARD_ADJ(state,dstate,Tstep,Nsteps,N,dstatea_out)
      ! Steps a number of states backward using the Lorenz 63 equation adjoint
      INTEGER N,Nsteps
      REAL(kind=8) :: Tstep
      REAL(kind=8), DIMENSION(N) :: dstate, dstatea_out
      REAL(kind=8), DIMENSION(N) :: dstate_dt_a
      REAL(kind=8), DIMENSION(Nsteps,N) :: state
      INTEGER :: i
      
!f2py intent(in) state, Tstep, Nsteps
!f2py integer intent(hide),depend(state) :: nsteps=shape(state,0)
!f2py integer intent(hide),depend(state) :: n=shape(state,1)
!f2py intent(out) dstatea_out
        
        ! First step
        dstatea_out = dstate
        
        do i = Nsteps-1, 1, -1
           ! Call the adjoint equations
           CALL LORENZ63_JACOBIAN_TRANSPOSE(state(i,:),&
    & dstatea_out,N,dstate_dt_a)
            ! Deploy a euler forward time difference to the adjoint
            dstatea_out = dstatea_out + dstate_dt_a * Tstep
        end do
      
      END SUBROUTINE RUN_BACKWARD_ADJ
      

END MODULE