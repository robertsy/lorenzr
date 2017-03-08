! Basic Lorenz63 implementation


SUBROUTINE l63integrate(state, kens, params, nsteps, deltat)
  USE utils

  IMPLICIT NONE

  INTEGER :: kens                                            !number of ensemble and 
  REAL(dp), DIMENSION(3*kens)                 :: state       ! state vector 
  REAL(dp), DIMENSION(3)                      :: params      ! forcing parameters (sigma, rho and beta)
  INTEGER                                     :: nsteps
  REAL(dp)                                    :: deltat      ! increments in time

  !! dummy variables:
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)       :: state_a     ! field array
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)       :: deltax      ! increments in x

  INTEGER                                     :: i,j
  REAL(dp)                                    :: sigma, rho, beta
  
  !--------------------------------------------------------------
  ! unpack parameters: 
  sigma = params(1)
  rho = params(2)
  beta = params(3)
  
  !sigma <- 10; rho <- 28; beta <- 8/3
  
  !--------------------------------------------------------------
  ! allocate state array 
  ALLOCATE( state_a( 3, kens ) )
  ALLOCATE( deltax ( 3, kens ) )
  
  ! stack vectors in array:
  DO j = 1,kens
      state_a(:,j) = state( ((j-1)*3 + 1) : j*3)
  END DO

  !--------------------------------------------------------------
  ! time integration  
  DO i=1,nsteps
    CALL l63diff(state_a, params, deltax) 
    state_a(:,:) = state_a(:,:) + deltat * deltax 
  END DO

  !--------------------------------------------------------------
  ! put back the array in big vectors:
  DO j = 1,kens
      state( ((j-1)*3 + 1) : j*3) = state_a(:,j)
  END DO
    
  

END SUBROUTINE l63integrate


