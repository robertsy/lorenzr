! Basic Lorenz96 implementation


SUBROUTINE l96integrate(state, kens, ndim, forcing, nsteps, deltat, forcingvec)
  USE utils

  IMPLICIT NONE

  INTEGER :: kens, ndim              !number of ensemble and size of the domain (in # of cells)

  REAL(dp), DIMENSION(ndim*kens)              :: state       ! field
  REAL(dp)                                    :: forcing     ! forcing parameter (typically 8)
  INTEGER                                     :: nsteps
  REAL(dp)                                    :: deltat      ! increments in time
  REAL(dp), DIMENSION(ndim)                   :: forcingvec  ! 

  !! dummy variables:
  INTEGER                                     :: buff_size=2 ! size of buffer for BC
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)       :: state_a     ! field array
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)       :: deltax      ! increments in x
  INTEGER                                     :: i,j
  INTEGER                                     :: ind(ndim)   ! indices of original grid points

  
  
  !--------------------------------------------------------------
  ! allocate state array with buffer on both sides for periodic BCs:
  !ALLOCATE( state_a( -(buff_size-1):(ndim+buff_size), kens ) )
  ALLOCATE( state_a( 2*buff_size+ndim, kens ) )

  
  ! indices of grid points (not including the BCs)
  ind = (/(i, i=buff_size+1, buff_size+ndim)/)
  
  ! stack vectors in array:
  DO j = 1,kens
      state_a(ind,j) = state( ((j-1)*ndim + 1) : j*ndim)
  END DO
    
  CALL ringbc( state_a, buff_size, ndim)
  
  
  !--------------------------------------------------------------
  ! allocate the array of increments deltax:
  ALLOCATE( deltax (ndim, kens ) )

  
  DO i=1,nsteps
    !CALL l96rk(state_a, forcing, ndim, kens, ind, buff_size, deltat, deltax)
    !CALL l96diff(state_a, forcing, ndim, ind, deltax)
    CALL l96diffvec(state_a, forcingvec, ndim, ind, deltax)
    state_a(ind,:) = state_a(ind,:) + deltat * deltax 
    CALL ringbc( state_a, buff_size, ndim)
  END DO

  !--------------------------------------------------------------
  ! put back the array in big vectors:
  DO j = 1,kens
      state( ((j-1)*ndim + 1) : j*ndim) = state_a(ind,j)
  END DO
    
  

END SUBROUTINE l96integrate


