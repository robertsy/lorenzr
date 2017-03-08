! Basic Lorenz96 implementation


SUBROUTINE lm2integrate(state, kens, ndim, kint, forcing, nsteps, deltat)
  USE utils

  IMPLICIT NONE

  INTEGER                                     :: kens, ndim  !number of ensemble and size of the domain (in # of cells)

  REAL(dp), DIMENSION(ndim*kens)              :: state       ! field
  INTEGER                                     :: kint        ! k, an even number, of [X,X]_Kn
  REAL(dp)                                    :: forcing     ! forcing parameter (typically 15)
  INTEGER                                     :: nsteps      ! number of time steps to integrate
  REAL(dp)                                    :: deltat      ! increments in time

  !! dummy variables:
  INTEGER                                     :: nj          ! used in summation for [X,X]_Kn
  INTEGER                                     :: buff_size   ! size of buffer for BC
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)       :: state_a     ! field array
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)       :: deltax      ! increments in x
  INTEGER                                     :: i,j
  INTEGER                                     :: ind(ndim)   ! indices of original grid points


  nj = kint/2
  buff_size = 5*nj
  
  !--------------------------------------------------------------
  ! allocate state array with buffer on both sides for periodic BCs:
  !ALLOCATE( state_a( -(buff_size-1):(ndim+buff_size), kens ) )
  ALLOCATE( state_a( 2*buff_size + ndim, kens ) )

  
  ! indices of grid points (not including the BCs)
  ind = (/(i, i=buff_size + 1, buff_size + ndim)/)
  
  ! stack vectors in array:
  DO j = 1,kens
      state_a(ind,j) = state( ((j-1)*ndim + 1) : j*ndim)
  END DO
    
  CALL ringbc( state_a, buff_size, ndim)
  
  
  !--------------------------------------------------------------
  ! allocate the array of increments deltax:
  ALLOCATE( deltax (ndim, kens ) )
  !CALL lm2rk(   state_a, kint, forcing, ndim, kens, ind, buff_size, deltat, deltax)
!  CALL lm2diff( state_a, kint, forcing, ndim, kens, ind, deltax, buff_size)
!  state_a(ind,:) = deltax
!  
  DO i=1,nsteps
    CALL lm2rk( state_a, kint, forcing, ndim, kens, ind, buff_size, deltat, deltax)
!    CALL lm2diff( state_a, kint, forcing, ndim, kens, ind, deltax)
    state_a(ind,:) = state_a(ind,:) + deltat * deltax 
    CALL ringbc( state_a, buff_size, ndim)
  END DO

  !--------------------------------------------------------------
  ! put back the array in big vectors:
  DO j = 1,kens
      state( ((j-1)*ndim + 1) : j*ndim) = state_a(ind,j)
  END DO
    
  

END SUBROUTINE lm2integrate


