! Basic Lorenz96 implementation
MODULE utils
  use, intrinsic ::  iso_fortran_env, dp=>real64

  IMPLICIT NONE
  !INTEGER, PARAMETER :: dp=real64!kind(0.d0)
  
  CONTAINS
  
  
  
  
  !! Periodic Boundary conditions:
  !! copy the elements from the right of state to the left padding, and vice-versa
  SUBROUTINE ringbc(state_a, buff_size, ndim)
    REAL(dp)           :: state_a(:,:)
    INTEGER            :: buff_size, ndim
  
    !! dummy variables:
    INTEGER            :: i
  
    !! Left BCs:
    DO i = 1,buff_size
      state_a(i,:) = state_a(ndim+i,:)
    END DO
  
    !! Right BCs:
    DO i = (ndim+1+buff_size), (ndim+2*buff_size)
      state_a(i,:) = state_a(i-ndim, :)
    END DO
  
  END SUBROUTINE ringbc

  
  ! compute dx/dt for lorenz63 model based on state_a, a matrix of ensemble states
  ! the returned value is in deltax
  SUBROUTINE l63diff(state_a, params, deltax)
    REAL(dp)           :: state_a(:,:)
    REAL(dp)           :: params(3)
    REAL(dp)           :: deltax(:,:)
    
    deltax(1,:) = params(1) * (state_a(2,:) - state_a(1,:))
    deltax(2,:) = state_a(1,:) * (params(2) - state_a(3,:)) - state_a(2,:)
    deltax(3,:) = state_a(1,:) * state_a(2,:) - params(3) * state_a(3,:)
  
  
  END SUBROUTINE l63diff
  
  
  ! compute dx/dt for lorenz96 model based on state_a, a matrix of ensemble states
  ! the returned value is in deltax
  SUBROUTINE l96diff(state_a, forcing, ndim, ind, deltax)
    REAL(dp)           :: state_a(:,:), forcing
    REAL(dp)           :: deltax(:,:)
    INTEGER            :: ndim
    INTEGER            :: ind(ndim)
    INTEGER            :: i
    
    DO i = 1,ndim
      deltax(i,:) = - state_a(ind(i)-2,:) * state_a(ind(i)-1,:) + &
         & state_a(ind(i)-1,:) * state_a(ind(i)+1,:) - &
         & state_a(ind(i),:) + forcing
    END DO
  
  END SUBROUTINE l96diff
  
  
  ! same as l96diff but forcing is a vector
  SUBROUTINE l96diffvec(state_a, forcingvec, ndim, ind, deltax)
    REAL(dp)           :: state_a(:,:), forcingvec(:)
    REAL(dp)           :: deltax(:,:)
    INTEGER            :: ndim
    INTEGER            :: ind(ndim)
    INTEGER            :: i
    
    DO i = 1,ndim
      deltax(i,:) = - state_a(ind(i)-2,:) * state_a(ind(i)-1,:) + &
         & state_a(ind(i)-1,:) * state_a(ind(i)+1,:) - &
         & state_a(ind(i),:) + forcingvec(i)
    END DO
  END SUBROUTINE l96diffvec
  
  ! compute the increment deltax of the lorenz96 model using
  ! runge-kutta. At the end one can do state_a + deltat*deltax
  ! it depends on l96diff
  SUBROUTINE l96rk(state_a, forcing, ndim, kens, ind, buff_size, deltat, deltax)
    REAL(dp)                        :: state_a(:,:), forcing
    REAL(dp), ALLOCATABLE           :: state_x(:,:)
    REAL(dp)                        :: deltax(:,:)
    REAL(dp)                        :: deltat
    INTEGER                         :: ndim,kens, buff_size
    INTEGER                         :: ind(ndim)
    REAL(dp), DIMENSION(ndim,kens)  :: k1, k2, k3, k4
    
    ALLOCATE(state_x(size(state_a,1), size(state_a,2)))
    
    !k1:
    CALL l96diff(state_a, forcing, ndim, ind, k1)
    
    !k2:
    state_x(ind,:) = state_a(ind,:) + deltat/2._dp * k1
    CALL ringbc( state_x, buff_size, ndim)
    CALL l96diff(state_x, forcing, ndim, ind, k2)
    
    !k3:
    state_x(ind,:) = state_a(ind,:) + deltat/2._dp * k2
    CALL ringbc( state_x, buff_size, ndim)
    CALL l96diff(state_x, forcing, ndim, ind, k3)
    
    !k4:
    state_x(ind,:) = state_a(ind,:) + deltat * k3
    CALL ringbc( state_x, buff_size, ndim)
    CALL l96diff(state_x, forcing, ndim, ind, k4)
    
    !all together:
    deltax = 1._dp/6._dp * ( k1 + 2*k2 + 2*k3 + k4 )
    
  END SUBROUTINE l96rk
    







  ! compute dx/dt for lorenz MII model based on state_a, a matrix of ensemble states
  ! the returned value is in deltax
  SUBROUTINE lm2diff(state_a, kint, forcing, ndim, kens, ind, deltax, buff_size)
    REAL(dp)           :: state_a(:,:), forcing
    INTEGER            :: kint
    REAL(dp)           :: deltax(:,:)
    INTEGER            :: buff_size
    INTEGER            :: ndim, kens
    INTEGER            :: ind(ndim)
    INTEGER            :: nj  ! used in summation for [X,X]_Kn
    INTEGER            :: i
    
    REAL(dp), ALLOCATABLE :: w(:,:)
    
    !! x is buffered on both side.
    !! dx needs w_n-2K and w_n-2K needs x_n-2K-J= X_n-5J
    !! but it is easier to have w the same size as x
    !! so:
    !! x = [ 5nj ndim 5nj] = 2*5*nj + ndim = 2*buff_size + ndim
    !! w = [ 5nj ndim 5nj] = 2*5*nj + ndim = 2*buff_size + ndim
    !! dx =[     ndim    ] = ndim

    nj = kint/2
    
    
    ALLOCATE( w( 2*buff_size + ndim, kens ) )
    
    
    !! initialize w1:
    w(ind(1),:) = 1._dp/kint * (                             &
                       & state_a( ind(1)+nj,:)/2._dp  +      &
                       & state_a( ind(1)-nj,:)/2._dp  +      &
                  & sum( state_a((ind(1)-nj+1):(ind(1)+nj-1),:), 1 ) )
    
    !! build w incrementally:
    DO i=1,ndim-1 ! size(w,1)-1
      w(ind(i+1),:) = w(ind(i),:) + 1._dp/kint/2._dp * (     &
                       & state_a(ind(i)+nj,:)   + &
                       & state_a(ind(i)+nj+1,:) - &
                       & state_a(ind(i)-nj,:)   - &
                       & state_a(ind(i)-nj+1,:) )
    END DO
    !! Repeat w periodically:
    CALL ringbc(w, buff_size, ndim)
    
    !! build deltax from w:
    DO i=1,ndim
      deltax(i,:) = - w(ind(i)-2*kint,:) * w(ind(i)-kint,:) + &
                    & w(ind(i)-kint-nj,:) * w(ind(i)+kint-nj,:)/2._dp/kint +&
                    & w(ind(i)-kint+nj,:) * w(ind(i)+kint+nj,:)/2._dp/kint + &
               & sum(w((ind(i)-kint-nj+1):(ind(i)-kint+nj-1),:)* &
                     w((ind(i)+kint-nj+1):(ind(i)+kint+nj-1),:), 1)/kint - &
                     state_a(ind(i),:) + forcing
    END DO

    
  
  END SUBROUTINE lm2diff



  ! compute the increment deltax of the lorenzM2 model using
  ! runge-kutta. At the end one can do state_a + deltat*deltax
  ! it depends on lm2diff
  SUBROUTINE lm2rk(state_a, kint, forcing, ndim, kens, ind, buff_size, deltat, deltax)
    REAL(dp)                        :: state_a(:,:), forcing
    INTEGER                         :: kint
    REAL(dp), ALLOCATABLE           :: state_x(:,:)
    REAL(dp)                        :: deltax(:,:)
    REAL(dp)                        :: deltat
    INTEGER                         :: ndim,kens, buff_size
    INTEGER                         :: ind(ndim)
    REAL(dp), DIMENSION(ndim,kens)  :: k1, k2, k3, k4
    
    ALLOCATE(state_x( ndim+2*buff_size, kens))
    
    !k1:
    CALL lm2diff(state_a, kint, forcing, ndim, kens, ind, k1, buff_size)
    
    !k2:
    state_x(ind,:) = state_a(ind,:) + deltat/2._dp * k1
    CALL ringbc( state_x, buff_size, ndim)
    CALL lm2diff(state_a, kint, forcing, ndim, kens, ind, k2, buff_size)
    
    !k3:
    state_x(ind,:) = state_a(ind,:) + deltat/2._dp * k2
    CALL ringbc( state_x, buff_size, ndim)
    CALL lm2diff(state_a, kint, forcing, ndim, kens, ind, k3, buff_size)
    
    !k4:
    state_x(ind,:) = state_a(ind,:) + deltat * k3
    CALL ringbc( state_x, buff_size, ndim)
    CALL lm2diff(state_a, kint, forcing, ndim, kens, ind, k4, buff_size)
    
    !all together:
    deltax = 1._dp/6._dp * ( k1 + 2*k2 + 2*k3 + k4 )
    
  END SUBROUTINE lm2rk

END MODULE utils

