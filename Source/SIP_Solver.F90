SUBROUTINE Solve_Equilibrium
    USE SIP_Module
    IMPLICIT NONE

    IF(poissonIsConverged == 0) THEN !Iterate until poisson solver has converged
        n = RhoN
        p = RhoP
        ro = n - p + doping
        
        PhiDelta=0D0
        
        !Set up contact boundary conditions
        Phi(-1:0,:,:) = mu(1)
        Phi(Npx-1:Npx,:,:) = mu(2)

        alpha = alphaArray(mod(ite, alphaArrayLength)+1)           !Update alpha
        
        !Solve poisson equation:
        CALL SIP_V1 !Run through SIP algorithm
    ENDIF       !end do loop, poisson equation solved

    IF (ite > 1) THEN
        WRITE(*,*) "Iteration =", ite
        WRITE(*,*) "Poisson Converged =", poissonIsConverged
    ENDIF

END SUBROUTINE Solve_Equilibrium


!     ======================================================================
!     ======================================================================
!                  IMPLEMENTATION OF THE STRONGLY IMPLICIT METHOD
!                    FOR THE SOLUTION OF THE 3D-POISSON EQUATION
!     ======================================================================

SUBROUTINE SIP_V1
    USE SIP_Module
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE  :: a, b, c, d, e, f, g
    
    REAL(8)           :: Evec_new
	
    ALLOCATE(R_vec((Npx-1),(Npy-1),(Npz-1)));              R_vec = 0D0
    ALLOCATE(V_vec(0:(Npx-1)+1,0:(Npy-1)+1,0:(Npz-1)+1));  V_vec = 0D0
    
    !Elements of Lower and Upper Triangular Matrices
    ALLOCATE(a(0:Npx,0:Npy,0:Npz));      a = 0D0
    ALLOCATE(b(0:Npx,0:Npy,0:Npz));      b = 0D0
    ALLOCATE(c(0:Npx,0:Npy,0:Npz));      c = 0D0
    ALLOCATE(d(0:Npx,0:Npy,0:Npz));      d = 0D0
    ALLOCATE(e(0:Npx,0:Npy,0:Npz));      e = 0D0
    ALLOCATE(f(0:Npx,0:Npy,0:Npz));      f = 0D0
    ALLOCATE(g(0:Npx,0:Npy,0:Npz));      g = 0D0
    
    !========================================
    !Solve Poisson Equation using SIP Method:
    !========================================
    !First Sweep - Elements of Upper and Lower Triangular Matrices:
    IF (mod(ite,2)==0) THEN
        DO kk = 1,(Npz-1)
            DO jj = 1,(Npy-1)
                DO ii = 1,(Npx-1)
                    IF((materialType(ii,jj,kk) == 1) .OR. (materialType(ii,jj,kk) == 2)) THEN   !BiSe Region
                        Evec_new = Evec(ii,jj,kk) - (ro(ii,jj,kk))
                    ELSE
                        Evec_new = Evec(ii,jj,kk)
                    ENDIF
                    a(ii,jj,kk) = Bvec(ii,jj,kk)/(1D0+alpha*(e(ii,jj,kk-1) + f(ii,jj,kk-1)))
                    IF (isnan(a(ii,jj,kk))) STOP "a is Nan"
                    
                    b(ii,jj,kk) = Cvec(ii,jj,kk)/(1D0+alpha*(e(ii,jj-1,kk) + g(ii,jj-1,kk)))
                    IF (isnan(b(ii,jj,kk))) STOP "b is Nan"
                    
                    c(ii,jj,kk) = Dvec(ii,jj,kk)/(1D0+alpha*(g(ii-1,jj,kk) + f(ii-1,jj,kk)))
                    IF (isnan(c(ii,jj,kk))) STOP "c is Nan"
                    
                    d(ii,jj,kk) = Evec_new - c(ii,jj,kk)*(e(ii-1,jj,kk) - alpha*g(ii-1,jj,kk) - alpha*f(ii-1,jj,kk)) &
                        - b(ii,jj,kk)*(f(ii,jj-1,kk) - alpha*e(ii,jj-1,kk) - alpha*g(ii,jj-1,kk))- a(ii,jj,kk)*(g(ii,jj,kk-1) &
                        - alpha*e(ii,jj,kk-1) - alpha*f(ii,jj,kk-1))
                    IF (isnan(d(ii,jj,kk))) STOP "d is Nan"
                    IF (d(ii,jj,kk) == 0) WRITE(*,*) "d(",ii,jj,kk,") is 0"
                    
                    e(ii,jj,kk) = (Fvec(ii,jj,kk) - alpha*(a(ii,jj,kk)*e(ii,jj,kk-1) + b(ii,jj,kk)*e(ii,jj-1,kk)))/d(ii,jj,kk)
                    IF (isnan(e(ii,jj,kk))) STOP "e is Nan"
                    
                    f(ii,jj,kk) = (Gvec(ii,jj,kk) - alpha*(c(ii,jj,kk)*f(ii-1,jj,kk) + a(ii,jj,kk)*f(ii,jj,kk-1)))/d(ii,jj,kk)
                    IF (isnan(f(ii,jj,kk))) STOP "f is Nan"
                    
                    g(ii,jj,kk) = (Hvec(ii,jj,kk) - alpha*(c(ii,jj,kk)*g(ii-1,jj,kk) + b(ii,jj,kk)*g(ii,jj-1,kk)))/d(ii,jj,kk)
                    IF (isnan(g(ii,jj,kk))) STOP "g is Nan"               
                ENDDO
            ENDDO
        ENDDO
        
        !Calculate Correction Term:
        DO kk = 1, (Npz-1)
            DO jj = 1,(Npy-1)
                DO ii = 1,(Npx-1)
                R_vec(ii,jj,kk) = ro(ii,jj,kk) - (&
                    Bvec(ii,jj,kk)*Phi(ii,jj,kk-1) + Cvec(ii,jj,kk)*Phi(ii,jj-1,kk) + Dvec(ii,jj,kk)*Phi(ii-1,jj,kk) + &
                    Evec(ii,jj,kk)*Phi(ii,jj,kk)   + Fvec(ii,jj,kk)*Phi(ii+1,jj,kk) + Gvec(ii,jj,kk)*Phi(ii,jj+1,kk) + Hvec(ii,jj,kk)*Phi(ii,jj,kk+1) )
                    
                V_vec(ii,jj,kk) = (R_vec(ii,jj,kk) - a(ii,jj,kk)*V_vec(ii,jj,kk-1) - b(ii,jj,kk)*V_vec(ii,jj-1,kk) &
                    - c(ii,jj,kk)*V_vec(ii-1,jj,kk))/d(ii,jj,kk)
                ENDDO
            ENDDO
        ENDDO 

        DO kk = (Npz-1), 1, -1
            DO jj = (Npy-1), 1, -1
                DO ii = (Npx-1), 1, -1
                    PhiDelta(ii,jj,kk) = V_vec(ii,jj,kk) - e(ii,jj,kk)*PhiDelta(ii+1,jj,kk) - &
                        f(ii,jj,kk)*PhiDelta(ii,jj+1,kk) - g(ii,jj,kk)*PhiDelta(ii,jj,kk+1)
                ENDDO
            ENDDO
        ENDDO  
    ENDIF
    !End for first-loop part
 
    !Second Sweep: 
    IF (mod(ite,2)==1) THEN 
        DO kk = (Npz-1), 1, -1
            DO jj = (Npy-1), 1, -1
                DO ii = 1, (Npx-1)
                    IF((materialType(ii,jj,kk) == 1) .OR. (materialType(ii,jj,kk) == 2)) THEN !GaAs Region
                        Evec_new = Evec(ii,jj,kk) - (n(ii,jj,kk) - p(ii,jj,kk))
                    ELSE
                        Evec_new = Evec(ii,jj,kk)
                    ENDIF
                    a(ii,jj,kk) = Hvec(ii,jj,kk)/(1D0+alpha*(e(ii,jj,kk+1)+f(ii,jj,kk+1)))
                    IF (isnan(a(ii,jj,kk))) STOP "a is Nan"
                    b(ii,jj,kk) = Gvec(ii,jj,kk)/(1D0+alpha*(e(ii,jj+1,kk)+g(ii,jj+1,kk)))
                    IF (isnan(b(ii,jj,kk))) STOP "b is Nan"
                    c(ii,jj,kk) = Dvec(ii,jj,kk)/(1D0+alpha*(g(ii-1,jj,kk)+f(ii-1,jj,kk)))
                    IF (isnan(c(ii,jj,kk))) STOP "c is Nan"
                    
                    d(ii,jj,kk) = Evec_new - a(ii,jj,kk)*(g(ii,jj,kk+1) - alpha*e(ii,jj,kk+1) - alpha*f(ii,jj,kk+1)) &
                        - b(ii,jj,kk)*(f(ii,jj+1,kk) - alpha*e(ii,jj+1,kk) - alpha*g(ii,jj+1,kk)) &
                        - c(ii,jj,kk)*(e(ii-1,jj,kk) - alpha*g(ii-1,jj,kk) - alpha*f(ii-1,jj,kk))
                    IF (isnan(d(ii,jj,kk))) THEN
                        WRITE(*,*) "ii,jj,kk are", ii,jj,kk
                        WRITE(*,*) "n, p are", n(ii,jj,kk),p(ii,jj,kk)
                        WRITE(*,*) "ro, Phi, mat type are ", ro(ii,jj,kk), Phi(ii,jj,kk), materialType(ii,jj,kk)
                        STOP "d is Nan"
                    ENDIF
                    
                    e(ii,jj,kk) = (Fvec(ii,jj,kk) - alpha*(a(ii,jj,kk)*e(ii,jj,kk+1) + b(ii,jj,kk)*e(ii,jj+1,kk)))/d(ii,jj,kk)
                    IF (isnan(e(ii,jj,kk))) STOP "e is Nan"
                    
                    f(ii,jj,kk) = (Cvec(ii,jj,kk) - alpha*(c(ii,jj,kk)*f(ii-1,jj,kk) + a(ii,jj,kk)*f(ii,jj,kk+1)))/d(ii,jj,kk)
                    IF (isnan(f(ii,jj,kk))) STOP "f is Nan"
                    
                    g(ii,jj,kk) = (Bvec(ii,jj,kk) - alpha*(c(ii,jj,kk)*g(ii-1,jj,kk) + b(ii,jj,kk)*g(ii,jj+1,kk)))/d(ii,jj,kk)
                    IF (isnan(g(ii,jj,kk))) STOP "g is Nan"
                ENDDO
            ENDDO
        ENDDO
                
        !Calculate Correction Term:   
        DO kk = (Npz-1), 1, -1
            DO jj = (Npy-1), 1, -1       
                DO ii = 1, (Npx-1)
                    R_vec(ii,jj,kk) = ro(ii,jj,kk) - (&
                         Bvec(ii,jj,kk)*Phi(ii,jj,kk-1) + Cvec(ii,jj,kk)*Phi(ii,jj-1,kk) + Dvec(ii,jj,kk)*Phi(ii-1,jj,kk) + &
                         Evec(ii,jj,kk)*Phi(ii,jj,kk)   + Fvec(ii,jj,kk)*Phi(ii+1,jj,kk) + Gvec(ii,jj,kk)*Phi(ii,jj+1,kk) + &
                         Hvec(ii,jj,kk)*Phi(ii,jj,kk+1) )
                    V_vec(ii,jj,kk) = (R_vec(ii,jj,kk) - a(ii,jj,kk)*V_vec(ii,jj,kk+1) - &
                         b(ii,jj,kk)*V_vec(ii,jj+1,kk) - c(ii,jj,kk)*V_vec(ii-1,jj,kk))/d(ii,jj,kk)
                ENDDO
            ENDDO
        ENDDO

        DO kk = 1, (Npz-1)
            DO jj = 1, (Npy-1)       
                DO ii = (Npx-1), 1, -1
                    PhiDelta(ii,jj,kk) = V_vec(ii,jj,kk) - e(ii,jj,kk)*PhiDelta(ii+1,jj,kk) - &
                        f(ii,jj,kk)*PhiDelta(ii,jj-1,kk) - g(ii,jj,kk)*PhiDelta(ii,jj,kk-1)
                ENDDO
            ENDDO
        ENDDO
    ENDIF
    !End of second-loop part
                
	!Check for convergence
	IF ( ite < 2 .OR. maxval(abs(PhiDelta)/(1D-10+abs(Phi))) > 1D-2 ) THEN
	    poissonIsConverged = 0
	ELSE
	    poissonIsConverged = 1
	ENDIF
	Phi = Phi + PhiDelta
	
	DEALLOCATE(a)
	DEALLOCATE(b)
	DEALLOCATE(c)
	DEALLOCATE(d)
	DEALLOCATE(e)
	DEALLOCATE(f)
	DEALLOCATE(g)
	DEALLOCATE(V_vec)
	DEALLOCATE(R_vec)
END SUBROUTINE SIP_V1