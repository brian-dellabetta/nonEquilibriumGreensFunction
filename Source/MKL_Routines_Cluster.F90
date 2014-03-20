!Get invG = (G)^-1
SUBROUTINE Invert_Matrix(invG, G, N, errorFlag)
    USE Constants
    IMPLICIT NONE
    
    !INCLUDE 'mkl.h' ! mkl_solver_ilp64_sequential.lib mkl_intel_ilp64.lib mkl_sequential.lib  mkl_core.lib'
    
    INTEGER(4)                            :: N, errorFlag  !invG, G must be NxN matrices
    COMPLEX(8), DIMENSION(:,:)            :: invG(N,N), G(N,N)
    COMPLEX(8), DIMENSION(:), ALLOCATABLE :: WORK
    INTEGER(4), DIMENSION(:), ALLOCATABLE :: IPIV
    INTEGER(4)                            :: LWORK = 64*8, INFO

    ALLOCATE(WORK(LWORK))               ;WORK = 0D0
    ALLOCATE(IPIV(N))                   ;IPIV = 0

    invG = G
    !Factor invG
    CALL ZGETRF(N, N, invG, N, IPIV, INFO)  !CALL ZHETRF('U', N, invG, N, IPIV, WORK, LWORK, INFO)
    
    IF (INFO == 0) THEN
        !Compute inverse of invG
        CALL ZGETRI(N, invG, N, IPIV, WORK, N, INFO)    !CALL ZHETRI('U', N, invG, N, IPIV, WORK, INFO)
        errorFlag = INFO
    ELSE
        WRITE(*,*) "Error in ZGETRF function:  INFO = ", INFO
        errorFlag = INFO
    ENDIF
    
    DEALLOCATE(WORK)
    DEALLOCATE(IPIV)
    
END SUBROUTINE invert_Matrix

!solve A*psi = eigval*psi
SUBROUTINE Get_Eigvals(Nx, Nz, Htotd, HopZ, eigvals)
    USE Constants
    IMPLICIT NONE
    INTEGER(4)  :: Nx, Nz, errorFlag, ii
    COMPLEX(8)  :: Htotd(Nx,Nx,Nz), HopZ(Nx,Nx)
    REAL(8)     :: eigvals(Nx*Nz)
    INTEGER(4)  :: N
    COMPLEX(16) :: WORK(Nx*Nz*(Nx*Nz+3))
    REAL(8)     :: RWORK(Nx*Nz*(2*Nx*Nz+6))
    INTEGER     :: IWORK(5*Nx*Nz+5)
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: HDMat
   
    errorFlag=0;
    N=Nx*Nz
    ALLOCATE(HDMat(N,N)); HDMat = c0
    
    DO ii = 1, Nz
        HDMat(1+(ii-1)*Nx:ii*Nx, 1+(ii-1)*Nx:ii*Nx) = Htotd(:,:,ii)
        
        IF (ii<Nz) THEN
            HDMat(1+(ii-1)*Nx:ii*Nx, 1+ii*Nx:(ii+1)*Nx) = HopZ
        
            HDMat(1+ii*Nx:(ii+1)*Nx, 1+(ii-1)*Nx:ii*Nx) = conjg(transpose(HopZ))
        ENDIF
    ENDDO

    CALL ZHEEVD('N', 'U', N, HDMat, N, eigvals, WORK, N*(N+3), RWORK, N*(2*N+6), IWORK, 5*N+5, errorFlag)
    IF (errorFlag /= 0) THEN
        WRITE(*,*) "Error in get_Eigvals computation, error in parameter", errorFlag
    ENDIF
END SUBROUTINE Get_Eigvals


!solve A*psi = eigval*B*psi
SUBROUTINE Get_CEigvals(N, A, B, eigval)
    USE Constants
    IMPLICIT NONE
    INTEGER(4)  :: N, errorFlag
    COMPLEX(8)  :: A(N,N), B(N,N), eigval(N), tau(N-1)
    COMPLEX(8)  :: C(N,N), WORK(4*N)
    errorFlag = 0
    !inv(B)*A*psi = eigval*psi
    
    CALL invert_Matrix(B, B, N, errorFlag)
    CALL zgemm('N', 'N', N, N, N, c1, B, N, A, N, c0, C, N)
    
    CALL ZGEHRD(N, 1, N, C, N, tau, WORK, 4*N, errorFlag)
    IF (errorFlag /= 0) THEN
        WRITE(*,*) "Error in Hess transformation, error in parameter", errorFlag
    ENDIF
    
    CALL ZHSEQR('E', 'N', N, 1, N, C, N, eigval, A*c0, N, WORK, 4*N, errorFlag)
    IF (errorFlag /= 0) THEN
        WRITE(*,*) "Error in get_Eigvals computation, error in parameter", errorFlag
    ENDIF
END SUBROUTINE Get_CEigvals


!Get invG = (G)^-1 for tridiagonal matrix G
SUBROUTINE Invert_Tridiag(invG, G, N, errorFlag)
    USE Constants
    IMPLICIT NONE
    
    !Recursive Green's functions for a tridiagonal N*N matrix G
    INTEGER(4)                            :: N
    COMPLEX(8)                            :: invG(N,N), G(N,N)
    INTEGER(4)                            :: errorFlag
    COMPLEX(8), DIMENSION(:), ALLOCATABLE :: gLR
    INTEGER(4)                            :: ii, jj
    
    ALLOCATE(gLR(1:N))          ;gLR = c0
    
    DO ii = 1,N
        DO jj = 1,N
            IF ((abs(ii-jj) >= 2) .AND. (G(ii,jj) /= 0D0)) THEN
                errorFlag = jj +N*(ii-1)
                RETURN
            ENDIF
        ENDDO
    ENDDO
    
    gLR(1) = c1/G(1,1)
    !forward recursion
    DO ii = 2, N-1
        gLR(ii) = c1/(G(ii,ii) - G(ii,ii-1)*gLR(ii-1)*G(ii-1,ii))
    ENDDO
    gLR(N) = c1/(G(N,N) - G(N,N-1)*gLR(N-1)*G(N-1,N))
    
    !backward recursion
    invG(N,N) = gLR(N)
    DO ii = N-1, 1, -1
        invG(ii,ii) = gLR(ii)*(1+G(ii,ii+1)*invG(ii+1,ii+1)*G(ii+1,ii)*gLR(ii))
    ENDDO
    
    !off diagonal coordinates
    DO ii = 0, N
        DO jj = N-1,1,-1
            IF (jj+1+ii <= N) THEN
                invG(jj,jj+1+ii) = -1D0*gLR(jj)*G(jj,jj+1)*invG(jj+1,jj+1+ii)
                invG(jj+1+ii,jj) = invG(jj,jj+1+ii)
            ENDIF
        ENDDO
    ENDDO
    
    !check if NaN's exist
    DO ii = 1, N
        DO jj = 1, N
            IF (isnan(abs(invG(ii,jj)))) THEN
                WRITE(*,*) "Error in invert_tridiag function..."
                WRITE(*,*) "ii, jj, invG(ii, jj) is", ii, jj, invG(ii,jj)
                errorFlag=1
            ENDIF
        ENDDO
    ENDDO

    DEALLOCATE(gLR)
    
END SUBROUTINE Invert_Tridiag