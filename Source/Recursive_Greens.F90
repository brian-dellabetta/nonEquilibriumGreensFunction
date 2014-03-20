!Get G = (H)^-1
SUBROUTINE Get_Retarded_Greens
    USE NEGF_Module
    IMPLICIT NONE
    
    !Assumes off-diagonal blocks are always Hermitian (Al = conjg(transpose(Au))) and invGRu is the same regardless of z index, just use (-HopX) instead
    COMPLEX(8), DIMENSION(:,:,:), ALLOCATABLE :: gLR, gRR
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE   :: cmat1, cmat2

    !For simplicity, store matrix in arrays of upper-,on-,and lower-diagonal components
    ALLOCATE(cmat1(NNz,NNz))    ;cmat1 = c0
    ALLOCATE(cmat2(NNz,NNz))    ;cmat2 = c0

    !Get gLR
    ALLOCATE(gLR(NNz,NNz,Npx))  ;gLR = c0
    DO ii = 1, Npx
        IF (ii>1) THEN
            CALL zgemm('C', 'N', NNz, NNz, NNz, c1, -HopX, NNz, gLR(:,:,ii-1), NNz, c0, cmat1, NNz)
            CALL zgemm('N', 'N', NNz, NNz, NNz, c1, cmat1, NNz, -HopX, NNz, c0, cmat2, NNz)
            cmat1 = InvGRd(:,:,ii) - cmat2
        ELSE
            cmat1 = InvGRd(:,:,ii)
        ENDIF

        !invert cmat1, store in gLR
        CALL Invert_Matrix(gLR(:,:,ii), cmat1, NNz, errorFlag)
        
        IF (errorFlag /= 0) THEN
            WRITE(*,*) "ERROR: Inverse of gLR is singular... Dumping data"
            WRITE(*,*) "INFO = ", errorFlag
            WRITE(*,*) "Energy = ", energy
            errorFlag=0
        ENDIF
    ENDDO
    
    !Backward Recursion to get GRd,GRu,GRl
    GRd(:,:,Npx) = gLR(:,:,Npx)
    DO ii = (Npx-1), 1, -1
        !GRl(:,:,ii)=-GRd(:,:,ii+1)*Al(:,:,ii)*gLR(:,:,ii)
        CALL zgemm('N', 'C', NNz, NNz, NNz, -c1, GRd(:,:,ii+1), NNz, -HopX, NNz, c0, cmat1, NNz)
        CALL zgemm('N', 'N', NNz, NNz, NNz, c1, cmat1, NNz, gLR(:,:,ii), NNz, c0, GRl(:,:,ii), NNz)
        !GRu(:,:,ii)=-gLR(:,:,ii)*-HopX(:,:,ii)*GRd(:,:,ii+1)
        CALL zgemm('N', 'N', NNz, NNz, NNz, -c1, gLR(:,:,ii), NNz, -HopX, NNz, c0, cmat1, NNz)
        CALL zgemm('N', 'N', NNz, NNz, NNz, c1, cmat1, NNz, GRd(:,:,ii+1), NNz, c0, GRu(:,:,ii), NNz)
        !GRd(:,:,ii)=gLR(:,:,ii)* ( eye(NNz)--HopX(:,:,ii)*GRl(:,:,ii) )
        CALL zgemm('N', 'N', NNz, NNz, NNz, -c1, -HopX, NNz, GRl(:,:,ii), NNz, c0, cmat1, NNz)
        DO jj = 1, NNz
            cmat1(jj,jj) = cmat1(jj,jj) + c1
        ENDDO
        CALL zgemm('N', 'N', NNz, NNz, NNz, c1, gLR(:,:,ii), NNz, cmat1, NNz, c0, GRd(:,:,ii), NNz)
    ENDDO
    
    GRRgt = c0
    GRRgt(1+(Npx-1)*NNz:Npx*NNz, :) = GRd(:,:,Npx)
    GRRgt(1+(Npx-2)*NNz:(Npx-1)*NNz, :) = GRu(:,:,Npx-1)

    !Expand to off diagonal parts GRRgt
    DO ii = 1, Npx-2
        jj = Npx-1-ii
        !GR(1+(jj-1)*NNz:jj*NNz,1+(jj+ii)*NNz:(jj+ii+1)*NNz) = -gLR(:,:,jj)*-HopX(:,:,jj)*GR(1+jj*NNz:(jj+1)*NNz,1+(jj+ii)*NNz:(jj+ii+1)*NNz);
        CALL zgemm('N', 'N', NNz, NNz, NNz, -c1, gLR(:,:,jj), NNz, -HopX(:,:), NNz, c0, cmat1, NNz)
        CALL zgemm('N', 'N', NNz, NNz, NNz, c1, cmat1, NNz, GRRgt(1+jj*NNz:(jj+1)*NNz, :), NNz, c0, GRRgt(1+(jj-1)*NNz:jj*NNz, :), NNz)
    ENDDO
    DEALLOCATE(gLR)
    
    !gRR needed for GRLef, GRTop
    ALLOCATE(gRR(NNz,NNz,Npx)); gRR = c0
    DO ii = Npx, 1, -1
        IF (ii<Npx) THEN
            CALL zgemm('N', 'N', NNz, NNz, NNz, c1, -HopX(:,:), NNz, gRR(:,:,ii+1), NNz, c0, cmat1, NNz)
            CALL zgemm('N', 'C', NNz, NNz, NNz, c1, cmat1, NNz, -HopX(:,:), NNz, c0, cmat2, NNz)
            cmat1 = InvGRd(:,:,ii) - cmat2
        ELSE
            cmat1 = InvGRd(:,:,ii)
        ENDIF
        
        !invert cmat1, store in gRR
        CALL Invert_Matrix(gRR(:,:,ii), cmat1, NNz, errorFlag)
        
        IF (errorFlag /= 0) THEN
            WRITE(*,*) "ERROR: Inverse of gRR is singular... Dumping data"
            WRITE(*,*) "INFO = ", errorFlag
            WRITE(*,*) "Energy = ", energy
            errorFlag=0
        ENDIF
    ENDDO

    GRLef = c0; 
    GRLef(1:NNz, :) = GRd(:,:,1)
    GRLef(1+NNz:2*NNz, :) = GRl(:,:,1)

    DO ii = 1, Npx-1
        !GR(1+(jj-1)*NNz:jj*NNz,1+(jj+ii)*NNz:(jj+ii+1)*NNz) = -gRR(:,:,jj)*-HopX(:,:,jj)*GR(1+jj*NNz:(jj+1)*NNz,1+(jj+ii)*NNz:(jj+ii+1)*NNz);
        CALL zgemm('N', 'C', NNz, NNz, NNz, -c1, gRR(:,:,ii+1), NNz, -HopX(:,:), NNz, c0, cmat1, NNz)
        CALL zgemm('N', 'N', NNz, NNz, NNz, c1, cmat1, NNz, GRLef(1+(ii-1)*NNz:ii*NNz, :), NNz, c0, GRLef(1+ii*NNz:(ii+1)*NNz, :), NNz)
    ENDDO
    DEALLOCATE(gRR)
    
    DEALLOCATE(cmat1)
    DEALLOCATE(cmat2)

    IF (errorFlag /= 0) THEN
        WRITE(*,*) "ERROR: Inverse of GR is singular... Dumping data"
        WRITE(*,*) "INFO = ", errorFlag
        WRITE(*,*) "Energy = ", energy
    ENDIF
END SUBROUTINE Get_Retarded_Greens