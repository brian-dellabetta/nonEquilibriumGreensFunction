SUBROUTINE Calculate_Gnp
    USE NEGF_Module
    IMPLICIT NONE
    INTEGER :: index
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: cmat1
    
    ALLOCATE(cmat1(Ntot, NNz))
    
    !!!!!!!!!!!!!Calculate Gn/Gp/A!!!!!!!!!
    cmat1 = c0; Gnd = c0; Gnu = c0; Gnl = c0
    !cmat1 = GR*Gamma(c1,c2)*f(1,2)
    DO ii = 1, Ntot
        DO jj = 1, NNz
            DO kk = 1, NNz
                cmat1(ii,jj) = cmat1(ii,jj) + f1*GRLef(ii, kk)*GammaL(kk,jj)
            ENDDO
        ENDDO
    ENDDO
    !Gn = cmat1*GR'
    DO xx = 1, Npx
        DO ii = 1, NNz
            DO jj = 1, NNz
                DO kk = 1, NNz
                    !Diagonal Block Gnd(x- and y-direction)
                    Gnd(ii, jj, xx) = Gnd(ii, jj, xx) + cmat1(ii+(xx-1)*NNz, kk)*conjg(GRLef(jj+(xx-1)*NNz, kk))
                    !Nearest neighbor Gn (depth z-direction)
                    IF (xx<Npx) THEN
                        Gnu(ii, jj, xx) = Gnu(ii, jj, xx) + cmat1(ii+(xx-1)*NNz, kk)*conjg(GRLef(jj+xx*NNz, kk))

                        Gnl(ii, jj, xx) = Gnl(ii, jj, xx) + cmat1(ii+xx*NNz, kk)*conjg(GRLef(jj+(xx-1)*NNz, kk))
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    cmat1 = c0; 
    DO ii = 1, Ntot
        DO jj = 1, NNz
            DO kk = 1, NNz
                cmat1(ii,jj) = cmat1(ii,jj) + f2*GRRgt(ii, kk)*GammaR(kk,jj)
            ENDDO
        ENDDO
    ENDDO
    DO xx = 1, Npx
        DO ii = 1, NNz
            DO jj = 1, NNz
                DO kk = 1, NNz
                    !Diagonal Block Gnd(x- and y-direction)
                    Gnd(ii, jj, xx) = Gnd(ii, jj, xx) + cmat1(ii+(xx-1)*NNz, kk)*conjg(GRRgt(jj+(xx-1)*NNz, kk))
                    
                    !Nearest neighbor Gn (depth z-direction)
                    IF (xx<Npx) THEN
                        Gnu(ii, jj, xx) = Gnu(ii, jj, xx) + cmat1(ii+(xx-1)*NNz, kk)*conjg(GRRgt(jj+xx*NNz, kk))

                        Gnl(ii, jj, xx) = Gnl(ii, jj, xx) + cmat1(ii+xx*NNz, kk)*conjg(GRRgt(jj+(xx-1)*NNz, kk))
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO

    cmat1 = c0; Gpd = c0; Gpu = c0; Gpl = c0
    !cmat1 = GR*Gamma(c1,c2)*f(1,2)v
    DO ii = 1, Ntot
        DO jj = 1, NNz
            DO kk = 1, NNz
                cmat1(ii,jj) = cmat1(ii,jj) + f1v*GRLef(ii, kk)*GammaL(kk,jj)
            ENDDO
        ENDDO
    ENDDO
    !Gp = cmat1*GR'
     DO xx = 1, Npx
        DO ii = 1, NNz
            DO jj = 1, NNz
                DO kk = 1, NNz
                    !Diagonal Block Gnd(x- and y-direction)
                    Gpd(ii, jj, xx) = Gpd(ii, jj, xx) + cmat1(ii+(xx-1)*NNz, kk)*conjg(GRLef(jj+(xx-1)*NNz, kk))
                    
                    !Nearest neighbor Gn (depth z-direction)
                    IF (xx<Npx) THEN
                        Gpu(ii, jj, xx) = Gpu(ii, jj, xx) + cmat1(ii+(xx-1)*NNz, kk)*conjg(GRLef(jj+xx*NNz, kk))

                        Gpl(ii, jj, xx) = Gpl(ii, jj, xx) + cmat1(ii+xx*NNz, kk)*conjg(GRLef(jj+(xx-1)*NNz, kk))
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    cmat1 = c0; 
    DO ii = 1, Ntot
        DO jj = 1, NNz
            DO kk = 1, NNz
                cmat1(ii,jj) = cmat1(ii,jj) + f2v*GRRgt(ii, kk)*GammaR(kk,jj)
            ENDDO
        ENDDO
    ENDDO
    !Gp = cmat1*GR'
     DO xx = 1, Npx
        DO ii = 1, NNz
            DO jj = 1, NNz
                DO kk = 1, NNz
                    !Diagonal Block Gnd(x- and y-direction)
                    Gpd(ii, jj, xx) = Gpd(ii, jj, xx) + cmat1(ii+(xx-1)*NNz, kk)*conjg(GRRgt(jj+(xx-1)*NNz, kk))
                    
                    !Nearest neighbor Gn (depth z-direction)
                    IF (xx<Npx) THEN
                        Gpu(ii, jj, xx) = Gpu(ii, jj, xx) + cmat1(ii+(xx-1)*NNz, kk)*conjg(GRRgt(jj+xx*NNz, kk))

                        Gpl(ii, jj, xx) = Gpl(ii, jj, xx) + cmat1(ii+xx*NNz, kk)*conjg(GRRgt(jj+(xx-1)*NNz, kk))
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    DEALLOCATE(cmat1)
    
    !!!!!!!on-site density matrix
    DO xx = 1, Npx
        DO yy = 1, Npy
            DO zz = 1, Npz
                DO ii = 1, Norb
                    index = ii+(yy-1)*Norb+(zz-1)*NN
                    RhoN(xx, yy, zz) = RhoN(xx, yy, zz) + real(Gnd(index, index, xx))*(deltaE/(2D0*pi))
                    RhoP(xx, yy, zz) = RhoP(xx, yy, zz) + real(Gpd(index, index, xx))*(deltaE/(2D0*pi))
                    IF (erCD == 1) THEN
                        ERRhon(xx, yy, zz, eCtr+1) = ERRhon(xx, yy, zz, eCtr+1) + real(Gnd(index, index, xx)) - real(Gpd(index, index, xx))
                    ENDIF
                    IF (orCD == 1) THEN
                        ORRhon(xx, yy, zz, ii) = ORRhon(xx, yy, zz, ii) + real(Gnd(index, index, xx))*(deltaE/(2D0*pi))
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    !!!!!!!on-site density matrix
END SUBROUTINE Calculate_Gnp
