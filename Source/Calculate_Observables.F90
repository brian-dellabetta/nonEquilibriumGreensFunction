!Calculate Transmission at each contact, modes, and current at each point
SUBROUTINE Calculate_Observables
    USE NEGF_Module
    IMPLICIT NONE
    INTEGER(4) :: index 
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: cmat1, cmat2
    REAL(4)    :: dumVar, sigmaZFactor
    
    !!!!!Spin expectation values
    IF (calcSigExp) THEN
        DO xx = 1, Npx
            DO yy = 1, Npy
                DO zz = 1, Npz
                    index = (yy-1)*Norb+(zz-1)*NN
                    DO ii = 1, Nh
                        DO jj = 1, Nh
                            SigXExp(xx,yy,zz) = SigXExp(xx,yy,zz) + deltaE/(2D0*pi) * &
                                    Gnd(index+ii, index+jj, xx)*gMx(jj,ii)
        
                            SigYExp(xx,yy,zz) = SigYExp(xx,yy,zz) + deltaE/(2D0*pi) * &
                                    Gnd(index+ii, index+jj, xx)*gMy(jj,ii)
        
                            SigZExp(xx,yy,zz) = SigZExp(xx,yy,zz) + deltaE/(2D0*pi) * &
                                    Gnd(index+ii, index+jj, xx)*gMz(jj,ii)
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDIF
    !!!!!Spin expectation values
    
    !!!!!!!!electron and hole DOS as a function of energy
    DO xx = 1, Npx
        DO yy = 1, Npy
            DO zz = 1, Npz
                DO ii = 1, Nh
                    EDensity(eCtr+1) = EDensity(eCtr+1) + real(Gnd(ii+(yy-1)*Norb+(zz-1)*NN, ii+(yy-1)*Norb+(zz-1)*NN, xx))/(2D0*pi)
                    HDensity(eCtr+1) = HDensity(eCtr+1) + real(Gpd(ii+(yy-1)*Norb+(zz-1)*NN, ii+(yy-1)*Norb+(zz-1)*NN, xx))/(2D0*pi)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    !!!!!!!!!!!!!!electron and hole DOS as a function of energy
    
    !!!!!!!!!!injected modes vs. energy
    !Modec1(eCtr+1) = trace(GammaL*(Gn+Gp)); Modec2(eCtr+1) = trace(GammaR*(Gn+Gp));
    DO ii = 1, NNz
        DO kk = 1, NNz
            Modec1(eCtr+1) = Modec1(eCtr+1) + real(GammaL(ii,kk)*(Gnd(kk,ii,1)+Gpd(kk,ii,1)))
            Modec2(eCtr+1) = Modec2(eCtr+1) + real(GammaR(ii,kk)*(Gnd(kk,ii,Npx)+Gpd(kk,ii,Npx)))
        ENDDO
    ENDDO
    !!!!!!!!!!injected modes vs. energy

    !!!!!!!!!transmissions as functions of energy, this routine calculates Josephson currents if SC contacts
    !T12(eCtr+1) = real(trace(GammaL*GR*GammaR*GR'));
    ALLOCATE(cmat1(NNz, NNz));  cmat1 = c0
    ALLOCATE(cmat2(NNz, NNz));  cmat2 = c0
    
    IF (isBdGHam) THEN
        !Supercurrent transmission from Equation 14, PRB 62 648 (2000):
        !T12(Ectr+1) = real(trace(SigmaLT * GR + Sigma * GLT)); SigmaLT = -Gamma, GLT = ci*Gn
        DO ii = 1, NNz
            DO jj = 1, NNz
                DO kk = 1, NNz
                    cmat1(ii,jj) = cmat1(ii,jj) + (ci*Gnd(ii, kk, 1))*conjg(SigmaL(jj, kk))

                    cmat2(ii,jj) = cmat2(ii,jj) + GRd(ii, kk, 1)*(-GammaL(kk, jj))
                ENDDO
            ENDDO
        ENDDO
        DO ii = 1, NNz
            IF (mod(ii-1, Norb) < Norb/2) THEN
                sigmaZFactor = 1D0
            ELSE
                sigmaZFactor = -1D0
            ENDIF
            !2-terminal device, T12 = T1 = -T2   
            T12(eCtr+1) = T12(eCtr+1) + sigmaZFactor*real(cmat1(ii, ii)+cmat2(ii, ii))
        ENDDO

        ContactCurrent(1) = ContactCurrent(1) + q*deltaE/h*( T12(eCtr+1));
        ContactCurrent(2) = ContactCurrent(2) + q*deltaE/h*(-T12(eCtr+1));
    ELSE
        DO ii = 1, NNz
            DO jj = 1, NNz
                DO kk = 1, NNz
                    !GammaL*GR
                    cmat1(ii,jj) = cmat1(ii,jj) + GammaL(ii, kk)*GRRgt(kk, jj)
                    !GammaR*GR'
                    cmat2(ii,jj) = cmat2(ii,jj) + GammaR(ii, kk)*conjg(GRRgt(jj, kk))
                ENDDO
            ENDDO
        ENDDO
        DO ii = 1, NNz
            DO jj = 1, NNz
                T12(eCtr+1) = T12(eCtr+1) + real(cmat1(ii, jj)*cmat2(jj,ii))
            ENDDO
        ENDDO

        ContactCurrent(1) = ContactCurrent(1) + q*deltaE/h*( T12(eCtr+1)*(f2-f1) );
        ContactCurrent(2) = ContactCurrent(2) + q*deltaE/h*( T12(eCtr+1)*(f1-f2) );
    ENDIF
    
    DEALLOCATE(cmat1)
    DEALLOCATE(cmat2)
    !!!!!!!!!transmissions as functions of energy

    !!!!!!!!!Spatially Resolved Current Profile
    DO xx = 1, Npx
        DO yy = 1, Npy
            DO zz = 1, Npz
                index = (yy-1)*Norb+(zz-1)*NN
                DO ii = 1, Nh
                    DO jj = 1, Nh
                        !y-direction
                        IF (yy < Npy) THEN
                            dumVar = real(ci*q*q/hbar* ( Htotd(index+ii, index+Norb+jj, xx)*Gnd(index+Norb+jj, index+ii, xx)   -&
                                                           Gnd(index+ii, index+Norb+jj, xx)*Htotd(index+Norb+jj, index+ii, xx) -&
                                                         Htotd(index+ii, index+Norb+jj, xx)*Gpd(index+Norb+jj, index+ii, xx) +&
                                                         Gpd(index+ii, index+Norb+jj, xx)*Htotd(index+Norb+jj, index+ii, xx)  ))
                            
                            IF (erCD == 1) THEN
                                ERJy(xx, yy, zz, eCtr+1) = ERJy(xx, yy, zz, eCtr+1) + dumVar
                            ENDIF
                            IF (orCD == 1) THEN
                                ORJy(xx, yy, zz, (ii-1)*Norb+jj) = ORJy(xx, yy, zz, (ii-1)*Norb+jj) + dumVar*deltaE/(2D0*pi)
                            ENDIF
                            Jy(xx,yy,zz) = Jy(xx,yy,zz) + dumVar*deltaE/(2D0*pi)
                        ENDIF
                           
                        !z-direction  
                        IF (zz < Npz) THEN
                            dumVar = real(ci*q*q/hbar* ( Htotd(index+ii, index+NN+jj, xx)*Gnd(index+NN+jj, index+ii, xx) -&
                                                         Gnd(index+ii, index+NN+jj, xx)*Htotd(index+NN+jj, index+ii, xx) -&
                                                         Htotd(index+ii, index+NN+jj, xx)*Gpd(index+NN+jj, index+ii, xx) +&
                                                         Gpd(index+ii, index+NN+jj, xx)*Htotd(index+NN+jj, index+ii, xx)  ))
                                                                       
                            IF (erCD == 1) THEN
                                ERJz(xx, yy, zz, eCtr+1) = ERJz(xx, yy, zz, eCtr+1) + dumVar
                            ENDIF
                            IF (orCD == 1) THEN
                                ORJz(xx, yy, zz, (ii-1)*Norb+jj) = ORJz(xx, yy, zz, (ii-1)*Norb+jj) + dumVar*deltaE/(2D0*pi)
                            ENDIF
                            Jz(xx,yy,zz) = Jz(xx,yy,zz) + dumVar*deltaE/(2D0*pi)
                        ENDIF
                        
                        !x-direction
                        IF (xx < Npx) THEN
                            dumVar = real(ci*q*q/hbar* ( HopX(index+ii, index+jj)*Gnl(index+jj, index+ii, xx)         -&
                                                         Gnu(index+ii, index+jj, xx)*conjg(HopX(index+ii, index+jj)) -&
                                                         HopX(index+ii, index+jj)*Gpl(index+jj, index+ii, xx)        +&
                                                         Gpu(index+ii, index+jj, xx)*conjg(HopX(index+ii, index+jj))  ))
                                                                    
                            IF (erCD == 1) THEN
                                ERJx(xx, yy, zz, eCtr+1) = ERJx(xx, yy, zz, eCtr+1) + dumVar
                            ENDIF
                            IF (orCD == 1) THEN
                                ORJx(xx, yy, zz, (ii-1)*Norb+jj) = ORJx(xx, yy, zz, (ii-1)*Norb+jj) + dumVar*deltaE/(2D0*pi)
                            ENDIF
                            Jx(xx,yy,zz) = Jx(xx,yy,zz) + dumVar*deltaE/(2D0*pi)
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    !!!!!!!!!Spatially Resolved Current Profile
END SUBROUTINE Calculate_Observables