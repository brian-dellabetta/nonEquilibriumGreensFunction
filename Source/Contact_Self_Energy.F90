SUBROUTINE Contact_Self_Energy
    USE NEGF_Module
    IMPLICIT NONE
    
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: alpha, beta
    REAL(8)                                 :: contactEnergy, tunnelingRate
    COMPLEX(8)                              :: BetaEps, dumVar1, dumVar2, malpha, dalpha

    tunnelingRate = 1D0             !phenomenological treatment of contacts (0<rate<=1)
    
    ALLOCATE(alpha(NNz,NNz));     alpha = c0
    ALLOCATE(beta(NNz,NNz));      beta = c0
    
    !include small imaginary term for smoother MBS/Transmission peaks
    contactEnergy = cmplx(energy, eta_sigma, 8)
    
    !Self energies from boundaries with contacts, NOT HERMITIAN!
    IF (contactType == 1) THEN  !Phenomenological, big contact
        DO ii = 1, NNz
            SigmaL(ii,ii) = -ci*tunnelingRate
        ENDDO
        SigmaR = SigmaL
    ELSEIF (contactType == 2) THEN  !Phenomenological, small contact
        DO ii = 1, NN
            SigmaL(ii,ii) = -ci*tunnelingRate
        ENDDO
        SigmaR = SigmaL
    ELSEIF (contactType == 3) THEN  !Superconducting contacts (From Sun et al., Phys. Rev. B 62, 648-660 (2000)). 
        !From PRB 62, 648 (2000)
        dumVar1 = DeltaSCL**2D0-contactEnergy**2D0
        dumVar1 = zsqrt( cmplx( real(dumVar1), imag(dumVar1), 8) )
        
        malpha = -tunnelingRate/2D0 * (contactEnergy/dumVar1)
        dalpha =  tunnelingRate/2D0 * (DeltaSCL/dumVar1)
        
        !Semi-infinite superconducting contacts, alpha = E*Id - Hsc
        SigmaL(1:Nh,1:Nh) = malpha*id(1:Nh,1:Nh)
        SigmaL(Nh+1:Norb,Nh+1:Norb) = malpha*id(1:Nh,1:Nh)
        
        SigmaL(1:Nh,Nh+1:Norb) = dalpha*zexp(cmplx(0D0, chiLR,8))*(gDS)
        SigmaL(Nh+1:Norb,1:Nh) = conjg(transpose(SigmaL(1:Nh,Nh+1:Norb)))     !dalpha*zexp(cmplx(0D0,-PhiL,8))*(gDA+gDB)

        
        dumVar1 = DeltaSCR**2D0-contactEnergy**2D0
        dumVar1 = zsqrt( cmplx( real(dumVar1), imag(dumVar1), 8) )
        
        malpha = -tunnelingRate/2D0 * (contactEnergy/dumVar1)
        dalpha =  tunnelingRate/2D0 * (DeltaSCR/dumVar1)
        
        SigmaR(1:Nh,1:Nh) = malpha*id(1:Nh,1:Nh)
        SigmaR(Nh+1:Norb,Nh+1:Norb) = malpha*id(1:Nh,1:Nh)
        
        SigmaR(1:Nh,Nh+1:Norb) = dalpha*zexp(cmplx(0D0, 0D0,8))*(gDS)
        SigmaR(Nh+1:Norb,1:Nh) = conjg(transpose(SigmaR(1:Nh,Nh+1:Norb)))     !dalpha*zexp(cmplx(0D0,-PhiR,8))*(gDA+gDB)
    
    ELSE                       !Contacts are a semi-infinite extension of the channel
        alpha = c0; beta = c0; errorFlag=0
        
        alpha = -Htotd(:,:,1)
        DO ii = 1, NNz
            alpha(ii,ii) = alpha(ii,ii)+contactEnergy
        ENDDO
        beta = -1D0*conjg(transpose(HopX))
        CALL Sancho_Rubio(NNz, alpha, beta, SigmaL, errorFlag)
        IF (errorFlag /= 0) THEN
            WRITE(*, *) "Error in SR algorithm, SigmaL not converging, Energy = ", energy
        ENDIF

        alpha = c0; beta = c0; errorFlag=0
        
        alpha = -Htotd(:,:,Npx)
        DO ii = 1, NNz
            alpha(ii,ii) = alpha(ii,ii)+contactEnergy
        ENDDO
        beta = -1D0*(HopX)
        CALL Sancho_Rubio(NNz, alpha, beta, SigmaR, errorFlag)
        IF (errorFlag /= 0) THEN
            WRITE(*, *) "Error in SR algorithm, SigmaR not converging, Energy = ", energy
        ENDIF
    ENDIF

     DO ii=1,NNz
        DOSc1(eCtr+1) = DOSc1(eCtr+1)-imag(SigmaL(ii,ii))/pi
        DOSc2(eCtr+1) = DOSc2(eCtr+1)-imag(SigmaR(ii,ii))/pi
    ENDDO

    !Compute gamma (broadening) terms needed to calculate transmission
    GammaL = ci*(SigmaL-conjg(transpose(SigmaL)))
    GammaR = ci*(SigmaR-conjg(transpose(SigmaR)))

    DEALLOCATE(alpha)
    DEALLOCATE(beta)

END SUBROUTINE Contact_Self_Energy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Sancho_Rubio(NN, alpha, beta, Sigma, errorFlag)
    USE Constants
    IMPLICIT NONE
    
    !From J. Phys. F Vol 14, 1984, 1205, M P Lopez Sancho, J. Rubio
    INTEGER(4)  :: NN
    COMPLEX(8)  :: alpha(NN,NN), beta(NN,NN), Sigma(NN,NN)
    INTEGER(4)  :: errorFlag
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: tmp, betad, t, tt, CapT, Toldt, IdMat, cmat1, cmat2, change, gfinal
    INTEGER(4)  :: counter=1, ii, jj
    

    ALLOCATE(gfinal(NN,NN))          ;gfinal = c0
    ALLOCATE(tmp(NN,NN))             ;tmp = c0
    ALLOCATE(betad(NN,NN))           ;betad = c0
    ALLOCATE(t(NN,NN))               ;t = c0
    ALLOCATE(tt(NN,NN))              ;tt = c0
    ALLOCATE(CapT(NN,NN))            ;CapT = c0
    ALLOCATE(Toldt(NN,NN))           ;Toldt = c0
    ALLOCATE(IdMat(NN,NN))           ;IdMat = c0
    DO ii = 1,NN
        IdMat(ii, ii) = c1
    ENDDO
    ALLOCATE(cmat1(NN,NN))           ;cmat1 = c0
    ALLOCATE(cmat2(NN,NN))           ;cmat2 = c0
    ALLOCATE(change(NN,NN))          ;change = c1
    
    !betad = conjugate transpose of beta
    betad = conjg(transpose(beta))

    !Eq. 8:  tmp = ginit = inv(alpha + i*eye*eta_sigma
    CALL Invert_Matrix(tmp, alpha + IdMat(1:NN, 1:NN)*cmplx(0D0, eta_sigma, 8), NN, errorFlag)
    IF (errorFlag /= 0) THEN
        WRITE(*, *) "Error in Invert_Matrix algorithm, Contact_Self_Energy routine"
    ENDIF
    
    !Eq. 8 (t = -tmp*betad)
    CALL zgemm('N', 'N', NN, NN, NN, -c1, tmp, NN, betad, NN, c0, t, NN)
    !Eq. 8 (tt = -tmp*beta)
    CALL zgemm('N', 'N', NN, NN, NN, -c1, tmp, NN, beta, NN, c0, tt, NN)
    !First term in Eq. 16
    CapT = t
    Toldt = IdMat

    DO WHILE ((counter <= 20000) .AND. (maxval(abs(change))>1D-8))  !Sancho Rubio converges when change matrix goes to zero
        counter = counter+1
        !Toldt = Toldt*tt   (Product of tilde t in subsequent terms in Eq. 16)
        CALL zgemm('N', 'N', NN, NN, NN, c1, Toldt, NN, tt, NN, c0, cmat1, NN)
        Toldt = cmat1
        !tmp = inv(Id - t*tt - tt*t)    Inverse part of Eq. 12
        cmat1 = IdMat
        CALL zgemm('N', 'N', NN, NN, NN, -c1, t, NN, tt, NN, c1, cmat1, NN)
        CALL zgemm('N', 'N', NN, NN, NN, -c1, tt, NN, t, NN, c1, cmat1, NN)
        CALL Invert_Matrix(tmp, cmat1, NN, errorFlag)
        IF (errorFlag) THEN
            WRITE(*, *) "Error in tmp = inv(Id - t*tt - tt*t) calculation, Sancho_Rubio routine, tmp = ", tmp
            STOP "Error in tmp = inv(Id - t*tt - tt*t) calculation, Sancho_Rubio routine"
        ENDIF
        !t = tmp*t*t    Eq.12 (t_i)
        cmat1 = c0
        CALL zgemm('N', 'N', NN, NN, NN, c1, tmp, NN, t, NN, c0, cmat1, NN)
        CALL zgemm('N', 'N', NN, NN, NN, c1, cmat1, NN, t, NN, c0, cmat2, NN)
        t = cmat2
        !tt = tmp*tt*tt     Eq.12 (t_i tilde)
        CALL zgemm('N', 'N', NN, NN, NN, c1, tmp, NN, tt, NN, c0, cmat1, NN)
        CALL zgemm('N', 'N', NN, NN, NN, c1, cmat1, NN, tt, NN, c0, cmat2, NN)
        tt = cmat2
        !change = Toldt*t   Next term of Eq. 16
        CALL zgemm('N', 'N', NN, NN, NN, c1, Toldt, NN, t, NN, c0, change, NN)
        !T=T+change         Add to T, Eq. 16
        CapT = CapT + change

        DO ii = 1,NN
            DO jj = 1,NN
                IF (isnan(abs(change(ii,jj))) == 1) THEN
                    WRITE(*, *) "Error in Sancho Rubio method, NaN found"
                    WRITE(*, *) "counter, ii, jj, change(ii,jj) are", counter, ii, jj, change(ii,jj)
                    errorFlag = 1
                    RETURN
                ENDIF
            ENDDO
        ENDDO
    ENDDO

    !g = inv(alpha+beta*T)      Eq. 18
    CALL zgemm('N', 'N', NN, NN, NN, c1, beta, NN, CapT, NN, c0, cmat1, NN)
    cmat1 = cmat1 + alpha
    CALL Invert_Matrix(gfinal, cmat1, NN, errorFlag)
    IF (errorFlag) THEN
        WRITE(*, *) "Error in g=inv(alpha+beta*T) algorithm, Sancho_Rubio routine"
    ENDIF
    
    IF (counter == 20000) THEN
        WRITE(*,*) "Error in SR method, does not converge"
        errorFlag = 1
    ENDIF
    
    !Sigma = beta*gfinal*beta';
    CALL zgemm('N', 'N', NN, NN, NN, c1, beta, NN, gfinal, NN, c0, cmat1, NN)
    CALL zgemm('N', 'C', NN, NN, NN, c1, cmat1, NN, beta, NN, c0, Sigma, NN)
    
    DEALLOCATE(gfinal)
    DEALLOCATE(tmp)
    DEALLOCATE(betad)
    DEALLOCATE(t)
    DEALLOCATE(tt)
    DEALLOCATE(CapT)
    DEALLOCATE(Toldt)
    DEALLOCATE(IdMat)
    DEALLOCATE(cmat1)
    DEALLOCATE(cmat2)
    DEALLOCATE(change)

END SUBROUTINE Sancho_Rubio


!!!!!!!!!!!!Recursive Inversion routine is necessary when contact is 2 unit cells (e.g. graphene)
SUBROUTINE Recursive_Inversion(NN, alpha1, beta1, alpha2, beta2, Sigma, errorFlag)
    USE Constants
    IMPLICIT NONE
    
    INTEGER(4)  :: NN
    COMPLEX(8)  :: alpha1(NN,NN), beta1(NN,NN), alpha2(NN,NN), beta2(NN,NN), Sigma(NN,NN)
    INTEGER(4)  :: errorFlag, ite, ii, jj
    REAL(8)     :: eps, epstmp
    REAL(8), PARAMETER                      :: damping = 5D-1
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: g_new1, g_old1, g_new2, g_old2, cmat1, cmat2, gfinal
    
    ALLOCATE(gfinal(NN,NN))    ;gfinal = c0
    ALLOCATE(g_new1(NN,NN))    ;g_new1 = c0
    ALLOCATE(g_old1(NN,NN))    ;g_old1 = c0
    ALLOCATE(g_new2(NN,NN))    ;g_new2 = c0
    ALLOCATE(g_old2(NN,NN))    ;g_old2 = c0
    ALLOCATE(cmat1(NN,NN))     ;cmat1 = c0
    ALLOCATE(cmat2(NN,NN))     ;cmat2 = c0
    
    CALL Invert_Matrix(g_old1, alpha1, NN, errorFlag)
    CALL Invert_Matrix(g_old2, alpha2, NN, errorFlag)
    ite = 1; eps = 1D0

    DO WHILE ((eps > 1D-6) .AND. (ite<=100000))
        cmat1 = c0; cmat2 = c0

        !g_new2 = inv(alpha2 - beta1'*g_old1*beta1);
        CALL zgemm('C', 'N', NN, NN, NN, c1, beta1, NN, g_old1, NN, c0, cmat1, NN)
        CALL zgemm('N', 'N', NN, NN, NN, c1, cmat1, NN, beta1, NN, c0, cmat2, NN)
        CALL Invert_Matrix(g_new2, alpha2-cmat2, NN, errorFlag)
        g_old2 = damping*g_new2 + (1D0-damping)*g_old2
        
        !g_new1 = inv(alpha1 - beta2'*g_old2*beta2);
        CALL zgemm('C', 'N', NN, NN, NN, c1, beta2, NN, g_old2, NN, c0, cmat1, NN)
        CALL zgemm('N', 'N', NN, NN, NN, c1, cmat1, NN, beta2, NN, c0, cmat2, NN)
        CALL Invert_Matrix(g_new1, alpha1-cmat2, NN, errorFlag)
    
        !eps is the maximum deviation in magnitude between like elements of g_new1 and g_old1
        IF (ite > 200) THEN
            eps = 0D0
            DO ii = 1,NN
                DO jj = 1,NN
                    epstmp = real(g_new1(ii, jj)-g_old1(ii,jj))/real(g_new1(ii,jj))
                    IF (abs(epstmp) > eps) eps = abs(epstmp)
                    epstmp = imag(g_new1(ii, jj)-g_old1(ii,jj))/imag(g_new1(ii,jj))
                    IF (abs(epstmp) > eps) eps = abs(epstmp)
                    epstmp = real(g_new2(ii, jj)-g_old2(ii,jj))/real(g_new2(ii,jj))
                    IF (abs(epstmp) > eps) eps = abs(epstmp)
                    epstmp = imag(g_new2(ii, jj)-g_old2(ii,jj))/imag(g_new2(ii,jj))
                    IF (abs(epstmp) > eps) eps = abs(epstmp)
                ENDDO
            ENDDO
        ENDIF

        g_old1 = damping*g_new1 + (1D0-damping)*g_old1
        ite = ite+1
        IF (ite == 100000) THEN
            errorFlag = 1
        ENDIF
    ENDDO
    gfinal = g_old1

    !Sigma = beta1*gfinal*beta1';
    CALL zgemm('N', 'N', NN, NN, NN, c1, beta1, NN, gfinal, NN, c0, cmat1, NN)
    CALL zgemm('N', 'C', NN, NN, NN, c1, cmat1, NN, beta1, NN, c0, Sigma, NN)
    
    DEALLOCATE(gfinal)
    DEALLOCATE(g_new1)
    DEALLOCATE(g_old1)
    DEALLOCATE(g_new2)
    DEALLOCATE(g_old2)
    DEALLOCATE(cmat1)
    DEALLOCATE(cmat2)
    
END SUBROUTINE Recursive_Inversion
