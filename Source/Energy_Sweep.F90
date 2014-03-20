SUBROUTINE Energy_Sweep
    USE SIP_Module
    IMPLICIT NONE
    
    IF (mpiRank == mpiRoot) WRITE(*,*) eCtr
    
    IF (eCtr /= (eStpnum-1)) THEN
        deltaE = EAxis(eCtr+2)-EAxis(eCtr+1)  !recalculate dE for integral in case of non-uniform mesh
    ENDIF

    IF (energy <= muT-mu(1)) THEN
        f1=1D0;
    ELSE
        f1=0D0;
    ENDIF

    IF (energy <= muT-mu(2)) THEN
        f2=1D0;
    ELSE
        f2=0D0;
    ENDIF
    f1v = 1D0-f1
    f2v = 1D0-f2

    GRd = c0; GRu = c0; GRl = c0; InvGRd = c0; GRLef = c0; GRRgt = c0; 
    SigmaL = c0; SigmaR = c0; GammaL = c0; GammaR = c0;

    CALL Contact_Self_Energy    !Calculate self-energy/broadening terms (Sigmas, Gammas)

    !invGR = (E+i*eta_sigma)*IdMat - Htot - (SigmaL+SigmaR+Sigmav1+Sigmav2)
    DO xx = 1, Npx
        InvGRd(:,:,xx) = -Htotd(:,:,xx)
        DO ii = 1, NNz
            InvGRd(ii, ii, xx) = InvGRd(ii, ii, xx) + cmplx(energy, eta_sigma, 8)
        ENDDO
    ENDDO
    InvGRd(:,:,1) = InvGRd(:,:,1) - SigmaL
    InvGRd(:,:,Npx) = InvGRd(:,:,Npx) - SigmaR

    CALL Get_Retarded_Greens    !Calculate GR = inv(invGR)

    CALL Calculate_Gnp  !Calculate Gn, Gp, rhon, rhop

    IF (ite==itemax .OR. poissonIsConverged==1) THEN
        CALL Calculate_Observables
    ENDIF
END SUBROUTINE Energy_Sweep