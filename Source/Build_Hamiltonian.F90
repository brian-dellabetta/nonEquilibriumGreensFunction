SUBROUTINE Build_Hamiltonian
    USE SIP_Module
    IMPLICIT NONE

    INTEGER  :: index
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE   :: HH2d, HH1d, HopY, HopZ

    ALLOCATE(HH2d(NNz, NNz))        ;HH2d = c0
    ALLOCATE(HH1d(NN, NN))          ;HH1d = c0
    ALLOCATE(HopY(Norb,Norb))       ;HopY = c0
    ALLOCATE(HopZ(NN,NN))           ;HopZ = c0

    IF (channelType == 0) THEN      !Topological Insulator Hamiltonian Unit Cell, assumes Norb=4
        MM = 15D-1
        AA(1:Nh,1:Nh) = (MM-3D0/(a0*a0))*g0 - muT*id(1:4,1:4)
        HopX(1:Nh,1:Nh) = 5D-1*(c1/(a0*a0)*g0 - ci/a0*g1)
        HopY(1:Nh,1:Nh) = 5D-1*(c1/(a0*a0)*g0 - ci/a0*g2)
        HopZ(1:Nh,1:Nh) = 5D-1*(c1/(a0*a0)*g0 - ci/a0*g3)
    ELSEIF (channelType == 1) THEN  !Bi2Se3 Hamiltonian Unit Cell, assumes Norb=4
        AA(1:Nh,1:Nh) = (mpC + (2D0*mpD1+4D0*mpD2)/(a0*a0))*id &
                      + (mpM - (2D0*mpB1+4D0*mpB2)/(a0*a0))*g0

        HopY(1:Nh,1:Nh) = mpA2/(2D0*a0)*ci*g2 &
                        - mpD2/(a0*a0)*id &
                        + mpB2/(a0*a0)*g0

        HopZ(1:Nh,1:Nh) = mpA1/(2D0*a0)*ci*g3 &
                        - mpD1/(a0*a0)*id &
                        + mpB1/(a0*a0)*g0

        HopX(1:Nh,1:Nh) = mpA2/(2D0*a0)*ci*g1 &
                        - mpD2/(a0*a0)*id &
                        + mpB2/(a0*a0)*g0
    ELSEIF (channelType == 2) THEN  !Single-band Metal Hamiltonian Unit Cell, assumes Norb=1
        AA(1:Nh,1:Nh) = 6D0*tMet
        HopY(1:Nh,1:Nh) = -tMet
        HopZ(1:Nh,1:Nh) = -tMet
        HopX(1:Nh,1:Nh) = -tMet
    ELSEIF (channelType == 3) THEN  !TI Surface state Hamiltonian, H(k) = k_x \sigma_y - k_y \sigma_x, assumbes Norb=2
        AA(1:Nh,1:Nh) = c0
        HopY(1:Nh,1:Nh) = 5D-1*ci*g1
        HopX(1:Nh,1:Nh) = -5D-1*ci*g2
        HopZ = c0
    ELSEIF (channelType == 4) THEN !1D Quantum wire with spin orbit coupling, H(kx) = k_x \sigma_y
        AA(1:Nh,1:Nh) = (-muT) * id(1:Nh,1:Nh)
        !x-direction, 1-d chain
        HopX(1:Nh,1:Nh) = -ci*5D-1*g2
    ELSE
        STOP "channelType not supported, must be >=0 and <=3"
    ENDIF
    
    IF (isBdGHam) THEN  !The 22 block of the BdG Hamiltonian is -conj(11 Block)
        AA(1+Nh:Norb,1+Nh:Norb)   = -conjg(AA(1:Nh,1:Nh))
        HopY(1+Nh:Norb,1+Nh:Norb) = -conjg(HopY(1:Nh,1:Nh))
        HopX(1+Nh:Norb,1+Nh:Norb) = -conjg(HopX(1:Nh,1:Nh))
        HopZ(1+Nh:Norb,1+Nh:Norb) = -conjg(HopZ(1:Nh,1:Nh))
    ENDIF

    DO ii = 1, Npy-1
        HopZ(ii*Norb+1:(ii+1)*Norb,ii*Norb+1:(ii+1)*Norb) = HopZ(1:Norb,1:Norb)
    ENDDO
    DO ii = 1, Npy*Npz-1
        HopX(ii*Norb+1:(ii+1)*Norb,ii*Norb+1:(ii+1)*Norb) = HopX(1:Norb,1:Norb)
    ENDDO

    !!!!!!!!!!Intraslice (y-direction) Interaction (1D)
    DO yy = 0, Npy-2
        HH1d(Norb*yy+1:Norb*(yy+1),Norb*yy+1:Norb*(yy+1)) = AA

        HH1d(Norb*yy+1:Norb*(yy+1),Norb*(yy+1)+1:Norb*(yy+2)) = HopY
        HH1d(Norb*(yy+1)+1:Norb*(yy+2),Norb*yy+1:Norb*(yy+1)) = conjg(transpose(HopY))
    ENDDO
    HH1d(NN-Norb+1:NN,NN-Norb+1:NN) = AA

    !Connect Npy -> 1 if periodic boundary conditions
    IF (bcy==1) THEN
        HH1d(1:Norb, NN-Norb+1:NN) = conjg(transpose(HopY))
        HH1d(NN-Norb+1:NN, 1:Norb) = HopY
    ENDIF
    !!!!!!!!!!Intraslice (y-direction) Interaction (1D)

    !!!!!!!!!!!!Interslice (z-direction) Interaction (2D)
    DO zz = 0,Npz-2
        HH2d(zz*NN+1:(zz+1)*NN,zz*NN+1:(zz+1)*NN) = HH1d
        HH2d(zz*NN+1:(zz+1)*NN,(zz+1)*NN+1:(zz+2)*NN)= HopZ
        HH2d((zz+1)*NN+1:(zz+2)*NN,zz*NN+1:(zz+1)*NN)= conjg(transpose(HopZ))
    ENDDO
    HH2d(NNz-NN+1:NNz,NNz-NN+1:NNz) = HH1d

    !Connect Npz -> 1 if periodic boundary conditions
    IF (bcz==1) THEN
        HH2d(1:NN,NN*(Npz-1)+1:NN*Npz) = conjg(transpose(HopZ))
        HH2d(NN*(Npz-1)+1:NN*Npz,1:NN) = HopZ
    ENDIF
    !!!!!!!!!!!!Interslice (z-direction) Interaction

    !!!!!!!!!!!!Interlayer (x-direction) Interaction (3D)
    DO xx = 1,Npx
        Htotd(:,:,xx) = HH2d
    ENDDO
    !!!!!!!!!!!!Interlayer (x-direction) Interaction

    !!!!!!!!!!!!Include electrostatics contribution from Poisson Equation
    IF (includePoisson==1) THEN
        DO xx = 1, Npx
            DO yy = 1, Npy
                DO zz = 1, Npz
                    CALL Insert_Impurity(xx,yy,zz,Phi(xx,yy,zz))
                ENDDO
            ENDDO
        ENDDO
    ENDIF
    !!!!!!!!!!!!Include electrostatics contribution from Poisson Equation

    !!!!!!!!!!!Add vector potential for B=B_x
    IF (abPhi > 1D-10) THEN
        CALL Add_AB_Phase
    ENDIF
    !!!!!!!!!!!Add vector potential for B=B_x

    !!!!!!!!!!!Add Correlated Disorder
    IF (dabs(disorderStrength) > 1D-20) THEN
        CALL Add_Correlated_Disorder
    ENDIF
    !!!!!!!!!!!Add Correlated Disorder
    
    !!!!!!!!!!!Add magnetic ring around surface for topological magneto-electric effect
    IF (magneticRingWidth > 0) THEN
        CALL Add_Magnetic_Ring
    ENDIF
    !!!!!!!!!!!Add magnetic ring around surface for topological magneto-electric effect

    DEALLOCATE(HH2d)
    DEALLOCATE(HH1d)
    DEALLOCATE(HopY)
    DEALLOCATE(HopZ)
END SUBROUTINE Build_Hamiltonian


!!!!!Raise on-site energy of vacancy (DOI 10.1007/s10825-006-0116-4)
SUBROUTINE Insert_Impurity(xPos, yPos, zPos, onsitePot)
    USE NEGF_Module
    IMPLICIT NONE
    INTEGER(4)   :: xPos, yPos, zPos, index
    REAL(8)      :: onsitePot

    index = (yPos-1)*Norb+(zPos-1)*NN
    Htotd(index+1:index+Nh, index+1:index+Nh, xPos) = Htotd(index+1:index+Nh, index+1:index+Nh, xPos) + onsitePot*id(1:Nh, 1:Nh)
    IF (isBdGHam) THEN
        Htotd(index+1+Nh:index+Norb, index+1+Nh:index+Norb, xPos) = -conjg(Htotd(index+1:index+Nh, index+1:index+Nh, xPos))
    ENDIF
END SUBROUTINE Insert_Impurity
!!!!!Raise on-site energy of vacancy (DOI 10.1007/s10825-006-0116-4)


SUBROUTINE Add_AB_Phase
    USE NEGF_Module
    IMPLICIT NONE
    INTEGER(4) :: index
    
    DO zz = 1, Npz
        DO xx = 1, Npx
            DO yy = 1, Npy-1
                index = (yy-1)*Norb+(zz-1)*NN

                Htotd(index+1:index+Nh,index+1+Norb:index+Norb+Nh,xx) = &
                        Htotd(index+1:index+Nh,index+1+Norb:index+Norb+Nh,xx)*zexp( ci*2D0*pi*real(zz-1)*abPhi/real(Npy-1)/real(Npz-1))
                        
                IF (isBdGHam) THEN
                    Htotd(index+1+Nh:index+Norb+Nh,index+1+Norb+Nh:index+Norb+2*Nh, xx) = -conjg(Htotd(index+1:index+Nh,index+1+Norb:index+Norb+Nh,xx))
                ENDIF

                Htotd(index+1+Norb:index+2*Norb,index+1:index+Norb,xx) = conjg(transpose(Htotd(index+1:index+Norb,index+1+Norb:index+2*Norb,xx)))
            ENDDO
        ENDDO
    ENDDO
END SUBROUTINE  Add_AB_Phase


!!!!!!!!!!!Add Correlated Disorder
SUBROUTINE Add_Correlated_Disorder
    USE SIP_Module
    IMPLICIT NONE
    REAL(8)  :: onsitePot, randVal
    INTEGER  :: index, kStpnum=100
    REAL(8)  :: qq, kmin, kmax, deltaK
    REAL(8), ALLOCATABLE :: KxAxis(:), KyAxis(:), KzAxis(:)
    COMPLEX(8), ALLOCATABLE :: DeltaQ(:,:,:)
    
    ALLOCATE(KxAxis(kStpnum))      ;KxAxis = 0D0
    ALLOCATE(KyAxis(kStpnum))      ;KyAxis = 0D0
    ALLOCATE(KzAxis(kStpnum))      ;KzAxis = 0D0
    ALLOCATE(DeltaQ(kStpnum,kStpnum,kStpnum)) ;DeltaQ = c0

    IF (dabs(disorderStrength) > 1D-20) THEN
        IF (ite == 1) THEN
            DeltaR = 0D0

            !Initialize PseudoRandom Number Generator for Correlated Random Impurities on top/bottom surfaces
            IF (randSeed == 1) THEN
                CALL RANDOM_SEED(PUT = (/1216546/)) !1 = 1216546, 2 =3598616, 3 = 9876248, 4 = 5248596
            ELSEIF (randSeed == 2) THEN
                CALL RANDOM_SEED(PUT = (/3598616/))
            ELSEIF (randSeed == 3) THEN
                CALL RANDOM_SEED(PUT = (/9876248/))
            ELSEIF (randSeed == 4) THEN
                CALL RANDOM_SEED(PUT = (/5248596/))
            ELSE
                randSeed = 5
                CALL RANDOM_SEED(PUT = (/6861389/))
            ENDIF

            IF (corrLength > 1D-12) THEN

                kmax = 1D1/corrLength; kmin = -kmax
                deltaK = (kmax-kmin)/(kStpnum-1);
                DO ii = 0, kStpnum-1
                    KxAxis(ii+1) = kmin + deltaK*ii
                    KyAxis(ii+1) = kmin + deltaK*ii
                    KzAxis(ii+1) = kmin + deltaK*ii
                ENDDO

                DO ii = 1, kStpnum
                    DO jj = 1, kStpnum
                        DO kk = 1, kStpnum
                            CALL RANDOM_NUMBER(randVal)
                            qq = dsqrt(KxAxis(ii)**2+KyAxis(jj)**2+KzAxis(kk)**2)
                            DeltaQ(ii,jj,kk) = zexp(-qq**2*corrLength**2/2D0 + ci*2D0*pi*randVal)
                        ENDDO
                    ENDDO
                ENDDO

                DO xx = 2, Npx-1
                  DO yy = 1, Npy
                    DO zz = 1, Npz
                      DO ii = 1, kStpnum
                        DO jj = 1, kStpnum
                          DO kk = 1, kStpnum
                              DeltaR(xx,yy,zz) = DeltaR(xx,yy,zz) + real(DeltaQ(ii,jj,kk)*zexp(-ci*(KxAxis(ii)*real(xx)+KyAxis(jj)*real(yy)+KzAxis(kk)*real(zz))))
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
                onsitePot = disorderStrength/maxval(maxval(maxval(abs(DeltaR),3),2),1)
                DeltaR = DeltaR * onsitePot
            ELSE
                DO xx = 2, Npx-1
                    DO yy = 1, Npy
                        DO zz = 1, Npz
                            CALL RANDOM_NUMBER(randVal)
                            DeltaR(xx,yy,zz) = (randVal-5D-1)*disorderStrength
                        ENDDO
                    ENDDO
                ENDDO
            ENDIF
        ENDIF

        DO xx = 2, Npx-1
            DO yy = 1, Npy
                DO zz = 1, Npz
                    CALL Insert_Impurity(xx,yy,zz,DeltaR(xx,yy,zz))
                ENDDO
            ENDDO
        ENDDO
    ENDIF
    
    DEALLOCATE(KxAxis)
    DEALLOCATE(KyAxis)
    DEALLOCATE(KzAxis)
    DEALLOCATE(DeltaQ)
END SUBROUTINE Add_Correlated_Disorder

!!!!!!!!!!!!Insert Ring of magnetic disorder into Hamiltonian
SUBROUTINE Add_Magnetic_Ring
    USE NEGF_Module
    IMPLICIT NONE
    INTEGER :: xMid
    xMid = int(real(Npx)/2D0)
    DO xx = xMid - int(magneticRingWidth/2), xMid + int(magneticRingWidth/2)
        DO yy = 1, Npy
            DO zz = 1, Npz
                IF (zz == 1) THEN
                    IF (((yy == 1) .OR. (yy == Npy)) .AND. (bcy == 0)) THEN
                        CALL Insert_Magnetic_Impurity(xx, yy, zz, 0D0, 0D0, magneticDisorderStrength/dsqrt(2D0))
                    ELSE
                        CALL Insert_Magnetic_Impurity(xx, yy, zz, 0D0, 0D0, magneticDisorderStrength)
                    ENDIF
                ELSEIF (zz == Npz) THEN
                    IF (((yy == 1) .OR. (yy == Npy)) .AND. (bcy == 0)) THEN
                        CALL Insert_Magnetic_Impurity(xx, yy, zz, 0D0, 0D0, -magneticDisorderStrength/dsqrt(2D0))
                    ELSE
                        CALL Insert_Magnetic_Impurity(xx, yy, zz, 0D0, 0D0, -magneticDisorderStrength)
                    ENDIF
                ENDIF

                IF ((yy == 1) .AND. (bcy == 0)) THEN
                    IF ((zz == 1) .OR. (zz == Npz)) THEN
                        CALL Insert_Magnetic_Impurity(xx, yy, zz, 0D0, magneticDisorderStrength/dsqrt(2D0), 0D0)
                    ELSE
                        CALL Insert_Magnetic_Impurity(xx, yy, zz, 0D0, magneticDisorderStrength, 0D0)
                    ENDIF
                ELSEIF ((yy == Npy) .AND. (bcy == 0)) THEN
                    IF ((zz == 1) .OR. (zz == Npz)) THEN
                        CALL Insert_Magnetic_Impurity(xx, yy, zz, 0D0, -magneticDisorderStrength/dsqrt(2D0), 0D0)
                    ELSE
                        CALL Insert_Magnetic_Impurity(xx, yy, zz, 0D0, -magneticDisorderStrength, 0D0)
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
    ENDDO
END SUBROUTINE Add_Magnetic_Ring


!!!!!Add magnetic disorder to location (xPos,yPos,zPos)
SUBROUTINE Insert_Magnetic_Impurity(xPos, yPos, zPos, magX, magY, magZ)
    USE SIP_Module
    IMPLICIT NONE
    INTEGER(4)            :: xPos, yPos, zPos, index
    REAL(8)               :: magX, magY, magZ

    index = (yPos-1)*Norb+(zPos-1)*NN
    Htotd(index+1:index+Nh, index+1:index+Nh, xPos) = Htotd(index+1:index+Nh, index+1:index+Nh, xPos) + magX*gMx + magY*gMy + magZ*gMz
    
    IF (isBdGHam) THEN
        Htotd(index+1+Nh:index+Norb, index+1+Nh:index+Norb, xPos) = -conjg(Htotd(index+1:index+Nh, index+1:index+Nh, xPos) )
    ENDIF
END SUBROUTINE Insert_Magnetic_Impurity
!!!!!Add magnetic disorder to location (xPos,yPos,zPos)


