SUBROUTINE Dump_Data
    USE SIP_Module
    IMPLICIT NONE
    
    !Open files to store scalar data
    OPEN(unit=13, file=trim(outputfiledir)//"/scalars.txt")  
    WRITE(13, *) Npx
    WRITE(13, *) Npy
    WRITE(13, *) Npz
    WRITE(13, *) Norb
    WRITE(13, 43) muT
    WRITE(13, 43) eMin
    WRITE(13, 43) eMax
    WRITE(13, *) eStpnum
    WRITE(13, 43) deltaE
    DO ii = 1,2
        WRITE(13, 43) mu(ii)
    ENDDO  
    WRITE(13, *) channelType
    IF (channelType == 1) THEN
        WRITE(13, *) a0
    ENDIF
    WRITE(13, *) abPhi
    WRITE(13, 43) corrLength
    WRITE(13, 43) disorderStrength
    CLOSE(13)
    
    OPEN(unit=13, file=trim(outputfiledir)//"/Greens.txt")  
    !Write matrices from NEGF formalism
    DO xx = 1, Npx
        DO yy = 1, Npy
            DO zz = 1, Npz
                WRITE(13, 43) real(RhoN(xx,yy,zz))
                WRITE(13, 43) real(RhoP(xx,yy,zz))
            ENDDO
        ENDDO
    ENDDO
    
    DO xx = 1, Npx
        DO yy = 1, Npy
            DO zz = 1, Npz
                WRITE(13, 43) DeltaR(xx,yy,zz)
            ENDDO
        ENDDO
    ENDDO
    CLOSE(13)
        
    OPEN(unit=13, file=trim(outputfiledir)//"/Transmissions.txt")  
    !electron/hole Density
    DO ii = 1, eStpnum
        WRITE(13, 43) EDensity(ii)
        WRITE(13, 43) HDensity(ii)
    ENDDO
    !Transmissions
    DO ii = 1, eStpnum
        WRITE(13, 43) T12(ii)
    ENDDO
    !Modec1, Modec2, Modev1, Modev2
    DO ii = 1, eStpnum
        WRITE(13, 43) Modec1(ii)
        WRITE(13, 43) Modec2(ii)
    ENDDO
    DO ii = 1, eStpnum
        WRITE(13, *) EAxis(ii)
    ENDDO
    DO ii = 1,2
        WRITE(13, *) ContactCurrent(ii)
    ENDDO
    DO xx = 1, Npx
        DO yy = 1, Npy
            DO zz = 1, Npz
                WRITE(13, *) Jx(xx, yy, zz)
                WRITE(13, *) Jy(xx, yy, zz)
                WRITE(13, *) Jz(xx, yy, zz)
            ENDDO
        ENDDO
    ENDDO
    DO ii = 1, eStpnum
        WRITE(13, 43) DOSc1(ii)
        WRITE(13, 43) DOSc2(ii)
    ENDDO
    CLOSE(13)
    
    IF (includePoisson) THEN
        OPEN(unit=13, file=trim(outputfiledir)//"/Poissons.txt")  
        !Values from poisson solver
        DO ii = 1, Npx
            DO jj = 1, Npy
                DO kk = 1, Npz
                    WRITE(13, 43) Phi(ii-1, jj-1, kk-1)

                    WRITE(13, 43) ro(ii-1, jj-1, kk-1)
                ENDDO
            ENDDO
        ENDDO
        CLOSE(13)   
    ENDIF
    
    
    IF (erCD == 1) THEN
        OPEN(unit=13, file=trim(outputfiledir)//"/ERCD.txt")  
        DO xx = 1, Npx
            DO yy = 1, Npy
                DO zz = 1, Npz
                    DO ii = 1, eStpnum
                        WRITE(13, 43) ERJx(xx, yy, zz, ii)
                        WRITE(13, 43) ERJy(xx, yy, zz, ii)
                        WRITE(13, 43) ERJz(xx, yy, zz, ii)
                        WRITE(13, 43) ERRhon(xx, yy, zz, ii)
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDIF
    
    IF (orCD == 1) THEN
        OPEN(unit=13, file=trim(outputfiledir)//"/ORCD.txt")  
        DO xx = 1, Npx
            DO yy = 1, Npy
                DO zz = 1, Npz
                    DO ii = 1, Norb*Norb
                        WRITE(13, 43) ORJx(xx, yy, zz, ii)
                        WRITE(13, 43) ORJy(xx, yy, zz, ii)
                        WRITE(13, 43) ORJz(xx, yy, zz, ii)
                    ENDDO
                    DO ii = 1, Norb
                        WRITE(13, 43) ORRhon(xx, yy, zz, ii)
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDIF
    
    IF (calcSigExp) THEN
        OPEN(unit=13, file=trim(outputfiledir)//"/SigExp.txt")  
        !Write spin densities
        DO xx = 1, Npx
            DO yy = 1, Npy
                DO zz = 1, Npz
                    WRITE(13, 43) real(SigXExp(xx,yy,zz))
                    WRITE(13, 43) real(SigYExp(xx,yy,zz))
                    WRITE(13, 43) real(SigZExp(xx,yy,zz))
                ENDDO
            ENDDO
        ENDDO
        CLOSE(13)
    ENDIF

    WRITE(*,*)"COMPLETE:  Data dumped to directory ", outputfiledir
    STOP "Run Complete!"

43 format(ES11.4)
END SUBROUTINE Dump_Data