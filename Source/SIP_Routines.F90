SUBROUTINE SIP_Main
    USE SIP_Module
    IMPLICIT NONE

    poissonIsConverged = 0  !reset poisson solver to rerun
    IF (ite == 1) THEN  !first ite, go through initialization
        CALL Initialize_Grid
        CALL Initialize_Doping
        CALL SIP_Coefficients
    ELSE
        CALL Solve_Equilibrium
    ENDIF
    
END SUBROUTINE SIP_Main

!!!!Initialize uniform mesh
SUBROUTINE Initialize_Grid
    USE SIP_Module
    IMPLICIT NONE
    
    ALLOCATE(Phi(-1:Npx,-1:Npy,-1:Npz));               Phi = 0D0
    ALLOCATE(PhiDelta(-1:Npx,-1:Npy,-1:Npz));          PhiDelta = 0D0
    ALLOCATE(ro(0:Npx-1,0:Npy-1,0:Npz-1));             ro = 0D0
    ALLOCATE(materialType(0:Npx-1,0:Npy-1,0:Npz-1));   materialType = 0
    ALLOCATE(epsMat(0:Npx-1,0:Npy-1,0:Npz-1));         epsMat = 0D0
    ALLOCATE(n(0:Npx-1,0:Npy-1,0:Npz-1));              n = 0D0
    ALLOCATE(p(0:Npx-1,0:Npy-1,0:Npz-1));              p = 0D0
        
    !Material Type:  1 = SiO2 Oxide, 2 = Bi2Se3, 3 = Left Gate, 4 = Right Gate
    
    !BULK REGION
    materialType(1:Npx-2,1:Npy-2,1:Npz-2) = 2
    epsMat(1:Npx-2,1:Npy-2,1:Npz-2) = epsBi/epsBi
    
    !left face contact (mu(1)):
    materialType(0,:,:) = 3
    epsMat(0,:,:) = epsGate/epsBi
    !right face contact (mu(2)):
    materialType(Npx-1,:,:) = 4
    epsMat(Npx-1,:,:) = epsGate/epsBi
    
    hx=a0
    hy=a0
    hz=a0
END SUBROUTINE Initialize_Grid
    

SUBROUTINE Initialize_Doping
    USE SIP_Module
    IMPLICIT NONE
    
    !doping array
    ALLOCATE(doping(0:Npx-1,0:Npy-1,0:Npz-1));           doping = 0D0
    ALLOCATE(NominalDoping(0:Npx-1,0:Npy-1,0:Npz-1));    NominalDoping = 0D0
    
    !System is undoped
    
    NominalDoping = doping  !uniform to be used for boundary conditions
END SUBROUTINE Initialize_Doping