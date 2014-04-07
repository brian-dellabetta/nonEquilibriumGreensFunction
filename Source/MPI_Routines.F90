SUBROUTINE MPI_Initialize
    USE NEGF_Module
    USE MPI
    IMPLICIT NONE
    
    IF (Nprocs > 1) THEN
        CALL MPI_Init( mpierror ); CALL MPI_ErrorCheck(mpierror, 1)
        CALL MPI_Comm_size( MPI_COMM_WORLD, mpisize, mpierror ); CALL MPI_ErrorCheck(mpierror, 1)
        CALL MPI_Comm_Rank( MPI_COMM_WORLD, mpirank, mpierror ); CALL MPI_ErrorCheck(mpierror, 1)
        ALLOCATE (mpistatus(MPI_STATUS_SIZE))
    ELSE
        mpisize = 1
        mpirank = 0
    ENDIF
    mpiroot = 0
END SUBROUTINE


SUBROUTINE MPI_Broadcast
    USE NEGF_Module
    USE MPI
    IMPLICIT NONE
    
    IF (Nprocs > 1) THEN
        CALL MPI_BARRIER(MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 2)
        
        CALL MPI_BCAST(Htotd, NNz*NNz*Npx, MPI_COMPLEX16, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 2)
        CALL MPI_BCAST(HopX, NNz*NNz, MPI_COMPLEX16, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 2)
    ENDIF
END SUBROUTINE MPI_Broadcast


SUBROUTINE MPI_Gather_Iteration
    USE NEGF_Module
    USE MPI
    IMPLICIT NONE
    
    REAL(8), ALLOCATABLE :: rmat1(:), rmat2(:,:), rmat3(:,:,:), rmat4(:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: cmat1(:), cmat2(:,:), cmat3(:,:,:), cmat4(:,:,:,:)
    
    IF (Nprocs > 1) THEN
        CALL MPI_BARRIER(MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)

        !Reduce observables to process mpiroot
        ALLOCATE(rmat3(Npx, Npy, Npz))
        rmat3 = RhoN
        CALL MPI_REDUCE(rmat3, RhoN, Npx*Npy*Npz, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
        rmat3 = RhoP
        CALL MPI_REDUCE(rmat3, RhoP, Npx*Npy*Npz, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
        DEALLOCATE(rmat3)
    ENDIF
END SUBROUTINE MPI_Gather_Iteration


SUBROUTINE MPI_Gather_Final
    USE NEGF_Module
    USE MPI
    
    REAL(8), ALLOCATABLE :: rmat1(:), rmat2(:,:), rmat3(:,:,:), rmat4(:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: cmat1(:), cmat2(:,:), cmat3(:,:,:), cmat4(:,:,:,:)

    IF (Nprocs > 1) THEN
        ALLOCATE(rmat1(eStpnum))
        rmat1 = EDensity
        CALL MPI_REDUCE(rmat1, EDensity, eStpnum, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
        rmat1 = HDensity
        CALL MPI_REDUCE(rmat1, HDensity, eStpnum, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)

        rmat1 = T12
        CALL MPI_REDUCE(rmat1, T12, eStpnum, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
        rmat1 = Modec1
        CALL MPI_REDUCE(rmat1, Modec1, eStpnum, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
        rmat1 = Modec2
        CALL MPI_REDUCE(rmat1, Modec2, eStpnum, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
        rmat1 = DOSc1
        CALL MPI_REDUCE(rmat1, DOSc1, eStpnum, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
        rmat1 = DOSc2
        CALL MPI_REDUCE(rmat1, DOSc2, eStpnum, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)

        rmat1(1:2) = contactCurrent
        CALL MPI_REDUCE(rmat1(1:2), contactCurrent, 2, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
        DEALLOCATE(rmat1)
        
        ALLOCATE(rmat3(Npx,Npy,Npz))
        rmat3 = Jx
        CALL MPI_REDUCE(rmat3, Jx, Npx*Npy*Npz, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
        rmat3 = Jy
        CALL MPI_REDUCE(rmat3, Jy, Npx*Npy*Npz, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
        rmat3 = Jz
        CALL MPI_REDUCE(rmat3, Jz, Npx*Npy*Npz, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
        DEALLOCATE(rmat3)
        
        IF (ERcd == 1) THEN
            ALLOCATE(rmat4(Npx,Npy,Npz,eStpnum))
            rmat4 = ERJx
            CALL MPI_REDUCE(rmat4, ERJx,   Npx*Npy*Npz*eStpnum, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
            rmat4 = ERJy
            CALL MPI_REDUCE(rmat4, ERJy,   Npx*Npy*Npz*eStpnum, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
            rmat4 = ERJz
            CALL MPI_REDUCE(rmat4, ERJz,   Npx*Npy*Npz*eStpnum, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
            rmat4 = ERRhon
            CALL MPI_REDUCE(rmat4, ERRhon, Npx*Npy*Npz*eStpnum, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
            DEALLOCATE(rmat4)
        ENDIF
        IF (ORcd == 1) THEN
            ALLOCATE(rmat4(Npx,Npy,Npz,Norb*Norb))
            rmat4 = ORJx
            CALL MPI_REDUCE(rmat4, ORJx,   Npx*Npy*Npz*Norb*Norb, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
            rmat4 = ORJy
            CALL MPI_REDUCE(rmat4, ORJy,   Npx*Npy*Npz*Norb*Norb, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
            rmat4 = ORJz
            CALL MPI_REDUCE(rmat4, ORJz,   Npx*Npy*Npz*Norb*Norb, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
            rmat4(:,:,:,1:Norb) = ORRhon
            CALL MPI_REDUCE(rmat4(:,:,:,1:Norb), ORRhon, Npx*Npy*Npz*Norb, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
            DEALLOCATE(rmat4)
        ENDIF
        IF (calcSigExp == 1) THEN
            ALLOCATE(rmat3(Npx, Npy, Npz))
            rmat3 = SigXExp
            CALL MPI_REDUCE(rmat3, SigXExp, Npx*Npy*Npz, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
            rmat3 = SigYExp
            CALL MPI_REDUCE(rmat3, SigYExp, Npx*Npy*Npz, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
            rmat3 = SigZExp
            CALL MPI_REDUCE(rmat3, SigZExp, Npx*Npy*Npz, MPI_DOUBLE_PRECISION, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 3)
            DEALLOCATE(rmat3)
        ENDIF
    
        CALL MPI_BARRIER(MPI_COMM_WORLD, mpierror); CALL MPI_ErrorCheck(mpierror, 4)
        
        CALL MPI_FINALIZE(mpierror); CALL MPI_ErrorCheck(mpierror, 4)
        
        IF (mpirank /= mpiroot) STOP
    ENDIF
END SUBROUTINE MPI_Gather_Final


SUBROUTINE MPI_ErrorCheck(mpierror, condition)
    IMPLICIT NONE
    INTEGER :: mpierror, condition

    IF (mpierror /= 0) THEN
        IF (condition == 1) THEN
            WRITE(*,*) "MPI Error in Initialize MPI Subroutine.  mpierror = ", mpierror
        ELSEIF (condition == 2) THEN
            WRITE(*,*) "MPI Error in Broadcast MPI Subroutine.  mpierror = ", mpierror
        ELSEIF (condition == 3) THEN
            WRITE(*,*) "MPI Error in GatherData MPI Subroutine.  mpierror = ", mpierror
        ELSEIF (condition == 4) THEN
            WRITE(*,*) "MPI Error in Finalize MPI Subroutine.  mpierror = ", mpierror
        ELSE
            WRITE(*,*) "MPI Error in unknown MPI Subroutine.  mpierror = ", mpierror
        ENDIF
    ENDIF  !If error occurred
END SUBROUTINE MPI_ErrorCheck