PROGRAM NEGF_Main
    USE SIP_Module
    IMPLICIT NONE

    CALL Read_In_Allocate
    
    DO ite = 1, iteMax
        Htotd = c0; HopX = c0;
        
        IF (mpiRank == mpiRoot) THEN
            IF (includePoisson==1) THEN
                CALL SIP_Main   !Calculate 3-D Voltage profile Phi as a function of RhoN, RhoP
            ENDIF

            CALL Build_Hamiltonian !Build Htotd, HopX
        ENDIF
        
        RhoN = c0; RhoP = c0;
        CALL MPI_Broadcast !Broadcast mpiRoot information to non-root processes for energy integration

        DO eCtr = 0, eStpnum-1
            !Split up energy integration over each mpi process
            IF ((mpiSize == 1) .OR. (mod(eCtr, mpiSize) == mpiRank)) THEN
                energy = EAxis(eCtr+1)
                CALL Energy_Sweep
            ENDIF
        ENDDO
        
        !Gather Data from each process, store in mpiRoot process
        CALL MPI_Gather_Iteration
        
        IF (ite==iteMax .OR. poissonIsConverged==1) THEN
            CALL MPI_Gather_Final

            IF (mpiRank==mpiRoot) THEN
                CALL Dump_Data !dump data to file and stop
            ENDIF
        ENDIF
    ENDDO
END PROGRAM NEGF_Main