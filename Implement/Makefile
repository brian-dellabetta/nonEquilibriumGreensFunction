FC=mpif90
FFLAGS=-O2
LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LDFLAGS=-L/site/local/pkg/intel/Compiler/11.1/056/mkl/lib/em64t
SUBFILES=Read_In_Allocate.o Build_Hamiltonian.o Contact_Self_Energy.o Recursive_Greens.o SIP_Functions.o SIP_Routines.o SIP_Solver.o SIP_Coefficients.o MKL_Routines_Cluster.o Calculate_Observables.o Calculate_Gnp.o Dump_Data.o MPI_Routines.o Energy_Sweep.o

nonEquilibriumGreensFunction.exe: SIP_Module.o $(SUBFILES) NEGF_Main.o
	$(FC) -heap-arrays $(LDFLAGS) -o nonEquilibriumGreensFunction.exe Constants.o NEGF_Module.o SIP_Module.o $(SUBFILES) NEGF_Main.o $(LIBS)

SIP_Module.o: NEGF_Module.o SIP_Module.F90
	$(FC) $(FFLAGS) -c SIP_Module.F90

NEGF_Module.o: Constants.o NEGF_Module.F90
	$(FC) $(FFLAGS) -c NEGF_Module.F90

Constants.o: Constants.F90
	$(FC) $(FFLAGS) -c Constants.F90

MPI_Routines.o: SIP_Module.o MPI_Routines.F90
	$(FC) $(FFLAGS) -c MPI_Routines.F90

Energy_Sweep.o: SIP_Module.o Energy_Sweep.F90
	$(FC) $(FFLAGS) -c Energy_Sweep.F90

Recursive_Greens.o: SIP_Module.o Recursive_Greens.F90
	$(FC) $(FFLAGS) -c Recursive_Greens.F90

MKL_Routines_Cluster.o: Constants.o MKL_Routines_Cluster.F90
	$(FC) $(FFLAGS) -c MKL_Routines_Cluster.F90

Read_In_Allocate.o:  SIP_Module.o Read_In_Allocate.F90
	 $(FC) $(FFLAGS) -c Read_In_Allocate.F90

Build_Hamiltonian.o: SIP_Module.o Build_Hamiltonian.F90
	$(FC) $(FFLAGS) -c Build_Hamiltonian.F90

Calculate_Gnp.o: SIP_Module.o Calculate_Gnp.F90
	$(FC) $(FFLAGS) -c Calculate_Gnp.F90

Calculate_Observables.o: SIP_Module.o Calculate_Observables.F90
	$(FC) $(FFLAGS) -c Calculate_Observables.F90

SIP_Functions.o: Constants.o SIP_Functions.F90
	$(FC) $(FFLAGS) -c SIP_Functions.F90

SIP_Routines.o: SIP_Module.o SIP_Routines.F90
	$(FC) $(FFLAGS) -c SIP_Routines.F90

SIP_Solver.o: SIP_Module.o SIP_Solver.F90
	$(FC) $(FFLAGS) -c SIP_Solver.F90

SIP_Coefficients.o: SIP_Module.o SIP_Coefficients.F90
	$(FC) $(FFLAGS) -c SIP_Coefficients.F90

Contact_Self_Energy.o: SIP_Module.o Contact_Self_Energy.F90
	$(FC) $(FFLAGS) -c Contact_Self_Energy.F90

Dump_Data.o: SIP_Module.o Dump_Data.F90
	$(FC) $(FFLAGS) -c Dump_Data.F90

NEGF_Main.o: SIP_Module.o $(SUBFILES) NEGF_Main.F90
	$(FC) $(FFLAGS) -c NEGF_Main.F90 

clean:
	rm -f *.mod *.o >& /dev/null
