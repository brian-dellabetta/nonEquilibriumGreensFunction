nonEquilibriumGreensFunction
============================


This project calculates quantum-level transport dynamics of a fully three-dimensional model Topolgoical Insulator (expandable and customizable to other materials in the Build_Hamiltonian file) channel in the Non-Equilibrium Green's Function Formalism. The code includes the capability to insert magnetic (both via Zeeman splitting and via a quantum phase brought about by the electromagnetic vector potential) and nonmagnetic disorder.  The code uses MPI to determine observables in a massively parallelized fashion, and uses a recursive Greens function routine to invert the block tridiagonal Hamiltonian at O(Npx*(Npy*Npz)^3) speed rather than O((Npx*Npy*Npz)^3).  The Math Kernel Library (MKL) is necessary to run this code.

The Makefile and example shell submission wrappers are included in the ./Implement folder. 
Fortran90 source code is located in the ./Source folder. 
A simple Matlab script to read in and analyze the data output by the Dump_Data routine is located in the ./Analyze folder.


==============================
./Implement Folder contains:
==============================
Multiple_Run.sh:  This file sets up the inidividual runs that are submitted to the cluster using the Sun Grid Engine.  Input parameters are set in this file

Single_Run.sh:  This file sets up the grid engine environment for each simulation, and is called by Multiple_Run.sh.  The parameters for the SGE environment are set in this file

Makefile:  Example makefile that must be customized for a given compiler to link correctly to MKL and MPI libraries.

==============================
./Source Folder contains:
==============================
NEGF_Main.F90:  Main program that calls all other subroutines in the formalism

Read_In_Allocate.F90:  Reads in command line parameters and allocates and prepares the necessary arrays

SIP_{Main,Functions,Routines,Solver}.F90:  Calculates the electrostatic profile brought about by the electron and hole concentrations via a Poisson Solver

Build_Hamiltonian.F90:  Builds the real-space Hamiltonian for a customizable channel.

Energy_Sweep.F90:  Each slave process calculates the Green's function response for a given energy, returning the results, which are all integrated over energy, to the master node.

Contact_Self_Energy.F90:  Sets up the self-energy terms for the customiable contacts. Subroutine includes iterative (Sancho Rubio and Recursive Inversion) procedures to calculate the surface Green's function for semi-infinite contacts.

Get_Retarded_Greens.F90:  Calculates the necessary blocks of the inverted matrix G = (invG)^-1 = ((E+i\eta)I-H-\Sigma)^-1.  Only certain components of G are calculated, as only a sparse section of it is necessary for calculation of the observables in a nearest-neighbor tight-binding model.

Calculate_Gnp.F90:  Calculates the electron- and hole-correlated Green's functions (Gn and Gp) from the Retarded Green's function.

Calculate_Observables.F90:  Uses Gn and Gp to calculate transport observables - transmission, density profiles, spin expectation values.

Dump_Data.F90:  Outputs all observables to a series of files that can then be analyzed with Matlab routines, an example of which is located in the ./Analyze folder.

MKL_Routines_{PC,Cluster}.F90:  Math Kernel Library routines called by other subroutines.  PC is used when running on a local, single-process machine for debugging purposes.

============================
./Analyze Folder contains:
============================
NEGF_ReadData.m:  Reads in data output by the Dump_Data.F90 subroutine, and plots a few example observables like transmission and spatially-resolved current density. 


===============================
This source code was used for the following publications:
S. Cho, B. Dellabetta et al., "The Aharonov-Bohm Effect in a 3D Topological Insulator Nanowire", in preparation (2014)
Abstract at <http://meetings.aps.org/Meeting/MAR14/Session/T43.1>

S. Cho, B. Dellabetta et al., "Symmetry Protected Josephson Supercurrents in Three-Dimensional Topological Insulators"
Nature Communications 4, 1689 (2013) 

B. Dellabetta et al., "Imaging Topologically Protected Transport with Quantum Degenerate Gases"
Phys. Rev. B 85, 205442 (2012)