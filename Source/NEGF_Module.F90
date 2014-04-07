MODULE  NEGF_Module
    USE Constants
    IMPLICIT NONE

    !Device parameters
    REAL(8)                     :: MM = 15D-1                   !Mass in dirac Hamiltonian (natural units)
    REAL(8)                     :: tMet = 1D0                   !hopping energy for metal (channel or contact)

    !Material Parameters Specific to Bi2Se3
    REAL(8), PARAMETER          :: mpA1 = 2.26                  !eV*Angstroms
    REAL(8), PARAMETER          :: mpA2 = 3.33                  !eV*Angstroms
    REAL(8), PARAMETER          :: mpC = -0.0083                !eV
    REAL(8), PARAMETER          :: mpD1 = 5.74                  !eV*Angstroms^2
    REAL(8), PARAMETER          :: mpD2 = 30.4                  !ev*Angstroms^2
    REAL(8), PARAMETER          :: mpM = 0.28                   !eV
    REAL(8), PARAMETER          :: mpB1 = 6.86                  !eV*Angstroms^2
    REAL(8), PARAMETER          :: mpB2 = 44.5                  !eV*Angstroms^2
    REAL(8)                     :: a0                           !1nm = 10 Angstroms (lattice constant)

    INTEGER                     :: channelType                  !0=natural units, 1=Bi2Se3, 2=metal, 3=QAHE
    INTEGER                     :: contactType                  !0=SI, 1=PBig, 2=PSmall

    !input parameters from command line
    INTEGER             :: nargs
    CHARACTER           :: arg*200, outputfiledir*200
    INTEGER             :: abStpnum, includePoisson
    INTEGER             :: Npx, Npy, Npz, Ntot, NN, NNx, NNz, Nh, Norb   !number of points per layer
    INTEGER             :: randSeed !Initial Seed for Random Impurity Disorder
    REAL                :: corrLength
    
    !iteration variables
    INTEGER             :: ite, iteMax

    !boundary conditions - 1 = periodic and 0 = open
    INTEGER             :: bcy, bcz

    !gamma matrices
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: g1,g2,g3,g0,gMx,gMy,gMz,id,AA,gDS

    !lead chemical potentials
    REAL(8), DIMENSION(2)       :: mu
    REAL(8)                     :: muT

    !Disorder variables
    REAL(8) :: disorderStrength   !Strength of nonmagnetic disorder added to on-site block
    REAL(8) :: magneticRingWidth
    REAL(8) :: magneticDisorderStrength
    INTEGER(4) :: isHollowRing

    !density matrix
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE    :: RhoN, RhoP          !rhon = integral of Gn(E), rhop = integral of Gp(E)  (from eMin to eMax)
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE    :: SigXExp, SigYExp, SigZExp, DeltaR
    INTEGER(4)                                :: calcSigExp

    !Hamiltonian matrices
    COMPLEX(8), DIMENSION(:,:,:), ALLOCATABLE :: Htotd
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE   :: HopX   !hopping matrices along inter-slize (x) direction
    
    !Superconductor Parameters (for SC contacts)
    INTEGER(4)                                :: isBdGHam=0   !BdG formalism necessary for superconducting contacts (double Hilbert space for inclusion of Nambu spinor)
    REAL(8)                                   :: chiLR              !Phase difference between sc contacts
    REAL(8)                                   :: DeltaSCL, DeltaSCR !Gap size of superconducting contacts
    
    !MPI variables
    INTEGER(4) :: Nprocs
    INTEGER(4) :: mpiRank, mpiSize, mpiRoot, mpierror
    INTEGER(4), ALLOCATABLE :: mpistatus(:)

    !Green's function matrices
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE   :: GRLef, GRRgt  !GRTop, GRBot
    COMPLEX(8), DIMENSION(:,:,:), ALLOCATABLE :: InvGRd, GRd, GRu, GRl, Gnd, Gpd, Gnu, Gpu, Gnl, Gpl

    !Self-energy/broadening matrices
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: GammaL, Gammav1, GammaR, Gammav2
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: SigmaL, SigmaR, Sigmav1, Sigmav2

    !Energy parameters
    INTEGER                                 :: eStpnum      !number of energy points (sum of biased/nonbiased energy sweeps)
    REAL(8)                                 :: energy, f1, f1v, f2, f2v
    REAL(8)                                 :: eMin, eMax, Emax0, deltaE, deltaE0
    REAL(8)                                 :: eAbsMin, eAbsMax
    REAL(8), DIMENSION(:), ALLOCATABLE      :: EAxis

    !Aharonov-Bohm  parameters
    REAL(8)                                 :: abPhi        !AB Phase from electromagnetic vector potential

    !direction-resolved current densities
    REAL(8), DIMENSION(:), ALLOCATABLE      :: ContactCurrent
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE  :: Jx, Jy, Jz

    !Transmission at each energy level from each contact to every other contact
    REAL(8), DIMENSION(:), ALLOCATABLE      :: T12
    REAL(8), DIMENSION(:), ALLOCATABLE      :: Modec1, Modec2
    REAL(8), DIMENSION(:), ALLOCATABLE      :: DOSc1, DOSc2

    !electron and hole density matrices, integral of rhon and rhop
    REAL(8), DIMENSION(:), ALLOCATABLE      :: EDensity, HDensity

    !Energy-resolved current densities
    REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE:: ERJx, ERJy, ERJz, ERRhon, ORJx, ORJy, ORJz, ORRhon
    INTEGER                                 :: erCD, orCD

    !temp/dummy variables
    INTEGER                                 :: eCtr, ii, jj, kk, xx, yy, zz
    INTEGER                                 :: errorFlag = 0
END MODULE NEGF_Module
