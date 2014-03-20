MODULE SIP_Module
    USE NEGF_Module
    IMPLICIT NONE
    
    !Temperature related quantities
    REAL(8), PARAMETER  :: T0=300D0             !Use room temperature for poisson solver

    !!!!Bi2Se3 Parameters
    REAL(8), PARAMETER  :: epsBi = 113D0*eps0   !Bi2Se3 dielectric constant
    
    !!!!SiO2 Parameters
    REAL(8), PARAMETER  :: epsGate = 1D0*eps0   !Gate dielectric (SiO2 dielectric constant)
    REAL(8), PARAMETER  :: epsBox = 1D0*eps0    !Oxide dielectric
    !!!!End SiO2 Parameters
    
    !!!!GaAs Parameters    (material params from: http://www.ioffe.rssi.ru/SVA/NSM/Semicond)
    !effective masses
    REAL(8), PARAMETER  :: ml_gaas = 0.063*m0, mt_gaas = 0.063*m0
    REAL(8), PARAMETER  :: mlh_gaas = 0.082*m0, mhh_gaas = 0.51*m0
    REAL(8), PARAMETER  :: mdh_gaas = ((mlh_gaas/m0)**1.5+(mhh_gaas/m0)**1.5)**(2D0/3D0)*m0
    REAL(8), PARAMETER  :: mde_gaas = ((ml_gaas/m0)*((mt_gaas/m0)**2D0))**(1D0/3D0)*m0
    REAL(8), PARAMETER  :: mc_gaas = 3D0/(1D0/ml_gaas+2D0/mt_gaas)
    !Band gap energy and density of states
    REAL(8), PARAMETER  :: Eg_gaas = 1.519-((5.405D-4)*(T0**2)/(T0+204D0))
    REAL(8), PARAMETER  :: Nc_gaas = 2*(2*pi*(mde_gaas/(hbar*2*pi))*(kB*T0/(hbar*2*pi)))**1.5	
    REAL(8), PARAMETER  :: Nv_gaas = 2*(2*pi*(mdh_gaas/(hbar*2*pi))*(kB*T0/(hbar*2*pi)))**1.5
    REAL(8)             :: delta_Ec = 0.237		!conduction band offset from Si to BOX (in eV)     
    !Intrinsic Fermi level w.r.t. Ec, normalized by Vt
    REAL(8)             :: dEc_gaas, Norm_Eg_gaas
    REAL(8)             :: Ncnorm_gaas, Nvnorm_gaas
    REAL(8), PARAMETER  :: ni_gaas = sqrt(Nc_gaas)*sqrt(Nv_gaas)*exp(-Eg_gaas*q/(2D0*kB*T0))
    !!!!End GaAs Parameters
    
    !!!!Mesh Parameters
    REAL(8), PARAMETER                      ::  VTg=0D0, VBg=-0D0            !Top, bottom gate voltages
    REAL(8), PARAMETER                      ::  b0=1D0                       !Lattice constant in z
    REAL(8)                                 ::  hx, hy, hz
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE  ::  Phi, ro, epsMat    !3-D voltage, density, and dielectric constant profiles
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE  ::  doping, NominalDoping       !3-D Doping matrices
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE  ::  n, p    !3-D electron and hole density matrices
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE  ::  materialType   !flag type:  1 = oxide, 2 = Bi2Se3, 3-8 = Contacts
    !!!!End Mesh Parameters
        
    !!!!Parameters for Stone's SIP algorithm
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE  :: Bvec, Cvec, Dvec, Evec, Fvec, Gvec, Hvec  !coefficient arrays
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE  :: R_vec, PhiDelta, V_vec
    INTEGER, PARAMETER                      :: alphaArrayLength=6
    REAL(8), DIMENSION(1:alphaArrayLength)  :: alphaArray
    DATA                                       alphaArray(1:alphaArrayLength) /.25,.35,.45,.35,.25,.15/
    REAL(8)                                 :: alpha
    INTEGER                                 :: poissonIsConverged = 0 !will be set to 1 if solution has converged
    !!!!End Parameters used for Stone's SIP algorithm

END MODULE SIP_Module