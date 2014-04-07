SUBROUTINE Read_In_Allocate
    USE NEGF_Module
    IMPLICIT NONE

    !Read in parameters
    nargs = iargc()
    IF (nargs == 19) THEN
        CALL getarg(1, arg); READ (arg, *) Nprocs
        CALL MPI_Initialize
        CALL getarg(2, arg); READ (arg, *) Npx
        CALL getarg(3, arg); READ (arg, *) Npy
        CALL getarg(4, arg); READ (arg, *) Npz
        CALL getarg(5, arg); READ (arg, *) bcy
        CALL getarg(6, arg); READ (arg, *) bcz
        CALL getarg(7, arg); READ (arg, *) mu(1)
        CALL getarg(8, arg); READ (arg, *) mu(2)
        CALL getarg(9, arg); READ (arg, *) muT
        CALL getarg(10, arg); READ (arg, *) channelType !0=toy TI, 1=Bi2Se3, 2=toy metal, 3=TI surface state, 4=1d wire with spin-orbit coupling
        CALL getarg(11, arg); READ (arg, *) contactType !0 = SI, 1 = PBig, 2 = PSmall, 3 = Superconducting
        CALL getarg(12, arg); READ (arg, *) a0
        CALL getarg(13, arg); READ (arg, *) abPhi
        CALL getarg(14, arg); READ (arg, *) disorderStrength
        CALL getarg(15, arg); READ (arg, *) corrLength
        CALL getarg(16, arg); READ (arg, *) randSeed
        CALL getarg(17, arg); READ (arg, *) includePoisson
        CALL getarg(18, arg); READ (arg, *) iteMax
        CALL getarg(19, outputfiledir)

    ELSE                    !default Values
        WRITE(*,*) "Insufficient command line arguments"
        WRITE(*,*) "nargs = ", nargs, "Using default values..."
        Nprocs = 1
        CALL MPI_Initialize

        Npx = 4; Npy = 3; Npz = 3
        bcy = 0; bcz = 0
        mu(1) = -1D-1
        mu(2) = 0D0
        muT = 0D0
        channelType = 0
        contactType = 0 !0
        a0 = 1D0
        disorderStrength = 0D0
        corrLength = 0D0
        randSeed = 1
        abPhi = 0D0
        outputfiledir = 'TRDebug/'
        includePoisson = 0
        iteMax = 10
    ENDIF
    IF (channelType == 2) THEN      !one-orbital metal with surface states
        Norb = 1
    ELSEIF (channelType == 3) THEN  !2D quantum hall effect
        Norb = 2
        Npz = 1  !no z hopping term
        bcz=0
    ELSEIF (channelType == 4) THEN  !1D Quantum wire, H(kx) = k_x*\sigma_y
        Norb = 2
        Npy = 1; Npz = 1
        bcy = 0; bcz = 0
    ELSEIF(channelType == 0 .OR. channelType == 1) THEN  !four-orbital 3D TI (toy model==0, Bi2Se3 model == 1)
        Norb = 4
    ELSE
        STOP "channelType must be >= 0 and <=4"
    ENDIF
    
    Nh = Norb
    
    IF (contactType == 3) THEN  !Superconducting contacts require BdG Hamiltonian
        isBdGHam = 1
        Norb = Norb*2
        chiLR = 8D-1       !Josephson phase between superconducting contacts
        mu = 0D0           !Formalism only set up to calculate Josephson current between SCs with equivalent chemical potential
        DeltaSCL = 1D-3
        DeltaSCR = 1D-3
    ENDIF
    
    IF (abPhi > 1D-10 .AND. (bcy .OR. bcz)) THEN
        WRITE(*,*) "Error - cannot include Aharanov-Bohm phase to periodic structure.  Setting abPhi=0..."
        abPhi=0D0
    ENDIF
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!Other variables to set
    eStpnum = 600       !number of discretized points in energy space
    
    magneticRingWidth = 0         !Creates a magnetic ring of width (magRing) in center of channel to study Topological Magnetoelectric Effect
    magneticDisorderStrength = 0D0
    
    isHollowRing = 0    !Flag to insert vacanacies in the bulk to model trivial surface state for AB effect
    
    erCD = 0            !Flag for energy resolved current density
    orCD = 0            !Flag for orbital resolved current densit
    calcSigExp = 0      !Flag to calculate spin densities
    !!!!!!!!!!!!!!!!!!!!!!!!!!Other variables to set
    
    NN  = Npy*Norb
    NNz = Npz*NN
    Ntot = Npx*Npy*Npz*Norb
    IF (includePoisson == 0) THEN
        iteMax = 1      !No self-consistent iterative procedure necessary
    ENDIF  

    !Allocate arrays...this is done to prevent stack overflow (automatic arrays are placed on the stack)
    ALLOCATE(Htotd(NNz,NNz,Npx))            ;Htotd = c0
    ALLOCATE(HopX(NNz,NNz))                 ;HopX = c0

    ALLOCATE(InvGRd(NNz,NNz,Npx))           ;InvGRd = c0

    ALLOCATE(GRd(NNz,NNz,Npx))              ;GRd = c0
    ALLOCATE(GRu(NNz,NNz,Npx-1))            ;GRu = c0
    ALLOCATE(GRl(NNz,NNz,Npx-1))            ;GRl = c0
    ALLOCATE(GRLef(Ntot, NNz))               ;GRLef = c0
    ALLOCATE(GRRgt(Ntot, NNz))               ;GRRgt = c0

    ALLOCATE(Gnd(NNz,NNz,Npx))              ;Gnd = c0
    ALLOCATE(Gpd(NNz,NNz,Npx))              ;Gpd = c0
    ALLOCATE(Gnu(NNz,NNz,Npx-1))            ;Gnu = c0
    ALLOCATE(Gpu(NNz,NNz,Npx-1))            ;Gpu = c0
    ALLOCATE(Gnl(NNz,NNz,Npx-1))            ;Gnl = c0
    ALLOCATE(Gpl(NNz,NNz,Npx-1))            ;Gpl = c0

    ALLOCATE(RhoN(Npx, Npy, Npz))           ;RhoN = 0D0
    ALLOCATE(RhoP(Npx, Npy, Npz))           ;RhoP = 0D0

    IF (calcSigExp) THEN
        ALLOCATE(SigXExp(Npx, Npy, Npz)) ;SigXExp = c0
        ALLOCATE(SigYExp(Npx, Npy, Npz)) ;SigYExp = c0
        ALLOCATE(SigZExp(Npx, Npy, Npz)) ;SigZExp = c0
    ENDIF

    ALLOCATE(GammaL(NNz,NNz))               ;GammaL = c0
    ALLOCATE(GammaR(NNz,NNz))               ;GammaR = c0
    ALLOCATE(SigmaL(NNz,NNz))               ;SigmaL = c0
    ALLOCATE(SigmaR(NNz,NNz))               ;SigmaR = c0

    ALLOCATE(EDensity(eStpnum))            ;EDensity = 0D0
    ALLOCATE(HDensity(eStpnum))            ;HDensity = 0D0
    ALLOCATE(T12(eStpnum))                 ;T12 = 0D0
    ALLOCATE(Modec1(eStpnum))              ;Modec1 = 0D0
    ALLOCATE(Modec2(eStpnum))              ;Modec2 = 0D0
    ALLOCATE(DOSc1(eStpnum))               ;DOSc1 = 0D0
    ALLOCATE(DOSc2(eStpnum))               ;DOSc2 = 0D0

    ALLOCATE(EAxis(eStpnum))                ;EAxis = 0D0
    ALLOCATE(Jx(Npx, Npy, Npz))             ;Jx = 0D0
    ALLOCATE(Jy(Npx, Npy, Npz))             ;Jy = 0D0
    ALLOCATE(Jz(Npx, Npy, Npz))             ;Jz = 0D0
    ALLOCATE(DeltaR(Npx,Npy,Npz))           ;DeltaR = 0D0
    ALLOCATE(ContactCurrent(2))             ;ContactCurrent = 0D0

    !Set up Gamma matrices
    ALLOCATE(g1(Nh,Nh)) ;g1=c0
    ALLOCATE(g2(Nh,Nh)) ;g2=c0
    ALLOCATE(g3(Nh,Nh)) ;g3=c0
    ALLOCATE(g0(Nh,Nh)) ;g0=c0
    ALLOCATE(gMx(Nh,Nh)) ;gMx=c0
    ALLOCATE(gMy(Nh,Nh)) ;gMy=c0
    ALLOCATE(gMz(Nh,Nh)) ;gMz=c0
    ALLOCATE(gDS(Nh,Nh)) ;gDS=c0

    IF (Nh==4) THEN
        !!!!!!!!!!!!!!!Chen's Basis for toy TI Model
        g1(1,:) = (/ c0,  c0,  c0, -ci/)
        g1(2,:) = (/ c0,  c0,  ci,  c0/)
        g1(3,:) = (/ c0, -ci,  c0,  c0/)
        g1(4,:) = (/ ci,  c0,  c0,  c0/)

        g2(1,:) = (/ c0,  c0,  c0, -c1/)
        g2(2,:) = (/ c0,  c0, -c1,  c0/)
        g2(3,:) = (/ c0, -c1,  c0,  c0/)
        g2(4,:) = (/-c1,  c0,  c0,  c0/)

        g3(1,:) = (/c0,  c0, -ci,  c0/)
        g3(2,:) = (/c0,  c0,  c0, -ci/)
        g3(3,:) = (/ci,  c0,  c0,  c0/)
        g3(4,:) = (/c0,  ci,  c0,  c0/)

        g0(1,:) = (/c1,  c0,  c0,  c0/)
        g0(2,:) = (/c0,  c1,  c0,  c0/)
        g0(3,:) = (/c0,  c0, -c1,  c0/)
        g0(4,:) = (/c0,  c0,  c0, -c1/)

        gMx(1,:) = (/c0,  c1,  c0,  c0/)
        gMx(2,:) = (/c1,  c0,  c0,  c0/)
        gMx(3,:) = (/c0,  c0,  c0,  c1/)
        gMx(4,:) = (/c0,  c0,  c1,  c0/)

        gMy(1,:) = (/c0, -ci,  c0,  c0/)
        gMy(2,:) = (/ci,  c0,  c0,  c0/)
        gMy(3,:) = (/c0,  c0,  c0, -ci/)
        gMy(4,:) = (/c0,  c0,  ci,  c0/)

        gMz(1,:) = (/c1,  c0,  c0,  c0/)
        gMz(2,:) = (/c0, -c1,  c0,  c0/)
        gMz(3,:) = (/c0,  c0,  c1,  c0/)
        gMz(4,:) = (/c0,  c0,  c0, -c1/)
        !!!!!!!!!!!!!!!Chen's Basis for toy TI Model
        
        gDS(1,:) = (/ c0, c1, c0, c0/)
        gDS(2,:) = (/-c1, c0, c0, c0/)
        gDS(3,:) = (/ c0, c0, c0, c1/)
        gDS(4,:) = (/ c0, c0,-c1, c0/)
    ELSEIF (Nh==2) THEN

        g1(1,:) = (/c0,  c1/)
        g1(2,:) = (/c1,  c0/)

        g2(1,:) = (/c0, -ci/)
        g2(2,:) = (/ci,  c0/)

        g3(1,:) = (/c1,  c0/)
        g3(2,:) = (/c0, -c1/)
        
        gDS(1,:) = (/ c0, c1/)
        gDS(2,:) = (/-c1, c0/)

        gMx = g1
        gMy = g2
        gMz = g3
    ENDIF

    ALLOCATE(AA(Norb,Norb))                 ;AA=c0
    ALLOCATE(id(Norb,Norb)) ;id=c0
    DO ii = 1, Norb
        id(ii,ii) = c1
    ENDDO

    IF (erCD == 1) THEN
        ALLOCATE(ERJx(Npx,Npy,Npz, eStpnum))       ;ERJx=0D0
        ALLOCATE(ERJy(Npx,Npy,Npz, eStpnum))       ;ERJy=0D0
        ALLOCATE(ERJz(Npx,Npy,Npz, eStpnum))       ;ERJz=0D0
        ALLOCATE(ERRhon(Npx,Npy,Npz, eStpnum))     ;ERRhon=0D0
    ENDIF
    IF (orCD == 1) THEN
        ALLOCATE(ORJx(Npx,Npy,Npz, Norb*Norb))       ;ORJx=0D0
        ALLOCATE(ORJy(Npx,Npy,Npz, Norb*Norb))       ;ORJy=0D0
        ALLOCATE(ORJz(Npx,Npy,Npz, Norb*Norb))       ;ORJz=0D0
        ALLOCATE(ORRhon(Npx,Npy,Npz, Norb))          ;ORRhon=0D0
    ENDIF

    !Set up energy range
    IF (includePoisson == 1) THEN
        eAbsMin = 0D0 !muT-(maxval(mu))+1D-10 !-1D0
        eAbsMax = muT-(minval(mu))-1D-10 !1D0
    ELSEIF (contactType == 3) THEN
        eAbsMin = -11D-1*dabs(DeltaSCR)
        eAbsMax =  11D-1*dabs(DeltaSCR)
    ELSE
        eAbsMin = muT-(maxval(mu))-1D-10 !-maxval(dabs(mu))+1D-10  !
        eAbsMax = muT-(minval(mu))+1D-10 !maxval(dabs(mu))-1D-10 !
    ENDIF

    eMin = eAbsMin
    eMax = eAbsMax
    deltaE = (eMax-eMin)/(eStpnum-1);
    DO eCtr = 0, eStpnum-1
        energy = eMin+deltaE*eCtr
        EAxis(eCtr+1) = energy
    ENDDO

END SUBROUTINE Read_In_Allocate
