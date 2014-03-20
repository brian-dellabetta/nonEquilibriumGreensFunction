MODULE Constants
    IMPLICIT NONE
    COMPLEX(8), PARAMETER :: ci = cmplx(0D0,1D0,8), c0 = cmplx(0D0,0D0,8), c1 = cmplx(1D0,0D0,8)
    REAL(8), PARAMETER    :: q = 1.602D-19, m0 = 9.1019D-31, pi = 4D0*atan(1D0)
    REAL(8), PARAMETER    :: eps0 = 8.85418782D-12, kB = 1.38066D-23
    REAL(8), PARAMETER    :: h = 6.626D-34, hbar = h/(2D0*pi)
    REAL(8), PARAMETER    :: eta_sigma = 1D-4
END MODULE Constants
