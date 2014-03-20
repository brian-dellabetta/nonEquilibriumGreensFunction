FUNCTION Fermi(X)
	REAL(8) :: x, pi = 4D0*atan(1D0)
	REAL(8) :: x2,rnu,zeta,Fermi

		x2=REAL(x,KIND(x2))
		IF((0.17*(x2+1)**2D0)<708D0) THEN
			rnu = x2**4D0 + 50D0 + 33.6*x2*(1D0-0.68*exp(-0.17*(x2+1)**2D0))
		ELSE
			rnu = x2**4D0 + 50D0
		ENDIF
		zeta = 3D0*SQRT(PI)/(4D0*rnu**0.375)
		Fermi = 1D0/(exp(-x) + zeta)
	RETURN
END FUNCTION Fermi
!*********************************************************************************

!*********************************************************************************
FUNCTION Ec_Ei(rn, rv, Eg)
      REAL*8 :: x,x1,x2,f1,f2,rtbis
	  REAL*8 :: Ec_Ei
	  DOUBLE PRECISION, EXTERNAL :: Fermi
	  INTEGER :: k
      LOGICAL I$FLAG
      REAL*8 FX(0:200)
      REAL*8 xMid, ddx
	  REAL*8 :: rn,rv,Eg,fmid

      x = -Eg
      ddx = 5.0
      fx(0) = rn*Fermi(-x) - rv*Fermi(x-Eg)
      I$FLAG = .FALSE.
      k = 0
      DO WHILE(.NOT.I$FLAG)
          k = k + 1
          IF(k > 100)THEN
              PRINT*,'Too small array!!!'
              PRINT*,'Function Ec_Ei'
              STOP
          ENDIF
          x = x + ddx
          fx(k) = rn*Fermi(-x) - rv*Fermi(x-Eg)
          IF(fx(k-1)*fx(k) < 0D0)THEN
              x1 = x - ddx
              x2 = x
              f1 = fx(k-1)
              f2 = fx(k)
              I$FLAG = .TRUE. 
          ENDIF
      ENDDO
      rtbis = x1
      I$FLAG = .FALSE.
      DO WHILE(.NOT.I$FLAG)
          ddx = 0.50*ddx
          xMid = rtbis + ddx
          fmid = rn*Fermi(-xMid) - rv*Fermi(xMid-Eg)
          IF(fmid > 0D0) rtbis = xMid
          IF((abs(ddx) < 1D-6) .OR. (fmid == 0D0)) I$FLAG = .TRUE.
      ENDDO
      Ec_Ei = xMid
   RETURN
END FUNCTION Ec_Ei
!*********************************************************************************

!*********************************************************************************
FUNCTION VOLTAGE(C, Ncnorm, Nvnorm, Eg, dEc)
!	called from init_potential -> "C" is at a given node
	IMPLICIT NONE
	LOGICAL I$FLAG
	REAL*8 :: FX(0:200)
	REAL*8 :: xMid, ddx
	REAL*8 :: Ncnorm, Nvnorm,C,Eg,dEc
	INTEGER :: k
	REAL*8 :: x,x1,x2,f1,f2,fmid,rtbis
	REAL*8 :: Voltage
	REAL*8, EXTERNAL :: Fermi

	x = -2D0*Eg

	ddx = 5D0
	fx(0) = Ncnorm*Fermi(x-dEc) - Nvnorm*Fermi(dEc-Eg-x) + C
	I$FLAG = .FALSE.
	k = 0
	DO WHILE(.NOT.I$FLAG)
	  k = k + 1
	  IF(k > 200) then
		PRINT*,'Too small array!!!'
		PRINT*,'Function Voltage'
		STOP
	  ENDIF
	  x = x + ddx
	  fx(k) = Ncnorm*Fermi(x-dEc) - Nvnorm*Fermi(dEc-Eg-x) + C
	  IF(fx(k-1)*fx(k) < 0)THEN
		 x1 = x - ddx
		 x2 = x
		 f1 = fx(k-1)
		 f2 = fx(k)
		 I$FLAG = .TRUE. 
	  ENDIF
	ENDDO
	rtbis = x1
	I$FLAG = .FALSE.
	DO WHILE(.NOT.I$FLAG)
	  ddx = 0.50*ddx
	  xMid = rtbis + ddx
	  fmid = Ncnorm*Fermi(xMid-dEc) - Nvnorm*Fermi(dEc-Eg-xMid) + C
	  IF(fx(0).lt.0)then
		IF(fmid < 0D0) rtbis = xMid
	  ELSEIF(fx(0) > 0D0)then
		IF(fmid > 0D0) rtbis = xMid
	  ENDIF
	  IF((abs(ddx) < 1D-12) .OR. (fmid == 0D0)) I$FLAG = .TRUE.
	ENDDO
	Voltage = xMid
	RETURN
END FUNCTION Voltage
