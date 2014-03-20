!Calculate coefficients for arrays of SIP algorithm
SUBROUTINE SIP_Coefficients
    USE SIP_Module
    IMPLICIT NONE

    REAL(8) :: epsPlusHalf, epsMinusHalf

    !Allocate coefficient arrays for SIP algorithm
    ALLOCATE(Bvec(1:Npx-1,1:Npy-1,1:Npz-1));  Bvec = 0D0
    ALLOCATE(Cvec(1:Npx-1,1:Npy-1,1:Npz-1));  Cvec = 0D0
    ALLOCATE(Dvec(1:Npx-1,1:Npy-1,1:Npz-1));  Dvec = 0D0
    ALLOCATE(Evec(1:Npx-1,1:Npy-1,1:Npz-1));  Evec = 0D0
    ALLOCATE(Fvec(1:Npx-1,1:Npy-1,1:Npz-1));  Fvec = 0D0
    ALLOCATE(Gvec(1:Npx-1,1:Npy-1,1:Npz-1));  Gvec = 0D0
    ALLOCATE(Hvec(1:Npx-1,1:Npy-1,1:Npz-1));  Hvec = 0D0
    
    !!!!Calculate coefficients B through H for GaAs layers
    DO ii = 1,Npx-1
        DO jj = 1,Npy-1
            DO kk = 1,Npz-1
                IF (ii == 1)	THEN
		            Dvec(ii,jj,kk)=0.0
		            Fvec(ii,jj,kk)=(2.0*epsMat(ii,jj,kk))/(hx*(hx+hx))
	            ELSEIF (ii == Npx-1)	THEN
		            Dvec(ii,jj,kk)=(2.0*epsMat(ii,jj,kk))/(hx*(hx+hx))
		            Fvec(ii,jj,kk)=0.0
	            ELSE		
		            epsPlusHalf=0.5*(epsMat(ii+1,jj,kk)+epsMat(ii,jj,kk))
		            epsMinusHalf=0.5*(epsMat(ii,jj,kk)+epsMat(ii-1,jj,kk))
		            Dvec(ii,jj,kk)=(2.0*epsMinusHalf)/(hx*(hx+hx))
		            Fvec(ii,jj,kk)=(2.0*epsPlusHalf)/(hx*(hx+hx))
	            ENDIF

	            IF (jj == 1)	THEN
		            Cvec(ii,jj,kk)=0.0
		            Gvec(ii,jj,kk)=(2.0*epsMat(ii,jj,kk))/(hy*(hy+hy))
	            ELSEIF (jj == Npy-1)	THEN
		            Cvec(ii,jj,kk)=(2.0*epsMat(ii,jj,kk))/(hy*(hy+hy))
		            Gvec(ii,jj,kk)=0.0
	            ELSE
		            epsPlusHalf=0.5*(epsMat(ii,jj+1,kk)+epsMat(ii,jj,kk))
		            epsMinusHalf=0.5*(epsMat(ii,jj,kk)+epsMat(ii,jj-1,kk))
		            Cvec(ii,jj,kk)=(2.0*epsMinusHalf)/(hy*(hy+hy)) 
		            Gvec(ii,jj,kk)=(2.0*epsPlusHalf)/(hy*(hy+hy))
	            ENDIF

	            IF (kk == 1)	THEN
		            Bvec(ii,jj,kk)=0.0
		            Hvec(ii,jj,kk)=(2.0*epsMat(ii,jj,kk))/(hz*(hz+hz))
	            ELSEIF (kk == Npz-1)	THEN
		            Bvec(ii,jj,kk)=(2.0*epsMat(ii,jj,kk))/(hz*(hz+hz)) 
		            Hvec(ii,jj,kk)=0.0
	            ELSE
		            epsPlusHalf=0.5*(epsMat(ii,jj,kk+1)+epsMat(ii,jj,kk))
		            epsMinusHalf=0.5*(epsMat(ii,jj,kk)+epsMat(ii,jj,kk-1))
		            Bvec(ii,jj,kk)=(2.0*epsPlusHalf)/(hz*(hz+hz))
		            Hvec(ii,jj,kk)=(2.0*epsMinusHalf)/(hz*(hz+hz))
	            ENDIF

	            IF((materialType(ii,jj,kk)>=3) .AND. (materialType(ii,jj,kk)<=4)) THEN  !If contact/gate
		            Cvec(ii,jj,kk)=0.0
		            Gvec(ii,jj,kk)=0.0
		            Bvec(ii,jj,kk)=0.0
		            Hvec(ii,jj,kk)=0.0
		            Dvec(ii,jj,kk)=0.0
		            Fvec(ii,jj,kk)=0.0

		            Evec(ii,jj,kk) = 1.0
	            ELSE
		            Evec(ii,jj,kk) = - (Cvec(ii,jj,kk)+Gvec(ii,jj,kk)+Bvec(ii,jj,kk)+Hvec(ii,jj,kk)+Dvec(ii,jj,kk)+Fvec(ii,jj,kk))
	            ENDIF
            ENDDO
        ENDDO
    ENDDO
END SUBROUTINE SIP_Coefficients