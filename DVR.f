C-------------------------------------------------------------------

C-------------------------------------------------------------------

      PROGRAM DVR


C

C This program uses the DVR described by D.T.Colber and W.H.Miller

C in J. Chem. Phys. 96, 1982 (1992).

C

C Writen by J.J.Soares Neto - Feb./1995

C-------------------------------------------------------------------

C Declare Variables.

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (NDIM=1000)

      DIMENSION W(3*NDIM*(NDIM+1)/2),IW(NDIM),EIGVCT(NDIM,NDIM),

     *          XVCT(NDIM,NDIM)
      call cpu_time ( t1 )
C-------------------------------------------------------------------

C Read in input data.

      CALL READ_I(A,B,NPOINT,AMASS,ELOW,EHIGH,NUMEIG)

C-------------------------------------------------------------------

C Build the Hamiltonian.

      CALL BUILD(A,B,NPOINT,AMASS,W,IW)

C-------------------------------------------------------------------

C Build the DVR From the Primitive Functions.

      CALL DVR_PR(A,B,NPOINT)

C-------------------------------------------------------------------

C Diagonilize the Hamiltonian.

      CALL EIGCLL(ELOW,EHIGH,NUMEIG,W,IW,NPOINT,NDIM,EIGVCT,

     *            NUMEV)



C-------------------------------------------------------------------

C Plot the Optimized Finite Basis Representation.

       CALL FBR(A,B,NPOINT,NDIM,EIGVCT)

C-------------------------------------------------------------------

C Generate the X Matrix.

      CALL OPTMZD(A,B,NDIM,EIGVCT,NPOINT,NUMEV,W,IW)

C-------------------------------------------------------------------

C Diagonilize the X Matrix and Obtain the Points of the

C Optimized Quadrature Rule.

!      CALL EIGCLL(A,B,NUMEV,W,IW,NUMEV,NDIM,XVCT,NUMEV)

C-------------------------------------------------------------------

C Build the DVR  From the Optimized Functions.

      CALL DVR_OP(A,B,NPOINT,NUMEV,W,NDIM,EIGVCT,XVCT)

C-------------------------------------------------------------------

CC      pause 'FIM'
      call cpu_time ( t2 )
      write ( *, * ) 'tempo ', t2 - t1

      READ *
      END





C-------------------------------------------------------------------
      SUBROUTINE BUILD(A,B,NPOINT,AMASS,W,IW)
C-------------------------------------------------------------------
C Declare Variables
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (PI=3.1415926535898D0)
      DIMENSION W(1),IW(1)
C-------------------------------------------------------------------
C Build the Hamiltonian
      FCT1=(PI**2)/(4.0D0*AMASS*(B-A)**2)
      FCT2=PI/(2.0D0*DBLE(NPOINT+1))
      FCT3=PI/DBLE(NPOINT+1)
      NCT=0
      DO 10 I=1,NPOINT
         DO 20 J=1,I-1
            NCT=NCT+1
            W(NCT)=(-1)**(I-J)*FCT1*(1.0D0/(SIN((I-J)*FCT2)**2)-
     *               1.0D0/(SIN((I+J)*FCT2)**2))
            W(NPOINT*(NPOINT+1)/2+NCT)=0.0D0
20       CONTINUE
         NCT=NCT+1
         XI=A+(B-A)*I/DBLE(NPOINT+1)
         W(NCT)=FCT1*((2.0D0*DBLE(NPOINT+1)**2+1.0D0)/3.0D0-
     *            1.0D0/DSIN(FCT3*I)**2)+V(XI)
         W(NPOINT*(NPOINT+1)/2+NCT)=1.0D0
         IW(I)=I*(I+1)/2
	
10    CONTINUE
C-------------------------------------------------------------------
      RETURN
      END
C-------------------------------------------------------------------

C-------------------------------------------------------------------
      SUBROUTINE DVR_OP(A,B,NPOINT,NUMEV,W,NDIM,EIGVCT,XVCT)
C-------------------------------------------------------------------
C Declare Variables.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EIGVCT(NDIM,NDIM),XVCT(NDIM,NDIM)
C-------------------------------------------------------------------
      PRINT*,'OPTIMIZED DVR'
      DO 10 K=1,1
         DO 20 I=1,NPOINT
            FXI=0.0D0
            DO 30 L=1,NUMEV
               FXI=FXI+EIGVCT(L,I)*XVCT(K,L)
30          CONTINUE
            XI=A+(B-A)*I/DBLE(NPOINT+1)

		  TMPkcm = dble(I)
            WRITE(11,*) XI,V(TMPkcm)
20       CONTINUE
10    CONTINUE
C-------------------------------------------------------------------
      RETURN
      END
C-------------------------------------------------------------------


C----------------------------------------------------------------

C-------------------------------------------------------------------
      SUBROUTINE DVR_PR(A,B,NPOINT)
C-------------------------------------------------------------------
C Declare Variables.
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (PI=3.1415926535898D0)
C-------------------------------------------------------------------
C Build the DVR Function.
      WRITE(10,*)'# PRIMITIVE DVR'
      DXPT=(B-A)/DBLE(NPOINT+1)
      XPT=A
      XI=A+1*(B-A)/(DBLE(NPOINT+1))
C While Loop.
100   CONTINUE
      IF(XPT.LT.B) THEN
         FXI=0.0D0
         DO 10 N=1,NPOINT
            XN=DBLE(N)
            FXI=FXI+DSIN(XN*PI*(XPT-A)/(B-A))*DSIN(XN*PI*(XI-A)/(B-A))
10       CONTINUE
         FXI=2.0D0*FXI/(DBLE(NPOINT+1))
c        FXI=2.0D0*FXI/(B-A)
         WRITE(10,*) XPT,FXI
         XPT=XPT+DXPT
      GOTO 100
      ENDIF
C-------------------------------------------------------------------
      RETURN
      END
C-------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE EIGCLL(ELOW,EHIGH,NUMEIG,W,IW,NPOINT,NDIM,EIGVCT,
     *                  NUMEV)
C
C   This routine sets the input and print out the output
C   of STLM.
C
C   Written by J. J. Soares Neto, July, 1991.
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C   Define the variables and arrays.
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER DAFILE,ERRNO,TCONV,X,TON,PMAX,PROFIL
      LOGICAL DIAGM
      DIMENSION W(6),IW(1),EIGVCT(NDIM,NDIM)
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C  Debugging part. If IEIGE=1 print out the Hamiltonian and the Over-
C  lap matrices.
      IEIGE=0
      IF(IEIGE.EQ.1) THEN
         NTOTAL=IW(NPOINT)
         DO 1 N=1,NtotaL
            PRINT*,N,W(N),W(N+NTOTAL)
1        CONTINUE
      ENDIF

C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C  Sets the input parameters for STLM. Read STLM User's Guide for de-
C  tails.
      A=ELOW
      B=EHIGH
      DIAGM=.FALSE.
      MAXW=820000
      MAXIW=16000
      MAXL=80
      MAXREC=550
      PROFIL=1
      PMAX=50
      MXREST=0
      MSGLVL=2
      MAXL=NUMEIG+NUMEIG
      DAFILE=19
      KFILE=18
      TON=8*NPOINT
      ISIZE=(NPOINT*MAXREC)/64 + 1
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C  Open Scratch files for STLM.
      CLOSE(DAFILE)
      OPEN(UNIT=DAFILE,ERR=271,STATUS='SCRATCH',ACCESS='DIRECT',
     *     RECL=TON)
      GOTO 272
271   STOP 'error opening dafile in EIGCAQ'
272   CONTINUE
      CLOSE(KFILE)
      OPEN(UNIT=KFILE,STATUS='SCRATCH',FORM='UNFORMATTED')
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C   Calculates Eigenvalues and Eigevectors using Lanczos procedure.
      CALL STLM(NPOINT,A,B,MAXL,PROFIL,PMAX,MXREST,MSGLVL,MAXW,
     *          MAXIW,DAFILE,MAXREC,KFILE,X,BG,TCONV,NLEFT,
     *          ERRNO,W,IW)
C-----------------------------------------------------------------------



C
C-----------------------------------------------------------------------
      NUMEV=MIN(NUMEIG,TCONV)
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
     

*	WRITE(4,*) 'ENERGIA ROVIBRACIONAL J=1(cm-1)'

	WRITE(4,*) 'ENERGIA VIBRACIONAL (cm-1)'

	WRITE(4,*) '----------------------------------------------'

      DO I=1,NUMEV

	   

	WRITE(4,*) I,219474.631*W(I)

	WRITE(4,*) '----------------------------------------------'

	END DO



	WRITE(4,*) '                                             '

	WRITE(4,*) '                                             '

*	WRITE(4,*) 'ESPECTRO ROVIBRACIONAL J=1 (cm-1)'

	WRITE(4,*) 'ESPECTRO VIBRACIONAL (cm-1)'
        WRITE(4,*) '----------------------------------------------'


      

	DO 30 I=1,NUMEV

	   WRITE(4,*) I,219474.631*(W(I+1)-W(1))

         
200      FORMAT(22H     Eigenvalue number   ,I3,3H = ,f25.9,
     *          9H hartree.)
         WRITE(4,*) '----------------------------------------------'
         

C------
C------
C Read the Eigenvectors from file 19 and write to 20.
!         CALL EXCHG(NDIM,EIGVCT,NPOINT,I)
C------
C------
30    CONTINUE  
      WRITE(4,*) '                                             '

	WRITE(4,*) '                                             '

*	WRITE(4,*) 'DIFERENÇA ENTRE OS NIVEIS ROVIBRACIONAIS J=1 (cm-1)'

      WRITE(4,*) 'DIFERENCA ENTRE OS NIVEIS VIBRACIONAIS (cm-1)'

	WRITE(4,*) '----------------------------------------------'

	DO I=1,NUMEV

	

         WRITE(4,*) I,219474.631*(W(I+1)-W(I))

	   WRITE(4,*) '----------------------------------------------'

	

	END DO

C     #aqui1	

      we=29.434588840876046
    
      wexe=1.99859909076558900E-002

      weye= 2.14735624358330610E-005
      alfae=1.0d0/8.0d0*(-12.0d0*(W(2)-W(1))*219474.631D0+4.0d0*

     *(W(3)-W(1))*219474.631D0+4.0d0*we-23.0d0*weye)

      gamae=1.0d0/4.0d0*(-2.0d0*(W(2)-W(1))*219474.631D0+(W(3)-W(1))*

     *219474.631D0+2.0d0*wexe-9.0d0*weye)


      we=1.0d0/24.0d0*(141.0d0*(W(2)-W(1))-93.0d0*(W(3)-W(1))+

     *23.0d0*(W(4)-W(1)))

      wexe=1.0d0/4.0d0*(13.0d0*(W(2)-W(1))-11.0d0*(W(3)-W(1))+

     *3.0d0*(W(4)-W(1)))

      weye=1.0d0/6.0d0*(3.0d0*(W(2)-W(1))-3.0d0*(W(3)-W(1))+(W(4)-W(1)))

      WRITE(97,*)we*219474.631D0

      

         WRITE(4,*) '                                              '

	   WRITE(4,*) '                                              '

	   WRITE(4,*) 'CONSTANTES ESPECTROSCÓPICAS'

	   WRITE(4,*) '----------------------------------------------'

	   WRITE(4,*) 'WE(cm-1)=',we*219474.631D0   !! s¢ imprimir essa

	   WRITE(4,*) 'ALFAE(cm-1)=',alfae

         WRITE(4,*) '----------------------------------------------'

	   WRITE(4,*) 'WEXE(cm-1)=',wexe*219474.631D0

	   WRITE(4,*) 'GAMAE(cm-1)=',gamae

	   WRITE(4,*) '----------------------------------------------'

	   WRITE(4,*) 'WEYE(cm-1)=',weye*219474.631D0

	   WRITE(4,*) '----------------------------------------------'

C-----------------------------------------------------------------------
C   Close scratch files used by STLM.
      CLOSE(UNIT=DAFILE)
      CLOSE(UNIT=KFILE)
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      RETURN
      END
C-----------------------------------------------------------------------

C--------------------------------------------------------------------
      SUBROUTINE EXCHG(NDIM,EIGVCT,NPOINT,I)
C--------------------------------------------------------------------
C Declare Variables.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EIGVCT(NDIM,NDIM)
C--------------------------------------------------------------------
C   Read the vectors which diagonilizes the Model Hamiltonian.
      READ(UNIT=19,REC=I) (EIGVCT(I,L),L=1,NPOINT)
C--------------------------------------------------------------------
      RETURN
      END
C--------------------------------------------------------------------

C-------------------------------------------------------------------
      SUBROUTINE FBR(A,B,NPOINT,NDIM,EIGVCT)
C-------------------------------------------------------------------
C Declare Variables.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EIGVCT(NDIM,NDIM)
C-------------------------------------------------------------------
      PRINT*,'FBR'
      DO 10 I=1,NPOINT
         XI=A+(B-A)*I/DBLE(NPOINT+1)
c       PRINT*,XI,EIGVCT(1,I)
10    CONTINUE
C-------------------------------------------------------------------
      RETURN
      END
C-------------------------------------------------------------------

C---------------------------------------------------------------------
      SUBROUTINE OPTMZD(A,B,NDIM,EIGVCT,NPOINT,NUMEV,W,IW)
C---------------------------------------------------------------------
C Declare Variables.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EIGVCT(NDIM,NDIM),W(1),IW(1)
C---------------------------------------------------------------------
C Build the Matrix X.
      NCT=0
      DO 10 I=1,NUMEV
         DO 20 J=1,I
            NCT=NCT+1
            W(NCT)=0.0D0
            DO 30 K=1,NPOINT
               XI=A+(B-A)*K/DBLE(NPOINT+1)
               W(NCT)=W(NCT)+EIGVCT(I,K)*XI*EIGVCT(J,K)	 
30    			CONTINUE
                  W(NUMEV*(NUMEV+1)/2+NCT)=0.0D0
            		  IF(I.EQ.J) W(NUMEV*(NUMEV+1)/2+NCT)=1.0D0
20          CONTINUE
          
	   IW(I)=I*(I+1)/2
	
10    CONTINUE
      
C---------------------------------------------------------------------
      RETURN
      END
C---------------------------------------------------------------------      

	REAL*8 FUNCTION RANNUM(IDUM)
         IMPLICIT REAL*8 (a-h,o-z)
         implicit real*8 (m)
         PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=2.5E-7)
c     PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RANNUM=MJ*FAC
      RETURN
      END




	real*8 function ranfx(dummy)
c
	integer dummy
	real*8 avrge, drand, rand
	logical first
c	external rand
	save first
	data first/.true./,avrge/1.5/
c
	if (first) then
        ranfx = rand(3)
	    first = .false.
	endif
c
	ranfx = rand(0) + rand(0) + rand(0) - avrge
	return
	end
c
	real*8 function second(dummy)
c
c	Suggested by Hans Jxrgen Aagaard Jensen 880126
c
	real*4 etime, tarr(2)
	integer dummy
		second = etime(tarr)
	return
	end
cC----------------------------------------------------------------
      SUBROUTINE READ_I(A,B,NPOINT,AMASS,ELOW,EHIGH,NUMEIG)
C----------------------------------------------------------------
C Declare Variables.
      IMPLICIT REAL*8 (A-H,O-Z)
C----------------------------------------------------------------
      READ(3,*) A,B,NPOINT,AMASS
      READ(3,*) ELOW,EHIGH,NUMEIG
C----------------------------------------------------------------
      RETURN
      END
C----------------------------------------------------------------CALLOC
      LOGICAL FUNCTION ALLOC(ID1, ID2, AD1, AD2, ACTION, W, ADDRSS)     L***
C
C **********************************************************************
C
C     PURPOSE - (VER = 1 OR 2)
C
C         THIS ROUTINE ESTABLISHES CONNECTIONS BETWEEN IDENTIFIERS
C         AND VECTORS OF LENGTH N IN W.
C         IF VEC(ID2) DOES NOT LIE IN W ALREADY IT IS
C         FETCHED FROM SECONDARY STORAGE WITH THE HELP OF IO. IF AD IS
C         THE ADDRESS IN W THEN
C         VEC(ID)=(W(AD), ..., W(AD+N-1)), FOR ID = ID1, ID2, AND
C         AD = AD1, AD2.
C
C     INPUT PARAMETERS -
C
C         ID1, ID2   IDENTIFIERS OF THE VECTORS
C         ACTION     THE CORRESPONDING ACTION
C         ADDRSS     ADDRESS VECTOR
C
C
C     OUTPUT PARAMETERS -
C
C         AD1, AD2   THE ADDRESSES CORRESPONDING TO ID1 AND ID2.
C
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          RDUMP, SECOND, ST, TIME, W(1)                           R MOD.
      INTEGER   ACTION, ACTIVE, AD(2), AD1, AD2, ADDRSS(1), COUNT,      I***
     +          DUMMY, ERRNO, FREE, FREPOS, HIT, I, ID, ID1, ID2, IDUMP,
     +          MV, MXNEW, MXOLD, MXRST, N, NBADMU, NIL, NMXRST, NOACTN,
     +          NOR, READID, READK, SAEVAL, SAVE, SAVFRE, SCPX, SOLCPX,
     +          V, WAD, WRITID, X
      LOGICAL   F, IO, T                                                L***
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMWH/ WAD(2), ACTIVE(2)
C
      DATA      NOR  / 1 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      ALLOC = T
C
      AD(1) = NIL
      AD(2) = NIL
      ID = ID1
C     *********************
C     LOOP FOR ID1 AND ID2.
C     *********************
      I = 1
C
10    FREPOS = 0
      IF(ACTIVE(1) .EQ. NIL) FREPOS = 1
      IF(ACTIVE(2) .EQ. NIL) FREPOS = 2
C
C     ******************************
C     CHECK IF ID ALREADY LIES IN W.
C     ******************************
      HIT = 0
      IF(ID .EQ. ACTIVE(1)) HIT = 1
      IF(ID .EQ. ACTIVE(2)) HIT = 2
C
      IF(.NOT. ADDRSS(ID) .LT. 0) GOTO 15
        AD(I) = -ADDRSS(ID)
        GOTO 40
C
15    IF(.NOT. HIT .NE. 0) GOTO 20
        AD(I) = WAD(HIT)
        GOTO 40
C
20      IF(.NOT. FREPOS .LE. 0) GOTO 30
C         ******************************************
C         NO ROOM LEFT. SHOULD NEVER HAPPEN, BUT ...
C         ******************************************
C
          IDUMP(1) = ID
          IDUMP(2) = ACTION
          CALL ERROR(NOR, 1)
          GOTO 8888
C
C     *******
C     UPDATE.
C     *******
30    ACTIVE(FREPOS) = ID
      AD(I) = WAD(FREPOS)
      IF(.NOT. I .EQ. 2) GOTO 40
C       ************
C       GET VEC(ID).
C       ************
        AD2 = AD(2)
        IF(.NOT. IO(W(AD2), ID, READID, N, ADDRSS)) GOTO 8888
C
40    IF(I .EQ. 2) GOTO 50
      I = 2
      ID = ID2
      IF(ID .EQ. NIL) GOTO 50
      GOTO 10
C
50    IF(ACTION .EQ. NOACTN .OR. ACTION .EQ. SAVE) GOTO 9999
C
C       **********
C       DELETE ID.
C       **********
        IF(ACTION .EQ. SAVFRE) ID = ID1
        IF(ID .EQ. ACTIVE(1)) ACTIVE(1) = NIL
        IF(ID .EQ. ACTIVE(2)) ACTIVE(2) = NIL
C
      GOTO 9999
C
8888  ALLOC = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      AD1 = AD(1)
      AD2 = AD(2)
C
      RETURN
      END
CALWAYS
      LOGICAL FUNCTION ALWAYS(W, IW)                                    
C
C **********************************************************************
C
C     PURPOSE -
C
C         HANDLES THE LANCZOS RUN FOR ONE SHIFT. COMPUTES EIGENVALUES
C         ON BOTH SIDES OF THE SHIFT.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  
     +          CK, CL, COEFF, RDUMP, SECOND, ST, TIME, TIMTQL, TLDL,   
     +          TNORM, TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1), WRR
      INTEGER   ALPHA, BETA, BETAPI, CNEG, CONV, COPT, COUNT, CPOS,     
     +          DUMMY, ERRNO, IDUMP, ITERNO, IW(1), LAMBDA, MV, MXNEG,
     +          MXNEW, MXOLD, MXREST, MXRST, N, NBADMU, NIL, NMXRST,
     +          NOR, NU, OLCPOS, P, PFCONV, PMAX, POINTR, POPT, REST,
     +          RNEW, ROLD, S, SCPX, SCR, SOLCPX, TCONV, V, WRI, X
      LOGICAL   CONTIN, CPALSO, F, FINAL, INITLA, LANCZO, MEQI, POSDOT, 
     +          POSNU, PREPSV, RANDVC, SAVEXL, T, TRIDIG, UPDATE, USEMX,
     +          VISITL, WRL, ZERBET, REACHB, USSMXR, USEDB
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
      COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMTS/ REACHB, USSMXR, USEDB
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
C
      DATA      NOR  / 2 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      ALWAYS = T
C
      VISITL = F
C     ******************************************************
C     GENERATE A RANDOMVECTOR AND ORTHOGONALIZE IT AGAINST
C     THE OLCPOS-VECTORS, PRODUCING A STARTINGVECTOR (LET US
C     CALL IT SV).
C     ******************************************************
10    IF(.NOT. RANDVC(V + 1, W, IW, .NOT. USEMX .AND. .NOT. MEQI))
     +                                             GOTO 8888
C
      IF(.NOT. PREPSV(SOLCPX, MXOLD, OLCPOS, W, IW)) GOTO 8888
      IF(.NOT. INITLA(CONTIN, FINAL, POSDOT, POSNU, W, IW)) GOTO 8888
C     ***********************
C     IS SV(T)*M*SV POSITIVE.
C     ***********************
      IF(POSDOT) GOTO 30
        IF(.NOT. VISITL) GOTO 20
          CALL ERROR(NOR, 1)
          GOTO 8888
C
20      VISITL = T
        GOTO 10
C
C     *******************************
C     RUN LANCZOS. CHECK CONVERGENCE.
C     *******************************
30    IF(.NOT. (CONTIN .AND. (.NOT. POSNU .OR.
     +         (CNEG .EQ. 0 .AND. REST .GT. 0)))) GOTO 40
C
        IF(USEDB .AND. CNEG .GE. REST) GOTO 50
C
        IF(.NOT. LANCZO(W(ALPHA), W(BETA), TNORM, W, IW)) GOTO 8888
        IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, T, W, IW)) GOTO 8888
        GOTO 30
C
40    IF(.NOT. (CONTIN .AND. P .LT. POPT .AND.
     +         (CNEG .LT. REST .OR. CNEG + CPOS .LT. COPT))) GOTO 45
C
        IF(USEDB .AND. CNEG .GE. REST) GOTO 50
C
        IF(.NOT. LANCZO(W(ALPHA), W(BETA), TNORM, W, IW)) GOTO 8888
        IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, T, W, IW)) GOTO 8888
        GOTO 40
C
45    IF(.NOT. (CONTIN .AND. CNEG + CPOS .EQ. 0)) GOTO 50
C
        IF(USEDB .AND. CNEG .GE. REST) GOTO 50
C
        IF(.NOT. LANCZO(W(ALPHA), W(BETA), TNORM, W, IW)) GOTO 8888
        IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, T, W, IW)) GOTO 8888
        GOTO 45
C
C     **********************
C     NO MORE LANCZOS STEPS.
C     **********************
50    FINAL = T
      IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, T, W, IW)) GOTO 8888
      CALL FREEID(MV + P + 1)
C
      IF(ITERNO .GT. 1) REST = REST - CNEG
      WRR(1) = TNORM
      CALL WRINFO(NOR, 1, W, IW)
      IF(.NOT. REST .LT. 0) GOTO 60
C       ****************************
C       SHOULD NEVER HAPPEN, BUT ...
C       ****************************
        CALL ERROR(NOR, 2)
        GOTO 8888
C
60    IF(.NOT. REST .EQ. 0) GOTO 80
C       ***************************
C       FINISHED WITH THE INTERVAL.
C       ***************************
        IF(.NOT. POSNU) GOTO 70
          CPALSO = T
          MXNEG = NIL
          SCPX = X + TCONV + CNEG
C         *****************************
C         COMPUTE AND STORE EIGENPAIRS.
C         *****************************
          IF(.NOT. SAVEXL(W(BETAPI), W(LAMBDA), W(NU), IW(POINTR), W(S),
     +                    MXNEG, PMAX, CPALSO, W, IW)) GOTO 8888
C         ******************
C         NEW POPT AND COPT.
C         ******************
          CALL NEWPC(W, IW)
          GOTO 9999
C
C       *****************************************
C       WE DO NOT HAVE A NEW SHIFT, BUT REST = 0.
C       COMPUTE AND STORE EIGENPAIRS.
C       *****************************************
70      CPALSO = F
        MXNEG = NIL
        IF(.NOT. SAVEXL(W(BETAPI), W(LAMBDA), W(NU), IW(POINTR), W(S),
     +                  MXNEG, PMAX, CPALSO, W, IW)) GOTO 8888
        CALL ERROR(NOR, 3)
        GOTO 8888
C
C     *************************************************
C     WE DO NOT HAVE A NEW SHIFT, AND REST IS POSITIVE.
C     *************************************************
80    IF(POSNU) GOTO 9999
        CALL ERROR(NOR, 4)
C
8888  ALWAYS = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CBLOCK
      BLOCK DATA
C
C **********************************************************************
C
C    THE IMPLEMENTOR SHOULD INITIALIZE THE VARIABLES DESCRIBED BELOW.
C
C
C
C    NERR   = STANDARD ERROR MESSAGE UNIT.
C             THE LOGICAL UNIT NUMBER FOR THE STANDARD ERROR
C             MESSAGE UNIT OF THE SYSTEM.
C
C    NOUT   = STANDARD OUTPUT UNIT.
C             THE LOGICAL UNIT NUMBER FOR THE STANDARD OUTPUT
C             UNIT OF THE SYSTEM.
C
C    SMALL  = 1.0/SRANGE, WHERE SRANGE IS THE LARGEST REAL NUMBER
C             X SUCH THAT THE ARITHMETIC OPERATIONS O ARE CORRECTLY
C             PERFORMED FOR ALL ELEMENTS A, B OF THE SYSTEM OF REAL
C             NUMBERS, PROVIDED THAT A, B AND THE EXACT MATHEMATICAL
C             RESULT OF A O B DO NOT HAVE AN ABSOLUTE VALUE OUTSIDE
C             THE RANGE (1/X, X).
C             (O STANDS FOR THE THE ARITHMETIC OPERATIONS, + - * /,
C             AND MONADIC -.)
C
C    VER    = VERSION OF STLM. PLEASE SEE THE INSTALLATION GUIDE.
C
C    SRELPR = RELATIVE PRECISION.
C             THE SMALLEST NUMBER X SUCH THAT
C             1.0-X .LT. 1.0 .LT. 1.0+X WHERE 1.0-X AND 1.0+X ARE
C             STORED VALUES OF THE COMPUTED RESULTS.
C
C    MORE INFORMATION ABOUT NERR, NOUT, SRANGE, AND SRELPR CAN BE
C    FOUND IN ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE, JUNE 1978,
C    VOLUME 4, NUMBER 2, 100-103.
C
C
C    NOTE. IF THE IMPLEMENTOR HAS USED THE TYPE CONVERTER TO CHANGE
C    THE TYPE REAL TO SOMETHING ELSE, THE CORRESPONDING CHANGES
C    SHOULD ALSO BE DONE TO SMALL AND SRELPR, SINCE THE TYPE CONVERTER
C    CAN NOT HANDLE THESE CHANGES.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  
     +          SMALL, SRELPR, RDUMMY                                   
      INTEGER   NERR, NOUT, VER, IDUMMY                                 
C
      COMMON   /STLMPR/ IDUMMY, NERR, NOUT
      COMMON   /STLMTL/ SRELPR, RDUMMY(9)
      COMMON   /STLMVR/ SMALL,  VER
C
      DATA              NERR,  NOUT   /6,       6/
      DATA              SMALL, VER    /1.0E-32, 1/
      DATA              SRELPR        /1.0E-14   /
C
      END
CCHECK
      SUBROUTINE CHECK(N, A, B, MAXL, PROFIL, PMAX, MXREST, MSGLVL,
     +                 MAXW, MAXIW, DAFILE, MAXREC, KFILE,
     +                 ERRNO, W, IW, LEN)
C
C **********************************************************************
C
C     PURPOSE -
C
C         THIS ROUTINE CHECKS INPUT.
C
C
C     OUTPUT PARAMETERS -
C
C         LEN =  0, NO ERROR HAS BEEN DETECTED.
C                1, ILLEGAL VER.
C                2, N LE 0.
C                3, ILLEGAL PROFIL.
C                4, MAXL LE 0.
C                5, MAXREC LE 0.
C                6, ERROR IN D (VER = 1 ONLY).
C                7, ERROR IN M (VER = 1 ONLY).
C                8, ILLEGAL MSGLVL.
C                9, ILLEGAL PMAX.
C               10, TOO SMALL MAXW.
C               11, TOO SMALL MAXIW.
C               12, B LE A.
C               13, ILLEGAL FILE UNIT(S).
C
C     PLEASE SEE THE USER GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  
     +          A, A1, B, B1, MU, RDUMP, SMALL, TOTALT, TPRED, W(1)     
      INTEGER   ADDRSS, CNEGF, COPT, CPOS, DAFIL1, DAFILE, ERRNO,       
     +          ERRNO1, I, IDUM, IDUMP, ITERNO, IW(1), KFILE, KFILE1,
     +          LEFTP, LEN, LENADR, LP, MAXIW, MAXIW1, MAXIWT, MAXL,
     +          MAXL1, MAXRE1, MAXREC, MAXW, MAXW1, MAXWT, MDIM,
     +          MSGLVL, MXREST, N, N1, N2, NREAD, NUMEIG,
     +          NUMVEC, NWRITE, PMAX, POPT, RFIRST, RIGHTC, RIGHTM,
     +          RIGHTP, RNEW, STADEW, TCONV, NERR  , NOUT , VER, X,
     +          ERRNO2, PROFI1, PROFIL, MXRES1, PMAX1, MSGLV1
      LOGICAL   LDUM, DIAGM, MEQI                                       
      COMMON   /STLMAD/ ADDRSS, DAFIL1, KFILE1, LP, MAXL1, NREAD, NWRITE
      COMMON   /STLMCT/ N2, IDUM(9), LDUM(2)
      COMMON   /STLMER/ RDUMP, ERRNO1, IDUMP(2)
      COMMON   /STLMEW/ LEFTP, LENADR, MAXRE1, NUMVEC, RIGHTC, RIGHTM,
     +                  RIGHTP, STADEW
      COMMON   /STLMIN/ A1, B1, NUMEIG, MAXW1, MAXIW1
      COMMON   /STLMPF/ PROFI1
      COMMON   /STLMPL/ PMAX1, POPT, COPT, MXRES1
      COMMON   /STLMPR/ MSGLV1, NERR  , NOUT
      COMMON   /STLMUI/ MU(3), TOTALT, TPRED, CNEGF, CPOS, ERRNO2,
     +                  ITERNO, N1, RFIRST, RNEW, TCONV, X
      COMMON   /STLMVR/ SMALL, VER
C
C
      DIAGM = PROFIL .GE. 2
      MEQI  = PROFIL .EQ. 3
C
      LEN = 1
      IF(VER .LT. 1 .OR. VER .GT. 3) GOTO 8888
C
      LEN = 2
      IF(N .LE. 0) GOTO 8888
C
      LEN = 3
      IF(PROFIL .LT. 1  .OR.  PROFIL .GT. 3) GOTO 8888
C
      LEN = 4
      IF(MAXL .LE. 0) GOTO 8888
C
      IF(.NOT. VER .LE. 2) GOTO 5
        LEN = 5
        IF(MAXREC .LE. 0) GOTO 8888
C
5     IF(.NOT. VER .EQ. 1) GOTO 10
        CALL DCHECK(IW, N, LEN, IDUMP(1))
        IF(LEN .EQ. 6) GOTO 8888
C
        I = IW(N) + 1
        CALL MCHECK(W(I), IW, N, DIAGM, MEQI, LEN, IDUMP(1))
        IF(LEN .EQ. 7) GOTO 8888
C
10    LEN = 8
      IF(MSGLVL .LT. 0 .OR. MSGLVL .GT. 4) GOTO 8888
C
      LEN = 9
      IF(PMAX .LT. 15) GOTO 8888
C
      IF(.NOT. VER .EQ. 1) GOTO 20
        MDIM = IW(N)
        IF(DIAGM) MDIM = N
        IF(MEQI)  MDIM = 0
        MAXWT = IW(N) + MDIM + PMAX * (5 + PMAX) + 2 * N + MAXL
        MAXIWT = N + 6 * PMAX + MAXL + 3
        IF(MEQI) MAXIWT = N + 2 * PMAX + MAXL + 2
        GOTO 40
C
20    IF(.NOT. VER .EQ. 2) GOTO 30
        MAXWT = PMAX * (5 + PMAX) + 2 * N + MAXL
        MAXIWT = 6 * PMAX + MAXL + 3
        IF(MEQI) MAXIWT = 2 * PMAX + MAXL + 2
        GOTO 40
C
30    MAXWT = PMAX * (5 + PMAX) + MAXL
      MAXIWT = PMAX + 1 + MAXL
C
40    LEN = 10
      IF(MAXW .LT. MAXWT) GOTO 8888
C
      LEN = 11
      IF(MAXIW .LT. MAXIWT) GOTO 8888
C
      LEN = 12
      IF(B .LE. A) GOTO 8888
C
      LEN = 13
      IF(.NOT. VER .LE. 2 ) GOTO 8000
        IF(DAFILE .EQ. NERR    .OR.  DAFILE .EQ. NOUT ) GOTO 8888
        IF(VER .EQ. 1  .AND.  (DAFILE .EQ. KFILE .OR.
     +     KFILE .EQ. NERR    .OR.  KFILE .EQ. NOUT ) ) GOTO 8888
C
C
8000  LEN = 0
      GOTO 9999
C
C     *************************************************
C     ERROR EXIT.
C     FILL STLMUI AND THE COMMON BLOCKS USED IN WRINFO.
C     *************************************************
8888  MU(1) = A
      MU(2) = A
      MU(3) = A
      TOTALT = 0.0D0                                                    
      TPRED = 0.0D0                                                     
      CNEGF = 0
      CPOS = 0
      ITERNO = 0
      N1 = N
      RFIRST = 0
      RNEW = 0
      TCONV = 0
      X = 0
C
      N2 = N
      A1 = A
      B1 = B
      MAXW1 = MAXW
      MAXIW1 = MAXIW
      MAXL1 = MAXL
      KFILE1 = KFILE
      DAFIL1 = DAFILE
      MAXRE1 = MAXREC
      PROFI1 = PROFIL
      PMAX1 = PMAX
      MXRES1 = MXREST
      MSGLV1 = 0
      IF(MSGLVL .EQ. 1) MSGLV1 = 1
      IF(MSGLVL .GT. 1 .AND. MSGLVL .LE. 4) MSGLV1 = 2
      CALL WRINFO(20, 1, W, IW)
C
      CALL ERROR(25, LEN)
      CALL ERROR(25, 25)
      CALL ERROR(20, 20)
      ERRNO = ERRNO1
      ERRNO2 = ERRNO1
C
C
9999  RETURN
      END
CCMPRSS
      SUBROUTINE CMPRSS(POINTR, PMAX)
C
C **********************************************************************
C
C     PURPOSE -
C
C         UPDATES POINTR, I.E. DELETES ELEMENTS POINTING TO FALSE
C         EIGENPAIRS.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      INTEGER   I, ISAVE, J, LEN, NHITS, PMAX, POINTR(1)                
C
      LEN = POINTR(PMAX + 1)
      IF(LEN .LE. 0) GOTO 9999
C     ********************************************
C     A FALSE ELEMENT IS .LE. 0. LOCATE THE FIRST.
C     ********************************************
      DO 10 I = 1, LEN
        IF(POINTR(I) .LE. 0) GOTO 20
10      CONTINUE
      GOTO 9999
C
20    NHITS = 0
      ISAVE = I
C     ****************************
C     WRITE OVER THE FALSE VALUES.
C     ****************************
      DO 40 J = ISAVE, LEN
        IF(.NOT. POINTR(J) .LE. 0) GOTO 30
          NHITS = NHITS + 1
          GOTO 40
C
30      POINTR(I) = POINTR(J)
        I = I + 1
C
40      CONTINUE
C
C     ******************
C     UPDATE THE LENGTH.
C     ******************
      POINTR(PMAX + 1) = LEN - NHITS
C
9999  RETURN
      END
CCOMPL
      SUBROUTINE COMPL(BETA, BETAPI, NU, LAMBDA, S, MAXNU, MINNU, PMAX,
     +                 LGEB)
C
C **********************************************************************
C
C     PURPOSE -
C
C         COMPUTES MAX AND MIN(NU(I)), I=1, ..., P, LAMBDA, AND BETAPI
C         FOR ALL NOT TOO SMALL, NU-VALUES.
C
C
C     OUTPUT PARAMETERS -
C
C         LGEB  = TRUE IF A CONVERGED EIGENVALUE IS .GE. B (SUCH
C                 EIGENVALUES ARE THROWN AWAY), AND FALSE OTHERWISE.
C                 IF LGEB = T, IT WILL CAUSE THE NEXT SHIFT TO BECOME B.
C         MAXNU = MAX(NU(1), ..., NU(P))
C         MINNU = MIN(NU(1), ..., NU(P))
C
C         IF ABS(NU(I)) .LE. TOLZNU THE VALUE IS NOT USED.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      INTEGER   PMAX                                                    I***
      DOUBLE PRECISION                                                  R INS.
     +          ALTMU, B, BETA(1), BETAPI(1), FACTOR, LAMBDA(1), SRELPR,R MOD.
     +          MAXNU, MINNU, MU, NEXTMU, NU(1), NUI, OLDMU, S(PMAX, 1),
     +          TOLBPI, TOLLDL, TOLPDM, TOLS1I, TOLZBT, TOLZNU,
     +          A, BB, CUT
      INTEGER   CNEG, CPOS, I, ITERNO, N, OLCPOS, P, REST, RNEW, ROLD,  I***
     +          TCONV, NUMEIG, MAXW, MAXIW
      LOGICAL   USEMX, ZERBET, REACHB, USSMXR, USEDB, T, F, LGEB        L***
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMIN/ A, B, NUMEIG, MAXW, MAXIW
      COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
      COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,
     +                  TOLZNU, TOLLDL(3)
      COMMON   /STLMTF/ T, F
      COMMON   /STLMTS/ REACHB, USSMXR, USEDB
C
C     ********************************************************
C     LOOK FOR A NU-VALUE WHERE 1/NU WILL BE SAFE TO COMPUTE.
C     INITIALIZE MINNU AND MAXNU.
C     ********************************************************
      DO 10 I = 1, P
        IF(.NOT. DABS(NU(I)) .GT. TOLZNU) GOTO 10                      
          MINNU = NU(I)
          MAXNU = MINNU
          GOTO 20
10      CONTINUE
C
      MINNU = -1.0D0                                             
      MAXNU = -1.0D0                                             
C
20    BB = BETA(P)
C     ****************************************************************
C     THE ONLY CASE WHEN CUT .NE. B, IS WHEN THE LDL(T) FACTORIZATION
C     FAILED FOR MU = B.
C     ****************************************************************
      CUT = B
      LGEB = F
      IF(REACHB) CUT = MU
      DO 40 I = 1, P
        NUI = NU(I)
        IF(.NOT. DABS(NUI) .GT. TOLZNU) GOTO 30                        
C         **********************************************************
C         FOR ALL NU-VALUES NOT NEAR ZERO, COMPUTE THE CORRESPONDING
C         LAMBDA AND CONVERGENCE VALUES.
C         UPDATE MINNU AND MAXNU.
C         **********************************************************
          LAMBDA(I) = MU + 1.0D0 / NUI                                 
          BETAPI(I) = BB * S(P, I) / NUI
          IF(.NOT. (LAMBDA(I) .GT. CUT  .AND.
     +             DABS(BETAPI(I)) .LE. TOLBPI)  ) GOTO 25             
C           **********************************************************
C           THROW AWAY CONVERGED EIGENVALUES GREATER THAN CUT, AND SET
C           LGEB.
C           **********************************************************
            BETAPI(I) = -1.0D20                                         
            LGEB      = T
C
25        MINNU = DMIN1(MINNU, NUI)                                     
          MAXNU = DMAX1(MAXNU, NUI)                                     
          GOTO 40
C
C       **************************************************
C       PUT IN VALUES SO THAT WRINFO WILL WORK.
C       BETAPI=-1E20 MEANS THAT IT WILL NEVER BE ACCEPTED.
C       **************************************************
30      LAMBDA(I) = -1.0D20                                            
        BETAPI(I) = -1.0D20                                            
C
40      CONTINUE
C
      RETURN
      END
CCONVER
      SUBROUTINE CONVER(BETAPI, LAMBDA, POINTR, S, PMAX)
C
C **********************************************************************
C
C     PURPOSE -
C
C         COMPUTE THE NUMBER OF EIGENVALUES IN THE INTERVAL
C         (OLDMU,MU) AND (OLDMU = -INFINITY FOR THE FIRST SHIFT) AND
C         THE NUMBER OF CONVERGED EIGENVALUES IN (MU,INFINITY).
C         COMPUTE CNEG AND CPOS, INITIALIZE POINTR.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      INTEGER   PMAX                                                    
      DOUBLE PRECISION                                                  
     +          ALTMU, BETAPI(1), CK, CL, COEFF, FACTOR, LAMBDA(1),     
     +          SRELPR, MU, NEXTMU, OLDMU, S(PMAX, 1), TIMTQL, TLDL,
     +          TOLBPI, TOLLDL, TOLPDM, TOLS1I, TOLZBT, TOLZNU, TOPINV,
     +          TOPM, TPRED, TSAVE, TVECOP
      INTEGER   CNEG, CNEGF, CONV, CPOS, I, ITERNO, LENP, N, OLCPOS, P, 
     +          PFCONV, POINTR(1), REST, RFIRST, RNEW, ROLD, TCONV
      LOGICAL   SAFRST, UPDATE, USEMX, ZERBET                           
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
      COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,
     +                  TOLZNU, TOLLDL(3)
C
      CNEG = 0
      CPOS = 0
      LENP = 0
      CONV = 0
C
      IF(.NOT. ITERNO .EQ. 1) GOTO 40
C       ***********************************************************
C       THE FIRST SHIFT IS SPECIAL, SINCE WE DO NOT HAVE ANY OLDMU.
C       ***********************************************************
        DO 30 I = 1, P
          IF(.NOT. DABS(BETAPI(I)) .LE. TOLBPI) GOTO 30                 
            CONV = CONV + 1
            IF(.NOT. MU .LT. LAMBDA(I)) GOTO 10
C             ************************************
C             CONVERGED TO THE RIGHT OF THE SHIFT.
C             ************************************
              CPOS = CPOS + 1
              GOTO 20
C
10          IF(.NOT. SAFRST) GOTO 30
C             ***********************************
C             CONVERGED TO THE LEFT OF THE SHIFT.
C             ***********************************
              CNEG = CNEG + 1
C
20          LENP = LENP + 1
C           ******************
C           INITIALIZE POINTR.
C           ******************
            POINTR(LENP) = I
C
30        CONTINUE
C       ***************************
C       LENGTH OF POINTR IS STORED.
C       ***************************
        POINTR(PMAX + 1) = LENP
        GOTO 9999
C
C     ********************************
C     THIS IS NOT THE FIRST ITERATION.
C     ********************************
40    DO 60 I = 1, P
        IF(.NOT. MU .LT. LAMBDA(I)) GOTO 50
          IF(.NOT. DABS(BETAPI(I)) .LE. TOLBPI) GOTO 60                 
C           *********
C           AS ABOVE.
C           *********
            CONV = CONV + 1
            CPOS = CPOS + 1
            LENP = LENP + 1
            POINTR(LENP) = I
            GOTO 60
C
50      IF(.NOT. (LAMBDA(I) .GT. OLDMU .AND. REST .GT. 0)) GOTO 60
C         **************************************************************
C         LAMBDA IS IN (OLDMU, MU).
C         IF REST = 0, NOTHING SHOULD CONVERGE, SO WE DO
C         NOT SAVE ANYTHING.
C         IF LAMBDA HAS NOT CONVERGED, POINTR IS GIVEN A NEGATIVE VALUE,
C         BUT CNEG IS NOT UPDATED IN THIS CASE.
C         **************************************************************
          LENP = LENP + 1
          POINTR(LENP) =  - I
          IF(.NOT. DABS(BETAPI(I)) .LE. TOLBPI) GOTO 60                 
            POINTR(LENP) = I
            CONV = CONV + 1
            CNEG = CNEG + 1
60      CONTINUE
C
      POINTR(PMAX + 1) = LENP
C     *******************************************************
C     IF CNEG =0, UPDATE POINTR, I.E. DELETE NEGATIVE VALUES.
C     *******************************************************
      IF(CNEG .EQ. 0) CALL CMPRSS(POINTR, PMAX)
C
9999  RETURN
      END
CDCHECK
      SUBROUTINE DCHECK(D, N, LEN, LOC)
C
C **********************************************************************
C
C     PURPOSE - (VER = 1)
C
C         THIS ROUTINE CHECKS THE COLUMN LENGTHS.
C
C     INPUT PARAMETERS -
C
C         D     = POINTER VECTOR TO DIAGONAL ELEMENTS OF K.
C         N     = DIMENSION OF M MATRIX.
C
C
C     OUTPUT PARAMETERS -
C
C         LEN = 0, NO ERROR HAS BEEN DETECTED.
C               6, A TOO SHORT, OR A TOO LONG COLUMN HAS BEEN DETECTED.
C         LOC = THE INDEX OF THE NUMBER OF THE COLUMN.
C
C **********************************************************************
C
      INTEGER   D(1), DI, DIM1, I, LEN, LENG, LOC, N                    
C
      DIM1 = 0
      DO 10 I = 1, N
        DI = D(I)
        LENG = DI - DIM1
        DIM1 = DI
C
        IF(.NOT. (LENG .LE. 0 .OR. LENG .GT. I)) GOTO 10
          LEN = 6
          LOC = I
          GOTO 8888
C
10      CONTINUE
C
8888  RETURN
      END
CDECOMP
      LOGICAL FUNCTION DECOMP(W, IW)                                    
C
C **********************************************************************
C
C     PURPOSE -
C
C         TO MAKE AN ACCEPTABLE LDL(T)-DECOMPOSITION OF K-MU*M. THE
C         DECOMPOSITION CAN BE REJECTED BECAUSE,
C
C         K-MU*M OR A PRINCIPAL LEADING SUBMATRIX (VER=1) OF IT IS
C         NEARLY SINGULAR (DECIDED IN LDLT, FLAG BADMU).
C
C         THE INTERVAL (OLDMU,MU) CONTAINS TOO MANY NOT YET COMPUTED
C         EIGENVALUES, I.E. REST .GT. MXREST.
C
C         SHOULD WE HAVE SUCH A MU AN OTHER IS CHOSEN AND WE TRY
C         AGAIN. WE MAKE ONLY ONE RETRY (THE VISIT FLAGS).
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  
     +          ALTMU, CK, CL, COEFF, MU, NEXTMU, OLDMU, RDUMP, SECOND, 
     +          ST, ST1, TIME, TIMTQL, TLDL, TMU, TOPINV, TOPM, TPRED,
     +          TSAVE, TVECOP, W(1), A, B
      INTEGER   CNEG, CONV, COPT, COUNT, CPOS, DUMMY, ERRNO, IDUMP,     
     +          ITERNO, IW(1), MXREST, N, NBADMU, NMXRST, NOR, OLCPOS,
     +          P, PFCONV, PMAX, POPT, REST, RNEW, ROLD, TCONV, TRNEW,
     +          NUMEIG, MAXW, MAXIW
      LOGICAL   BADMU, F, LDLT, T, UPDATE, USEMX, VISITD, VISITM,       
     +          ZERBET, REACHB, USSMXR, USEDB
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMIN/ A, B, NUMEIG, MAXW, MAXIW
      COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTS/ REACHB, USSMXR, USEDB
      COMMON   /STLMTF/ T, F
C
      DATA      NOR  / 3 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      DECOMP = T
C
C     ****************
C     SET VISIT FLAGS.
C     ****************
      VISITD = F
      VISITM = F
C
10    ST1 = SECOND()
      COUNT(11) = COUNT(11) + 1
C     ************************
C     MAKE LDLT-DECOMPOSITION.
C     ************************
      TMU = MU
      IF(.NOT. LDLT(TMU, TRNEW, BADMU, W, IW)) GOTO 8888
      RNEW = TRNEW
      TIME(11) = TIME(11) + SECOND() - ST1
      CALL WRINFO(NOR, 1, W, IW)
C     ***********************************
C     CHECK IF IT IS AN ACCEPTABLE SHIFT.
C     ***********************************
      IF(.NOT. BADMU) GOTO 40
        NBADMU = NBADMU + 1
        IF(.NOT. VISITD) GOTO 20
          CALL ERROR(NOR, 1)
          GOTO 8888
C
C       ***********************
C       NO. TRY AN OTHER SHIFT.
C       ***********************
20      VISITD = T
        IF(.NOT. ITERNO .EQ. 1) GOTO 30
          IF(MU .EQ. 0.0D0) MU = -1.0D-2                                
          MU = MU - 1.0D-2 * DABS(MU)                                   
          TLDL = TIME(11)
          GOTO 10
C
30      MU = MU + 1.0D-1 * (MU - OLDMU)                                 
C       ***************************
C       CHECK IF WE HAVE REACHED B.
C       ***************************
        IF(.NOT. MU .GE. B) GOTO 10
          REACHB = T
          USEDB  = T
          GOTO 10
C
40    VISITD = F
C
      IF(.NOT. ITERNO .GT. 1) GOTO 9999
C       ****************
C       HOW MANY REMAIN.
C       ****************
        REST = (RNEW - ROLD) - OLCPOS
        IF(.NOT. REST .GT. MXREST) GOTO 9999
          NMXRST = NMXRST + 1
          IF(.NOT. VISITM) GOTO 50
            CALL ERROR(NOR, 2)
            GOTO 8888
C
C         *****************************
C         TOO MANY, TRY AN OTHER SHIFT.
C         *****************************
50        MU = ALTMU
          VISITM = T
C         *****************************
C         CHECK IF WE BACK AWAY FROM B.
C         *****************************
          IF(.NOT. MU .LT. B) GOTO 10
            REACHB = F
            USEDB  = F
            GOTO 10
C
8888  DECOMP = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CDELDUP
      SUBROUTINE DELDUP(NU, POINTR, S, MAXNU, MINNU, NUMDEL, PMAX,
     +                  FINAL)
C
C **********************************************************************
C
C     PURPOSE -
C
C         DELETES DUPLICATE EIGENPAIRS TO THE LEFT OF THE SHIFT. THERE
C         ARE TWO CASES, VIZ.
C
C         T HAS NUMERICALLY MULTIPLE EIGENVALUES. LIKELY TO HAPPEN IF
C         WE ARE CLOSE TO AN EIGENVALUE OF (HIGH) MULTIPLICITY.
C
C         THE EIGENVECTOR HAS A SMALL COMPONENT IN THE DIRECTION OF
C         THE STARTINGVECTOR. LIKELY TO HAPPEN IF MU IS CLOSE TO THE
C         CORRESPONDING EIGENVALUE AND THE PAIR IS ALREADY COMPUTED.
C         WE CHECK ON S(1,I).
C
C
C     INPUT PARAMETERS -
C
C         MAXNU  = MAX(NU(I)), I=1, ..., P.
C         MINNU  = MIN(NU(I)), I=1, ..., P.
C         NUMDEL = NUMBER OF DELETED VECTORS IN CASE 1.
C         FINAL  = TRUE IF THE LAST LANCZOS STEP HAS BEEN TAKEN.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      INTEGER   PMAX                                                    
      DOUBLE PRECISION                                                  
     +          ADNU, FACTOR, SRELPR, MAXNU, MINNU, NU(1), NU1, NU2, S1,
     +          S2, S(PMAX, 1), TOL, TOLBPI, TOLDEL, TOLLDL, TOLPDM,
     +          TOLS1I, TOLZBT, TOLZNU, TSMAX, WRR
      INTEGER   CNEG, CPOS, I, ITERNO, LEN, N, NUMDEL, OLCPOS, P,       
     +          POINT, POINTR(1), REST, RNEW, ROLD, S1IDEL, TCONV, WRI
      LOGICAL   FINAL, USEMX, WRL, ZERBET                               
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,
     +                  TOLZNU, TOLLDL(3)
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
C
C     ************************************************
C     NUMDEL = NUMBER OF DELETED PAIRS, FIRST REASON.
C     S1IDEL =                        , SECOND REASON.
C     ************************************************
      NUMDEL = 0
      S1IDEL = 0
      IF(CNEG .EQ. 0 .OR. ITERNO .EQ. 1) GOTO 9999
C     **********************
C     COMPUTE THE TOLERANCE.
C     **********************
      TOLDEL = (MAXNU - MINNU) * FACTOR
      TOL = TOLDEL * 1.0D-4                                             D MOD.
      IF(MAXNU .EQ. MINNU) TOLDEL = 1.0D0                               D MOD.
      I = 1
C     ***********************************
C     IABS, SINCE POINTR MAY BE NEGATIVE.
C     ***********************************
      POINT = IABS(POINTR(I))
      NU2 = NU(POINT)
C     *****************************************************************
C     S1 AND S2 ARE THE CORRESPONDING TOP ELEMENTS OF THE EIGENVECTORS.
C     *****************************************************************
      S2 = DABS(S(1, POINT))                                            D MOD.
C
C     *******************************************
C     LEN = THE NUMBER OF LAMBDAS IN (OLDMU, MU).
C     *******************************************
      LEN = POINTR(PMAX + 1) - CPOS
10    IF(.NOT. I .LT. LEN) GOTO 40
        NU1 = NU2
        S1 = S2
        I = I + 1
        POINT = IABS(POINTR(I))
        NU2 = NU(POINT)
        S2 = DABS(S(1, POINT))                                          D MOD.
C
        ADNU = DABS(NU1 - NU2)                                          D MOD.
C       ********************************************************
C       THE CONTROL CONSISTS OF TWO STEPS.
C         IS ADNU SMALL ENOUGH (WE MUST HAVE THIS CONTROL SINCE
C         THE SECOND MAY BE TRUE IF ADNU IS BIG).
C
C         CHECK ON ADNU*(THE QUOTIENT BETWEEN THE TOP ELEMENTS).
C
C       TO AVOID DIVISION WITH ZERO, WE DO NOT FORM THE
C       QUOTIENT EXPLICITLY.
C       ********************************************************
        IF(.NOT. ADNU .LE. TOLDEL) GOTO 10
C
          TSMAX = DMAX1(S1, S2)                                         D MOD.
          IF(.NOT. ADNU * DMIN1(S1, S2) .LE. TSMAX * TOL) GOTO 10       D MOD.
            IF(TSMAX .EQ. 0.0D0) GOTO 10                                D MOD.
C
C           ******************************************************
C           CHECK SO THAT WE DO NOT HAVE TWO NOT CONVERGED VALUES.
C           ******************************************************
            IF(POINTR(I - 1) .LT. 0 .AND. POINTR(I) .LT. 0) GOTO 10
C
C           ******************
C           MAKE THE ROTATION.
C           ******************
            CALL ROTATE(S, POINTR(I - 1), FINAL, PMAX, NUMDEL)
C
            I = I + 1
            IF(I .GE. LEN) GOTO 40
            POINT = IABS(POINTR(I))
            NU2 = NU(POINT)
            S2 = DABS(S(1, POINT))                                      D MOD.
            GOTO 10
C
C     ********************************
C     DELETE THOSE WITH SMALL S(1, I).
C     ********************************
40    DO 50 I = 1, LEN
        POINT = POINTR(I)
        IF(POINT .LE. 0) GOTO 50
        IF(.NOT. DABS(S(1, POINT)) .LE. TOLS1I) GOTO 50                 D MOD.
          POINTR(I) = 0
          S1IDEL = S1IDEL + 1
50      CONTINUE
C
C     **************
C     UPDATE POINTR.
C     **************
      IF(FINAL) CALL CMPRSS(POINTR, PMAX)
C
C
9999  WRI(2) = S1IDEL
C     ************
C     ADJUST CNEG.
C     ************
      CNEG = CNEG - NUMDEL - S1IDEL
C
      RETURN
      END
CERROR
      SUBROUTINE ERROR(NOR,LEN)
C
C **********************************************************************
C
C     PURPOSE -
C
C         ERROR HANDLING. THE FIRST TIME ERROR IS CALLED AN ERROR
C         MESSAGE IS WRITTEN. ON THE SECOND CALL A HEADING IS WRITTEN.
C         THE REMAINING CALLS GIVE A TRACEBACK.
C
C
C     INPUT PARAMETERS -
C
C         NOR    = NUMBER OF THE CALLING ROUTINE.
C         LEN    = LOCAL ERROR NUMBER IN THAT ROUTINE (ONLY USED IN
C                  THE FIRST CALL).
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          ALTMU, MU, NEXTMU, OLDMU, RDUMP, TIME, SMALL            R MOD.
      INTEGER   CNEG, COUNT, CPOS, DUMMY, ENOR, ERRNO, IDUMP, ITERNO,   I***
     +          LEN, MSGLVL, N, NBADMU, NMXRST, NOR, OLCPOS, P, REST,
     +          RNEW, ROLD, TCONV, TL, TN, NERR  , NOUT , UNOR, VER
      LOGICAL   FIRST, USEMX, ZERBET                                    L***
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
      COMMON   /STLMPR/ MSGLVL, NERR  , NOUT
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMVR/ SMALL, VER
C
      DATA      ENOR /4/
C
      COUNT(ENOR)=COUNT(ENOR)+1
      IF(COUNT(ENOR) .EQ. 2) GOTO 9999
      TN=NOR
      TL=LEN
C
C     *********************************************
C     CHECK IF LEN LIES IN THE CORRECT INTERVAL.
C     *********************************************
      IF(TL .LE. 0 .OR. TL .GT. 99) TL=0
      FIRST=COUNT(ENOR) .EQ. 1
      IF(.NOT. FIRST) GOTO 5
C       *************************
C       COMPUTE THE ERROR NUMBER.
C       *************************
        ERRNO=100*TN+TL
5     IF(MSGLVL .EQ. 0) GOTO 9999
C
      IF(.NOT. FIRST) GOTO 15
C       ******************
C       WRITE THE HEADING.
C       ******************
        IF(MSGLVL .LE. 3) WRITE(NERR  ,10)
10      FORMAT(1H1///1X,60(1H*),2(/1X,1H*,58X,1H*)/1X,1H*,10X,
     +         38H E R R O R        I N        S T L M  ,
     +         10X,1H*,2(/1X,1H*,58X,1H*)/1X,60(1H*))
C
        IF(MSGLVL .EQ. 4) WRITE(NERR  ,11)
11      FORMAT(1H1///1X,120(1H*)/1X,120(1H*)/1X,120(1H*)/
     +         2(1X,3H***,114X,3H***/),1X,3H***,
     +         38X,38H E R R O R        I N        S T L M  ,38X,3H***
     +         /2(1X,3H***,114X,3H***/),1X,120(1H*)/1X,120(1H*)/
     +         1X,120(1H*))
C
        WRITE(NERR  ,12) ERRNO,NOR,LEN
12      FORMAT(//4X,8HERRNO  =,I18/4X,8HNOR    =,I18/
     +         4X,8HLEN    =,I18//4X,17HERROR OCCURRED IN)
C
C     ****************************
C     WRITE THE TRACEBACK HEADING.
C     ****************************
15    IF(COUNT(ENOR) .EQ. 3) WRITE(NERR  ,20)
20    FORMAT(/4X,9HTRACEBACK/4X,9(1H-))
C
      IF(.NOT. TN .EQ. 0) GOTO 40
        WRITE(NERR  ,30)
30      FORMAT(4X,13H***** UNKNOWN)
        GOTO 9999
C
C     **************************************************
C     WRITE INFORMATION FOR THE ROUTINE WITH NUMBER NOR.
C     AN ERROR MESSAGE OR THE ROUTINE NAME IS WRITTEN.
C     **************************************************
40    IF(NOR .GE. 100) GOTO 5000
      GOTO(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,
     +     1400,1500,1600,1700,1800,1900,2000,2100,2200,2300, 9999,
     +     2500), NOR
C
100   WRITE(NERR  ,110)
110   FORMAT(4X,5HALLOC)
      IF(.NOT. FIRST) GOTO 9999
      WRITE(NERR  ,120) IDUMP
120   FORMAT(4X,28HALL WORKING VECTORS OCCUPIED/
     +       4X,8HID     =,I18/4X,8HACTION =,I18)
      GOTO 9999
C
200   WRITE(NERR  ,210)
210   FORMAT(4X,6HALWAYS)
      IF(.NOT. FIRST) GOTO 9999
      IF(LEN .EQ. 1) WRITE(NERR  ,220) RDUMP
      IF(LEN .EQ. 2) WRITE(NERR  ,230) REST
      IF(LEN .EQ. 3) WRITE(NERR  ,240)
      IF(LEN .EQ. 3) WRITE(NERR  ,250)
      IF(LEN .EQ. 4) WRITE(NERR  ,260)
      IF(LEN .EQ. 4) WRITE(NERR  ,240)
      IF(LEN .EQ. 4) WRITE(NERR  ,270) REST
220   FORMAT(4X,38HTHE M-MATRIX IS PRESUMABLY INDEFINITE./
     +       4X,8HDOT    =,1P1D18.7)                                    F MOD.
230   FORMAT(4X,17HREST IS NEGATIVE./4X,8HREST   =,I18)
240   FORMAT(4X,28HCAN NOT COMPUTE A NEW SHIFT.)
250   FORMAT(4X,22HCNEG PAIRS ARE STORED.)
260   FORMAT(4X,20HREST IS POSITIVE AND)
270   FORMAT(4X,26HCNEG PAIRS ARE NOT STORED./4X,8HREST   =,I18)
      GOTO 9999
C
300   WRITE(NERR  ,310)
310   FORMAT(4X,6HDECOMP)
      IF(.NOT. FIRST) GOTO 9999
      IF(.NOT. LEN .EQ. 1) GOTO 330
        WRITE(NERR  ,350) MU
        IF(.NOT. IDUMP(2) .EQ. 1) GOTO 9999
          WRITE(NERR  ,370) RDUMP,IDUMP(1)
          GOTO 9999
330   WRITE(NERR  ,380) MU,REST
      GOTO 9999
350   FORMAT(4X,39HLDLT FAILED FOR TWO CONSECUTIVE SHIFTS./
     +       4X,8HMU     =,1P1D18.7)                                    F MOD.
370   FORMAT(4X,8HRDUMP  =,1P1D18.7/4X,8HI      =,I18)                  F MOD.
380   FORMAT(4X,52HREST GREATER THAN MXREST FOR TWO CONSECUTIVE SHIFTS./
     +       4X,8HMU     =,1P1D18.7/4X,8HREST   =,I18)                  F MOD.
C
400   GOTO 9999
C
500   WRITE(NERR  ,510)
510   FORMAT(4X,6HFRSTIT)
      IF(.NOT. FIRST) GOTO 9999
      WRITE(NERR  ,520) MU
520   FORMAT(4X,35HFIRST SHIFT GREATER THAN LAMBDA(N).
     +       /4X,8HMU     =,1P1D18.7)                                   F MOD.
      GOTO 9999
C
600   WRITE(NERR  ,610)
610   FORMAT(4X,6HIMTQL2)
      IF(.NOT. FIRST) GOTO 9999
      WRITE(NERR  ,620) IDUMP(1)
620   FORMAT(4X,40HMORE THEN 30 ITERATIONS ARE REQUIRED TO
     +       /4X,30HDETERMINE NU-EIGENVALUE NUMBER,I10)
      GOTO 9999
C
700   WRITE(NERR  ,710)
710   FORMAT(4X,6HINITLA)
      GOTO 9999
C
800   WRITE(NERR  ,810)
810   FORMAT(4X,5HINITU)
      GOTO 9999
C
900   WRITE(NERR  ,910)
910   FORMAT(4X,2HIO)
      IF(.NOT. FIRST) GOTO 9999
      IF(LEN .EQ. 1) WRITE(NERR  , 920)
920   FORMAT(4X,15HDAFILE IS FULL.)
      IF(LEN .EQ. 2) WRITE(NERR  , 930) IDUMP(1)
930   FORMAT(4X,22HTRIED TO READ RECORD =, I6)
      IF(LEN .EQ. 3) WRITE(NERR  , 940) IDUMP(1)
940   FORMAT(4X,23HTRIED TO WRITE RECORD =, I6)
      GOTO 9999
C
1000  WRITE(NERR  ,1010)
1010  FORMAT(4X,6HLANCZO)
      IF(.NOT. FIRST) GOTO 9999
      WRITE(NERR  ,220) RDUMP
      GOTO 9999
C
1100  WRITE(NERR  ,1110)
1110  FORMAT(4X,4HLDLT)
      GOTO 9999
C
1200  WRITE(NERR  ,1210)
1210  FORMAT(4X,6HMULVEC)
      GOTO 9999
C
1300  WRITE(NERR  ,1310)
1310  FORMAT(4X,5HOPINV)
      GOTO 9999
C
1400  WRITE(NERR  ,1410)
1410  FORMAT(4X,3HOPM)
      GOTO 9999
C
1500  WRITE(NERR  ,1510)
1510  FORMAT(4X,6HPREPSV)
      GOTO 9999
C
1600  WRITE(NERR  ,1610)
1610  FORMAT(4X,6HRANDVC)
      GOTO 9999
C
1700  WRITE(NERR  ,1710)
1710  FORMAT(4X,6HSAVEXL)
      IF(.NOT. FIRST) GOTO 9999
      WRITE(NERR  , 1720) RDUMP
1720  FORMAT(4X,44HINSUFFICIENT STORAGE TO STORE AN EIGENVALUE./
     +       4X,8HLAMBDA =, 1P1D18.7)                                   F MOD.
      GOTO 9999
C
1800  WRITE(NERR  ,1810)
1810  FORMAT(4X,6HSCPROD)
      GOTO 9999
C
1900  WRITE(NERR  ,1910)
1910  FORMAT(4X,6HSELDOM)
      IF(.NOT. FIRST) GOTO 9999
      IF(LEN .EQ. 1) WRITE(NERR  ,220) RDUMP
      IF(LEN .EQ. 2) WRITE(NERR  ,1920) REST
      IF(LEN .EQ. 3) WRITE(NERR  ,230) REST
1920  FORMAT(4X,43HEIGENPAIRS ARE MISSING, BUT NO CONVERGENCE./
     +       4X,8HREST   =,I18)
      GOTO 9999
C
2000  WRITE(NERR  ,2010)
2010  FORMAT(4X,4HSTLM//4X,20HTRACEBACK COMPLETED.///
     +       1X,13H*** ERROR ***)
      GOTO 9999
C
2100  WRITE(NERR  ,2110)
2110  FORMAT(4X,6HSUBVEC)
      GOTO 9999
C
2200  WRITE(NERR  ,2210)
2210  FORMAT(4X,6HTRANSF)
      GOTO 9999
C
2300  WRITE(NERR  ,2310)
2310  FORMAT(4X,6HTRIDIG)
      GOTO 9999
C
2500  WRITE(NERR  , 2510)
2510  FORMAT(4X,5HINITD)
      IF(.NOT. FIRST) GOTO 9999
      IF(LEN .EQ. 1) WRITE(NERR  , 2515)
2515  FORMAT(4X, 12HILLEGAL VER.)
      IF(LEN .EQ. 2) WRITE(NERR  , 2520)
2520  FORMAT(4X, 15HN LESS EQUAL 0.)
      IF(LEN .EQ. 3) WRITE(NERR  , 2525)
2525  FORMAT(4X, 15HILLEGAL PROFIL.)
      IF(LEN .EQ. 4) WRITE(NERR  , 2530)
2530  FORMAT(4X, 18HMAXL LESS EQUAL 0.)
      IF(LEN .EQ. 5) WRITE(NERR  , 2535)
2535  FORMAT(4X, 20HMAXREC LESS EQUAL 0.)
      IF(LEN .EQ. 6) WRITE(NERR  , 2540) IDUMP(1)
2540  FORMAT(4X, 11HERROR IN D., I10)
      IF(LEN .EQ. 7) WRITE(NERR  , 2545) IDUMP(1)
2545  FORMAT(4X, 11HERROR IN M., I10)
      IF(LEN .EQ. 8) WRITE(NERR  , 2550)
2550  FORMAT(4X, 15HILLEGAL MSGLVL.)
      IF(LEN .EQ. 9) WRITE(NERR  , 2555)
2555  FORMAT(4X, 13HILLEGAL PMAX.)
      IF(LEN .EQ. 10) WRITE(NERR  , 2560)
2560  FORMAT(4X, 15HTOO SMALL MAXW.)
      IF(LEN .EQ. 11) WRITE(NERR  , 2565)
2565  FORMAT(4X, 16HTOO SMALL MAXIW.)
      IF(LEN .EQ. 12) WRITE(NERR  , 2570)
2570  FORMAT(4X, 15HB LESS EQUAL A.)
      IF(LEN .EQ. 13) WRITE(NERR  , 2575)
2575  FORMAT(4X, 19HILLEGAL FILE UNITS.)
      GOTO 9999
C
C     ***********
C     USER CALLS.
C     ***********
5000  UNOR = NOR - 99
      IF(.NOT. FIRST) GOTO 9999
      GOTO(5100, 5200, 5300, 5400, 5500, 5600, 5700), UNOR
5100  WRITE(NERR  , 5110) VER
5110  FORMAT(4X, 3HLDL, I1)
      GOTO 9999
5200  WRITE(NERR  , 5210) VER
5210  FORMAT(4X, 3HOPM, I1)
      GOTO 9999
5300  WRITE(NERR  , 5310) VER
5310  FORMAT(4X, 3HSOL, I1)
      GOTO 9999
5400  WRITE(NERR  , 5410) VER
5410  FORMAT(4X, 3HMUL, I1)
      GOTO 9999
5500  WRITE(NERR  , 5510) VER
5510  FORMAT(4X, 3HRAN, I1)
      GOTO 9999
5600  WRITE(NERR  , 5610) VER
5610  FORMAT(4X, 3HSCP, I1)
      GOTO 9999
5700  WRITE(NERR  , 5710) VER
5710  FORMAT(4X, 3HSUB, I1)
      GOTO 9999
C
9999  RETURN
      END
CFREEID
      SUBROUTINE FREEID(ID)
C
C **********************************************************************
C
C     PURPOSE -
C
C         RELEASES VEC(ID) FROM W.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      INTEGER   ACTIVE, ID, MV, MXNEW, MXOLD, MXRST, NIL,               I***
     +          SCPX, SOLCPX, V, WAD, X
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMWH/ WAD(2), ACTIVE(2)
C
      IF(ID .EQ. NIL) GOTO 9999
      IF(ID .EQ. ACTIVE(1)) ACTIVE(1) = NIL
      IF(ID .EQ. ACTIVE(2)) ACTIVE(2) = NIL
C
9999  RETURN
      END
CFRSTIT
      LOGICAL FUNCTION FRSTIT(W, IW)                                    L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         DETERMINES WHAT METHOD SHOULD BE USED FOR THE
C         ORTHOGONALIZATION OF THE STARTINGVECTOR. INITIALIZES POPT,
C         COPT, AND CONSTANTS USED IN LATER UPDATES. GIVES VALUES TO
C         /STLMID/.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          CK, CL, COEFF, DOT, SECOND, ST, TCK, TCL, TEMP, TIME,   R MOD.
     +          TIMTQL, TLDL, TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1)
      INTEGER   ADR, ALPHA, BETA, BETAPI, CNEG, CNEGF, CONV, COPT,      I***
     +          COUNT, CPOS, DUMMY, FREE, I, ITERNO, IW(1), L, LAMBDA,
     +          LIM, MV, MXNEW, MXOLD, MXREST, MXRST, N, NBADMU, NIL,
     +          NMXRST, NOACTN, NOR, NU, OLCPOS, P, PFCONV, PMAX,
     +          POINTR, POPT, REST, RFIRST, RNEW, ROLD, S, SAVE, SAVFRE,
     +          SCPX, SCR, SOLCPX, TCONV, TCOPT, TPOPT, V, X
      LOGICAL   F, FINAL, IMTQL2, MEQI, MULVEC, OPINV, OPM, RANDVC,     L***
     +          SAFRST, SCPROD, SUBVEC, T, UPDATE, USEID, USEMX, USEMX1,
     +          REACHB, USSMXR, USEDB
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, USEMX1
      COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
      COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTS/ REACHB, USSMXR, USEDB
      COMMON   /STLMTF/ T, F
C
      DATA      NOR  / 5 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      FRSTIT = T
C
      RFIRST = RNEW
      IF(.NOT. (RFIRST .EQ. N .AND. .NOT. SAFRST)) GOTO 10
C       *********************************************************
C       THE FIRST SHIFT IS GREATER THAN THE LAST EIGENVALUE, AND
C       WE DO NOT WISH TO SAVE THOSE THAT MIGHT CONVERGE TO THE
C       LEFT OF THE SHIFT. HE MUST IN THE TOMS VERSION.
C       *********************************************************
        CALL ERROR(NOR, 1)
        GOTO 8888
C
10    USEID = T
C     *************************************
C     CHECK IF THE USER HAS GIVEN /STLMID/.
C     HE MAY NOT IN THE TOMS VERSION.
C     *************************************
      IF(.NOT. V .EQ. X) GOTO 15
        USEID = F
        MV = 0
        V = PMAX + 1
        X = V + 2
        IF(MEQI) V = MV
C
C     *********************************************
C     TLDL = TIME FOR THE FIRST LDLT-DECOMPOSITION.
C     *********************************************
15    TLDL = TIME(3) - TLDL
C     ****************************************************************
C     IF( DEFAULT MX OR UPDATE OF POPT, COPT ) THEN ... ELSE GOTO 100.
C     ****************************************************************
      IF(.NOT. ((USEMX .AND. .NOT. USEMX1) .OR. UPDATE)) GOTO 100
        IF(.NOT. RANDVC(V + 1, W, IW, F)) GOTO 8888
        IF(.NOT. TOPM .LE. 0.0D0) GOTO 20                               D MOD.
          TOPM = 0.0D0                                                  D MOD.
          IF(MEQI) GOTO 20
C         ************************************************************
C         TOPM = TIME FOR ONE Y = M * X. THE CONTROL ABOVE IS
C         MADE TO SEE IF THE USER HAS SUPPLIED A VALUE ON TOPM OR NOT.
C         HE MAY NOT SUPPLY VALUES IN THE TOMS VERSION.
C         ************************************************************
          TOPM = SECOND()
          IF(.NOT. OPM(MV + 1, V + 1, W, IW)) GOTO 8888
          TOPM = SECOND() - TOPM
C
20      IF(.NOT. TVECOP .LE. 0.0D0) GOTO 40                             D MOD.
          TEMP = 3 * N
C         *********************************************
C         HOW MANY LOOPS. ADAPTED FOR A MACHINE WHERE A
C         MULTIPLICATION TAKES ABOUT 1E-5 SECONDS.
C         *********************************************
          L = DMAX1(1.5D0, 5.0D-1 / (TEMP * 1.0D-5))                    D MOD.
          IF(L .GT. 300) L = 300
          TVECOP = SECOND()
C         ***************************
C         TIME THE VECTOR OPERATIONS.
C         ***************************
          ADR = SAVE
          IF(MEQI) ADR = NOACTN
          DO 30 I = 1, L
            IF(.NOT. MULVEC(MV + 1, V + 1, 0.4D0, W, IW, ADR, 1))       D MOD.
     +                 GOTO 8888
            IF(.NOT. SUBVEC(V + 1, MV + 1, 0.4D0, W, IW, SAVE, 1))      D MOD.
     +                 GOTO 8888
            IF(.NOT. SCPROD(MV + 1, V + 1, DOT, W, IW, SAVE, 1))
     +                 GOTO 8888
30          CONTINUE
          TEMP = 3 * L
C         ****************************************
C         COMPUTE THE MEAN VALUE FOR N OPERATIONS.
C         ****************************************
          TVECOP = (SECOND() - TVECOP) / TEMP
C
C       ****************************************************
C       IF ONLY MX-DECISION, WE DO NOT NEED THE TIMES BELOW.
C       ****************************************************
40      IF(.NOT. UPDATE) GOTO 90
          IF(.NOT. TOPINV .LE. 0.0D0) GOTO 50                           D MOD.
C           ********************************************
C           TOPINV = TIME FOR ONE X =(LDL(T)**(-1)) * B.
C           ********************************************
            TOPINV = SECOND()
            IF(.NOT. OPINV(MV + 1, V + 1, W, IW)) GOTO 8888
            TOPINV = SECOND() - TOPINV
C
50        IF(.NOT. TIMTQL .LE. 0.0D0) GOTO 70                           D MOD.
C           *********************************************
C           SIMULATE THE TIME FOR ONE COMPUTATION OF THE
C           EIGENPAIRS OF THE TRIDIAGONAL MATRIX T.
C           *********************************************
            L = MIN0(PMAX, 20)
            FINAL = F
            TIMTQL = SECOND()
            DO 60 I = 1, L
              ADR = ALPHA + I - 1
              W(ADR) = I
              ADR = BETA + I - 1
              W(ADR) = 1.0D0                                            D MOD.
              IF(I .EQ. L) FINAL = T
              IF(.NOT. IMTQL2(W(ALPHA), W(BETA), W(NU), W(SCR), W(S),
     +                        I, PMAX, FINAL)) GOTO 8888
60            CONTINUE
            TIMTQL = SECOND() - TIMTQL
C
          TEMP = L**3
C         ********************************************
C         COMPUTE C WHERE TIMTQL = C * (DIM. OF T)**3.
C         COMPUTE COEFFICIENTS USED IN NEWPC AND HALF.
C         ********************************************
          TIMTQL = TIMTQL / TEMP
70        COEFF(1) = TLDL + 2.0D0 * (TOPM + TVECOP)                     D MOD.
          IF(MEQI) COEFF(1) = COEFF(1) - TVECOP
          COEFF(2) = TOPINV + 7.0D0 * TVECOP + TOPM                     D MOD.
          IF(MEQI) COEFF(2) = COEFF(2) - TVECOP
          COEFF(3) = 2.0D0 * TVECOP                                     D MOD.
          COEFF(4) = TVECOP
C
C         ******************************************
C         COMPUTE POPT AND COPT, AND CHECK IF OK.
C         THE TEMPORARY VARIABLES ARE TO PASS PFORT.
C         ******************************************
          TCK = CK
          TCL = CL
          TCOPT = COPT
          TPOPT = POPT
          CALL HALF(TCK, TCL, TPOPT, TCOPT, LIM)
          COPT = TCOPT
          POPT = TPOPT
          IF(.NOT. LIM .LE. 3) GOTO 80
            POPT = MIN0(N, PMAX)
            TEMP = PMAX - 10
            COPT = TEMP * 0.5D0 + 0.5D0                                 D MOD.
            IF(COPT .LT. 1) COPT = 1
C
80        IF(.NOT. USSMXR) MXREST = COPT
C
90      CALL FREEID(V + 1)
        CALL FREEID(MV + 1)
C
100   TEMP = POPT + 1
C     ******************************
C     SHOULD WE SAVE THE MV-VECTORS.
C     ******************************
      IF(USEMX .AND. .NOT. USEMX1) USEMX = TEMP * TVECOP .LT. TOPM
      IF(USEMX .AND. MEQI) USEMX = F
      IF(USEID) GOTO 110
C       *************************************************
C       INITIALIZE /STLMID/ IF THE USER HAS NOT DONE IT.
C       HE MAY NOT DO IT IN THE TOMS VERSION.
C       *************************************************
        X = V + PMAX + 1
        IF(.NOT. USEMX) GOTO 110
          MXNEW = V + PMAX + 1
          MXOLD = MXNEW + PMAX
          MXRST = MXOLD + PMAX
          X = MXRST + PMAX
C
110   IF(.NOT. UPDATE) GOTO 9999
C       ******************
C       CORRECT FOR USEMX.
C       ******************
        I = 0
        IF(USEMX) I = 1
        TEMP = I
        COEFF(3) = COEFF(3) + (1.0D0 - TEMP) * TOPM * 0.5D0 +           D MOD.
     +             TEMP * TVECOP * 0.5D0                                D MOD.
        COEFF(4) = COEFF(4) + TEMP * TVECOP * 0.5D0                     D MOD.
        GOTO 9999
C
8888  FRSTIT = F
      CALL ERROR(NOR, NOR)
C
9999  IF(FRSTIT) CALL WRINFO(NOR, 1, W, IW)
      TIME(NOR) = TIME(NOR) + SECOND() - ST
      TSAVE = TSAVE + TIME(NOR)
      RETURN
      END
CHALF
      SUBROUTINE HALF(K, L, TP, TC, LIM)
C
C **********************************************************************
C
C     PURPOSE -
C
C         COMPUTES COEFFICIENTS IN THE POLYNOMIAL WHOSE ZERO GIVES
C         COPT. COMPUTES THE ZERO AND POPT.
C
C
C     INPUT PARAMETERS -
C
C         K,L =  WE ASSUME THAT P=K*C+L, WHERE C = NUMBER OF CONVERGED
C                EIGENPAIRS FOR P LANCZOS STEPS.
C
C
C     OUTPUT PARAMETERS -
C
C         TP  = TEMPORARY POPT.
C         TC  = TEMPORARY COPT.
C         LIM = 6 NORMAL TERMINATION.
C               1 K OR L NEGATIVE.
C               2 CAN NOT FIND COPT.
C               3 TP OR TL IS ABSURD.
C               4 PMAX IS TOO SMALL COMPARED WITH TP.
C               5 TP IS GREATER THAN N.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          A, B, C, CK, CL, COEFF, FM, FN, FP, K, L, MEAN, NEG,    R MOD.
     +          POS, TEMP, TIMTQL, TLDL, TOPINV, TOPM, TPRED, TSAVE,
     +          TVECOP, WRR
      INTEGER   CNEG, CONV, COPT, CPOS, ITERNO, LIM, LOOP, MXREST, N,   I***
     +          OLCPOS, P, PFCONV, PMAX, POPT, REST, RNEW, ROLD, TC,
     +          TCONV, TP, WRI
      LOGICAL   MEQI, UPDATE, USEMX, WRL, ZERBET                        L***
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
C
      LIM = 6
      IF(.NOT. (K .LE. 0.0D0 .OR. L .LE. 0.0D0)) GOTO 10                D MOD.
C       **************************
C       SHOULD NOT HAPPEN, BUT ...
C       **************************
        LIM = 1
        GOTO 9999
C
C     **************************************************
C     COMPUTE COEFFICIENTS IN TIME OPTIMIZATION FORMULA.
C     **************************************************
10    A = 2.0D0 * TIMTQL * K**3                                         D MOD.
      B = K * COEFF(4) + TVECOP * K * K + 3.0D0 * TIMTQL * K * K * L    D MOD.
      C = COEFF(1) + L * (COEFF(2) + L * (TVECOP + TIMTQL * L))
      WRR(1) = A
      WRR(2) = B
      WRR(3) = C
      POS = DSQRT(C / B)                                                D MOD.
      FP = A * POS**3
      NEG = 0.0D0                                                       D MOD.
      FN =  - C
      LOOP = 0
C
C     ********************************************************
C     SEEK AN APPROXIMATION TO A ZERO TO
C     A * X**3 + B * X**2 - C = 0, BY SUCCESSIVE BISECTIONS.
C     LOOP IS USED TO AVOID AN INFINITE LOOP IN CASE SOMETHING
C     SHOULD GO WRONG. (LOOP = 5 IS NORMAL).
C     ********************************************************
20    IF(.NOT. (POS - NEG .GT. 4.5D-1 / K .AND. LOOP .LE. 100)) GOTO 50 D MOD.
        MEAN = (NEG + POS) * 0.5D0                                      D MOD.
        FM = MEAN * MEAN * (A * MEAN + B) - C
        IF(FM .EQ. 0.0D0) GOTO 70                                       D MOD.
        IF(.NOT. FM .GT. 0.0D0) GOTO 30                                 D MOD.
          POS = MEAN
        GOTO 40
30        NEG = MEAN
40      LOOP = LOOP + 1
        GOTO 20
C
50    IF(.NOT. LOOP .GT. 100) GOTO 60
        LIM = 2
        GOTO 9999
C
60    MEAN = (NEG + POS) * 0.5D0                                        D MOD.
C     ********************************
C     COMPUTE TEMPORARY POPT AND COPT.
C     ********************************
70    TC = MEAN + 0.5D0                                                 D MOD.
      TEMP = TC
      TP = K * TEMP + L + 0.5D0                                         D MOD.
C     ********************************************************
C     EXPERIENCE HAS SHOWN THAT A TP .LT. 20 CAN CAUSE TROUBLE
C     WITH THE SHIFTING.
C     ********************************************************
      IF(.NOT. TP .LT. 20) GOTO 75
        TP = 20
        TEMP = TP
        TC = DMAX1(1.5D0, (TEMP - L) / K +0.5D0)                        D MOD.
C
75    WRI(1) = TC
      WRI(2) = TP
      WRI(3) = LOOP
C
      IF(.NOT. (TC .LE. 0 .OR. TP .LE. 0)) GOTO 80
C       ****************************
C       SHOULD NEVER HAPPEN, BUT ...
C       ****************************
        LIM = 3
        GOTO 9999
C
80    IF(.NOT. TP .GT. PMAX) GOTO 90
C       ******************************************
C       WE CAN NOT TAKE SO MANY STEPS DUE TO PMAX.
C       ******************************************
        LIM = 4
        TP = PMAX
        TEMP = TP
        TC = DMAX1(1.5D0, (TEMP - L) / K + 0.5D0)                       D MOD.
C
90    IF(.NOT. TP .GT. N) GOTO 9999
C       ***************************************
C       THE PROBLEM IS VERY SMALL, LIMITS POPT.
C       ***************************************
        LIM = 5
        TP = N
        TEMP = TP
        TC = DMAX1(1.5D0, (TEMP - L) / K + 0.5D0)                       D MOD.
C
9999  WRI(4) = LIM
      RETURN
      END
CIMTQL2
      LOGICAL FUNCTION IMTQL2(ALPHA, BETA, D, E, Z, N, NM, FINAL)       L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         A SOMEWHAT MODIFIED VERSION OF THE EISPACK SUBROUTINE
C         IMTQL2. SEE B.S. GARBOW AND J.J. DONGARRA, PATH CHART AND
C         DOCUMENTATION FOR THE EISPACK PACKAGE OF MATRIX EIGENSYSTEM
C         ROUTINES, ARGONNE NATIONAL LABORATORY, 1975.
C
C         THE MAJOR MODIFICATIONS ARE,
C
C                 DIFFERENT NUMBER OF DUMMY PARAMETERS.
C                 INITIALIZATION OF D, E, AND, Z INSIDE THE PROCEDURE.
C                 POSSIBILTY TO COMPUTE EIGENVALUES AND ONLY TOP AND
C                 BOTTOM ELEMENTS OF Z (THE EIGENVECTORS).
C                 NO SORTING IS MADE.
C                 SINGLE PRECISION.
C
C         THE LOOP WHERE THE COMPUTATIONS ARE MADE (DO 240 ...)
C         IS HOWEVER NEARLY INTACT. IN THIS CASE WE WILL DESCRIBE
C         THE PARAMETERS SINCE THE NAMES ARE DIFFERENT FROM WHAT WE
C         ARE USED TO.
C
C         IMTQL2 IS USED AS FOLLOWS IN THE REFERENCING ROUTINE,
C
C
C         IF(.NOT. IMTQL2(W(ALPHA), W(BETA), W(NU), W(SCR), W(S),
C        +                P, PMAX, FINAL)) GOTO 8888
C
C
C     INPUT PARAMETERS -
C
C         ALPHA = DIAGONAL IN THE TRIDIAGONAL MATRIX T.
C         BETA  = SUBDIAGONAL IN T.
C         E     = IS A WORKING AREA.
C         N     = THE DIMENSION OF T.
C         NM    = ROW DIMENSION OF Z, THE EIGENVECTOR MATRIX.
C         FINAL = IF TRUE, COMPUTE THE WHOLE EIGENVECTORS,
C                 IF FALSE, COMPUTE ONLY THE TOP AND BOTTOM ELEMENTS.
C
C
C     OUTPUT PARAMETERS -
C
C         D = CONTAINS THE EIGENVALUES.
C         Z = CONTAINS THE EIGENVECTORS OR THE END ELEMENTS.
C
C
C     NOTE. IN THE ORIGINAL VERSION, IMTQL2 OVERWRITES ALPHA WITH
C     THE EIGENVALUES AND DESTROYS BETA. Z SHOULD FURTHERMORE BE SET TO
C     THE IDENTITY MATRIX. WE THEREFORE COPY ALPHA TO D, BETA TO E,
C     AND Z IS SET TO I (OR THE FIRST AND LAST ROWS OF I) BEFORE
C     ANY COMPUTATION TAKES PLACE.
C
C **********************************************************************
C
C
      INTEGER   NM                                                      I***
      DOUBLE PRECISION                                                  R INS.
     +          ALPHA(1), B, BETA(1), C, D(1), E(1), F, FACTOR, G,      R MOD.
     +          SRELPR, P, R, RDUMP, S, SECOND, ST, TIME, TOLBPI,
     +          TOLLDL, TOLPDM, TOLS1I, TOLZBT, TOLZNU, Z(NM, 1)
      INTEGER   COUNT, DUMMY, ERRNO, I, IDUMP, II, J, K, L, M, MML, N,  I***
     +          NBADMU, NMXRST, NOR, STEP
      LOGICAL   FA, FINAL, T                                            L***
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, FA
      COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,
     +                  TOLZNU, TOLLDL(3)
C
      DATA      NOR  / 6 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      IMTQL2 = T
C
C     *****
C     COPY.
C     *****
      DO 10 I = 1, N
        D(I) = ALPHA(I)
        E(I) = BETA(I)
10      CONTINUE
      E(N) = 0.0D0                                                      D MOD.
C
C     *************
C     Z = IDENTITY.
C     *************
      IF(.NOT. FINAL) GOTO 40
        STEP = 1
        DO 30 I = 1, N
          DO 20 J = 1, N
            Z(I, J) = 0.0D0                                             D MOD.
20          CONTINUE
          Z(I, I) = 1.0D0                                               D MOD.
30        CONTINUE
        GOTO 60
C
C     **********************************
C     Z = IDENTITY, FIRST AND LAST ROWS.
C     **********************************
40    STEP = N - 1
      DO 50 I = 1, N
        Z(1, I) = 0.0D0                                                 D MOD.
        Z(N, I) = 0.0D0                                                 D MOD.
50      CONTINUE
      Z(1, 1) = 1.0D0                                                   D MOD.
      Z(N, N) = 1.0D0                                                   D MOD.
C
60    IF(N .EQ. 1) GOTO 9999
C
      DO 240 L = 1, N
         J = 0
C        ***************************************
C        LOOK FOR SMALL SUB - DIAGONAL ELEMENTS.
C        ***************************************
105      DO 110 M = L, N
            IF(M .EQ. N) GOTO 120
            IF(DABS(E(M)) .LE. SRELPR * (DABS(D(M)) + DABS(D(M + 1))))  D MOD.
     +                                          GOTO 120
110      CONTINUE
C
120      P = D(L)
         IF(M .EQ. L) GOTO 240
         IF(.NOT. J .EQ. 30) GOTO 121
C        ***************
C        ERROR HANDLING.
C        ***************
          IDUMP(1) = L
          CALL ERROR(NOR, 1)
          GOTO 8888
C
121      J = J + 1
C        ***********
C        FORM SHIFT.
C        ***********
         G = (D(L + 1) - P) / (2.0D0 * E(L))                            D MOD.
         R = DSQRT(G * G + 1.0D0)                                       D MOD.
         G = D(M) - P + E(L) / (G + DSIGN(R, G))                        D MOD.
         S = 1.0D0                                                      D MOD.
         C = 1.0D0                                                      D MOD.
         P = 0.0D0                                                      D MOD.
         MML = M - L
C        **********************************
C        FOR I = M - 1 STEP  - 1 DO  -  -
C        **********************************
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF(DABS(F) .LT. DABS(G)) GOTO 150                           D MOD.
            C = G / F
            R = DSQRT(C * C + 1.0D0)                                    D MOD.
            E(I + 1) = F * R
            S = 1.0D0 / R                                               D MOD.
            C = C * S
            GOTO 160
150         S = F / G
            R = DSQRT(S * S + 1.0D0)                                    D MOD.
            E(I + 1) = G * R
            C = 1.0D0 / R                                               D MOD.
            S = S * C
160         G = D(I + 1) - P
            R = (D(I) - G) * S + 2.0D0 * C * B                          D MOD.
            P = S * R
            D(I + 1) = G + P
            G = C * R - B
C           *****************************************************
C           FORM VECTOR.
C           STEP = 1, COMPUTE THE WHOLE VECTORS.
C           STEP = N-1, COMPUTE ONLY THE TOP AND BOTTOM ELEMENTS.
C           *****************************************************
            DO 180 K = 1, N, STEP
               F = Z(K, I + 1)
               Z(K, I + 1) = S * Z(K, I) + C * F
               Z(K, I) = C * Z(K, I) - S * F
180         CONTINUE
C
200      CONTINUE
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0D0                                                   D MOD.
         GOTO 105
240   CONTINUE
C
      GOTO 9999
8888  IMTQL2 = FA
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CINITAD
      SUBROUTINE INITAD(NUM, BASEID, ADR, N, ADDRSS)
C
C **********************************************************************
C
C     PURPOSE - (VER = 1 OR 2)
C
C         INITIALIZES THE ADDRESS VECTOR.
C
C     INPUT PARAMETERS -
C
C         NUM    = THE NUMBER OF VECTORS IN ONE IDENTIFIER GROUP.
C         BASEID = THE IDENTIFIER OF THE FIRST VECTOR IN THE GROUP - 1.
C         ADR    = THE ADDRESS IN THE W-ARRAY.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      INTEGER   ADDRSS(1), ADR, BASEID, I, ID, N, NUM                   I***
C
      IF(NUM .EQ. 0) GOTO 9999
      DO 10 I = 1, NUM
        ID = BASEID + I
        ADDRSS(ID) = ADR
        ADR = ADR - N
10      CONTINUE
C
9999  RETURN
      END
CINITD
      SUBROUTINE INITD(N, A, B, MAXL, PROFIL, PMAX, MXREST, MSGLVL,
     +                 MAXW, MAXIW, DAFILE, MAXREC, KFILE, X, BG,
     +                 TCONV, NLEFT, ERRNO, W, IW, LEN)
C
C **********************************************************************
C
C     PURPOSE -
C
C         INITALIZES ALL THE COMMON BLOCKS. SOME INITIALIZATIONS ARE
C         MADE ONLY SO THAT A PREMATURE TERMINATION WILL NOT CAUSE A
C         SUBSEQUENT ABORT WHEN TRYING TO WRITE THE VALUE OF AN
C         A NOT INITIALIZED VARIABLE. OTHERS ARE MADE BECAUSE WRINF4
C         WOULD OTHERWISE SOMETIMES WRITE VALUES NOT INITIALIZED.
C
C
C     OUTPUT PARAMETERS -
C
C         LEN HAS THE VALUE COMPUTED IN CHECK.
C
C
C         THE REASON FOR THE VARIABLE NAMES BELOW (E.G. N1) IS THAT
C         IT IS FORBIDDEN TO HAVE DUMMY PARAMETERS IN COMMON.
C
C     PLEASE SEE THE USER GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          A, A1, ALTMU, B, B1, CK, CL, COEFF, FACTOR, SRELPR, MU, R MOD.
     +          NEXTMU, OLDMU, RDUMP, SMALL, TEMP, TIME, TIMTQL,
     +          TLDL, TOLBPI, TOLLDL, TOLPDM, TOLS1I, TOLZBT, TOLZNU,
     +          TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1), WRR, BG
      INTEGER   ACTIVE, ADDRSS, ALPHA, BETA, BETAPI, CNEG, CNEGF, CONV, I***
     +          COPT, COUNT, CPOS, D, DAFIL1, DAFILE, DUMMY, ERRNO,
     +          FREE, I, I1, IDUMP, ITERNO, IW(1), K, KFILE, KFILE1,
     +          LAMBDA, LEFTP, LEN, LENADR, LP, M, MAXIW, MAXIW1, MAXL,
     +          MAXL1, MAXRE1, MAXREC, MAXW, MAXW1, MV, MXNEW,
     +          MXOLD, MXREST, MXRST, N, N1, N2, NBADMU, NIL, NMXRST,
     +          NOACTN, NREAD, NU, NUMEIG, NUMVEC, NWRITE,
     +          OLCPOS, P, PFCONV, PMAX, POINTR, POPT, PV1, PV2,
     +          READID, READK, REST, RFIRST, RIGHTC, RIGHTM, RIGHTP,
     +          RNEW, ROLD, S, SAEVAL, SAVE, SAVFRE, SCPX, SCR, SOLCPX,
     +          STADEW, TCONV, V, VER, VTEMP, WAD, WRI,
     +          WRITID, X, TCONV1, NLEFT, ERRNO1, MSGLV1, MSGLVL, PMAX1,
     +          PROFIL, PROFI1, MXRES1, NERR  , NOUT , X1
      LOGICAL   DIAGM, F, MEQI, SAFRST, T, UPDATE,                      L***
     +          USEMX, WRL, ZERBET, REACHB, USSMXR, USEDB
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMAD/ ADDRSS, DAFIL1, KFILE1, LP, MAXL1, NREAD, NWRITE
      COMMON   /STLMCT/ N1, ITERNO, TCONV1, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMDS/ K, M, D, DIAGM
      COMMON   /STLMER/ RDUMP, ERRNO1, IDUMP(2)
      COMMON   /STLMEW/ LEFTP, LENADR, MAXRE1, NUMVEC, RIGHTC, RIGHTM,
     +                  RIGHTP, STADEW
      COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX,
     +                  X1
      COMMON   /STLMIN/ A1, B1, NUMEIG, MAXW1, MAXIW1
      COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N2
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMPL/ PMAX1, POPT, COPT, MXRES1
      COMMON   /STLMPF/ PROFI1
      COMMON   /STLMPR/ MSGLV1, NERR  , NOUT
      COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,
     +                  TOLZNU, TOLLDL(3)
      COMMON   /STLMTS/ REACHB, USSMXR, USEDB
      COMMON   /STLMVR/ SMALL, VER
      COMMON   /STLMWH/ WAD(2), ACTIVE(2)
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
C
C     *****************************
C     INITIALIZE OUTPUT PARAMETERS.
C     *****************************
      X = 0
      BG = A
      WRR(1) = A
      TCONV = 0
      NLEFT = 0
      ERRNO = 0
      WRI(1) = 0
C     ************
C     CHECK INPUT.
C     ************
      COUNT(4) = 0
      CALL CHECK(N, A, B, MAXL, PROFIL, PMAX, MXREST, MSGLVL,
     +           MAXW, MAXIW, DAFILE, MAXREC, KFILE,
     +           ERRNO, W, IW, LEN)
      IF(LEN .GT. 0) GOTO 9999
C
C
C     *************************
C     INITIALIZE COMMON BLOCKS.
C     *************************
C
C
C     **********
C     / STLMAC /
C     **********
      NOACTN = 0
      FREE = 1
      SAVE = 2
      SAVFRE = 3
C     **********
C     / STLMCT /
C     **********
      N1 = N
      ITERNO = 0
      TCONV1 = 0
      CNEG = 0
      CPOS = 0
      OLCPOS = 0
      RNEW = 0
      ROLD = 0
      REST = 0
      P = 0
      USEMX = .TRUE.
      ZERBET = .FALSE.
C     **********
C     / STLMDS /
C     **********
      K = 1
      M = 1
      D = 1
      DIAGM = PROFIL .GE. 2
C     **********
C     / STLMMI /
C     **********
      MEQI = PROFIL .EQ. 3
C     ****************************************
C     PV1 POINTS TO THE START OF ALPHA IN W.
C     PV2 POINTS TO THE START OF POINTR IN IW.
C     ****************************************
      PV1 = 1
      PV2 = 1
      IF(.NOT. VER .EQ. 1) GOTO 10
        M = IW(N) + 1
        IF(MEQI) PV1 = M
        IF(DIAGM .AND. .NOT. MEQI) PV1 = M + N
        IF(.NOT. DIAGM) PV1 = M + IW(N)
        PV2 = N + 1
C       *******
C       SAVE K.
C       *******
        REWIND KFILE
        I1 = IW(N)
        WRITE(KFILE) (W(I), I=1, I1)
C     **********
C     / STLMER /
C     **********
10    RDUMP = 0.0D0                                                     D MOD.
      ERRNO1 = 0
      IDUMP(1) = 0
      IDUMP(2) = 0
C     **********
C     / STLMFT /
C     **********
      CNEGF = 0
      RFIRST = 0
      SAFRST = .FALSE.
C     **********
C     / STLMID /
C     **********
      NIL = 0
      MV = 0
      V = 0
      MXNEW = 0
      MXOLD = 0
      MXRST = 0
      SCPX = 0
      SOLCPX = 0
      X1 = 0
C     **********
C     / STLMIO /
C     **********
      SAEVAL = 1
      READID = 2
      WRITID = 3
      READK = 4
      N2 = N
C     **********
C     / STLMMU /
C     **********
      MU = A
      OLDMU = 0.0D0                                                     D MOD.
      NEXTMU = 0.0D0                                                    D MOD.
      ALTMU = 0.0D0                                                     D MOD.
C     **********
C     / STLMOP /
C     **********
      TOPINV =  -1.0D0                                                  D MOD.
      TOPM =  -1.0D0                                                    D MOD.
      TIMTQL =  -1.0D0                                                  D MOD.
      TVECOP =  -1.0D0                                                  D MOD.
      TPRED = 0.0D0                                                     D MOD.
      TSAVE = 0.0D0                                                     D MOD.
      DO 20 I = 1, 4
        COEFF(I) = 0.0D0                                                D MOD.
20      CONTINUE
      CK = 2.0D0                                                        D MOD.
      CL = 1.0D1                                                        D MOD.
      CONV = 0
      PFCONV = 0
      UPDATE = .TRUE.
      TLDL = 0.0D0                                                      D MOD.
C     ***********
C     / STLMPL /
C     **********
      POPT = 0
      COPT = 0
      MXRES1 = MXREST
      PMAX1 = PMAX
C     **************************************
C     / STLMPR / SEE BLOCKDATA PROGRAM ALSO.
C     **************************************
      MSGLV1 = MSGLVL
C     **********
C     / STLMPF /
C     **********
      PROFI1 = PROFIL
C     **********
C     / STLMPV /
C     **********
      ALPHA = PV1
      BETA = ALPHA + PMAX
      BETAPI = BETA + PMAX
      LAMBDA = BETAPI + PMAX
      NU = LAMBDA + PMAX
      POINTR = PV2
      S = NU + PMAX
      SCR = BETAPI
C     **********
C     / STLMST /
C     **********
      DO 30 I = 1, 24
        TIME(I) = 0.0D0                                                 D MOD.
        COUNT(I) = 0
30      CONTINUE
      NBADMU = 0
      NMXRST = 0
C     **********
C     / STLMTF /
C     **********
      T = .TRUE.
      F = .FALSE.
C     **************************************
C     / STLMTL / SEE BLOCKDATA PROGRAM ALSO.
C     **************************************
C
      TOLBPI = 1.0D-7                                                   D MOD.
      TEMP = N
      TOLS1I = 1.0D-2 / DSQRT(TEMP)                                     D MOD.
      FACTOR = 1.0D-6                                                   D MOD.
      TOLPDM = SMALL * TEMP
      TOLZBT = 1.0D-5                                                   D MOD.
      TOLZNU = SMALL * 1.0D5                                            D MOD.
      TOLLDL(1) = SMALL * 1.0D5                                         D MOD.
      TOLLDL(2) = 1.0D5                                                 D MOD.
      TOLLDL(3) = -1.0D0                                                D MOD.
C     **********
C     / STLMTS /
C     **********
      REACHB = .FALSE.
      USSMXR = MXREST .GE. 1
      USEDB  = .FALSE.
C     **********
C     / STLMUI /
C     **********
C     **********
C     / STLMWH /
C     **********
      ACTIVE(1) = NIL
      ACTIVE(2) = NIL
      WAD(1) = S + PMAX**2
      WAD(2) = WAD(1) + N
C     **********
C     / STLMWR /
C     **********
      DO 50 I = 1, 5
        WRR(I) = 0.0D0                                                  D MOD.
        WRI(I) = 0
        WRL(I) = .TRUE.
50      CONTINUE
C
C
C     **********
C     / STLMAD /
C     **********
      ADDRSS = POINTR + PMAX + 1
      IF(VER .EQ. 3) ADDRSS = 0
      DAFIL1 = DAFILE
      IF(VER. EQ. 3) DAFIL1 = 0
      KFILE1 = KFILE
      IF(VER .GE. 2) KFILE1 = 0
      LP = WAD(2) + N
      IF(VER .EQ. 3) LP = WAD(1)
      MAXL1 = MAXL
      NREAD = 0
      NWRITE = 0
C
C     **********
C     / STLMEW /
C     ***********
      STADEW = LP + MAXL
      IF(.NOT. VER .EQ. 3) GOTO 60
        STADEW = 0
        NUMVEC = 0
        LENADR = 0
        MAXRE1 = 0
        GOTO 80
C
C
C
60    NUMVEC = (MAXW - STADEW + 1) / N
      IF(NUMVEC .LT. 0) NUMVEC = 0
      LENADR = 5 * PMAX + 2 + MAXL
      IF(MEQI) LENADR = PMAX + 1 + MAXL
C
      I1 = ADDRSS + LENADR - 1
      DO 70 I = ADDRSS, I1
        IW(I) = 0
70      CONTINUE
C
C     *****************************************************
C     SINCE USEMX ETC. ARE NOT KNOWN AT THIS STAGE, WE ONLY
C     GIVE ADDRESSES TO V + 1 (AND MV + 1), SINCE THEY MAY
C     BE USED IN FRSTIT.
C     *****************************************************
      VTEMP = 0
      IF(.NOT. MEQI) VTEMP = PMAX + 1
      I1 = ADDRSS + VTEMP
      IF(NUMVEC .GE. 1) IW(I1) = -STADEW
      I1 = ADDRSS + MV
      IF(NUMVEC .GE. 2 .AND. .NOT. MEQI) IW(I1) = -STADEW - N
C
      LEFTP = 1
      RIGHTP = MAXREC
      RIGHTC = MAXREC + 1
      RIGHTM = RIGHTC
      MAXRE1 = MAXREC
C
C
C     **********
C     / STLMIN /
C     **********
80    A1 = A
      B1 = B
      NUMEIG = 0
      MAXW1 = MAXW
      MAXIW1 = MAXIW
C
C     *********************************
C     / STLMVR / SEE BLOCKDATA PROGRAM.
C     *********************************
C
C
9999  RETURN
      END
CINITEW
      SUBROUTINE INITEW(ADDRSS)
C
C **********************************************************************
C
C     PURPOSE - (VER = 1 OR 2)
C
C         THIS ROUTINE UPDATES THE ADDRESS VECTOR TO THE
C         DIRECT ACCESS FILE.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          FC, FN                                                  R MOD.
      INTEGER   ADDRSS(1), ADR, CNEG, COPT, CPOS, I, I1, ITERNO, LEFTP, I***
     +          LENADR, MAXREC, MV, MXNEW, MXOLD, MXREST, MXRST, N, NIL,
     +          NUMMV, NUMMXR, NUMV, NUMVEC, OLCPOS, P, PMAX, POPT,
     +          REMAIN, REST, RIGHTC, RIGHTM, RIGHTP, RNEW, ROLD, SCPX,
     +          SOLCPX, STADEW, TCONV, V, X
      LOGICAL   MEQI, USEMX, ZERBET                                     L***
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMEW/ LEFTP, LENADR, MAXREC, NUMVEC, RIGHTC, RIGHTM,
     +                  RIGHTP, STADEW
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
C
C
      IF(.NOT. USEMX) GOTO 25
C       *****************************
C       SET IF OLCPOS SHOULD BE ZERO.
C       *****************************
        RIGHTC = MAXREC + 1
        RIGHTM = RIGHTC
        IF(.NOT. OLCPOS .GT. 0) GOTO 25
C
C         **************************************
C         ERASE ALL BUT THE X AND MXOLD VECTORS.
C         **************************************
          DO 21 I = 1, MXOLD
            ADDRSS(I) = 0
21          CONTINUE
C
          I1 = MXOLD + OLCPOS + 1
          DO 22 I = I1, X
            ADDRSS(I) = 0
22          CONTINUE
C
C         ***********************************************
C         RIGHTC POINTS TO THE  RECORD HOLDING  THE FIRST
C         VECTOR IN MXOLD. RIGHTM POINTS TO THE LEFT OF
C         THE LAST VECTOR IN MXOLD.
C         ***********************************************
          RIGHTM = ADDRSS(I1 - 1) - 1
          RIGHTC = ADDRSS(MXOLD + 1)
C
          GOTO 35
C
C
C     ************************************
C     WE COME HERE IF USEMX = F. ERASE ALL
C     BUT THE X-VECTORS.
C     ************************************
25    DO 30 I = 1, X
        ADDRSS(I) = 0
30      CONTINUE
C     *************
C     RESET RIGHTP.
C     *************
35    RIGHTP = MAXREC
C
C     *************************************
C     CHECK IF EVERYTHING LIES ON THE FILE.
C     *************************************
      IF(NUMVEC .EQ. 0) GOTO 9999
C
C     *********************************************
C     MAKE AN OPTIMAL DIVISION OF THE STORAGE IN W.
C     *********************************************
      FN = NUMVEC
      FC = COPT
      IF(MEQI) NUMMV = 0
      IF(.NOT. MEQI .AND. .NOT. USEMX) NUMMV = (FN - FC) * 0.5D0        D MOD.
      IF(.NOT. MEQI .AND.       USEMX) NUMMV = (FN - FC * 0.5D0) * 0.5D0D MOD.
C
      IF(NUMMV .LT. 0) NUMMV = 0
C     ************************************
C     COMPUTE OPTIMAL NUMBER OF V-VECTORS.
C     CHECK IF GREATER THAN PMAX+1.
C     ************************************
      NUMV = NUMVEC - NUMMV
      I = PMAX + 1
      IF(NUMV .GT. I) NUMV = I
      IF(NUMMV .GT. I) NUMMV = I
C
      REMAIN = NUMVEC - (NUMV + NUMMV)
      IF(.NOT. USEMX) GOTO 50
C       ***********************************************
C       IF THERE REMAIN VECTORS WE PLACE THEM IN MXRST,
C       IF USEMX = T.
C       ***********************************************
        IF(REMAIN .GT. PMAX) REMAIN = PMAX
        NUMMXR = REMAIN
C
C
50    ADR = - STADEW
C     ***************************
C     HERE ADDRSS IS INITIALIZED.
C     ***************************
      CALL INITAD(NUMMV, MV, ADR, N, ADDRSS)
      CALL INITAD(NUMV, V, ADR, N, ADDRSS)
      IF(USEMX) CALL INITAD(NUMMXR, MXRST, ADR, N, ADDRSS)
C
C
9999  RETURN
      END
CINITLA
      LOGICAL FUNCTION INITLA(CONTIN, FINAL, POSDOT, POSNU, W, IW)      L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         INITIALIZES PARAMETERS USED IN CONNECTION WITH THE LANCZOS
C         STEPS. NORMALIZES THE STATRINGVECTOR (UNLESS  THE M-SCALAR
C         PRODUCT IS TOO SMALL).
C
C
C     OUTPUT PARAMETERS -
C
C         CONTIN = USED TO DECIDE IF AN OTHER LANCZOS STEP SHOULD BE
C                  TAKEN.
C         FINAL  = TRUE IF THE LAST STEP HAS BEEN TAKEN.
C         POSDOT = TRUE IF THE M-SCALAR PRODUCT IS GREATER THAN
C                  TOLPDM.
C         POSNU  = TRUE IF IT EXISTS A NU-EIGENVALUE GREATER THAN
C                  TOLZNU.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          CK, CL, COEFF, DOT, FACTOR, SRELPR, RDUMP, SECOND, ST,  R MOD.
     +          TIME, TIMTQL, TLDL, TOLBPI, TOLLDL, TOLPDM, TOLS1I,
     +          TOLZBT, TOLZNU, TOPINV, TOPM, TPRED, TSAVE, TVECOP,
     +          W(1), WRR
      INTEGER   ACTION, CNEG, CONV, COUNT, CPOS, DUMMY, ERRNO, FREE,    I***
     +          IDUMP, ITERNO, IW(1), MV, MXNEW, MXOLD, MXRST, N,
     +          NBADMU, NIL, NMXRST, NOACTN, NOR, OLCPOS, P, PFCONV,
     +          REST, RNEW, ROLD, SAVE, SAVFRE, SCPX, SOLCPX, TCONV, V,
     +          WRI, X
      LOGICAL   CONTIN, F, FINAL, MEQI, MULVEC, OPM, POSDOT, POSNU,     L***
     +          SCPROD, T, UPDATE, USEMX, WRL, ZERBET
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,
     +                  TOLZNU, TOLLDL(3)
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
C
      DATA      NOR  / 7 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      INITLA = T
C
C     ********************************************************
C     COMPUTE STARTINGVECTOR(T) * M * STARTINGVECTOR AND CHECK
C     IF IT IS GREAT ENOUGH.
C     ********************************************************
      IF(MEQI) GOTO 10
      IF(.NOT. OPM(MV + 1, V + 1, W, IW)) GOTO 8888
10    IF(.NOT. SCPROD(MV + 1, V + 1, DOT, W, IW, NOACTN, 1)) GOTO 8888
C
      POSDOT = DOT .GT. TOLPDM
      WRL(2) = POSDOT
      RDUMP = DOT
C
      IF(.NOT. POSDOT) GOTO 30
C       ***********************************
C       NORMALIZE V (AND MV IF IT IS USED).
C       ***********************************
        DOT = 1.0D0 / DSQRT(DOT)                                        D MOD.
        ACTION = SAVFRE
        IF(MEQI) ACTION = SAVE
        IF(.NOT. MULVEC(V + 1, V + 1, DOT, W, IW, ACTION, 1)) GOTO 8888
        IF(MEQI) GOTO 20
        IF(.NOT. MULVEC(MV + 1, MV + 1, DOT, W, IW, SAVE, 1)) GOTO 8888
C
C       **************************
C       MAKE SOME INITIALIZATIONS.
C       **************************
20      PFCONV = 0
        CONTIN = T
        FINAL = F
        ZERBET = F
        POSNU = F
        P = 0
        CNEG = 0
        GOTO 9999
C
30    IF(.NOT. MEQI) CALL FREEID(MV + 1)
      GOTO 9999
C
8888  INITLA = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CIO
      LOGICAL FUNCTION IO(VEC, ID, CASE, NUMEL, ADDRSS)                 L***
C
C **********************************************************************
C
C     PURPOSE - (VER = 1 OR 2)
C
C         THIS ROUTINE READS OR WRITES A VECTOR OF LENGTH N
C         (VER = 1 OR 2).
C         IF VER = 1, IT ALSO READS THE K-MATRIX.
C
C     INPUT PARAMETERS -
C
C         VEC    = IS THE VECTOR TO BE WRITTEN.
C         ID     = IS THE IDENTIFIER OF THE VECTOR.
C         NUMEL  = THE LENGTH OF THE VECTOR.
C         CASE   = 1, NOT USED IN THE TOMS VERSION.
C                  2, READ THE VECTOR.
C                  3, WRITE THE VECTOR.
C                  4, READ THE K-MATRIX.
C         ADDRSS = ADDRESS VECTOR FOR THE N-VECTORS
C
C     OUTPUT PARAMETERS -
C
C         VEC    = IS THE VECTOR TO BE WRITTEN.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      INTEGER   NUMEL                                                   I***
      DOUBLE PRECISION                                                  R INS.
     +          RDUMP, SECOND, ST, TIME, VEC(NUMEL)                     R MOD.
      INTEGER   ADDRS1, ADDRSS(1), CASE, COUNT, DAFILE, DUMMY, ERRNO,   I***
     +          ID, IDUMP, KFILE, LEFTP, LENADR, LP, MAXL, MAXREC, MV,
     +          MXNEW, MXOLD, MXRST, NBADMU, NIL, NMXRST, NOR, NREAD,
     +          NUMVEC, NWRITE, RIGHTC, RIGHTM, RIGHTP, SCPX, SOLCPX,
     +          STADEW, V, X
      LOGICAL   F, T                                                    L***
      COMMON   /STLMAD/ ADDRS1, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMEW/ LEFTP, LENADR, MAXREC, NUMVEC, RIGHTC, RIGHTM,
     +                  RIGHTP, STADEW
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      DATA      NOR  / 9 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      IO = T
C
      GOTO(9999, 100, 200, 300), CASE
C     *************
C     READ VEC(ID).
C     *************
100   IF(ADDRSS(ID) .LT. 0) GOTO 9999
C     ********************************************************
C     THE IMPLEMENTOR SHOULD REPLACE THIS FORTRAN77 STATEMENT,
C     IF NECESSARY.
C     THE STATEMENT SHOULD READ VEC OF LENGTH NUMEL (=N), WITH
C     INDEX ADDRSS(ID) ON THE DIRECT ACCESS FILE DAFILE.
C     IN CASE OF ERROR, GOTO 110.
C     ********************************************************
      READ(DAFILE, REC = ADDRSS(ID), ERR = 110) VEC
      NREAD = NREAD + 1
      GOTO 9999
C
110   IDUMP(1) = ADDRSS(ID)
      CALL ERROR(NOR, 2)
      GOTO 8888
C     **************
C     WRITE VEC(ID).
C     **************
200   IF(ADDRSS(ID) .LT. 0) GOTO 9999
      IF(.NOT. ADDRSS(ID) .EQ. 0) GOTO 202
        IF(RIGHTP .EQ. RIGHTC) RIGHTP = RIGHTM
        IF(RIGHTP .LT. LEFTP .OR. RIGHTM .LT. LEFTP) GOTO 220
        IF(.NOT. ID .LE. X) GOTO 201
          ADDRSS(ID) = RIGHTP
          RIGHTP = RIGHTP - 1
          GOTO 202
C
201     ADDRSS(ID) = LEFTP
        LEFTP = LEFTP + 1
C
C     *****************************************************************
C     THE IMPLEMENTOR SHOULD REPLACE THE FOLLOWING FORTRAN77 STATEMENT,
C     IF NECESSARY.
C     THE STATEMENT SHOULD WRITE VEC OF LENGTH NUMEL (=N), WITH INDEX
C     ADDRSS(ID) ON THE DIRECT ACCESS FILE DAFILE.
C     IN CASE OF ERROR GOTO 210.
C     *****************************************************************
202   CONTINUE
      WRITE(DAFILE, REC = ADDRSS(ID), ERR = 210) VEC
      NWRITE = NWRITE + 1
      GOTO 9999
C
210   IDUMP(1) = ADDRSS(ID)
      CALL ERROR(NOR, 3)
      GOTO 8888
220   CALL ERROR(NOR, 1)
      GOTO 8888
C     **************
C     READ K-MATRIX.
C     **************
300   REWIND KFILE
      READ(KFILE) VEC
      GOTO 9999
C
8888  IO = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CLANCZO
      LOGICAL FUNCTION LANCZO(ALPHA, BETA, TNORM, W, IW)                L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         TAKES ONE LANCZOS STEP.
C
C
C     INPUT PARAMETERS -
C
C         TNORM = MAXIMUM NORM OF THE TRIDIAGONAL T-MATRIX OF
C         DIMENSION P (.GT. 1).
C
C
C     OUTPUT PARAMETERS -
C
C         TNORM = UPDATED TNORM, I.E. NORM OF T OF ORDER P+1.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          ALPHA(1), BETA(1), DOT, FACTOR, SRELPR, RDUMP, SECOND,  R MOD.
     +          ST, ST1, TIME, TNORM, TOLBPI, TOLLDL, TOLPDM, TOLS1I,
     +          TOLZBT, TOLZNU, W(1), XX
      INTEGER   ACTION, CNEG, COUNT, CPOS, DUMMY, ERRNO, FREE, I, IDUMP,I***
     +          ITERNO, IW(1), MV, MXNEW, MXOLD, MXRST, N, NBADMU,
     +          NEXTV, NIL, NMXRST, NOACTN, NOR, OLCPOS, P, REST, RNEW,
     +          ROLD, SAVE, SAVFRE, SCPX, SOLCPX, TCONV, V, X
      LOGICAL   F, MEQI, MULVEC, OPINV, OPM, SCPROD, SUBVEC, T, USEMX,  L***
     +          ZERBET
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,
     +                  TOLZNU, TOLLDL(3)
C
      DATA      NOR  / 10 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      LANCZO = T
C
      P = P + 1
C     ********************************
C     NEXTV = ID TO THE NEXT V-VECTOR.
C     ********************************
      NEXTV = V + P + 1
      ACTION = FREE
      IF(MEQI) ACTION = NOACTN
C
      IF(.NOT. OPINV(NEXTV, MV + P, W, IW)) GOTO 8888
C
      IF(.NOT. P .GT. 1) GOTO 10
        CALL FREEID(MV + P)
        IF(.NOT. SUBVEC(NEXTV, V + P - 1, BETA(P - 1), W, IW, FREE, 1))
     +                                  GOTO 8888
C
10    IF(.NOT. SCPROD(NEXTV, MV + P, ALPHA(P), W, IW, ACTION, 1))
     +                    GOTO 8888
C
      IF(.NOT. SUBVEC(NEXTV, V + P, ALPHA(P), W, IW, FREE, 1)) GOTO 8888
C
C     ********************
C     REORTHOGONALIZATION.
C     ********************
      ST1 = SECOND()
      DO 20 I = 1, P
        COUNT(24) = COUNT(24) + 1
        IF(.NOT. SCPROD(NEXTV, MV + I, DOT, W, IW, ACTION, 2)) GOTO 8888
        IF(.NOT. SUBVEC(NEXTV, V + I, DOT, W, IW, FREE, 3)) GOTO 8888
20      CONTINUE
      TIME(24) = TIME(24) + SECOND() - ST1
C
      IF(MEQI) GOTO 25
      IF(.NOT. OPM(MV + P + 1, NEXTV, W, IW)) GOTO 8888
25    IF(.NOT. SCPROD(MV + P + 1, NEXTV, DOT, W, IW, NOACTN, 1))
     +               GOTO 8888
C
C
      IF(.NOT. P .GT. 2) GOTO 30
C       *************
C       UPDATE TNORM.
C       *************
        TNORM = DMAX1(TNORM, BETA(P - 2) + DABS(ALPHA(P - 1)) +         D MOD.
     +                       BETA(P - 1))
        GOTO 40
C
30    IF(P .EQ. 1) TNORM = DABS(ALPHA(1))                               D MOD.
      IF(P .EQ. 2) TNORM = TNORM + BETA(1)
C
40    BETA(P) =  DSIGN(DSQRT(DABS(DOT)), DOT)                           D MOD.
C
      IF(.NOT. DABS(BETA(P)) .LE. TNORM * TOLZBT) GOTO 50               D MOD.
C       ***********************************************
C       BETA(P) IS CONSIDERED TO BE ZERO. HALT LANCZOS.
C       ***********************************************
        XX = 1.0D0                                                      D MOD.
        ZERBET = T
        GOTO 60
C
50    IF(.NOT. BETA(P) .LT.  (- TNORM * TOLZBT) ) GOTO 60
C       *********************************
C       M IS CONSIDERED TO BE INDEFINITE.
C       *********************************
        RDUMP = DOT
        CALL ERROR(NOR, 1)
        GOTO 8888
C
60    IF(.NOT. ZERBET) XX = 1.0D0 / BETA(P)                             D MOD.
      ACTION = SAVFRE
      IF(MEQI) ACTION = SAVE
      IF(.NOT. MULVEC(NEXTV, NEXTV, XX, W, IW, ACTION, 1)) GOTO 8888
      IF(MEQI) GOTO 9999
      IF(.NOT. MULVEC(MV + P + 1, MV + P + 1, XX, W, IW, SAVE, 1))
     +                 GOTO 8888
C
      GOTO 9999
C
8888  LANCZO = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CLDL2
      SUBROUTINE LDL2(W, IW, X, Y, MU, N, RNEW, PROFIL, LEN)
C
C **********************************************************************
C
C     PURPOSE (VER = 2)
C
C     DUMMY ROUTINE. PLEASE SEE THE INSTALLATION GUIDE FOR
C     MORE DETAILS.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +         W(1), X(1), Y(1), MU                                     R MOD.
      INTEGER  IW(1), N, RNEW, LEN, PROFIL                              I***
C
      RETURN
      END
CLDL3
      SUBROUTINE LDL3(W, IW, MU, N, RNEW, PROFIL, LEN)
C
C **********************************************************************
C
C     PURPOSE (VER = 3)
C
C     DUMMY ROUTINE. PLEASE SEE THE INSTALLATION GUIDE FOR
C     MORE DETAILS.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +         W(1), MU                                                 R MOD.
      INTEGER  IW(1), N, RNEW, LEN, PROFIL                              I***
C
      RETURN
      END
CLDLSUB
      SUBROUTINE LDLSUB(K, M, D, MXNORM, ONENRM, MU, N, RNEW,
     +                  DIAGM, MEQI, LEN)
C
C **********************************************************************
C
C     PURPOSE - (VER = 1)
C
C         COMPUTES K-MU*M AND MAKES AN LDL(T)-DECOMPOSITION OF IT.
C         COMPUTES THE NUMBER OF NEGATIVE D(I,I)-ELEMENTS.
C
C
C     INPUT PARAMETERS -
C
C         D       THE POINTER VECTOR TO THE DIAGONAL ELEMENTS IN K.
C         DIAGM   PROFIL .LE. 2.
C         K       K MATRIX (STORED USING PROFILE STORAGE).
C         LEN     LOCAL ERROR NUMBER.
C         M       M MATRIX (STORED ACCORDING TO THE USER GUIDE).
C         MEQI    PROFIL .EQ. 3
C         N       DIMENSION OF PROBLEM.
C
C         MXNORM AND ONENRM ARE SCRATCH VECTORS OF LENGTH N.
C
C     OUTPUT PARAMETERS -
C
C         K       IS OVERWRITTEN AND CONTAINS L AND D, IN THE LDL(T)
C                 FACTORIZATION OF K-MU*M.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          DIAGEL, FACTOR, K(1), LIJ, M(1), SRELPR, MU, MXNORM(1), R MOD.
     +          NORM, ONENRM(1), RDUMP, SCPR, SIJ, TOL, TOLBPI, TOLLDL,
     +          TOLPDM, TOLS1I, TOLZBT, TOLZNU
      INTEGER   D(1), DI, DIM1, DIMI, DJ, ERRNO, I, IDUMP, IJAD, IJAD1, I***
     +          IJAD2, J, LEN, LENM1, LOOPL, MXQIQJ, N, NM1, QI, RNEW,
     +          STARTL, STARTS, TOP
      LOGICAL   DIAGM, F, MEQI, T                                       L***
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMTF/ T, F
      COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,
     +                  TOLZNU, TOLLDL(3)
C
      RNEW = 0
      LEN = 0
C
C     ********************
C     FORM K = K - MU * M.
C     ********************
      IF(.NOT. DIAGM) GOTO 40
        IF(.NOT. MEQI) GOTO 20
C         ******
C         M = I.
C         ******
          DO 10 I = 1, N
            DI = D(I)
            K(DI) = K(DI) - MU
10          CONTINUE
          GOTO 50
C
C       **********************************
C       M IS DIAGONAL, BUT NOT EQUAL TO I.
C       **********************************
20      DO 30 I = 1, N
          DI = D(I)
          K(DI) = K(DI) - MU * M(I)
30        CONTINUE
        GOTO 50
C
C     ****************************
C     M HAS THE SAME PROFILE AS K.
C     ****************************
40    CALL SUBV(K, M, MU, D(N))
C
C     *****************************
C     COMPUTE MAX NORM(K - MU * M),
C     AND THE FIRST TOLERANCE.
C     *****************************
50    CALL MAXNRM(K, MXNORM, D, N, NORM)
      IF(NORM .EQ. 0.0D0) NORM = 1.0D0                                  D MOD.
      TOL = TOLLDL(1) * NORM
C
      DO 60 I = 1, N
        MXNORM(I) = 0.0D0                                               D MOD.
        ONENRM(I) = 0.0D0                                               D MOD.
60      CONTINUE
C
      DIM1 = 0
C
C     *********************
C     FOR I = 1 TO N DO ...
C     *********************
      DO 120 I = 1, N
        TOP = DIM1 + 1
        DI = D(I)
        LENM1 = DI - TOP
C       *****************************************
C       IF (QI .LE. I - 1) THEN ... ELSE GOTO 110
C       *****************************************
        IF(.NOT. LENM1 .GE. 1) GOTO 110
          QI = I - LENM1
          J = QI
          IJAD2 = DI - 1
C         ****************************************
C         IF (QI .LT. I - 1) THEN ... ELSE GOTO 80
C         ****************************************
          IF(.NOT. LENM1 .GT. 1) GOTO 80
            IJAD1 = TOP + 1
            DIMI = DI - I
C           ******************************
C           FOR J = QI + 1 TO I - 1 DO ...
C           ******************************
            DO 70 IJAD = IJAD1, IJAD2
              J = J + 1
              DJ = D(J)
C             **********************
C             MXQIQJ = MAX (QI, QJ).
C             **********************
              MXQIQJ = MAX0(QI, J - DJ + D(J - 1) + 1)
              LOOPL = J - MXQIQJ
              IF(.NOT. LOOPL .GE. 1) GOTO 70
C               *****************************
C               COMPUTE SUM L(J, U) * S(I, U),
C               FOR U = MAX(QI, QJ) TO J - 1.
C               *****************************
                STARTL = DJ - LOOPL
                STARTS = DIMI + MXQIQJ
                K(IJAD) = K(IJAD) - SCPR(K(STARTL), K(STARTS), LOOPL)
70            CONTINUE
C
          J = QI
80        DIAGEL = K(DI)
C         **************************
C         FOR J = QI TO I - 1 DO ...
C         **************************
          DO 100 IJAD = TOP, IJAD2
            DJ = D(J)
            SIJ = K(IJAD)
            LIJ = SIJ / K(DJ)
            ONENRM(J) = ONENRM(J) + DABS(LIJ)                           D MOD.
            MXNORM(J) = DMAX1(MXNORM(J), DABS(SIJ))                     D MOD.
            DIAGEL = DIAGEL - LIJ * SIJ
            K(IJAD) = LIJ
            J = J + 1
100         CONTINUE
C
          K(DI) = DIAGEL
C
110     DIAGEL = K(DI)
C       ************
C       UPDATE RNEW.
C       ************
        IF(DIAGEL .LT. 0.0D0) RNEW = RNEW + 1                           D MOD.
        DIM1 = DI
C
C       ***************************
C       IF D(J) IS TOO SMALL, STOP.
C       ***************************
        IF(.NOT. DABS(DIAGEL) .LE. TOL) GOTO 120                        D MOD.
          LEN = -1
          IDUMP(1) = I
          IDUMP(2) = 1
          RDUMP = DIAGEL
          GOTO 9999
C
120     CONTINUE
C
      IF(N .EQ. 1) GOTO 9999
C     ****************************
C     COMPUTE THE OTHER TOLERANCE.
C     ****************************
      TOL = TOLLDL(2) * NORM
      NM1 = N - 1
C     *************************
C     FOR I = 1 TO N - 1 DO ...
C     *************************
      DO 130 I = 1, NM1
        IF(.NOT. ONENRM(I) * MXNORM(I) .GT. TOL) GOTO 130
          LEN = -1
          IDUMP(1) = I
          IDUMP(2) = 1
          RDUMP = ONENRM(I) * MXNORM(I)
          GOTO 9999
C
130     CONTINUE
C
9999  RETURN
      END
CLDLT
      LOGICAL FUNCTION LDLT(MU, RNEW, BADMU, W, IW)                     L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         DRIVER ROUTINE FOR LDL(T)-FACTORIZATION ROUTINES.
C
C
C     OUTPUT PARAMETERS -
C
C         BADMU = TRUE, IF MU IS CONSIDERED TO BE UNSUITABLE.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          MU, SMALL, W(1)                                         R MOD.
      INTEGER   ACTIVE, AD, ADDRSS, D, DAFILE, IW(1), K, KFILE, LEN, LP,I***
     +          M, MAXL, N, NOR, NREAD, NUMEL, NWRITE, READID, READK,
     +          RNEW, SAEVAL, VER, WAD1, WAD2, WRITID, PROFIL
      LOGICAL   BADMU, DIAGM, F, IO, MEQI, T                            L***
      COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMDS/ K, M, D, DIAGM
      COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMPF/ PROFIL
      COMMON   /STLMTF/ T, F
      COMMON   /STLMVR/ SMALL, VER
      COMMON   /STLMWH/ WAD1, WAD2, ACTIVE(2)
C
      DATA      NOR  / 11 /
C
      LDLT = T
C
      IF(.NOT. VER .EQ. 1) GOTO 10
        AD = D + N - 1
C       *******************************
C       FETCH K FROM SECONDARY STORAGE.
C       *******************************
        NUMEL = IW(AD)
        IF(.NOT. IO(W(K), 1, READK, NUMEL, IW(ADDRSS) )) GOTO 8888
C
C       *****************************************************
C       MAKE DECOMPOSITION. WAD1 AND WAD2 ARE ADDRESSES TO
C       N-VECTORS IN W. THEY ARE USED AS SCRATCH VECTORS AND
C       DO NOT HAVE ANY IDENTIFIERS. (WE DO NOT NEED ALLOC IN
C       THIS SPECIAL CASE).
C       *****************************************************
        CALL LDLSUB(W(K), W(M), IW(D), W(WAD1), W(WAD2),
     +              MU, N, RNEW, DIAGM, MEQI, LEN)
C
        BADMU = F
        IF(LEN .EQ. (-1)) BADMU = T
        GOTO 9999
C
10    IF(VER .EQ. 3) CALL LDL3(W, IW, MU, N, RNEW, PROFIL, LEN)
C
      IF(VER .EQ. 2) CALL LDL2(W, IW, W(WAD1), W(WAD2), MU, N,
     +                         RNEW, PROFIL, LEN)
C
      BADMU = F
      IF(LEN .EQ. (-1)) BADMU = T
      IF(.NOT. LEN .GT. 0) GOTO 9999
        CALL ERROR(100, LEN)
        CALL ERROR(100, 100)
C
8888  LDLT = F
      CALL ERROR(NOR, NOR)
C
9999  RETURN
      END
CMAXNRM
      SUBROUTINE MAXNRM(A, X, D, N, NORM)
C
C **********************************************************************
C
C     PURPOSE - (VER = 1)
C
C         COMPUTES THE MAXIMUM NORM OF A MATRIX STORED WITH THE
C         HELP OF THE D POINTER VECTOR.
C
C
C     INPUT PARAMETERS -
C
C         A = THE MATRIX.
C         D = POINTER VECTOR TO DIAGONAL ELEMENTS IN A.
C         N = DIMENSION OF A MATRIX.
C         X = A SCRATCH ARRAY OF LENGTH N.
C
C
C     OUTPUT PARAMETERS -
C
C         NORM = MAX ABS(A(I,1))+, ..., +ABS(A(I,N)).
C                 I
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          A(1), DOT, NORM, X(1)                                   R MOD.
      INTEGER   D(1), DI, DIM1, I, J, J2, LENM1, N, ROWNO, TOP          I***
C
C     *******************************************
C     FOR I = 1 TO N DO
C       X(I) = ABS(A(I, 1)) + ... + ABS(A(I, N)).
C     *******************************************
      DIM1 = 0
      DO 30 I = 1, N
        TOP = DIM1 + 1
        DI = D(I)
        LENM1 = DI - TOP
C
        DOT = 0.0D0                                                     D MOD.
        DO 10 J = TOP, DI
          DOT = DOT + DABS(A(J))                                        D MOD.
10        CONTINUE
        X(I) = DOT
C
        DIM1 = DI
        IF(.NOT. LENM1 .GT. 0) GOTO 30
          ROWNO = I - LENM1
          J2 = I - 1
          DO 20 J = ROWNO, J2
            X(J) = X(J) + DABS(A(TOP))                                  D MOD.
            TOP = TOP + 1
20          CONTINUE
C
30      CONTINUE
C
C     *********************************
C     COMPUTE MAX(X(I)), I = 1, ..., N.
C     *********************************
      NORM = X(1)
      DO 40 I = 1, N
        NORM = DMAX1(NORM, X(I))                                        D MOD.
40      CONTINUE
C
      RETURN
      END
CMCHECK
      SUBROUTINE MCHECK(M, D, N, DIAGM, MEQI, LEN, LOC)
C
C **********************************************************************
C
C     PURPOSE - (VER = 1)
C
C         THIS ROUTINE CHECKS THE DIAGONAL ELEMENTS IN M.
C
C     INPUT PARAMETERS -
C
C         D     = POINTER VECTOR TO DIAGONAL ELEMENTS OF M (IF
C                 M IS NOT DIAGONAL).
C         DIAGM = PROFIL .LE. 2
C         M     = M MATRIX STORED ACCORDING TO THE USERGUIDE.
C         MEQI  = PROFIL .EQ. 1
C         N     = DIMENSION OF M MATRIX.
C         X     = VECTOR TO BE MULTIPLIED WITH M.
C
C
C     OUTPUT PARAMETERS -
C
C         LEN  = 0, NO ERROR HAS BEEN DETECTED.
C                7, A NEGATIVE DIAGONAL ELEMENT HAS BEEN DETECTED.
C         LOC  = IS THE LOCATION OF THE NEGATIVE ELEMENT.
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          M(1)                                                    R MOD.
      INTEGER   D(1), DI, I, LEN, LOC, N                                I***
      LOGICAL   DIAGM, MEQI                                             L***
C
      LEN = 0
      IF(MEQI) GOTO 8888
C
      IF(.NOT. DIAGM) GOTO 30
        DO 20 I = 1, N
          IF(.NOT. M(I) .LT. 0.0D0) GOTO 20                             D MOD.
            LEN = 7
            LOC = I
            GOTO 8888
20        CONTINUE
C
        GOTO 8888
C
30    DO 40 I = 1, N
        DI = D(I)
        IF(.NOT. M(DI) .LT. 0.0D0) GOTO 40                              D MOD.
          LEN = 7
          LOC = I
          GOTO 8888
40      CONTINUE
C
8888  RETURN
      END
CMUL3
      SUBROUTINE MUL3(ID1, ID2, C, N, W, IW, LEN)
C
C **********************************************************************
C
C     PURPOSE (VER = 3)
C
C     DUMMY ROUTINE. PLEASE SEE THE INSTALLATION GUIDE FOR
C     MORE DETAILS.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +         C, W(1)                                                  R MOD.
      INTEGER  ID1, ID2, N, IW(1), LEN                                  I***
C
      RETURN
      END
CMULV
      SUBROUTINE MULV(X, Y, C, N)
C
C **********************************************************************
C
C     PURPOSE - (VER = 1 OR 2)
C
C         COMPUTES X=C*Y.
C
C
C     INPUT PARAMETERS -
C
C         Y = INCOMING VECTOR.
C         C = REAL CONSTANT.
C         N = DIMENSION OF THE VECTOR (NOT NECESSARILY OF THE PROBLEM).
C
C
C     OUTPUT PARAMETERS -
C
C         X = RESULTING VECTOR, C*Y.
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          C, X(1), Y(1)                                           R MOD.
      INTEGER   I, N                                                    I***
C
      DO 10 I = 1, N
        X(I) = C * Y(I)
10      CONTINUE
C
      RETURN
      END
CMULVEC
      LOGICAL FUNCTION MULVEC(ID1, ID2, C, W, IW, ACTION, CASE)         L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         DRIVER ROUTINE FOR MULV. COMPUTES VEC(ID1)=C*VEC(ID2).
C
C     INPUT PARAMETERS -
C
C         ID1    = IDENTIFIER OF THE RESULTING VECTOR.
C         ID2    = IDENTIFIER OF THE INCOMING VECTOR.
C         C      = REAL CONSTANT.
C         ACTION = DETERMINES IF VEC(ID1) OR VEC(ID2) SHOULD
C                  BE SAVED OR NOT.
C         CASE   = NOT USED IN THIS VERSION.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          C, SECOND, ST, TIME, W(1), SMALL                        R MOD.
      INTEGER   ACTION, AD1, AD2, ADDRSS, CASE, COUNT, DAFILE, DUMMY,   I***
     +          FREE, ID1, ID2, IW(1), KFILE,  LP, MAXL, N, NBADMU,
     +          NMXRST, NOACTN, NOR, NREAD, NWRITE, READID, READK,
     +          SAEVAL, SAVE, SAVFRE, WRITID, VER, LEN
      LOGICAL   ALLOC, F, IO, T                                         L***
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMVR/ SMALL, VER
C
      DATA      NOR  / 12 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      MULVEC = T
C
      IF(.NOT. VER .LE. 2) GOTO 10
C       *************
C       VER = 1 OR 2.
C       *************
        IF(.NOT. ALLOC(ID1, ID2, AD1, AD2, ACTION, W, IW(ADDRSS) ))
     +                 GOTO 8888
C
C
        CALL MULV(W(AD1), W(AD2), C, N)
C
        IF(.NOT. (ACTION .EQ. SAVE .OR. ACTION .EQ. SAVFRE)) GOTO 9999
          IF(.NOT. IO(W(AD1), ID1, WRITID, N, IW(ADDRSS)) ) GOTO 8888
          GOTO 9999
C
C     ********
C     VER = 3.
C     ********
10    CALL MUL3(ID1, ID2, C, N, W, IW, LEN)
C
      IF(.NOT. LEN .NE. 0) GOTO 9999
        CALL ERROR(103, LEN)
        CALL ERROR(103, 103)
C
8888  MULVEC = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CNEWMU
      SUBROUTINE NEWMU(LAMBDA, POINTR, MAXNU, LGEB)
C
C **********************************************************************
C
C     PURPOSE -
C
C         COMPUTES THE NEXT SHIFT. THERE ARE THREE CASES.
C
C         EIGENVALUES .GE. B HAVE CONVERGED. SET NEXT SHIFT TO B.
C
C         WE HAVE AT LEAST ONE CONVERGED EIGENVALUE GREATER
C         THE SHIFT.
C
C         WE HAVE ONLY NOT CONVERGED EIGENVALUES GREATER THAN THE SHIFT.
C
C
C     INPUT PARAMETERS -
C
C         LGEB  = T IF A CONVERGED EIGENVALUE .GE. B HAS BEEN COMPUTED AND
C                 THROWN AWAY. FALSE OTHERWISE.
C         MAXNU = MAX(NU(I)), I=1, ..., P.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          ALTMU, LAMBDA(1), MAXNU, MU, NEXTMU, OLDMU, A, B        R MOD.
      INTEGER   CNEG, CPOS, ITERNO, N, OLCPOS, P, POINT, POINTR(1),     I***
     +          REST, RNEW, ROLD, TCONV, NUMEIG, MAXW, MAXIW
      LOGICAL   USEMX, ZERBET, T, F, REACHB, USSMXR, USEDB, LGEB        L***
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMIN/ A, B, NUMEIG, MAXW, MAXIW
      COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
      COMMON   /STLMTF/ T, F
      COMMON   /STLMTS/ REACHB, USSMXR, USEDB
C
C     ***************************************************
C     CONVERGED EIGENVALUES GREATER THAN B HAVE OCCURRED.
C     SET THE NEXT SHIFT TO B.
C     ***************************************************
      IF(.NOT. LGEB) GOTO 5
        REACHB = T
        NEXTMU = B
        ALTMU  = MU + 0.1D0 * (B - MU)                                  D MOD.
        GOTO 9999
C
5     IF(.NOT. CPOS .GT. 0) GOTO 10
C       ***********
C       FIRST CASE.
C       ***********
        POINT = CNEG + CPOS
        POINT = POINTR(POINT)
C       *****************************************************
C       THE GREATEST LAMBDA IS BISECTING THE DISTANCE BETWEEN
C       THE NEXT SHIFT (NEXTMU) AND MU.
C       ALTMU IS USED IN CONNECTION WITH MXREST.
C       *****************************************************
        NEXTMU = 2.0D0 * LAMBDA(POINT) - MU                             D MOD.
        ALTMU = 1.1D0 * LAMBDA(POINT) - 0.1D0 * MU                      D MOD.
C       ***************************
C       CHECK IF WE HAVE REACHED B.
C       ***************************
        IF(.NOT. NEXTMU .GE. B) GOTO 9999
          REACHB = T
          NEXTMU = B
          GOTO 9999
C
C     ************
C     SECOND CASE.
C     ************
10    NEXTMU = MU + 0.9D0 / MAXNU                                       D MOD.
C     ***************************
C     CHECK IF WE HAVE REACHED B.
C     ***************************
      IF(.NOT. NEXTMU .GE. B) GOTO 20
        REACHB = T
        NEXTMU = B
20    ALTMU = 0.9D0 * MU + 0.1D0 * NEXTMU                               D MOD.
C
9999  RETURN
      END
CNEWPC
      SUBROUTINE NEWPC(W, IW)
C
C **********************************************************************
C
C     PURPOSE -
C
C         UPDATES POPT, COPT, CK, AND CL. GIVES A PREDICTION ON
C         HOW MUCH TIME THAT WILL BE USED FOR THE NEXT SHIFT.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          CK, CL, COEFF, FCOPT, FPOPT, K, L, SECOND, TEMP, TIMTQL,R MOD.
     +          TLDL, TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1), WEIGHT,
     +          WRR
      INTEGER   CNEG, CONV, COPT, CPOS, DUMMY, ITERNO, IW(1), LIM,      I***
     +          MXREST, N, NCOPT, NPOPT, OLCPOS, P, PFCONV, PMAX, POPT,
     +          REST, RNEW, ROLD, TCONV, WRI
      LOGICAL   MEQI, UPDATE, USEMX, WRL, ZERBET, REACHB, USSMXR, USEDB L***
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
      COMMON   /STLMTS/ REACHB, USSMXR, USEDB
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
C
      IF(.NOT. UPDATE) GOTO 9999
C
      WEIGHT = SECOND() - TSAVE
      FPOPT = P
      FCOPT = CONV
C     ***********************************************************
C     LET T(P, C) = TIME USED (ACCORDING TO OUR FORMULA) FOR ONE
C     SHIFT WHEN LANCZOS IS RUN P STEPS, AND C EIGENPAIRS HAVE
C     CONVERGED.
C     THEN WEIGHT = TIME USED SO FAR FOR THIS SHIFT / T(P, CONV).
C     WEIGHT SHOULD LIE AROUND 1.0, WHICH IS USUALLY DOES.
C     ***********************************************************
      WEIGHT = WEIGHT / (COEFF(1) + FPOPT * (COEFF(2) + FCOPT * COEFF(4)
     +         + FPOPT * (TVECOP + FPOPT * TIMTQL)) + FCOPT * COEFF(3))
C
      IF(.NOT. CONV .GT. 0) GOTO 10
C       *********************************************
C       COMPUTE NEW K AND L.
C       (P = K * NUMBER OF CONVERGED EIGENPAIRS + L).
C       *********************************************
        L = PFCONV - 1
        TEMP = P
        K = TEMP - L
        TEMP = CONV
        K = K / TEMP
C
C       **************************
C       COMPUTE NEW POPT AND COPT.
C       **************************
        CALL HALF(K, L, NPOPT, NCOPT, LIM)
C
        WRR(4) = WEIGHT
        CALL WRINFO(25, 1, W, IW)
        WRR(1) = L
        WRR(2) = K
C
C       *************************************
C       CHECK THAT THE VALUES ARE REASONABLE.
C       *************************************
        IF(.NOT. (2 * NPOPT .GE. POPT .AND. NPOPT .LE. 2 * POPT
     +            .AND. 2 * NCOPT .GE. COPT .AND. NCOPT .LE. 2 * COPT
     +            .AND. LIM .GE. 4)) GOTO 10
C
C         **************************************
C         YES, THEY ARE. REPLACE THE OLD VALUES.
C         **************************************
          POPT = NPOPT
          COPT = NCOPT
          IF(.NOT. USSMXR) MXREST = COPT
          CK = K
          CL = L
C
10    FPOPT = POPT
      FCOPT = COPT
C     ***********************************************************
C     PREDICTED TIME FOR THE NEXT SHIFT = WEIGHT * T(POPT, COPT).
C     ***********************************************************
      TPRED = WEIGHT * (COEFF(1) + FPOPT * (COEFF(2) + FCOPT * COEFF(4)
     +        + FPOPT * (TVECOP + FPOPT * TIMTQL)) + FCOPT * COEFF(3))
C
      CALL WRINFO(25, 2, W, IW)
C
9999  RETURN
      END
COPINV
      LOGICAL FUNCTION OPINV(ID1, ID2, W, IW)                           L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         DRIVER ROUTINE FOR THE SOLVER ROUTINE. SOLVES
C         (LDL(T))*VEC(ID)=VEC(ID2).
C
C     INPUT PARAMETERS -
C
C         ID1    = IDENTIFIER OF THE RESULTING VECTOR.
C         ID2    = IDENTIFIER OF THE INCOMING VECTOR.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          SECOND, SMALL, ST, TIME, W(1)                           R MOD.
      INTEGER   AD1, AD2, ADDRSS, COUNT, D, DAFILE, DUMMY, FREE, ID1,   I***
     +          ID2, IW(1), K, KFILE, LEN, LP, M, MAXL, N, NBADMU,
     +          NMXRST, NOACTN, NOR, NREAD, NWRITE, READID, READK,
     +          SAEVAL, SAVE, SAVFRE, VER, WRITID
      LOGICAL   ALLOC, DIAGM, F, MEQI, T                                L***
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMDS/ K, M, D, DIAGM
      COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMVR/ SMALL, VER
C
      DATA      NOR  / 13 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      OPINV = T
C
      IF(VER .EQ. 3) GOTO 10
      IF(.NOT. ALLOC(ID1, ID2, AD1, AD2, NOACTN, W, IW(ADDRSS)))
     +            GOTO 8888
C
      IF(VER .EQ. 2) GOTO 20
      CALL SOLVE(W(AD1), W(AD2), W(K), IW(D), N)
      GOTO 9999
C
10    CALL SOL3(ID1, ID2, W, IW, N, LEN)
20    IF(VER .EQ. 2) CALL SOL2(W(AD1), W(AD2),  W, IW, N, LEN)
      IF(.NOT. LEN .GT. 0) GOTO 9999
        CALL ERROR(102, LEN)
        CALL ERROR(102, 102)
C
8888  OPINV = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
COPM
      LOGICAL FUNCTION OPM(ID1, ID2, W, IW)                             L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         DRIVER ROUTINE FOR THE M-MULTIPLICATION ROUTINE. COMPUTES
C         VEC(ID1)=M*VEC(ID2).
C
C     INPUT PARAMETERS -
C
C         ID1    = IDENTIFIER OF THE RESULTING VECTOR.
C         ID2    = IDENTIFIER OF THE INCOMING VECTOR.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          SECOND, SMALL, ST, TIME, W(1)                           R MOD.
      INTEGER   AD1, AD2, ADDRSS, COUNT, D, DAFILE, DUMMY, FREE, ID1,   I***
     +          ID2, IW(1), K, KFILE, LEN, LP, M, MAXL, N, NBADMU,
     +          NMXRST, NOACTN, NOR, NREAD, NWRITE, READID, READK,
     +          SAEVAL, SAVE, SAVFRE, VER, WRITID, PROFIL
      LOGICAL   ALLOC, DIAGM, F, T                                      L***
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMDS/ K, M, D, DIAGM
      COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMPF/ PROFIL
      COMMON   /STLMTF/ T, F
      COMMON   /STLMVR/ SMALL, VER
C
      DATA      NOR  / 14 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      OPM = T
C
      IF(VER .EQ. 3) GOTO 10
      IF(.NOT. ALLOC(ID1, ID2, AD1, AD2, NOACTN, W, IW(ADDRSS)))
     +                GOTO 8888
C
      IF(VER .EQ. 2) GOTO 20
      CALL OPMSUB(W(AD1), W(AD2), W(M), IW(D), N, DIAGM)
      GOTO 9999
C
10    CALL OPM3(ID1, ID2, W, IW, N, PROFIL, LEN)
20    IF(VER .EQ. 2) CALL OPM2(W(AD1), W(AD2), W, IW, N, PROFIL, LEN)
C
      IF(.NOT. LEN .GT. 0) GOTO 9999
        CALL ERROR(101, LEN)
        CALL ERROR(101, 101)
C
8888  OPM = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
COPM2
      SUBROUTINE OPM2(MX, X, W, IW, N, PROFIL, LEN)
C
C **********************************************************************
C
C     PURPOSE (VER = 2)
C
C     DUMMY ROUTINE. PLEASE SEE THE INSTALLATION GUIDE FOR
C     MORE DETAILS.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +         MX(1), X(1), W(1)                                        R MOD.
      INTEGER  IW(1), N, LEN, PROFIL                                    I***
C
      RETURN
      END
COPM3
      SUBROUTINE OPM3(IDMX, IDX, W, IW, N, PROFIL, LEN)
C
C **********************************************************************
C
C     PURPOSE (VER = 3)
C
C     DUMMY ROUTINE. PLEASE SEE THE INSTALLATION GUIDE FOR
C     MORE DETAILS.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +         W(1)                                                     R MOD.
      INTEGER  IW(1), N, LEN, IDMX, IDX, PROFIL                         I***
C
      RETURN
      END
COPMSUB
      SUBROUTINE OPMSUB(MX, X, M, D, N, DIAGM)
C
C **********************************************************************
C
C     PURPOSE - (VER = 1)
C
C         MX=M*X IS COMPUTED.
C
C
C     INPUT PARAMETERS -
C
C         D     = POINTER VECTOR TO DIAGONAL ELEMENTS OF M (IF
C                 M IS NOT DIAGONAL).
C         DIAGM = PROFIL .LE. 2
C         M     = M MATRIX STORED ACCORDING TO THE USER GUIDE.
C         N     = DIMENSION OF M MATRIX.
C         X     = VECTOR TO BE MULTIPLIED BY M.
C
C
C     OUTPUT PARAMETERS -
C
C         MX = RESULTING VECTOR.
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          M(1), MX(1), SCPR, X(1)                                 R MOD.
      INTEGER   D(1), DI, DIM1, I, LENM1, N, ROWNO, TOP                 I***
      LOGICAL   DIAGM                                                   L***
C
      IF(.NOT. DIAGM) GOTO 40
C       ***********************************
C       M IS DIAGONAL (BUT NOT EQUAL TO I).
C       ***********************************
        DO 30 I = 1, N
          MX(I) = M(I) * X(I)
30        CONTINUE
        GOTO 9999
C
C     ******************
C     M IS NON DIAGONAL.
C     ******************
40    DIM1 = 0
C
      DO 50 I = 1, N
        TOP = DIM1 + 1
        DI = D(I)
        LENM1 = DI - TOP
        ROWNO = I - LENM1
        MX(I) = SCPR(X(ROWNO), M(TOP), LENM1 + 1)
        IF(LENM1 .GT. 0) CALL SUBV(MX(ROWNO), M(TOP), - X(I), LENM1)
        DIM1 = DI
50      CONTINUE
C
9999  RETURN
      END
CPREPSV
      LOGICAL FUNCTION PREPSV(IDX, IDMX, NUMBER, W, IW)                 L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         PREPARES THE STARTINGVECTOR, I.E. ORTHOGONALIZES IT AGAINST
C         SOME ALREADY COMPUTED EIGENVECTORS.
C
C
C
C     INPUT PARAMETERS -
C
C         IDX    = IDX+1, ..., IDX+NUMBER ARE THE IDENTIFIERS FOR THE
C                  X-VECTORS.
C         IDMX   = IDMX+1, ..., IDMX+NUMBER ARE THE IDENTIFIERS
C                  FOR THE MX-VECTORS (IF USED).
C         NUMBER = IS THE NUMBER OF VECTORS.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          DOT, SECOND, ST, TIME, W(1)                             R MOD.
      INTEGER   CNEG, COUNT, CPOS, DUMMY, FREE, I, IDMX, IDX, ITERNO,   I***
     +          IW(1), MV, MXNEW, MXOLD, MXRST, N, NBADMU, NIL, NMXRST,
     +          NOACTN, NOR, NUMBER, OLCPOS, P, REST, RNEW, ROLD, SAVE,
     +          SAVFRE, SCPX, SOLCPX, TCONV, V, X
      LOGICAL   F, MEQI, OPM, SCPROD, SUBVEC, T, USEMX, ZERBET          L***
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
C
      DATA      NOR  / 15 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      PREPSV = T
C
      IF(NUMBER .LE. 0) GOTO 9999
C
      IF(.NOT. USEMX) GOTO 15
C       ************************
C       THE MX-VECTORS ARE USED.
C       ************************
        DO 10 I = 1, NUMBER
C         *****************************
C         V1 = V1 - (V1(T) * (MX)) * X.
C         *****************************
          IF(.NOT. SCPROD(V + 1, IDMX + I, DOT, W, IW, FREE, 2))
     +      GOTO 8888
C
          IF(.NOT. SUBVEC(V + 1, IDX + I, DOT, W, IW, FREE, 3))
     +      GOTO 8888
C
10        CONTINUE
        GOTO 9999
C
15    IF(.NOT. MEQI) GOTO 30
C       ******
C       M = I.
C       ******
        DO 20 I = 1, NUMBER
C         **************************
C         V1 = V1 - (V1(T) * X) * X.
C         **************************
          IF(.NOT. SCPROD(V + 1, IDX + I, DOT, W, IW, NOACTN, 2))
     +      GOTO 8888
C
          IF(.NOT. SUBVEC(V + 1, IDX + I, DOT, W, IW, FREE, 3))
     +      GOTO 8888
C
20      CONTINUE
      GOTO 9999
C
C     ********************
C     M IS NOT EQUAL TO I.
C     ********************
30    CALL FREEID(V + 1)
      DO 40 I = 1, NUMBER
C       ****************************************************
C       V2 = M * X, (V2 = V + 2 IS USED AS A SCRATCH VECTOR)
C       V1 = V1 - (V1(T) * V2) * X.
C       ****************************************************
        IF(.NOT. OPM(V + 2, IDX + I, W, IW)) GOTO 8888
        CALL FREEID(IDX + I)
        IF(.NOT. SCPROD(V + 2, V + 1, DOT, W, IW, NOACTN, 2)) GOTO 8888
        CALL FREEID(V + 2)
        IF(.NOT. SUBVEC(V + 1, IDX + I, DOT, W, IW, SAVFRE, 2))
     +      GOTO 8888
C
        CALL FREEID(IDX + I)
40      CONTINUE
      GOTO 9999
C
8888  PREPSV = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CRAN3
      SUBROUTINE RAN3(ID, N, W, IW, LEN)
C
C **********************************************************************
C
C     PURPOSE (VER = 3)
C
C     DUMMY ROUTINE. PLEASE SEE THE INSTALLATION GUIDE FOR
C     MORE DETAILS.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +         W(1)                                                     R MOD.
      INTEGER  ID, N, IW(1), LEN                                        I***
C
      RETURN
      END
CRANDVC
      LOGICAL FUNCTION RANDVC(ID, W, IW, SAVE1)                         L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         VEC(ID) = (K - MU * M)**(-1) * M * RANDOMVECTOR.
C
C
C     INPUT PARAMETERS -
C
C         ID    = IDENTIFIER OF THE RESULTING VECTOR.
C         SAVE1 = TRUE, VEC(ID) IS STORED AFTER INITIALIZATION.
C                 FALSE, VEC(ID) IS NOT STORED.
C                 SAVE1 IS ONLY USED IF VER = 1 OR 2.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          RANFX, SECOND, ST, TIME, W(1), SMALL                     R MOD.
      INTEGER   AD1, AD2, ADDRSS, ADSAVE, COUNT, DAFILE, DUMMY, FREE, I,I***
     +          ID, ID1, IW(1), KFILE, LP, MAXL, MV, MXNEW, MXOLD,
     +          MXRST, N, NBADMU, NIL, NMXRST, NOACTN, NOR, NREAD,
     +          NWRITE, READID, READK, SAEVAL, SAVE, SAVFRE, SCPX,
     +          SOLCPX, V, WRITID, X, VER, ID2, LEN
      LOGICAL   ALLOC, F, IO, MEQI, OPM, SAVE1, T, OPINV                L***
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMVR/ SMALL, VER
C
      DATA      NOR  / 16 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      RANDVC = T
C
      ID1 = ID
      ID2 = V + 2
      IF(.NOT. MEQI) GOTO 5
        ID2 = ID
        ID1 = V + 2
C
5     IF(.NOT. VER .LE. 2) GOTO 20
C     *************
C     VER = 1 OR 2.
C     *************
        IF(.NOT. ALLOC(ID1, NIL, AD1, AD2, NOACTN, W, IW(ADDRSS) ))
     +      GOTO 8888
C
        ADSAVE = AD1
        AD1    = AD1 - 1
C
        DO 10 I = 1, N
          AD1 = AD1 + 1
C         ******************************************
C         NOTE. RANF SHOULD BE SUPPLIED BY THE USER.
C         ******************************************
          W(AD1) = RANFX(DUMMY)
10        CONTINUE
C
        GOTO 30
C
C     ********
C     VER = 3.
C     ********
20    CALL RAN3(ID1, N, W, IW, LEN)
      IF(LEN .NE. 0) GOTO 8888
C
30    IF(.NOT. MEQI) GOTO 40
C       *****************
C       TRANSFORM VECTOR.
C       *****************
        IF(.NOT. OPINV(ID2, ID1, W, IW)) GOTO 8888
        CALL FREEID(ID1)
        GOTO 9999
C
C     *****************
C     TRANSFORM VECTOR.
C     *****************
40    IF(.NOT. OPM(ID2, ID1, W, IW)) GOTO 8888
      IF(.NOT. OPINV(ID1, ID2, W, IW)) GOTO 8888
      CALL FREEID(ID2)
C
      IF(.NOT. VER .LE. 2) GOTO 9999
        IF(.NOT. SAVE1) GOTO 9999
C         ************
C         SAVE VECTOR.
C         ************
          IF(.NOT. IO(W(ADSAVE), ID, WRITID, N, IW(ADDRSS) )) GOTO 8888
C
          GOTO 9999
C
8888  IF(.NOT. VER .EQ. 3) GOTO 8889
        CALL ERROR(104, LEN)
        CALL ERROR(104, 104)
C
8889  RANDVC = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CREPL
      SUBROUTINE REPL
C
C **********************************************************************
C
C     PURPOSE -
C
C         UPDATES CERTAIN CONTROL VARIABLES.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          ALTMU, MU, NEXTMU, OLDMU                                R MOD.
      INTEGER   BUFF, CNEG, CPOS, ITERNO, MV, MXNEW, MXOLD, MXRST, N,   I***
     +          NIL, OLCPOS, P, REST, RNEW, ROLD, SCPX, SOLCPX, TCONV,
     +          V, X
      LOGICAL   USEMX, ZERBET                                           L***
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
C
      OLDMU = MU
      MU = NEXTMU
      ROLD = RNEW
      OLCPOS = CPOS
C     ***********************
C     INTERCHANGE MX-BUFFERS.
C     ***********************
      BUFF = MXNEW
      MXNEW = MXOLD
      MXOLD = BUFF
      SOLCPX = SCPX
C
      RETURN
      END
CROTATE
      SUBROUTINE ROTATE(S, POINTR, FINAL, PMAX, NUMDEL)
C
C **********************************************************************
C
C     PURPOSE -
C
C         GIVEN TWO S-VECTORS CORRESPONDING TO NUMERICALLY INDENTICAL
C         NU-EIGENVALUES THIS ROUTINE CONSTRUCTS ONE GOOD VECTOR. THIS
C         IS MADE WITH A ROTATION MATRIX.
C
C
C     INPUT PARAMETERS -
C
C         POINTR = A PART OF THE POINTR VECTOR.
C         FINAL  = TRUE IF THE LAST LANCZOS STEP HAS BEEN TAKEN.
C
C
C     OUTPUT PARAMETERS -
C
C         NUMDEL = UPDATED VALUE ON THE NUMBER OF DELETED AND
C                  CONVERGED VECTORS.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C
C **********************************************************************
C
      INTEGER   PMAX                                                    I***
      DOUBLE PRECISION                                                  R INS.
     +          H1, H2, Q, S(PMAX, 2)                                   R MOD.
      INTEGER   BAD, CNEG, CPOS, GOOD, I, ITERNO, N, NUMDEL,            I***
     +          OLCPOS, P, POINTR(2), REST, RNEW, ROLD, STEP, TCONV
      LOGICAL   FINAL, USEMX, ZERBET                                    L***
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
C
C     *****************************************************
C     IF( ... ) THEN WE HAVE DELETED A CONVERGED EIGENPAIR.
C     *****************************************************
      IF(POINTR(1) .GT. 0 .AND. POINTR(2) .GT. 0) NUMDEL = NUMDEL + 1
      I = 1
C     *******************************************************
C     THE BETTER VECTOR (GOOD) HAS A GREATER FIRST COMPONENT.
C     *******************************************************
      GOOD = IABS(POINTR(1))
      BAD = IABS(POINTR(2))
      IF(.NOT. DABS(S(1, GOOD)) .LT. DABS(S(1, BAD))) GOTO 10           D MOD.
        I = GOOD
        GOOD = BAD
        BAD = I
        I = 2
10    POINTR(I) = GOOD
      I = 3 - I
C     *********************
C     CLEAR THE BAD POINTR.
C     *********************
      POINTR(I) = 0
C
C     ******************************
C     CONSTRUCT THE ROTATION MATRIX.
C     ******************************
      Q = S(1, BAD) / S(1, GOOD)
      H2 = DSQRT(1.0D0  + Q * Q)                                        D MOD.
      H1 = 1.0D0  / H2                                                  D MOD.
      H2 = Q / H2
C
C     ******************************************************
C     IF NOT FINAL, CHANGE ONLY THE TOP AND BOTTOM ELEMENTS.
C     ******************************************************
      STEP = P - 1
      IF(FINAL) STEP = 1
C
C     *******
C     ROTATE.
C     *******
      DO 20 I = 1, P, STEP
        S(I, GOOD) = H1 * S(I, GOOD) + H2 * S(I, BAD)
20      CONTINUE
C
      RETURN
      END
CSAVEXL
      LOGICAL FUNCTION SAVEXL(BETAPI, LAMBDA, NU, POINTR, S, MXNEG,     L***
     +                        PMAX, CPALSO, W, IW)
C
C **********************************************************************
C
C     PURPOSE -
C
C         COMPUTES AND STORES EIGENVECTORS (AND MX-VECTORS IF USEMX).
C         SAVES EIGENVALUSES.
C
C
C     INPUT PARAMETERS -
C
C         MXNEG  = IF THE MX(CNEG) VECTORS ARE COMPUTED THEY GET
C                  IDENTIFIERS, MXNEG+1, ..., MXNEG+CNEG.
C         CPALSO = FALSE, COMPUTE AND SAVE CNEG-PAIRS. IF TRUE, HANDLE
C                  THE CPOS-PAIRS ALSO.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      INTEGER   PMAX                                                    I***
      DOUBLE PRECISION                                                  R INS.
     +          ADD, BETAPI(1), LAMBDA(1), NU(1), RDUMP, S(PMAX,1),     R MOD.
     +          SECOND, ST, TIME, W(1)
      INTEGER   ADDRSS, CNEG, COUNT, CPOS, DAFILE, DUMMY, ERRNO, I, I1, I***
     +          I2, IDUMP, ITERNO, IW(1), J, J1, KFILE, LP, MAXL, MV,
     +          MXNEG, MXNEW, MXOLD, MXRST, MXS, N, NALSO, NBADMU, NIL,
     +          NMXRST, NOR, NREAD, NWRITE, OLCPOS, P, POINT, POINTR(1),
     +          READID, READK, REST, RNEW, ROLD, SAEVAL, SCPX,
     +          SOLCPX, TCONV, V, WRITID, X
      LOGICAL   CPALSO, F, T, TRANSF, USEMX, ZERBET                     L***
      COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, NALSO
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
C
      DATA      NOR  / 17 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      SAVEXL = T
C
C     ***********************
C     PREPARE FOR CNEG PAIRS.
C     ***********************
      MXS = MXNEG
      I1 = 1
      I2 = CNEG
      IF(CNEG .EQ. 0) GOTO 30
C
10    J = 1
C     *************************************
C     CHECK IF WE SHALL COMPUTE MX-VECTORS.
C     *************************************
      IF(.NOT. USEMX) MXS = NIL
C
      DO 20 I = I1, I2
        POINT = POINTR(I)
C       ********************
C       UPDATE TCONV.
C       COMPUTE EIGENVECTOR.
C       SAVE EIGENVALUE.
C       ********************
        TCONV = TCONV + 1
C       **********************************
C       CHECK IF ROOM TO STORE EIGENVALUE.
C       **********************************
        IF(.NOT. TCONV .GT. MAXL) GOTO 15
          TCONV = TCONV - 1
          RDUMP = LAMBDA(POINT)
          CALL ERROR(NOR, 1)
          GOTO 8888
C
C       ****************
C       SAVE EIGENVALUE.
C       ****************
15      J1 = LP + TCONV - 1
        W(J1) = LAMBDA(POINT)
C
        ADD = -BETAPI(POINT)
        IF(ZERBET) ADD = -S(P, POINT) / NU(POINT)
        IF(.NOT. TRANSF(ADD, S(1, POINT), X + TCONV, V, W, IW))
     +                   GOTO 8888
C
        IF(.NOT. MXS .NE. NIL) GOTO 20
C         ******************
C         COMPUTE MX-VECTOR.
C         ******************
          IF(.NOT. TRANSF(ADD, S(1, POINT), MXS + J, MV, W, IW))
     +                     GOTO 8888
C
C
          J = J + 1
20      CONTINUE
C
C     ***************
C     CHECK IF READY.
C     ***************
30    IF(.NOT. (CPALSO .AND. I2 .EQ. CNEG .AND. CPOS .GT. 0)) GOTO 9999
C       *************************
C       PREPARE FOR CPOS VECTORS.
C       *************************
        I1 = CNEG + 1
        I2 = CNEG + CPOS
        MXS = MXNEW
        GOTO 10
C
8888  SAVEXL = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CSCP3
      SUBROUTINE SCP3(ID1, ID2, C, N, W, IW, LEN)
C
C **********************************************************************
C
C     PURPOSE (VER = 3)
C
C     DUMMY ROUTINE. PLEASE SEE THE INSTALLATION GUIDE FOR
C     MORE DETAILS.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +         C, W(1)                                                  R MOD.
      INTEGER  ID1, ID2, N, IW(1), LEN                                  I***
C
      RETURN
      END
CSCPR
      DOUBLE PRECISION                                                  R INS.
     +     FUNCTION SCPR(X, Y, N)                                       R MOD.
C
C **********************************************************************
C
C     PURPOSE - (VER = 1 OR 2)
C
C         COMPUTES SCPR=X(T)*Y.
C
C
C     INPUT PARAMETERS -
C
C         X AND Y = VECTORS PARTICIPATING IN THE COMPUTATION.
C         N       = THE DIMENSION OF THE VECTORS (NOT NECESSARILY OF THE
C                   PROBLEM).
C
C
C     OUTPUT PARAMETERS -
C
C         SCPR = THE SCALARPRODUCT BETWEEN THE VECTORS.
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          X(1), Y(1)                                              R MOD.
      INTEGER   I, N                                                    I***
C
      SCPR = 0.0D0                                                      D MOD.
      DO 10 I = 1, N
        SCPR = SCPR + X(I) * Y(I)
10      CONTINUE
C
      RETURN
      END
CSCPROD
      LOGICAL FUNCTION SCPROD(ID1, ID2, DOT, W, IW, ACTION, CASE)       L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         DRIVER ROUTINE  FOR THE SCALAR PRODUCT ROUTINE. COMPUTES
C         DOT=VEC(ID1)(T)*VEC(ID2).
C
C     INPUT PARAMETERS -
C
C         ID1    = IDENTIFIER OF INCOMING VECTOR.
C         ID2    = IDENTIFIER OF INCOMING VECTOR.
C         ACTION = DETERMINES IF VEC(ID1) OR VEC(ID2) SHOULD
C                  BE SAVED OR NOT.
C         CASE   = NOT USED IN THIS VERSION.
C
C     OUTPUT PARAMETERS -
C
C         DOT   = THE SCALARPRODUCT.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          DOT, SCPR, SECOND, ST, TIME, W(1), SMALL, C             R MOD.
      INTEGER   ACTION, AD1, AD2, ADDRSS, CASE, COUNT, DAFILE, DUMMY,   I***
     +          FREE, ID1, ID2, IW(1), KFILE, LP, MAXL, N, NBADMU,
     +          NMXRST, NOACTN, NOR, NREAD, NWRITE, READID, READK,
     +          SAEVAL, SAVE, SAVFRE, WRITID, VER, LEN
      LOGICAL   ALLOC, F, IO, T                                         L***
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMVR/ SMALL, VER
C
      DATA      NOR  / 18 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      SCPROD = T
C
      IF(.NOT. VER .LE. 2) GOTO 10
C       *************
C       VER = 1 OR 2.
C       *************
        IF(.NOT. ALLOC(ID1, ID2, AD1, AD2, ACTION, W, IW(ADDRSS) ))
     +                 GOTO 8888
C
        DOT = SCPR(W(AD1), W(AD2), N)
C
        IF(.NOT. (ACTION .EQ. SAVE .OR. ACTION .EQ. SAVFRE)) GOTO 9999
          IF(.NOT. IO(W(AD1), ID1, WRITID, N, IW(ADDRSS) )) GOTO 8888
          GOTO 9999
C
C     ********
C     VER = 3.
C     ********
10    CALL SCP3(ID1, ID2, C, N, W, IW, LEN)
C
      IF(LEN .EQ. 0) DOT = C
      IF(.NOT. LEN .NE. 0) GOTO 9999
        CALL ERROR(105, LEN)
        CALL ERROR(105, 105)
C
8888  SCPROD = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CSELDOM
      LOGICAL FUNCTION SELDOM(W, IW)                                    L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         HANDLES THE CASE WHEN REST IS POSITIVE AFTER ALWAYS, I.E.
C         THERE ARE MISSING EIGENVALUES IN THE INTERVAL. THE ALGORITHM
C         IS AS FOLLOWS,
C
C               10 ORTHOGONALIZE THE STARTINGVECTOR AGAINST THE VECTORS
C                  THAT HAVE CONVERGED, AND WHOSE EIGENVALUES LIE
C                  IN (OLDMU,MU) (AND AGAINST THOSE TO THE RIGHT OF
C                  MU THE FIRST TIME).
C
C                  RUN LANCZOS AS USUAL.
C
C                  SAVE THE NEW VECTORS WITH LAMBDA IN THE INTERVAL.
C
C                  IF STILL MISSING, GOTO 10.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          RDUMP, SECOND, ST, TIME, TNORM, W(1), WRR               R MOD.
      INTEGER   ALPHA, BETA, BETAPI, CNEG, COPT, COUNT, CPOS, DUMMY,    I***
     +          ERRNO, IDUMP, ITERNO, IW(1), LAMBDA, MV, MXNEG, MXNEW,
     +          MXOLD, MXREST, MXRST, N, NBADMU, NIL, NMXRST, NOR, NU,
     +          NUMR, OLCPOS, P, PMAX, POINTR, POPT, REST, RNEW, ROLD,
     +          S, SACNEG, SACPOS, SCNEGX, SCPOSX, SCPX, SCR, SOLCPX,
     +          SRESTX, TCONV, V, WRI, X
      LOGICAL   CONTIN, CPALSO, F, FINAL, INITLA, LANCZO, MEQI, POSDOT, L***
     +          POSNU, PREPSV, RANDVC, SAVEXL, T, TRIDIG, USEMX, VISITL,
     +          WRL, ZERBET
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
      COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
C
      DATA      NOR  / 19 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      SELDOM = T
C
C     ************************************************************
C     SET UP POINTER INFORMATION.
C       SCNEGX = START OF CNEG X-VECTORS.
C       SCPOSX =          CPOS
C       SRESTX =          REST
C
C     SAVE CNEG AND CPOS, SINCE THEY WOULD BE DESTROYED OTHERWISE.
C     NUMR = NUMBER OF REST VECTORS.
C     ************************************************************
      SCNEGX = X + TCONV
      SCPOSX = SCNEGX + CNEG
      SRESTX = SCPOSX + CPOS
      SACNEG = CNEG
      SACPOS = CPOS
      NUMR = 0
C
      CPALSO = T
C     ********************************************************
C     MXNEG = START OF CNEG MX-VECTORS.
C     NOTE, THERE IS ROOM FOR THE CPOS MX-VECTORS TO THE LEFT.
C     ********************************************************
      MXNEG = MXNEW + CPOS
      SCPX = SCPOSX
C     *********************************
C     COMPUTE AND STORE THE EIGENPAIRS.
C     *********************************
      IF(.NOT. SAVEXL(W(BETAPI), W(LAMBDA), W(NU), IW(POINTR), W(S),
     +                MXNEG, PMAX, CPALSO, W, IW)) GOTO 8888
C
      VISITL = F
C
C     *************************************************
C     COMPUTE RANDOMVECTOR AND ORTHOGONALIZE IT AGAINST
C       CPOS
C       OLD CPOS
C       CNEG
C       REST X-VECTORS.
C     *************************************************
10    IF(.NOT. RANDVC(V + 1, W, IW, .NOT. USEMX .AND. .NOT. MEQI))
     +                                GOTO 8888
C
      IF(.NOT. PREPSV(SCPOSX, MXNEW, SACPOS, W, IW)) GOTO 8888
      IF(.NOT. PREPSV(SOLCPX, MXOLD, OLCPOS, W, IW)) GOTO 8888
      IF(.NOT. PREPSV(SCNEGX, MXNEW + SACPOS, SACNEG, W, IW)) GOTO 8888
      IF(.NOT. PREPSV(SRESTX, MXRST, NUMR, W, IW)) GOTO 8888
      IF(.NOT. INITLA(CONTIN, FINAL, POSDOT, POSNU, W, IW)) GOTO 8888
C
      IF(POSDOT) GOTO 30
C       ************************************
C       M INDEFINITE, OR BAD STARTINGVECTOR.
C       ************************************
        IF(.NOT. VISITL) GOTO 20
          CALL ERROR(NOR, 1)
          GOTO 8888
20      VISITL = T
        GOTO 10
C
C     *******************************
C     RUN LANCZOS. CHECK CONVERGENCE.
C     *******************************
30    IF(.NOT. (CONTIN .AND. P .LT. POPT .AND. CNEG .LT. REST)) GOTO 40
        IF(.NOT. LANCZO(W(ALPHA), W(BETA), TNORM, W, IW)) GOTO 8888
        IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, F, W, IW)) GOTO 8888
        GOTO 30
C
40    IF(.NOT. (CONTIN .AND. CNEG .EQ. 0)) GOTO 50
        IF(.NOT. LANCZO(W(ALPHA), W(BETA), TNORM, W, IW)) GOTO 8888
        IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, F, W, IW)) GOTO 8888
        GOTO 40
C
50    IF(.NOT. (.NOT. CONTIN .AND. CNEG .EQ. 0)) GOTO 60
C       *************************************
C       CAN NOT CONTINUE, AND NO CONVERGENCE.
C       *************************************
        CALL ERROR(NOR, 2)
        GOTO 8888
C
C     **********************
C     NO MORE LANCZOS STEPS.
C     **********************
60    FINAL = T
      IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, F, W, IW)) GOTO 8888
      CALL FREEID(MV + P + 1)
C
C     **************************************************
C     NEW REST. REST .LT. 0 SHOULD NEVER HAPPEN, BUT ...
C     **************************************************
      REST = REST - CNEG
      IF(.NOT. REST .LT. 0) GOTO 70
        CALL ERROR(NOR, 3)
        GOTO 8888
C
70    CPALSO = F
C     *******************************************************
C     MXNEG = START OF REST MX-VECTORS.
C     IF REST = 0, WE ARE READY, AND DO NOT HAVE TO SAVE ANY.
C     COMPUTE AND SAVE REST EIGENPAIRS. NOTE CPALSO = F.
C     UPDATE NUMR.
C     *******************************************************
      MXNEG = MXRST + NUMR
      IF(REST .EQ. 0) MXNEG = NIL
C
      IF(.NOT. SAVEXL(W(BETAPI), W(LAMBDA), W(NU), IW(POINTR), W(S),
     +                MXNEG, PMAX, CPALSO, W, IW)) GOTO 8888
C
      NUMR = NUMR + CNEG
      WRI(1) = NUMR
      WRR(1) = TNORM
      CALL WRINFO(NOR, 1, W, IW)
      IF(REST .GT. 0) GOTO 10
      GOTO 9999
C
8888  SELDOM = F
      CALL ERROR(NOR, NOR)
C
9999  CNEG = SACNEG
      CPOS = SACPOS
      TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CSOL2
      SUBROUTINE SOL2(X, B, W, IW, N, LEN)
C
C **********************************************************************
C
C     PURPOSE (VER = 2)
C
C     DUMMY ROUTINE. PLEASE SEE THE INSTALLATION GUIDE FOR
C     MORE DETAILS.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +         X(1), B(1), W(1)                                         R MOD.
      INTEGER  IW(1), N, LEN                                            I***
C
      RETURN
      END
CSOL3
      SUBROUTINE SOL3(IDX, IDB, W, IW, N, LEN)
C
C **********************************************************************
C
C     PURPOSE (VER = 3)
C
C     DUMMY ROUTINE. PLEASE SEE THE INSTALLATION GUIDE FOR
C     MORE DETAILS.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +         W(1)                                                     R MOD.
      INTEGER  IW(1), N, LEN, IDX, IDB                                  I***
C
      RETURN
      END
CSOLVE
      SUBROUTINE SOLVE(X, B, K, D, N)
C
C **********************************************************************
C
C     PURPOSE - (VER = 1)
C
C         SOLVES (LDL(T))*X=B WITH THE USUAL METHOD, VIZ.
C
C         1) SOLVE L*Z=B      (FORWARD SUBSTITUTION)
C         2) COMPUTE U=(D**(-1))*Z
C         3) SOLVE L(T)*X=U   (BACK SUBSTITUTION)
C
C         IN THIS ROUTINE Z,U, AND X OCCUPY THE SAME MEMORY.
C
C
C     INPUT PARAMETERS -
C
C         B = KNOWN VECTOR.
C         D = POINTER VECTOR TO THE DIAGONAL ELEMENTS IN L.
C         K = CONTAINS THE STRICTLY LOWER TRIANGLE OF L, AND
C             THE DIAGONAL, D, OF LDL(T) (STORED IN THE DIAGONAL OF K).
C         N = DIMENSION OF L MATRIX.
C
C
C     OUTPUT PARAMETERS -
C
C         X = SOLUTION.
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          B(1), K(1), SCPR, X(1), XNDECR                          R MOD.
      INTEGER   D(1), DI, DIM1, I, LENM1, N, NDECR, ROWNO, TOP          I***
C
C     *********************
C     FORWARD SUBSTITUTION.
C     *********************
      DIM1 = 0
      DO 20 I = 1, N
        DI = D(I)
        TOP = DIM1 + 1
        LENM1 = DI - TOP
        DIM1 = DI
        ROWNO = I - LENM1
        IF(.NOT. LENM1 .GT. 0) GOTO 10
          X(I) = B(I) - SCPR(K(TOP), X(ROWNO), LENM1)
          GOTO 20
10      X(I) = B(I)
20      CONTINUE
C
C     **********
C     D-INVERSE.
C     **********
      DO 40 I = 1, N
        DI = D(I)
        X(I) = X(I) / K(DI)
40      CONTINUE
C
      IF(N .EQ. 1) GOTO 9999
C
C     ******************
C     BACK SUBSTITUTION.
C     ******************
      NDECR = N
      DO 50 I = 2, N
        TOP = D(NDECR - 1) + 1
        LENM1 = D(NDECR) - TOP
        ROWNO = NDECR - LENM1
        XNDECR = X(NDECR)
        IF(LENM1 .GT. 0) CALL SUBV(X(ROWNO), K(TOP), XNDECR, LENM1)
        NDECR = NDECR - 1
50      CONTINUE
C
9999  RETURN
      END
CSORTP
      SUBROUTINE SORTP(POINTR, LAMBDA, LAST)
C
C **********************************************************************
C
C     PURPOSE -
C
C         THIS IS A USUAL SELECTION SORT WITH ONE DIFFERENCE. IT
C         SORTS THE POINTR VECTOR WITH RESPECT TO LAMBDA (I.E. WE
C         CHANGE POINTR AND NOT LAMBDA).
C
C
C     INPUT PARAMETERS -
C
C         LAMBDA       = LAMBDA EIGENVALUES.
C         POINTR(LAST) = NUMBER OF ELEMENTS TO BE SORTED.
C
C
C     OUTPUT PARAMETERS -
C
C         POINTR IS SORTED SO THAT  (IF P =POINTR(LAST))
C         LAMBDA(IABS(POINTR(1))), ..., LAMBDA(IABS(POINTR(P))) ARE IN
C         ASCENDING ORDER.
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          LAMBDA(1), TLAM                                         R MOD.
      INTEGER   I, IP1, J, K, LAST, P, PM1, POINT, POINTR(1)            I***
C
      P = POINTR(LAST)
      PM1 = P - 1
      IF(PM1 .LE. 0) GOTO 9999
C
      DO 20 I = 1, PM1
        K = I
        POINT = IABS(POINTR(K))
        TLAM = LAMBDA(POINT)
        IP1 = I + 1
C
        DO 10 J = IP1, P
          POINT = IABS(POINTR(J))
C         ********
C         COMPARE.
C         ********
          IF(.NOT. LAMBDA(POINT) .LT. TLAM) GOTO 10
            K = J
            POINT = IABS(POINTR(K))
            TLAM = LAMBDA(POINT)
10        CONTINUE
C
C       ***************************
C       INTERCHANGE POINTER VALUES.
C       ***************************
        IF(.NOT. I .NE. K) GOTO 20
          POINT = POINTR(I)
          POINTR(I) = POINTR(K)
          POINTR(K) = POINT
C
20      CONTINUE
C
9999  RETURN
      END
CSTLM
      SUBROUTINE STLM(N, A, B, MAXL, PROFIL, PMAX, MXREST, MSGLVL,
     +                MAXW, MAXIW, DAFILE, MAXREC, KFILE, X, BG,
     +                TCONV, NLEFT, ERRNO, W, IW)
C
C **********************************************************************
C
C     PURPOSE -
C
C         MAIN ROUTINE FOR STLM.
C
C         PLEASE SEE THE USER GUIDE FOR FURTHER COMMENTS.
C
C
C **********************************************************************
C
C
      DOUBLE PRECISION                                                  R INS.
     +          A, A1, ALTMU, B, B1, CK, CL, COEFF, FACTOR, SRELPR, MU, R MOD.
     +          NEXTMU, OLDMU, RDUMP, RV, SECOND, SMALL, TIME, TIMTQL,
     +          TLDL, TOLBPI, TOLLDL, TOLPDM, TOLS1I, TOLZBT, TOLZNU,
     +          TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1), WRR, BG
      INTEGER   ACTIVE, ADDRSS, ALPHA, BETA, BETAPI, CNEG, CNEGF, CONV, I***
     +          COPT, COUNT, CPOS, D, DAFIL1, DAFILE, DUMMY, ERRNO,
     +          FREE, IDUMP, ITERNO, IV, IW(1), K, KFILE, KFILE1,
     +          LAMBDA, LEFTP, LEN, LENADR, LP, M, MAXIW, MAXIW1, MAXL,
     +          MAXL1, MAXRE1, MAXREC, MAXW, MAXW1, MSGLVL, MV, MXNEW,
     +          MXOLD, MXREST, MXRST, N, N1, N2, NBADMU, NIL, NMXRST,
     +          NOACTN, NOR, NREAD, NU, NUMEIG, NUMVEC, NWRITE,
     +          OLCPOS, P, PFCONV, PMAX, POINTR, POPT, READID, READK,
     +          REST, RFIRST, RIGHTC, RIGHTM, RIGHTP, RNEW, ROLD, S,
     +          SAEVAL, SAVE, SAVFRE, SCPX, SCR, SOLCPX, STADEW, TCONV,
     +          NERR  , NOUT , V, VER, WAD, WRI, WRITID, X, X1,
     +          TCONV1, NLEFT, ERRNO1, PMAX1, MXRES1, MSGLV1, PROFIL,
     +          PROFI1
      LOGICAL   ALWAYS, DECOMP, DIAGM, F, FRSTIT, MEQI,                 L***
     +          SAFRST, SELDOM, T, TERMIN, UPDATE, USEMX, WRL, ZERBET,
     +          REACHB, USSMXR, USEDB
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMAD/ ADDRSS, DAFIL1, KFILE1, LP, MAXL1, NREAD, NWRITE
      COMMON   /STLMCT/ N1, ITERNO, TCONV1, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMDS/ K, M, D, DIAGM
      COMMON   /STLMER/ RDUMP, ERRNO1, IDUMP(2)
      COMMON   /STLMEW/ LEFTP, LENADR, MAXRE1, NUMVEC, RIGHTC, RIGHTM,
     +                  RIGHTP, STADEW
      COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX,
     +                  X1
      COMMON   /STLMIN/ A1, B1, NUMEIG, MAXW1, MAXIW1
      COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N2
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMPL/ PMAX1, POPT, COPT, MXRES1
      COMMON   /STLMPF/ PROFI1
      COMMON   /STLMPR/ MSGLV1, NERR  , NOUT
      COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,
     +                  TOLZNU, TOLLDL(3)
      COMMON   /STLMTS/ REACHB, USSMXR, USEDB
      COMMON   /STLMUI/ RV(5), IV(9)
      COMMON   /STLMVR/ SMALL, VER
      COMMON   /STLMWH/ WAD(2), ACTIVE(2)
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
C
      DATA      NOR  / 20 /
C
C     ********************
C     INITIALIZATION PART.
C     ********************
      CALL INITD(N, A, B, MAXL, PROFIL, PMAX, MXREST, MSGLVL,
     +           MAXW, MAXIW, DAFILE, MAXREC, KFILE, X, BG,
     +           TCONV, NLEFT, ERRNO, W, IW, LEN)
      IF(LEN .GT. 0) RETURN
C
C     ***********************
C     WRITE INPUT PARAMETERS.
C     ***********************
      CALL WRINFO(NOR, 1, W, IW)
C
      TIME(NOR) = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
C
C     **********
C     MAIN LOOP.
C     **********
10      ITERNO = ITERNO + 1
        TSAVE = SECOND()
        CALL WRINFO(NOR, 2, W, IW)
C       *******************************************
C       COMPUTE LDL(T)-DECOMPOSITION OF K - MU * M.
C       *******************************************
        IF(.NOT. DECOMP(W, IW)) GOTO 8888
C
        CALL WRINFO(NOR, 3, W, IW)
        IF(USEDB .AND. REST .EQ. 0) GOTO 6666
C
        IF(.NOT. ITERNO .EQ. 1) GOTO 20
C         ****************************************
C         INITIALIZATIONS FOR THE FIRST ITERATION.
C         ****************************************
          RV(1) = MU
          IF(.NOT. FRSTIT(W, IW)) GOTO 8888
C
C       *************************************
C       HANDLES THE COMPUTATIONS FOR A SHIFT.
C       *************************************
20      IF(VER .LT. 3) CALL INITEW(IW(ADDRSS))
        IF(.NOT. ALWAYS(W, IW)) GOTO 8888
C
        IF(.NOT. REST .GT. 0) GOTO 30
C         ********************
C         MISSING EIGENVALUES.
C         ********************
          IF(.NOT. SELDOM(W, IW)) GOTO 8888
C
30      CALL WRINFO(NOR, 4, W, IW)
C       *********************************************
C       CHECK IF WE HAVE REACHED THE LAST EIGENVALUE.
C       *********************************************
        IF(CPOS + RNEW .EQ. N) GOTO 7777
C
C       ************************************
C       DECIDE IF WE SHOULD CONTINUE OR NOT.
C       ************************************
6666    IF(TERMIN(1, W, IW, BG, TCONV, NLEFT, ERRNO)) GOTO 9999
C       *********************************
C       UPDATE CERTAIN CONTROL VARIABLES.
C       *********************************
        CALL REPL
C
        GOTO 10
C
C     ************************************
C     WE HAVE REACHED THE LAST EIGENVALUE.
C     ************************************
7777  IF(TERMIN(2, W, IW, BG, TCONV, NLEFT, ERRNO)) GOTO 9999
      GOTO 9999
C
C     ***********************************
C     ERROR EXIT, UNLESS ERRORS IN INPUT.
C     ***********************************
8888  CALL ERROR(NOR, NOR)
      IF(TERMIN(3, W, IW, BG, TCONV, NLEFT, ERRNO)) GOTO 9999
C
9999  TIME(NOR) = SECOND() - TIME(NOR)
C
C
C     *************************
C     WRITE OUTPUT INFORMATION.
C     *************************
      CALL WRINFO(NOR, 5, W, IW)
C
      X = X1
C
      RETURN
      END
CSUB3
      SUBROUTINE SUB3(ID1, ID2, C, N, W, IW, LEN)
C
C **********************************************************************
C
C     PURPOSE (VER = 3)
C
C     DUMMY ROUTINE. PLEASE SEE THE INSTALLATION GUIDE FOR
C     MORE DETAILS.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +         C, W(1)                                                  R MOD.
      INTEGER  ID1, ID2, N, IW(1), LEN                                  I***
C
      RETURN
      END
CSUBV
      SUBROUTINE SUBV(X, Y, C, N)
C
C **********************************************************************
C
C     PURPOSE - (VER = 1 OR 2)
C
C         COMPUTES X=X-C*Y.
C
C
C     INPUT PARAMETERS -
C
C         X AND Y = N-VECTORS PARTICIPATING IN THE COMPUTATION.
C         C       = REAL CONSTANT.
C         N       = DIMENSION OF THE VECTORS (NOT NECESSARILY OF THE
C                   PROBLEM).
C
C
C     OUTPUT PARAMETERS -
C
C         X = UPDATED X.
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          C, X(1), Y(1)                                           R MOD.
      INTEGER   I, N                                                    I***
C
      DO 10 I = 1, N
        X(I) = X(I) - C * Y(I)
10      CONTINUE
C
      RETURN
      END
CSUBVEC
      LOGICAL FUNCTION SUBVEC(ID1, ID2, C, W, IW, ACTION, CASE)         L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         DRIVER ROUTINE FOR SUBV. COMPUTES
C         VEC(ID1)=VEC(ID1)-C*VEC(ID2).
C
C
C     INPUT PARAMETERS -
C
C         ID1    = IDENTIFIER OF THE RESULTING VECTOR.
C         ID2    = IDENTIFIER OF THE INCOMING VECTOR.
C         C      = REAL CONSTANT.
C         ACTION = DETERMINES IF VEC(ID1) OR VEC(ID2) SHOULD
C                  BE SAVED OR NOT.
C         CASE   = NOT USED IN THIS VERSION.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          C, SECOND, ST, TIME, W(1), SMALL                        R MOD.
      INTEGER   ACTION, AD1, AD2, ADDRSS, CASE, COUNT, DAFILE, DUMMY,   I***
     +          FREE, ID1, ID2, IW(1), KFILE, LP, MAXL, N, NBADMU,
     +          NMXRST, NOACTN, NOR, NREAD, NWRITE, READID, READK,
     +          SAEVAL, SAVE, SAVFRE, WRITID, VER, LEN
      LOGICAL   ALLOC, F, IO, T                                         L***
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
      COMMON   /STLMVR/ SMALL, VER
C
      DATA      NOR  / 21 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      SUBVEC = T
C
      IF(.NOT. VER .LE. 2) GOTO 10
C       *************
C       VER = 1 OR 2.
C       *************
        IF(.NOT. ALLOC(ID1, ID2, AD1, AD2, ACTION, W, IW(ADDRSS) ))
     +                 GOTO 8888
C
        CALL SUBV(W(AD1), W(AD2), C, N)
C
        IF(.NOT. (ACTION .EQ. SAVE .OR. ACTION .EQ. SAVFRE)) GOTO 9999
          IF(.NOT. IO(W(AD1), ID1, WRITID, N, IW(ADDRSS) )) GOTO 8888
          GOTO 9999
C
C     ********
C     VER = 3.
C     ********
10    CALL SUB3(ID1, ID2, C, N, W, IW, LEN)
C
      IF(.NOT. LEN .NE. 0) GOTO 9999
        CALL ERROR(106, LEN)
        CALL ERROR(106, 106)
C
8888  SUBVEC = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CTERMIN
      LOGICAL FUNCTION TERMIN(CAUSE, W, IW, BG, TCONV, NLEFT, ERRNO)    L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         THIS ROUTINE CHECKS IF WE HAVE COMPUTED THE
C         REQUESTED EIGENPAIRS OR NOT.
C
C     INPUT PARAMETERS -
C
C         CAUSE = 1, DECIDE IF WE SHOULD CONTINUE OR NOT.
C                 2, WE HAVE REACHED THE LAST EIGENVALUE.
C                 3, AN ERROR HAS OCCURRED.
C
C     OUTPUT PARAMETERS -
C
C         TERMIN  = FALSE, IF WE SHOULD TAKE ANOTHER SHIFT,
C                   AND TRUE OTHERWISE.
C
C     PLEASE SEE THE PROGRAMMERS AND USER GUIDES FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          A, B, CK, CL, COEFF, MU, SMALL, TIMTQL, TLDL, TOPINV,   R MOD.
     +          TOPM, TOTALT, TPRED, TPRED1, TSAVE, TVECOP, W(1), BG,
     +          TW, WRR
      INTEGER   ADDRSS, CAUSE, CNEGF, CONV, COPT, CPOS, DAFILE, ERRNO,  I***
     +          I, ITERNO, IW(1), KFILE, LP, MAXIW, MAXL, MAXW, MXREST,
     +          N, NREAD, NUMEIG, NWRITE, PFCONV, PMAX, POPT, RFIRST,
     +          RNEW, TCONV, VER, X, TM1, TCSAV, IP1, J, K, TIW, ISTART,
     +          WRI, TCONV1, NLEFT, ERRNO1
      LOGICAL   F, T, UPDATE, REACHB, USSMXR, USEDB, WRL                L***
      COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMIN/ A, B, NUMEIG, MAXW, MAXIW
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED1,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
      COMMON   /STLMTF/ T, F
      COMMON   /STLMTS/ REACHB, USSMXR, USEDB
      COMMON   /STLMUI/ MU(3), TOTALT, TPRED, CNEGF, CPOS, ERRNO1,
     +                  ITERNO, N, RFIRST, RNEW, TCONV1, X
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
      COMMON   /STLMVR/ SMALL, VER
C
C
C     ***********
C     INITIALIZE.
C     ***********
      CALL UI
      NLEFT = RFIRST
      ERRNO = ERRNO1
      TCONV = TCONV1
C
C     *********************************************************
C     TERMIN = T, IF WE HAVE REACHED B, THE LAST EIGENVALUE, OR
C                 IF AN ERROR HAS OCCURRED.
C     **********************************************************
      TERMIN = (REACHB .AND. USEDB) .OR. CAUSE .GT. 1
      IF(REACHB) USEDB = T
      IF(.NOT. TERMIN) GOTO 9999
C     ********************************************************
C     BG = LAST SHIFT.
C     IF (REACHED THE LAST EIGENVALUE) BG = B
C     IF (ERROR) BG = THE PREVIOUS SHIFT
C     IF (FIRST ITERATION AND ERROR) BG = A
C     IF (NOT FIRST ITERATION AND CAN NOT COMPUTE A NEW SHIFT)
C           BG = LAST SHIFT
C     ********************************************************
      BG = MU(3)
      IF(CAUSE .EQ. 2) BG = B
      IF(CAUSE .EQ. 3) BG = MU(2)
      IF(ITERNO .EQ. 1  .AND. CAUSE .EQ. 3) BG = A
      IF(ITERNO .GT. 1  .AND.  ERRNO .EQ. 203) BG = MU(3)
      IF(BG .GT. B) BG = B
      IF(BG .LT. A) BG = A
C
C
      IF(TCONV .LE. 0) GOTO 9999
C
C     *********************************************************
C     CONSTRUCT POINTER VECTOR IW, AND MOVE EIGENVALUES TO THE
C     BEGINNING OF W. RESULTS IN W(IW(1)), ..., W(IW(TCONV)) IN
C     ASCENDING ORDER.
C     *********************************************************
      DO 5 I = 1, TCONV
        IW(I) = I
        W(I)  = W(LP)
        LP    = LP + 1
5       CONTINUE
C
      TM1 = TCONV - 1
      IF(TM1 .LE. 0) GOTO 30
      DO 20 I = 1, TM1
        K   = I
        IP1 = I + 1
C
        DO 10 J = IP1, TCONV
          IF(W(J) .LT. W(K)) K = J
10        CONTINUE
C
        IF(.NOT. I .NE. K) GOTO 20
          TIW   = IW(I)
          IW(I) = IW(K)
          IW(K) = TIW
C
          TW    = W(I)
          W(I)  = W(K)
          W(K)  = TW
20        CONTINUE
C
30    IF(.NOT. ERRNO .NE. 1701) GOTO 55
C       ***********************************************************
C       IF NOT ERROR DUE TO MAXL TOO SHORT, FIND LARGEST EIGENVALUE
C       SMALLER THAN OR EQUAL TO BG.
C       ***********************************************************
        I = TCONV + 1
        DO 40 J = 1, TCONV
          I = I - 1
          IF(W(I) .LE. BG) GOTO 50
40        CONTINUE
        I = 0
C
C     **************************************
C     THROW ALL EIGENVALUES GREATER THAN BG.
C     **************************************
50    TCONV = I
C
C
55    IF(.NOT. (MU(1) .LT. A  .AND. TCONV .GT. 0) ) GOTO 9999
        TCSAV = TCONV
C       ****************************************************
C       FIND SMALLEST EIGENVALUE GREATER THAN OR EQUAL TO A.
C       ****************************************************
        DO 60 I = 1, TCONV
          IF(W(I) .GE. A) GOTO 70
60        CONTINUE
C
        I = TCONV + 1
C
C       *********************************
C       THROW THOSE THAT ARE LESS THAN A.
C       *********************************
70      TCONV = TCONV - (I - 1)
C
C
        ISTART = I
C       ************************************************
C       ADJUST NLEFT, AND MOVE EIGENVALUES AND POINTERS.
C       ************************************************
        NLEFT  = RFIRST + ISTART - 1
        IF(.NOT. (TCONV .GT. 0  .AND.  TCSAV .NE. TCONV) ) GOTO 9999
          J = 1
          DO 80 I = ISTART, TCSAV
            W(J)  = W(I)
            IW(J) = IW(I)
            J     = J + 1
80          CONTINUE
C
9999  WRI(1) = TCONV
      WRI(2) = NLEFT
      WRR(1) = BG
C
C
      RETURN
      END
CTRANSF
      LOGICAL FUNCTION TRANSF(ADD, SCOL, XORMX, VORMV, W, IW)           L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         FOR A GIVEN S-VECTOR THIS ROUTINE COMPUTES X=V*S-ADD*V(P+1)
C         OR MX=MV*S-ADD*MV(P+1), WHERE ADD=-S(P)*BETA(P)/NU,
C         (ADD=-S(P)/NU IF ZERBET IS TRUE). V AND MV ARE MATRICES, SO
C         THE ROUTINE PERFORMS ESSENTIALLY A MATRIX VECTOR
C         MULTIPLICATION.
C
C
C     INPUT PARAMETERS -
C
C         ADD    = -S(P,I)*BETA(P)/NU(I) FOR AN I, IF ZERBET IS FALSE.
C                = -S(P, I)/NU(I), IF ZERBET = TRUE.
C         SCOL   = (S(1,I), ..., S(P,I))(T), I.E. COLUMN NUMBER I IN S.
C         XORMX  = THE COMPUTED X (OR MX) VECTOR GETS THE IDENTIFIER
C                  XORMX+1.
C         VORMV  = THE COLUMNS IN THE V (OR MV) MATRIX HAVE THE
C                  IDENTIFIERS, VORMV+1, ..., VORM+P, (VORMV+P+1).
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          ADD, SCOL(1), SECOND, ST, TIME, W(1)                    R MOD.
      INTEGER   CNEG, COUNT, CPOS, DUMMY, FREE, I, ITERNO,              I***
     +          IW(1), N, NBADMU, NMXRST, NOACTN, NOR, OLCPOS, P, REST,
     +          RNEW, ROLD, SAVE, SAVFRE, TCONV, VORMV, XORMX
      LOGICAL   F, MULVEC, SUBVEC, T, USEMX, ZERBET                     L***
      COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTF/ T, F
C
      DATA      NOR  / 22 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      TRANSF = T
C     ************************
C     X = V1 * S1 (OR MX, MV).
C     ************************
      IF(.NOT. MULVEC(XORMX, VORMV + 1, SCOL(1), W, IW, FREE, 1))
     +      GOTO 8888
C
      IF(.NOT. P .GT. 1) GOTO 50
        DO 40 I = 2, P
C         ********************************************
C         X = X + V2 * S2 + ... + VP * SP (OR MX, MV).
C         ********************************************
          IF(.NOT. SUBVEC(XORMX, VORMV + I, - SCOL(I), W, IW, FREE, 1))
     +                                      GOTO 8888
40        CONTINUE
C
C     *********************************
C     X = X - ADD * V(P+1) (OR MX, MV).
C     *********************************
50    IF(.NOT. SUBVEC(XORMX, VORMV + P + 1, ADD, W, IW, SAVFRE, 1))
     +                                     GOTO 8888
      CALL FREEID(VORMV + P + 1)
      GOTO 9999
C
8888  TRANSF = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CTRIDIG
      LOGICAL FUNCTION TRIDIG(CONTIN, FINAL, POSNU, UPDMU, W, IW)       L***
C
C **********************************************************************
C
C     PURPOSE -
C
C         THIS ROUTINE HANDLES ALL COMPUTATIONS CONNECTED WITH THE
C         EIGENSYSTEM OF T (THE TRIDIAGONAL P X P-MATRIX).
C
C
C     INPUT PARAMETERS -
C
C         FINAL = TRUE IF WE HAVE TAKEN THE LAST LANCZOS STEP, FALSE
C                 OTHERWISE.
C         UPDMU = TRUE, IF NEXTMU SHOULD BE COMPUTED, FALSE
C                 OTHERWISE.
C
C
C     OUTPUT PARAMETERS -
C
C         CONTIN = TRUE IF WE CAN TAKE AN OTHER STEP, I.E. IF THERE IS
C                  ROOM LEFT (P .LT. PMAX), VECTORS HAVE NOT BEGUN TO
C                  SPLIT (NUMDEL .EQ. 0), AND BETA(P) IS LARGE ENOUGH
C                  (.NOT. ZERBET). FALSE OTHERWISE.
C         POSNU  = TRUE IF WE HAVE A POSITIVE NU EIGENVALUE OR IF MU IS
C                  GREATER THAN LAMBDA(N), OR IF WE HAVE USED
C                  B AS SHIFT. FALSE OTHERWISE.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          CK, CL, COEFF, MAXNU, MINNU, SECOND, ST, TIME, TIMTQL,  R MOD.
     +          TLDL, TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1), WRR
      INTEGER   ALPHA, BETA, BETAPI, CNEG, CNEGF, CONV, COPT, COUNT,    I***
     +          CPOS, DUMMY, ITERNO, IW(1), LAMBDA, MXREST, N, NBADMU,
     +          NMXRST, NOR, NU, NUMDEL, OLCPOS, P, PFCONV, PMAX,
     +          POINTR, POPT, REST, RFIRST, RNEW, ROLD, S, SCR, TCONV,
     +          WRI
      LOGICAL   CONTIN, F, FINAL, IMTQL2, POSNU, SAFRST, T, UPDATE,     L***
     +          UPDMU, USEMX, WRL, ZERBET, REACHB, USSMXR, USEDB,
     +          LGEB
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
      COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTS/ REACHB, USSMXR, USEDB
      COMMON   /STLMTF/ T, F
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
C
      DATA      NOR  / 23 /
C
      ST = SECOND()
      COUNT(NOR) = COUNT(NOR) + 1
      TRIDIG = T
C
C     *******************************
C     SOLVE TRIDIAGONAL EIGENPROBLEM.
C     *******************************
      IF(.NOT. IMTQL2(W(ALPHA), W(BETA), W(NU), W(SCR), W(S), P, PMAX,
     +                FINAL)) GOTO 8888
C
C     **************************
C     COMPUTE LAMBDA AND BETAPI.
C     **************************
      CALL COMPL(W(BETA), W(BETAPI), W(NU), W(LAMBDA), W(S),
     +           MAXNU, MINNU, PMAX, LGEB)
C
C     ******************
C     CHECK CONVERGENCE.
C     ******************
      CALL CONVER(W(BETAPI), W(LAMBDA), IW(POINTR), W(S), PMAX)
C
C     ************
C     SORT POINTR.
C     ************
      CALL SORTP(IW(POINTR), W(LAMBDA), PMAX + 1)
C
      IF(FINAL) CALL WRINFO(NOR, 1, W, IW)
C
C     ***********************
C     DELETE DUPLICATE PAIRS.
C     ***********************
      CALL DELDUP(W(NU), IW(POINTR), W(S), MAXNU, MINNU, NUMDEL, PMAX,
     +            FINAL)
C
      CONV = CONV - NUMDEL
      IF(CONV .GT. 0 .AND. PFCONV .EQ. 0) PFCONV = P
C     ******************************************************
C     NOTE RNEW .EQ. N. WE CAN NOT EXPECT TO GET A NEW SHIFT
C     IF WE ARE TO THE RIGHT OF LAMBDA(N).
C     ******************************************************
      POSNU = MAXNU .GT. 0.0D0 .OR. RNEW .EQ. N .OR. USEDB              D MOD.
C
C     **********************************************
C     CHECK IF WE SHOULD TAKE AN OTHER LANCZOS STEP.
C     **********************************************
      CONTIN = P .LT. N .AND. P .LT. PMAX .AND. NUMDEL .EQ. 0 .AND.
     +         .NOT. ZERBET
      IF(.NOT. FINAL) GOTO 9999
C       ***************************************************
C       NO MORE LANCZOS STEPS. RNEW .EQ. N, SEE NOTE ABOVE.
C       COMPUTE THE NEW SHIFT.
C       ***************************************************
        IF(RNEW .EQ. N) CPOS = 0
        IF(UPDMU .AND. MAXNU .GT. 0.0D0 .AND. RNEW .LT. N               D MOD.
     +     .AND. (.NOT. USEDB) )
     +     CALL NEWMU(W(LAMBDA), IW(POINTR), MAXNU, LGEB)
C
        IF(ITERNO .EQ. 1) CNEGF = CNEG
C
        WRI(1) = NUMDEL
        WRL(1) = POSNU
        WRR(1) = MINNU
        WRR(2) = MAXNU
        CALL WRINFO(NOR, 2, W, IW)
C
        GOTO 9999
C
8888  TRIDIG = F
      CALL ERROR(NOR, NOR)
C
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
      RETURN
      END
CUI
      SUBROUTINE UI
C
C **********************************************************************
C
C     PURPOSE -
C
C         WHEN CALLED THIS ROUTINE FILLS THE COMMONBLOCK /STLMUI/ WITH
C         INFORMATION.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          ALTMU, CK, CL, COEFF, MU, NEXTMU, OLDMU, R, RDUMP,      R MOD.
     +          SECOND, TIME, TIMTQL, TLDL, TOPINV, TOPM, TPRED, TSAVE,
     +          TVECOP
      INTEGER   CNEG, CNEGF, CONV, COUNT, CPOS, DUMMY, ERRNO, I, IDUMP, I***
     +          ITERNO, MV, MXNEW, MXOLD, MXRST, N, NBADMU, NIL, NMXRST,
     +          OLCPOS, P, PFCONV, REST, RFIRST, RNEW, ROLD, SCPX,
     +          SOLCPX, TCONV, V, X
      LOGICAL   SAFRST, UPDATE, USEMX, ZERBET                           L***
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMUI/ R(5), I(9)
C
      R(2) = OLDMU
      R(3) = MU
      R(4) = SECOND() - TIME(20)
      R(5) = TPRED
C
      I(1) = CNEGF
      I(2) = CPOS
      I(3) = ERRNO
      I(4) = ITERNO
      I(5) = N
      I(6) = RFIRST
      I(7) = RNEW
      I(8) = TCONV
      I(9) = X
C
      RETURN
      END
CWRINF3
      SUBROUTINE WRINF3(NOR,LOC,LAMBDA,POINTR)
C
C **********************************************************************
C
C     PURPOSE -
C
C         WRITES INFORMATION FOR MSGLVL=2 AND 3.
C
C     INPUT PARAMETERS -
C
C         NOR    = THE NUMBER OF THE ROUTINE CALLING WRINFO.
C         LOC    = THE LOCATION IN THAT ROUTINE.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          A, ALTMU, B, CK, CL, COEFF, FACTOR, LAMBDA(1), SRELPR,  R MOD.
     +          MU, NEXTMU, OLDMU, RDUMP, SECOND, SMALL, TEMP, TIME,
     +          TIMTQL, TLDL, TOLBPI, TOLLDL, TOLPDM, TOLS1I, TOLZBT,
     +          TOLZNU, TOPINV, TOPM, TPRED, TSAVE, TVECOP, WRR
      INTEGER   ADDRSS, CNEG, CNEGF, CONV, COPT, COUNT, CPOS, D, DAFILE,I***
     +          DUMMY, ERRNO, I, IDUMP, IMAX, ITERNO, K, KFILE, LEFTP,
     +          LENADR, LOC, LP, M, MAXIW, MAXL, MAXREC, MAXW, MSGLVL,
     +          MXREST, N, NBADMU, NMXRST, NOR, NREAD, NUMEIG, NUMVEC,
     +          NWRITE, OLCPOS, P, PFCONV, PMAX, POINTN, POINTP,
     +          POINTR(1), POPT, REST, RFIRST, RIGHTC, RIGHTM, RIGHTP,
     +          RNEW, ROLD, SACNEG, SACPOS, STADEW, TCONV, NERR  ,
     +          NOUT , VER, PROFIL, WRI
      LOGICAL   DIAGM, MEQI, SAFRST, UPDATE, USEMX, ZERBET, WRL         L***
      COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMDS/ K, M, D, DIAGM
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMEW/ LEFTP, LENADR, MAXREC, NUMVEC, RIGHTC, RIGHTM,
     +                  RIGHTP, STADEW
      COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
      COMMON   /STLMIN/ A, B, NUMEIG, MAXW, MAXIW
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
      COMMON   /STLMPF/ PROFIL
      COMMON   /STLMPR/ MSGLVL, NERR  , NOUT
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,
     +                  TOLZNU, TOLLDL(3)
      COMMON   /STLMVR/ SMALL, VER
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
C
C
      IF(.NOT. NOR .EQ. 20) GOTO 100
      GOTO(10,15,20,30,60),LOC
C
C     *******************************
C     ROUTINE STLM. INPUT PARAMETERS.
C     *******************************
10    WRITE(NOUT ,1000) VER,N,A,B,MAXL,PROFIL,PMAX,MXREST,
     +                  MSGLVL, MAXW, MAXIW, NERR  , NOUT
C
1000  FORMAT(1H1///1X,60(1H*),2(/1X,1H*,58X,1H*)/1X,1H*,10X,
     +       38H I N P U T        T O        S T L M  ,
     +       10X,1H*,2(/1X,1H*,58X,1H*)/1X,60(1H*)///
     +       4X, 3HVER, I54/
     +       4X, 1HN, I56/
     +       4X, 1HA, 1P1D56.7/                                         F MOD.
     +       4X, 1HB, 1P1D56.7/                                         F MOD.
     +       4X, 4HMAXL, I53/
     +       4X, 6HPROFIL, I51/
     +       4X, 4HPMAX, I53/
     +       4X, 6HMXREST, I51/
     +       4X, 6HMSGLVL, I51/
     +       4X, 4HMAXW, I53/
     +       4X, 5HMAXIW, I52/
     +       4X, 6HNERR  , I51/
     +       4X, 5HNOUT , I52)
C
      IF(VER .LE. 2) WRITE(NOUT , 1002) DAFILE, MAXREC
1002  FORMAT(4X, 6HDAFILE, I51/ 4X, 6HMAXREC, I51)
C
      IF(VER .EQ. 1) WRITE(NOUT , 1001) KFILE
1001  FORMAT(4X, 5HKFILE, I52)
C
      WRITE(NOUT , 1003)
1003  FORMAT(//1X,13H*** INPUT ***)
C
      GOTO 9999
C
C     *******************
C     START OF ITERATION.
C     *******************
15    IF(ITERNO .EQ. 1) WRITE(NOUT ,1010)
1010  FORMAT(1H1)
C
      WRITE(NOUT ,1020) ITERNO
1020  FORMAT(1X,60(1H=)///4X,16HITERATION NUMBER,I41)
C
      GOTO 9999
C
C     *************
C     AFTER DECOMP.
C     *************
20    WRITE(NOUT ,1025) MU,RNEW,REST
1025  FORMAT(/4X,5HSHIFT,1P1D52.7/                                      F MOD.
     +       4X,38HNO. OF EIGENVALUES LESS THAN THE SHIFT,I19/
     +       4X,44HNO. OF REMAINING EIGENVALUES IN THE INTERVAL,I13)
C
      GOTO 9999
C
C     ****************************
C     BEFORE END OF THE ITERATION.
C     ****************************
30    TEMP=SECOND()-TSAVE
      WRITE(NOUT ,1030) TEMP
1030  FORMAT(4X,19HTHIS ITERATION TOOK,1P1D30.1,8H SECONDS)             F MOD.
C
      IF(.NOT. CPOS+RNEW .EQ. N) GOTO 40
        WRITE(NOUT ,1050)
1050    FORMAT(4X,36HWE HAVE REACHED THE LAST EIGENVALUE.)
C
        GOTO 50
C
40    IF(UPDATE) WRITE(NOUT ,1060) TPRED
1060  FORMAT(4X,37HPREDICTED TIME FOR THE NEXT ITERATION,
     +       1P1D12.1,8H SECONDS)                                       F MOD.
C
50    WRITE(NOUT ,1080)
1080  FORMAT(/1X,13H*** SHIFT ***)
C
      GOTO 9999
C
C     ******************
C     OUTPUT PARAMETERS.
C     ******************
60    WRITE(NOUT ,1090) ITERNO,WRI(1),WRR(1)
1090  FORMAT(1H1///1X,60(1H*),2(/1X,1H*,58X,1H*)/1X,1H*,7X,
     +       44H O U T P U T        F R O M        S T L M  ,
     +       7X,1H*,2(/1X,1H*,58X,1H*)/1X,60(1H*)///
     +       4X,20HNUMBER OF ITERATIONS,I37/
     +       4X,35HTOTAL NUMBER OF COMPUTED EIGENPAIRS,I22/
     +       4X,10HBG        ,1P1D47.7)                                 F MOD.
C
      IF(ERRNO .LT. 2501 .OR. ERRNO .GT. 2513) WRITE(NOUT ,1100) WRI(2)
1100  FORMAT(/4X,33HNUMBER OF EIGENVALUES LESS THAN A,I24)
C
      IF(SAFRST) WRITE(NOUT ,1110) CNEGF
1110  FORMAT(4X,8HOF THESE,I5,15H WERE ACCEPTED.)
C
      IF(ERRNO .NE. 0) WRITE(NOUT ,1120) ERRNO
1120  FORMAT(/4X,31HAN ERROR OCCURRED, ERROR NUMBER,I26)
C
      WRITE(NOUT ,1130) TIME(20)
1130  FORMAT(/4X,27HTOTAL TIME USED FOR THE RUN,1P1D22.1,8H SECONDS)    F MOD.
C
      IF(.NOT. TCONV .GT. 0) GOTO 70
        TEMP=TCONV
        TEMP=TIME(20)/TEMP
        WRITE(NOUT ,1140) TEMP
1140    FORMAT(4X,14HTIME/EIGENPAIR,1P1D35.1,8H SECONDS)                F MOD.
C
70    WRITE(NOUT ,1150) TLDL
1150  FORMAT(4X,33HTHE FIRST LDLT-DECOMPOSITION TOOK,1P1D16.1,          F MOD.
     +       8H SECONDS)
C
      TEMP=N
      TEMP=TVECOP/TEMP
      IF(UPDATE) WRITE(NOUT ,1160) TOPINV,TOPM,TEMP
1160  FORMAT(4X,24HOPINV  (SOLVES LDLT*X=B),1P1D25.1,8H SECONDS/        F MOD.
     +       4X,14HOPM    (Y=M*X),1P1D35.1,8H SECONDS/                  F MOD.
     +       4X,16HVECTOROPERATIONS,1P1D18.1,23H SECONDS/MULTIPLICATION)F MOD.
C
      WRITE(NOUT ,1170)
1170  FORMAT(//1X,14H*** OUTPUT ***//1X,12H*** STLM ***)
C
      GOTO 9999
C
C     ***************
C     ROUTINE ALWAYS.
C     ***************
100   IF(.NOT. NOR .EQ. 2) GOTO 200
C
      WRITE(NOUT ,2000)
2000  FORMAT(4X,40HTHE FOLLOWING EIGENVALUES WERE ACCEPTED.)
C
      IF(.NOT. CNEG+CPOS .GT. 0) GOTO 140
        WRITE(NOUT ,2010)
2010    FORMAT(/6X,20HLESS THAN THE SHIFT.,
     +         9X,23HGREATER THAN THE SHIFT.)
C
        SACNEG=CNEG
        SACPOS=CPOS
        IMAX=MAX0(CNEG,CPOS)
C
        DO 130 I=1,IMAX
          IF(.NOT. (SACNEG .GT. 0 .AND. SACPOS .GT. 0)) GOTO 110
            POINTN=POINTR(I)
            POINTP=I+CNEG
            POINTP=POINTR(POINTP)
C
            WRITE(NOUT ,2020) LAMBDA(POINTN),LAMBDA(POINTP)
2020        FORMAT(1X,1P1D25.7,1P1D30.7)                                F MOD.
C
            SACNEG=SACNEG-1
            SACPOS=SACPOS-1
C
            GOTO 130
C
110       IF(.NOT. SACNEG .GT. 0) GOTO 120
            POINTN=POINTR(I)
            WRITE(NOUT ,2030) LAMBDA(POINTN)
2030        FORMAT(1X,1P1D25.7)                                         F MOD.
C
            GOTO 130
C
120       POINTP=I+CNEG
          POINTP=POINTR(POINTP)
C
          WRITE(NOUT ,2040) LAMBDA(POINTP)
2040      FORMAT(1X,1P1D55.7)                                           F MOD.
C
130       CONTINUE
C
140   IF(REST .EQ. 0) WRITE(NOUT ,2050)
2050  FORMAT(/4X,25HWE HAVE NOW COMPUTED ALL ,
     +       28HEIGENVALUES IN THE INTERVAL.)
C
      IF(.NOT. REST .GT. 0) GOTO 9999
        WRITE(NOUT ,2060) REST
2060    FORMAT(/4X,15HTHERE REMAIN(S),I5,
     +         31H EIGENVALUE(S) IN THE INTERVAL./
     +         4X,29HFUNCTION SELDOM WILL BE USED.)
C
        WRITE(NOUT ,2000)
        WRITE(NOUT ,2010)
C
        GOTO 9999
C
C     ***************
C     ROUTINE SELDOM.
C     ***************
200   IF(.NOT. NOR .EQ. 19) GOTO 9999
      IF(.NOT. CNEG .GT. 0) GOTO 9999
        DO 210 I=1,CNEG
          POINTN=POINTR(I)
          WRITE(NOUT ,2030) LAMBDA(POINTN)
210       CONTINUE
C
        IF(REST .EQ. 0) WRITE(NOUT ,2070)
2070    FORMAT(/4X,40HWE HAVE NOW COMPUTED THE REMAINING ONES.)
C
C
9999  RETURN
      END
CWRINF4
      SUBROUTINE WRINF4(NOR,LOC,ALPHA,BETA,BETAPI,LAMBDA,
     +                  NU,POINTR,S,PMAX)
C
C **********************************************************************
C
C     PURPOSE -
C
C         WRITES INFORMATION FOR MSGLVL=4.
C
C     INPUT PARAMETERS -
C
C         NOR = NUMBER OF THE ROUTINE CALLING WRINFO.
C         LOC = THE LOCATION IN THAT ROUTINE.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      INTEGER   PMAX                                                    I***
      DOUBLE PRECISION                                                  R INS.
     +          STIME1(8), STIME2(8), STIME3(8),                        R MOD.
     +          A, ALPHA(1), ALTMU, B, BETA(1), BETAPI(1), CK, CL,
     +          COEFF, FACTOR, FTCONV, LAMBDA(1), SRELPR, MU, NEXTMU,
     +          NU(1), OLDMU, RDUMP, S(PMAX,1), SECOND, SMALL, TEMP,
     +          TIME, TIMTQL, TLDL, TOLBPI, TOLLDL, TOLPDM, TOLS1I,
     +          TOLZBT, TOLZNU, TOPINV, TOPM, TPRED, TSAVE, TSTLM,
     +          TVECOP, WRR
      INTEGER   ACTIVE, ADDRSS, CNEG, CNEGF,                            I***
     +          CONV, COPT, COUNT, CPOS, D, DAFILE, DUMMY, ERRNO, I,
     +          IDUMP, ITERNO, J, J1, J2, J3, K, KFILE, LEFTP, LENADR,
     +          LOC, LP, M, MAXIW, MAXL, MAXREC, MAXW, MSGLVL, MV,
     +          MXNEW, MXOLD, MXREST, MXRST, N, NBADMU, NIL, NMXRST,
     +          NOR, NREAD, NUMEIG, NUMVEC, NWRITE, OLCPOS, P,
     +          PERCNT(8), PFCONV, PMAX1, POINTR(1), POPT, PVADR, REST,
     +          RFIRST, RIGHTC, RIGHTM, RIGHTP, RNEW, ROLD, SCPX,
     +          SOLCPX, STADEW, TCONV, NERR  , NOUT , V, VER, WAD, WRI,
     +          X
      LOGICAL   DIAGM, MEQI, SAFRST, UPDATE, USEMX, WRL, ZERBET         L***
      COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
      COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,
     +                  ROLD, REST, P, USEMX, ZERBET
      COMMON   /STLMDS/ K, M, D, DIAGM
      COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
      COMMON   /STLMEW/ LEFTP, LENADR, MAXREC, NUMVEC, RIGHTC, RIGHTM,
     +                  RIGHTP, STADEW
      COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
      COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
      COMMON   /STLMIN/ A, B, NUMEIG, MAXW, MAXIW
      COMMON   /STLMMI/ MEQI
      COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
      COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,
     +                  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
      COMMON   /STLMPL/ PMAX1, POPT, COPT, MXREST
      COMMON   /STLMPR/ MSGLVL, NERR  , NOUT
      COMMON   /STLMPV/ PVADR(8)
      COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
      COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,
     +                  TOLZNU, TOLLDL(3)
      COMMON   /STLMVR/ SMALL, VER
      COMMON   /STLMWH/ WAD(2), ACTIVE(2)
      COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
C
CC      DATA      S1I,BLANK,BPI /1HS,1H ,1HB/

      CHARACTER(LEN=1) COMM1, COMM2
      CHARACTER(LEN=1), PARAMETER:: s1I="S", BLANK=" ", BPI="B"
C
C
C
      IF(.NOT. NOR .EQ. 20) GOTO 100
      GOTO(10,20,25,30,40),LOC
C
C     *******************************
C     ROUTINE STLM. INPUT PARAMETERS.
C     *******************************
10    WRITE(NOUT ,1000)
1000  FORMAT(1H1///1X,120(1H*)/1X,120(1H*)/1X,120(1H*)/
     +       2(1X,3H***,114X,3H***/),1X,3H***,
     +       38X,38H I N P U T        T O        S T L M  ,38X,3H***
     +       /2(1X,3H***,114X,3H***/),1X,120(1H*)/1X,120(1H*)/
     +       1X,120(1H*))
C
      WRITE(NOUT ,1010) A,ADDRSS,PVADR(1),B,PVADR(2),PVADR(3),
     +                  COPT,D,DAFILE,DIAGM,FACTOR,K,KFILE,PVADR(4),
     +                  LENADR,LP,M,SRELPR,MAXIW,MAXL,MAXREC,MAXW,MEQI,
     +                  MSGLVL
C
      WRITE(NOUT ,1020) MXREST,N,PVADR(5),NUMEIG,NUMVEC,PMAX,
     +                  PVADR(6),POPT,PVADR(7),PVADR(8),SMALL,STADEW,
     +                  TOLBPI,TOLLDL(1),TOLLDL(2),TOLPDM,TOLS1I,
     +                  TOLZBT,TOLZNU,NOUT ,NERR  ,VER,WAD(1),WAD(2)
C
1010  FORMAT(///
     +       6X,8HA      = ,1P1D17.7,5X,8HADDRSS = ,I17     ,           F MOD.
     +       5X,8HALPHA  = ,I17     ,5X,8HB      = ,1P1D17.7//          F MOD.
     +       6X,8HBETA   = ,I17     ,5X,8HBETAPI = ,I17     ,
     +       5X,8HCOPT   = ,I17     ,5X,8HD      = ,I17     //
     +       6X,8HDAFILE = ,I17     ,5X,8HDIAGM  = ,L17     ,
     +       5X,8HFACTOR = ,1P1D17.7,5X,8HK      = ,I17     //          F MOD.
     +       6X,8HKFILE  = ,I17     ,5X,8HLAMBDA = ,I17     ,
     +       5X,8HLENADR = ,I17     ,5X,8HLP     = ,I17     //
     +       6X,8HM      = ,I17     ,5X,8HSRELPR = ,1P1D17.7,           F MOD.
     +       5X,8HMAXIW  = ,I17     ,5X,8HMAXL   = ,I17     //
     +       6X,8HMAXREC = ,I17     ,5X,8HMAXW   = ,I17     ,
     +       5X,8HMEQI   = ,L17     ,5X,8HMSGLVL = ,I17     /)
C
1020  FORMAT(6X,8HMXREST = ,I17     ,5X,8HN      = ,I17     ,
     +       5X,8HNU     = ,I17     ,5X,8HNUMEIG = ,I17     //
     +       6X,8HNUMVEC = ,I17     ,5X,8HPMAX   = ,I17     ,
     +       5X,8HPOINTR = ,I17     ,5X,8HPOPT   = ,I17     //
     +       6X,8HS      = ,I17     ,5X,8HSCR    = ,I17     ,
     +       5X,8HSMALL  = ,1P1D17.7,5X,8HSTADEW = ,I17     //          F MOD.
     +       6X,8HTOLBPI = ,1P1D17.7,5X,8HTOLLDL1= ,1P1D17.7,           F MOD.
     +       5X,8HTOLLDL2= ,1P1D17.7,5X,8HTOLPDM = ,1P1D17.7//          F MOD.
     +       6X,8HTOLS1I = ,1P1D17.7,5X,8HTOLZBT = ,1P1D17.7,           F MOD.
     +       5X,8HTOLZNU = ,1P1D17.7,5X,8HNOUT   = ,I17     //          F MOD.
     +       6X,8HNERR   = ,I17     ,5X,8HVER    = ,I17     ,
     +       5X,8HWAD(1) = ,I17     ,5X,8HWAD(2) = ,I17
     +       //1X,13H*** INPUT ***)
C
      GOTO 9999
C
C     *******************
C     START OF ITERATION.
C     *******************
20    IF(ITERNO .EQ. 1) WRITE(NOUT ,1110)
1110  FORMAT(1H1)
C
      WRITE(NOUT ,1120) ITERNO
1120  FORMAT(1X,120(1H=)/1X,56(1H=),8H ITERNO ,56(1H=)/1X,
     +       120(1H=)/1X,I61)
C
      IF(ITERNO .GT. 1) WRITE(NOUT ,1121) TCONV,OLCPOS,ROLD,OLDMU
1121  FORMAT(//16X,5HTCONV,14X,6HOLCPOS,16X,4HROLD,15X,5HOLDMU/
     +       1X,3I20,1P1D20.7)                                          F MOD.
C
      GOTO 9999
C
C     *************
C     AFTER DECOMP.
C     *************
25    WRITE(NOUT ,1122) MU,REST,RNEW
1122  FORMAT(/1X,4HSTLM,14X,2HMU,16X,4HREST,16X,4HRNEW/
     +       1X,1P1D20.7,2I20)                                          F MOD.
C
      IF(VER .LT. 3) WRITE(NOUT , 1271) LEFTP, RIGHTC, RIGHTM, RIGHTP
C
      GOTO 9999
C
C     *******************************
C     BEFORE THE END OF AN ITERATION.
C     *******************************
30    TEMP=SECOND()-TSAVE
      WRITE(NOUT ,1130) SCPX,TCONV,TEMP
1130  FORMAT(/1X,4HSTLM,12X,4HSCPX,15X,5HTCONV,8X,12HSECOND-TSAVE/
     +       1X,2I20,1P1D20.7//1X,13H*** SHIFT ***)                     F MOD.
C
      GOTO 9999
C
C     ******************
C     OUTPUT PARAMETERS.
C     ******************
40    WRITE(NOUT ,1140)
1140  FORMAT(1H1///1X,120(1H*)/1X,120(1H*)/1X,120(1H*)/
     +       2(1X,3H***,114X,3H***/),1X,3H***,
     +       35X,44H O U T P U T        F R O M        S T L M  ,35X,
     +       3H***/2(1X,3H***,114X,3H***/),1X,120(1H*)/1X,120(1H*)/
     +       1X,120(1H*))
C
      WRITE(NOUT ,1150) ITERNO,TCONV,CNEG,CPOS,OLCPOS,RNEW,
     +                  ROLD,REST,P,USEMX,ZERBET
1150  FORMAT(///1X,8H/STLMCT//15X,6HITERNO,15X,5HTCONV,
     +       16X,4HCNEG,16X,4HCPOS,14X,6HOLCPOS,16X,4HRNEW/1X,
     +       6I20//17X,4HROLD,16X,4HREST,19X,1HP,15X,5HUSEMX,
     +       14X,6HZERBET/1X,3I20,2L20)
C
      WRITE(NOUT ,1160) ERRNO
1160  FORMAT(/1X,8H/STLMER//16X,5HERRNO/1X,I20)
C
      IF(VER .LT. 3) WRITE(NOUT , 1165) LEFTP, RIGHTC, RIGHTM, RIGHTP
1165  FORMAT(/1X,8H/STLMEW//16X,5HLEFTP,14X,6HRIGHTC,14X,6HRIGHTM,
     +        14X,6HRIGHTP/1X,4I20)
C
      WRITE(NOUT ,1170) CNEGF,RFIRST
1170  FORMAT(/1X,8H/STLMFT//16X,5HCNEGF,14X,6HRFIRST/1X,2I20)
C
      WRITE(NOUT ,1175) SCPX,SOLCPX
1175  FORMAT(/1X,8H/STLMID//17X,4HSCPX,14X,6HSOLCPX/1X,2I20)
C
      WRITE(NOUT ,1180) MU,OLDMU,NEXTMU,ALTMU
1180  FORMAT(/1X,8H/STLMMU//19X,2HMU,15X,5HOLDMU,
     +       14X,6HNEXTMU,15X,5HALTMU/1X,1P4D20.7)                      F MOD.
C
      WRITE(NOUT ,1190) ACTIVE
1190  FORMAT(/1X,8H/STLMWH//12X,9HACTIVE(1),11X,9HACTIVE(2)/1X,2I20)
C
      WRITE(NOUT ,1200)
1200  FORMAT(1H1////1X,25HINFORMATION FROM /STLMST/)
C
      TSTLM=TIME(20)
      FTCONV=TCONV
      IF(TCONV .EQ. 0) FTCONV=1.0D0                                     D MOD.
      IF(.NOT. UPDATE) GOTO 45
C       ************************************************************
C       IF THE ROUTINE IS ONLY TIMED IN FRSTIT,
C       TOTAL TIME FOR THE ROUTINE = COUNT * TIME FOR ONE REFERENCE.
C       ************************************************************
        TEMP=COUNT(12)
        IF(TIME(12) .LE. 1.0D-20) TIME(12)=TEMP*TVECOP                  D MOD.
        TEMP=COUNT(13)
        IF(TIME(13) .LE. 1.0D-20) TIME(13)=TEMP*TOPINV                  D MOD.
        TEMP=COUNT(14)
        IF(TIME(14) .LE. 1.0D-20) TIME(14)=TEMP*TOPM                    D MOD.
        TEMP=COUNT(18)
        IF(TIME(18) .LE. 1.0D-20) TIME(18)=TEMP*TVECOP                  D MOD.
        TEMP=COUNT(21)
        IF(TIME(21) .LE. 1.0D-20) TIME(21)=TEMP*TVECOP                  D MOD.
45    CONTINUE
C
      DO 60 I=1,3
C       *************************
C       WRITE IN GROUPS OF EIGHT.
C       *************************
        J1=1+(I-1)*8
        J2=J1+7
C
        IF(I .EQ. 1) WRITE(NOUT ,1210)
1210    FORMAT(32X,5HALLOC,6X,6HALWAYS,6X,6HDECOMP,
     +         7X,5HERROR,6X,6HFRSTIT,6X,6HIMTQL2,
     +         6X,6HINITLA,7X,5HINITU)
C
        IF(I .EQ. 2) WRITE(NOUT ,1220)
1220    FORMAT(//35X,2HIO,6X,6HLANCZO,8X,4HLDLT,6X,6HMULVEC,
     +         7X,5HOPINV,9X,3HOPM,6X,6HPREPSV,6X,6HRANDVC)
C
        IF(I .EQ. 3) WRITE(NOUT ,1230)
1230    FORMAT(//31X,6HSAVEXL,6X,6HSCPROD,6X,6HSELDOM,8X,4HSTLM,
     +         6X,6HSUBVEC,6X,6HTRANSF,6X,6HTRIDIG,2X,10H(REORTHO.))
C
        WRITE(NOUT ,1240) (COUNT(J),J=J1,J2)
1240    FORMAT(1X,5HCOUNT,19X,8I12)
C
C       ************************
C       HERE WE COMPUTE
C       COUNT = NUMBER OF CALLS.
C       TIME  = USED TIME.
C       TIME / COUNT.
C       TIME / TIME FOR STLM.
C       TIME / TCONV.
C
C       REORTHO. IS THE
C       REORTHOGONALIZATION IN
C       LANCZOS.
C       ************************
        DO 50 J=J1,J2
          J3=J-J1+1
          STIME1(J3)=TIME(J)
          TEMP=COUNT(J)
          STIME2(J3)=0.0D0                                              D MOD.
          IF(COUNT(J) .NE. 0) STIME2(J3)=TIME(J)/TEMP
          PERCNT(J3)=(TIME(J)/TSTLM)*1.0D2+0.5D0                        D MOD.
          STIME3(J3)=TIME(J)/FTCONV
50        CONTINUE
C
        WRITE(NOUT ,1250) STIME1,STIME2,PERCNT,STIME3
1250    FORMAT(1X,4HTIME,20X,1P8D12.1/1X,10HTIME/COUNT,14X,1P8D12.1/    F MOD.
     +         1X,17HTIME/STLM IN P.C.,7X,8I12/
     +         1X,10HTIME/TCONV,14X,1P8D12.1)                           F MOD.
C
60      CONTINUE
C
      WRITE(NOUT ,1260) NBADMU,NMXRST,NREAD,NWRITE,WRI(1),WRI(2),WRR(1)
1260  FORMAT(/15X,6HNBADMU,14X,6HNMXRST,
     +       15X,5HNREAD,14X,6HNWRITE/1X,4I20//
     +       15X,6HNTCONV,15X,5HNLEFT,18X,2HBG/
     +       1X,2I20,1P1D20.7/                                          F MOD.
     +       //1X,14H*** OUTPUT ***//1X,12H*** STLM ***)
C
      GOTO 9999
C
C     ***************
C     ROUTINE ALWAYS.
C     ***************
100   IF(.NOT. NOR .EQ. 2) GOTO 200
C
      WRITE(NOUT ,1270) RDUMP,WRL(2),REST,WRR(1)
1270  FORMAT(/1X,6HALWAYS,11X,3HDOT,14X,6HPOSDOT,
     +       16X,4HREST,15X,5HTNORM/1X,1P1D20.7,L20,I20,1P1D20.7)       F MOD.
C
      IF(VER .LT. 3) WRITE(NOUT , 1271) LEFTP, RIGHTC, RIGHTM, RIGHTP
1271  FORMAT(/16X,5HLEFTP,14X,6HRIGHTC,14X,6HRIGHTM,
     +        14X,6HRIGHTP/1X,4I20)
C
      GOTO 9999
C
C     ***************
C     ROUTINE DECOMP.
C     ***************
200   IF(.NOT. NOR .EQ. 3) GOTO 300
C
      WRITE(NOUT ,1280) MU,RNEW,NBADMU,NMXRST,REST
1280  FORMAT(/1X,6HDECOMP,12X,2HMU,16X,4HRNEW,14X,6HNBADMU,
     +       14X,6HNMXRST,16X,4HREST/1X,1P1D20.7,4I20)                  F MOD.
C
      GOTO 9999
C
C     ***************
C     ROUTINE FRSTIT.
C     ***************
300   IF(.NOT. NOR .EQ. 5) GOTO 400
C
      WRITE(NOUT ,1290) NIL,MV,V,MXNEW,MXOLD,MXRST,X,POPT,COPT,MXREST
1290  FORMAT(/1X,6HFRSTIT,11X,3HNIL,18X,2HMV,19X,1HV,15X,5HMXNEW,
     +       15X,5HMXOLD,15X,5HMXRST/1X,6I20//20X,1HX,16X,4HPOPT,
     +       16X,4HCOPT,14X,6HMXREST/1X,4I20)
C
      WRITE(NOUT ,1300) TLDL,TOPINV,TOPM,TIMTQL,TVECOP,TPRED,
     +                  TSAVE,COEFF,CK,CL,UPDATE,USEMX
1300  FORMAT(/17X,4HTLDL,14X,6HTOPINV,16X,4HTOPM,14X,6HTIMTQL,
     +       14X,6HTVECOP,15X,5HTPRED/1X,1P6D20.7//16X,5HTSAVE,         F MOD.
     +       12X,8HCOEFF(1),12X,8HCOEFF(2),12X,8HCOEFF(3),
     +       12X,8HCOEFF(4)/1X,1P5D20.7//19X,2HCK,18X,2HCL,             F MOD.
     +       14X,6HUPDATE,15X,5HUSEMX/1X,1P2D20.7,2L20)                 F MOD.
C
      WRITE(NOUT ,1305) (WRR(I),I=1,3),(WRI(I),I=1,4)
1305  FORMAT(/20X,1HA,19X,1HB,19X,1HC,
     +       18X,2HTC,18X,2HTP,16X,4HLOOP/1X,1P3D20.7,3I20//18X,3HLIM/  F MOD.
     +       1X,I20)
C
      GOTO 9999
C
C     ***************
C     ROUTINE TRIDIG.
C     ***************
400   IF(.NOT. NOR .EQ. 23) GOTO 500
      GOTO(410,430),LOC
C
410   WRITE(NOUT ,1310) CNEG,CPOS,CONV,PFCONV,P,ZERBET
1310  FORMAT(/1X,6HTRIDIG,10X,4HCNEG,16X,4HCPOS,16X,4HCONV,
     +       14X,6HPFCONV,19X,1HP,14X,6HZERBET/1X,5I20,L20//
     +       5X,1HI,1X,4HCOMM,11X,9HLAMBDA(I),1X,9HBETAPI(I),
     +       15X,5HNU(I),4X,6HS(1,I),4X,6HS(P,I),12X,8HALPHA(I),
     +       13X,7HBETA(I))
C
      WRITE(NOUT ,1330)
C
      DO 420 I=1,P
        COMM1=BLANK
        COMM2=BLANK
        IF(DABS(BETAPI(I)) .LE. TOLBPI) COMM1=BPI                       D MOD.
        IF(DABS(S(1,I)) .LE. TOLS1I) COMM2=S1I                          D MOD.
C
        WRITE(NOUT ,1320) I,COMM1,COMM2,LAMBDA(I),BETAPI(I),NU(I),
     +                    S(1,I),S(P,I),ALPHA(I),BETA(I)
1320    FORMAT(1X,I5,3X,2A1,2(1P1D20.7,1P1D10.2),1P1D10.2,1P2D20.7)     F MOD.
C
        IF(MOD(I,5) .EQ. 0) WRITE(NOUT ,1330)
1330    FORMAT(1X,120(1H-))
C
420     CONTINUE
C
      GOTO 9999
C
430   WRITE(NOUT ,1340) CONV,CNEG,CPOS,WRI(1),WRI(2),
     +                  WRR(1),WRR(2),WRL(1),NEXTMU,ALTMU
1340  FORMAT(/1X,6HTRIDIG,10X,4HCONV,16X,4HCNEG,16X,4HCPOS,
     +       14X,6HNUMDEL,14X,6HS1IDEL/1X,5I20//16X,5HMINNU,
     +       15X,5HMAXNU,15X,5HPOSNU,14X,6HNEXTMU,15X,5HALTMU/
     +       1X,1P2D20.7,L20,1P2D20.7)                                  F MOD.
C
      I=CNEG+CPOS
C
      IF(I .GT. 0) WRITE(NOUT ,1350) (POINTR(J),J=1,I)
1350  FORMAT(/1X,6HPOINTR/(1X,30I4))
C
      GOTO 9999
C
C     **************
C     ROUTINE NEWPC.
C     **************
500   IF(.NOT. NOR .EQ. 25) GOTO 600
      GOTO(510,520),LOC
C
510   WRITE(NOUT ,1360) (WRR(I),I=1,4),(WRI(I),I=1,4)
1360  FORMAT(/1X,5HNEWPC,14X,1HA,19X,1HB,19X,1HC,14X,6HWEIGHT,
     +       18X,2HTC,18X,2HTP/1X,1P4D20.7,2I20//17X,4HLOOP,17X,3HLIM/  F MOD.
     +       1X,2I20)
C
      GOTO 9999
C
520   TEMP=SECOND()-TSAVE
      WRITE(NOUT ,1370) (WRR(I),I=1,2),COPT,POPT,MXREST,CL,CK,TPRED,
     +                  TEMP
1370  FORMAT(/20X,1HL,19X,1HK,16X,4HCOPT,16X,4HPOPT,14X,6HMXREST,
     +       18X,2HCL/1X,1P2D20.7,3I20,1P1D20.7//19X,2HCK,15X,5HTPRED,  F MOD.
     +       8X,12HSECOND-TSAVE/1X,1P3D20.7)                            F MOD.
C
      GOTO 9999
C
C     ***************
C     ROUTINE SELDOM.
C     ***************
600   IF(.NOT. NOR .EQ. 19) GOTO 9999
C
      WRITE(NOUT ,1380) REST,WRI(1),RDUMP,WRL(2),WRR(1)
1380  FORMAT(/1X,6HSELDOM,10X,4HREST,16X,4HNUMR,
     +       17X,3HDOT,14X,6HPOSDOT,15X,5HTNORM/
     +       1X,2I20,1P1D20.7,L20,1P1D20.7)                             F MOD.
C
C
9999  RETURN
      END
CWRINFO
      SUBROUTINE WRINFO(NOR, LOC, W, IW)
C
C **********************************************************************
C
C     PURPOSE -
C
C         DRIVER ROUTINE FOR THE ROUTINES THAT WRITE INFORMATION.
C
C
C     INPUT PARAMETERS -
C
C         NOR = NUMBER OF THE CALLING ROUTINE.
C         LOC = LOCATION IN THAT ROUTINE.
C
C     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
C     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
C     THE FUNCTION OF THE ROUTINE.
C
C
C **********************************************************************
C
      DOUBLE PRECISION                                                  R INS.
     +          W(1)                                                    R MOD.
      INTEGER   ALPHA, BETA, BETAPI, COPT, IW(1), LAMBDA, LOC, MSGLVL,  I***
     +          MXREST, NOR, NU, PMAX, POINTR, POPT, S, SCR, NERR  ,
     +          NOUT
      COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
      COMMON   /STLMPR/ MSGLVL, NERR  , NOUT
      COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
C
      IF(MSGLVL .LE. 1) GOTO 9999
C
      IF(.NOT. (MSGLVL .EQ. 3 .OR. (MSGLVL .EQ. 2 .AND.
     +   NOR .EQ. 20 .AND. (LOC .EQ. 1 .OR. LOC .EQ. 5)))) GOTO 10
C
        IF(NOR .EQ. 20 .OR. NOR .EQ. 2 .OR. NOR .EQ. 19)
     +                   CALL WRINF3(NOR, LOC, W(LAMBDA), IW(POINTR))
C
10    IF(MSGLVL .EQ. 4)
     +  CALL WRINF4(NOR, LOC, W(ALPHA), W(BETA), W(BETAPI),
     +              W(LAMBDA), W(NU), IW(POINTR), W(S), PMAX)
C
9999  RETURN

      END





      REAL*8 FUNCTION V(XI)

C----------------------------------------------------------------

C Declare Variables

      IMPLICIT REAL*8 (A-H,O-Z)

C----------------------------------------------------------------
C     #aqui2

      INCLUDE 'Constantes2.txt'
      XE= R0
      r_m = r_m*1822.88839
      xx=XI-XE
      
      V=-De*(1+c1*xx+c2*xx**2+c3*xx**3+c4*xx**4+c5*xx**5+c6*xx**6)*

     >exp(-c1*xx)+ De

     >+(J*(J+1.0d0))/(2*r_m*(XI**2))

      if (mod (XI, 1.).NE.0) write (98,*)XI,V
      
      write(99,*)XI,V

      RETURN

      END
