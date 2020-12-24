C $Header: /u/gcmpack/MITgcm/model/src/remove_mean.F,v 1.3 2009/04/28 18:30:33 jmc Exp $
C $Name: checkpoint63a $

#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: REMOVE_MEAN_RL
C     !INTERFACE:
      SUBROUTINE REMOVE_MEAN_RL(
     I                myNr, arr, arrMask, arrhFac, arrArea, arrDr,
     I                arrName, myTime,
     I                myThid )
C     !DESCRIPTION: \bv
C     /==========================================================\
C     | SUBROUTINE REMOVE_MEAN_RL                                |
C     | o Calculate mean of global array "_RL arr" and substract |
C     |   it from the array                                      |
C     \==========================================================/
C     \ev

      IMPLICIT NONE

C     === Global data ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C     === Functions ====
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE

C     === Routine arguments ===
      INTEGER myNr
      _RL arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,myNr,nSx,nSy)
      _RS arrMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,myNr,nSx,nSy)
      _RS arrhFac(1-OLx:sNx+OLx,1-OLy:sNy+OLy,myNr,nSx,nSy)
      _RS arrArea(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS arrDr(myNr)
      CHARACTER*(*) arrName
      _RL myTime
      INTEGER myThid

C#ifdef ALLOW_BALANCE_FLUXES
C     === Local variables ====
      INTEGER bi,bj,I,J,K
      _RL tmpVal
      _RL theMean
      _RL theVol
      _RL tmpVol
      CHARACTER*(max_len_mbuf) msgbuf
CEOP

      theMean=0.
      theVol=0.

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO K=1,myNr
         DO J=1,sNy
          DO I=1,sNx
           tmpVal=arr(I,J,K,bi,bj)
           IF (arrMask(I,J,K,bi,bj).NE.0.) THEN
            tmpVol = arrArea(I,J,bi,bj)*arrhFac(I,J,K,bi,bj)*arrDr(K)
            theVol = theVol   + tmpVol
            theMean = theMean + tmpVol*tmpVal
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      _GLOBAL_SUM_RL(theVol,myThid)
      _GLOBAL_SUM_RL(theMean,myThid)

      IF (theVol.GT.0.) THEN
       theMean=theMean/theVol

       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO K=1,myNr
          DO J=1,sNy
           DO I=1,sNx
            IF (arrMask(I,J,K,bi,bj).NE.0.) THEN
             arr(I,J,K,bi,bj) = arr(I,J,K,bi,bj) - theMean
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO

      ENDIF

C     Print the global mean to standard output, this is a measure for
C     the drift of the array arr
      IF ( balancePrintMean ) THEN
       _BEGIN_MASTER( myThid )
       WRITE(msgbuf,'(a,a,a,e24.17)')
     &      'REMOVE_MEAN_RL: Global mean of ',
     &      arrName, ' = ', theMean
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &      SQUEEZE_RIGHT , 1)
       _END_MASTER( myThid )
      ENDIF

C#endif /* ALLOW_BALANCE_FLUXES */
      RETURN
      END

CBOP
C     !ROUTINE: REMOVE_MEAN_RS
C     !INTERFACE:
      SUBROUTINE REMOVE_MEAN_RS(
     I                myNr, arr, arrMask, arrhFac, arrArea, arrDr,
     I                arrName, myTime,
     I                myThid )
C     !DESCRIPTION: \bv
C     /==========================================================\
C     | SUBROUTINE REMOVE_MEAN_RS                                |
C     | o Calculate mean of global array "_RS arr" and substract |
C     |   it from the array                                      |
C     \==========================================================/
C     \ev

      IMPLICIT NONE

C     === Global data ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C     === Functions ====
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE

C     === Routine arguments ===
      INTEGER myNr
      _RS arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,myNr,nSx,nSy)
      _RS arrMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,myNr,nSx,nSy)
      _RS arrhFac(1-OLx:sNx+OLx,1-OLy:sNy+OLy,myNr,nSx,nSy)
      _RS arrArea(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS arrDr(myNr)
      CHARACTER*(*) arrName
      _RL myTime
      INTEGER myThid

C#ifdef ALLOW_BALANCE_FLUXES
C     === Local variables ====
      INTEGER bi,bj,I,J,K
      _RS tmpVal
      _RS theMean
      _RS theVol

      _RS theMean_A !Atlantic BOX
      _RS theMean_M !Mediterranean
      _RS theVol_A
      _RS theVol_M


      _RS tmpVol
      CHARACTER*(max_len_mbuf) msgbuf
cgmBALANCE(
      INTEGER iG, jG, iGo, jGo
cgmBALANCE)

CEOP

      theMean=0.
      theVol=0.

      theMean_A=0.
      theVol_A=0.

      theMean_M=0.
      theVol_M=0.


cgmOriginal(
c      DO bj=myByLo(myThid),myByHi(myThid)
c       DO bi=myBxLo(myThid),myBxHi(myThid)
c        DO K=1,myNr
c         DO J=1,sNy
c          DO I=1,sNx
c           tmpVal=arr(I,J,K,bi,bj)
c           IF (arrMask(I,J,K,bi,bj).NE.0.) THEN
c            tmpVol = arrArea(I,J,bi,bj)*arrhFac(I,J,K,bi,bj)*arrDr(K)
c            theVol = theVol   + tmpVol
c            theMean = theMean + tmpVol*tmpVal
c           ENDIF
c          ENDDO
c         ENDDO
c        ENDDO
c       ENDDO
c      ENDDO
cgmOriginal)
cgmBALANCE(
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO K=1,myNr
         DO J=1,sNy
          DO I=1,sNx
C convert from local y index J to global y index jG
c                jG = myYGlobalLo-1+(bj-1)*sNy+J
                iG = myXGlobalLo-1+(bi-1)*sNx+I

cgm(it needs to be modified also few lines down from here
         IF (iG.le.67) then  !1/12
cgm)
           tmpVal=arr(I,J,K,bi,bj)
           IF (arrMask(I,J,K,bi,bj).NE.0.) THEN
            tmpVol = arrArea(I,J,bi,bj)*arrhFac(I,J,K,bi,bj)*arrDr(K)
            theVol_A = theVol_A   + tmpVol
            theMean_A = theMean_A + tmpVol*tmpVal
           ENDIF
	ENDIF
c------
         IF (iG.gt.67) then  !1/12
cgm)
           tmpVal=arr(I,J,K,bi,bj)
           IF (arrMask(I,J,K,bi,bj).NE.0.) THEN
            tmpVol = arrArea(I,J,bi,bj)*arrhFac(I,J,K,bi,bj)*arrDr(K)
            theVol_M = theVol_M   + tmpVol
            theMean_M = theMean_M + tmpVol*tmpVal
           ENDIF
	ENDIF
c------

          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
cgmBALANCE)
      _GLOBAL_SUM_RS(theVol_A,myThid)
      _GLOBAL_SUM_RS(theMean_A,myThid)
      _GLOBAL_SUM_RS(theVol_M,myThid)
      _GLOBAL_SUM_RS(theMean_M,myThid)

c      _GLOBAL_SUM_RS(theVol,myThid)
c      _GLOBAL_SUM_RS(theMean,myThid)

	
      IF (theVol_A.GT.0.) THEN
       theMean=(theMean_A+ theMEAN_M)/theVol_A

cgmBalance(
cgmOriginal(
c       DO bj=myByLo(myThid),myByHi(myThid)
c        DO bi=myBxLo(myThid),myBxHi(myThid)
c         DO K=1,myNr
c          DO J=1,sNy
c           DO I=1,sNx
c            IF (arrMask(I,J,K,bi,bj).NE.0.) THEN
c             arr(I,J,K,bi,bj) = arr(I,J,K,bi,bj) - theMean
c            ENDIF
c           ENDDO
c          ENDDO
c         ENDDO
c        ENDDO
c       ENDDO
cgmOriginal)
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO K=1,myNr
          DO J=1,sNy
           DO I=1,sNx
C convert from local y index J to global y index jG
c                jG = myYGlobalLo-1+(bj-1)*sNy+J
                iG = myXGlobalLo-1+(bi-1)*sNx+I
cgm(it needs to be modified also few lines upper from here
         IF (iG.le.67) then  !1/12
cgm)
            IF (arrMask(I,J,K,bi,bj).NE.0.) THEN
             arr(I,J,K,bi,bj) = arr(I,J,K,bi,bj) - theMean


		ENDIF
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO

cgm.....ricalcolo per verificare che sia zero 
c      DO bj=myByLo(myThid),myByHi(myThid)
c       DO bi=myBxLo(myThid),myBxHi(myThid)
c        DO K=1,myNr
c         DO J=1,sNy
c          DO I=1,sNx
c           tmpVal=arr(I,J,K,bi,bj)
c           IF (arrMask(I,J,K,bi,bj).NE.0.) THEN
c            tmpVol = arrArea(I,J,bi,bj)*arrhFac(I,J,K,bi,bj)*arrDr(K)
c            theVol = theVol   + tmpVol
c            theMean = theMean + tmpVol*tmpVal
c           ENDIF
c          ENDDO
c         ENDDO
c        ENDDO
c       ENDDO
c      ENDDO

c      _GLOBAL_SUM_RS(theVol,myThid)
c      _GLOBAL_SUM_RS(theMean,myThid)

c      IF (theVol.GT.0.) THEN
c       theMean=theMean/theVol
c	ENDIF
c	print*,'CONFERMA.....:',theMean
c----------------------
cgmBALANCE)

      ENDIF

C     Print the global mean to standard output, this is a measure for
C     the drift of the array arr
      IF ( balancePrintMean ) THEN
       _BEGIN_MASTER( myThid )
       WRITE(msgbuf,'(a,a,a,e24.17)')
     &      'REMOVE_MEAN_RS: Global mean of ',
     &      arrName, ' = ', theMean
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &      SQUEEZE_RIGHT , 1)
       _END_MASTER( myThid )
      ENDIF

C#endif /* ALLOW_BALANCE_FLUXES */
      RETURN
      END
