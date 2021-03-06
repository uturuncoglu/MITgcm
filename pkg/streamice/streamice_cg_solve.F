C $Header: /u/gcmpack/MITgcm/pkg/streamice/streamice_cg_solve.F,v 1.5 2014/07/09 16:09:23 dgoldberg Exp $
C $Name: checkpoint65 $

#include "STREAMICE_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
      SUBROUTINE STREAMICE_CG_SOLVE( 
     U                               cg_Uin,     ! x-velocities
     U                               cg_Vin,     ! y-velocities
     I                               cg_Bu,      ! force in x dir
     I                               cg_Bv,      ! force in y dir
     I                               A_uu,       ! section of matrix that multiplies u and projects on u
     I                               A_uv,       ! section of matrix that multiplies v and projects on u
     I                               A_vu,       ! section of matrix that multiplies u and projects on v
     I                               A_vv,       ! section of matrix that multiplies v and projects on v
     I                               tolerance, 
     O                               iters,
     I                               myThid )
C     /============================================================\
C     | SUBROUTINE                                                 |   
C     | o                                                          |
C     |============================================================|
C     |                                                            |
C     \============================================================/
      IMPLICIT NONE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "STREAMICE.h"
#include "STREAMICE_CG.h"



!#ifdef ALLOW_PETSC
!#include "finclude/petsc.h"
! UNCOMMENT IF V3.0
!#include "finclude/petscvec.h"
!#include "finclude/petscmat.h"
!#include "finclude/petscksp.h"
!#include "finclude/petscpc.h"
!#endif
C     === Global variables ===

      
C     !INPUT/OUTPUT ARGUMENTS
C     cg_Uin, cg_Vin - input and output velocities
C     cg_Bu, cg_Bv - driving stress
      INTEGER myThid
      INTEGER iters
      _RL tolerance
      _RL cg_Uin (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL cg_Vin (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL cg_Bu (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL cg_Bv (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL 
     & A_uu (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,-1:1,-1:1),
     & A_vu (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,-1:1,-1:1),
     & A_uv (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,-1:1,-1:1),
     & A_vv (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,-1:1,-1:1)

C     LOCAL VARIABLES
      INTEGER i, j, bi, bj, cg_halo, conv_flag
      INTEGER iter, is, js, ie, je, colx, coly, k
      _RL dot_p1, dot_p2, alpha_k, beta_k, resid, resid_0
      _RL dot_p1_tile (nSx,nSy)
      _RL dot_p2_tile (nSx,nSy)
      CHARACTER*(MAX_LEN_MBUF) msgBuf


!#ifdef ALLOW_PETSC
!      INTEGER indices(2*(snx*nsx*sny*nsy))
!      INTEGER n_dofs_cum_sum (0:nPx*nPy-1), idx(1)
!      _RL rhs_values(2*(snx*nsx*sny*nsy))
!      _RL solution_values(2*(snx*nsx*sny*nsy))
!      _RL mat_values (2*Nx*Ny,2*(snx*nsx*sny*nsy))
!      _RL mat_values (18,1), mat_val_return(1)
!      INTEGER indices_col(18)
!      INTEGER local_dofs, global_dofs, dof_index, dof_index_col
!      INTEGER local_offset
!      Mat matrix
!      KSP ksp
!      PC  pc
!      Vec rhs
!      Vec solution
!      PetscErrorCode ierr
!#ifdef ALLOW_USE_MPI
!      integer mpiRC, mpiMyWid
!#endif 
!#endif


#ifdef ALLOW_STREAMICE



      CALL TIMER_START ('STREAMICE_CG_SOLVE',myThid)
#ifndef STREAMICE_SERIAL_TRISOLVE

#ifdef ALLOW_PETSC

      CALL STREAMICE_CG_SOLVE_PETSC(
     U         cg_Uin,     ! x-velocities
     U         cg_Vin,     ! y-velocities
     I         cg_Bu,      ! force in x dir
     I         cg_Bv,      ! force in y dir
     I         A_uu,       ! section of matrix that multiplies u and projects on u
     I         A_uv,       ! section of matrix that multiplies v and projects on u
     I         A_vu,       ! section of matrix that multiplies u and projects on v
     I         A_vv,       ! section of matrix that multiplies v and projects on v
     I         tolerance,
     O         iters,
     I         myThid )


#else  /* ALLOW_PETSC */


      iters = streamice_max_cg_iter
      conv_flag = 0


      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          Zu_SI (i,j,bi,bj) = 0. _d 0
          Zv_SI (i,j,bi,bj) = 0. _d 0
          Ru_SI (i,j,bi,bj) = 0. _d 0
          Rv_SI (i,j,bi,bj) = 0. _d 0
          Au_SI (i,j,bi,bj) = 0. _d 0
          Av_SI (i,j,bi,bj) = 0. _d 0
          Du_SI (i,j,bi,bj) = 0. _d 0
          Dv_SI (i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      
C     FIND INITIAL RESIDUAL, and initialize r

! #ifdef STREAMICE_CONSTRUCT_MATRIX

        

        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j=1,sNy
           DO i=1,sNx
            DO colx=-1,1
             DO coly=-1,1
              Au_SI(i,j,bi,bj) = Au_SI(i,j,bi,bj) + 
     &         A_uu(i,j,bi,bj,colx,coly)*
     &         cg_Uin(i+colx,j+coly,bi,bj)+
     &         A_uv(i,j,bi,bj,colx,coly)*     
     &         cg_Vin(i+colx,j+coly,bi,bj)


              Av_SI(i,j,bi,bj) = Av_SI(i,j,bi,bj) + 
     &         A_vu(i,j,bi,bj,colx,coly)*
     &         cg_Uin(i+colx,j+coly,bi,bj)+
     &         A_vv(i,j,bi,bj,colx,coly)*     
     &         cg_Vin(i+colx,j+coly,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

             
      _EXCH_XY_RL( Au_SI, myThid )
      _EXCH_XY_RL( Av_SI, myThid )


      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          Ru_SI(i,j,bi,bj)=cg_Bu(i,j,bi,bj)-
     &     Au_SI(i,j,bi,bj)
          Rv_SI(i,j,bi,bj)=cg_Bv(i,j,bi,bj)-
     &     Av_SI(i,j,bi,bj)
         ENDDO
        ENDDO
        dot_p1_tile(bi,bj) = 0. _d 0
        dot_p2_tile(bi,bj) = 0. _d 0  
       ENDDO
      ENDDO

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          IF (STREAMICE_umask(i,j,bi,bj).eq.1.0)
     &      dot_p1_tile(bi,bj)=dot_p1_tile(bi,bj)+Ru_SI(i,j,bi,bj)**2
          IF (STREAMICE_vmask(i,j,bi,bj).eq.1.0)
     &      dot_p1_tile(bi,bj)=dot_p1_tile(bi,bj)+Rv_SI(i,j,bi,bj)**2
         ENDDO
        ENDDO
       ENDDO
      ENDDO


      CALL GLOBAL_SUM_TILE_RL( dot_p1_tile, dot_p1, myThid ) 
      resid_0 = sqrt(dot_p1)

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

       WRITE(msgBuf,'(A,I1,I1,E14.7)') 'CONJ GRAD INIT RESID LOCAL, ',
     &         bi,bj, dot_p1_tile(bi,bj)
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT , 1)

       enddo
      enddo

      WRITE(msgBuf,'(A,E14.7)') 'CONJ GRAD INIT RESID, ',
     &         resid_0
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT , 1)

C    CCCCCCCCCCCCCCCCCCCC

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          IF (STREAMICE_umask(i,j,bi,bj).eq.1.0)
     &      Zu_SI(i,j,bi,bj)=Ru_SI(i,j,bi,bj) / DIAGu_SI(i,j,bi,bj)
          IF (STREAMICE_vmask(i,j,bi,bj).eq.1.0)
     &      Zv_SI(i,j,bi,bj)=Rv_SI(i,j,bi,bj) / DIAGv_SI(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      cg_halo = min(OLx-1,OLy-1)
      conv_flag = 0

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          Du_SI(i,j,bi,bj)=Zu_SI(i,j,bi,bj) 
          Dv_SI(i,j,bi,bj)=Zv_SI(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      resid = resid_0
      iters = 0
      
c  !!!!!!!!!!!!!!!!!!
c  !!              !!
c  !! MAIN CG LOOP !!
c  !!              !!
c  !!!!!!!!!!!!!!!!!!


     
  
c  ! initially, b-grid data is valid up to 3 halo nodes out -- right? (check for MITgcm!!)

      WRITE(msgBuf,'(A)') 'BEGINNING MAIN CG LOOP'
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT , 1)

!       IF(STREAMICE_construct_matrix) CALL STREAMICE_CG_MAKE_A(myThid)


      do iter = 1, streamice_max_cg_iter
       if (resid .gt. tolerance*resid_0) then

c      to avoid using "exit"
       iters = iters + 1

       is = 1 - cg_halo 
       ie = sNx + cg_halo 
       js = 1 - cg_halo 
       je = sNy + cg_halo
      
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx 
           Au_SI(i,j,bi,bj) = 0. _d 0
           Av_SI(i,j,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO

!        IF (STREAMICE_construct_matrix) THEN

! #ifdef STREAMICE_CONSTRUCT_MATRIX

        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j=js,je
           DO i=is,ie
            DO colx=-1,1
             DO coly=-1,1
              Au_SI(i,j,bi,bj) = Au_SI(i,j,bi,bj) + 
     &         A_uu(i,j,bi,bj,colx,coly)*
     &         Du_SI(i+colx,j+coly,bi,bj)+
     &         A_uv(i,j,bi,bj,colx,coly)*     
     &         Dv_SI(i+colx,j+coly,bi,bj)
              Av_SI(i,j,bi,bj) = Av_SI(i,j,bi,bj) + 
     &         A_vu(i,j,bi,bj,colx,coly)*
     &         Du_SI(i+colx,j+coly,bi,bj)+
     &         A_vv(i,j,bi,bj,colx,coly)*     
     &         Dv_SI(i+colx,j+coly,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

!        else
! #else
! 
!         CALL STREAMICE_CG_ACTION( myThid, 
!      O     Au_SI, 
!      O     Av_SI, 
!      I     Du_SI, 
!      I     Dv_SI, 
!      I     is,ie,js,je)
! 
! !        ENDIF 
! 
! #endif
       

       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         dot_p1_tile(bi,bj) = 0. _d 0
         dot_p2_tile(bi,bj) = 0. _d 0  
        ENDDO
       ENDDO

       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1,sNy
          DO i=1,sNx
           IF (STREAMICE_umask(i,j,bi,bj).eq.1.0) THEN
           dot_p1_tile(bi,bj)=dot_p1_tile(bi,bj)+Zu_SI(i,j,bi,bj)*
     &        Ru_SI(i,j,bi,bj)
            dot_p2_tile(bi,bj)=dot_p2_tile(bi,bj)+Du_SI(i,j,bi,bj)*
     &        Au_SI(i,j,bi,bj)
           ENDIF
           IF (STREAMICE_vmask(i,j,bi,bj).eq.1.0) THEN
            dot_p1_tile(bi,bj)=dot_p1_tile(bi,bj)+Zv_SI(i,j,bi,bj)*
     &        Rv_SI(i,j,bi,bj)
            dot_p2_tile(bi,bj)=dot_p2_tile(bi,bj)+Dv_SI(i,j,bi,bj)*
     &        Av_SI(i,j,bi,bj)
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO

       CALL GLOBAL_SUM_TILE_RL( dot_p1_tile, dot_p1, myThid )
       CALL GLOBAL_SUM_TILE_RL( dot_p2_tile, dot_p2, myThid )            
       alpha_k = dot_p1/dot_p2

       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx

           IF (STREAMICE_umask(i,j,bi,bj).eq.1.0) THEN
            cg_Uin(i,j,bi,bj)=cg_Uin(i,j,bi,bj)+
     &       alpha_k*Du_SI(i,j,bi,bj)
            Ru_old_SI(i,j,bi,bj) = Ru_SI(i,j,bi,bj)
            Zu_old_SI(i,j,bi,bj) = Zu_SI(i,j,bi,bj)
            Ru_SI(i,j,bi,bj) = Ru_SI(i,j,bi,bj)-
     &       alpha_k*Au_SI(i,j,bi,bj)
            Zu_SI(i,j,bi,bj) = Ru_SI(i,j,bi,bj) / 
     &       DIAGu_SI(i,j,bi,bj)
           ENDIF

           IF (STREAMICE_vmask(i,j,bi,bj).eq.1.0) THEN
            cg_Vin(i,j,bi,bj)=cg_Vin(i,j,bi,bj)+
     &       alpha_k*Dv_SI(i,j,bi,bj)
            Rv_old_SI(i,j,bi,bj) = Rv_SI(i,j,bi,bj)
            Zv_old_SI(i,j,bi,bj) = Zv_SI(i,j,bi,bj)
            Rv_SI(i,j,bi,bj) = Rv_SI(i,j,bi,bj)-
     &       alpha_k*Av_SI(i,j,bi,bj)
            Zv_SI(i,j,bi,bj) = Rv_SI(i,j,bi,bj) / 
     &       DIAGv_SI(i,j,bi,bj)

           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO

       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         dot_p1_tile(bi,bj) = 0. _d 0
         dot_p2_tile(bi,bj) = 0. _d 0  
        ENDDO
       ENDDO

       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1,sNy
          DO i=1,sNx

           IF (STREAMICE_umask(i,j,bi,bj).eq.1.0) THEN
            dot_p1_tile(bi,bj)=dot_p1_tile(bi,bj)+Zu_SI(i,j,bi,bj)*
     &        Ru_SI(i,j,bi,bj)
            dot_p2_tile(bi,bj)=dot_p2_tile(bi,bj)+Zu_old_SI(i,j,bi,bj)*
     &        Ru_old_SI(i,j,bi,bj)
           ENDIF

           IF (STREAMICE_vmask(i,j,bi,bj).eq.1.0) THEN
            dot_p1_tile(bi,bj)=dot_p1_tile(bi,bj)+Zv_SI(i,j,bi,bj)*
     &        Rv_SI(i,j,bi,bj)
            dot_p2_tile(bi,bj)=dot_p2_tile(bi,bj)+Zv_old_SI(i,j,bi,bj)*
     &        Rv_old_SI(i,j,bi,bj)
           ENDIF

          ENDDO
         ENDDO
        ENDDO
       ENDDO

       CALL GLOBAL_SUM_TILE_RL( dot_p1_tile, dot_p1, myThid )
       CALL GLOBAL_SUM_TILE_RL( dot_p2_tile, dot_p2, myThid )
                   
       beta_k = dot_p1/dot_p2

       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           IF (STREAMICE_umask(i,j,bi,bj).eq.1.0) 
     &      Du_SI(i,j,bi,bj)=beta_k*Du_SI(i,j,bi,bj)+
     &      Zu_SI(i,j,bi,bj)
           IF (STREAMICE_vmask(i,j,bi,bj).eq.1.0) 
     &      Dv_SI(i,j,bi,bj)=beta_k*Dv_SI(i,j,bi,bj)+
     &      Zv_SI(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO

       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         dot_p1_tile(bi,bj) = 0. _d 0
        ENDDO
       ENDDO

       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1,sNy
          DO i=1,sNx
           IF (STREAMICE_umask(i,j,bi,bj).eq.1.0)
     &      dot_p1_tile(bi,bj)=dot_p1_tile(bi,bj)+Ru_SI(i,j,bi,bj)**2
           IF (STREAMICE_vmask(i,j,bi,bj).eq.1.0)
     &      dot_p1_tile(bi,bj)=dot_p1_tile(bi,bj)+Rv_SI(i,j,bi,bj)**2
          ENDDO
         ENDDO
        ENDDO
       ENDDO

       CALL GLOBAL_SUM_TILE_RL( dot_p1_tile, dot_p1, myThid ) 
       resid = sqrt(dot_p1)

!        IF (iter .eq. 1) then
!         print *, alpha_k, beta_k, resid
!        ENDIF

       cg_halo = cg_halo - 1

       if (cg_halo .eq. 0) then
        cg_halo = min(OLx-1,OLy-1)
        _EXCH_XY_RL( Du_SI, myThid )
        _EXCH_XY_RL( Dv_SI, myThid )
        _EXCH_XY_RL( Ru_SI, myThid )
        _EXCH_XY_RL( Rv_SI, myThid )
        _EXCH_XY_RL( cg_Uin, myThid )
        _EXCH_XY_RL( cg_Vin, myThid )
       endif


       endif
      enddo ! end of CG loop
      
c     to avoid using "exit"
c     if iters has reached max_iters there is no convergence
      
      IF (iters .lt. streamice_max_cg_iter) THEN
       conv_flag = 1
      ENDIF

!       DO bj = myByLo(myThid), myByHi(myThid)
!        DO bi = myBxLo(myThid), myBxHi(myThid)
!         DO j=1-OLy,sNy+OLy
!          DO i=1-OLy,sNx+OLy
!           IF (STREAMICE_umask(i,j,bi,bj).eq.3.0)
!      &     cg_Uin(i,j,bi,bj)=u_bdry_values_SI(i,j,bi,bj)
!           IF (STREAMICE_vmask(i,j,bi,bj).eq.3.0)
!      &     cg_Vin(i,j,bi,bj)=v_bdry_values_SI(i,j,bi,bj)
!          ENDDO
!         ENDDO
!        ENDDO
!       ENDDO       
! 
!       _EXCH_XY_RL( cg_Uin, myThid )
!       _EXCH_XY_RL( cg_Vin, myThid )        

#endif /* ifndef ALLOW_PETSC */

#else /* STREAMICE_SERIAL_TRISOLVE */

      iters = 0

      CALL STREAMICE_TRIDIAG_SOLVE( 
     U                               cg_Uin,     ! x-velocities
     U                               cg_Vin,
     U                               cg_Bu,      ! force in x dir
     I                               A_uu,       ! section of matrix that multiplies u and projects on u
     I                               STREAMICE_umask,
     I                               myThid )

#endif

      CALL TIMER_STOP ('STREAMICE_CG_SOLVE',myThid)


#endif
      RETURN
      END
