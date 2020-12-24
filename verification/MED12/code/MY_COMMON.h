C $Header: /u/gcmpack/MITgcm/model/inc/SURFACE.h,v 1.20 2014/04/29 21:03:00 jmc Exp $
C $Name: checkpoint65 $
C
CBOP
C    !ROUTINE: MY_COMMON.h
C    !INTERFACE:
C    include MY_COMMON.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | SURFACE.h
C     | o Header file defining surface-related model variables
C     *==========================================================*
C     | Contains variables relative to the surface position
C     | that are held fixed in linear free-surface formulation
C     | but can vary with time with a non-linear free-surface.
C     *==========================================================*
C     \ev
CEOP


cgmMASK(
C     my_2D_mask    :: 2D Mask for heat and salt fluxes
C     SurRelMsk_T   :: Surface Relaxation Mask x Temperature
C     SurRelMsk_S   :: Surface Relaxation Mask x Salinity
C
C     EmPmR_msk     :: 2D Mask where EmPmR modify the salinity
C     EmPmR_SALT_msk:: 2D Mask where is specified the minimum salinity allowed
cgmMASK)

cgmMASK(
      COMMON /SURF_MASK/ SurRelMsk_T,SurRelMsk_S,my_2D_mask,
     &EmPmR_msk,EmPmR_SALT_msk
      _RL  SurRelMsk_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  SurRelMsk_S(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL   my_2D_mask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL        EmPmR_msk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL   EmPmR_SALT_msk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
cgmMASK)
