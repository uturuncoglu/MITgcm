C $Header: /u/gcmpack/MITgcm/pkg/mom_common/MOM_COMMON_OPTIONS.h,v 1.4 2013/07/30 19:05:54 jmc Exp $
C $Name: checkpoint65 $

C CPP options file for mom_common package
C Use this file for selecting CPP options within the mom_common package

#ifndef MOM_COMMON_OPTIONS_H
#define MOM_COMMON_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_MOM_COMMON
C     Package-specific options go here

C allow isotropic 3-D Smagorinsky viscosity
#undef ALLOW_SMAG_3D

C allow full 3D specification of horizontal Laplacian Viscosity
#undef ALLOW_3D_VISCAH

C allow full 3D specification of horizontal Biharmonic Viscosity
#undef ALLOW_3D_VISCA4

#endif /* ALLOW_MOM_COMMON */
#endif /* MOM_COMMON_OPTIONS_H */
