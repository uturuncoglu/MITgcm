C $Header: /u/gcmpack/MITgcm/pkg/cheapaml/cheapaml_diagnostics_init.F,v 1.7 2013/06/11 01:42:50 jmc Exp $
C $Name: checkpoint65 $

#include "CHEAPAML_OPTIONS.h"
CBOP
C     !ROUTINE: CHEAPAML_DIAGNOSTICS_INIT
C     !INTERFACE:
      SUBROUTINE CHEAPAML_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE CHEAPAML_DIAGNOSTICS_INIT
C     | o Routine to initialize list of all available diagnostics
C     |   for CHEAPAML package
C     *==========================================================*
C     \ev
C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "CHEAPAML.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     === Local variables ===
C     diagNum   :: diagnostics number in the (long) list of available diag.
C     diagMate  :: diag. mate number in the (long) list of available diag.
C     diagName  :: local short name (8c) of a diagnostics
C     diagCode  :: local parser field with characteristics of the diagnostics
C              cf head of S/R DIAGNOSTICS_INIT_EARLY or DIAGNOSTICS_MAIN_INIT
C     diagUnits :: local string (16c): physical units of a diagnostic field
C     diagTitle :: local string (80c): description of field in diagnostic
      INTEGER       diagNum
      INTEGER       diagMate
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


      diagName  = 'CH_TAIR '
      diagTitle = 'cheapAML Surface Air Temp'
      diagUnits = 'degres C        '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I            diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CH_QNET '
      diagTitle = 'cheapAML Surface heating '
      diagUnits = 'W/m2            '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I            diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CH_QAIR '
      diagTitle = 'cheapAML Specific Humidity'
      diagUnits = 'kg/kg           '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I            diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CH_EmP  '
      diagTitle = 'cheapAML Evap-Precip     '
      diagUnits = 'kg/m2/s         '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I            diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CH_Trace'
      diagTitle = 'cheapAML Passive Tracer  '
      diagUnits = 'Depend on Tracer '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I            diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CH_Uwind'
      diagTitle = 'cheapAML Zonal Wind      '
      diagUnits = 'm/s           '
      diagCode  = 'UU      L1      '
      diagmate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I      diagName, diagCode, diagUnits, diagTitle, diagmate, myThid )

      diagName  = 'CH_Vwind'
      diagTitle = 'cheapAML Meridional Wind '
      diagUnits = 'm/s             '
      diagCode  = 'VV      L1      '
      diagmate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I      diagName, diagCode, diagUnits, diagTitle, diagmate, myThid )

      diagName  = 'CH_LH   '
      diagTitle = 'cheapAML latent heat     '
      diagUnits = 'W/m2            '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I            diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CH_SH   '
      diagTitle = 'cheapAML sensible heat   '
      diagUnits = 'W/m2            '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I            diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CH_Prec '
      diagTitle = 'cheapAML precipitation   '
      diagUnits = 'kg/m2/s         '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I            diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CH_q100 '
      diagTitle = 'cheapAML upper Specific Humidity'
      diagUnits = 'kg/kg           '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I            diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CH_ssqt '
      diagTitle = 'cheapAML surf sat Specific Humid.'
      diagUnits = 'kg/kg           '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I            diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
