C $Header: /u/gcmpack/MITgcm/pkg/ecco/ecco_readparms.F,v 1.70 2014/06/27 14:23:19 gforget Exp $
C $Name: checkpoint65 $

#include "ECCO_OPTIONS.h"

      subroutine ecco_readparms( mythid )

c     ==================================================================
c     SUBROUTINE ecco_readparms
c     ==================================================================
c
c     o This routine initialises the package cost.
c     started: Ralf Giering 18-Jan-2001
c
c     ==================================================================
c     SUBROUTINE ecco_readparms
c     ==================================================================

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "PARAMS.h"

#include "ecco_cost.h"

c     == routine arguments ==

      integer mythid

c     == local variables ==
C     msgBuf      - Informational/error meesage buffer
C     iUnit       - Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER k, iUnit, num_file, num_var
#ifdef ALLOW_GENCOST_CONTRIBUTION
      INTEGER IL
#endif

Catn-- retired parameters
      character*(MAX_LEN_FNAM) topexmeanfile
      integer nRetired

c     == external functions ==

      integer  ilnblnk
      external ilnblnk

c     == end of interface ==

c--   Read the namelist input.
      namelist /ecco_cost_nml/
     &                   data_errfile,
     &                   tbarfile, sbarfile,sstbarfile,
     &                   psbarfile, bpbarfile,iestaubarfile,
     &                   ubarfile, vbarfile, wbarfile,
     &                   tauxbarfile, tauybarfile,
     &                   hfluxmeanbarfile, sfluxmeanbarfile,
     &                   costTranspDataFile, mdtdatfile,
     &                   mdtstartdate1, mdtstartdate2,
     &                   mdtenddate1, mdtenddate2,
     &                   using_cost_altim, 
     &                   using_cost_bp, using_cost_scat, 
     &                   using_cost_sst, using_cost_seaice,
     &                   using_topex, using_ers, using_gfo,
     &                   topexstartdate1, topexstartdate2,
     &                   topexperiod, topexfile, tpTimeMaskFile,
     &                   topexintercept, topexslope, ersTimeMaskFile,
     &                   ersstartdate1, ersstartdate2, gfoTimeMaskFile,
     &                   ersperiod, ersfile,
     &                   ersintercept, ersslope,
     &                   gfostartdate1, gfostartdate2,
     &                   gfoperiod, gfofile,
     &                   gfointercept, gfoslope,
     &                   scatstartdate1, scatstartdate2,
     &                   scatperiod, scatxdatfile, scatydatfile,
     &                   scatx_errfile, scaty_errfile,
     &                   ssh_errfile, sshv4cost_scalefile,
     &                   sshv4cost_errfile, sshv4cost_errfactor,
     &                   tp_errfile, ers_errfile, gfo_errfile,
     &                   tmistartdate1, tmistartdate2, tmidatfile,
     &                   sststartdate1, sststartdate2, sstdatfile,
     &                   sssstartdate1, sssstartdate2, sssdatfile,
     &                   bpstartdate1, bpstartdate2, bpdatfile,
     &                   iesstartdate1, iesstartdate2, iesdatfile,
     &                   tdatfile, sdatfile, ctdtfile, ctdsfile,
     &                   ctdtclimfile, ctdsclimfile, xbtfile,
     &                   argotstartdate1, argotstartdate2, argotfile,
     &                   argosstartdate1, argosstartdate2, argosfile,
     &                   udriftfile, vdriftfile,
     &                   udrifterrfile, vdrifterrfile,
     &                   curmtrufile,curmtrvfile,
     &                   salterrfile, temperrfile, velerrfile,
     &                   salt0errfile, temp0errfile, vel0errfile,
     &                   etan0errfile, uvel0errfile, vvel0errfile,
     &                   ssterrfile, ssserrfile, bperrfile,ieserrfile,
     &                   hflux_errfile, sflux_errfile,
     &                   tauu_errfile, tauv_errfile,
     &                   atemp_errfile, aqh_errfile,
     &                   precip_errfile, swflux_errfile, swdown_errfile,
     &                   lwflux_errfile, lwdown_errfile, evap_errfile,
     &                   snowprecip_errfile, apressure_errfile,
     &                   uwind_errfile, vwind_errfile,
     &                   geoid_errfile, geoid_covariancefile,
     &                   bottomdrag_errfile,edtau_errfile,
     &                   kapgm_errfile,kapredi_errfile,diffkr_errfile,
     &                   usercost_datafile, usercost_errfile,
     &                   mult_hflux, mult_sflux, mult_tauu, mult_tauv,
     &                   mult_hfluxmm, mult_sfluxmm,
     &                   mult_atemp, mult_aqh,
     &                   mult_precip, mult_swflux, mult_swdown,
     &                   mult_lwflux, mult_lwdown, mult_evap,
     &                   mult_snowprecip, mult_apressure,
     &                   mult_uwind, mult_vwind,
     &                   mult_hmean, mult_h, mult_sshv4cost,
     &                   mult_tp, mult_ers, mult_gfo,
     &                   mult_temp, mult_salt,
     &                   mult_tmi, mult_sst, mult_sss, mult_bp,
     &                   mult_ctdt, mult_ctds, mult_ies,
     &                   mult_ctdtclim, mult_ctdsclim,
     &                   mult_xbt, mult_drift,
     &                   mult_argot, mult_argos,
     &                   mult_scatx, mult_scaty,
     &                   mult_sdrift, mult_tdrift, mult_wdrift,
     &                   mult_temp0, mult_salt0,
     &                   mult_etan0, mult_uvel0, mult_vvel0,
     &                   mult_obcsn, mult_obcss,
     &                   mult_obcsw, mult_obcse, mult_obcsvol,
     &                   mult_curmtr, mult_ageos,
     &                   mult_kapgm, mult_kapredi, mult_diffkr,
     &                   mult_edtau, mult_bottomdrag,
     &                   mult_usercost, mult_transp,
     &                   mult_smooth_ic, mult_smooth_bc,
     &                   whflux0, wsflux0, wtau0,
     &                   wbottomdrag0,
     &                   watemp0, waqh0,
     &                   wprecip0, wswflux0, wswdown0, wwind0,
     &                   wsnowprecip0, wlwflux0, wlwdown0, wevap0,
     &                   wapressure0, wdiffkr0, wkapgm0, wedtau0,
     &                   wkapredi0, wmean_hflux, wmean_sflux, wmean_tau,
     &                   wmean_atemp, wmean_aqh, wmean_precip,
     &                   wmean_swflux, wmean_swdown, wmean_wind,
     &                   wmean_lwflux, wmean_lwdown, wmean_evap,
     &                   wmean_snowprecip, wmean_apressure,
     &                   cost_iprec, cost_yftype
     &                   ,topexmeanfile

#ifdef ALLOW_GENCOST_CONTRIBUTION
      namelist /ecco_gencost_nml/
     &         using_gencost, 
     &         gencost_barfile,
     &         gencost_datafile,
     &         gencost_name,
     &         gencost_scalefile,
     &         gencost_errfile,
     &         gencost_mask,
     &         gencost_spmin,
     &         gencost_spmax,
     &         gencost_spzero,
     &         gencost_avgperiod,
#ifdef ALLOW_SMOOTH
     &         gencost_smooth2Ddiffnbt,
#endif /* ALLOW_SMOOTH */
#ifdef ALLOW_GENCOST_TIMEVARY_WEIGHT
     &         gencost_timevaryweight,
#endif /* ALLOW_GENCOST_TIMEVARY_WEIGHT */
     &         mult_gencost
#endif

      _BEGIN_MASTER( myThid )

c--     Set default values.

c run time switches for main cost terms
#ifdef ALLOW_SSH_COST_CONTRIBUTION
      using_cost_altim=.TRUE.
#else
      using_cost_altim=.FALSE.
#endif
#ifdef ALLOW_SEAICE_COST_CONTRIBUTION
      using_cost_seaice=.TRUE.
#else
      using_cost_seaice=.FALSE.
#endif
#ifdef ALLOW_BP_COST_CONTRIBUTION
      using_cost_bp=.TRUE.
#else
      using_cost_bp=.FALSE.
#endif
#ifdef ALLOW_SST_COST_CONTRIBUTION
      using_cost_sst=.TRUE.
#else
      using_cost_sst=.FALSE.
#endif
#if (defined (ALLOW_SCAT_COST_CONTRIBUTION) || \
     defined (ALLOW_DAILYSCAT_COST_CONTRIBUTION) )
      using_cost_scat=.TRUE.
#else
      using_cost_scat=.FALSE.
#endif

c other run time parameters
        data_errfile          =    ' '
        tbarfile              =    'tbar'
        sbarfile              =    'sbar'
        sstbarfile            =    'sstbar'
        psbarfile             =    'psbar'
        bpbarfile             =    'bpbar'
        iestaubarfile         =    'iestaubar'
        ubarfile              =    'ubar'
        vbarfile              =    'vbar'
        wbarfile              =    'wbar'
        tauxbarfile           =    'tauxbar'
        tauybarfile           =    'tauybar'
        hfluxmeanbarfile      =    'hfluxmeanbar'
        sfluxmeanbarfile      =    'sfluxmeanbar'
        costTranspDataFile    =    ' '
#ifdef ALLOW_SSH_TPANOM_COST_CONTRIBUTION
        using_topex           = .true.
#else
        using_topex           = .false.
#endif
        topexstartdate1       =      0
        topexstartdate2       =      0
        topexperiod           =      0. _d 0
        topexintercept        =      0. _d 0
        topexslope            =      0. _d 0
        topexfile             =    ' '
catn: defaulting mdt[start,end]date1 to pavlis so that
catn  eccov4 will run without needing to read in these params
        mdtdatfile            =    ' '
        mdtstartdate1         =      19930101
        mdtstartdate2         =      0
        mdtenddate1           =      20041231
        mdtenddate2           =      0
#ifdef ALLOW_SSH_ERSANOM_COST_CONTRIBUTION
        using_ers             = .true.
#else
        using_ers             = .false.
#endif
        ersstartdate1         =      0
        ersstartdate2         =      0
        ersperiod             =      0. _d 0
        ersintercept          =      0. _d 0
        ersslope              =      0. _d 0
        ersfile               =    ' '
#ifdef ALLOW_SSH_GFOANOM_COST_CONTRIBUTION
        using_gfo             = .true.
#else
        using_gfo             = .false.
#endif
        gfostartdate1         =      0
        gfostartdate2         =      0
        gfoperiod             =      0. _d 0
        gfointercept          =      0. _d 0
        gfoslope              =      0. _d 0
        gfofile               =    ' '
        scatstartdate1        =      0
        scatstartdate2        =      0
        scatperiod            =      0. _d 0
        scatxdatfile          =    ' '
        scatydatfile          =    ' '
        ssh_errfile           =    ' '
        tp_errfile            =    ' '
        ers_errfile           =    ' '
        gfo_errfile           =    ' '
        scatx_errfile         =    ' '
        scaty_errfile         =    ' '
        tmistartdate1         =      0
        tmistartdate2         =      0
        tmidatfile            =    ' '
        sststartdate1         =      0
        sststartdate2         =      0
        sstdatfile            =    ' '
        sssstartdate1         =      0
        sssstartdate2         =      0
        sssdatfile            =    ' '
        bpstartdate1          =      0
        bpstartdate2          =      0
        bpdatfile             =    ' '
        iesstartdate1         =      0
        iesstartdate2         =      0
        iesdatfile            =    ' '
        tdatfile              =    ' '
        sdatfile              =    ' '
        ctdtfile              =    ' '
        ctdsfile              =    ' '
        ctdtclimfile          =    ' '
        ctdsclimfile          =    ' '
        curmtrufile           =    ' '
        curmtrvfile           =    ' '
        xbtfile               =    ' '
        argotstartdate1       =      0
        argotstartdate2       =      0
        argotfile             =    ' '
        argosstartdate1       =      0
        argosstartdate2       =      0
        argosfile             =    ' '
        udriftfile            =    ' '
        vdriftfile            =    ' '
        udrifterrfile         =    ' '
        vdrifterrfile         =    ' '
        salterrfile           =    ' '
        temperrfile           =    ' '
        velerrfile            =    ' '
        salt0errfile          =    ' '
        temp0errfile          =    ' '
        etan0errfile          =    ' '
        uvel0errfile          =    ' '
        vvel0errfile          =    ' '
        vel0errfile           =    ' '
        ssterrfile            =    ' '
        ssserrfile            =    ' '
        bperrfile             =    ' '
        ieserrfile            =    ' '
        hflux_errfile         =    ' '
        sflux_errfile         =    ' '
        tauu_errfile          =    ' '
        tauv_errfile          =    ' '
        atemp_errfile         =    ' '
        aqh_errfile           =    ' '
        precip_errfile        =    ' '
        swflux_errfile        =    ' '
        swdown_errfile        =    ' '
        lwflux_errfile        =    ' '
        lwdown_errfile        =    ' '
        evap_errfile          =    ' '
        snowprecip_errfile    =    ' '
        apressure_errfile     =    ' '
        uwind_errfile         =    ' '
        vwind_errfile         =    ' '
        geoid_errfile         =    ' '
        geoid_covariancefile  =    ' '
        bottomdrag_errfile    =    ' '
        edtau_errfile         =    ' '
        kapgm_errfile         =    ' '
        kapredi_errfile       =    ' '
        diffkr_errfile        =    ' '
        do k=1,NUSERCOST
         usercost_datafile(k)  =    ' '
         usercost_errfile(k)   =    ' '
        enddo
        do k=1,NSSHV4COST
         sshv4cost_scalefile(k) =    ' '
         sshv4cost_errfile(k)   =    ' '
         sshv4cost_errfactor(k) =     1. _d 0
        enddo
        tpTimeMaskFile        =    ' '
        ersTimeMaskFile       =    ' '
        gfoTimeMaskFile       =    ' '

        mult_hflux            =      0. _d 0
        mult_sflux            =      0. _d 0
        mult_hfluxmm          =      0. _d 0
        mult_sfluxmm          =      0. _d 0
        mult_tauu             =      0. _d 0
        mult_tauv             =      0. _d 0
        mult_atemp            =      0. _d 0
        mult_aqh              =      0. _d 0
        mult_precip           =      0. _d 0
        mult_swflux           =      0. _d 0
        mult_swdown           =      0. _d 0
        mult_lwflux           =      0. _d 0
        mult_lwdown           =      0. _d 0
        mult_evap             =      0. _d 0
        mult_snowprecip       =      0. _d 0
        mult_apressure        =      0. _d 0
        mult_uwind            =      0. _d 0
        mult_vwind            =      0. _d 0
        mult_hmean            =      0. _d 0
        mult_h                =      1. _d 0
        mult_tp               =      0. _d 0
        mult_ers              =      0. _d 0
        mult_gfo              =      0. _d 0
        mult_temp             =      0. _d 0
        mult_salt             =      0. _d 0
        mult_temp0            =      0. _d 0
        mult_salt0            =      0. _d 0
        mult_etan0            =      0. _d 0
        mult_uvel0            =      0. _d 0
        mult_vvel0            =      0. _d 0
        mult_tmi              =      0. _d 0
        mult_sst              =      0. _d 0
        mult_sss              =      0. _d 0
        mult_bp               =      0. _d 0
        mult_ies              =      0. _d 0
        mult_ctdt             =      0. _d 0
        mult_ctds             =      0. _d 0
        mult_ctdtclim         =      0. _d 0
        mult_ctdsclim         =      0. _d 0
        mult_xbt              =      0. _d 0
        mult_argot            =      0. _d 0
        mult_argos            =      0. _d 0
        mult_drift            =      0. _d 0
        mult_tdrift           =      0. _d 0
        mult_sdrift           =      0. _d 0
        mult_wdrift           =      0. _d 0
        mult_scatx            =      0. _d 0
        mult_scaty            =      0. _d 0
        mult_obcsn            =      0. _d 0
        mult_obcss            =      0. _d 0
        mult_obcsw            =      0. _d 0
        mult_obcse            =      0. _d 0
        mult_obcsvol          =      0. _d 0
        mult_ageos            =      0. _d 0
        mult_curmtr           =      0. _d 0
        mult_kapgm            =      0. _d 0
        mult_kapredi            =      0. _d 0
        mult_diffkr           =      0. _d 0
        mult_edtau            =      0. _d 0
        mult_bottomdrag       =      0. _d 0
        do k=1,NUSERCOST
        mult_usercost(k)       =      0. _d 0
        enddo
        do k=1,NSSHV4COST
        mult_sshv4cost(k)       =      0. _d 0
        enddo
        mult_smooth_ic        =      0. _d 0
        mult_smooth_bc        =      0. _d 0
        mult_transp           =      0. _d 0
        watemp0               =      0. _d 0
        waqh0                 =      0. _d 0
        wprecip0              =      0. _d 0
        wswflux0              =      0. _d 0
        wswdown0              =      0. _d 0
        wlwflux0              =      0. _d 0
        wlwdown0              =      0. _d 0
        wevap0                =      0. _d 0
        wsnowprecip0          =      0. _d 0
        wapressure0           =      0. _d 0
        wwind0                =      0. _d 0
        wdiffkr0              =      1. _d 0
        wkapgm0               =      1. _d 0
        wkapredi0               =      1. _d 0
        wedtau0               =      1. _d 0
        whflux0               =      0. _d 0
        wsflux0               =      0. _d 0
        wtau0                 =      0. _d 0
        wbottomdrag0          =      1. _d 0
        wmean_atemp           =      0. _d 0
        wmean_aqh             =      0. _d 0
        wmean_precip          =      0. _d 0
        wmean_swflux          =      0. _d 0
        wmean_swdown          =      0. _d 0
        wmean_lwflux          =      0. _d 0
        wmean_lwdown          =      0. _d 0
        wmean_evap            =      0. _d 0
        wmean_snowprecip      =      0. _d 0
        wmean_apressure       =      0. _d 0
        wmean_wind            =      0. _d 0
        wmean_hflux           =      0. _d 0
        wmean_sflux           =      0. _d 0
        wmean_tau             =      0. _d 0
        cost_iprec  = precFloat32
        cost_yftype = 'RL'

Catn-- retired parameters
        topexmeanfile         =    ' '
        nRetired              =      0

c--     Next, read the data.ecco file.
        WRITE(msgBuf,'(A)') 'ECCO_READPARMS: opening data.ecco'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)

        CALL OPEN_COPY_DATA_FILE(
     I                          'data.ecco', 'ECCO_READPARMS',
     O                          iUnit,
     I                          myThid )

        READ(  iUnit, nml = ecco_cost_nml )

        WRITE(msgBuf,'(A)')
     &       'ECCO_READPARMS: finished reading #1: ecco_cost_nml'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &       SQUEEZE_RIGHT , 1)

C--   Check parameters:
        IF ( cost_yftype.NE.'RL' ) THEN
          STOP 'ABNORMAL END: S/R ECCO_READPARMS: invalid cost_yftype'
        ENDIF
Catn -- retired parameters
      IF ( topexmeanfile .NE. ' ' ) THEN
       nRetired = nRetired + 1
       WRITE(msgBuf,'(A,A)')
     &  'S/R ECCO_READPARMS: "topexmeanfile" ',
     &  'is no longer allowed in file "data.ecco"'
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(A)')
     &  'S/R ECCO_READPARMS: set "mdtdatfile" instead'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R ECCO_READPARMS'
      ENDIF
      IF ( nRetired .GT. 0 ) THEN
       WRITE(msgBuf,'(2A)') 'S/R ECCO_READPARMS: ',
     &  'Error reading parameter file "data.ecco"'
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(A)')
     &  'some out of date parameters were found in the namelist'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R ECCO_READPARMS'
      ENDIF

      _END_MASTER( myThid )

#ifdef ALLOW_GENCOST_CONTRIBUTION

      _BEGIN_MASTER( myThid )

c--   Set default values.
      do k=1,NGENCOST
         using_gencost(k)      = .FALSE.
         gencost_flag(k)       = 0
         gencost_avgperiod(k)  = '     '
         gencost_datafile(k)   = ' '
         gencost_name(k)       = 'gencost'
         gencost_scalefile(k)  = ' '
         gencost_errfile(k)    = ' '
         gencost_mask(k)       = 'c'
cgf         write(gencost_barfile(k)(1:8),'(a,i2.2)')
cgf     &        'genbar', k
         gencost_barfile(k)    = ' '
         gencost_spmin(k)      = 0. _d 0
         gencost_spmax(k)      = 0. _d 0
         gencost_spzero(k)     = 9876. _d 0
         mult_gencost(k)       = 0. _d 0
#ifdef ALLOW_SMOOTH
         gencost_smooth2Ddiffnbt(k) = 300
#endif /* ALLOW_SMOOTH */
#ifdef ALLOW_GENCOST_TIMEVARY_WEIGHT
         gencost_timevaryweight(k) = .FALSE.
#endif /* ALLOW_GENCOST_TIMEVARY_WEIGHT */
      enddo

      READ(  iUnit, nml = ecco_gencost_nml )

      WRITE(msgBuf,'(A)')
     &     'ECCO_READPARMS: finished reading #2: ecco_gencost_nml'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &     SQUEEZE_RIGHT , 1)

      do k=1,NGENCOST
         if ( (gencost_name(k).EQ.'sshv4-mdt').OR.
     &        (gencost_name(k).EQ.'sshv4-tp').OR.
     &        (gencost_name(k).EQ.'sshv4-ers').OR.
     &        (gencost_name(k).EQ.'sshv4-gfo').OR.
     &        (gencost_name(k).EQ.'sshv4-lsc').OR.
     &        (gencost_name(k).EQ.'sshv4-gmsl').OR.
     &        (gencost_name(k).EQ.'boxmean').OR.
     &        (gencost_name(k).EQ.'sstv4-amsre').OR.
     &        (gencost_name(k).EQ.'sstv4-amsre-lsc') ) then
c pre-specified generic cost terms using otherwise created barfiles
           gencost_flag(k)=-1
           using_gencost(k)=.TRUE.
#ifdef ALLOW_GENCOST_SEAICEV4
         elseif ( 
     &        (gencost_name(k).EQ.'siv4-conc').OR.
     &        (gencost_name(k).EQ.'siv4-deconc').OR.
     &        (gencost_name(k).EQ.'siv4-exconc') ) then
c pre-specified generic cost terms that create their own barfile
           gencost_flag(k)=2
           using_gencost(k)=.TRUE.
#endif
         elseif ( gencost_datafile(k) .NE. ' ' ) then
c fully interactive generic cost terms
           gencost_flag(k)=1
           using_gencost(k)=.TRUE.
         endif
      enddo

C--   Check parameters:
      do k=1,NGENCOST

        if ( (using_gencost(k)).AND.( gencost_flag(k).GE.1 ) ) then

          if ( gencost_name(k) .EQ. 'gencost' )
     &        gencost_name(k) = gencost_datafile(k)

          if ( gencost_avgperiod(k) .NE. 'day  ' .AND.
     &        gencost_avgperiod(k) .NE. 'DAY  ' .AND.
     &        gencost_avgperiod(k) .NE. 'month' .AND.
     &        gencost_avgperiod(k) .NE. 'MONTH' .AND.
     &        gencost_avgperiod(k) .NE. 'year ' .AND.
     &        gencost_avgperiod(k) .NE. 'YEAR ' ) then

            il=ilnblnk(gencost_name(k))
            WRITE(msgBuf,'(A,I2,2A)')
     &          'ERROR in ECCO_READPARMS: for gencost',k,
     &          '  ',gencost_name(k)(1:il)
            CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &          SQUEEZE_RIGHT , 1)
            STOP
     &           'ECCO_READPARMS: gencost_avgperiod not properly set'
          end if

          if ( gencost_spmin(k) .EQ. 0. _d 0 .AND.
     &        gencost_spmax(k) .EQ. 0. _d 0 ) then
            il=ilnblnk(gencost_name(k))
            WRITE(msgBuf,'(A,I2,2A)')
     &          'ERROR in ECCO_READPARMS: for gencost',k,
     &          '  ',gencost_name(k)(1:il)
            CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &          SQUEEZE_RIGHT , 1)
            STOP
     &       'ECCO_READPARMS: gencost_spmin, gencost_spmax not set'
          endif

          if ( gencost_spzero(k) .EQ. 9876. _d 0 ) then
            il=ilnblnk(gencost_name(k))
            WRITE(msgBuf,'(A,I2,2A)')
     &          'ERROR in ECCO_READPARMS: for gencost',k,
     &          '  ',gencost_name(k)(1:il)
            CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &          SQUEEZE_RIGHT , 1)
            STOP
     &       'ECCO_READPARMS: gencost_spzero not set'
          endif

        endif

      end do

      _END_MASTER( myThid )

#endif /* ALLOW_GENCOST_CONTRIBUTION */

      _BEGIN_MASTER( myThid )

      CLOSE( iUnit )

      _END_MASTER( myThid )

      WRITE(msgBuf,'(A)')
     &     'ECCO_READPARMS: done'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &     SQUEEZE_RIGHT , 1)

      _BARRIER

      end
