C $Header: /u/gcmpack/MITgcm/pkg/ecco/cost_gencost_all.F,v 1.13 2014/06/27 14:23:19 gforget Exp $
C $Name: checkpoint65 $

#include "ECCO_OPTIONS.h"


      subroutine cost_gencost_all( myiter, mytime, mythid )

c     ==================================================================
c     SUBROUTINE cost_gencost_all
c     ==================================================================
c
c     o Evaluate user defined penalty terms
c
c     ==================================================================
c     SUBROUTINE cost_gencost_all
c     ==================================================================

      implicit none

c     == global variables ==

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"

#include "cal.h"
#include "CTRL_SIZE.h"
#include "ctrl.h"
#include "ctrl_dummy.h"
#include "ecco_cost.h"

c     == routine arguments ==

      integer mythid
      integer myiter
      _RL     mytime

#ifdef ALLOW_GENCOST_CONTRIBUTION
c     == local variables ==

      integer i,j,k
      character*(max_len_mbuf) msgbuf

      integer nnzbar
      integer nnzobs
      integer nrecloc

c     == external functions ==

c     == end of interface ==

      do k = 1, NGENCOST

         if ( (using_gencost(k)).AND.( gencost_flag(k).EQ.1 ) ) then

         write(msgbuf,'(a,x,i2.2)') 'ph-cost gencost #', k
         call print_message( msgbuf, standardmessageunit,
     &        SQUEEZE_RIGHT , mythid)
         nnzbar = 1
         nnzobs = 1
         call cost_generic(
     &        nnzbar, gencost_barfile(k),
     &        gencost_barfld(1-Olx,1-Oly,1,1,k),
     &        xx_genbar_dummy(k),
     &        nnzobs, gencost_datafile(k), mult_gencost(k),
     &        gencost_nrec(k), modelstartdate, gencost_period(k),
     &        gencost_mask(k), gencost_weight(1-Olx,1-Oly,1,1,k),
     &        gencost_spmin(k), gencost_spmax(k),
     &        gencost_spzero(k),
     &        objf_gencost(1,1,k), num_gencost(1,1,k),
     &        myiter, mytime, mythid )

         endif

      end do

      call cost_gencost_boxmean( mythid )

# if (defined(ALLOW_GENCOST_SSHV4))
      IF (using_cost_altim) 
     &    call cost_gencost_sshv4( myiter, mytime, mythid )
# endif /* ALLOW_GENCOST_SSHV4 */

# if (defined(ALLOW_GENCOST_SSTV4))
      IF (using_cost_sst)
     &    call cost_gencost_sstv4( myiter, mytime, mythid )
# endif /* ALLOW_GENCOST_SSTV4 */

#ifdef ALLOW_GENCOST_SEAICEV4
      IF (using_cost_seaice)
     &    call cost_gencost_seaicev4( mythid )
# endif /* ALLOW_GENCOST_SEAICEV4 */

#endif /* ALLOW_GENCOST_CONTRIBUTION */

      end
