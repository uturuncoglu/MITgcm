C $Header: /u/gcmpack/MITgcm/pkg/ecco/cost_etan0.F,v 1.5 2012/09/04 15:04:51 gforget Exp $
C $Name: checkpoint65 $

#include "ECCO_OPTIONS.h"


      subroutine cost_etan0(
     I                           myiter,
     I                           mytime,
     I                           mythid
     &                         )

c     ==================================================================
c     SUBROUTINE cost_zonstress
c     ==================================================================
c
c     o Calculate the etan0 contribution to the cost function.
c
c     started: Christian Eckert eckert@mit.edu 30-Jun-1999
c
c     changed: Christian Eckert eckert@mit.edu 25-Feb-2000
c
c              - Restructured the code in order to create a package
c                for the MITgcmUV.
c
c     ==================================================================
c     SUBROUTINE cost_zonstress
c     ==================================================================

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"

#include "ecco_cost.h"
#include "CTRL_SIZE.h"
#include "ctrl.h"
#include "ctrl_dummy.h"
#include "optim.h"

c     == routine arguments ==

      integer myiter
      _RL     mytime
      integer mythid

c     == local variables ==

      integer bi,bj
      integer i,j,k
      integer itlo,ithi
      integer jtlo,jthi
      integer jmin,jmax
      integer imin,imax
      integer nrec
      integer irec
      integer ilfld

      _RL fctile
      _RL tmpx
      _RL lengthscale

      logical doglobalread
      logical ladinit

      character*(80) fnamefld

      character*(MAX_LEN_MBUF) msgbuf

c     == external functions ==

      integer  ilnblnk
      external ilnblnk

c     == end of interface ==

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)
      jmin = 1
      jmax = sny
      imin = 1
      imax = snx

      lengthscale = 1. _d 0

c--   Read state record from global file.
      doglobalread = .false.
      ladinit      = .false.

      irec = 1

#ifdef ALLOW_ETAN0_COST_CONTRIBUTION

      if (optimcycle .ge. 0) then
        ilfld = ilnblnk( xx_etan_file )
        write(fnamefld(1:80),'(2a,i10.10)')
     &       xx_etan_file(1:ilfld),'.',optimcycle
      endif

      call active_read_xy(
     &     fnamefld, tmpfld2d, irec, doglobalread,
     &     ladinit, optimcycle, mythid, xx_etan_dummy )

c--     Loop over this thread s tiles.
        do bj = jtlo,jthi
          do bi = itlo,ithi

c--         Determine the weights to be used.

            fctile = 0. _d 0
             do j = jmin,jmax
              do i = imin,imax
                if (_hFacC(i,j,1,bi,bj) .ne. 0.) then
                  tmpx = tmpfld2d(i,j,bi,bj)
                IF ( .NOT.ctrlSmoothCorrel2D ) THEN
                  fctile = fctile
     &                 + wetan(i,j,bi,bj)*cosphi(i,j,bi,bj)
     &                 *tmpx*tmpx
                ELSE !IF ( .NOT.ctrlSmoothCorrel2D ) THEN
                  fctile = fctile + tmpx*tmpx
                ENDIF !IF ( .NOT.ctrlSmoothCorrel2D ) THEN
                  if ( wetan(i,j,bi,bj)*cosphi(i,j,bi,bj).ne.0. )
     &                 num_etan0(bi,bj) = num_etan0(bi,bj) + 1. _d 0
                endif
             enddo
            enddo

            objf_etan0(bi,bj) = objf_etan0(bi,bj) + fctile

          enddo
        enddo

#endif

      return
      end


