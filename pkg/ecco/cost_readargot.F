C $Header: /u/gcmpack/MITgcm/pkg/ecco/cost_readargot.F,v 1.3 2012/08/10 19:45:26 jmc Exp $
C $Name: checkpoint65 $

#include "ECCO_OPTIONS.h"


      subroutine cost_Readargot(
     I                               irec,
     I                               mythid
     &                             )

c     ==================================================================
c     SUBROUTINE cost_Readargot
c     ==================================================================
c
c     o Read a given record of the
c
c     started:Armin Koehl akoehl@ucsd.edu
c              - Restructured the code in order to create a package
c                for the MITgcmUV.
c
c     ==================================================================
c     SUBROUTINE cost_Readargot
c     ==================================================================

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"

#include "cal.h"
#include "ecco_cost.h"

c     == routine arguments ==

      integer irec
      integer mythid

c     == local variables ==

      integer bi,bj
      integer i,j,k
      integer itlo,ithi
      integer jtlo,jthi
      integer jmin,jmax
      integer imin,imax
      integer nobs
      integer argotrec
      integer beginargot
      integer beginrun

      _RL spval
      _RL vartile

c     == end of interface ==

      parameter (spval = -9999. )
ce    --> there is certainly a better place for this.

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)
      jmin = 1
      jmax = sny
      imin = 1
      imax = snx

      beginargot = argotstartdate(1)/10000
      beginrun = modelstartdate(1)/10000
      if ( beginargot .eq. beginrun ) then
        argotrec = mod(modelstartdate(1)/100,100) -
     &           mod(argotstartdate(1)/100,100) + irec
      else
        argotrec = ( beginrun - beginargot - 1)*nmonthyear +
     &           (nmonthyear - mod(argotstartdate(1)/100,100) +
     &            1) + mod(modelstartdate(1)/100,100) - 1 + irec
      endif
      if (argotrec.gt.0) then
      call mdsreadfield( argotfile, cost_iprec, cost_yftype, nr,
     &                   argotobs,argotrec, mythid )
      else
         do bj = jtlo,jthi
            do bi = itlo,ithi
               do k = 1,nr
               do j = jmin,jmax
                  do i = imin,imax
                     argotobs(i,j,k,bi,bj)=spval
                  enddo
               enddo
               enddo
            enddo
         enddo
      endif
      return
      end

