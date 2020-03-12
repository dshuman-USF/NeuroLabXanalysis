c Copyright 2005-2020 Kendall F. Morris

c This file is part of the Xanalysis software suite.

c     The Xanalysis software suite is free software: you can redistribute
c     it and/or modify it under the terms of the GNU General Public
c     License as published by the Free Software Foundation, either
c     version 3 of the License, or (at your option) any later version.

c     The suite is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c     GNU General Public License for more details.

c     You should have received a copy of the GNU General Public License
c     along with the suite.  If not, see <https://www.gnu.org/licenses/>.

      module mod_calc_ach_v4
      contains
*       *****************************************************************
*       *                                                               *
*       *       Filename = calc_ach_v4.f                                *
*       *       date of last revision:   27-oct-2004    lss             *
*       *                                                               *
*       *       mar-2000        lss                                     *
*       *       modified to allow ID codes to range from 1 to 999;      *       
*       *        maximum number of units allowed still = 120;           *       
*       *        *.defs files inserted                                  *
*       *        (MAX_NUM_CODES=999; MAX_NUM_CHAN=120)                  *
*       *                                                               *
*       *       link with x2000_v2 code                                 *
*       *                                                               *
*       *       This subroutine of Xanalysis calculates and builds      *
*       *       autocorrelograms.                                       *
*       *                                                               *
*       *****************************************************************
*       A variable pointer (ITB) is used to enhance the efficiency of this subroutine by preventing
*       the comparison of every target event to each reference event.
*
*       The following is a stab at an explanation of why I have monkeyed around with the spike times to create
*       temporary arrays (REFtmp and TARtmp):

*               Double precision arithmetic is a treacherous thing.  DP floats may not be represented exactly within
*               the binary depths of the computer.  So ... it behooves one to do as much integer arithmetic as
*               possible.  Since our spike times are initially expressed as milliseconds (xxx.x) in SPIKETIMES(),
*               we can get rid of the decimal point by multiplying by 10 (*10D0) and then rounding to the nearest 
*               integer (IDNINT):  
*                       TIME = IDNINT(TIME * 10D0)
*               Now each time is expressed in tenths of a millisecond, and we can do "integer" math with it.
*       
*               OK, now for the binwidth thing and how to know which bin of the histogram to increment:
*               Our binwidths are also expressed as milliseconds (xx.x).  That's not OK to use with the IDNINT'ed
*               spike times, so we have to do BINW = IDNINT(BINW*10D0).  

*               Now we can calculate the bin by: BIN# = (DELT/BINW_rnd) + 1     [BIN# is defined as LAG in this subroutine]
*

      SUBROUTINE calcACH(SPIKETIMES,ITAL,BINW,IHIST,LAST_REF,
     +     LAST_TAR,IRPT,ITPT,ZOFST,ISTOPT,IFO,
     +     IDs)

      include 'x2000parameter.defs'

      double precision
     +     DELT,TREF,TTAR,WIND
      double precision
     +     SPIKETIMES(:,:)
      double precision, allocatable ::
     +     REFtmp(:),
     +     TARtmp(:)

      integer IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN),IHIST(NUM_BINS),
     +     LAG,R,T,ITB,LAST_REF,LAST_TAR
      real ZOFST                !offset for each REF event - currently set to 0.0

      if(.false.)print *,IFO    !suppress unused variable warning
      allocate (reftmp (size (spiketimes, 1)))
      allocate (tartmp (size (spiketimes, 1)))
      
      do i = 1, NUM_BINS
         IHIST(i)=0             !initialize the histogram array
      end do
      LAST_REF=0
      LAST_TAR=0
*
      WIND=BINW*100.0D0*10D0    !WIND = total histogram time "window" (ACHs have 100 bins)
      ITB=1                     !initial pointer value for TARGET_EVENT do loop.

      do i = 1, ITAL(IDs(IRPT))
         REFtmp(i) = IDNINT((SPIKETIMES(i,IDs(IRPT))+ZOFST)*10D0)
      end do
      do i = 1, ITAL(IDs(ITPT))
         TARtmp(i) = IDNINT((SPIKETIMES(i,IDs(ITPT)))*10D0)
      end do


*       
      REF_EVENT: DO R=1,ITAL(IDs(IRPT)) !get a reference event
      TREF=REFtmp(R)            !apply the offset to the ref event time
      LAST_REF=LAST_REF+1       !LAST_REF = a counter of # of ref events used in this ACH
c           TARGET_EVENT: DO T=ITB,ITAL(IDs(ITPT))      !get a target event
      TARGET_EVENT: DO T = R+1, ITAL(IDs(ITPT)) !get a target event
      TTAR=TARtmp(T)
c             if(TTAR.lt.TREF)cycle TARGET_EVENT
      DELT = TTAR - TREF        !DELT=time between ref and target events (msec.)
*                                                        (only want target events that occur after 
*                                                        reference events for an ACH)
*                                                        so ... get next target event
      IF(DELT.GT.WIND) CYCLE REF_EVENT !DELT>WIND --> target occurs too much later
*                                                        than reference event so...get next reference event
*
      IF(DELT.eq.0.0)cycle TARGET_EVENT
c             IF((DELT.eq.0.0).AND.(IRPT.EQ.ITPT)) then !ref and target events are same - don't include in ACH and...
c                ITB=T                                  !update the pointer value for TARGET_EVENT do loop and...
c                CYCLE TARGET_EVENT                     !get next target event
c             end if

*       ***** at this point, have a valid target event that occurs within the histogram window ***

      LAST_TAR=LAST_TAR+1       !LAST_TAR = a counter of # of target events used

      IF((ISTOPT.NE.0).AND.(LAST_REF.GT.ISTOPT))EXIT REF_EVENT !have used up allotted # stimuli so this ACH is 
!completed; leave subroutine (ISTOPT=0-->use all stimuli)
      LAG = DELT/IDNINT(BINW * 10D0) + 1 !DON'T FORGET THAT THIS WILL TRUNCATE BECAUSE ASSIGNING TO INT VALUE!!!
      IHIST(LAG)=IHIST(LAG)+1   !increment that location in IHIST()
c             IF(IFO.NE.0)CYCLE REF_EVENT               ! currently first order ACH not used
      end do TARGET_EVENT       !get another target event
      end do REF_EVENT

      IHIST(NUM_BINS)=0
      RETURN
      END

      end module mod_calc_ach_v4
