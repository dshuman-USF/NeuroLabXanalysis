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

      module mod_intvar4_v2a
      contains
*       lss     6-apr-2004      version = intvar4_v2.f
*
*       link with x2000_v2
*
*       increased max # spike trains and # events
*
*       ZMAX = max allowed interval (in msec fp)
*       FOR COEFFICIENT OF VARIANCE FOR CUMSUM SIG. CALC
*
*       ****************************************************************
*       ****************************************************************
*
      subroutine INTVAR(SPIKETIMES,NCELL,ICELPT,ZMAX,COEFVR,
     +     IDs)
*
*
      include 'x2000parameter.defs'
      PARAMETER (ICT=100)

      double precision SPIKETIMES(:,:)
*
      integer IDs(MAX_NUM_CODES)
      double precision TRV,SUMTRV
*
*
      DIMENSION RQ(ICT)
*
*
*
      RQ = 0.0                  !initialize array
*
*
*
*       ***** CONVERT SPIKE TIMES INTERVALS *****
*
      NINT=0
      SUMTRV=0.0
      SUMSQR=0.0
      COEFVR=0.0
*
      DO IR=1,(NCELL-1)         !look at each pair of spikes in the train
         TRV=SPIKETIMES(IR+1,IDs(ICELPT)) -
     +        SPIKETIMES(IR,IDs(ICELPT)) !TRV = interspike interval
         IF(TRV.LE.ZMAX)NINT=NINT+1 !this interval is OK, so count it
         IF(TRV.GT.ZMAX)cycle   !this interval too great - reject
         RQ(NINT)=TRV           !load RQ() with durations of acceptable intervals
         SUMTRV=SUMTRV+TRV      !running total of acceptable interval durations
         if (NINT.ge.ICT)goto 2001 !don't accumulate more than 100 acceptable intervals
      end do
*
      if((NINT.eq.0).or.(NINT.eq.1).or.(SUMTRV.eq.0.0))then !in the event of "no accepted cycles..."
         COEFVR = 0.0
         return
      end if
*
*       ***** CALCULATE MEAN,S.D., COEF.OF VARIANCE: *****
*
 2001 ZMEAN=SUMTRV/NINT         !ZMEAN = mean duration of this sample of acceptable intervals
*
      DO i=1,NINT                     
         SUMSQR=((RQ(i)-ZMEAN)**2)+SUMSQR
      end do
*
      STDDEV=SQRT(SUMSQR/(NINT-1)) !compute standard deviation
      COEFVR=STDDEV/ZMEAN       !compute coefficient of variance
*
*       *************************************************************
*
*       PRINT '(2X,I5,'' INTERVALS USED'',2X,/,
*     +         ''MAX INTERVAL ACCEPTED= '',F10.1,//,2X,
*     +         ''cell ID = '',I3,/,2x,
*     +         ''SUMTRV = '',f10.5,2x,
*     +         ''MEAN= '',F10.1,2X,''S.D.= +/-'',F10.1,2X,
*     +         ''COEF. OF VAR.='',F10.4)', 
*     +    NINT,ZMAX,IDs(ICELPT),ZMEAN,STDDEV,COEFVR
*       read ('(A)')
*
      return
      END
      end module mod_intvar4_v2a
