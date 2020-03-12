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

      module mod_calc_cardCCH
      contains
*       *****************************************************************
*       *                                                               *
*       *       Filename = calc_cardCCH.f                               *
*       *         (derived from calc_cardCTH_v2a.f  22-May-2002 lss)    *
*       *                                                               *
*       *       date of last revision:   04-Nov-2004    lss             *
*       *                                                               *
*       *       may-2002        lss                                     *
*       *        calculate CCH (no offset) instead of CTH               *
*       *                                                               *
*       *       apr-2001        lss                                     *
*       *        only those cardiac events which occur within the       *
*       *         respiratory cycles used to calculate resp CTHs        *
*       *         are used to create cardiac CTHs                       *
*       *                                                               *
*       *       mar-2000        lss                                     *
*       *       modified to allow ID codes to range from 1 to 999;      *       
*       *        maximum number of units allowed still = 120;           *       
*       *        *.defs files inserted                                  *
*       *        (MAX_NUM_CODES=999; MAX_NUM_CHAN=120)                  *
*       *                                                               *
*       *       link with x2002 code                                    *
*       *                                                               *
*       *       This subroutine of x2002 calculates and builds          *
*       *       cardiac CCHs.                                           *
*       *                                                               *
*       *****************************************************************
*
*
      SUBROUTINE calc_cardCCH(SPIKETIMES,ITAL,BINW,IHIST,
     +     LAST_REF,LAST_TAR,IRPT,ITPT,
     +     IDs,selected_C,
     +     total_sel_C_pulses)
*
      INCLUDE 'x2000parameter.defs'

      DOUBLE PRECISION
     +     SPIKETIMES(:,:),
     +     selected_C(:)
      DOUBLE PRECISION, allocatable ::
     +     selected_C_tmp(:),
     +     TARtmp(:)
*
      integer IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN)
      dimension ILH(IHALF),IRH(IHALF),IHIST(NUM_BINS)
*
      DOUBLE PRECISION
     +     DELT,TREF,TTAR,
     +     WL,WR,DL,HALFB,BINWtmp

*
      INTEGER LAG,R,T,ITB,LAST_REF,LAST_TAR,
     +     total_sel_C_pulses

      SAVE selected_C_tmp
      max_num_events = size (spiketimes, 1)
      if(.not.allocated(selected_C_tmp)) allocate (selected_C_tmp(MAX_NUM_EVENTS))
      allocate (TARtmp(MAX_NUM_EVENTS))
*
*
*       ***** initialize the arrays and variables *****
*
      do i = 1, NUM_BINS
         IHIST(i)=0
      end do
      do i = 1, IHALF
         ILH(i)=0
         IRH(i)=0
      end do
      IZER=0
      LAST_REF=0
      LAST_TAR=0
      BINWtmp=IDNINT(BINW*10D0)
      HALFB=(BINWtmp/2.0D0)
c       WR=(BINW*50.0)+HALFB            !WL and WR = left and right "edges" of the histogram
      WR=(BINWtmp*float(IHALF))+HALFB !WL and WR = left and right "edges" of the histogram
      WL=-(WR)
c       WINDH=1.0/BINW                  !WINDH = fraction of total histogram
*                                          window occupied by one binwidth
      ITB=1                     !initial pointer value for GET_TAR do loop.
!  A variable pointer is used to enhance the
!  efficiency of this subroutine by preventing
!  the comparison of every target event to each
!  reference event.

*
*       ***** if calculating the cardiac pulse overlay (cardiac pulse vs. cardiac pulse): *****
*
      if (IRPT.eq.ITPT) then  
         do i = 1, total_sel_C_pulses
            selected_C_tmp(i)=IDNINT(selected_C(i)*10D0)
         end do
         GET_REF: do R=1,total_sel_C_pulses
c          TREF = IDNINT(selected_C(R)*10D0)
         TREF = selected_C_tmp(R)
         LAST_REF=LAST_REF+1    !LAST_REF = a counter of # of ref events used in this CCH
         GET_TAR: do T=ITB,total_sel_C_pulses
c           TTAR = IDNINT(selected_C(T)*10D0)
         TTAR = selected_C_tmp(T)
         DELT = TTAR-TREF       !DELT = time elasped between REF and TAR event
*           if(DELT.eq.0)cycle GET_TAR          !same event - do not use this target event; get next one
         IF(DELT .LE. WL) then  !target event occurs "too early", so
            ITB=T               !  reset the do loop variable and
            cycle GET_TAR       !  get another TAR event
         end if
         IF(DELT .GE. WR) cycle GET_REF !target event occurs "too late" - get next REF
         LAST_TAR=LAST_TAR+1    !have a valid event - update counter
         DL=ABS(DELT)-HALFB
         IF(DL .lt. 0) then
            IZER=IZER+1         !an event which "occurs at 0"
            cycle GET_TAR       !get next TAR event
         end if
c           LAG = (DL*WINDH)                    !calculate which bin of the histogram
c           LAG = LAG + 1                       !  should be increased
         LAG = (DL/BINWtmp)+1   !calculate which bin of the histogram
         IF(DELT .GT. 0.0) then
            IRH(LAG) = IRH(LAG) + 1
            cycle GET_TAR
         end if
         IF(DELT.lt.0.0)then
            ILH(LAG) = ILH(LAG) + 1
            cycle GET_TAR
         end if
         
      end do GET_TAR
      end do GET_REF
*
c        IHIST(IHALF:1:-1)=ILH(1:IHALF) !load the histogram array
c        IHIST(IHALFplus1+1:NUM_BINS)=IRH(1:IHALF)

      do i = 1, 50
         IHIST(i) = ILH(IHALF-i+1)
         IHIST(i+IHALFplus1) = IRH(i)
      end do
      IHIST(IHALFplus1)=IZER
*
      RETURN
      end if
      
      
*       **** if calculating cardiac pulse vs. a unit: *****

      if(IRPT.ne.ITPT)then
         do i = 1, ITAL(IDs(ITPT))
            TARtmp(i)=IDNINT(SPIKETIMES(i,IDs(ITPT))*10D0)
         end do

         GET_REF1: DO R=1,total_sel_C_pulses !get a cardiac pulse and ...
         TREF=selected_C_tmp(R) !TREF = the time of the cardiac ref event
         LAST_REF=LAST_REF+1    !LAST_REF = a counter of # of ref events used in this CCH
         GET_TAR1: DO T=ITB,ITAL(IDs(ITPT)) !get a target event
c           TTAR = SPIKETIMES(T,IDs(ITPT))
         TTAR = TARtmp(T)
         DELT = TTAR-TREF       !DELT = time elasped between REF and TAR event
         IF(DELT .LE. WL) then  !target event occurs "too early", so
            ITB=T               !  reset the do loop variable and
            cycle GET_TAR1      !  get another TAR event
         end if
         IF(DELT .GE. WR) cycle GET_REF1 !target event occurs "too late" - get next REF
         LAST_TAR=LAST_TAR+1    !have a valid event - update counter
         DL=ABS(DELT)-HALFB
         IF(DL .lt. 0) then
            IZER=IZER+1         !an event which "occurs at 0"
            cycle GET_TAR1      !get next TAR event
         end if
c           LAG = (DL*WINDH)                    !calculate which bin of the histogram
c           LAG = LAG + 1                       !  should be increased
         LAG = (DL/BINWtmp)+1   !calculate which bin of the histogram
         IF(DELT .GT. 0) then
            IRH(LAG) = IRH(LAG) + 1
            cycle GET_TAR1
         end if
         IF(DELT.lt.0)then
            ILH(LAG) = ILH(LAG) + 1
            cycle GET_TAR1
         end if
         
      end do GET_TAR1
      end do GET_REF1
*
c        IHIST(50:1:-1)=ILH(1:50)       !load the histogram array
c        IHIST(52:NUM_BINS)=IRH(1:50)

      do i = 1, 50
         IHIST(i) = ILH(IHALF-i+1)
         IHIST(i+IHALFplus1) = IRH(i)
      end do
      IHIST(IHALFplus1)=IZER

      IHIST(51)=IZER
*
      RETURN
      end if
      END
      end module mod_calc_cardCCH
