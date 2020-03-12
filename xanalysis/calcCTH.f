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

      module mod_calcCTH
      contains
*       *************************************************************************
*       *                                                                       *
*       *       filename = calcCTH.f                                            *
*       *       lss 04-Nov-2004                                                 *
*       *mar-2000       lss
*        modified to allow ID codes to range from 1 to 999;     
*        maximum number of units allowed still = 120;           
*        *.defs files inserted                                  
*        (MAX_NUM_CODES=999; MAX_NUM_CHAN=120)                  
*       *** INDIRECT POINTERS are now used to access unit data
*       ***  within the following arrays:
*       ***     SPIKETIMES, ITAL, card_type, coef, coefnum, TALLY_NORM,
*       ***     ETA2_*, zmodsig_*, zmodsig2_*, analyzed_cells, analyzed_pairs,
*       ***     resp_type, card, CELL_NAMES
*       ***     [these arrays are now dimensioned to MAX_NUM_CHAN;
*       ***      the location of a cell's information (user-assigned ID code = i)
*       ***      within these arrays is array(IDs(i))]
*       ***
*       *** DIRECT POINTERS are used with these arrays: exclude, included
*       ***     (these arrays are dimensioned to MAX_NUM_CODES)
**      *                                                                       *
*       *       This subroutine of Xanalysis calculates and builds              *
*       *       cycle-triggered histograms.  Only those cycles which have       *
*       *       been previously selected by the user will be used in these      *
*       *       calculations.  A temporary array is created that contains only  *
*       *       events which occur in selected cycles.  In the event that       *
*       *       selected cycles are not sequential, event times are adjusted    *
*       *       to create an artificially sequential data flow with no gaps     *
*       *       between selected cycles.  This will prevent the inclusion of    *
*       *       events from non-selected cycles.                                *
*       *                                                                       *
*       *       The user has the option to override the recommended binwidth    *
*       *         (=((mean_resp_cycle+2*std_dev)+(1/2)(mean_I+2*std_dev_I))/100)*
*       *       and the recommended offset (=-(mean_I+(2*std_dev_I))).          *
*       *                                                                       *
*       *************************************************************************
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

      SUBROUTINE calcCTH (SPIKETIMES,ITAL,BINW,IHIST,LAST_REF,
     +     LAST_TAR,REF,TAR,ZOFST,ISTOPT,FIRST_ORDER,
     +     begin_E,end_E,TOTAL_SELECTED_CYCLES,IDs)
*
*
      INCLUDE 'x2000parameter.defs'
*
      double precision SPIKETIMES(:,:)
*
      integer IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN)
      DIMENSION IHIST(NUM_BINS)
      real ZOFST
*
      DOUBLE PRECISION, allocatable ::
     +     SELECTED_DATA(:)

*
      INTEGER LAST_REF,LAST_TAR,FIRST_ORDER,REF,TAR,
     +     TOTAL_SELECTED_CYCLES
      integer*4 LAG
      DOUBLE PRECISION E_OFFSET,
     +     begin_E(:), end_E(:),
     +     WIND, begin_E_tmp,end_E_tmp,TARtmp,ZOFSTtmp
      DOUBLE PRECISION, allocatable :: 
     +     new_begin_E(:),
     +     new_end_E(:)
*
*
*
*
*       ***** suppress unused variable warnings *****
      if(.false.)print *,FIRST_ORDER
      if(.false.)print *,REF
      allocate (new_begin_E(size(begin_E)))
      allocate (new_end_E(size(end_E)))
      allocate (SELECTED_DATA(size(spiketimes, 1)))

      WIND=IDNINT(BINW*10D0)*100.0D0 !WIND = total CTH time "window" (CTHs have 100 bins)

*       WINDH=1.0/BINW                  !WINDH = fraction of total histogram window occupied by one bin
      do i = 1, NUM_BINS
         IHIST=0                !initialize the histogram
      end do
      LAST_REF=0                !initialize these counters
      LAST_TAR=0
      ZOFSTtmp = IDNINT(ZOFST*10D0)
      SELECTED_DATA = 0.0    !initialize the temporary arrays
      new_begin_E=0.0
      new_end_E=0.0
*
*       *****  create a temporary array to hold events that occur within  *****
*       *****   selected cycles                                           *****
*
      begin_E_tmp = IDNINT(begin_E(1)*10D0)
      end_E_tmp = IDNINT(end_E(1)*10D0)
      j=0
      do i= 1, ITAL(IDs(TAR))   !transfer data from first
         TARtmp = IDNINT(SPIKETIMES(i,IDs(TAR))*10D0)
         if(TARtmp.ge.end_E_tmp)exit ! selected cycle to the
         if((TARtmp.ge.begin_E_tmp).and.(TARtmp.lt.end_E_tmp))then ! temporary array
            j=j+1
            SELECTED_DATA(j)=TARtmp
         end if
      end do
      new_begin_E(1) = begin_E_tmp
      new_end_E(1) = end_E_tmp

*
      do i = 2, TOTAL_SELECTED_CYCLES !now take care of the remaining cycles
         begin_E_tmp = IDNINT(begin_E(i)*10D0)
         end_E_tmp = IDNINT(end_E(i)*10D0)
         do k = 1, ITAL(IDs(TAR))
            TARtmp = IDNINT(SPIKETIMES(k,IDs(TAR))*10D0)
            if((TARtmp.ge.begin_E_tmp).and. !check to see if event occurs
     +           (TARtmp.lt.end_E_tmp))then !  within a selected cycle
               j=j+1
               SELECTED_DATA(j)=TARtmp-(begin_E_tmp-new_end_E(i-1))
            end if
         end do
         new_begin_E(i)=new_end_E(i-1)
         new_end_E(i)=new_begin_E(i)+(end_E_tmp-begin_E_tmp)
      end do
*
*       *****  Now have an array filled only with events from selected cycles   *****
*       *****     and any gaps which may have been in the data due to picking   *****
*       *****     and choosing respiratory cycles have been removed.            *****

      GET_E:  do i = 2, TOTAL_SELECTED_CYCLES
      E_OFFSET = new_begin_E(i)+ZOFSTtmp !add offset to this selected E_pulse 
      if((ISTOPT.ne.0).and.(i.ge.ISTOPT))goto 5 !have used specified # of cycles
*
*       ***** now have a "good" E-pulse, so... *****
*
      BUILD_CTH:  do j=1,ITAL(IDs(TAR)) !get target unit spike 
      if(SELECTED_DATA(j).lt.E_OFFSET)cycle BUILD_CTH !unit spike < E pulse
      if((SELECTED_DATA(j)-E_OFFSET).gt.WIND)cycle GET_E !spike occurs too long after
! the (possibly offset) E pulse
! to be included in CTH
*
*       ***** now have a good E-pulse and a good target event, so... *****
*
      LAST_TAR=LAST_TAR+1       !LAST_TAR = # of unit events used
      LAG=((SELECTED_DATA(j)-E_OFFSET)/IDNINT(BINW*10D0))+1 !determine which bin will be incremented and...
      IHIST(LAG)=IHIST(LAG)+1   !  increment it
c            IF(FIRST_ORDER.NE.0)cycle GET_REF ! currently first order not used
      end do BUILD_CTH
      end do GET_E
 5    LAST_REF=i                !total # of E pulses used to make this CTH
      IHIST(NUM_BINS)=0
      RETURN
      END
      end module mod_calcCTH
