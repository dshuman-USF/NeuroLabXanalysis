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

      module mod_respsig6_v4
      contains
*       filename = respsig6_v4.f        
*
*
*       september-2001  lss
*         calculation of ETA**2, ANOVA, and BINARY tests revised
*
*       mar-2000        lss
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
*
*       link with x2000_v* code
*
*       This subroutine of x2000_v* returns an indication of the significance
*               of the respiratory modulation of a cell's activity pattern.
*
*       1.  The ETA**2 value and ANOVA test will be calculated using the 
*           first 50 acceptable cycles of the control period.
*           An acceptable cycle is one during which the cell
*           is active; all zero-activity cycles will be eliminated.
*           If 50 acceptable cycles are not available, the user will be
*           informed that the test has been performed using less than the
*           desired number of respiratory cycles.  (isubETA = # of cycles used
*           to calculate the ETA**2 value and ANOVA test; isubETA is < or = 50)
*            
*       2.  The BINARY test splits the respiratory cycle into 2 halves
*           according to the cell's activity pattern so that the difference
*           between the two parts is maximal (the cycle is originally divided
*           into 20 "aisles" [=itrt]; each half is comprised of 10 contiguous aisles).
*           This means that the next cycle may be required to "donate" some aisles. For that
*           reason, we want to look at groups of 50+1 cycles when testing the hypothesis; the
*           remainder of the aisles in the 51st cycle are ignored, resulting in the effective
*           "dropping" of the 51st cycle from binary testing.
*
*           The BINARY test is most robust if 50 cycles are
*           used to form the hypothesis of a respiratory pattern, with the next at least
*           51 cycles (preferably more) used to confirm or deny that hypothesis.
*           The BINARY test will compare all the remaining cycles (preferably in groups of 51)
*           to the first group of 50 cycles.  
*           Again, only those cycles during which the cell fires are used. 
*           If at least 101 acceptable cycles are not available, the BINARY test will
*           be performed using half of the available cycles to form the hypothesis, and
*           the remaining half to test the hypothesis.  In this case, the user will be
*           warned that 101 cycles were not available for a robust BINARY test.
*           (isubBIN = # of cycles used to form the hypothesis for the BINARY test;
*            isubBIN is < or = 50)
*
*               
*       original code by km
*
*       Respiratory significance is evaluated under several sets 
*         of circumstances:
*           1. using the entire data set (case 1)
*               a. the user eliminate cycles whose duration is
*                  greater than a user-defined value (=MAXCYCLE) (case 2)
*
*           2. using only those cycles within a control period of the recording:
*               (control period may be defined with a boundary code, as the first
*                 x number of cycles in a recording, or a period of time within
*                 the recording)
*               a. use all events (mean_control_cycle = mean control cycle duration)
*                   (case 3)
*               b. cycles whose duration is less than the mean + 2 SD of 
*                  all control cycles (=control_max_value) (case 4)
*               c. the user may eliminate cycles whose duration is greater than
*                  a user-defined value (= user_max_value) (case 5)
*               (NOTE: if a boundary code has not been entered, 
*                      cases 2a, 2b, and 2c will be omitted)
*
*           3. using only those cycles actually used to compute CTHs (case 6)
*
*       ONLY THOSE CYCLES DURING WHICH THE CELL IS ACTIVE ARE USED TO COMPUTE SIGNIFICANCE
*
*       itrt = # of "aisles" into which each cycle will be subdivided 
*       **********************************************************************************
*       **********************************************************************************
*
*
      subroutine respsig5(SPIKETIMES,CELL,E_pulse,
     +     mean_control_cycle,MAXCYCLE,
     +     control_max_value,user_max_value,
     +     mean_resp_cycle,begin_E,end_E,
     +     BNDRY_START_TIME,BNDRY_END_TIME,
     +     ETA2_1,zmodsig_1,zmodsig2_1,
     +     ETA2_2,zmodsig_2,zmodsig2_2,
     +     ETA2_3,zmodsig_3,zmodsig2_3,
     +     ETA2_4,zmodsig_4,zmodsig2_4,
     +     ETA2_5,zmodsig_5,zmodsig2_5,
     +     ETA2_6,zmodsig_6,zmodsig2_6,
     +     IDs,mode)
*
*
      use mod_betai
      INCLUDE 'x2000parameter.defs'

      double precision SPIKETIMES(:,:)
*
      integer IDs(MAX_NUM_CODES)
      PARAMETER (SIGLVL=.05,ZSIGLVL=1.645,itrt=20,
     +     isubmax=51,isub=50)
*
      DIMENSION IY(itrt,isubmax), JSUM(itrt)

      DOUBLE PRECISION
     +     begin_E(:),end_E(:)
      DOUBLE PRECISION
     +     cycle_start_time,cycle_end_time,CLLK,
     +     BNDRY_START_TIME,BNDRY_END_TIME
      DOUBLE PRECISION, allocatable ::
     +     temp_begin_E(:),
     +     temp_end_E(:)

      REAL*4 PROB,mean_control_cycle,
     +     control_max_value,user_max_value,
     +     z,ETA2,MAXCYCLE,
     +     mean_resp_cycle
*
      integer E_pulse,E_pointer,CELL,CELL_pointer,
     +     cycle,aisle
*       
      character*5 ETA2_1(MAX_NUM_CHAN),
     +     ETA2_2(MAX_NUM_CHAN),
     +     ETA2_3(MAX_NUM_CHAN),
     +     ETA2_4(MAX_NUM_CHAN),
     +     ETA2_5(MAX_NUM_CHAN),
     +     ETA2_6(MAX_NUM_CHAN)
*
      character*3 zmodsig_1(MAX_NUM_CHAN),
     +     zmodsig_2(MAX_NUM_CHAN),
     +     zmodsig_3(MAX_NUM_CHAN),
     +     zmodsig_4(MAX_NUM_CHAN),
     +     zmodsig_5(MAX_NUM_CHAN),
     +     zmodsig_6(MAX_NUM_CHAN),
     +     zmodsig2_1(MAX_NUM_CHAN),
     +     zmodsig2_2(MAX_NUM_CHAN),
     +     zmodsig2_3(MAX_NUM_CHAN),
     +     zmodsig2_4(MAX_NUM_CHAN),
     +     zmodsig2_5(MAX_NUM_CHAN),
     +     zmodsig2_6(MAX_NUM_CHAN)

      character*2 mode
      integer :: isubmaxBIN = 0   

*       ***** suppress unused variable warnings *****
      if(.false.)print *,mean_resp_cycle
      if(.false.)print *,mean_control_cycle

      max_num_events = size (spiketimes, 1)
      allocate (temp_begin_E(MAX_NUM_EVENTS))
      allocate (temp_end_E(MAX_NUM_EVENTS))
      max_cycles = size(end_E)
      ETA2_1(IDs(CELL))='     ' !initialize values for this cell only
      ETA2_2(IDs(CELL))='     '
      ETA2_3(IDs(CELL))='     '
      ETA2_4(IDs(CELL))='     '
      ETA2_5(IDs(CELL))='     '
      ETA2_6(IDs(CELL))='     '
      zmodsig_1(IDs(CELL))='  '
      zmodsig_2(IDs(CELL))='  '
      zmodsig_3(IDs(CELL))='  '
      zmodsig_4(IDs(CELL))='  '
      zmodsig_5(IDs(CELL))='  '
      zmodsig_6(IDs(CELL))='  '
      zmodsig2_1(IDs(CELL))='  '
      zmodsig2_2(IDs(CELL))='  '
      zmodsig2_3(IDs(CELL))='  '
      zmodsig2_4(IDs(CELL))='  '
      zmodsig2_5(IDs(CELL))='  '
      zmodsig2_6(IDs(CELL))='  '
      temp_begin_E = 0.0
      temp_end_E = 0.0
      
      if(mode.eq.'vt')then
         ETA2_1 = 'n/a'
         ETA2_2 = 'n/a'
         ETA2_3 = 'n/a'
         ETA2_4 = 'n/a'
         ETA2_5 = 'n/a'
         zmodsig_1(IDs(CELL)) = 'n/a'
         zmodsig_2(IDs(CELL)) = 'n/a'
         zmodsig_3(IDs(CELL)) = 'n/a'
         zmodsig_4(IDs(CELL)) = 'n/a'
         zmodsig_5(IDs(CELL)) = 'n/a'
         zmodsig2_1(IDs(CELL)) = 'n/a'
         zmodsig2_2(IDs(CELL)) = 'n/a'
         zmodsig2_3(IDs(CELL)) = 'n/a'
         zmodsig2_4(IDs(CELL)) = 'n/a'
         zmodsig2_5(IDs(CELL)) = 'n/a'

         goto 602               !view-only mode - skip most resp significance tests
      end if
*
*
*
*       ***** This subroutine calculates a total of 6 significance values each time it is called *****
*       *****  in order to cover all possible circumstances as listed above.  *****
*
*
*       *********************************************************************
*       *****                                                           *****
*       *****           CASE 1 -- use all cycles                        *****
*       *****                                                           *****
*       *********************************************************************
*
*                                       !calculate respiratory significance using data
! from entire recording; no cycles eliminated
! due to cycle length - only if cell is not
! active during that cycle
*
*       ***** check to see how many cycles are acceptable:    *****
*       *****  (want min of 101 for most robust BINARY test)  *****
*
      j=0
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
      isubmaxBIN = 0            
*
 105  if(E_pointer.lt.MAX_NUM_EVENTS)then
         cycle_start_time =
     +        SPIKETIMES(E_pointer,IDs(E_pulse))
         cycle_end_time =
     +        SPIKETIMES(E_pointer+1,IDs(E_pulse))
         if(cycle_end_time.eq.0)goto 120 !look at all cycles
*
 106     if(CELL_pointer.le.MAX_NUM_EVENTS)then
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
*
            if (CLLK.lt.cycle_start_time) then
               CELL_pointer=CELL_pointer+1
               goto 106         !loop until find the FIRST valid cell event
! that occurs after the start of the cycle
            endif
*
            if (CLLK.gt.cycle_end_time) then !if the event occurs after the
               E_pointer=E_pointer+1 ! end of the cycle,
               goto 105         ! go to the next cycle
            endif
*
*       ***** at this point, we have detected the first event in the first *****
*       *****  acceptable cycle, so increment the counter:  *****
*
            j=j+1               !keep track of number of acceptable cycles
            E_pointer=E_pointer+1 !point to the next respiratory cycle
            goto 105            !get next cycle
         end if
*
      end if
*
 120  if(j.eq.0)then            !cannot do any significance testing for this cell --
*         ETA2=0.0              ! no cycles are acceptable!
*         PROB=1.0
*         goto __________               !????????????????????????????????
      end if
*
      if(j.ge.50)then           !50 acceptable cycles are available to calculate
         isubETA=50             ! the ETA**2 value and perform the ANOVA test
      else                      
         isubETA=j              !have to use less than 50 cycles - warn the user
      end if
*
      if(j.ge.101)then          !at least 101 acceptable cycles available to
         isubBIN=50             ! perform the BINARY test - good
         isubmaxBIN=isubBIN + 1 
      else                      !have to use less than 101 cycles - warn the user
         isubBIN=INT(j/2)
         isubmaxBIN=isubBIN + 1 
      end if     
      df1=float(itrt-1)                         
      df2=float(itrt*isubETA-itrt)      
*
*       ***** Calculate the ETA**2 value and the ANOVA statistic (F test): *****
*
*
      ETA2=0.0                  !default value
      PROB=1.0                  !default value
      z=0.0
      RMSB=0
      RMSW=0
      IY = 0                    !IY(itrt,isub) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=1                   ! = respiratory cycle designation
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Find the FIRST respiratory cycle in the recording: *****
*
      cycle_start_time=
     +     SPIKETIMES(E_pointer,IDs(E_pulse))   
      cycle_end_time=
     +     SPIKETIMES(E_pointer+1,IDs(E_pulse))
*
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
         ETA2=0.0
         PROB=1.0
         goto 160
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 100  CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 2'')'
*          ETA2=0.0             !set default values
*           PROB=1.0
            goto 160
         ENDIF
         GOTO 100               !loop until find the FIRST valid cell event
! that occurs after the start of this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 110  FIND_THE_AISLE_1: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 111     if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 3'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 150
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES ETA2 4'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 150         !go to next resp sig test
            endif
*
            GOTO 110            !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 111 !itrt = total # aisles per cycle
*
         ENDIF
      end if FIND_THE_AISLE_1
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubETA)goto 150
      E_pointer=E_pointer+1     !point to next E pulse
*
      cycle_start_time=
     +     SPIKETIMES(E_pointer,IDs(E_pulse)) !read times of new beginning and 
      cycle_end_time=
     +     SPIKETIMES(E_pointer+1,IDs(E_pulse)) ! ending E pulses
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''NOT ENOUGH CYCLES ETA2 5'')' 
*           ETA2=0.0            !set default values
*            PROB=1.0
         goto 150
      end if
*
 140  IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 140
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubETA) GOTO 110 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubETA (optimally 50) cycles,               *****
*       *****  divided each cycle into 20 aisles, and counted the number of         *****
*       *****  cell events that occur within each aisle;                            *****
*       ***** now we're ready for some statistical analysis of respiratory pattern: *****
*
 150  AVGY=float(IYSUM)/float(itrt*isubETA) !total # cell events in the 1st isubETA cycles /
! total # aisles in the 1st isubETA cycles =
!  # cell events expected in each aisle if
! the cell's activity is NOT respiratory modulated
      DO aisle=1,itrt
         RMSB=RMSB+(float(JSUM(aisle))
     +        /float(isubETA)-AVGY)**2
      end do
      RMSB=isubETA*RMSB/df1     !RMSB = between group sum of squares (?)
*
      DO aisle=1,itrt
         DO cycle=1,isubETA
            RMSW=RMSW+(float(IY(aisle,cycle))
     +           -float(JSUM(aisle))/float(isubETA))**2
         end do
      end do
*
      RMSW=RMSW/df2             !RMSW = within group sum of squares (?)
*
      if(RMSW.NE.0) then
         F=RMSB/RMSW
      else
         F=1.0
      end if
*
      if(RMSW.ne.0)then
         ETA2=df1*F/(df1*F+df2) !ETA squared statistic
      else
         ETA2=0.0
      end if
*
      PROB = BETAI(0.5*df2,0.5*df1,df2/(df2+df1*F))
 160  write (ETA2_1(IDs(CELL)),('(F5.2)')) ETA2
      if(PROB.lt.SIGLVL)then
         zmodsig_1(IDs(CELL)) = 'R' !cell's activity is respiratory modulated
      else if((PROB.ge.SIGLVL).and.(isubETA.lt.50))then
         zmodsig_1(IDs(CELL)) = 'NED' !ned = not enough data to confidently say cell is not resp. mod.
      else if((PROB.ge.SIGLVL).and.(isubETA.ge.50))then
         zmodsig_1(IDs(CELL)) = 'N' !n = not respiratory modulated
      else
      end if
*
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*         <<<<< PROCEED WITH BINARY TEST for case 1: >>>>>
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*
*
      z=0.0
      ICNT1=0
      ICNT2=0
      ITEMPSUM1=0
      ITEMPSUM2=0
      IY = 0                    !IY(itrt,isub) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=0                   ! = respiratory cycle designation
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Find the first selected respiratory cycle in the recording: *****
*
      cycle_start_time=
     +     SPIKETIMES(E_pointer,IDs(E_pulse)) !read times of new beginning and 
      cycle_end_time=
     +     SPIKETIMES(E_pointer+1,IDs(E_pulse)) ! ending E pulses
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
*           PRINT '(''NOT ENOUGH CYCLES 1'')' 
         goto 1099
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 1000 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES 2'')'
            goto 1099
         ENDIF
         GOTO 1000              !loop until find the FIRST valid cell event
! in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 1010 FIND_THE_AISLE_1_2: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 1011    if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN !ran out of cell events
               goto 1099
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES 4'')' 
               goto 1099        !go to next resp sig test
            endif
*
            GOTO 1010           !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 1011 !itrt = total # aisles per cycle
*
         end if
      END IF FIND_THE_AISLE_1_2
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubBIN)goto 1041 !proceed to BINARY analysis
      E_pointer=E_pointer+1     !point to next E pulse
*
      cycle_start_time=
     +     SPIKETIMES(E_pointer,IDs(E_pulse)) !read times of new beginning and 
      cycle_end_time=
     +     SPIKETIMES(E_pointer+1,IDs(E_pulse)) ! ending E pulses
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''NOT ENOUGH CYCLES 5'')' 
         goto 1099
      end if
*
 1040 IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 1040
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubBIN) GOTO 1010 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubBIN (optimally 50) cycles,                   *****
*       *****  divided each one into 20 aisles,                                         *****
*       *****  and counted the number of cell events that occur within each aisle;      *****
*
*       ***** make 2 groups of "contiguous" aisles                      *****
*       ***** (ie., 1-10 & 11-20; 2-11 & 12-20 + 1; etc.);              *****
*       ***** Examine all possible groupings until find the case in     *****
*       *****  which the difference in total # cell events between the  *****
*       *****  two groups is maximal:                                   *****
*
 1041 IPT1=(itrt/2) - 1         ! = # of first aisle in the 1st group
      IPT2=itrt - 1             ! = # of first aisle in the 2nd group
      IOLDDIFF=0
      IBIGPT=1
      DO I=1,itrt/2             !update pointers to first aisle in each group
         IPT1=IPT1 - itrt/2 + 1
         IPT2=IPT2 - itrt/2 + 1
         IF (IPT2.LE.0) IPT2=IPT2 + itrt
         ITEMPSUM1=0            !ITEMPSUM1 = total # cell events in 1st group
         ITEMPSUM2=0            !ITEMPSUM2 = total # cell events in 2nd group
*
         DO J=1,itrt/2
            IPT1=IPT1+1         !point to next aisle in each group
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) IPT2=IPT2-itrt
            ITEMPSUM1=ITEMPSUM1+JSUM(IPT1)
            ITEMPSUM2=ITEMPSUM2+JSUM(IPT2)
         end do
*
         IF (ITEMPSUM1.GT.ITEMPSUM2) THEN !compare the #s of events which occurred in the
            INEWDIFF=ITEMPSUM1-ITEMPSUM2 ! 1st group vs. those in the 2nd group
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT1 - itrt/2 + 1 !store # of 1st aisle in group 1
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ELSE
            INEWDIFF=ITEMPSUM2-ITEMPSUM1 !ITEMPSUM2 > ITEMPSUM1
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT2 - itrt/2 + 1 !store # of 1st aisle in group 2
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ENDIF
*
      end do
*
      totalBINARY=isubBIN       !keep track of total # cycles used in BINARY testing
*
*       ***** Now have formed a hypothesis and gathered data based on the first isubBIN cycles *****
*       *****  about whether or not the cell's activity is greater in one half of the          *****
*       *****  cycle than in the other half.  Now need to look at remaining                       *****
*       *****  cycles to see if the relationship/pattern holds true:                           *****
*
*
      ICNT1=0                   != # times ITEMPSUM1 > ITEMPSUM2
      ICNT2=0                   != # times ITEMPSUM1 <> ITEMPSUM2
 1045 IY = 0
      cycle=0
      aisle=1
*
*      ***** move to the next set of E pulses (hence, the next cycle): *****
*      
      E_pointer=E_pointer+1
*
      cycle_start_time=
     +     SPIKETIMES(E_pointer,IDs(E_pulse)) !read times of new beginning and 
      cycle_end_time=
     +     SPIKETIMES(E_pointer+1,IDs(E_pulse)) ! ending E pulses
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*       PRINT '(''OUT OF CYCLES 6'')' 
         goto 1099
      end if
*
*      ***** find the first cell event after the first cycle start: *****
*
 1060 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN !ran out of cell events
*         PRINT '(''OUT OF CYCLES 7'')'
            GOTO 1099
         ENDIF
         GOTO 1060              !loop until find 1st cell event in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 1070 FIND_THE_AISLE_1_3: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
         if (CLLK.LT.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/itrt))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts events in array subscripted
*                                                       by aisle # and cycle #
            CELL_pointer=CELL_pointer+1
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*         PRINT '(''OUT OF CYCLES 8'')'
               GOTO 1098
            ENDIF
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event
*
            IF (CLLK.EQ.0) THEN
*         PRINT '(''OUT OF CYCLES 9'')'
               GOTO 1098
            ENDIF
*
            GOTO 1070           !loop back to analyze the new event time
*
         ELSE                   !cell event does not occur in current aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 1070 !look at the next aisle (or go to the next cycle)
*
         end if
      END IF FIND_THE_AISLE_1_3
*
      if(cycle.eq.isubmaxBIN)goto 1098
      E_pointer=E_pointer+1     !proceed to the next resp cycle
*
      cycle_start_time=
     +     SPIKETIMES(E_pointer,IDs(E_pulse)) !read times of new beginning and 
      cycle_end_time=
     +     SPIKETIMES(E_pointer+1,IDs(E_pulse)) ! ending E pulses
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''OUT OF CYCLES 10'')' 
         goto 1098
      end if
*
 1076 IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
         IF (CLLK.EQ.0) GOTO 1098
         GOTO 1076              !loop until read a valid cell event time
      ENDIF
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubmaxBIN) GOTO 1070 !for the next cycle
*
*       ***** look at cell's activity within this group of cycles to see if   *****
*       *****  the relationship/pattern detected with the 1st isubBIN cycles holds true: *****
*
 1098 DO icycle=1,cycle-1       !for each cycle in the group:
         totalBINARY = totalBINARY+1 !keep a total of # of cycles used in BINARY testing
         ITEMPSUM1=0
         ITEMPSUM2=0
         IPT1=IBIGPT-1
         IPT2=IBIGPT+itrt/2-1
         IF (IPT2.GT.itrt) IPT2=IPT2-itrt
         DO J=1,itrt/2
            IPT1=IPT1+1
            IF (IPT1.GT.itrt) THEN
               ITEMPSUM1=ITEMPSUM1+IY(IPT1-itrt,icycle+1)
            ELSE
               ITEMPSUM1=ITEMPSUM1+IY(IPT1,icycle)
            ENDIF
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) THEN
               ITEMPSUM2=ITEMPSUM2+IY(IPT2-itrt,icycle+1)
            ELSE
               ITEMPSUM2=ITEMPSUM2+IY(IPT2,icycle)
            ENDIF
         end do
*
         IF (ITEMPSUM1.GT.ITEMPSUM2) ICNT1=ICNT1+1 !ICNT1=# cycles which support the hypothesis
         IF (ITEMPSUM1.NE.ITEMPSUM2) ICNT2=ICNT2+1 !ICNT2=# cycles in which the cell's activity in
*                                                       ! one half of the cycle does not equal the
*                                                       ! activity in the other half of the cycle
      end do
*
      goto 1045                 
*
 1099 z=(float(ICNT1)-float(ICNT2)/2.)/(.5*SQRT(FLOAT(ICNT2)))
      if(z.gt.ZSIGLVL)then
         zmodsig2_1(IDs(CELL)) = 'R' ! bsig = Z value
      else if((z.le.ZSIGLVL).and.(isubBIN.lt.50))then
         zmodsig2_1(IDs(CELL)) = 'NED' !ned = not enough data for a good test
      else if((z.le.ZSIGLVL).and.(isubBIN.ge.50))then
         zmodsig2_1(IDs(CELL)) = 'N'
      else
      end if

*       print '(''case 1 '',I3,'': ETA2  = '',A,''   ANOVA = '',A1,
*     +         ''   BINARY = '',A1)',
*     +    CELL,ETA2_1(IDs(CELL)),zmodsig_1(IDs(CELL)),
*     +    zmodsig2_1(IDs(CELL))        
*
*
*       *********************************************************************
*       *****                                                           *****
*       ***** CASE 2 -- use entire recording, but reject a cycle if it  *****
*       *****   is longer than MAXCYCLE (default: MAXCYCLE=10,000.0 ms) *****
*       *****                                                           *****
*       *********************************************************************
*
*
c202    if(MAXCYCLE.eq.1.0)then         !user does not want to use this feature
c           ETA2_2 = 'n/a'              !  enter "not applicable" for these variables
c           zmodsig_2(IDs(CELL)) = 'n/a'
c           zmodsig2_2(IDs(CELL)) = 'n/a'
c           goto 302
c           end if
*
*       ***** check to see how many cycles are acceptable:    *****
*       *****  (want min of 101 for most robust BINARY test)  *****
*
      j=0
      temp_begin_E = 0.0
      temp_end_E = 0.0
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
      isubmaxBIN = 0            
*
 205  if(E_pointer.lt.MAX_NUM_EVENTS)then
         cycle_start_time =
     +        SPIKETIMES(E_pointer,IDs(E_pulse))
         cycle_end_time =
     +        SPIKETIMES(E_pointer+1,IDs(E_pulse))
         if(cycle_end_time-cycle_start_time.gt.MAXCYCLE)then
            E_pointer=E_pointer+1 !this cycle too long, so get next one
            goto 205
         end if
         if(cycle_end_time.eq.0)goto 220 !look at all cycles
*
 206     if(CELL_pointer.le.MAX_NUM_EVENTS)then
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
*
            if (CLLK.lt.cycle_start_time) then
               CELL_pointer=CELL_pointer+1
               goto 206         !loop until find the FIRST valid cell event
! that occurs after the start of the cycle
            endif
*
            if (CLLK.gt.cycle_end_time) then !if the event occurs after the
               E_pointer=E_pointer+1 ! end of the cycle,
               goto 205         ! go to the next cycle
            endif
*
*       ***** at this point, we have detected the first event in the first *****
*       *****  acceptable cycle, so increment the counter:  *****
*
            j=j+1               !keep track of number of acceptable cycles
            temp_begin_E(j)=cycle_start_time !load array with times of acceptable E pulses
            temp_end_E(j)=cycle_end_time !   for case 2
            E_pointer=E_pointer+1 !point to the next respiratory cycle
            goto 205            !get next cycle
         end if
*
      end if
*
 220  if(j.eq.0)then            !cannot do any significance testing for this cell --
*         ETA2=0.0              ! no cycles are acceptable!
*         PROB=1.0
*         goto _____________            !????????????????????????????????
      end if
*
      if(j.ge.50)then           !50 acceptable cycles are available to calculate
         isubETA=50             ! the ETA**2 value and perform the ANOVA test
      else                      
         isubETA=j              !have to use less than 50 cycles - warn the user
      end if
*
      if(j.ge.101)then          !at least 101 acceptable cycles available to
         isubBIN=50             ! perform the BINARY test - good
         isubmaxBIN=isubBIN + 1 
      else                      !have to use less than 101 cycles - warn the user
         isubBIN=INT(j/2)
         isubmaxBIN=isubBIN + 1 
      end if     
      df1=float(itrt-1)                         
      df2=float(itrt*isubETA-itrt)      
*
*       ***** Calculate the ETA**2 value and the ANOVA statistic (F test): *****
*
*
      ETA2=0.0                  !default value
      PROB=1.0                  !default value
      z=0.0
      RMSB=0
      RMSW=0
      IY = 0                    !IY(itrt,isub) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=1                   ! = respiratory cycle designation
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Get the FIRST acceptable respiratory cycle in the recording: *****
*
*201    cycle_start_time=
*     +         SPIKETIMES(E_pointer,IDs(E_pulse))      
*       cycle_end_time=
*     +         SPIKETIMES(E_pointer+1,IDs(E_pulse))
      cycle_start_time=temp_begin_E(E_pointer)
      cycle_end_time=temp_end_E(E_pointer)
*
*       if(cycle_end_time-cycle_start_time.gt.MAXCYCLE)then
*           E_pointer=E_pointer+1               !this cycle too long, so get next one
*           goto 201
*           end if
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
         ETA2=0.0
         PROB=1.0
         goto 260
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 200  CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 2'')'
*          ETA2=0.0             !set default values
*           PROB=1.0
            goto 260
         ENDIF
         GOTO 200               !loop until find the FIRST valid cell event
! that occurs after the start of this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 210  FIND_THE_AISLE_2: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 211     if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 3'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 250
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES ETA2 4'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 250         !go to next resp sig test
            endif
*
            GOTO 210            !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 211 !itrt = total # aisles per cycle
*
         ENDIF
      end if FIND_THE_AISLE_2
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubETA)goto 250
      E_pointer=E_pointer+1     !point to next E pulse
*
*       cycle_start_time=
*     +         SPIKETIMES(E_pointer,IDs(E_pulse))      !read times of new beginning and 
*       cycle_end_time=
*     +         SPIKETIMES(E_pointer+1,IDs(E_pulse))    ! ending E pulses
*       if(cycle_end_time-cycle_start_time.gt.MAXCYCLE)then     !this cycle too long, so get next one
*           goto 230
*           end if
      cycle_start_time=temp_begin_E(E_pointer)
      cycle_end_time=temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''NOT ENOUGH CYCLES ETA2 5'')' 
*           ETA2=0.0            !set default values
*            PROB=1.0
         goto 250
      end if
*
 240  IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 240
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubETA) GOTO 210 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubETA (optimally 50) cycles,               *****
*       *****  divided each cycle into 20 aisles, and counted the number of         *****
*       *****  cell events that occur within each aisle;                            *****
*       ***** now we're ready for some statistical analysis of respiratory pattern: *****
*
 250  AVGY=float(IYSUM)/float(itrt*isubETA) !total # cell events in the 1st isubETA cycles /
! total # aisles in the 1st isubETA cycles =
!  # cell events expected in each aisle if
! the cell's activity is NOT respiratory modulated
      DO aisle=1,itrt
         RMSB=RMSB+(float(JSUM(aisle))
     +        /float(isubETA)-AVGY)**2
      end do
      RMSB=isubETA*RMSB/df1     !RMSB = between group sum of squares (?)
*
      DO aisle=1,itrt
         DO cycle=1,isubETA
            RMSW=RMSW+(float(IY(aisle,cycle))
     +           -float(JSUM(aisle))/float(isubETA))**2
         end do
      end do
*
      RMSW=RMSW/df2             !RMSW = within group sum of squares (?)
*
      if(RMSW.NE.0) then
         F=RMSB/RMSW
      else
         F=1.0
      end if
*
      if(RMSW.ne.0)then
         ETA2=df1*F/(df1*F+df2) !ETA squared statistic
      else
         ETA2=0.0
      end if
*
      PROB = BETAI(0.5*df2,0.5*df1,df2/(df2+df1*F))
 260  write (ETA2_2(IDs(CELL)),('(F5.2)')) ETA2
      if(PROB.lt.SIGLVL)then
         zmodsig_2(IDs(CELL)) = 'R' !cell's activity is respiratory modulated
      else if((PROB.ge.SIGLVL).and.(isubETA.lt.50))then
         zmodsig_2(IDs(CELL)) = 'NED' !ned = not enough data
      else if((PROB.ge.SIGLVL).and.(isubETA.ge.50))then
         zmodsig_2(IDs(CELL)) = 'N' !n = not respiratory modulated
      else
      end if
c       if(PROB.lt.SIGLVL)zmodsig_2(IDs(CELL)) = 'R'
*
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*         <<<<< PROCEED WITH BINARY TEST for case 2: >>>>>
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*
*
      z=0.0
      ICNT1=0
      ICNT2=0
      ITEMPSUM1=0
      ITEMPSUM2=0
      IY = 0                    !IY(itrt,isub) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=0                   ! = respiratory cycle designation
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Get the first acceptable respiratory cycle in the recording: *****
*
*2095   cycle_start_time=
*     +         SPIKETIMES(E_pointer,IDs(E_pulse))      !read times of new beginning and 
*       cycle_end_time=
*     +         SPIKETIMES(E_pointer+1,IDs(E_pulse))    ! ending E pulses
*       if(cycle_end_time-cycle_start_time.gt.MAXCYCLE)then
*           E_pointer=E_pointer+1               !this cycle too long, so get next one
*           goto 2095
*           end if
      cycle_start_time=temp_begin_E(E_pointer)
      cycle_end_time=temp_end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
*           PRINT '(''NOT ENOUGH CYCLES 1'')' 
         goto 2099
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 2000 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES 2'')'
            goto 2099
         ENDIF
         GOTO 2000              !loop until find the FIRST valid cell event
! in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 2010 FIND_THE_AISLE_2_2: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 2011    if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN !ran out of cell events
               goto 2099
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES 4'')' 
               goto 2099        !go to next resp sig test
            endif
*
            GOTO 2010           !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 2011 !itrt = total # aisles per cycle
*
         end if
      END IF FIND_THE_AISLE_2_2
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubBIN)goto 2041 !proceed to BINARY analysis
      E_pointer=E_pointer+1     !point to next E pulse
*
*       cycle_start_time=
*     +         SPIKETIMES(E_pointer,IDs(E_pulse))      !read times of new beginning and 
*       cycle_end_time=
*     +         SPIKETIMES(E_pointer+1,IDs(E_pulse))    ! ending E pulses
*
*       if(cycle_end_time-cycle_start_time.gt.MAXCYCLE)then     !this cycle too long, so get next one
*           goto 2030
*           end if
      cycle_start_time=temp_begin_E(E_pointer)
      cycle_end_time=temp_end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''NOT ENOUGH CYCLES 5'')' 
         goto 2099
      end if
*
 2040 IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 2040
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubBIN) GOTO 2010 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubBIN (optimally 50) cycles,                   *****
*       *****  divided each one into 20 aisles,                                         *****
*       *****  and counted the number of cell events that occur within each aisle;      *****
*
*       ***** make 2 groups of "contiguous" aisles                      *****
*       ***** (ie., 1-10 & 11-20; 2-11 & 12-20 + 1; etc.);              *****
*       ***** Examine all possible groupings until find the case in     *****
*       *****  which the difference in total # cell events between the  *****
*       *****  two groups is maximal:                                   *****
*
 2041 IPT1=(itrt/2) - 1         ! = # of first aisle in the 1st group
      IPT2=itrt - 1             ! = # of first aisle in the 2nd group
      IOLDDIFF=0
      IBIGPT=1
      DO I=1,itrt/2             !update pointers to first aisle in each group
         IPT1=IPT1 - itrt/2 + 1
         IPT2=IPT2 - itrt/2 + 1
         IF (IPT2.LE.0) IPT2=IPT2 + itrt
         ITEMPSUM1=0            !ITEMPSUM1 = total # cell events in 1st group
         ITEMPSUM2=0            !ITEMPSUM2 = total # cell events in 2nd group
*
         DO J=1,itrt/2
            IPT1=IPT1+1         !point to next aisle in each group
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) IPT2=IPT2-itrt
            ITEMPSUM1=ITEMPSUM1+JSUM(IPT1)
            ITEMPSUM2=ITEMPSUM2+JSUM(IPT2)
         end do
*
         IF (ITEMPSUM1.GT.ITEMPSUM2) THEN !compare the #s of events which occurred in the
            INEWDIFF=ITEMPSUM1-ITEMPSUM2 ! 1st group vs. those in the 2nd group
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT1 - itrt/2 + 1 !store # of 1st aisle in group 1
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ELSE
            INEWDIFF=ITEMPSUM2-ITEMPSUM1 !ITEMPSUM2 > ITEMPSUM1
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT2 - itrt/2 + 1 !store # of 1st aisle in group 2
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ENDIF
*
      end do
*
      totalBINARY=isubBIN       !keep track of total # cycles used in BINARY testing
*
*       ***** Now have formed a hypothesis and gathered data based on the first isubBIN cycles *****
*       *****  about whether or not the cell's activity is greater in one half of the          *****
*       *****  cycle than in the other half.  Now need to look at remaining                       *****
*       *****  cycles to see if the relationship/pattern holds true:                           *****
*
*
      ICNT1=0                   != # times ITEMPSUM1 > ITEMPSUM2
      ICNT2=0                   != # times ITEMPSUM1 <> ITEMPSUM2
 2045 IY = 0
      cycle=0
      aisle=1
*
*      ***** move to the next set of E pulses (hence, the next cycle): *****
*      
      E_pointer=E_pointer+1
*
*2050   cycle_start_time=
*     +         SPIKETIMES(E_pointer,IDs(E_pulse))      !read times of new beginning and 
*       cycle_end_time=
*     +         SPIKETIMES(E_pointer+1,IDs(E_pulse))    ! ending E pulses
*       if(cycle_end_time-cycle_start_time.gt.MAXCYCLE)then
*           E_pointer=E_pointer+1               !this cycle too long, so get next one
*           goto 2050
*           end if
      cycle_start_time=temp_begin_E(E_pointer)
      cycle_end_time=temp_end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*       PRINT '(''OUT OF CYCLES 6'')' 
         goto 2099
      end if
*
*      ***** find the first cell event after the first cycle start: *****
*
 2060 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN !ran out of cell events
*         PRINT '(''OUT OF CYCLES 7'')'
            GOTO 2099
         ENDIF
         GOTO 2060              !loop until find 1st cell event in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 2070 FIND_THE_AISLE_2_3: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
         if (CLLK.LT.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/itrt))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts events in array subscripted
*                                                       by aisle # and cycle #
            CELL_pointer=CELL_pointer+1
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*         PRINT '(''OUT OF CYCLES 8'')'
               GOTO 2098
            ENDIF
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event
*
            IF (CLLK.EQ.0) THEN
*         PRINT '(''OUT OF CYCLES 9'')'
               GOTO 2098
            ENDIF
*
            GOTO 2070           !loop back to analyze the new event time
*
         ELSE                   !cell event does not occur in current aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 2070 !look at the next aisle (or go to the next cycle)
*
         end if
      END IF FIND_THE_AISLE_2_3
*
      if(cycle.eq.isubmaxBIN)goto 2098
      E_pointer=E_pointer+1     !proceed to the next resp cycle
*
*       cycle_start_time=
*     +         SPIKETIMES(E_pointer,IDs(E_pulse))      !read times of new beginning and 
*       cycle_end_time=
*     +         SPIKETIMES(E_pointer+1,IDs(E_pulse))    ! ending E pulses
*
*       if(cycle_end_time-cycle_start_time.gt.MAXCYCLE)then !this cycle too long, so get next one
*           goto 2075
*           end if
      cycle_start_time=temp_begin_E(E_pointer)
      cycle_end_time=temp_end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''OUT OF CYCLES 10'')' 
         goto 2098
      end if
*
 2076 IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
         IF (CLLK.EQ.0) GOTO 2098
         GOTO 2076              !loop until read a valid cell event time
      ENDIF
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubmaxBIN) GOTO 2070 !for the next cycle
*
*       ***** look at cell's activity within this group of cycles to see if   *****
*       *****  the relationship/pattern detected with the 1st isubBIN cycles holds true: *****
*
 2098 DO icycle=1,cycle-1       !for each cycle in the group:
         totalBINARY = totalBINARY+1 !keep a total of # of cycles used in BINARY testing
         ITEMPSUM1=0
         ITEMPSUM2=0
         IPT1=IBIGPT-1
         IPT2=IBIGPT+itrt/2-1
         IF (IPT2.GT.itrt) IPT2=IPT2-itrt
         DO J=1,itrt/2
            IPT1=IPT1+1
            IF (IPT1.GT.itrt) THEN
               ITEMPSUM1=ITEMPSUM1+IY(IPT1-itrt,icycle+1)
            ELSE
               ITEMPSUM1=ITEMPSUM1+IY(IPT1,icycle)
            ENDIF
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) THEN
               ITEMPSUM2=ITEMPSUM2+IY(IPT2-itrt,icycle+1)
            ELSE
               ITEMPSUM2=ITEMPSUM2+IY(IPT2,icycle)
            ENDIF
         end do
*
         IF (ITEMPSUM1.GT.ITEMPSUM2) ICNT1=ICNT1+1 !ICNT1=# cycles which support the hypothesis
         IF (ITEMPSUM1.NE.ITEMPSUM2) ICNT2=ICNT2+1 !ICNT2=# cycles in which the cell's activity in
*                                                       ! one half of the cycle does not equal the
*                                                       ! activity in the other half of the cycle
      end do
*
      goto 2045                 
*
 2099 z=(float(ICNT1)-float(ICNT2)/2.)/(.5*SQRT(FLOAT(ICNT2)))
      if(z.gt.ZSIGLVL)then
         zmodsig2_2(IDs(CELL)) = 'R' ! bsig = Z value
      else if((z.le.ZSIGLVL).and.(isubBIN.lt.50))then
         zmodsig2_2(IDs(CELL)) = 'NED' !ned = not enough data for a good test
      else if((z.le.ZSIGLVL).and.(isubBIN.ge.50))then
         zmodsig2_2(IDs(CELL)) = 'N'
      else
      end if
*       print '(''case 2 '',I3,'': ETA2  = '',A,''   ANOVA = '',A,
*     +         ''   BINARY = '',A)',
*     +    CELL,ETA2_2(IDs(CELL)),zmodsig_2(IDs(CELL)),
*     +    zmodsig2_2(IDs(CELL))        
*
*       *********************************************************************
*       *****                                                           *****
*       *****           CASE 3 -- use all control period cycles         *****
*       *****                                                           *****
*       *********************************************************************
*
*
      if((BNDRY_START_TIME.eq.0.0).and. !user has not designated a control period,
     +     (BNDRY_END_TIME.eq.0.0))then ! so cases 3, 4, and 5 do not occur
         ETA2_3 = 'n/a'
         ETA2_4 = 'n/a'
         ETA2_5 = 'n/a'
         zmodsig_3(IDs(CELL)) = 'n/a'
         zmodsig_4(IDs(CELL)) = 'n/a'
         zmodsig_5(IDs(CELL)) = 'n/a'
         zmodsig2_3(IDs(CELL)) = 'n/a'
         zmodsig2_4(IDs(CELL)) = 'n/a'
         zmodsig2_5(IDs(CELL)) = 'na'
         goto 602               !proceed with case 6
      end if
*
*       ***** check to see how many cycles are acceptable:    *****
*       *****  (want min of 101 for most robust BINARY test)  *****
*
      j=0
      temp_begin_E = 0.0
      temp_end_E = 0.0
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
      isubmaxBIN = 0            
*
 305  if(E_pointer.lt.MAX_NUM_EVENTS)then
         cycle_start_time =
     +        SPIKETIMES(E_pointer,IDs(E_pulse))
         cycle_end_time =
     +        SPIKETIMES(E_pointer+1,IDs(E_pulse))
         if(cycle_start_time.lt.BNDRY_START_TIME)then !be sure that both E pulses for
            E_pointer=E_pointer+1 ! this cycle occur within the
            goto 305            ! control period
         end if
*
         if((cycle_start_time.ge.BNDRY_END_TIME).or.
     +        (cycle_end_time.gt.BNDRY_END_TIME))then
            goto 320
         end if
         if(cycle_end_time.eq.0)goto 320 !look at all cycles
*
 306     if(CELL_pointer.le.MAX_NUM_EVENTS)then
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
*
            if (CLLK.lt.cycle_start_time) then
               CELL_pointer=CELL_pointer+1
               goto 306         !loop until find the FIRST valid cell event
! that occurs after the start of the cycle
            endif
*
            if (CLLK.gt.cycle_end_time) then !if the event occurs after the
               E_pointer=E_pointer+1 ! end of the cycle,
               goto 305         ! go to the next cycle
            endif
*
*       ***** at this point, we have detected the first event in the first *****
*       *****  acceptable cycle, so increment the counter:  *****
*
            j=j+1               !keep track of number of acceptable cycles
            temp_begin_E(j)=cycle_start_time !load array with times of acceptable cycles
            temp_end_E(j)=cycle_end_time                
            E_pointer=E_pointer+1 !point to the next respiratory cycle
            goto 305            !get next cycle
         end if
*
      end if
*
 320  if(j.eq.0)then            !cannot do any significance testing for this cell 
         ETA2_3 = 'n/a'         !  because no cycles are acceptable;
         ETA2_4 = 'n/a'         !  enter "not applicable" for these variables
         ETA2_5 = 'n/a'
         zmodsig_3(IDs(CELL)) = 'n/a'
         zmodsig_4(IDs(CELL)) = 'n/a'
         zmodsig_5(IDs(CELL)) = 'n/a'
         zmodsig2_3(IDs(CELL)) = 'n/a'
         zmodsig2_4(IDs(CELL)) = 'n/a'
         zmodsig2_5(IDs(CELL)) = 'n/a'
         goto 602               !proceed to case 6 testing
      end if
*
      if(j.ge.50)then           !50 acceptable cycles are available to calculate
         isubETA=50             ! the ETA**2 value and perform the ANOVA test
      else                      
         isubETA=j              !have to use less than 50 cycles - warn the user
      end if
*
      if(j.ge.101)then          !at least 101 acceptable cycles available to
         isubBIN=50             ! perform the BINARY test - good
         isubmaxBIN=isubBIN + 1 
      else                      !have to use less than 101 cycles - warn the user
         isubBIN=INT(j/2)
         isubmaxBIN=isubBIN + 1 
      end if     
      df1=float(itrt-1)                         
      df2=float(itrt*isubETA-itrt)      
*
*       ***** Calculate the ETA**2 value and the ANOVA statistic (F test): *****
*
*
      ETA2=0.0                  !default value
      PROB=1.0                  !default value
      z=0.0
      RMSB=0
      RMSW=0
      IY = 0                    !IY(itrt,isub) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=1                   ! = respiratory cycle designation
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Get the FIRST acceptable respiratory cycle in the control period: *****
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
         ETA2=0.0
         PROB=1.0
         goto 360
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 300  CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 2'')'
*          ETA2=0.0             !set default values
*           PROB=1.0
            goto 360
         ENDIF
         GOTO 300               !loop until find the FIRST valid cell event
! that occurs after the start of this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 310  FIND_THE_AISLE_3: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 311     if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 3'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 350
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES ETA2 4'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 350         !go to next resp sig test
            endif
*
            GOTO 310            !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 311 !itrt = total # aisles per cycle
*
         ENDIF
      end if FIND_THE_AISLE_3
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubETA)goto 350
      E_pointer=E_pointer+1     !point to next E pulse
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''NOT ENOUGH CYCLES ETA2 5'')' 
*           ETA2=0.0            !set default values
*            PROB=1.0
         goto 350
      end if
*
 340  IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 340
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubETA) GOTO 310 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubETA (optimally 50) cycles,               *****
*       *****  divided each cycle into 20 aisles, and counted the number of         *****
*       *****  cell events that occur within each aisle;                            *****
*       ***** now we're ready for some statistical analysis of respiratory pattern: *****
*
 350  AVGY=float(IYSUM)/float(itrt*isubETA) !total # cell events in the 1st isubETA cycles /
! total # aisles in the 1st isubETA cycles =
!  # cell events expected in each aisle if
! the cell's activity is NOT respiratory modulated
      DO aisle=1,itrt
         RMSB=RMSB+(float(JSUM(aisle))
     +        /float(isubETA)-AVGY)**2
      end do
      RMSB=isubETA*RMSB/df1     !RMSB = between group sum of squares (?)
*
      DO aisle=1,itrt
         DO cycle=1,isubETA
            RMSW=RMSW+(float(IY(aisle,cycle))
     +           -float(JSUM(aisle))/float(isubETA))**2
         end do
      end do
*
      RMSW=RMSW/df2             !RMSW = within group sum of squares (?)
*
      if(RMSW.NE.0) then
         F=RMSB/RMSW
      else
         F=1.0
      end if
*
      if(RMSW.ne.0)then
         ETA2=df1*F/(df1*F+df2) !ETA squared statistic
      else
         ETA2=0.0
      end if
*
      PROB = BETAI(0.5*df2,0.5*df1,df2/(df2+df1*F))
 360  write (ETA2_3(IDs(CELL)),('(F5.2)')) ETA2
c       if(PROB.lt.SIGLVL)zmodsig_3(IDs(CELL)) = 'R'
      if(PROB.lt.SIGLVL)then
         zmodsig_3(IDs(CELL)) = 'R' !cell's activity is respiratory modulated
      else if((PROB.ge.SIGLVL).and.(isubETA.lt.50))then
         zmodsig_3(IDs(CELL)) = 'NED' !ned = not enough data
      else if((PROB.ge.SIGLVL).and.(isubETA.ge.50))then
         zmodsig_3(IDs(CELL)) = 'N' !n = not respiratory modulated
      else
      end if
*
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*         <<<<< PROCEED WITH BINARY TEST for case 3: >>>>>
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*
*
      z=0.0
      ICNT1=0
      ICNT2=0
      ITEMPSUM1=0
      ITEMPSUM2=0
      IY = 0                    !IY(itrt,isub) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=0                   ! = respiratory cycle designation
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Find the first selected respiratory cycle in the recording: *****
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
*           PRINT '(''NOT ENOUGH CYCLES 1'')' 
         goto 3099
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 3000 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES 2'')'
            goto 3099
         ENDIF
         GOTO 3000              !loop until find the FIRST valid cell event
! in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 3010 FIND_THE_AISLE_3_2: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 3011    if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN !ran out of cell events
               goto 3099
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES 4'')' 
               goto 3099        !go to next resp sig test
            endif
*
            GOTO 3010           !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 3011 !itrt = total # aisles per cycle
*
         end if
      END IF FIND_THE_AISLE_3_2
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubBIN)goto 3041 !proceed to BINARY analysis
      E_pointer=E_pointer+1     !point to next E pulse
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or. !check for legitimate E pulses
     +     (cycle_end_time.eq.0.0))then         
*           PRINT '(''NOT ENOUGH CYCLES 5'')' 
         goto 3099
      end if
*
 3040 IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 3040
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubBIN) GOTO 3010 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubBIN (optimally 50) cycles,                   *****
*       *****  divided each one into 20 aisles,                                         *****
*       *****  and counted the number of cell events that occur within each aisle;      *****
*
*       ***** make 2 groups of "contiguous" aisles                      *****
*       ***** (ie., 1-10 & 11-20; 2-11 & 12-20 + 1; etc.);              *****
*       ***** Examine all possible groupings until find the case in     *****
*       *****  which the difference in total # cell events between the  *****
*       *****  two groups is maximal:                                   *****
*
 3041 IPT1=(itrt/2) - 1         ! = # of first aisle in the 1st group
      IPT2=itrt - 1             ! = # of first aisle in the 2nd group
      IOLDDIFF=0
      IBIGPT=1
      DO I=1,itrt/2             !update pointers to first aisle in each group
         IPT1=IPT1 - itrt/2 + 1
         IPT2=IPT2 - itrt/2 + 1
         IF (IPT2.LE.0) IPT2=IPT2 + itrt
         ITEMPSUM1=0            !ITEMPSUM1 = total # cell events in 1st group
         ITEMPSUM2=0            !ITEMPSUM2 = total # cell events in 2nd group
*
         DO J=1,itrt/2
            IPT1=IPT1+1         !point to next aisle in each group
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) IPT2=IPT2-itrt
            ITEMPSUM1=ITEMPSUM1+JSUM(IPT1)
            ITEMPSUM2=ITEMPSUM2+JSUM(IPT2)
         end do
*
         IF (ITEMPSUM1.GT.ITEMPSUM2) THEN !compare the #s of events which occurred in the
            INEWDIFF=ITEMPSUM1-ITEMPSUM2 ! 1st group vs. those in the 2nd group
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT1 - itrt/2 + 1 !store # of 1st aisle in group 1
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ELSE
            INEWDIFF=ITEMPSUM2-ITEMPSUM1 !ITEMPSUM2 > ITEMPSUM1
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT2 - itrt/2 + 1 !store # of 1st aisle in group 2
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ENDIF
*
      end do
*
      totalBINARY=isubBIN       !keep track of total # cycles used in BINARY testing
*
*       ***** Now have formed a hypothesis and gathered data based on the first isubBIN cycles *****
*       *****  about whether or not the cell's activity is greater in one half of the          *****
*       *****  cycle than in the other half.  Now need to look at remaining                       *****
*       *****  cycles to see if the relationship/pattern holds true:                           *****
*
*
      ICNT1=0                   != # times ITEMPSUM1 > ITEMPSUM2
      ICNT2=0                   != # times ITEMPSUM1 <> ITEMPSUM2
 3045 IY = 0
      cycle=0
      aisle=1
*
*      ***** move to the next set of E pulses (hence, the next cycle): *****
*      
      E_pointer=E_pointer+1
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*       PRINT '(''OUT OF CYCLES 6'')' 
         goto 3099
      end if
*
*      ***** find the first cell event after the first cycle start: *****
*
 3060 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN !ran out of cell events
*         PRINT '(''OUT OF CYCLES 7'')'
            GOTO 3099
         ENDIF
         GOTO 3060              !loop until find 1st cell event in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 3070 FIND_THE_AISLE_3_3: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
         if (CLLK.LT.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/itrt))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts events in array subscripted
*                                                       by aisle # and cycle #
            CELL_pointer=CELL_pointer+1
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*         PRINT '(''OUT OF CYCLES 8'')'
               GOTO 3098
            ENDIF
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event
*
            IF (CLLK.EQ.0) THEN
*         PRINT '(''OUT OF CYCLES 9'')'
               GOTO 3098
            ENDIF
*
            GOTO 3070           !loop back to analyze the new event time
*
         ELSE                   !cell event does not occur in current aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 3070 !look at the next aisle (or go to the next cycle)
*
         end if
      END IF FIND_THE_AISLE_3_3
*
      if(cycle.eq.isubmaxBIN)goto 3098
      E_pointer=E_pointer+1     !proceed to the next resp cycle
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''OUT OF CYCLES 10'')' 
         goto 3098
      end if
*
 3076 IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
         IF (CLLK.EQ.0) GOTO 3098
         GOTO 3076              !loop until read a valid cell event time
      ENDIF
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubmaxBIN) GOTO 3070 !for the next cycle
*
*       ***** look at cell's activity within this group of cycles to see if   *****
*       *****  the relationship/pattern detected with the 1st isubBIN cycles holds true: *****
*
 3098 DO icycle=1,cycle-1       !for each cycle in the group:
         totalBINARY = totalBINARY+1 !keep a total of # of cycles used in BINARY testing
         ITEMPSUM1=0
         ITEMPSUM2=0
         IPT1=IBIGPT-1
         IPT2=IBIGPT+itrt/2-1
         IF (IPT2.GT.itrt) IPT2=IPT2-itrt
         DO J=1,itrt/2
            IPT1=IPT1+1
            IF (IPT1.GT.itrt) THEN
               ITEMPSUM1=ITEMPSUM1+IY(IPT1-itrt,icycle+1)
            ELSE
               ITEMPSUM1=ITEMPSUM1+IY(IPT1,icycle)
            ENDIF
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) THEN
               ITEMPSUM2=ITEMPSUM2+IY(IPT2-itrt,icycle+1)
            ELSE
               ITEMPSUM2=ITEMPSUM2+IY(IPT2,icycle)
            ENDIF
         end do
*
         IF (ITEMPSUM1.GT.ITEMPSUM2) ICNT1=ICNT1+1 !ICNT1=# cycles which support the hypothesis
         IF (ITEMPSUM1.NE.ITEMPSUM2) ICNT2=ICNT2+1 !ICNT2=# cycles in which the cell's activity in
*                                                       ! one half of the cycle does not equal the
*                                                       ! activity in the other half of the cycle
      end do
*
      goto 3045                 
*
 3099 z=(float(ICNT1)-float(ICNT2)/2.)/(.5*SQRT(FLOAT(ICNT2)))
      if(z.gt.ZSIGLVL)then
         zmodsig2_3(IDs(CELL)) = 'R' ! bsig = Z value
      else if((z.le.ZSIGLVL).and.(isubBIN.lt.50))then
         zmodsig2_3(IDs(CELL)) = 'NED' !ned = not enough data for a good test
      else if((z.le.ZSIGLVL).and.(isubBIN.ge.50))then
         zmodsig2_3(IDs(CELL)) = 'N'
      else
      end if
*       print '(''case 3 '',I3,'': ETA2  = '',A,''   ANOVA = '',A1,
*     +         ''   BINARY = '',A1)',
*     +    CELL,ETA2_3(IDs(CELL)),zmodsig_3(IDs(CELL)),
*     +    zmodsig2_3(IDs(CELL))        
*
*       *********************************************************************
*       *****                                                           *****
*       *****   CASE 4 -- reject control cycles > control_max_value     *****
*       *****                   (= mean + 2 sd of all control cycles)   *****
*       *****                                                           *****
*       *********************************************************************
*
*
*
*       ***** check to see how many cycles are acceptable:    *****
*       *****  (want min of 101 for most robust BINARY test)  *****
*
      j=0
      temp_begin_E = 0.0
      temp_end_E = 0.0
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
      isubmaxBIN = 0            
*
 405  if(E_pointer.lt.MAX_NUM_EVENTS)then
         cycle_start_time =
     +        SPIKETIMES(E_pointer,IDs(E_pulse))
         cycle_end_time =
     +        SPIKETIMES(E_pointer+1,IDs(E_pulse))
         if(cycle_start_time.lt.BNDRY_START_TIME)then !be sure that both E pulses for
            E_pointer=E_pointer+1 ! this cycle occur within the
            goto 405            ! control period
         end if
*
         if((cycle_start_time.ge.BNDRY_END_TIME).or.
     +        (cycle_end_time.gt.BNDRY_END_TIME))then
            goto 420
         end if
*
         if(cycle_end_time.eq.0)goto 420 !look at all cycles
*
         if(cycle_end_time-cycle_start_time.gt.
     +        control_max_value) then
            E_pointer=E_pointer+1 !this cycle too long, so get next one
            goto 405
         end if
*
 406     if(CELL_pointer.le.MAX_NUM_EVENTS)then
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
*
            if (CLLK.lt.cycle_start_time) then
               CELL_pointer=CELL_pointer+1
               goto 406         !loop until find the FIRST valid cell event
! that occurs after the start of the cycle
            endif
*
            if (CLLK.gt.cycle_end_time) then !if the event occurs after the
               E_pointer=E_pointer+1 ! end of the cycle,
               goto 405         ! go to the next cycle
            endif
*
*       ***** at this point, we have detected the first event in the first *****
*       *****  acceptable cycle, so increment the counter:  *****
*
            j=j+1               !keep track of number of acceptable cycles
            temp_begin_E(j)=cycle_start_time !load array with times of acceptable cycles
            temp_end_E(j)=cycle_end_time                
            E_pointer=E_pointer+1 !point to the next respiratory cycle
            goto 405            !get next cycle
         end if
*
      end if
*
 420  if(j.eq.0)then            !cannot do any significance testing for this cell 
!  in this case because no cycles are acceptable;
         ETA2_4 = 'n/a'         !  enter "not applicable" for these variables
         zmodsig_4(IDs(CELL)) = 'n/a'
         zmodsig2_4(IDs(CELL)) = 'n/a'
         goto 502               !proceed to case 5 testing
      end if
*
      if(j.ge.50)then           !50 acceptable cycles are available to calculate
         isubETA=50             ! the ETA**2 value and perform the ANOVA test
      else                      
         isubETA=j              !have to use less than 50 cycles - warn the user
      end if
*
      if(j.ge.101)then          !at least 101 acceptable cycles available to
         isubBIN=50             ! perform the BINARY test - good
         isubmaxBIN=isubBIN + 1 
      else                      !have to use less than 101 cycles - warn the user
         isubBIN=INT(j/2)
         isubmaxBIN=isubBIN + 1 
      end if     
      df1=float(itrt-1)                         
      df2=float(itrt*isubETA-itrt)      
*
*       ***** Calculate the ETA**2 value and the ANOVA statistic (F test): *****
*
*
      ETA2=0.0                  !default value
      PROB=1.0                  !default value
      z=0.0
      RMSB=0
      RMSW=0
      IY = 0                    !IY(itrt,isub) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=1                   ! = respiratory cycle designation
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Get the FIRST acceptable respiratory cycle in the control period: *****
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
         ETA2=0.0
         PROB=1.0
         goto 460
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 400  CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 2'')'
*          ETA2=0.0             !set default values
*           PROB=1.0
            goto 460
         ENDIF
         GOTO 400               !loop until find the FIRST valid cell event
! that occurs after the start of this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 410  FIND_THE_AISLE_4: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 411     if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 3'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 450
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES ETA2 4'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 450         !go to next resp sig test
            endif
*
            GOTO 410            !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 411 !itrt = total # aisles per cycle
*
         ENDIF
      end if FIND_THE_AISLE_4
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubETA)goto 450
      E_pointer=E_pointer+1     !point to next E pulse
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''NOT ENOUGH CYCLES ETA2 5'')' 
*           ETA2=0.0            !set default values
*            PROB=1.0
         goto 450
      end if
*
 440  IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 440
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubETA) GOTO 410 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubETA (optimally 50) cycles,               *****
*       *****  divided each cycle into 20 aisles, and counted the number of         *****
*       *****  cell events that occur within each aisle;                            *****
*       ***** now we're ready for some statistical analysis of respiratory pattern: *****
*
 450  AVGY=float(IYSUM)/float(itrt*isubETA) !total # cell events in the 1st isubETA cycles /
! total # aisles in the 1st isubETA cycles =
!  # cell events expected in each aisle if
! the cell's activity is NOT respiratory modulated
      DO aisle=1,itrt
         RMSB=RMSB+(float(JSUM(aisle))
     +        /float(isubETA)-AVGY)**2
      end do
      RMSB=isubETA*RMSB/df1     !RMSB = between group sum of squares (?)
*
      DO aisle=1,itrt
         DO cycle=1,isubETA
            RMSW=RMSW+(float(IY(aisle,cycle))
     +           -float(JSUM(aisle))/float(isubETA))**2
         end do
      end do
*
      RMSW=RMSW/df2             !RMSW = within group sum of squares (?)
*
      if(RMSW.NE.0) then
         F=RMSB/RMSW
      else
         F=1.0
      end if
*
      if(RMSW.ne.0)then
         ETA2=df1*F/(df1*F+df2) !ETA squared statistic
      else
         ETA2=0.0
      end if
*
      PROB = BETAI(0.5*df2,0.5*df1,df2/(df2+df1*F))
 460  write (ETA2_4(IDs(CELL)),('(F5.2)')) ETA2
c       if(PROB.lt.SIGLVL)zmodsig_4(IDs(CELL)) = 'R'
      if(PROB.lt.SIGLVL)then
         zmodsig_4(IDs(CELL)) = 'R' !cell's activity is respiratory modulated
      else if((PROB.ge.SIGLVL).and.(isubETA.lt.50))then
         zmodsig_4(IDs(CELL)) = 'NED' !ned = not enough data
      else if((PROB.ge.SIGLVL).and.(isubETA.ge.50))then
         zmodsig_4(IDs(CELL)) = 'N' !n = not respiratory modulated
      else
      end if
*
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*         <<<<< PROCEED WITH BINARY TEST for case 4: >>>>>
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*
*
      z=0.0
      ICNT1=0
      ICNT2=0
      ITEMPSUM1=0
      ITEMPSUM2=0
      IY = 0                    !IY(itrt,isub) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=0                   ! = respiratory cycle designation
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Find the first selected respiratory cycle in the recording: *****
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
*           PRINT '(''NOT ENOUGH CYCLES 1'')' 
         goto 4099
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 4000 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES 2'')'
            goto 4099
         ENDIF
         GOTO 4000              !loop until find the FIRST valid cell event
! in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 4010 FIND_THE_AISLE_4_2: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 4011    if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN !ran out of cell events
               goto 4099
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES 4'')' 
               goto 4099        !go to next resp sig test
            endif
*
            GOTO 4010           !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 4011 !itrt = total # aisles per cycle
*
         end if
      END IF FIND_THE_AISLE_4_2
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubBIN)goto 4041 !proceed to BINARY analysis
      E_pointer=E_pointer+1     !point to next E pulse
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or. !check for legitimate E pulses
     +     (cycle_end_time.eq.0.0))then         
*           PRINT '(''NOT ENOUGH CYCLES 5'')' 
         goto 4099
      end if
*
 4040 IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 4040
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubBIN) GOTO 4010 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubBIN (optimally 50) cycles,                   *****
*       *****  divided each one into 20 aisles,                                         *****
*       *****  and counted the number of cell events that occur within each aisle;      *****
*
*       ***** make 2 groups of "contiguous" aisles                      *****
*       ***** (ie., 1-10 & 11-20; 2-11 & 12-20 + 1; etc.);              *****
*       ***** Examine all possible groupings until find the case in     *****
*       *****  which the difference in total # cell events between the  *****
*       *****  two groups is maximal:                                   *****
*
 4041 IPT1=(itrt/2) - 1         ! = # of first aisle in the 1st group
      IPT2=itrt - 1             ! = # of first aisle in the 2nd group
      IOLDDIFF=0
      IBIGPT=1
      DO I=1,itrt/2             !update pointers to first aisle in each group
         IPT1=IPT1 - itrt/2 + 1
         IPT2=IPT2 - itrt/2 + 1
         IF (IPT2.LE.0) IPT2=IPT2 + itrt
         ITEMPSUM1=0            !ITEMPSUM1 = total # cell events in 1st group
         ITEMPSUM2=0            !ITEMPSUM2 = total # cell events in 2nd group
*
         DO J=1,itrt/2
            IPT1=IPT1+1         !point to next aisle in each group
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) IPT2=IPT2-itrt
            ITEMPSUM1=ITEMPSUM1+JSUM(IPT1)
            ITEMPSUM2=ITEMPSUM2+JSUM(IPT2)
         end do
*
         IF (ITEMPSUM1.GT.ITEMPSUM2) THEN !compare the #s of events which occurred in the
            INEWDIFF=ITEMPSUM1-ITEMPSUM2 ! 1st group vs. those in the 2nd group
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT1 - itrt/2 + 1 !store # of 1st aisle in group 1
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ELSE
            INEWDIFF=ITEMPSUM2-ITEMPSUM1 !ITEMPSUM2 > ITEMPSUM1
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT2 - itrt/2 + 1 !store # of 1st aisle in group 2
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ENDIF
*
      end do
*
      totalBINARY=isubBIN       !keep track of total # cycles used in BINARY testing
*
*       ***** Now have formed a hypothesis and gathered data based on the first isubBIN cycles *****
*       *****  about whether or not the cell's activity is greater in one half of the          *****
*       *****  cycle than in the other half.  Now need to look at remaining                       *****
*       *****  cycles to see if the relationship/pattern holds true:                           *****
*
*
      ICNT1=0                   != # times ITEMPSUM1 > ITEMPSUM2
      ICNT2=0                   != # times ITEMPSUM1 <> ITEMPSUM2
 4045 IY = 0
      cycle=0
      aisle=1
*
*      ***** move to the next set of E pulses (hence, the next cycle): *****
*      
      E_pointer=E_pointer+1
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*       PRINT '(''OUT OF CYCLES 6'')' 
         goto 4099
      end if
*
*      ***** find the first cell event after the first cycle start: *****
*
 4060 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN !ran out of cell events
*         PRINT '(''OUT OF CYCLES 7'')'
            GOTO 4099
         ENDIF
         GOTO 4060              !loop until find 1st cell event in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 4070 FIND_THE_AISLE_4_3: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
         if (CLLK.LT.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/itrt))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts events in array subscripted
*                                                       by aisle # and cycle #
            CELL_pointer=CELL_pointer+1
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*         PRINT '(''OUT OF CYCLES 8'')'
               GOTO 4098
            ENDIF
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event
*
            IF (CLLK.EQ.0) THEN
*         PRINT '(''OUT OF CYCLES 9'')'
               GOTO 4098
            ENDIF
*
            GOTO 4070           !loop back to analyze the new event time
*
         ELSE                   !cell event does not occur in current aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 4070 !look at the next aisle (or go to the next cycle)
*
         end if
      END IF FIND_THE_AISLE_4_3
*
      if(cycle.eq.isubmaxBIN)goto 4098
      E_pointer=E_pointer+1     !proceed to the next resp cycle
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''OUT OF CYCLES 10'')' 
         goto 4098
      end if
*
 4076 IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
         IF (CLLK.EQ.0) GOTO 4098
         GOTO 4076              !loop until read a valid cell event time
      ENDIF
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubmaxBIN) GOTO 4070 !for the next cycle
*
*       ***** look at cell's activity within this group of  cycles to see if   *****
*       *****  the relationship/pattern detected with the 1st isubBIN cycles holds true: *****
*
 4098 DO icycle=1,cycle-1       !for each cycle in the group:
         totalBINARY = totalBINARY+1 !keep a total of # of cycles used in BINARY testing
         ITEMPSUM1=0
         ITEMPSUM2=0
         IPT1=IBIGPT-1
         IPT2=IBIGPT+itrt/2-1
         IF (IPT2.GT.itrt) IPT2=IPT2-itrt
         DO J=1,itrt/2
            IPT1=IPT1+1
            IF (IPT1.GT.itrt) THEN
               ITEMPSUM1=ITEMPSUM1+IY(IPT1-itrt,icycle+1)
            ELSE
               ITEMPSUM1=ITEMPSUM1+IY(IPT1,icycle)
            ENDIF
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) THEN
               ITEMPSUM2=ITEMPSUM2+IY(IPT2-itrt,icycle+1)
            ELSE
               ITEMPSUM2=ITEMPSUM2+IY(IPT2,icycle)
            ENDIF
         end do
*
         IF (ITEMPSUM1.GT.ITEMPSUM2) ICNT1=ICNT1+1 !ICNT1=# cycles which support the hypothesis
         IF (ITEMPSUM1.NE.ITEMPSUM2) ICNT2=ICNT2+1 !ICNT2=# cycles in which the cell's activity in
*                                                       ! one half of the cycle does not equal the
*                                                       ! activity in the other half of the cycle
      end do
*
      goto 4045                 
*
 4099 z=(float(ICNT1)-float(ICNT2)/2.)/(.5*SQRT(FLOAT(ICNT2)))
      if(z.gt.ZSIGLVL)then
         zmodsig2_4(IDs(CELL)) = 'R' ! bsig = Z value
      else if((z.le.ZSIGLVL).and.(isubBIN.lt.50))then
         zmodsig2_4(IDs(CELL)) = 'NED' !ned = not enough data for a good test
      else if((z.le.ZSIGLVL).and.(isubBIN.ge.50))then
         zmodsig2_4(IDs(CELL)) = 'N'
      else
      end if
*       print '(''case 4 '',I3,'': ETA2  = '',A,''   ANOVA = '',A1,
*     +         ''   BINARY = '',A1)',
*     +    CELL,ETA2_4(IDs(CELL)),zmodsig_4(IDs(CELL)),
*     +    zmodsig2_4(IDs(CELL))        
*
*       *********************************************************************
*       *****                                                           *****
*       *****   CASE 5 -- reject control cycles > user_max_value        *****
*       *****                   (default user_max_value=10,000.0 ms     *****
*       *****                                                           *****
*       *********************************************************************
*
*
c502    if(user_max_value.eq.1.0)then           !user does not want this feature
c           ETA2_5 = 'n/a'                      !  enter "not applicable" for these variables
c           zmodsig_5(IDs(CELL)) = 'n/a'
c           zmodsig2_5(IDs(CELL)) = 'n/a'
c           goto 602                            !proceed to case 6
c         end if
*
*       ***** check to see how many cycles are acceptable:    *****
*       *****  (want min of 101 for most robust BINARY test)  *****
*
 502  j=0
      temp_begin_E = 0.0
      temp_end_E = 0.0
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
      isubmaxBIN = 0            
*
 505  if(E_pointer.lt.MAX_NUM_EVENTS)then
         cycle_start_time =
     +        SPIKETIMES(E_pointer,IDs(E_pulse))
         cycle_end_time =
     +        SPIKETIMES(E_pointer+1,IDs(E_pulse))
         if(cycle_start_time.lt.BNDRY_START_TIME)then !be sure that both E pulses for
            E_pointer=E_pointer+1 ! this cycle occur within the
            goto 505            ! control period
         end if
*
         if((cycle_start_time.ge.BNDRY_END_TIME).or.
     +        (cycle_end_time.gt.BNDRY_END_TIME))then
            goto 520
         end if
*
         if(cycle_end_time.eq.0)goto 520 !look at all cycles
*
         if(cycle_end_time-cycle_start_time.gt.
     +        user_max_value) then
            E_pointer=E_pointer+1 !this cycle too long, so get next one
            goto 505
         end if
*
 506     if(CELL_pointer.le.MAX_NUM_EVENTS)then
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
*
            if (CLLK.lt.cycle_start_time) then
               CELL_pointer=CELL_pointer+1
               goto 506         !loop until find the FIRST valid cell event
! that occurs after the start of the cycle
            endif
*
            if (CLLK.gt.cycle_end_time) then !if the event occurs after the
               E_pointer=E_pointer+1 ! end of the cycle,
               goto 505         ! go to the next cycle
            endif
*
*       ***** at this point, we have detected the first event in the first *****
*       *****  acceptable cycle, so increment the counter:  *****
*
            j=j+1               !keep track of number of acceptable cycles
            temp_begin_E(j)=cycle_start_time !load array with times of acceptable cycles
            temp_end_E(j)=cycle_end_time                
            E_pointer=E_pointer+1 !point to the next respiratory cycle
            goto 505            !get next cycle
         end if
*
      end if
*
 520  if(j.eq.0)then            !cannot do any significance testing for this cell 
!  in this case because no cycles are acceptable;
         ETA2_5 = 'n/a'         !  enter "not applicable" for these variables
         zmodsig_5(IDs(CELL)) = 'n/a'
         zmodsig2_5(IDs(CELL)) = 'n/a'
         goto 602               !proceed to case 6 testing
      end if
*
      if(j.ge.50)then           !50 acceptable cycles are available to calculate
         isubETA=50             ! the ETA**2 value and perform the ANOVA test
      else                      
         isubETA=j              !have to use less than 50 cycles - warn the user
      end if
*
      if(j.ge.101)then          !at least 101 acceptable cycles available to
         isubBIN=50             ! perform the BINARY test - good
         isubmaxBIN=isubBIN + 1 
      else                      !have to use less than 101 cycles - warn the user
         isubBIN=INT(j/2)
         isubmaxBIN=isubBIN + 1 
      end if     
      df1=float(itrt-1)                         
      df2=float(itrt*isubETA-itrt)      
*
*       ***** Calculate the ETA**2 value and the ANOVA statistic (F test): *****
*
*
      ETA2=0.0                  !default value
      PROB=1.0                  !default value
      z=0.0
      RMSB=0
      RMSW=0
      IY = 0                    !IY(itrt,isub) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=1                   ! = respiratory cycle designation
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Get the FIRST acceptable respiratory cycle in the control period: *****
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
         ETA2=0.0
         PROB=1.0
         goto 560
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 500  CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 2'')'
*          ETA2=0.0             !set default values
*           PROB=1.0
            goto 560
         ENDIF
         GOTO 500               !loop until find the FIRST valid cell event
! that occurs after the start of this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 510  FIND_THE_AISLE_5: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 511     if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 3'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 550
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES ETA2 5'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 550         !go to next resp sig test
            endif
*
            GOTO 510            !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 511 !itrt = total # aisles per cycle
*
         ENDIF
      end if FIND_THE_AISLE_5
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubETA)goto 550
      E_pointer=E_pointer+1     !point to next E pulse
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''NOT ENOUGH CYCLES ETA2 5'')' 
*           ETA2=0.0            !set default values
*            PROB=1.0
         goto 550
      end if
*
 540  IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 540
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubETA) GOTO 510 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubETA (optimally 50) cycles,               *****
*       *****  divided each cycle into 20 aisles, and counted the number of         *****
*       *****  cell events that occur within each aisle;                            *****
*       ***** now we're ready for some statistical analysis of respiratory pattern: *****
*
 550  AVGY=float(IYSUM)/float(itrt*isubETA) !total # cell events in the 1st isubETA cycles /
! total # aisles in the 1st isubETA cycles =
!  # cell events expected in each aisle if
! the cell's activity is NOT respiratory modulated
      DO aisle=1,itrt
         RMSB=RMSB+(float(JSUM(aisle))
     +        /float(isubETA)-AVGY)**2
      end do
      RMSB=isubETA*RMSB/df1     !RMSB = between group sum of squares (?)
*
      DO aisle=1,itrt
         DO cycle=1,isubETA
            RMSW=RMSW+(float(IY(aisle,cycle))
     +           -float(JSUM(aisle))/float(isubETA))**2
         end do
      end do
*
      RMSW=RMSW/df2             !RMSW = within group sum of squares (?)
*
      if(RMSW.NE.0) then
         F=RMSB/RMSW
      else
         F=1.0
      end if
*
      if(RMSW.ne.0)then
         ETA2=df1*F/(df1*F+df2) !ETA squared statistic
      else
         ETA2=0.0
      end if
*
      PROB = BETAI(0.5*df2,0.5*df1,df2/(df2+df1*F))
 560  write (ETA2_5(IDs(CELL)),('(F5.2)')) ETA2
c       if(PROB.lt.SIGLVL)zmodsig_5(IDs(CELL)) = 'R'
      if(PROB.lt.SIGLVL)then
         zmodsig_5(IDs(CELL)) = 'R' !cell's activity is respiratory modulated
      else if((PROB.ge.SIGLVL).and.(isubETA.lt.50))then
         zmodsig_5(IDs(CELL)) = 'NED' !ned = not enough data
      else if((PROB.ge.SIGLVL).and.(isubETA.ge.50))then
         zmodsig_5(IDs(CELL)) = 'N' !n = not respiratory modulated
      else
      end if
*
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*         <<<<< PROCEED WITH BINARY TEST for case 5: >>>>>
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*
*
      z=0.0
      ICNT1=0
      ICNT2=0
      ITEMPSUM1=0
      ITEMPSUM2=0
      IY = 0                    !IY(itrt,isub) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=0                   ! = respiratory cycle designation
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Find the first selected respiratory cycle in the recording: *****
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
*           PRINT '(''NOT ENOUGH CYCLES 1'')' 
         goto 5099
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 5000 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES 2'')'
            goto 5099
         ENDIF
         GOTO 5000              !loop until find the FIRST valid cell event
! in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 5010 FIND_THE_AISLE_5_2: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 5011    if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN !ran out of cell events
               goto 5099
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES 4'')' 
               goto 5099        !go to next resp sig test
            endif
*
            GOTO 5010           !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 5011 !itrt = total # aisles per cycle
*
         end if
      END IF FIND_THE_AISLE_5_2
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubBIN)goto 5041 !proceed to BINARY analysis
      E_pointer=E_pointer+1     !point to next E pulse
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or. !check for legitimate E pulses
     +     (cycle_end_time.eq.0.0))then         
*           PRINT '(''NOT ENOUGH CYCLES 5'')' 
         goto 5099
      end if
*
 5040 IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 5040
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubBIN) GOTO 5010 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubBIN (optimally 50) cycles,                   *****
*       *****  divided each one into 20 aisles,                                         *****
*       *****  and counted the number of cell events that occur within each aisle;      *****
*
*       ***** make 2 groups of "contiguous" aisles                      *****
*       ***** (ie., 1-10 & 11-20; 2-11 & 12-20 + 1; etc.);              *****
*       ***** Examine all possible groupings until find the case in     *****
*       *****  which the difference in total # cell events between the  *****
*       *****  two groups is maximal:                                   *****
*
 5041 IPT1=(itrt/2) - 1         ! = # of first aisle in the 1st group
      IPT2=itrt - 1             ! = # of first aisle in the 2nd group
      IOLDDIFF=0
      IBIGPT=1
      DO I=1,itrt/2             !update pointers to first aisle in each group
         IPT1=IPT1 - itrt/2 + 1
         IPT2=IPT2 - itrt/2 + 1
         IF (IPT2.LE.0) IPT2=IPT2 + itrt
         ITEMPSUM1=0            !ITEMPSUM1 = total # cell events in 1st group
         ITEMPSUM2=0            !ITEMPSUM2 = total # cell events in 2nd group
*
         DO J=1,itrt/2
            IPT1=IPT1+1         !point to next aisle in each group
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) IPT2=IPT2-itrt
            ITEMPSUM1=ITEMPSUM1+JSUM(IPT1)
            ITEMPSUM2=ITEMPSUM2+JSUM(IPT2)
         end do
*
         IF (ITEMPSUM1.GT.ITEMPSUM2) THEN !compare the #s of events which occurred in the
            INEWDIFF=ITEMPSUM1-ITEMPSUM2 ! 1st group vs. those in the 2nd group
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT1 - itrt/2 + 1 !store # of 1st aisle in group 1
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ELSE
            INEWDIFF=ITEMPSUM2-ITEMPSUM1 !ITEMPSUM2 > ITEMPSUM1
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT2 - itrt/2 + 1 !store # of 1st aisle in group 2
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ENDIF
*
      end do
*
      totalBINARY=isubBIN       !keep track of total # cycles used in BINARY testing
*
*       ***** Now have formed a hypothesis and gathered data based on the first isubBIN cycles *****
*       *****  about whether or not the cell's activity is greater in one half of the          *****
*       *****  cycle than in the other half.  Now need to look at remaining                       *****
*       *****  cycles to see if the relationship/pattern holds true:                           *****
*
*
      ICNT1=0                   != # times ITEMPSUM1 > ITEMPSUM2
      ICNT2=0                   != # times ITEMPSUM1 <> ITEMPSUM2
 5045 IY = 0
      cycle=0
      aisle=1
*
*      ***** move to the next set of E pulses (hence, the next cycle): *****
*      
      E_pointer=E_pointer+1
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*       PRINT '(''OUT OF CYCLES 6'')' 
         goto 5099
      end if
*
*      ***** find the first cell event after the first cycle start: *****
*
 5060 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN !ran out of cell events
*         PRINT '(''OUT OF CYCLES 7'')'
            GOTO 5099
         ENDIF
         GOTO 5060              !loop until find 1st cell event in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 5070 FIND_THE_AISLE_5_3: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
         if (CLLK.LT.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/itrt))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts events in array subscripted
*                                                       by aisle # and cycle #
            CELL_pointer=CELL_pointer+1
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*         PRINT '(''OUT OF CYCLES 8'')'
               GOTO 5098
            ENDIF
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event
*
            IF (CLLK.EQ.0) THEN
*         PRINT '(''OUT OF CYCLES 9'')'
               GOTO 5098
            ENDIF
*
            GOTO 5070           !loop back to analyze the new event time
*
         ELSE                   !cell event does not occur in current aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 5070 !look at the next aisle (or go to the next cycle)
*
         end if
      END IF FIND_THE_AISLE_5_3
*
      if(cycle.eq.isubmaxBIN)goto 5098
      E_pointer=E_pointer+1     !proceed to the next resp cycle
*
      cycle_start_time = temp_begin_E(E_pointer)
      cycle_end_time = temp_end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''OUT OF CYCLES 10'')' 
         goto 5098
      end if
*
 5076 IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
         IF (CLLK.EQ.0) GOTO 5098
         GOTO 5076              !loop until read a valid cell event time
      ENDIF
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubmaxBIN) GOTO 5070 !for the next cycle
*
*       ***** look at cell's activity within this group of  cycles to see if   *****
*       *****  the relationship/pattern detected with the 1st isubBIN cycles holds true: *****
*
 5098 DO icycle=1,cycle-1       !for each cycle in the group:
         totalBINARY = totalBINARY+1 !keep a total of # of cycles used in BINARY testing
         ITEMPSUM1=0
         ITEMPSUM2=0
         IPT1=IBIGPT-1
         IPT2=IBIGPT+itrt/2-1
         IF (IPT2.GT.itrt) IPT2=IPT2-itrt
         DO J=1,itrt/2
            IPT1=IPT1+1
            IF (IPT1.GT.itrt) THEN
               ITEMPSUM1=ITEMPSUM1+IY(IPT1-itrt,icycle+1)
            ELSE
               ITEMPSUM1=ITEMPSUM1+IY(IPT1,icycle)
            ENDIF
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) THEN
               ITEMPSUM2=ITEMPSUM2+IY(IPT2-itrt,icycle+1)
            ELSE
               ITEMPSUM2=ITEMPSUM2+IY(IPT2,icycle)
            ENDIF
         end do
*
         IF (ITEMPSUM1.GT.ITEMPSUM2) ICNT1=ICNT1+1 !ICNT1=# cycles which support the hypothesis
         IF (ITEMPSUM1.NE.ITEMPSUM2) ICNT2=ICNT2+1 !ICNT2=# cycles in which the cell's activity in
*                                                       ! one half of the cycle does not equal the
*                                                       ! activity in the other half of the cycle
      end do
*
      goto 5045                 
*
 5099 z=(float(ICNT1)-float(ICNT2)/2.)/(.5*SQRT(FLOAT(ICNT2)))
      if(z.gt.ZSIGLVL)then
         zmodsig2_5(IDs(CELL)) = 'R' ! bsig = Z value
      else if((z.le.ZSIGLVL).and.(isubBIN.lt.50))then
         zmodsig2_5(IDs(CELL)) = 'NED' !ned = not enough data for a good test
      else if((z.le.ZSIGLVL).and.(isubBIN.ge.50))then
         zmodsig2_5(IDs(CELL)) = 'N'
      else
      end if
*       print '(''case 5 '',I3,'': ETA2  = '',A,''   ANOVA = '',A,
*     +         ''   BINARY = '',A)',
*     +    CELL,ETA2_5(IDs(CELL)),zmodsig_5(IDs(CELL)),
*     +    zmodsig2_5(IDs(CELL))        
*
*       *********************************************************************
*       *****                                                           *****
*       *****   CASE 6 -- use same control cycles as in CTH             *****
*       *****                                                           *****
*       *********************************************************************
*
*       *********************************************************************
*       *****                                                           *****
*       *****   calculate the ETA**2 and ANOVA statistic using the      *****
*       *****   first 50 acceptable control cycles of cell activity:    *****
*       *****                                                           *****
*       ***** (acceptable cycles are those during which the cell fires) *****
*       *****   If 50 acceptable cycles are not available, the user is  *****
*       *****   notified and calculation proceeds using available cycles*****
*       *****                                                           *****
*       *********************************************************************
*
*                       
*       ***** check to see how many cycles are acceptable:    *****
*       *****  (want min of 101 for most robust BINARY test)  *****
*
 602  j=0
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
      isubmaxBIN = 0            
*
 605  if(E_pointer.le.max_cycles)then
         cycle_start_time=begin_E(E_pointer)
         cycle_end_time=end_E(E_pointer)
         if(cycle_end_time.eq.0)goto 620 !look at all cycles
*
 606     if(CELL_pointer.le.MAX_NUM_EVENTS)then
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
*
            if (CLLK.lt.cycle_start_time) then
               CELL_pointer=CELL_pointer+1
               goto 606         !loop until find the FIRST valid cell event
! that occurs after the start of the cycle
            endif
*
            if (CLLK.gt.cycle_end_time) then !if the event occurs after the
               E_pointer=E_pointer+1 ! end of the cycle,
               goto 605         ! go to the next cycle
            endif
*
*       ***** at this point, we have detected the first event in the first *****
*       *****  acceptable cycle, so increment the counter:  *****
*
            j=j+1               !keep track of number of acceptable cycles
            E_pointer=E_pointer+1 !point to the next respiratory cycle
            goto 605            !get next cycle
         end if
*
      end if
*
 620  if(j.eq.0)then            !cannot do any significance testing for this cell --
*         ETA2=0.0              ! no cycles are acceptable!
*         PROB=1.0
         return                 !????????????????????????????????
      end if
*
      if(j.ge.50)then           !50 acceptable cycles are available to calculate
         isubETA=50             ! the ETA**2 value and perform the ANOVA test
      else                      
         isubETA=j              !have to use less than 50 cycles - warn the user
      end if
*
      if(j.ge.101)then          !at least 101 acceptable cycles available to
         isubBIN=50             ! perform the BINARY test - good
         isubmaxBIN=isubBIN + 1 
      else                      !have to use less than 101 cycles - warn the user
         isubBIN=INT(j/2)
         isubmaxBIN=isubBIN + 1 
      end if     
      df1=float(itrt-1)                         
      df2=float(itrt*isubETA-itrt) !?????????????????? isubETA OK??
*
*       ***** Calculate the ETA**2 value and the ANOVA statistic (F test): *****
*
*
      ICNT1=0
      ICNT2=0
      ETA2=0.0                  !default value
      PROB=1.0                  !default value
      z=0.0
      RMSB=0
      RMSW=0
      ITEMPSUM1=0
      ITEMPSUM2=0
      IY = 0                    !IY(itrt,isubmax) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=0                   ! = # of acceptable respiratory cycles used so far
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Find the first selected respiratory cycle in the recording: *****
*
      cycle_start_time=begin_E(E_pointer) !start at the first control cycle
      cycle_end_time=end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
*           PRINT '(''NOT ENOUGH CYCLES ETA2 1'')' 
*           ETA2=0.0            !set default values
*            PROB=1.0
         goto 660
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 600  CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 2'')'
*          ETA2=0.0             !set default values
*           PROB=1.0
            goto 660
         ENDIF
         GOTO 600               !loop until find the FIRST valid cell event
! that occurs after the start of this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 610  FIND_THE_AISLE_6: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 611     if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*           PRINT '(''NOT ENOUGH CYCLES ETA2 3'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 650
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES ETA2 4'')' 
*              ETA2=0.0
*               PROB=1.0
               goto 650         !go to next resp sig test
            endif
*
            GOTO 610            !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 611 !itrt = total # aisles per cycle
*
         ENDIF
      end if FIND_THE_AISLE_6
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubETA)goto 650
      E_pointer=E_pointer+1     !point to next E pulse
*
      cycle_start_time=begin_E(E_pointer) !read times of new beginning and 
      cycle_end_time=end_E(E_pointer) ! ending E pulses
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''NOT ENOUGH CYCLES ETA2 5'')' 
*           ETA2=0.0            !set default values
*            PROB=1.0
         goto 650
      end if
*
 640  IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 640
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubETA) GOTO 610 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubETA (optimally 50) cycles,               *****
*       *****  divided each cycle into 20 aisles, and counted the number of         *****
*       *****  cell events that occur within each aisle;                            *****
*       ***** now we're ready for some statistical analysis of respiratory pattern: *****
*
 650  AVGY=float(IYSUM)/float(itrt*isubETA) !total # cell events in the 1st isubETA cycles /
! total # aisles in the 1st isubETA cycles =
!  # cell events expected in each aisle if
! the cell's activity is NOT respiratory modulated
      DO aisle=1,itrt
         RMSB=RMSB+(float(JSUM(aisle))
     +        /float(isubETA)-AVGY)**2
      end do
      RMSB=isubETA*RMSB/df1     !RMSB = between group sum of squares (?)
*
      DO aisle=1,itrt
         DO cycle=1,isubETA
            RMSW=RMSW+(float(IY(aisle,cycle))
     +           -float(JSUM(aisle))/float(isubETA))**2
         end do
      end do
*
      RMSW=RMSW/df2             !RMSW = within group sum of squares (?)
*
      if(RMSW.NE.0) then
         F=RMSB/RMSW
      else
         F=1.0
      end if
*
      if(RMSW.ne.0)then
         ETA2=df1*F/(df1*F+df2) !ETA squared statistic
      else
         ETA2=0.0
      end if
*
      PROB = BETAI(0.5*df2,0.5*df1,df2/(df2+df1*F))
 660  write (ETA2_6(IDs(CELL)),('(F5.2)')) ETA2
c       if(PROB.lt.SIGLVL)zmodsig_6(IDs(CELL)) = 'R'
      if(PROB.lt.SIGLVL)then
         zmodsig_6(IDs(CELL)) = 'R' !cell's activity is respiratory modulated
      else if((PROB.ge.SIGLVL).and.(isubETA.lt.50))then
         zmodsig_6(IDs(CELL)) = 'NED' !ned = not enough data
      else if((PROB.ge.SIGLVL).and.(isubETA.ge.50))then
         zmodsig_6(IDs(CELL)) = 'N' !n = not respiratory modulated
      else
      end if
*
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>
*         <<<<< PROCEED WITH BINARY TEST: >>>>>
*       <<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>
*
*
      z=0.0
      ICNT1=0
      ICNT2=0
*          RMSB=0
*          RMSW=0
      IY = 0                    !IY(itrt,isub) holds #s of events according to aisle & cycle
      JSUM = 0                  !JSUM(itrt) holds the sum of events in each aisle
      IYSUM=0                   != total of cell events 
      cycle=0                   ! = respiratory cycle designation
      aisle=1                   ! = aisle designation
      CELL_pointer= 1           !CELL_pointer = pointer for cell data
      E_pointer= 1              !E_pointer = pointer for E_pulse data
*
*
*       ***** Find the first selected respiratory cycle in the recording: *****
*
      cycle_start_time=begin_E(E_pointer)       
      cycle_end_time=end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or.
     +     (cycle_end_time.eq.0.0))then
*           PRINT '(''NOT ENOUGH CYCLES 1'')' 
         goto 6099
      end if
*
*
*      ***** find the first cell event after the first cycle start: *****
*
 6000 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
*           PRINT '(''NOT ENOUGH CYCLES 2'')'
            goto 6099
         ENDIF
         GOTO 6000              !loop until find the FIRST valid cell event
! in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 6010 FIND_THE_AISLE_6_2: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
 6011    if(CLLK.lt.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/float(itrt)))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts cell events in array subscripted
! by aisle and cycle 
            IYSUM=IYSUM + 1     !sum of all cell events which occur in 1st 50 cycles
            JSUM(aisle)=JSUM(aisle) + 1 !sum of events in each 1/itrt th
            CELL_pointer=CELL_pointer+1 !point to next cell event
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN !ran out of cell events
               goto 6099
            endif
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event time
*
            IF (CLLK.EQ.0) THEN !no more data for this cell
*           PRINT '(''NOT ENOUGH CYCLES 4'')' 
               goto 6099        !go to next resp sig test
            endif
*
            GOTO 6010           !loop to analyze the new event time
*
         ELSE                   !event does not occur within this aisle, so go to the next aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 6011 !itrt = total # aisles per cycle
*
         end if
      END IF FIND_THE_AISLE_6_2
*
*
*       ***** Now finished with this respiratory cycle, so move to next one: *****
*
      if(cycle.eq.isubBIN)goto 6041 !proceed to ETA2 and ANOVA analysis
      E_pointer=E_pointer+1     !point to next E pulse
*
      cycle_start_time=begin_E(E_pointer) !read times of new beginning and 
      cycle_end_time=end_E(E_pointer) ! ending E pulses
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''NOT ENOUGH CYCLES 5'')' 
         goto 6099
      end if
*
 6040 IF (CLLK.LT.cycle_start_time) THEN !cell event occurs prior to the start
         CELL_pointer=CELL_pointer+1 ! of this cycle, so read next time; loop
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) ! until obtain a valid time
         goto 6040
      end if
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubBIN) GOTO 6010 !see where this cell event occurs within the cycle
*
*
*       ***** now have looked at first isubBIN (optimally 50) cycles,                   *****
*       *****  divided each one into 20 aisles,                                         *****
*       *****  and counted the number of cell events that occur within each aisle;      *****
*
*       ***** make 2 groups of "contiguous" aisles                      *****
*       ***** (ie., 1-10 & 11-20; 2-11 & 12-20 + 1; etc.);              *****
*       ***** Examine all possible groupings until find the case in     *****
*       *****  which the difference in total # cell events between the  *****
*       *****  two groups is maximal:                                   *****
*
 6041 IPT1=(itrt/2) - 1         ! = # of first aisle in the 1st group
      IPT2=itrt - 1             ! = # of first aisle in the 2nd group
      IOLDDIFF=0
      IBIGPT=1
      DO I=1,itrt/2             !update pointers to first aisle in each group
         IPT1=IPT1 - itrt/2 + 1
         IPT2=IPT2 - itrt/2 + 1
         IF (IPT2.LE.0) IPT2=IPT2 + itrt
         ITEMPSUM1=0            !ITEMPSUM1 = total # cell events in 1st group
         ITEMPSUM2=0            !ITEMPSUM2 = total # cell events in 2nd group
*
         DO J=1,itrt/2
            IPT1=IPT1+1         !point to next aisle in each group
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) IPT2=IPT2-itrt
            ITEMPSUM1=ITEMPSUM1+JSUM(IPT1)
            ITEMPSUM2=ITEMPSUM2+JSUM(IPT2)
         end do
*
         IF (ITEMPSUM1.GT.ITEMPSUM2) THEN !compare the #s of events which occurred in the
            INEWDIFF=ITEMPSUM1-ITEMPSUM2 ! 1st group vs. those in the 2nd group
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT1 - itrt/2 + 1 !store # of 1st aisle in group 1
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ELSE
            INEWDIFF=ITEMPSUM2-ITEMPSUM1 !ITEMPSUM2 > ITEMPSUM1
            IF (INEWDIFF.GT.IOLDDIFF) THEN
               IBIGPT=IPT2 - itrt/2 + 1 !store # of 1st aisle in group 2
               IF (IBIGPT.LE.0) IBIGPT=IBIGPT+itrt
               IOLDDIFF=INEWDIFF
            ENDIF
         ENDIF
*
      end do
*
      totalBINARY=isubBIN       !keep track of total # cycles used in BINARY testing
*
*       ***** Now have formed a hypothesis and gathered data based on the first isubBIN cycles *****
*       *****  about whether or not the cell's activity is greater in one half of the          *****
*       *****  cycle than in the other half.  Now need to look at remaining                       *****
*       *****  cycles to see if the relationship/pattern holds true:                           *****
*
*
      ICNT1=0                   != # times ITEMPSUM1 > ITEMPSUM2
      ICNT2=0                   != # times ITEMPSUM1 <> ITEMPSUM2
 6045 IY = 0
      cycle=0
      aisle=1
*
*      ***** move to the next set of E pulses (hence, the next cycle): *****
*      
      E_pointer=E_pointer+1
*
      if (E_pointer.gt.MAX_CYCLES)goto 6099
      cycle_start_time=begin_E(E_pointer)
      cycle_end_time=end_E(E_pointer)
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*       PRINT '(''OUT OF CYCLES 6'')' 
         goto 6099
      end if
*
*      ***** find the first cell event after the first cycle start: *****
*
 6060 CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
      IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         IF ((CELL_pointer.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN !ran out of cell events
*         PRINT '(''OUT OF CYCLES 7'')'
            GOTO 6099
         ENDIF
         GOTO 6060              !loop until find 1st cell event in this cycle
      ENDIF
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))cycle=cycle+1 ! the cycle, count it as an
*                                                       ! acceptable cycle
*
*       ***** determine in which aisle of the cycle the event occurs and *****
*       *****  count the cell events in each aisle of the cycle: *****
*
 6070 FIND_THE_AISLE_6_3: if(CLLK.lt.cycle_end_time)then !if the event occurs within the cycle:
         if (CLLK.LT.(cycle_start_time+
     +        aisle*(cycle_end_time-cycle_start_time)/itrt))then
*
            IY(aisle,cycle)=IY(aisle,cycle)+1 !counts events in array subscripted
*                                                       by aisle # and cycle #
            CELL_pointer=CELL_pointer+1
*
            IF (CELL_pointer.GE.MAX_NUM_EVENTS) THEN
*         PRINT '(''OUT OF CYCLES 8'')'
               GOTO 6098
            ENDIF
*
            CLLK=SPIKETIMES(CELL_pointer,IDs(CELL)) !read next cell event
*
            IF (CLLK.EQ.0) THEN
*         PRINT '(''OUT OF CYCLES 9'')'
               GOTO 6098
            ENDIF
*
            GOTO 6070           !loop back to analyze the new event time
*
         ELSE                   !cell event does not occur in current aisle
            aisle=aisle+1
            IF (aisle.LE.itrt) GOTO 6070 !look at the next aisle (or go to the next cycle)
*
         end if
      END IF FIND_THE_AISLE_6_3
*
      if(cycle.eq.isubmaxBIN)goto 6098
      E_pointer=E_pointer+1     !proceed to the next resp cycle
*
      cycle_start_time=begin_E(E_pointer)
      cycle_end_time=end_E(E_pointer)
*
      if((cycle_start_time.eq.0.0).or. !be sure that both E pulses for this
     +     (cycle_end_time.eq.0.0))then ! cycle occur within the control period
*           PRINT '(''OUT OF CYCLES 10'')' 
         goto 6098
      end if
*
 6076 IF (CLLK.LT.cycle_start_time) THEN
         CELL_pointer=CELL_pointer+1
         CLLK=SPIKETIMES(CELL_pointer,IDs(CELL))
         IF (CLLK.EQ.0) GOTO 6098
         GOTO 6076              !loop until read a valid cell event time
      ENDIF
*
*       ***** now have a cell event time which may occur within this new cycle: *****
*
      if((CLLK.ge.cycle_start_time).and. !if the event occurs within
     +     (CLLK.lt.cycle_end_time))then ! the cycle, count it as an
         cycle=cycle+1          ! acceptable cycle and
         aisle = 1              ! reset / update array pointers
      end if
*
      IF (cycle.LE.isubmaxBIN) GOTO 6070 !for the next cycle
*
*       ***** look at cell's activity within this group of  cycles to see if     *****
*       *****  the relationship/pattern detected with the 1st isubBIN cycles holds true: *****
*
 6098 DO icycle=1,cycle-1       !for each cycle in the group:
         totalBINARY = totalBINARY+1 !keep a total of # of cycles used in BINARY testing
         ITEMPSUM1=0
         ITEMPSUM2=0
         IPT1=IBIGPT-1
         IPT2=IBIGPT+itrt/2-1
         IF (IPT2.GT.itrt) IPT2=IPT2-itrt
         DO J=1,itrt/2
            IPT1=IPT1+1
            IF (IPT1.GT.itrt) THEN
               ITEMPSUM1=ITEMPSUM1+IY(IPT1-itrt,icycle+1)
            ELSE
               ITEMPSUM1=ITEMPSUM1+IY(IPT1,icycle)
            ENDIF
            IPT2=IPT2+1
            IF (IPT2.GT.itrt) THEN
               ITEMPSUM2=ITEMPSUM2+IY(IPT2-itrt,icycle+1)
            ELSE
               ITEMPSUM2=ITEMPSUM2+IY(IPT2,icycle)
            ENDIF
         end do
*

         IF (ITEMPSUM1.GT.ITEMPSUM2) ICNT1=ICNT1+1 !ICNT1=# cycles which support the hypothesis
         If (ITEMPSUM1.NE.ITEMPSUM2) ICNT2=ICNT2+1 !ICNT2=# cycles in which the cell's activity in
*                                                       ! one half of the cycle does not equal the
*                                                       ! activity in the other half of the cycle
      end do
*
      goto 6045                 
*
 6099 z=(float(ICNT1)-float(ICNT2)/2.)/(.5*SQRT(FLOAT(ICNT2)))
*         print '(''cell '',I4,'' ICNT1 = '',I5,'' ICNT2 = '',I5,
*     +     '' z = '',f15.5)',CELL,ICNT1,ICNT2,z
*         read '(A)'
      if(z.gt.ZSIGLVL)then
         zmodsig2_6(IDs(CELL)) = 'R' ! bsig = Z value
      else if((z.le.ZSIGLVL).and.(isubBIN.lt.50))then
         zmodsig2_6(IDs(CELL)) = 'NED' !ned = not enough data for a good test
      else if((z.le.ZSIGLVL).and.(isubBIN.ge.50))then
         zmodsig2_6(IDs(CELL)) = 'N'
      else
      end if
*
*
*
*
*       **********************************************************************
*       **********************************************************************
*
*
*
      return
      END 
      end module mod_respsig6_v4
