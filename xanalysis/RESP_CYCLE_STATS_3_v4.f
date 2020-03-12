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

      module mod_RESP_CYCLE_STATS_3_v4
      contains
*       filename = RESP_CYCLE_STATS_3_v4.f
*
*       date of latest revision = 05-oct-2005   lss 
*
*       apr-2001        lss
*        BOUNDARY CODES ARE NO LONGER REQUIRED FOR CTH GENERATION!
*        modified to allow user to select respiratory cycles for CTH calculation:
*               A -- use all cycles in the recording
*               F -- use first 10000 (=MAX_CYCLES) cycles in the recording
*               B -- use all cycles within a boundary coded section of the recording
*               C -- use the first x cycles in the recording
*               T -- user has designated a period of time within the recording
*
*       mar-2000        lss
*        modified to allow ID codes to range from 1 to 999;     
*        maximum number of units allowed still = 120;           
*        *.defs files inserted                                  
*        (MAX_NUM_CODES=999; MAX_NUM_CHAN=120)                  
*
*
*       ***** link with main program x2000_v5 *****
*
*       This subroutine allows the user to select the respiratory cycles which will be used
*       to calculate cycle-triggered histograms, both regular and normalized.
*       Three criteria may be used to select respiratory cycles:
*
*       1.  The user must select which cycles are to be considered for inclusion 
*               (see modifications for apr-2001, above).  THIS CRITERION MUST BE MET.
*
*       2.  The cycles must be of an acceptable duration.  The user may, if desired,
*               eliminate cycles that are too short or too long.  This is of use when there
*               is jitter in the E pulse record or when the control period has been edited
*               to exclude undesirable cycles, thereby artificially creating cycles of long
*               duration.
*
*       3.  The user may, if desired, select cycles based upon certain characteristics
*               of the inspiratory and/or expiratory period(s).  For example, if the data
*               indicates that there are periods of apneusis within the control period, 
*               the user may elect to eliminate those respiratory cycles whose duration is
*               greater than a user-entered value.
*
*
*       This subroutine also derives values that will be used later to look at the
*        significance of a cell's respiratory firing pattern.  Respiratory significance
*        is evaluated under several sets of circumstances:
*         1. using the entire data set
*               a. the user may, if desired, eliminate cycles whose duration is
*                  greater than a user-defined value
*         2. using events which occur only within the selected control period
*               a. use all events (mean_control_cycle = mean control cycle duration)
*               b. the user may, if desired, eliminate cycles whose duration is
*                  greater than the mean + 2 SD of control cycles 
*               c. OR the user may eliminate cycles whose duration is greater than
*                  a user-defined value
*         3. using only those control cycles actually used to calculate CTHs
*
*       These calculations are made before the user has a chance to eliminate respiratory
*       cycles from the control period.
*
*       **********************************************************************************
*       **********************************************************************************
*
*

      function z(s)
      character(len=*) :: s
      character(len=len_trim(s) + 1) :: z
      n = len_trim(s) + 1
      z = s
      z(n:n) = char (0)
      end function

      subroutine RESP_CYCLE_STATS(SPIKETIMES,I_pulse,E_pulse,
     +     BOUNDARY,icycles,
     +     BNDRY_START_TIME,BNDRY_END_TIME,phrenic,
     +     ITAL,IHIST,
     +     fildes,mouse,mean_resp_cycle,
     +     std_dev,IMAXCYC,begin_E,end_E,
     +     TOTAL_SELECTED_CYCLES,cycle_durations,
     +     rec_BW,rec_offset,mean_control_cycle,
     +     control_max_value,IDs,mean_resp_cycle_all,
     +     std_dev_all,longest_control,shortest_control,
     +     mean_E,shortest,longest,DB_FILES,middle_I)
      use mod_clear_routines
      use mod_locate_region2
      use mod_miscellaneous_subroutines
      use mod_new_draw_button
      use mod_new_plot

      type minmax
         real MIN_CYCLE
         real MAX_CYCLE
         real I_MIN
         real I_MAX
         real E_MIN
         real E_MAX
      end type minmax

      type (minmax) :: MM

      INCLUDE 'gopen_type.defs'
*
*
      include 'x2000parameter.defs'
*
      double precision SPIKETIMES(:,:)
*
      integer IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN)
*
      integer E_pulse,I_pulse,phrenic,TOTAL_SELECT_E,
     +     TOTAL_NUM_CYCLES,resp_offset,IHIST(101),
     +     TOTAL_SELECTED_CYCLES,I_BIN,E_BIN,ORIG_Ti_HIST(101),
     +     ORIG_Te_HIST(101),rep,ORIG_DIST_HIST(101),
     +     TOTAL_NUM_CONTROL_CYCLES,IHIST1(101),criterion
*
      integer*4 mouse,region1,ICN0(6)
*
      DOUBLE PRECISION
     +     BNDRY_START_TIME,BNDRY_END_TIME

      DOUBLE PRECISION
     +     begin_E(:),
     +     end_E(:),
     +     middle_I(:)
      REAL
     +     cycle_durations(:)

      DOUBLE PRECISION, allocatable ::
     +     control_end_E(:),
     +     control_begin_E(:),
     +     control_middle_I(:),
     +     TEMP_I(:),
     +     TEMP_E(:)
      REAL, allocatable ::
     +     init_I_durations(:),
     +     init_E_durations(:),
     +     I_durations(:),
     +     E_durations(:),
     +     control_cycle_durations(:)

      real mean_resp_cycle, MIN_CYCLE,MAX_CYCLE,
     +     shortest,longest,std_dev,sum_of_I,sum_of_e,I_MIN,
     +     I_MAX,E_MIN,E_MAX,longest_I,shortest_I,longest_E,
     +     shortest_E,longest_I_all,longest_E_all,
     +     longest_control,shortest_control,rec_BW,rec_offset,
     +     mean_I,std_dev_I,mean_E,mean_I_all,mean_E_all,
     +     mean_control_cycle,mean_resp_cycle_all,
     +     control_max_value,x0,y0,height,width,
     +     pixels_per_bin,std_dev_all
*
      character*(*) DB_FILES
      character*120 text
      character*20 c_short_selected,c_short_all,c_long_selected,
     +     c_long_all,c_mean_selected,
     +     c_mean_all,c_std,c_std_all,c_bw
      character*10 c_selected,c_total
*
      character*12 again
      character*3 DIS,OK
*
      character*1 info,BOUNDARY

      double precision :: sum_of_squares_E = 0.0
      double precision :: sum_of_squares_I = 0.0
      logical debug, use_minmax
      real :: zero = 0
*       ***** suppress unused variable warnings *****

      if(.false.)print *,BOUNDARY
      if(.false.)print *,icycles
      if(.false.)print *,mean_control_cycle
      if(.false.)print *,mouse
      if(.false.)print *,phrenic
      if(.false.)print *,TOTAL_NUM_CYCLES
      if(.false.)print *,fildes2 !suppress unused variable warning

      max_cycles = size(end_E)
      
      allocate (control_end_E(max_cycles))
      allocate (control_begin_E(max_cycles))
      allocate (control_middle_I(max_cycles))
      allocate (TEMP_I(max_cycles+1))
      allocate (TEMP_E(max_cycles+1))
      allocate (init_I_durations(max_cycles))
      allocate (init_E_durations(max_cycles))
      allocate (I_durations(max_cycles))
      allocate (E_durations(max_cycles))
      allocate (control_cycle_durations(max_cycles))

      ICN0 = 0
      MM%I_MIN = 0
      MM%I_MAX = 0
      MM%E_MIN = 0
      MM%E_MAX = 0
      MM%MIN_CYCLE = 0
      MM%MAX_CYCLE = 0

C     The calculation of rec_BW in the first call to CALCULATE depends
C     on these before they've been calculated, but that value of rec_BW
C     is never used. We initialize them here to NaN so we will know if the
C     bad value is ever used.
      mean_I = 0/zero
      std_dev_I = 0/zero

C     The last test of region1 in this subroutine can happen without
C     region1 having been set elsewhere, so we initialize it here.
      region1 = 0

      call getenv ("XANALYSIS_DEBUG",text)
      debug = text.eq.'on'
      call getenv ("USE_MINMAX",text)
      use_minmax = text.eq.'on'
*
*
*       *****   set screen and window parameters  *****
*
      call vdc_extent(fildes,0.0,0.0,0.0,1.0,1.0,0.0)
      call mapping_mode(fildes,1)
      call view_port(fildes,0.,0.,1.,1.)
      call view_window(fildes,1.,1.,1500.,600.)
*
*
*       ***** initialize variables and define array usage: *****
*
      DIS = 'DIS'
      OK = 'OK'
      again = 'SELECT RANGE'
*
      x0 = 100.                 !define physical parameters for plot of
      y0 = 200.                 !  respiratory duration distribution, 
      x0_I = 100.               !  I phase duration distribution, and
      y0_I = 300.
      x0_E = 100.               !  E phase duration distribution histograms
      y0_E = 21.
      height = 250.                   
      width = 500.
*
      TOTAL_SELECT_E = 0
      I_BIN = 0
      E_BIN = 0
      do i = 1, 101
         IHIST(i) = 0
         ORIG_DIST_HIST(i) = 0  !holds cycle duration distribution histogram computed
                                !  using all cycles in the control period.  Histograms
                                !  computed for subsequently "selected" cycles will
                                !  be compared to this original histogram
         ORIG_Ti_HIST(i) = 0
         ORIG_Te_HIST(i) = 0
      end do

      control_begin_E = 0.0     !holds beginning E pulse for each selected cycle
      control_end_E = 0.0       !holds ending E pulse for each selected cycle
      control_middle_I = 0.0    !holds the I pulse included in each selected cycle
      begin_E = 0.0             !holding arrays
      end_E = 0.0
      middle_I = 0.0

      sum_of_I = 0.0
      sum_of_E = 0.0
      resp_offset = 0
      MIN_CYCLE = 0.0
      MAX_CYCLE = 0.0
      I_MIN = 0.0
      I_MAX = 0.0
      E_MIN = 0.0
      E_MAX = 0.0
      bin_duration = 0.0        !bin width used to display selected cycles
      TOTAL_NUM_CYCLES = 0      !total # of cycles in entire record
      TOTAL_SELECTED_CYCLES = 0 !total # of cycles that meet all selection criteria
      TOTAL_NUM_CONTROL_CYCLES = 0 !total # of cycles in control period

      TEMP_E = 0.0              !holds times of E pulses for selected cycles
      TEMP_I = 0.0              !holds times of I pulses for selected cycles

      I_durations = 0.0         !holds durations of Ti for all selected cycles
      E_durations = 0.0         !holds durations of Te for all selected cycles
      init_I_durations = 0.0
      init_E_durations = 0.0

      rep = 0
      pixels_per_bin = width/100.
*
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                                               *
*       *       CRITERION #1:  BASIC SELECTION ROUTINE (NOT OPTIONAL)                   *
*       *                                                                               *
*       *       A -- use all cycles in the recording                                    *
*       *       F -- use first 10000 (=MAX_CYCLES) cycles in the recording              *
*       *       B -- use all cycles within a boundary-coded section of the recording    *
*       *       C -- use the first x cycles in the recording                            *       
*       *       T -- user has designated a period of time within the recording          *
*       *                                                                               *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
      criterion = 1             !flag to subroutine CALCULATE
      j=0
      do i=1,ITAL(IDs(E_pulse))
         if((SPIKETIMES(i,IDs(E_pulse)).ge.
     +        BNDRY_START_TIME).and. !be sure that selected E pulses
     +        (SPIKETIMES(i,IDs(E_pulse)).le. ! occur w/in the control period
     +        BNDRY_END_TIME))then                    
            j=j+1               !update the pointer
            TEMP_E(j)=SPIKETIMES(i,IDs(E_pulse))
         end if
         if(SPIKETIMES(i,IDs(E_pulse)).gt.BNDRY_END_TIME)
     +        goto 50
      end do

 50   TOTAL_SELECT_E = j        !TOTAL_SELECT_E = # of E pulses in TEMP_E
*
*
*       *****   find the I pulses associated with the above E pulses:   *****
*
      FILL_TEMP_I: do i = 1,(TOTAL_SELECT_E-1) !for every E pulse in the control
      do j = 1,ITAL(IDs(I_pulse)) !  period (except the last one),
         if((SPIKETIMES(j,IDs(I_pulse)).gt.TEMP_E(i)).and. !  find the I pulse that
     +        (SPIKETIMES(j,IDs(I_pulse)).lt.TEMP_E(i+1)))then !  occurs immediately after it but
            TEMP_I(i)=SPIKETIMES(j,IDs(I_pulse)) !  before the next E pulse
            cycle FILL_TEMP_I
         end if
      end do
      end do FILL_TEMP_I
*
*       *****   The data now contained in TEMP_E and TEMP_I reflect the following:      *****
*       *****                                                                           *****
*       *****           E----I----E----I----E---I----E....I----E                        *****
*       *****                                                                           *****
*       *****           # of E = n and # E-E cycles = (n-1)                             *****
*       *****                                                                           *****
*       *****           so ...............                                              *****
*
      TOTAL_NUM_CONTROL_CYCLES=TOTAL_SELECT_E-1 != total number of cycles included
                                !  in control period  
      TOTAL_NUM_CYCLES=ITAL(IDs(E_pulse)) - 1 != total number of cycles in entire 
                                !  recording
      TOTAL_SELECTED_CYCLES=TOTAL_SELECT_E-1 != number of cycles selected thus far
*
*       ***** calculate the mean cycle duration for ALL cycles in the CONTROL period *****
*       *****   = mean_control_cycle   (control_max_value = mean + 2 s.d.)           *****
*
      call read_minmax (Z(DB_FILES), MM, MMOK)
 100  call CALCULATE (criterion,shortest,longest,TEMP_E,
     +     mean_resp_cycle,
     +     std_dev,IMAXCYC,MAX_CYCLE,MIN_CYCLE,
     +     TEMP_I,TOTAL_NUM_CYCLES,TOTAL_NUM_CONTROL_CYCLES,
     +     cycle_durations,TOTAL_SELECTED_CYCLES,begin_E,
     +     end_E,middle_I,ORIG_DIST_HIST,
     +     control_cycle_durations,control_begin_E,
     +     control_end_E,control_middle_I,fildes,
     +     longest_control,shortest_control,rec_BW,
     +     rec_offset,mean_I,std_dev_I,x0,y0,
     +     width,bin_duration,control_max_value,
     +     mean_resp_cycle_all,std_dev_all)

*
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                                               *
*       *       CRITERION #2:   cycle durations must be within a minimum and a maximum  *
*       *                       value as specified by the user  (OPTIONAL)              *
*       *                                                                               *
*       *               MIN_CYCLE <= x >= MAX_CYCLE                                     *
*       *                                                                               *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
      criterion = 2
*
*       ***** allow user to enter new values for MIN and MAX_CYCLE, if desired: *****
*
*
      if (MMOK.eq.1) then
         MIN_CYCLE = MM%MIN_CYCLE
         MAX_CYCLE = MM%MAX_CYCLE
         MMOK = 2
         if(MIN_CYCLE.ne.0.or.MAX_CYCLE.ne.0)goto 100
      else
         MIN_CYCLE=0.0          !reset variable values
         MAX_CYCLE=0.0
      end if
*
*
*       call clear_all(fildes)
      text = 'DONE'
      call draw_button (fildes,1400.,110.,1450.,160.,text//char(0),    
     +     1405.,85.," ",0.,0.)                               
      text='You may choose the range of cycle durations that will '
     +     //'be used to calculate the respiratory CTHs.'
      call text2d(fildes,610.,575.,text//char(0),
     +     ANNOTATION_TEXT,0)
      text='Keep choosing until you get it the way you want it.'
     +     //'  Always choose the lower end of the range first.'
      call text2d(fildes,610.,555.,text//char(0),
     +     ANNOTATION_TEXT,0)
      text='Click DONE when you are satisfied with your choice.'
      call text2d(fildes,610.,530.,text//char(0),
     +     ANNOTATION_TEXT,0)
      call make_picture_current(fildes)
      if(use_minmax.and.MMOK.eq.2) then
         region1 = 91
         goto 226
      endif

 225  call request_locator
     +     (fildes,1,2e9,valid,x_dc,y_dc,z_dc)
      call vdc_to_wc(fildes,x_dc,y_dc,z_dc,
     +     xloc,yloc,zloc)      !translate into fp coords.
      call locate_region2(xloc,yloc,region1)
*
 226  if(region1.eq.91)then     !DONE - go on to criterion #3
         call text2d(fildes,900.,750.,text//char(0),
     +        ANNOTATION_TEXT,0)
         MIN_CYCLE = shortest_control
         MAX_CYCLE = longest_control
         call clear_all(fildes)
         goto 290
      else                      !select range of bins
         call clear_bottom_CTH(fildes)
         text = 'DONE'
         call draw_button (fildes,1400.,110.,1450.,160.,text//char(0),  
     +        1405.,85.," ",0.,0.)                             
         min_selected=INT((xloc-x0)/pixels_per_bin)+1 !first bin of selected cycles
         if(min_selected.lt.1)min_selected = 1 !allow a pick to the left of the histogram
         if(min_selected.gt.100)then !oops - picked to the right of the histogram
            call text_color(fildes,1.,0.,0.) !red text
            text='Select the MINIMUM cycle duration first.'
            call text2d(fildes,x0+100.,55.,text//char(0),
     +           ANNOTATION_TEXT,0)
            call make_picture_current (fildes)
            call text_color(fildes,0.,0.,0.) !black text
            goto 225
         end if
         if(debug)min_selected = 1
         MIN_CYCLE=shortest_control+(min_selected-1)*bin_duration !minimum allowed cycle time
         call clear_bottom_CTH(fildes)
         call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
         call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc)
         max_selected=INT((xloc-x0)/pixels_per_bin)+1 !last bin of selected cycles
         if(max_selected.lt.min_selected)then !can't have a MAX less than the MIN!
            call text_color(fildes,1.,0.,0.) !red text
            text='No, no, no.  Bad user.  Start over and remember '
     +           //'to select the MINIMUM duration first.'
            call text2d(fildes,x0+100.,55.,text//char(0),
     +           ANNOTATION_TEXT,0)
            call make_picture_current (fildes)
            call text_color(fildes,0.,0.,0.) !black text
            goto 225
         end if
         if(max_selected.gt.100)max_selected=100
         if(debug)max_selected = 6
         MAX_CYCLE=shortest_control+(max_selected*bin_duration) !maximum allowed cycle time
*         pause
         call clear_bottom_CTH(fildes) !clear bottom of window
         MM%MIN_CYCLE = MIN_CYCLE
         MM%MAX_CYCLE = MAX_CYCLE
         goto 100
      end if
*
*
*       ***** Selected respiratory cycles have met criteria #1 and #2: (1) each         *****
*       *****   cycle occurs within the control period and (2) the duration of each     *****
*       *****   cycle is within the limits as set by MIN_CYCLE and MAX_CYCLE.           *****
*       *****                                                                           *****
*       ***** Re-load "control_" arrays with durations of and I and E pulses assoc-     *****
*       *****   iated with selected cycles.  Use of "control_" and "init_" arrays       *****
*       *****   allows the user to enter any desired values to further select           *****
*       *****   cycles based upon times of inspiration and expiration.                  *****
*
 290  control_begin_E = begin_E
      control_end_E = end_E
      control_middle_I = middle_I
      control_cycle_durations = cycle_durations
      CRITERION2_CYCLES = TOTAL_SELECTED_CYCLES !store the number of cycles
                                ! that have met criteria #1 & #2
*
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                                               *
*       *       CRITERION #3:   select cycles based upon certain characteristics of the *
*       *                       inspiratory and/or expiratory cycle(s):  (OPTIONAL)     *
*       *                                                                               *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
      criterion = 3
*
      do i = 1,TOTAL_SELECTED_CYCLES !calculate Ti and Te for each selected
         init_I_durations(i)=end_E(i)-middle_I(i) !  cycle and store the values in
         init_E_durations(i)=middle_I(i)-begin_E(i) !  init_I and init_E_durations()
         sum_of_I = sum_of_I + init_I_durations(i)
         sum_of_E = sum_of_E + init_E_durations(i)
      end do
      I_durations = init_I_durations !fill the holding arrays
      E_durations = init_E_durations
      mean_I = sum_of_I / float(TOTAL_SELECTED_CYCLES)
      mean_E = sum_of_E / float(TOTAL_SELECTED_CYCLES)
      mean_I_all = mean_I       !store these values for later display
      mean_E_all = mean_E
*
*       *****  calculate standard deviation of mean Ti and mean Te: *****
*
      do i=1,TOTAL_SELECTED_CYCLES  
         sum_of_squares_I = sum_of_squares_I + 
     +        ((I_durations(i)-mean_I)**2)
         sum_of_squares_E = sum_of_squares_E + 
     +        ((E_durations(i)-mean_E)**2)
      end do
      std_dev_I=SQRT(sum_of_squares_I/
     +     float(TOTAL_SELECTED_CYCLES-1))
      std_dev_I_all = std_dev_I !store this value for later display
      std_dev_E=SQRT(sum_of_squares_E/
     +     float(TOTAL_SELECTED_CYCLES-1))
      std_dev_E_all = std_dev_E !store this value for later display
*
*
*
*       *****   allow user to select cycles based on duration(s) of     *****
*       *****           inspiration (Ti) and/or expiration (Te):        *****
*
      longest_I = I_durations(1) !load comparison variables
      shortest_I = I_durations(1)
      longest_E = E_durations(1)
      shortest_E = E_durations(1)
*
      do i = 1,TOTAL_SELECTED_CYCLES !find longest and 
         longest_I = AMAX1(I_durations(i),longest_I) !  shortest Ti and Te
         shortest_I = AMIN1(I_durations(i),shortest_I)
         longest_E = AMAX1(E_durations(i),longest_E)
         shortest_E = AMIN1(E_durations(i),shortest_E)
      end do
*
      I_MIN = shortest_I
      I_MAX = longest_I                       
      E_MIN = shortest_E
      E_MAX = longest_E
*
      shortest_I_all= shortest_I !store these parameters for later display
      longest_I_all = longest_I ! and calculations
      shortest_E_all = shortest_E
      longest_E_all = longest_E
*
      binwidth_I=(longest_I+.1-shortest_I)/100.0 !these binwidths are constant and will be
      binwidth_E=(longest_E+.1-shortest_E)/100.0 !  used to display Ti and Te distribution
                                !  histograms - allows overlay of cycles
                                !  selected according to criterion #3 over 
                                !  all cycles meeting #1 and #2
      call CALCULATE (criterion,shortest,longest,TEMP_E,
     +     mean_resp_cycle,
     +     std_dev,IMAXCYC,MAX_CYCLE,MIN_CYCLE,
     +     TEMP_I,TOTAL_NUM_CYCLES,TOTAL_NUM_CONTROL_CYCLES,
     +     cycle_durations,TOTAL_SELECTED_CYCLES,begin_E,
     +     end_E,middle_I,ORIG_DIST_HIST,
     +     control_cycle_durations,control_begin_E,
     +     control_end_E,control_middle_I,fildes,
     +     longest_control,shortest_control,rec_BW,
     +     rec_offset,mean_I,std_dev_I,x0,y0,
     +     width,bin_duration,control_max_value,
     +     mean_resp_cycle_all,std_dev_all)
*
*       *****   construct the Ti and Te distribution histograms using cycles    *****
*       *****           that have been selected thus far:                       *****
*
 305  rep = rep + 1             !keep track of how many times this section of code has
                                !  been entered
      IHIST = 0                 !initialize the histogram array
      do i = 1,TOTAL_SELECTED_CYCLES
         I_BIN=((I_durations(i)-shortest_I_all)/binwidth_I)+1
         IHIST(I_BIN) = IHIST(I_BIN)+1
      end do
      IHIST(101)=0
      if(rep.eq.1)ORIG_Ti_HIST=IHIST !save the original Ti distribution 
*                                               !  histogram 
      NHW=(binwidth_I*50.)
      info='n'
      call new_plot(fildes,100,ORIG_Ti_HIST,IHIST,IHIST1,DIS,x0_I,y0_I,
     +     height,width,info,'','','','',NHW,0,'',0,'',ICN0,0,'',
     +     '','',[integer::],'','','',0,0.,[integer::],[integer::])
      call character_height(fildes,.060)
      call character_width (fildes,.015)
      write (c_selected,'(I5)') TOTAL_SELECTED_CYCLES
      call remove_leading_blanks(c_selected,LEN(c_selected))
      call strlength(c_selected,LEN(c_selected),l_s)
      write (c_total,'(I5)') TOTAL_NUM_CONTROL_CYCLES
      call remove_leading_blanks(c_total,LEN(c_total))
      call strlength(c_total,LEN(c_total),l_t)
      text=c_selected(1:l_s)//' of '//c_total(1:l_t)//
     +     ' control cycles selected'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0_I+200.,y0_I+height+75.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      call character_height(fildes,.060)
      call character_width (fildes,.020)
      call text_orientation2d(fildes,-1.,0.,0.,1.) !-sind(90.),cosd(90.), cosd(90.),sind(90.)
      text='I PHASE DURATIONS'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0_I-30.,y0_I,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call character_height(fildes,.020)
      call character_width(fildes,.015)
      write (c_bw,'(f7.1)') binwidth_I
      call remove_leading_blanks(c_bw,LEN(c_bw))
      call strlength(c_bw,LEN(c_bw),l)
      text='binwidth = '//c_bw(1:l)//' msec.'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0_I-10.,y0_I+10.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call text_orientation2d(fildes,0.,1.,1.,0.)
      call character_height(fildes,.040)
      call character_width(fildes,.010)
      write (c_bw,'(f7.1)') binwidth_I
      call remove_leading_blanks(c_bw,LEN(c_bw))
      call strlength(c_bw,LEN(c_bw),l_bw)
      write (c_short_selected,'(f7.1)') shortest_I
      call remove_leading_blanks(c_short_selected,
     +     LEN(c_short_selected))
      call strlength(c_short_selected,LEN(c_short_selected),l_ss)
      write (c_short_all,'(f7.1)') shortest_I_all
      call remove_leading_blanks(c_short_all,LEN(c_short_all))
      call strlength(c_short_all,LEN(c_short_all),l_sa)
      write (c_long_selected,'(f7.1)') longest_I
      call remove_leading_blanks(c_long_selected,
     +     LEN(c_long_selected))
      call strlength(c_long_selected,LEN(c_long_selected),l_ls)
      write (c_long_all,'(f7.1)') longest_I_all
      call remove_leading_blanks(c_long_all,LEN(c_long_all))
      call strlength(c_long_all,LEN(c_long_all),l_la)
      write (c_mean_selected,'(f7.1)') mean_I
      call remove_leading_blanks(c_mean_selected,
     +     LEN(c_mean_selected))
      call strlength(c_mean_selected,LEN(c_mean_selected),l_ms)
      write (c_mean_all,'(f7.1)') mean_I_all
      call remove_leading_blanks(c_mean_all,LEN(c_mean_all))
      call strlength(c_mean_all,LEN(c_mean_all),l_ma)
      write (c_std,'(f7.1)') std_dev_I
      call remove_leading_blanks(c_std,LEN(c_std))
      call strlength(c_std,LEN(c_std),l_std)
      write (c_std_all,'(f7.1)') std_dev_I_all
      call remove_leading_blanks(c_std_all,LEN(c_std_all))
      call strlength(c_std_all,LEN(c_std_all),l_std_a)
      call text2d(fildes,x0_I-30.,y0_I-20.,
     +     c_short_all(1:l_sa)//char(0),ANNOTATION_TEXT,0)
      call joinstrings(c_long_all(1:l_la),' ms',text,l)
      call text2d(fildes,x0_I+width-50.,y0_I-20.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      text='shortest Ti:'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+150.,y0+height-130.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-130.,
     +     c_short_selected(1:l_ss)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-130.,
     +     c_short_all(1:l_sa)//char(0),ANNOTATION_TEXT,0)
      text='longest Ti:'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+150.,y0+height-150.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-150.,
     +     c_long_selected(1:l_ls)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-150.,
     +     c_long_all(1:l_la)//char(0),ANNOTATION_TEXT,0)

      text='mean Ti:'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+150.,y0+height-170.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-170.,
     +     c_mean_selected(1:l_ms)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-170.,
     +     c_mean_all(1:l_ma)//char(0),ANNOTATION_TEXT,0)

      text='standard dev:  +/-'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+150.,y0+height-190.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-190.,
     +     c_std(1:l_std)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-190.,
     +     c_std_all(1:l_std_a)//char(0),ANNOTATION_TEXT,0)

*
      IHIST = 0                 !initialize the histogram array
      do i = 1,TOTAL_SELECTED_CYCLES
         E_BIN=((E_durations(i)-shortest_E_all)/binwidth_E)+1
         IHIST(E_BIN) = IHIST(E_BIN)+1
      end do
      IHIST(101)=0
      if(rep.eq.1)ORIG_Te_HIST=IHIST !save the original Te distribution 
*                                               !  histogram
      if (MMOK.eq.2) then
         I_MIN = MM%I_MIN
         I_MAX = MM%I_MAX
         E_MIN = MM%E_MIN
         E_MAX = MM%E_MAX
         MMOK = 3
         region1 = 0
         goto 311
      end if

      NHW = (binwidth_E*50.)
      call new_plot(fildes,100,ORIG_Te_HIST,IHIST,IHIST1,DIS,x0_E,y0_E,
     +     height,width,info,'','','','',NHW,0,'',0,'',ICN0,0,'',
     +     '','',[integer::],'','','',0,0.,[integer::],[integer::])
*
      call character_height(fildes,.060)
      call character_width (fildes,.020)
      call text_orientation2d(fildes,-1.,0.,0.,1.) !-sind(90.),cosd(90.), cosd(90.),sind(90.)
      text='E PHASE DURATIONS'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0_E-30.,y0_E,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call character_height(fildes,.020)
      call character_width(fildes,.015)
      write (c_bw,'(f7.1)') binwidth_E
      call remove_leading_blanks(c_bw,LEN(c_bw))
      call strlength(c_bw,LEN(c_bw),l)
      text='binwidth = '//c_bw(1:l)//' msec.'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0_E-10.,y0_E+10.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call text_orientation2d(fildes,0.,1.,1.,0.)
*
      write (c_short_selected,'(f7.1)') shortest_E
      call remove_leading_blanks(c_short_selected,
     +     LEN(c_short_selected))
      call strlength(c_short_selected,LEN(c_short_selected),l_ss)
      write (c_short_all,'(f7.1)') shortest_E_all
      call remove_leading_blanks(c_short_all,LEN(c_short_all))
      call strlength(c_short_all,LEN(c_short_all),l_sa)
      write (c_long_selected,'(f7.1)') longest_E
      call remove_leading_blanks(c_long_selected,
     +     LEN(c_long_selected))
      call strlength(c_long_selected,LEN(c_long_selected),l_ls)
      write (c_long_all,'(f7.1)') longest_E_all
      call remove_leading_blanks(c_long_all,LEN(c_long_all))
      call strlength(c_long_all,LEN(c_long_all),l_la)
      write (c_mean_selected,'(f7.1)') mean_E
      call remove_leading_blanks(c_mean_selected,
     +     LEN(c_mean_selected))
      call strlength(c_mean_selected,LEN(c_mean_selected),l_ms)
      write (c_mean_all,'(f7.1)') mean_E_all
      call remove_leading_blanks(c_mean_all,LEN(c_mean_all))
      call strlength(c_mean_all,LEN(c_mean_all),l_ma)
      write (c_std,'(f7.1)') std_dev_E
      call remove_leading_blanks(c_std,LEN(c_std))
      call strlength(c_std,LEN(c_std),l_std)
      write (c_std_all,'(f7.1)') std_dev_E_all
      call remove_leading_blanks(c_std_all,LEN(c_std_all))
      call strlength(c_std_all,LEN(c_std_all),l_std_a)
      call character_height(fildes,.040)
      call character_width(fildes,.010)
      call text2d(fildes,x0_E-30.,y0_E-20.,
     +     c_short_all(1:l_sa)//char(0),
     +     ANNOTATION_TEXT,0)
      call joinstrings(c_long_all(1:l_la),' ms',text,l)
      call text2d(fildes,x0_E+width-50.,y0_E-20.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)

      text='shortest Te:'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+150.,y0+height-220.,
     +     text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-220.,
     +     c_short_selected(1:l_ss)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-220.,
     +     c_short_all(1:l_sa)//char(0),ANNOTATION_TEXT,0)
      text='longest Te:'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+150.,y0+height-240.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-240.,
     +     c_long_selected(1:l_ls)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-240.,
     +     c_long_all(1:l_la)//char(0),ANNOTATION_TEXT,0)
      text='mean Te:'
      call text2d(fildes,x0+width+150.,y0+height-260.,text//char(0),
     +     ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-260.,
     +     c_mean_selected(1:l_ms)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-260.,
     +     c_mean_all(1:l_ma)//char(0),ANNOTATION_TEXT,0)
      text='standard dev:  +/-'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+150.,y0+height-280.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-280.,
     +     c_std(1:l_std)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-280.,
     +     c_std_all(1:l_std_a)//char(0),ANNOTATION_TEXT,0)

      text='You may choose the range of Ti''s and/or Te''s '
     +     //'that will be used to calculate the respiratory CTHs.'
      call text2d(fildes,600.,575.,text//char(0),
     +     ANNOTATION_TEXT,0)
      text='You may use either histogram to make your choice.'
     +     //'  Keep choosing until you get it the way you want it.'
      call text2d(fildes,600.,555.,text//char(0),
     +     ANNOTATION_TEXT,0)
      text='Click DONE when you are satisfied with your choice.'
      call text2d(fildes,600.,530.,text//char(0),
     +     ANNOTATION_TEXT,0)
      call make_picture_current(fildes)
 310  text = 'DONE'
      call draw_button (fildes,1400.,110.,1450.,160.,text//char(0),     
     +     1405.,85.," ",0.,0.)
      call make_picture_current(fildes)

      if(use_minmax.and.MMOK.eq.3) then
         region1 = 91
      else
         call request_locator
     +        (fildes,1,2e9,valid,x_dc,y_dc,z_dc)
         call vdc_to_wc(fildes,x_dc,y_dc,z_dc,
     +        xloc,yloc,zloc)   !translate into fp coords.
         call locate_region2(xloc,yloc,region1)
      endif
*
      if(region1.eq.91)then
         call clear_all(fildes)
         goto 400
      else                      !select a range of bins
         call clear_quarter(fildes)
*
         if(((yloc.ge.y0_I).and.(yloc.le.y0_I+height)).and.
     +        ((xloc.ge.x0_I-100.).and.(xloc.le.x0_I+width+100.)))then !1st click is in Ti histogram:  change Ti selection range
            min_selected=INT((xloc-x0_I)/pixels_per_bin)+1 !first bin of selected cycles
            if(min_selected.lt.1)min_selected = 1 !allow a pick to the left of the histogram
            if(min_selected.gt.100)then !oops - 1st pick is to the right of the histogram
               call text_color(fildes,1.,0.,0.) !red text
               text='Select the MINIMUM cycle duration first.'
               call text2d(fildes,x0+width+100.,55.,text//char(0),
     +              ANNOTATION_TEXT,0)
               call make_picture_current (fildes)
               call text_color(fildes,0.,0.,0.) !black text
               goto 310
            end if
c             I_MIN=shortest_I_all+(min_selected-1)*binwidth_I          !minimum allowed Ti
            call clear_quarter(fildes) !clear bottom of window
            call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
            call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc)
            if(((yloc.ge.y0_I).and.(yloc.le.y0_I+height)).and. !is 2nd click in the Ti histogram?
     +           ((xloc.ge.x0_I-100.).and.(xloc.le.x0_I+width+100.)))
     +           then
               max_selected=INT((xloc-x0_I)/pixels_per_bin)+1 !yes - last bin of selected cycles
               if(max_selected.lt.min_selected)then
                  call text_color(fildes,1.,0.,0.) !red text
                  text='No, no, no.  Bad user.  Start over and ' //
     +                 'remember to select the MINIMUM duration first.'
                  call text2d(fildes,x0+width+100.,55.,text//char(0),
     +                 ANNOTATION_TEXT,0)
                  call make_picture_current (fildes)
                  call text_color(fildes,0.,0.,0.) !black text
                  goto 310
               end if
               if(max_selected.gt.100)max_selected=100
               if(debug)min_selected = 32
               if(debug)max_selected = 91
               I_MIN=shortest_I_all+(min_selected-1)*binwidth_I !minimum allowed Ti
               I_MAX=shortest_I_all+(max_selected*binwidth_I) !maximum allowed Ti
            else                !2nd click not in same histogram - start over
               call text_color(fildes,1.,0.,0.) !red text
               text='You may use only one histogram at a time to '
     +              //'pick your range. Start over.'
               call text2d(fildes,x0+width+100.,55.,text//char(0),
     +              ANNOTATION_TEXT,0)
               call make_picture_current (fildes)
               call text_color(fildes,0.,0.,0.) !black text
               goto 310
            end if
         else if(((yloc.ge.y0_E).and.(yloc.le.y0_E+height)).and. !1st click is in Te histogram:  change Te selection range
     +           ((xloc.ge.x0_E-100.).and.(xloc.le.x0_E+width+100.)))
     +           then
            min_selected=INT((xloc-x0_E)/pixels_per_bin)+1 !first bin of selected cycles
            if(min_selected.lt.1)min_selected=1 !allow a pick to the left of the histogram
            if(min_selected.gt.100)then !oops - 1st pick is to the right of the histogram
               call text_color(fildes,1.,0.,0.) !red text
               text='Select the MINIMUM cycle duration first.'
               call text2d(fildes,x0+width+100.,55.,text//char(0),
     +              ANNOTATION_TEXT,0)
               call make_picture_current (fildes)
               call text_color(fildes,0.,0.,0.) !black text
               goto 310
            end if
            call clear_quarter(fildes) !clear bottom of window
            call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
            call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc)
            if(((yloc.ge.y0_E).and.(yloc.le.y0_E+height)).and. !is 2nd click in the Te histogram?
     +           ((xloc.ge.x0_E-100.).and.(xloc.le.x0_E+width+100.)))
     +           then
               max_selected=INT((xloc-x0_E)/pixels_per_bin)+1 !yes - last bin of selected cycles
               if(max_selected.lt.min_selected)then
                  call text_color(fildes,1.,0.,0.) !red text
                  text='No, no, no.  Bad user.  Start over and ' //
     +                 'remember to select the MINIMUM duration first.'
                  call text2d(fildes,x0+width+100.,55.,text//char(0),
     +                 ANNOTATION_TEXT,0)
                  call make_picture_current (fildes)
                  call text_color(fildes,0.,0.,0.) !black text
                  goto 310
               end if
               if(max_selected.gt.100)max_selected=100
               if(debug)min_selected = 32
               if(debug)max_selected = 85
               E_MIN=shortest_E_all+(min_selected-1)*binwidth_E !minimum allowed Te
               E_MAX=shortest_E_all+(max_selected*binwidth_E) !maximum allowed Te
            else                !2nd click not in same histogram - start over
               call text_color(fildes,1.,0.,0.) !red text
               text='You may use only one histogram at a time to '
     +              //'pick your range. Start over.'
               call text2d(fildes,x0+width+100.,55.,text//char(0),
     +              ANNOTATION_TEXT,0)
               call make_picture_current (fildes)
               call text_color(fildes,0.,0.,0.) !black text
               goto 310
            end if
         else
            goto 310
         end if
      end if
*
*
*
*       *****   include cycles if the Ti for that cycle falls within the        *****
*       *****           user-entered range:                                     *****
*
 311  begin_E = 0.0             !re-initialize the holding arrays:
      end_E = 0.0
      middle_I = 0.0
      I_durations = 0.0
      E_durations = 0.0
      cycle_durations = 0.0
      sum_of_squares_I=0.0
      sum_of_squares_E=0.0
      sum_of_I = 0.0
      sum_of_E = 0.0
      j = 0                     !initialize the array pointer
      do i=1,INT(CRITERION2_CYCLES)  !for all cycles meeting criteria 1 & 2,
         if((init_I_durations(i).ge.I_MIN).and. !  test for selection based on Ti/Te
     +        (init_I_durations(i).le.I_MAX).and.
     +        (init_E_durations(i).ge.E_MIN).and.
     +        (init_E_durations(i).le.E_MAX))then
            j = j+1
            I_durations(j)=init_I_durations(i) !update holding arrays to 
            E_durations(j)=init_E_durations(i) !  include data pertaining 
            cycle_durations(j)=control_cycle_durations(i) !  only to cycles that
            begin_E(j)=control_begin_E(i) !  meet criterion #3
            end_E(j)=control_end_E(i)
            middle_I(j)=control_middle_I(i)             
         end if
      end do
      TOTAL_SELECTED_CYCLES = j !update # of selected cycles
      longest_I = I_durations(1) !load comparison variables
      shortest_I = I_durations(1)
      longest_E = E_durations(1)
      shortest_E = E_durations(1)
      do i = 1,TOTAL_SELECTED_CYCLES !find longest and shortest Ti 
         longest_I = AMAX1(I_durations(i),longest_I) !  and Te
         shortest_I = AMIN1(I_durations(i),shortest_I)
         sum_of_I = sum_of_I + I_durations(i)
         longest_E = AMAX1(E_durations(i),longest_E)   
         shortest_E = AMIN1(E_durations(i),shortest_E)
         sum_of_E = sum_of_E + E_durations(i)
      end do
      mean_I = sum_of_I / float(TOTAL_SELECTED_CYCLES)
      mean_E = sum_of_E / float(TOTAL_SELECTED_CYCLES)
*
*       *****  calculate standard deviation of mean Ti *****
*
      do i=1,TOTAL_SELECTED_CYCLES  
         sum_of_squares_I = sum_of_squares_I + 
     +        ((I_durations(i)-mean_I)**2)
         sum_of_squares_E = sum_of_squares_E + 
     +        ((E_durations(i)-mean_E)**2)
      end do
      std_dev_I=SQRT(sum_of_squares_I/float(TOTAL_SELECTED_CYCLES-1))
      std_dev_E=SQRT(sum_of_squares_E/float(TOTAL_SELECTED_CYCLES-1))
*
*
*
*
*       *****  Go back to CALCULATE to get new values of mean respiratory cycle *****  
*       *****   duration, standard deviation, etc.:                             *****
*
 400  call clear_all(fildes)
      call CALCULATE (criterion,shortest,longest,TEMP_E,
     +     mean_resp_cycle,
     +     std_dev,IMAXCYC,MAX_CYCLE,MIN_CYCLE,
     +     TEMP_I,TOTAL_NUM_CYCLES,TOTAL_NUM_CONTROL_CYCLES,
     +     cycle_durations,TOTAL_SELECTED_CYCLES,begin_E,
     +     end_E,middle_I,ORIG_DIST_HIST,
     +     control_cycle_durations,control_begin_E,
     +     control_end_E,control_middle_I,
     +     fildes,longest_control,shortest_control,
     +     rec_BW,rec_offset,mean_I,std_dev_I,
     +     x0,y0,
     +     width,bin_duration,control_max_value,
     +     mean_resp_cycle_all,std_dev_all)
*
*
      if(region1.eq.91)then
         call clear1(fildes,0)  !clear the screen
         MM%I_MIN = I_MIN
         MM%I_MAX = I_MAX
         MM%E_MIN = E_MIN
         MM%E_MAX = E_MAX
         call write_minmax (Z(DB_FILES), MM)
         return
      end if
*
      goto 305
*
*
*       ***** RETURN   *****    RETURN  *****   RETURN  *****   RETURN  *****   RETURN  *****
*
*
*       return
      end
*
*
*
*
*
*       **************************************************************************************
*       **************************************************************************************
*
*       SUBROUTINE CALCULATE:
*
*       1st time into this subroutine, all cycles within the boundaried control period
*               are used (criterion #1 has been meet); MIN_CYCLE and MAX_CYCLE = 0.0.
*
*       In subsequent passes, cycles will be selected according to values of MIN_CYCLE
*               and MAX_CYCLE.
*
*       The "control_" arrays hold the initial values of duration, begin_E, end_E,
*               and middle_I (i.e., the cycles that have passed criterion #1 by virtue of
*               their inclusion in the boundaried control period).  This arrangement allows
*               the user complete freedom in selection of MIN and MAX durations (i.e., the
*               selection of cycles is not an ever-tightening noose!).
*
*
      subroutine CALCULATE (criterion,shortest,longest,TEMP_E,
     +     mean_resp_cycle,
     +     std_dev,IMAXCYC,MAX_CYCLE,MIN_CYCLE,
     +     TEMP_I,TOTAL_NUM_CYCLES,
     +     TOTAL_NUM_CONTROL_CYCLES,
     +     cycle_durations,TOTAL_SELECTED_CYCLES,
     +     begin_E,end_E,middle_I,ORIG_DIST_HIST,
     +     control_cycle_durations,control_begin_E,
     +     control_end_E,control_middle_I,
     +     fildes,longest_control,
     +     shortest_control,rec_BW,rec_offset,mean_I,
     +     std_dev_I,x0,y0,
     +     width,bin_duration,control_max_value,
     +     mean_resp_cycle_all,std_dev_all)
*
      use mod_clear_routines
      use mod_miscellaneous_subroutines
      use mod_new_plot
      include 'x2000parameter.defs'
*
      integer TOTAL_NUM_CYCLES,ORIG_DIST_HIST(101),IHIST(101),
     +     TOTAL_SELECTED_CYCLES,TOTAL_NUM_CONTROL_CYCLES,
     +     IHIST1(101),criterion,resp_offset,ICN0(6)
*
      INCLUDE 'gopen_type.defs'
*
      DOUBLE PRECISION
     +     begin_E(:),
     +     end_E(:),
     +     middle_I(:),
     +     control_end_E(:),
     +     control_begin_E(:),
     +     control_middle_I(:),
     +     TEMP_E(:),
     +     TEMP_I(:)
      REAL
     +     cycle_durations(:),
     +     control_cycle_durations(:)

      real longest,shortest,bin_duration,mean_resp_cycle,
     +     mean_resp_cycle_all,std_dev,std_dev_all,
     +     sum_of_durations,sum_of_squares,
     +     MIN_CYCLE,
     +     MAX_CYCLE,longest_control,shortest_control,
     +     rec_BW,rec_offset,mean_I,std_dev_I,x0,y0,height,
     +     width,
     +     control_max_value
* 
      character*120 text
      character*20 c_short_selected,c_short_all,c_long_selected,
     +     c_long_all,c_mean_selected,
     +     c_mean_all,c_std,c_std_all,c_bw
      character*10 c_selected,c_total
      character*3 DIS
      if(.false.)print *,fildes2 !suppress unused variable warning

      ICN0 = 0
      sum_of_squares = 0.0
      sum_of_squares_e = 0.0
      sum_of_squares_i = 0.0

      if(.false.)print *,TOTAL_NUM_CYCLES !suppress unused variable warning
*
*
      DIS = 'DIS'
      height=300.
*
      if(criterion.eq.3)goto 100 !coming to CALCULATE only for statistical analysis
*
*
*       ****************************************************
*       ***** for the first time into this subroutine: *****
*       ****************************************************
*
      if(criterion.eq.1)then
         j = 0                  !initialize the array pointer
         control_cycle_durations = 0.0
         cycle_durations = 0.0  !initialize variable arrays
         ORIG_DIST_HIST = 0
         sum_of_durations = 0.0
*         do i = 1, TOTAL_NUM_CONTROL_CYCLES-1
         do i = 1, TOTAL_NUM_CONTROL_CYCLES
            if(TEMP_I(i).eq.0.0)cycle !eliminate cycles that don't
                                !  have an intermediate I pulse
            j = j + 1           !update the pointer
            control_cycle_durations(j)=TEMP_E(i+1)-TEMP_E(i) !fill control_cycle_durations()
                                ! with the durations of ALL
                                ! cycles contained within the 
                                ! boundaried control period
            control_begin_E(j)=TEMP_E(i) !also, fill arrays with I and E
            control_end_E(j)=TEMP_E(i+1) !  pulses associated with each 
            control_middle_I(j)=TEMP_I(i) !  cycle
         end do
         TOTAL_NUM_CONTROL_CYCLES = j !total number of complete control cycles
                                !  (i.e., have 2 E pulses and 1 I pulse)
         TOTAL_SELECTED_CYCLES = j !total number of cycles selected thus 
         cycle_durations = control_cycle_durations !  far
         begin_E = control_begin_E                             
         end_E = control_end_E                                 
         middle_I = control_middle_I
         longest = cycle_durations(1) !load variables with initial values
         shortest = cycle_durations(1) ! for the sake of comparison
         do i = 1,TOTAL_NUM_CONTROL_CYCLES
            longest=AMAX1(cycle_durations(i),longest) !find the longest cycle duration
            shortest=AMIN1(cycle_durations(i),shortest) !find the shortest cycle duration
            sum_of_durations = sum_of_durations + !calculate total of
     +           cycle_durations(i) ! durations
         end do
         longest_control=longest !store the values of greatest and least
         shortest_control=shortest ! durations of control cycle durations
         bin_duration=(longest+.1-shortest)/100. !this binwidth will be used for the display of
                                ! all respiratory cycle duration distribution 
                                ! histograms -- will allow display of selected
                                ! cycles contrasted with all control cycles
         mean_resp_cycle_all = sum_of_durations /
     +        FLOAT(TOTAL_SELECTED_CYCLES) !mean duration using all control cycles
         do i = 1,TOTAL_NUM_CONTROL_CYCLES
            sum_of_squares=sum_of_squares +
     +           ((cycle_durations(i)-mean_resp_cycle_all)**2)
         end do
         std_dev_all=SQRT(sum_of_squares/
     +        float(TOTAL_NUM_CONTROL_CYCLES-1))
         control_max_value = 
     +        mean_resp_cycle_all + (2*(std_dev_all)) !max cycle length used for case 4 sig. testing
      end if
*
*       *************************************************************
*       *****   for subsequent forays into this subroutine:     *****
*       *************************************************************
*
      if(criterion.eq.2)then
         j=0                    !  so......
         cycle_durations = 0.0  !initialize holding arrays
         begin_E = 0.0
         end_E = 0.0
         middle_I = 0.0
*
*       ***** Go through all cycles that have met criterion #1 and select those *****
*       *****   whose durations fall between MIN and MAX (criterion #2).  Store *****
*       *****   the durations and pulses associated with these selected cycles  *****
*       *****   in holding arrays:                                              *****
*
         CHECK_DURATIONS:  do i=1,TOTAL_NUM_CONTROL_CYCLES     
         if((control_cycle_durations(i).ge.MIN_CYCLE).and. 
     +        (control_cycle_durations(i).le.MAX_CYCLE))then
            j=j+1
            cycle_durations(j)=control_cycle_durations(i)
            begin_E(j)=control_begin_E(i)
            end_E(j)=control_end_E(i)
            middle_I(j)=control_middle_I(i)
         end if
      end do CHECK_DURATIONS
      TOTAL_SELECTED_CYCLES = j !update tally of cycles selected
      end if                    !  thus far
*
*
*       ***** contruct and display the distribution histogram *****
*
*
      IHIST = 0
      BUILD_DISTRIBUTION_HISTOGRAM: 
     +     do i = 1,TOTAL_SELECTED_CYCLES
      ibin=(((cycle_durations(i)-shortest_control)/
     +     bin_duration))+1
      IHIST(ibin)=IHIST(ibin)+1
      end do BUILD_DISTRIBUTION_HISTOGRAM
      IHIST(101)=0
*
      if((MIN_CYCLE.eq.0.0).and.(MAX_CYCLE.eq.0.0))
     +     ORIG_DIST_HIST=IHIST 
      NHW = (bin_duration*100.)
      call clear_all(fildes)
      call new_plot(fildes,100,ORIG_DIST_HIST,IHIST,IHIST1,DIS,x0,y0,   
     +     height,width,'','','','','',NHW,0,'',0,'',ICN0,0,'',
     +     '','',[integer::],'','','',0,0.,[integer::],[integer::])
      write (c_bw,'(f7.1)') bin_duration
      call remove_leading_blanks(c_bw,LEN(c_bw))
      call strlength(c_bw,LEN(c_bw),l_bw)
      write (c_short_all,'(f7.1)') shortest_control
      call remove_leading_blanks(c_short_all,LEN(c_short_all))
      call strlength(c_short_all,LEN(c_short_all),l_sa)
      write (c_long_all,'(f7.1)') longest_control
c      call remove_leading_blanks(c_long_selected,
c     +     LEN(c_long_selected))
c      call strlength(c_long_selected,LEN(c_long_selected),l_ls)
      call character_height(fildes,.040)
      call character_width(fildes,.010)
      call text2d(fildes,x0,y0-20.,c_short_all(1:l_sa)//char(0),
     +     ANNOTATION_TEXT,0)
      call strlength(c_long_all,LEN(c_long_all),l_la)
      call joinstrings(c_long_all(1:l_la),' msec',text,l)
      call text2d(fildes,x0+width-50.,y0-20.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      text='binwidth = '//c_bw(1:l_bw)//' msec.'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0,y0-50.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
*
*
*       ***** proceed with statistical analysis of "good data" *****
*
*
 100  sum_of_durations = 0.0    !initialize variables
      mean_resp_cycle = 0.0
      std_dev = 0.0
      sum_of_squares = 0.0
*
      longest=cycle_durations(1) !load longest & shortest with initial
      shortest=cycle_durations(1) ! values (for comparison)
*
      do i = 1,TOTAL_SELECTED_CYCLES
         longest=AMAX1(cycle_durations(i),longest) !find the longest and shortest
         shortest=AMIN1(cycle_durations(i),shortest) ! durations
         sum_of_durations = sum_of_durations + !calculate total of
     +        cycle_durations(i) ! durations
      end do
*
      if(criterion.eq.1)then
         MIN_CYCLE = shortest
         MAX_CYCLE = longest
      end if
*
      mean_resp_cycle = sum_of_durations /
     +     FLOAT(TOTAL_SELECTED_CYCLES) !mean duration (fp format)
      resp_offset=INT(mean_resp_cycle) !mean duration (integer)
      do i = 1,TOTAL_SELECTED_CYCLES
         sum_of_squares=sum_of_squares +
     +        ((cycle_durations(i)-mean_resp_cycle)**2)
      end do
      std_dev=SQRT(sum_of_squares/float(TOTAL_SELECTED_CYCLES-1))
      IMAXCYC = resp_offset + (2*(std_dev)) !max cycle length used for ANOVA
*
      call character_height(fildes,.060)
      call character_width(fildes,.015)
      write (c_selected,'(I5)') TOTAL_SELECTED_CYCLES
      call remove_all_blanks(c_selected,LEN(c_selected))
      call strlength(c_selected,LEN(c_selected),l)
      write (c_total,'(I5)') TOTAL_NUM_CONTROL_CYCLES
      call remove_all_blanks(c_total,LEN(c_total))
      call strlength(c_total,LEN(c_total),m)
      text=c_selected(1:l)//' of '//c_total(1:m)//
     +     ' control cycles have been selected'
      call text2d(fildes,x0,y0+height+50.,text//char(0),
     +     ANNOTATION_TEXT,0)
      call character_height(fildes,.040)
      call character_width(fildes,.010)
      write (c_bw,'(f7.1)') bin_duration
      call remove_leading_blanks(c_bw,LEN(c_bw))
      call strlength(c_bw,LEN(c_bw),l_bw)
      write (c_short_selected,'(f7.1)') shortest
      call remove_leading_blanks(c_short_selected,
     +     LEN(c_short_selected))
      call strlength(c_short_selected,LEN(c_short_selected),l_ss)
      write (c_short_all,'(f7.1)') shortest_control
      call remove_leading_blanks(c_short_all,LEN(c_short_all))
      call strlength(c_short_all,LEN(c_short_all),l_sa)
      write (c_long_selected,'(f0.1)') longest
      call remove_leading_blanks(c_long_selected,
     +     LEN(c_long_selected))
      call strlength(c_long_selected,LEN(c_long_selected),l_ls)
      write (c_long_all,'(f0.1)') longest_control
      call remove_leading_blanks(c_long_all,LEN(c_long_all))
      call strlength(c_long_all,LEN(c_long_all),l_la)
      write (c_mean_selected,'(f7.1)') mean_resp_cycle
      call remove_leading_blanks(c_mean_selected,
     +     LEN(c_mean_selected))
      call strlength(c_mean_selected,LEN(c_mean_selected),l_ms)
      write (c_mean_all,'(f7.1)') mean_resp_cycle_all
      call remove_leading_blanks(c_mean_all,LEN(c_mean_all))
      call strlength(c_mean_all,LEN(c_mean_all),l_ma)
      write (c_std,'(f7.1)') std_dev
      call remove_leading_blanks(c_std,LEN(c_std))
      call strlength(c_std,LEN(c_std),l_std)
      write (c_std_all,'(f7.1)') std_dev_all
      call remove_leading_blanks(c_std_all,LEN(c_std_all))
      call strlength(c_std_all,LEN(c_std_all),l_std_a)

      text='SELECTED'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+350.,y0+height-15.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      text='CYCLES'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+350.,y0+height-35.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      text='ALL'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+500.,y0+height-15.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      text='CYCLES'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+500.,y0+height-35.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      text='(all times in milliseconds)'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+620.,y0+height-35.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)

      text='shortest cycle length:'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+150.,y0+height-60.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-60.,
     +     c_short_selected(1:l_ss)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-60.,
     +     c_short_all(1:l_sa)//char(0),ANNOTATION_TEXT,0)
      text='longest cycle length:'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+150.,y0+height-80.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-80.,
     +     c_long_selected(1:l_ls)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-80.,
     +     c_long_all(1:l_la)//char(0),ANNOTATION_TEXT,0)

      text='mean cycle length:'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+150.,y0+height-100.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-100.,
     +     c_mean_selected(1:l_ms)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-100.,
     +     c_mean_all(1:l_ma)//char(0),ANNOTATION_TEXT,0)

      text='standard dev:  +/-'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x0+width+150.,y0+height-120.,
     +     text(1:l)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+350.,y0+height-120.,
     +     c_std(1:l_std)//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,x0+width+500.,y0+height-120.,
     +     c_std_all(1:l_std_a)//char(0),ANNOTATION_TEXT,0)

      call make_picture_current (fildes)
*
*
*
*
*       ***** Calculate recommended binwidth and offset for resp CTHs:  *****
*
      calc_BW=((mean_resp_cycle+(2*std_dev)) +
     +     ((mean_I+(2*std_dev_I))/2))/100
*
*               ***** make sure that rec_BW is in *.0 format and is a   *****
*               *****   multiple of 5:                                  *****
* 
      if(INT(AMOD(calc_BW,5.0)).eq.1)then !INT(AMOD(calc_BW,5.0))=1
         rec_BW=FLOAT(INT(calc_BW+4.0)) 
      elseif(INT(AMOD(calc_BW,5.0)).eq.2)then !INT(AMOD(calc_BW,5.0))=2
         rec_BW=FLOAT(INT(calc_BW+3.0)) 
      elseif(INT(AMOD(calc_BW,5.0)).eq.3)then !INT(AMOD(calc_BW,5.0))=3
         rec_BW=FLOAT(INT(calc_BW+2.0)) 
      elseif(INT(AMOD(calc_BW,5.0)).eq.4)then !INT(AMOD(calc_BW,5.0))=4
         rec_BW=FLOAT(INT(calc_BW+1.0)) 
      elseif(INT(AMOD(calc_BW,5.0)).eq.0)then !INT(AMOD(calc_BW,5.0))=0
         rec_BW=calc_BW
      end if
*
      rec_offset = -(mean_I+(2*std_dev_I))
      rec_offset = ANINT (rec_offset * 2.) / 2.
*
*
      RETURN
      END
      end module mod_RESP_CYCLE_STATS_3_v4
