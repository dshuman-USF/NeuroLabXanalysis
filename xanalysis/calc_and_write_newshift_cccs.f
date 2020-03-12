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


      module mod_calc_and_write_newshift_cccs
      contains
*       filename = calc_and_write_newshift_cccs.f

*       19-Feb-2004      lss
*               Shift-predictor (SP) control CCHs now calculated as follows:
*                       the reference spike train is shifted, in turn, into 20 different cycles,
*                       providing 20 SP histograms which will then be averaged into the SP for
*                       that pair of cells at that binwidth.  
*               Events which occur within a non-acceptable respiratory cycle (as determined by
*                       a user-entered maximum allowable cycle length) are not used for CCHs or SPs

*       29-may-2003     lss
*               version cccs of calc_and _write_newshift uses a subset of the data file
*               to compute the CCH and CONTROL histograms for each cell pair; the subset
*               is defined as those spikes which occur during acceptable respiratory cycles
*
*       may-2003        lss
*               new method of shifting REF train for calculation of shift-control:
*               shift is now proportional, by one respiratory cycle (as defined by E pulses)
*
*       may-2002        lss
*               change cardiac CTHs to CCHs; no offset
*
*       17-may-2002     lss
*               calculate and display rates for CTHs
*
*       apr-2002        lss
*               allow import of cell coordinate data
*               start preparing for ddt files!
*
*       apr-2001        lss
*        modified to allow user to select respiratory cycles for CTH calculation:
*               A -- use all cycles in the recording
*               F -- use first 10000 (=MAX_CYCLES) cycles in the recording
*               B -- use all cycles within a boundary coded section of the recording
*               C -- use the first x cycles in the recording
*               T -- designate a period of time within the recording
*
*       21-mar-2000     lss
*        version numbers included in *.qdt file to avoid
*          analysis using a version different from the one
*          used to generate the histograms
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
*       *** DIRECT POINTERS are used with these arrays: excluded, included
*       ***     (these arrays are dimensioned to MAX_NUM_CODES)
*
*       link with x2002 code
*
*       NOTE:  THIS VERSION CONTAINS DIRECT ACCESS PARAMETER FOR QDT FILE
*       
*       This file contains code to be used to create the .qdt file.  The user is asked
*               to enter date of experiment, recording #, binwidths, etc.  
*               Database and text files are created with a different subroutine.
*
*       Format of .qdt file:    [# cells = N; # pairs = (N(N-1)/2)]
*
*               histogram #
*                       1       phrenic overlay                 total # overlays = 3
*                       2       normalized phrenic overlay
*                       3       cardiac overlay
*               -----------------------------------------------------
*                       4       CTHs for 1st cell:      respiratory CTH
*                       5                               norm. resp. CTH
*                       6                               cardiac CTH
*                       .
*                       .                                       total # CTHs = 3N
*                       .
*                    3N+3       cardiac CTH for cell N
*               -------------------------------------------------------
*                    3N+4       ACHs for 1st cell:      @ binwidth 1
*                    3N+5                               @ binwidth 2
*                    3N+6                               @ binwidth 3
*                    3N+7                               @ binwidth 4
*                       .
*                       .                                       total # ACHs = 4N
*                       .
*                    7N+3       ACH for cell N @ binwidth 4
*               -------------------------------------------------------
*                    7N+4       CCHs for 1st pair:      @ binwidth 1
*                    7N+5                               @ binwidth 2
*                    7N+6                               @ binwidth 3
*                    7N+7                               @ binwidth 4
*                       .
*                       .                                       total # CCHs = 2(N**2)-2N
*                       .
*            2(N**2)+5N+3       CCH for last pair @ binwidth 4
*               -------------------------------------------------------
*            2(N**2)+5N+4  1-shiftcontrol CCHs for 1st pair:    @ binwidth 1  (REF train shifted into the next cycle)
*            2(N**2)+5N+5                                       @ binwidth 2
*            2(N**2)+5N+6                                       @ binwidth 3
*            2(N**2)+5N+7                                       @ binwidth 4
*                       .
*                       .                                       total # 1-shift control CCHs = 2(N**2)-2N
*                       .
*            4(N**2)+3N+3       control CCH for last pair @ binwidth 4
*               -------------------------------------------------------
*            4(N**2)+3N+4  2-shiftcontrol CCHs for 1st pair:    @ binwidth 1  (REF train shifted by 2 cycles)
*            4(N**2)+3N+5                                       @ binwidth 2
*            4(N**2)+3N+6                                       @ binwidth 3
*            4(N**2)+3N+7                                       @ binwidth 4
*                       .
*                       .                                       total # 2-shift control CCHs = 2(N**2)-2N
*                       .
*            6(N**2)+N+3        control CCH for last pair @ binwidth 4
*               -------------------------------------------------------
*            6(N**2)+N+4  20-shiftcontrol CCHs for 1st pair:    @ binwidth 1  (save the average of these 20 shifts)
*            6(N**2)+N+5                                        @ binwidth 2
*            6(N**2)+N+6                                        @ binwidth 3
*            6(N**2)+N+7                                        @ binwidth 4
*                       .
*                       .                                       total # 20-shift control CCHs = 2(N**2)-2N
*                       .
*            8(N**2)-N+3        control CCH for last pair @ binwidth 4
*               -------------------------------------------------------
*                                
*               total # of histograms written to .qdt = 8(N**2) - N + 3 = total_histograms
*
*
*

      subroutine calculate_histograms(version,mode,
     +     BDT_FILE,
     +     qdt_files,SPIKETIMES,IDs,ITAL,abort,c_format,
     +     total_num_qdts,BWs,date_exp,
     +     recording,protocol,exp_name,DB_FILES)
*
      use mod_RESP_CYCLE_STATS_3_v4
      use mod_calcCCH_newshift_2
      use mod_calcCTH
      use mod_calc_ach_v4
      use mod_calc_cardCCH
      use mod_cardsig_delta2
      use mod_clear_routines
      use mod_intvar4_v2a
      use mod_mean_and_sd_routines
      use mod_miscellaneous_subroutines
      use mod_qdttxt_write
      use mod_read_and_write_DBSAV
      use mod_respsig6_v4
      use mod_scalcn
      use mod_serotonergic_discr_func
      use mod_showCTHs_v3
      include 'x2000parameter.defs'
*
      integer ICN(6),jTOTAL_SHIFTS
*
      real, allocatable ::
     +     cycle_durations(:)
      double precision, allocatable :: 
     +     begin_E(:),
     +     norm_begin_E(:),
     +     end_E(:),
     +     middle_I(:)

      double precision
     +     REFtime,diff,shift1,shift2,
     +     BNDRY_START_TIME,BNDRY_END_TIME,first_E,last_E,
     +     MAX_INT,factor,
     +     first_spike_time,last_spike_time,elapsed_time,
     +     shifted,SJOFSET,CSTDDEV
      double precision
     +     SPIKETIMES(:,:)
      double precision, allocatable ::
     +     TEMP_durations(:),
     +     E_begin(:),
     +     E_end(:),
     +     SPIKESHIFT(:,:),
     +     cycle_length(:),
     +     I_assoc(:),
     +     Te(:),
     +     Ti(:),
     +     selected_C(:),
     +     RTAL(:)
      double precision SURR_TIMES(size(spiketimes,1),MAX_NUM_CHAN)

      POINTER(STP,SURR_TIMES)
      integer IHIST(NUM_BINS),excluded(MAX_NUM_CODES),
     +     CONTROL(MAX_NUM_QDTS,4,NUM_BINS),
     +     included(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN),
     +     IDs(MAX_NUM_CODES),SURR_ITAL(MAX_NUM_CHAN),
     +     TOTAL_SELECTED_CYCLES,
     +     phrenic_hist(NUM_BINS),
     +     norm_phrenic_hist(NUM_BINS),cardiac_hist(NUM_BINS),
     +     width,height,STATUS,total_num_qdts,pair_num,
     +     SP(NUM_BINS),
     +     sp_per_cycle(MAX_NUM_CHAN,MAX_NUM_ACC_CYCLES),
     +     Q_pos,rel_loc,
     +     REC_num_control_1,REC_num_control_2,REC_num_control_avg,
     +     SPC_array,set_A(4),set_B(4),savedITAL(MAX_NUM_CODES),
     +     doomed(MAX_NUM_CODES),ebins,ibins,spike_count,
     +     cells(MAX_NUM_CODES),markers(MAX_NUM_CODES)

      integer recnum(total_num_qdts)
*
      double precision
     +     mean_resp_cycle_all,std_dev_all

      real BINW,BINW2,BINWa,BINWx,BINWs,coefvr,STIM_OFFSET,NORM_BW,
     +     mean_resp_cycle,
     +     rec_BW,rec_offset,NORM_OFFSET,BINW_1,BINW_2,BINW_3,BINW_4,
     +     coefnum(MAX_NUM_CHAN),
     +     mean_control_cycle,MAXCYCLE,user_max_value,
     +     MAXduration,BWs(MAX_NUM_QDTS,4),
     +     mean_resp_cycle_control,
     +     longest_control,shortest_control,shortest,longest,
     +     std_dev_control,longest_all,shortest_all,mean_E

      character*(*) version,mode,BDT_FILE,qdt_files(MAX_NUM_QDTS),
     +     date_exp,recording,protocol,exp_name,DB_FILES

      character*10 fiveHT(MAX_NUM_CHAN),meanISI(MAX_NUM_CHAN),
     +     sdISI(MAX_NUM_CHAN),mean_rISI(MAX_NUM_CHAN),
     +     sd_rISI(MAX_NUM_CHAN),num_rej_ISI(MAX_NUM_CHAN),
     +     num_rej_rISI(MAX_NUM_CHAN),c_MAX_INT
*
      integer REF,TAR,E_pulse, cardiac_pls, phrenic, BNDRY, 
     +     flag,change_it,cardiac_offset,LAST_REF,
     +     LAST_TAR,histogram_number,total_histograms,
     +     flagA,flagB,
     +     qdttxt_count,total_num_cells,prevREFcode,
     +     cell,previous_percent,count, total_sel_C_pulses,
     +     spikecount(iTOTAL_SHIFTS+1)
*
      integer*4 mouse
      INCLUDE 'gopen_type.defs'
*
*
      character*250 qdtstring,MEAN__E
      character*120 text
      character*300 text300
      character*60 QDTSAV,QDTTXT
      character*50 QDT_FILENAME
      character*30 perimeter
*
      CHARACTER*20 text1,c_format
*
      character*30 TODAY
      character*20 cntrl
      logical cycle_shift, cth_cch, surrogate
*
*
*
      character*8 CSTDDEV_text,c_mean_all,c_std_all,
     +     c_short_all,c_long_all,c_mean_con,c_std_con,
     +     c_short_con,c_long_con
*
      character*12 bwtext1,bwtext2,bwtext3
      character*7bwtext
*
      character*5 coef(MAX_NUM_CHAN),
     +     ETA2_1(MAX_NUM_CHAN),
     +     ETA2_2(MAX_NUM_CHAN),
     +     ETA2_3(MAX_NUM_CHAN),
     +     ETA2_4(MAX_NUM_CHAN),
     +     ETA2_5(MAX_NUM_CHAN),
     +     ETA2_6(MAX_NUM_CHAN),
     +     cardiac_offset_text,
     +     ICMAXCYC_text,DELTA2(MAX_NUM_CHAN)
*
      CHARACTER*3 card_type(MAX_NUM_CHAN)
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
     +     zmodsig2_6(MAX_NUM_CHAN),file_ext

      character*2 tedfactor(MAX_NUM_CHAN),envtxt
*
      CHARACTER*1 OK, last_one,
     +     BOUNDARY,abort,calc_shift,C_pulse
*
c        SAVE included, excluded
      if(.false.)print *,fildes2 !suppress unused variable warning

      max_num_events = size (spiketimes, 1)

      allocate (TEMP_durations(MAX_NUM_EVENTS))
      allocate (E_begin(MAX_NUM_EVENTS))
      allocate (E_end(MAX_NUM_EVENTS))
      allocate (SPIKESHIFT(MAX_NUM_EVENTS,iTOTAL_SHIFTS+1))
      allocate (cycle_length(MAX_NUM_EVENTS))
      allocate (I_assoc(MAX_NUM_EVENTS))
      allocate (Te(MAX_NUM_EVENTS))
      allocate (Ti(MAX_NUM_EVENTS))
      allocate (selected_C(MAX_NUM_EVENTS))
      allocate (RTAL(MAX_NUM_EVENTS))
*       
*
*
*       *****************************************************************
*       *****************************************************************
*
*
c$$$    character*40 geometry
c$$$    character*6 c_width,c_height
c$$$    call GET LOG(USER)
c$$$    WINDOW1='XAnalysis'
*
*       ***** create the graphics display window: *****
*

      call screensize(width,height,STATUS)
c$$$     if(STATUS.eq.0) then
c$$$       print '(''cannot open display'')'
c$$$       return
c$$$       end if
c$$$     write (c_width,'(I6)') width
c$$$         height=height*0.5
c$$$     write (c_height,'(I6)') height
c$$$     geometry=' '
c$$$     geometry=c_width//'x'//c_height//'-700-5'
c$$$     ito=0
c$$$     do i=LEN(geometry),1,-1                        !remove blanks
c$$$       if(geometry(i:i).ne.' ')exit
c$$$       end do
c$$$     do ifrom=1,i
c$$$       if(geometry(ifrom:ifrom).ne.' ')then
c$$$         ito=ito+1
c$$$         geometry(ito:ito)=geometry(ifrom:ifrom)
c$$$         end if
c$$$       end do
c$$$     geometry((ito+1):LEN(geometry))=' '
c$$$     do i=LEN(geometry),1,-1
c$$$       if(geometry(i:i).ne.' ')exit
c$$$        end do
c$$$        call strlength(USER,LEN(USER),l)
c$$$        isys=SYSTEM ('xwcreate -wmdir /dev/screen/'//USER(1:l)//
c$$$     +          ' -title '//WINDOW1//
c$$$     +          ' -geometry ='//geometry(1:i)//char(0))
c$$$
c$$$        fildes=gopen ('/dev/screen/'//USER(1:l)//'/'//
c$$$     +              WINDOW1//char(0),
c$$$     +      OUTINDEV,'sox11'//char(0),INIT)
      shift2 = 0
      fildes=gopen(width,height/2,-700,-5,'XAnalysis'//char(0))
*
*
*       *****   set screen and window parameters  *****
*
      call vdc_extent(fildes,0.0,0.0,0.0,1.0,1.0,0.0)
      call mapping_mode(fildes,1)
      call shade_mode(fildes,IOR(INIT,CMAP_NORMAL),0)
      call background_color(fildes,1.,1.,1.) !window background = white
      call clear(fildes)
      call view_port(fildes,0.,0.,1.,1.)
      call view_window(fildes,1.,1.,1700.,600.)
      call text_font_index(fildes,6)
      call text_color(fildes,0.,0.,0.)
      call line_color(fildes,0.,0.,0.)
      call perimeter_color(fildes,0.,0.,0.)
*
      call vdc_extent(fildes,0.0,0.0,0.0,1.0,1.0,0.0)
      call view_port(fildes,0.,0.,1.,1.)
      call track(fildes,fildes,1)
      call echo_type(fildes,0,1,.625,.1,0.0)
      call set_locator(fildes,1,.625,.1,0.0)

c        do i = 1, 1000000000
c           call clear_quarter(fildes)
c        end do

*
*
*       * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                               *
*       *       initialize variables and arrays:        *
*       *                                               *
*       * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
      do j = 1, MAX_NUM_ACC_CYCLES
         do i = 1, MAX_NUM_CHAN
            sp_per_cycle (i,j)= 0
         end do
      end do

      SPIKESHIFT = 0.0

      jTOTAL_SHIFTS = 0
      BINW=0.0
      BINW_1=0.0
      BINW_2=0.0
      BINW_3=0.0
      BINW_4=0.0
      rec_BW=0.0
      rec_offset=0.0
      do i = 1, NUM_BINS
         IHIST(i)=0             !array that holds bin values of histograms from subs
      end do
      NUM_STIM =0               !number of stimuli used to compute CTHs
      STIM_OFFSET=0.0           !offset applied to resp CTHs
      NORM_OFFSET=0.0           !offset applied to normalized resp CTHs
      mean_control_cycle=0.0
      histogram_number=0
      total_histograms=0
      previous_percent = 0
      MAXCYCLE=0.0
      MAX_INT = 0.0
      user_max_value=0.0
      BNDRY_START_TIME=0.0
      BNDRY_END_TIME=0.0
      perimeter='perimeter'
      change_it=0
      do i = 1, MAX_NUM_CODES
         excluded(i)=0
         included(i)=0
         savedITAL(i)=0
      end do
      qdttxt_count=0
      OK=' '
      abort='n'
      do i = 1, MAX_NUM_CHAN
         DELTA2(i)= ' '
         tedfactor(i)=' '
         card_type(i) = ' '
         fiveHT (i)= ' '
         meanISI(i) = ' '
         sdISI(i) = ' '
         mean_rISI(i) = ' '
         sd_rISI(i) = ' '
         num_rej_ISI(i) = ' '
         num_rej_rISI(i) = ' '
         coef(i) = ' '
         coefnum(i) = 0.0
      end do
*     ***** set defaults for CTH parameter sets: *****
      set_A(1) = 97             !I pulse
      set_A(2) = 98             !E pulse
      set_A(3) = 89             !phrenic
      set_A(4) = 99             !cardiac pulse
      set_B(1) = 97             !I pulse
      set_B(2) = 98             !E pulse
      set_B(3) = 90             !phrenic
      set_B(4) = 99             !cardiac pulse
      flagA = 0
      flagB = 0
      do i = 1, 4
         call CHECK_CODE(set_A(i),flagA,IDs,'n')
         if(flagA.eq.1)exit     !if flagA = 1 --> one of the IDcodes for set A does not exist
      end do
      do i = 1, 4
         call CHECK_CODE(set_B(i),flagB,IDs,'n')
         if(flagB.eq.1)exit     !if flagB = 1 --> one of the IDcodes for set B does not exist
      end do
*
*
*
*               ***** open "gamesave" file: *****
*

      QDT_FILENAME = qdt_files(1)


*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                       *
*       *   OBTAIN DATA THAT IS COMMON TO ALL CELLS IN GROUP:   * 
*       *                                                       *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
 50   if((flagA.ne.1).and.(flagB.eq.1))
     +     print '(30(/),
     +''Available parameter set for CTH calculation:'',
     +//,                T24,''  A'', 
     +/,                T24,''----'',
     +//,T10,''I pulse'',T24,I3, 
     +/,T10,''E pulse'',T24,I3, 
     +/,T10,''phrenic'',T24,I3, 
     +/,T10,''cardiac pulse'',T24,I3,
     +///,T5,''CHOOSE:'',T15,''A -- set A'',
     +/,               T15,''C -- Custom parameter set'',
     +//,T18,''>> '',$)',(set_A(i),i=1,4)
      if((flagA.eq.1).and.(flagB.ne.1))
     +     print '(30(/),
     +''Available parameter set for CTH calculation:'',
     +//,                T24,''  B'',
     +/,                T24,''----'',
     +//,T10,''I pulse'',T24,I3, 
     +/,T10,''E pulse'',T24,I3, 
     +/,T10,''phrenic'',T24,I3, 
     +/,T10,''cardiac pulse'',T24,I3,
     +///,T5,''CHOOSE:'',T15,''B -- set B'',
     +/,               T15,''C -- Custom parameter set'',
     +//,T18,''>> '',$)',(set_B(i),i=1,4)
      if((flagA.ne.1).and.(flagB.ne.1))
     +     print '(30(/),
     +''Available parameter sets for CTH calculation:'',
     +//,                T25,''A'',  T35,''B'',
     +/,                T24,''---'',T34,''---'',
     +//,T10,''I pulse'',T25,''97'', T35,''97'',
     +/,T10,''E pulse'',T25,''98'', T35,''98'',
     +/,T10,''phrenic'',T25,''89'', T35,''90'',
     +/,T10,''cardiac pulse'',T25,''99'', T35,''99'',
     +///,T5,''CHOOSE:'',T15,''A -- set A'',
     +/,               T15,''B -- set B'',
     +/,               T15,''C -- Custom parameter set'',
     +//,T18,''>> '',$)'
      if((flagA.eq.1).and.(flagB.eq.1))then !neither parameter set exist for the data - user must enter
         OK = 'C'
         goto 54
      end if

      OK = ' '
      read (*,fmt='(A1)',err=50) OK
      call upper_case(OK,LEN(OK))
 54   if(OK.eq.'A')then
         if(flagA.eq.1)goto 50
         I_pulse = set_A(1)
         E_pulse = set_A(2)
         phrenic = set_A(3)
         cardiac_pls = set_A(4)
         C_pulse = 'y'
      else if(OK.eq.'B')then
         if(flagB.eq.1)goto 50
         I_pulse = set_B(1)
         E_pulse = set_B(2)
         phrenic = set_B(3)
         cardiac_pls = set_B(4)
         C_pulse = 'y'
      else if(OK.eq.'C')then
         print '(30(/),T5,
     +''Please enter the following information for ''
     +''CTH calculation:'')'
 55      PRINT '(/,T10,''I pulse ID code  >> '',$)'
         read (*,fmt='(I5)',err=55) I_pulse
         call CHECK_CODE(I_pulse,flag,IDs,'y')
         if(flag.eq.1)goto 55
 60      PRINT '(T10,''E pulse ID code  >> '',$)'
         read (*,fmt='(I5)',err=60) E_pulse
         call CHECK_CODE(E_pulse,flag,IDs,'y')
         if(flag.eq.1)goto 60
 70      PRINT '(T10,''Phrenic ID code  >> '',$)'
         read (*,fmt='(I5)',err=70) phrenic
         call CHECK_CODE(phrenic,flag,IDs,'y')
         if(flag.eq.1)goto 70
 72      PRINT '(T10,''Cardiac pulse ID code ..or.. ''
     +''<cr> if cardiac pulse is not available  >> '',$)'
         cardiac_pls = 0
         flag = 0
         read (*,fmt='(I5)',err=72) cardiac_pls
         if(cardiac_pls.eq.0)then
            C_pulse = 'n'
         else
            C_pulse = 'y'
            call CHECK_CODE(cardiac_pls,flag,IDs,'y')
         end if
         if(flag.eq.1)goto 72
      else
         goto 50
      end if
*
      max_cycles = ITAL(IDs(E_pulse))-1

      allocate (cycle_durations(max_cycles))
      allocate (begin_E(max_cycles))
      allocate (norm_begin_E(max_cycles))
      allocate (end_E(max_cycles))
      allocate (middle_I(max_cycles))

      cycle_durations=0.0
      norm_begin_E=0.0
      begin_E=0.0
      end_E=0.0

 80   icycles=0
      start_time=0.0
      end_time=0.0
      BNDRY_START_TIME=0.0
      BNDRY_END_TIME=0.0
      BNDRY = 0
      print '(///,''Xanalysis can handle a maximum of '',I6,
     +'' cycles for calculation of CTHs.  '',//,
     +''There are '',I6,'' cycles in this file.'',///,
     +T5,''Select method of defining a control period.'',//,
     +T8,''ENTER:'',//,
     +T15,''A --> use All respiratory cycles'',/,
     +T15,''F --> use First '',I6,'' respiratory cycles'',/,
     +T15,''B --> use a Boundary code'',/,
     +T15,''C --> designate the number of Cycles to be used'',/,
     +T15,''T --> designate the period of Time to be used'',//,
     +T15,''X --> eXit to main menu'',//,
     +T20,''>>  '',$)',
     +     MAX_CYCLES,ITAL(IDs(E_pulse))-1,MAX_CYCLES
      read (*,'(A)',err=80),BOUNDARY
      call upper_case(BOUNDARY,LEN(BOUNDARY))
*
 804  if(BOUNDARY.eq.'X')then
         abort='y'
         call clear(fildes)
         goto 500
*
      else if(BOUNDARY.eq.'A')then
 805     icycles=ITAL(IDs(E_pulse))-1 !total # cycles = (# E pulses - 1)
         if(icycles.gt.MAX_CYCLES)then
            print '(/,T5,''Error: Xanalysis can handle a ''
     +''maximum of '',
     +I6,'' cycles for CTHs!  This file contains ''
     +''approximately'',I6,
     +'' cycles.'',
     +/,T5,''Use first '',I6,'' cycles? (y/n)  >> '',$)',
     +           MAX_CYCLES,ITAL(IDs(E_pulse))-1,MAX_CYCLES
            read (*,'(A)',err=805),OK
            call upper_case(OK,LEN(OK))
            if(OK.eq.'Y')then
               BOUNDARY='F'
               goto 804
            else if(OK.eq.'N')then
               goto 80
            else
               goto 805
            end if
         end if
         BNDRY_START_TIME=SPIKETIMES(1,IDs(E_pulse)) !define boundary times 
         BNDRY_END_TIME=SPIKETIMES(icycles+1,IDs(E_pulse)) !define boundary times 
*
      else if(BOUNDARY.eq.'F')then
 807     icycles=MAX_CYCLES
         if(icycles.gt.ITAL(IDs(E_pulse))-1)then
            print '(/,T5,''Warning: this file does not contain '',
     +''that many cycles!  Use all '',I6,'' cycles'',
     +'' in this file? (y/n)  >> '',$)',
     +           ITAL(IDs(E_pulse))-1
            read (*,'(A)',err=807),OK
            call upper_case(OK,LEN(OK))
            if(OK.eq.'Y')then
               BOUNDARY='A'
               goto 804
            else if(OK.eq.'N')then
               goto 80
            else
               goto 807
            end if
         end if
         BNDRY_START_TIME=SPIKETIMES(1,IDs(E_pulse)) !define boundary times 
         BNDRY_END_TIME=SPIKETIMES(icycles+1,IDs(E_pulse)) !define boundary times 
*
      else if(BOUNDARY.eq.'B')then
 81      print '(/,T10,''Enter boundary ID code ..or.. ''
     +''<cr> to return to previous menu  >> '',$)'
         BNDRY_START_TIME=0.0
         BNDRY_END_TIME=0.0
         read (*,fmt='(I5)',err=81) BNDRY
         if(BNDRY.eq.0)goto 80
         call CHECK_CODE(BNDRY,flag,IDs,'y')
         if(flag.eq.1)goto 81
         BNDRY_START_TIME=SPIKETIMES(1,IDs(BNDRY)) !define boundary times for
         BNDRY_END_TIME=SPIKETIMES(2,IDs(BNDRY)) ! respiratory CTHs
         icycles=0
         do i = 1,ITAL(IDs(E_pulse)) !count the # of cycles encompassed
            if((SPIKETIMES(i,IDs(E_pulse)).ge.BNDRY_START_TIME)
     +           .and.          !  by the boundary code
     +           (SPIKETIMES(i,IDs(E_pulse)).le.BNDRY_END_TIME))
     +           icycles=icycles+1
            if(SPIKETIMES(i,IDs(E_pulse)).gt.BNDRY_END_TIME)exit
         end do
         if(icycles.gt.MAX_CYCLES)then
            print '(//,T5,''Boundary code '',I3,
     +'' includes more than '',I6,
     +'' cycles.'')',BNDRY,MAX_CYCLES
            goto 80
         end if
*
      else if(BOUNDARY.eq.'C')then
 82      print '(//,T5,''Total number of cycles in file = '',I6)',
     +        ITAL(IDs(E_pulse))-1
         print '(/,T5,''Enter number of cycles to be used ..or.. ''
     +''<cr> to return  >> '',$)'
         read (*,'(I4)',err=82),icycles
         if(icycles.eq.0)goto 80 !return to previous menu
         if(icycles.gt.MAX_CYCLES)then
            if(ITAL(IDs(E_pulse))-1.gt.MAX_CYCLES)then
 801           print '(/,T5,''Error: Xanalysis can handle a ''
     +''maximum of '',
     +I6,'' cycles!  This file contains '',
     +I6,''cycles.'',
     +/,T12,''Use first '',I6,
     +'' cycles? (y/n)  >> '',$)',
     +              MAX_CYCLES,ITAL(IDs(E_pulse))-1,MAX_CYCLES
               read (*,'(A)',err=801),OK
               call upper_case(OK,LEN(OK))
               if(.false.)then
*
               else if(OK.eq.'Y')then
*                icycles=MAX_CYCLES
                  BOUNDARY='F'
                  goto 804
               else if(OK.eq.'N')then
                  goto 80
               else
                  goto 801
               end if
            end if
            if(ITAL(IDs(E_pulse))-1.lt.MAX_CYCLES)then
 802           print '(/,T5,''Error: Xanalysis can handle a ''
     +''maximum of '',
     +I6,'' cycles!  This file contains '',
     +I6,'' cycles.'',
     +/,T12,''Use all '',I6,'' cycles? (y/n)  >> '',$)',
     +              MAX_CYCLES,ITAL(IDs(E_pulse))-1,ITAL(IDs(E_pulse))-1
               read (*,'(A)',err=802),OK
               call upper_case(OK,LEN(OK))
               if(.false.)then
*
               else if(OK.eq.'Y')then
                  BOUNDARY='A'
                  goto 804
               else if(OK.eq.'N')then
                  goto 80
               else
                  goto 802
               end if
            end if
         end if
         if(icycles.gt.ITAL(IDs(E_pulse))-1)then
 800        print '(/,T5,''Warning: this file does not contain '',
     +''that many cycles!  Use all '',I6,'' cycles'',
     +'' in this file? (y/n)  >> '',$)',
     +           ITAL(IDs(E_pulse))-1
            read (*,'(A)',err=82),OK
            call upper_case(OK,LEN(OK))
            if(.false.)then
*
            else if(OK.eq.'Y')then
               BOUNDARY='A'
               goto 804
            else if(OK.eq.'N')then
               goto 80
            else
               goto 800
            end if
         end if
         BNDRY_START_TIME=SPIKETIMES(1,IDs(E_pulse)) !define boundary times 
         BNDRY_END_TIME=SPIKETIMES(icycles+1,IDs(E_pulse)) !define boundary times 
*
      else if(BOUNDARY.eq.'T')then
         first_spike_time=SPIKETIMES(1,IDs(E_pulse))/60000.
         last_spike_time=
     +        SPIKETIMES(ITAL(IDs(E_pulse)),IDs(E_pulse))/60000.
         elapsed_time=last_spike_time - first_spike_time
 83      print '(///,T5,''Total elapsed time of file = '',F7.3,
     +'' minutes.'',/,T8,''(first spike time: '',f7.3,
     +''; last spike time: '',f7.3,'' min.)'')',
     +        elapsed_time,first_spike_time,last_spike_time
         print '(/,T5,''Enter time period to be used ''
     +''(minutes; fmt=f5.1)  >> '',
     +//,T8,''start time  >> '',$)'
         read (*,'(f5.1)',err=83),start_time
         if(start_time.lt.first_spike_time)then
            print '(20(/)T5,''You can''''t start before the first ''
     +''time!  Please re-enter your start time.'')'
            goto 83
         end if
         if(start_time.ge.last_spike_time)then
            print '(20(/)T5,''You can''''t start after the finish ''
     +''line!  Please re-enter your start time.'')'
            goto 83
         end if
         print '(/,T8,''end time  >> '',$)'
         read (*,'(f5.1)',err=83),end_time
         if(end_time.eq.0.0)goto 80 !return to previous menu
         if(end_time.lt.start_time)goto 83 !end time must be greater than start time
         if(end_time.gt.last_spike_time)then
            print '(20(/),T5,''Error: time entered is greater ''
     +''than the file length!  Try again.'')'
            goto 83
         end if
         BNDRY_START_TIME=start_time*60000. !convert min --> msec
         BNDRY_END_TIME=end_time*60000. !convert min --> msec
         do i = 1,ITAL(IDs(E_pulse)) !count the # of cycles encompassed
            if((SPIKETIMES(i,IDs(E_pulse)).ge.BNDRY_START_TIME).and. !  by the time period
     +           (SPIKETIMES(i,IDs(E_pulse)).le.BNDRY_END_TIME))
     +           icycles=icycles+1
            if(SPIKETIMES(i,IDs(E_pulse)).gt.BNDRY_END_TIME)exit
         end do
         if(icycles.gt.MAX_CYCLES)then
            print '(/,T5,''The time period you entered includes ''
     +''more than '',I6,'' events.'')',MAX_CYCLES
            goto 80
         end if
*
      else
         goto 80
*
      end if
*
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                               *
*       *  Go to subroutine RESP_CYCLE_STATS to select which            *
*       *   respiratory cycles will be used to compute regular CTHs     *
*       *   and normalized CTHs.  Selection is based upon whether       *
*       *   a cycle occurs within a user-selected period (of time or    *
*       *   number of cycles) (REQUIRED),                               *
*       *   its duration is neither too short nor too long (OPTIONAL),  *
*       *   and the duration of its inspiratory and expiratory          *
*       *   phases are within acceptable limits (OPTIONAL).             *
*       *                                                               *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
      call RESP_CYCLE_STATS(SPIKETIMES,I_pulse,E_pulse,
     +     BOUNDARY,icycles,
     +     BNDRY_START_TIME,
     +     BNDRY_END_TIME,phrenic,ITAL,IHIST,
     +     fildes,mouse,mean_resp_cycle,std_dev,IMAXCYC,
     +     begin_E,end_E,TOTAL_SELECTED_CYCLES,
     +     cycle_durations,rec_BW,rec_offset,mean_control_cycle,
     +     control_max_value,IDs,mean_resp_cycle_control,
     +     std_dev_control,longest_control,shortest_control,mean_E,
     +     shortest,longest,DB_FILES,middle_I)


      rec_offset = rec_offset-100.
      BINW=rec_BW               !default values for binwidths and offsets
      STIM_OFFSET = rec_offset
      NORM_BW=rec_BW            !the recommended BW and offset are used to calculate the normalized CTHs
*      NORM_OFFSET=rec_offset
      NORM_OFFSET= ANINT (rec_offset / NORM_BW) * NORM_BW

      call phase_bins (begin_E, middle_I, end_E, NORM_BW,
     +     TOTAL_SELECTED_CYCLES, ebins, ibins)
      mean_E = ebins * NORM_BW

*
*
*       *****  Display the recommended binwidth and offset for CTHs as  *****
*       *****   calculated in RESP_CYCLE_STATS.                         *****
*       *****   The user may override these recommended values when     *****
*       *****   calculating "regular" CTHs.                             *****
*       *****                                                           *****
*       *****   Normalized CTHs will always use the recommended values  *****
*       *****    for binwidth and offset.                               *****
*
*
*
*       ***** Calculate and display phrenic overlay using calculated parameters.*****
*       *****   Ask user to verify calculated parameters or to select others.   *****
*
 84   CALL calcCTH(SPIKETIMES,ITAL,BINW,phrenic_hist,LAST_REF,
     +     LAST_TAR,E_pulse,phrenic,STIM_OFFSET,NUM_STIM,0,
     +     begin_E,end_E,TOTAL_SELECTED_CYCLES,IDs)
*
      call showCTHs(2,c_format,fildes,phrenic_hist,phrenic_hist,       
     +     norm_phrenic_hist,norm_phrenic_hist,cardiac_hist,  
     +     cardiac_hist,0.,0.,0.,bwtext1,bwtext2,bwtext3,'','',    
     +     E_pulse,LAST_REF,0,0,phrenic,LAST_TAR,BDT_FILE,[integer::],
     +     [integer::],'',QDT_FILENAME,'','','','', '','','','','','',
     +     '','','','','','','','', '','','','','','','','','','','','',
     +     0.,0.,0,0,[integer::],[integer::],cardiac_pls, STIM_OFFSET,
     +     NORM_OFFSET,0.,'','n')
      
      call character_height(fildes,.050)
      call character_width(fildes,.012)
      write (text1,'(f7.1)') BINW
      call strlength(text1,LEN(text1),l)
      text='binwidth = '//text1(1:l)//' msec.'
      call strlength(text,LEN(text),l)
      call text2d(fildes,950.,375.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      write (text1,'(f9.1)') STIM_OFFSET
      call strlength(text1,LEN(text1),l)
      text='offset = '//text1(1:l)//' msec.'
      call strlength(text,LEN(text),l)
      call text2d(fildes,950.,350.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call make_picture_current(fildes)
*
 85   print '(20(/),T5,''ENTER a value for the binwidth ''
     +''(fmt=f7.1 msec) ..or.. '',
     +''<cr> to accept this binwidth >> '',$)'
      BINW_tmp = 0.0
      read (*,fmt='(F7.1)',err=85),BINW_tmp
      if(BINW_tmp.ne.0.0)then
         BINW = BINW_tmp
         goto 84                !recalculate phrenic CTH with new binwidth
      else
         goto 90
      end if
*
 192  CALL calcCTH(SPIKETIMES,ITAL,BINW,phrenic_hist,LAST_REF,
     +     LAST_TAR, E_pulse,phrenic,STIM_OFFSET,NUM_STIM,0,
     +     begin_E,end_E,TOTAL_SELECTED_CYCLES,IDs)
*
      call showCTHs(2,c_format,fildes,phrenic_hist,phrenic_hist,     
     +     norm_phrenic_hist,norm_phrenic_hist,cardiac_hist,
     +     cardiac_hist,0.,0.,0.,bwtext1,bwtext2,bwtext3,'','',
     +     E_pulse,LAST_REF,0,0,phrenic,LAST_TAR,BDT_FILE,[integer::],
     +     [integer::],'',QDT_FILENAME,'','','','','','','','','','','',
     +     '','','','','','','','','','','','','','','','','','','',
     +     0.,0.,0,0,[integer::],[integer::],         
     +     cardiac_pls,STIM_OFFSET,NORM_OFFSET,mean_E,
     +     '','n')

      call character_height(fildes,.050)
      call character_width(fildes,.012)
      write (text1,'(f7.1)') BINW
      call strlength(text1,LEN(text1),l)
      text='binwidth = '//text1(1:l)//' msec.'
      call strlength(text,LEN(text),l)
      call text2d(fildes,950.,375.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      write (text1,'(f9.1)') STIM_OFFSET
      call strlength(text1,LEN(text1),l)
      text='offset = '//text1(1:l)//' msec.'
      call strlength(text,LEN(text),l)
      call text2d(fildes,950.,350.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call make_picture_current(fildes)
 90   print '(20(/),T5,''ENTER a value for the offset ''
     +''(fmt=+/- f9.1 msec) ..or.. ''
     +''<cr> to accept this offset  >> '',$)'
      STIM_OFFSET_tmp=0.0
      read (*,fmt='(F9.1)',err=90),STIM_OFFSET_tmp
      if(STIM_OFFSET_tmp.ne.0.0)then
         STIM_OFFSET = STIM_OFFSET_tmp
         goto 192
      else
         call clear_right_half(fildes)
      end if
*
      if(NUM_STIM.eq.0)NUM_STIM=TOTAL_SELECTED_CYCLES !use all cycles

*
*       ***** Calculate and display phrenic overlay using selected parameters.  *****
*       *****   Ask user to verify selected parameters or to re-enter them.     *****
*

      call character_height(fildes,.050)
      call character_width(fildes,.012)
      write (text1,'(I5)') I_pulse
      call strlength(text1,LEN(text1),l)
      text='I pulse = '//text1(1:l)
      call text2d(fildes,950.,300.,text//char(0),
     +     ANNOTATION_TEXT,0)
      write (text1,'(I5)') E_pulse
      call strlength(text1,LEN(text1),l)
      text='E pulse = '//text1(1:l)
      call strlength(text1,LEN(text1),l)
      call text2d(fildes,950.,275.,text//char(0),
     +     ANNOTATION_TEXT,0)
      write (text1,'(I5)') phrenic
      call strlength(text1,LEN(text1),l)
      text='phrenic = '//text1(1:l)
      call text2d(fildes,950.,250.,text//char(0),
     +     ANNOTATION_TEXT,0)
      if(BNDRY.ne.0)then
         write (text1,'(I5)') BNDRY
         call strlength(text1,LEN(text1),l)
         text='boundary code = '//text1(1:l)
         call text2d(fildes,950.,225.,text//char(0),
     +        ANNOTATION_TEXT,0)
      end if
      if(C_pulse.eq.'y')then
         write (text1,'(I5)') cardiac_pls
         call strlength(text1,LEN(text1),l)
         text='cardiac pulse = '//text1(1:l)
      elseif(C_pulse.eq.'n')then
         text='no cardiac pulse'
      else
         text='problem with C pulse'
      end if
      call strlength(text,LEN(text),l)
      call text2d(fildes,950.,200.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      write (text1,'(f7.1)') BINW
      call strlength(text1,LEN(text1),l)
      text='binwidth (msec)= '//text1(1:l)
      call text2d(fildes,950.,150.,text//char(0),
     +     ANNOTATION_TEXT,0)
      write (text1,'(f9.1)') STIM_OFFSET
      call strlength(text1,LEN(text1),l)
      text='offset (msec)= '//text1(1:l)
      call text2d(fildes,950.,125.,text//char(0),
     +     ANNOTATION_TEXT,0)
      write (text1,'(I5)') NUM_STIM
      call strlength(text1,LEN(text1),l)
      text='# E pulses used = '//text1(1:l)
      call text2d(fildes,950.,100.,text//char(0),
     +     ANNOTATION_TEXT,0)
c       text='Significance testing of respiratory modulation:'
c          call text2d(fildes,450.,50.,text//char(0),
c     +        ANNOTATION_TEXT,0)
c           write (text1,'(f9.1)') MAXCYCLE
c            call strlength(text1,LEN(text1),l)
c           text='max. allowed cycle = '//text1(1:l)//' msec.'
c           call text2d(fildes,950.,50.,text//char(0),
c     +        ANNOTATION_TEXT,0)
c          write (text1,'(f9.1)') user_max_value
c            call strlength(text1,LEN(text1),l)
c           text='max. control cycle = '//text1(1:l)//' msec.'
c           call text2d(fildes,950.,25.,text//char(0),
c     +        ANNOTATION_TEXT,0)
      call make_picture_current(fildes)
*
 95   print '(20(/),T5,
     +''Accept these values for respiratory CTHs (y/n)?  >> '',$)'
*
      OK=' '
      read (*,fmt='(A)',err=95), OK
      call upper_case(OK,LEN(OK))
*
      if(OK.eq.'Y')then
      else if(OK.eq.'N')then    !re-enter CTH parameters
         I_pulse=0
         E_pulse=0
         phrenic=0
         BNDRY=0
         BINW=0.0
         STIM_OFFSET=0.0
         NUM_STIM=0
         cardiac_pls=0
         goto 50
      else
         goto 95
      end if

      call clear (fildes)
      if(NUM_STIM.eq.TOTAL_SELECTED_CYCLES)NUM_STIM=0 !reset the flag
*
*
*       
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                       *
*       *       get information about cardiac cycle times..     *
*       *                                                       *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*       ***** use only those cardiac pulses which occur within *****
*       *****  selected respiratory cycles to calculate        *****
*       *****  cardiac CTHs:                                   *****
*
      BINW2 = 0
      if(C_pulse.eq.'y')then    !IF_THERE_IS_A_C_PULSE: 
         istart=1
         selected_C=0.0
         k=0
         do i = 1,TOTAL_SELECTED_CYCLES !SELECT_C_PULSES
            do j = istart,MAX_NUM_EVENTS
               if(SPIKETIMES(j,IDs(cardiac_pls)).lt.begin_E(i))cycle !get another C pulse
*
               if(SPIKETIMES(j,IDs(cardiac_pls)).ge.end_E(i))then
                  istart=j      !C pulse occurs after this resp cycle, so ...
                  exit          !cycle SELECT_C_PULSES               ! get the next selected resp cycle
               end if
*
               if((SPIKETIMES(j,IDs(cardiac_pls)).ge.begin_E(i)).and.
     +              (SPIKETIMES(j,IDs(cardiac_pls)).lt.end_E(i)))then
                  k=k+1
                  selected_C(k)=SPIKETIMES(j,IDs(cardiac_pls))
               end if
*
            end do
         end do                 !SELECT_C_PULSES
         total_sel_C_pulses=k
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *       cardiac_offset = average time between cardiac   * 
*       *             timing pulses...                          *
*       *                                                       * 
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*       
*      define cardiac_offset ...(use all cardiac cycles)
*
         itmp = ITAL(IDs(cardiac_pls)) - 1
         sumdif=0.0
         do i=1,itmp
            diff=SPIKETIMES(i+1,IDs(cardiac_pls)) -
     +           SPIKETIMES(i,IDs(cardiac_pls))
            RTAL(i) = diff
            sumdif=sumdif+diff
         end do
         SJOFSET = sumdif/itmp  !mean
         SUMSQR=0.0
         do I=1,itmp
            SUMSQR=((RTAL(i)-SJOFSET)**2)+SUMSQR
         end do
*
*       calculate mean and standard deviation of cardiac cycle length:
*
         CSTDDEV=SQRT(SUMSQR/(itmp-1))
         cardiac_offset = INT(SJOFSET) ! integer default mean
*       max cycle length used for anova test for card modulation
         ICMAXCYC = cardiac_offset + 2*(CSTDDEV) 
         call dp_real_mean_and_sd(RTAL,ITAL(IDs(cardiac_pls))-1,
     +        SJOFSET,CSTDDEV)
         write (cardiac_offset_text,'(I5)') cardiac_offset
         write (CSTDDEV_text,'(f8.1)') CSTDDEV
         write (ICMAXCYC_text,'(I5)') ICMAXCYC
*
*       ***** Display cardiac pulse overlay using default binwidth=10 ms: *****
*
         BINW2 = 10.0           !default BW for cardiac CCHs is 10 msec.
 99      CALL calc_cardCCH(SPIKETIMES,ITAL,BINW2,cardiac_hist,
     +        LAST_REF, LAST_TAR,cardiac_pls,
     +        cardiac_pls,IDs,
     +        selected_C,total_sel_C_pulses)
         call showCTHs(3,c_format,fildes,phrenic_hist,phrenic_hist,       
     +        norm_phrenic_hist,norm_phrenic_hist,cardiac_hist,  
     +        cardiac_hist,0.,0.,0.,bwtext1,bwtext2,bwtext3,'','',    
     +        cardiac_pls,LAST_REF,0,0,cardiac_pls,LAST_TAR,     
     +        BDT_FILE,[integer::],[integer::],'',QDT_FILENAME,
     +        '','','','','','','','','','','','','','','',
     +        '','','','','','','','','','','','','','','',
     +        0.,0.,0,0,[integer::],[integer::], 
     +        cardiac_pls,STIM_OFFSET,NORM_OFFSET,0.,'','n')
*
         call character_height(fildes,.050)
         call character_width(fildes,.012)
         text='mean cardiac cycle = '//cardiac_offset_text//
     +        ' +/- '//CSTDDEV_text//' msec.'
         call text2d(fildes,325.,350.,text//char(0),
     +        ANNOTATION_TEXT,0)
         text='maximum cardiac cycle length for ANOVA = '
     +        //ICMAXCYC_text//' msec.'
         call text2d(fildes,325.,300.,text//char(0),
     +        ANNOTATION_TEXT,0)
         write (bwtext,'(f4.1)') BINW2
c         bwtext3=bwtext//' ms.'        
         text='binwidth = '//bwtext//' msec.'
         call text2d(fildes,325.,250.,text//char(0),
     +        ANNOTATION_TEXT,0)
         call make_picture_current(fildes)
*
*
*       ***** display cardiac CTH parameters and ask for verification: *****
*
 111     print '(20(/),T5,''ENTER a binwidth for the cardiac CCHs''
     +'' (fmt=f7.1 msec) ..or.. '',
     +''<cr> to accept this binwidth >> '',$)'
         BINW2_tmp = 0.0
         read (*,'(f7.1)',err=111) BINW2_tmp
         if(BINW2_tmp.ne.0.0)then
            BINW2 = BINW2_tmp
            goto 99
         end if
         call clear(fildes)

      end if                    !IF_THERE_IS_A_C_PULSE

*     ************************************************
*     ***** re-display respiratory cycle values: *****
*     ************************************************

*     SKIP THIS CODE FOR NOW - IT DOESN'T WORK CORRECTLY

c     re-enable part of this code - see later goto
c      if (.true.) goto 92

c       mean_resp_cycle_all = 0.0
c        sum_of_squares = 0.0
c        sum_of_durations = 0.0
c        std_dev_all = 0.0
      longest_all = 0.0
      shortest_all = 0.0
      do i = 1, ITAL(IDs(E_pulse))-1
         TEMP_durations(i) = 
     +        SPIKETIMES(i+1,IDs(E_pulse))-SPIKETIMES(i,IDs(E_pulse))
c           sum_of_durations = sum_of_durations + TEMP_durations(i)
      end do

      call dp_real_mean_and_sd(TEMP_durations,ITAL(IDs(E_pulse))-1,
     +     mean_resp_cycle_all,std_dev_all)
      print '(''mean = '',f7.1)',mean_resp_cycle_all
      print '(''sd = '',f7.1)',std_dev_all

c        mean_resp_cycle_all = sum_of_durations/FLOAT(ITAL(IDs(E_pulse)))
c        do i=1,ITAL(IDs(E_pulse))-1
c           sum_of_squares = sum_of_squares +
c     +                  ((TEMP_durations(i)-mean_resp_cycle-all)**2)
c        end do
c       std_dev_all=SQRT(sum_of_squares/float(ITAL(IDs(E_pulse)-1)))
      longest_all = TEMP_durations(1) !load longest & shortest with initial
      shortest_all = TEMP_durations(1) ! values (for comparison)
      do i = 1, ITAL(IDs(E_pulse))-1
         longest_all=AMAX1(real(TEMP_durations(i)),longest_all) !find the longest and shortest
         shortest_all=AMIN1(real(TEMP_durations(i)),shortest_all) ! durations
      end do

c     code still disabled from here - see SKIP THIS CODE, above
      goto 92

      text='ALL'
      call text2d(fildes,500.,400.,text//char(0),
     +     ANNOTATION_TEXT,0)
      text='CONTROL CYCLES (all times in milliseconds)'
      call text2d(fildes,700.,400.,text//char(0),
     +     ANNOTATION_TEXT,0)

      text='shortest cycle length:'
      call text2d(fildes,150.,360.,text//char(0),
     +     ANNOTATION_TEXT,0)
      write (c_short_all,'(f7.1)') shortest_all
      call remove_all_blanks(c_short_all,LEN(c_short_all))
      call strlength(c_short_all,LEN(c_short_all),l)
      call text2d(fildes,500.,360.,
     +     c_short_all(1:l)//char(0),ANNOTATION_TEXT,0)
      write (c_short_con,'(f7.1)') shortest_control
      call remove_all_blanks(c_short_con,LEN(c_short_con))
      call strlength(c_short_con,LEN(c_short_con),l)
      call text2d(fildes,700.,360.,
     +     c_short_con(1:l)//char(0),ANNOTATION_TEXT,0)

      text='longest cycle length:'
      call text2d(fildes,150.,340.,text//char(0),
     +     ANNOTATION_TEXT,0)
      write (c_long_all,'(f7.1)') longest_all
      call remove_all_blanks(c_long_all,LEN(c_long_all))
      call strlength(c_long_all,LEN(c_long_all),l)
      call text2d(fildes,500.,340.,
     +     c_long_all(1:l)//char(0),ANNOTATION_TEXT,0)
      write (c_long_con,'(f7.1)') longest_control
      call remove_all_blanks(c_long_con,LEN(c_long_con))
      call strlength(c_long_con,LEN(c_long_con),l)
      call text2d(fildes,700.,340.,
     +     c_long_con(1:l)//char(0),ANNOTATION_TEXT,0)

      text='mean cycle length:'
      call text2d(fildes,150.,320.,text//char(0),
     +     ANNOTATION_TEXT,0)
      write (c_mean_all,'(f7.1)') mean_resp_cycle_all
      call remove_all_blanks(c_mean_all,LEN(c_mean_all))
      call strlength(c_mean_all,LEN(c_mean_all),l)
      call text2d(fildes,500.,320.,
     +     c_mean_all(1:l)//char(0),ANNOTATION_TEXT,0)
      write (c_mean_con,'(f7.1)') mean_resp_cycle_control
      call remove_all_blanks(c_mean_con,LEN(c_mean_con))
      call strlength(c_mean_con,LEN(c_mean_con),l)
      call text2d(fildes,700.,320.,
     +     c_mean_con(1:l)//char(0),ANNOTATION_TEXT,0)

      text='standard dev:  +/-'
      call text2d(fildes,150.,300.,text//char(0),
     +     ANNOTATION_TEXT,0)
      write (c_std_all,'(f7.1)') std_dev_all
      call remove_all_blanks(c_std_all,LEN(c_std_all))
      call strlength(c_std_all,LEN(c_std_all),l)
      call text2d(fildes,500.,300.,
     +     c_std_all(1:l)//char(0),ANNOTATION_TEXT,0)
      write (c_std_con,'(f7.1)') std_dev_control
      call remove_all_blanks(c_std_con,LEN(c_std_con))
      call strlength(c_std_con,LEN(c_std_con),l)
      call text2d(fildes,700.,300.,
     +     c_std_con(1:l)//char(0),ANNOTATION_TEXT,0)

      call make_picture_current(fildes)

 92   print '(20(/),''The following are default values:''
c MAXCYCLE (case 2 in resp_sig)
     +//,T5,''1.  10,000.0 msec'',T30,''The longest respiratory ''
     +''cycle length that will be included ''
     +/,T30 ''when calculating the significance of ''
     +''a cell''''s respiratory modulation ''
     +/,T30,''over the entire recording.''
c user_max_value (case 5 in resp_sig)
     +//,T5,''2.  10,000.0 msec'',T30,''The longest ''
     +''respiratory cycle length that will be included ''
     +/,T30,''when calculating the significance of ''
     +''a cell''''s respiratory modulation ''
     +/,T30,''over the cycles within the control period.''
     +//,T5,''3.  10,000.0 msec'',T30,''The longest respiratory ''
c MAXduration
     +''cycle length that will be included''
     +/,T30,''when calculating ACHs and CCHs.''
c MAX_INT
     +//,T5,''4.  10,000.0 msec'',T30,''The maximum interspike ''
     +''interval used when testing for ''
     +/,T30,''serotonergic-like properties of each cell.''
     +//,T5,''NOTE: the longest selected cycle for CTHs = '',f0.1,
     +'' msec.''
     +/,T5,''      the longest control cycle           = '',f0.1,
     +'' msec.''
     +/,T5,''      the longest cycle in the file       = '',f0.1,
     +'' msec.''
     +///,T15,''Accept these default values ''
     +''[y/n/i(information)]?  >> '',$)',
     +     longest,longest_control,longest_all
      
      read (*,'(A1)') OK
      call upper_case(OK,LEN(OK))
      if(OK.eq.'Y')then
         MAXCYCLE = 10000.0
         user_max_value = 10000.0
         MAX_INT = 10000.0
         MAXduration = 10000.0
      elseif((OK.eq.'N').or.(OK.eq.'I'))then
         print '(20(/),T5,''1.  Respiratory cycles can be excluded ''
     +''from significance testing of their ''
     +''respiratory pattern ''
     +/,T9,''based upon their duration.''
     +/,T9,''Spikes occuring during a cycle within the ''
     +''entire recording that is ''
     +''longer than this value ''
     +/,T9,''will not be used for ''
     +''significance testing of a cell''''s ''
     +''respiratory modulation. (case 2)''
     +/,T9,''Spikes from such a cycle will, however, be ''
     +''included in the respiratory CTH ''
     +''if you earlier chose to ''
     +/,T9,''include cycles longer ''
     +''than this value.'')'
         if(OK.eq.'N')then
 93         print '(/,T9,''ENTER this maximum cycle length in ''
     +''milliseconds (f9.1; must have decimal point)  ''
     +/,T12,''..or.. <cr> for default value of ''
     +''10,000.0 msec  >> '',$)'
            read (*,fmt='(F9.1)',err=93), MAXCYCLE
            if(MAXCYCLE.eq.0.0)MAXCYCLE=10000.0 !default value = 10 sec max. duration
         end if
         print '(/,T5,''2.  Respiratory cycles can be excluded ''
     +''from significance testing of their ''
     +''respiratory pattern ''
     +/,T9,''based upon their duration.''
     +/,T9,''Spikes occuring during a cycle within the ''
     +''control period that is ''
     +''longer than this value ''
     +/,T9,''will not be used for ''
     +''significance testing of a cell''''s ''
     +''respiratory modulation. (case 5)''
     +/,T9,''Spikes from such a cycle will, however, be ''
     +''included in the respiratory CTH ''
     +''if you earlier chose to ''
     +/,T9,''include cycles longer ''
     +''than this value.'')'
         if(OK.eq.'N')then
 94         print '(/,T9,''ENTER this maximum cycle length in ''
     +''milliseconds (f9.1; must have decimal point)  ''
     +/,T12,''..or.. <cr> for default value of ''
     +''10,000.0 msec  >> '',$)'
            read (*,fmt='(F9.1)',err=94), user_max_value
            if(user_max_value.eq.0.0)user_max_value=10000.0 !default value = 10 sec max. duration
         end if
         print '(/,T5,''3.  Spikes occuring during cycles ''
     +''longer than this value will not be included in ''
     +/,T9,''ACHs, CCHs, or shift-control CCHs.'')'
         if(OK.eq.'N')then
 112        print '(/,T9,''ENTER this maximum cycle length in ''
     +''milliseconds (f9.1; must have decimal point)  ''
     +/,T12,''..or.. <cr> for default value of ''
     +''10,000.0 msec  >> '',$)'
            read (*,fmt='(F9.1)',err=112), MAXduration
            if(MAXduration.eq.0.0)MAXduration=10000.0 !default value = 10 sec max. duration
         end if
         
         print '(/,T5,''4.  Interspike intervals greater ''
     +''than this value will not be used to calculate ''
     +''discriminant scores ''
     +/,T9,''for evaluation of serotonergic-like ''
     +''properties of each neuron. '')'
         if(OK.eq.'N')then
 1121       print '(/,T9,''ENTER the maximum interspike interval in ''
     +''milliseconds (f9.1; must have decimal point)  ''
     +/,T12,''..or.. <cr> for default value of ''
     +''10,000.0 msec  >> '',$)'
            read (*,fmt='(f8.1)',err=1121), MAX_INT
            if(MAX_INT.eq.0.0)MAX_INT=10000.0 !default value = 10 sec max. duration
         end if
         if(OK.eq.'I')then
            print '(//,T10,''Press <cr> to continue  >> '',$)'
            read (*,'(A1)')
            goto 92             !ask again
         end if
      else
         goto 92                !force an appropriate answer
      end if
      write (c_MAX_INT,'(f8.1)') MAX_INT !write integer value as text
      
*
*       ***** load array [excluded()] that will hold ID codes of cells which will *****
*       *****           NOT be included in CTHs                                  *****
*


      excluded = 0
      print '(15(/),''ENTER digital codes to be EXCLUDED''
     +'' from all calculations:''
     +/,'' (NOTE: I pulse, E pulse,''
     +'' boundary code, and phrenic have''
     +'' already been excluded.)'')'
      call doomed_cells (SPIKETIMES, size(spiketimes,1), ITAL, IDs,
     +     MAXduration, excluded,e_pulse, i_pulse, doomed)
      print '(/)'
*     +          /,''       Enter codes one at a time.'',/)'
      excluded(E_pulse)=1       !load flags for IDs of: E pulse, 
      excluded(I_pulse)=1       !                       I pulse, 
      excluded(phrenic)=1       !                       phrenic,
      if(C_pulse.eq.'y')excluded(cardiac_pls)=1 !    cardiac pulse, (if there is one)
      if(BNDRY.ne.0)excluded(BNDRY)=1 !                       and boundary code (if there is one)
      excluded(89)=1            !exclude both possible IDcodes for phrenic (89 and 90)
      excluded(90)=1
!the above ID codes are "automatically" excluded from
!  calculation of CCHs and CTHs
*
*       ***** make the user verify code status *****
*
 130  print '(T15,''INCLUDED CODES'',T75,''EXCLUDED CODES'',$)'
      call exclude_doomed (excluded, doomed)
      j = 0
      j1 = 0
      n = 1
      m = 1
      k = 0
      i = 0
      j2 = 0
 134  if(i.eq.MAX_NUM_CODES)then
         j2=0
         print '(/,$)'
         goto 135
      end if
      print '(/,$)'
      do i = n, MAX_NUM_CODES         
         if((IDs(i).ne.0).and.(excluded(i).eq.0))then
            print '(I3,''  '',$)',i ! print INcluded codes one by one
            j = j+1
         end if   
         if((j.eq.10).or.(i.eq.MAX_NUM_CODES))then
            j2 = j              !store the value of j
            j = 0
            n = i + 1
            goto 135
         end if
      end do
*
 135  if(k.eq.MAX_NUM_CODES)goto 136 !all EXcluded codes have been printed
      do k=1,5*(10-j2)+10
         print '('' '',$)'
      end do
*
      do k = m, MAX_NUM_CODES
         if((IDs(k).ne.0).and.(excluded(k).ne.0))then
            print '(I3,''  '',$)',k !print EXcluded codes one by one
            j1 = j1+1
         end if
         if((j1.eq.10).or.(k.eq.MAX_NUM_CODES))then
            m = k + 1
            j1 = 0
            goto 134
         end if
      end do
 136  if(i.lt.MAX_NUM_CODES)goto 134
*
      change_it=0
 140  print '(//,''ENTER the ID code you wish to change ..or.. ''
     +''-1 to exclude all ..or.. <cr> to continue  >> '',$)'
      read (*,fmt='(I5)',err=140) change_it
      if(change_it.eq.-1)then   !move all "included" to "excluded"
         do kk = 1, MAX_NUM_CODES
            if(IDs(kk).ne.0)excluded(kk)=1
         end do
         goto 130
      end if
      if(change_it.ne.0)then
         if(IDs(change_it).eq.0)then
            print '(''not a valid ID code'')'
            print '(30(/))'
            goto 130
         end if
         if(excluded(change_it).eq.0)then
            excluded(change_it)=1
            print '(30(/))'
            goto 130
         end if
         if(excluded(change_it).eq.1)then
            excluded(change_it)=0
            print '(30(/))'
            goto 130
         end if
      end if

      markers=excluded         
      cells=0
      print '(15(/),''Please indicate which of the excluded''
     +'' codes are markers and which are cells.'',//)'
1130  print '(T15,''MARKER CODES'',T75,''CELL CODES'',$)'
      j = 0
      j1 = 0
      n = 1
      m = 1
      k = 0
      i = 0
      j2 = 0
 1134  if(i.eq.MAX_NUM_CODES)then
         j2=0
         print '(/,$)'
         goto 1135
      end if
      print '(/,$)'
      do i = n, MAX_NUM_CODES         
         if((IDs(i).ne.0).and.(markers(i).eq.1))then
            print '(I3,''  '',$)',i ! print excluded/marker codes one by one
            j = j+1
         end if   
         if((j.eq.10).or.(i.eq.MAX_NUM_CODES))then
            j2 = j              !store the value of j
            j = 0
            n = i + 1
            goto 1135
         end if
      end do
*
 1135  if(k.eq.MAX_NUM_CODES)goto 1136 !all EXcluded/marker codes have been printed
      do k=1,5*(10-j2)+10
         print '('' '',$)'
      end do
*
      do k = m, MAX_NUM_CODES
         if((IDs(k).ne.0).and.(cells(k).ne.0))then
            print '(I3,''  '',$)',k !print cell codes one by one
            j1 = j1+1
         end if
         if((j1.eq.10).or.(k.eq.MAX_NUM_CODES))then
            m = k + 1
            j1 = 0
            goto 1134
         end if
      end do
 1136  if(i.lt.MAX_NUM_CODES)goto 1134
*
      change_it=0
 1140  print '(//,''ENTER the ID code you wish to change ..or.. ''
     +''-1 to indicate that all these codes are cells ''
     +''..or.. <cr> to continue  >> '',$)'
      read (*,fmt='(I5)',err=1140) change_it
      if(change_it.eq.-1)then   !move all "markers" to "cells"
         do kk = 1, MAX_NUM_CODES
            if(IDs(kk).ne.0)cells(kk)=1
         end do
         goto 1130
      end if
      if(change_it.ne.0)then
         if(IDs(change_it).eq.0)then
            print '(''not a valid ID code'')'
            print '(30(/))'
            goto 1130
         end if
         if(cells(change_it).eq.0)then
            cells(change_it)=1
            markers(change_it)=0
            print '(30(/))'
            goto 1130
         elseif(cells(change_it).eq.1)then
            cells(change_it)=0
            markers(change_it)=0
            print '(30(/))'
            goto 1130
         end if
      end if


*
*
*
*
*       ***** load array [excluded()] that will hold ID codes of cells  *****
*       *****           which will NOT be included in CTHs and CCHs     *****
*
*
*       ***** remove excluded cells from the ID list:  *****
*
      do i=1,MAX_NUM_CODES
c         if((excluded(i).eq.1).and.(i.ne.E_pulse).and.
c     +        (i.ne.I_pulse).and.(i.ne.phrenic).and.
c     +        (i.ne.cardiac_pls).and.(i.ne.BNDRY))then

         if(cells(i).eq.1)then
            IDs(i) = 0
            excluded(i) = 0
         end if
      end do                    !now the excluded cell (which has not been
!  used in any calculations) is no longer
!  included in the ID list.  It is as if
!  it no longer exists. However, marker ID codes (stim pulses, etc.) are still in the excluded list.
*
*
*       ***** Create a new array (included()) that will allow calculation of            *****
*       *****   the location of the histograms (= the record number)                    *****
*       *****   for cell n within the .qdt file:        [where N=total # cells]         *****
*       *****                                                                           *****
*       *****           location of resp CTH = (((included(n)*3)+1)                     *****
*       *****           location of ACH @ binwidth 1 = (((included(n)*4)+3N             *****
*       *****                                                                           *****
*       *****  To determine the location (record number) of the CCH @ binwidth 1:       *****
*       *****                                                                           *****
*       *****   record number = 3 + #CTHs + #ACHs + relative position of CCH            *****
*       *****                                           w/in all CCHs                   *****
*       *****                                                                           *****
*       *****   1. determine the position of the cell pair w/in the pair queue:         *****
*       *****      (sample queue(REF,TAR): 1,2  1,3  1,4  2,3  2,4  3,4)                *****
*       *****           N=4 in this case                                                *****
*       *****                                                                           *****
*       *****    1st occurrence of REF cell in the queue:                               *****
*       *****           =(N-1)+(N-2)+..+(N-included(REF)-1))+1                          *****
*       *****                   Q_pos = 0                                               *****
*       *****                   if(included(REF).eq.1) goto >A<                         *****
*       *****                   do i=1,included(REF)-1                                  *****
*       *****                     Q_pos = Q_pos + (N-i)                                 *****
*       *****                     end do                                                *****
*       *****                                                                           *****
*       *****    relative position of TAR w/in the group of REF cells:                  *****
*       *****           =included(TAR)-included(REF)                                    *****
*       *****                                                                           *****
*       *****       >A<         Q_pos = Q_pos + (included(TAR)-included(REF))           *****
*       *****                                                                           *****
*       *****     so Q_pos now equals the position of the pair within the queue         *****
*       *****                                                                           *****
*       *****   2. now use the queue position to determine the relative position        *****
*       *****           of the 1st CCH for a pair among all CCHs:                       *****
*       *****                                                                           *****
*       *****           rel_loc = (Q_pos *4)-3                                          *****
*       *****                                                                           *****
*       *****   SO..........                                                            *****
*       *****      record number = 3 + 3N + 4N + rel_loc = 3 + 7N + rel_loc             *****
*       *****                                                                           *****
*       *****  To determine the record number of a control CCH @ binwidth 1:            *****
*       *****                                                                           *****
*       *****   1. find Q_pos as above                                                  *****
*       *****   2. calculate rel_loc as above                                           *****
*       *****   3. record number=3+#CTHs+#ACHs+#CCHs+relative position of control CCH   *****
*       *****      record number=3+3N+4N+(2(N**2)-2N)+rel_loc                           *****
*       *****                   = 3 + 5N + 2(N**2) + rel-loc                            *****
*       *****                                                                           *****
*


      E_begin = 0.0
      E_end = 0.0
      cycle_length = 0.0
      Te = 0
      Ti = 0

      j=0
      do i = 1, ITAL(IDs(E_pulse))-1
         if(SPIKETIMES(i+1,IDs(E_pulse))-SPIKETIMES(i,IDs(E_pulse)) !want to use cycles of appropriate length
     +        .le.MAXduration)then
            do k = 1, ITAL(IDs(I_pulse)) !find the I pulse between these 2 E pulses
               if((SPIKETIMES(k,IDs(I_pulse)).gt.
     +              SPIKETIMES(i,IDs(E_pulse))).and.
     +              (SPIKETIMES(k,IDs(I_pulse)).lt.
     +              SPIKETIMES(i+1,IDs(E_pulse)))

     +              .and.(k.eq.ITAL(IDs(I_pulse)).or. !only allow one I_pulse
     +              SPIKETIMES(k+1,IDs(I_pulse)).gt.
     +              SPIKETIMES(i+1,IDs(E_pulse)))

     +              )then
                  j=j+1         !at the end, j = # acceptable cycles for shift
                  E_begin(j) = SPIKETIMES(i,IDs(E_pulse)) !fill the arrays with data for acceptable cycles
                  E_end(j) = SPIKETIMES(i+1,IDs(E_pulse))
                  I_assoc(j) = SPIKETIMES(k,IDs(I_pulse))
                  cycle_length(j) = SPIKETIMES(i+1,IDs(E_pulse))-
     +                 SPIKETIMES(i,IDs(E_pulse))
                  Te(j) = I_assoc(j) - E_begin(j) !store the durations of E phases of acceptable cycles
                  Ti(j) = E_end(j) - I_assoc(j) !store the durations of I phases of acceptable cycles
               end if
            end do
         end if
      end do
      num_acc_cycles = j
      if(num_acc_cycles.gt.MAX_NUM_ACC_CYCLES)
     +     num_acc_cycles=MAX_NUM_ACC_CYCLES

      IF(num_acc_cycles.le.iTOTAL_SHIFTS+1)then
         print '(/,T5,''***** WARNING: Only '',I2,
     +'' acceptable cycles.'',//,T5,''Continue with this''
     +''less-than-optimal number for average ''
     +''shift-predictor control calculation? ''
     +''(y/n)  >> '',$)',num_acc_cycles
 380     read (*,'(A)'),OK
         call upper_case(OK,LEN(OK))
         jTOTAL_SHIFTS = iTOTAL_SHIFTS
         if(OK.eq.'Y')then
            jTOTAL_SHIFTS = num_acc_cycles-2
         elseif(OK.eq.'N')then
            goto 500
         else
            goto 380
         end if
      else
         jTOTAL_SHIFTS = iTOTAL_SHIFTS
      end if
      do i = 1,total_num_qdts
         call strlength(qdt_files(i),LEN(qdt_files(i)),l)
         OPEN(UNIT=i+200,FILE=qdt_files(i)(1:l),
     +        FORM='FORMATTED',
     +        ACCESS='DIRECT',RECL=QDT_RECL)
         recnum(i) = 1
      end do


*       ***** last U-turn before qdt file generation: *****

      call getenv ("XANALYSIS_MAKE_SCRIPT",envtxt)
 150  print '(//,T5,''LAST U-TURN BEFORE HISTOGRAM FILE IS ''
     +''CREATED.  CONTINUE? (y/n)  >> '',$)'
      read (*,'(A)',err=150),OK
      call upper_case(OK,LEN(OK))
      if(OK.eq.'N')then
         if (envtxt.eq.'on') stop
         abort='y'
         call clear(fildes)
         goto 500  
      else if(OK.eq.'Y')then
         if (envtxt.eq.'on') stop
      else
         goto 150
      end if

      call character_height(fildes,.060)
      call character_width(fildes,.015)
      call FDATE(TODAY)
      text='Start time: '//TODAY
      call strlength(text,LEN(text),l)
      call text2d(fildes,30.,180.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)

      qdtstring = ' '
      do i = 1, total_num_qdts
         call strlength(qdtstring,LEN(qdtstring),l)
         qdtstring = qdtstring(1:l)//qdt_files(i)//','
      end do
      call remove_all_blanks(qdtstring,LEN(qdtstring))
      call strlength(qdtstring,LEN(qdtstring),l)
      qdtstring = qdtstring(1:l-1)
      if(total_num_qdts.eq.1)then
         text300 = 'Please wait . . . creating file '//qdtstring
      else
         text300 = 'Please wait . . . creating files '//qdtstring
      end if
      call text2d(fildes,30.,150.,text300,
     +     ANNOTATION_TEXT,0)

      
      j=0
      do i=1,MAX_NUM_CODES      !COUNT_THE_INCLUDED_CODES: 
         if((IDs(i).ne.0).and.(excluded(i).eq.0))then
            j=j+1
            included(i)=j
         end if
      end do                    !COUNT_THE_INCLUDED_CODES
      total_num_cells = j       !total number of cells used to generate histograms
*
*
      print *,'total_num_cells: ',total_num_cells

      total_histograms = 8*(j**2) - j + 3 !total_histograms = # histograms in each qdt file
      print '(30(/),T15,''Please wait while histograms are '',
     +''calculated ... '')'
      call make_picture_current(fildes)

      call concat(DB_FILES,'_mean_E.sav',MEAN__E,l)
      OPEN (UNIT=27,FILE=MEAN__E,FORM='UNFORMATTED')
      write (27) mean_E
      close (27)
*

*
*       ***** Go through spike trains of all cells that will be used *****
*       *****  to calculate CTHs, ACHs, CCHs, and shift-control CCHs *****
*       *****  and discard events that occur prior to the first or   *****
*       *****  after the last respiratory cycle in the recording:    *****
*
*
      first_E = SPIKETIMES(1,IDs(E_pulse))
      last_E = SPIKETIMES(ITAL(IDs(E_pulse)),IDs(E_pulse))
*
c       CELLS: do m = 1,MAX_NUM_CODES
c        if(IDs(m).eq.0)cycle CELLS                     !not a valid code
c        if(excluded(m).eq.1)cycle CELLS                        !this code has been excluded
c        count = 0
c*
c         savedITAL(m) = ITAL(IDs(m))
c        do j = 1,ITAL(IDs(m))
c         if(SPIKETIMES(j,IDs(m)).eq.0.0)cycle CELLS            !reached end of data for this cell
c         if(SPIKETIMES(j,IDs(m)).lt.first_E)cycle              !get another cell event time
c*
c         if(SPIKETIMES(j,IDs(m)).ge.last_E)then                !discard all events that occur after
c          SPIKETIMES(count+1:ITAL(IDs(m)),IDs(m))=0.0                 ! the last cycle
c          cycle CELLS
c          end if
c*
c         count=count+1
c         SPIKETIMES(count,IDs(m))=SPIKETIMES(j,IDs(m))
c         ITAL(IDs(m)) = count
c*
c         end do
c        end do CELLS
      do m = 1, MAX_NUM_CODES   !CELLS1: 
         if(IDs(m).eq.0)cycle   !CELLS1                 !not a valid code
         if(excluded(m).eq.1)cycle !CELLS1             !this code has been excluded
         icount = 0
         do j = 1,ITAL(IDs(m))
c            do j = 1,savedITAL(m)
            if((SPIKETIMES(j,IDs(m)).ge.first_E).and.
     +           (SPIKETIMES(j,IDs(m)).lt.last_E).and.
     +           (SPIKETIMES(j,IDs(m)).ne.0))then
               icount = icount+1
               SPIKETIMES(icount,IDs(m))=SPIKETIMES(j,IDs(m))
            end if
         end do
         if(icount.ne.ITAL(IDs(m)))then
c            if(icount.ne.savedITAL(m))then
            SPIKETIMES(icount+1:ITAL(IDs(m)),IDs(m))=0.0 !discard all events that occur after the last cycle
c               SPIKETIMES(icount+1:savedITAL(m),IDs(m))=0.0      !discard all events that occur after the last cycle
            ITAL(IDs(m)) = icount
            newITAL = icount
         end if
         if(ITAL(IDs(m)).eq.0)then !exclude this cell - no events left
c           if(newITAL.eq.0)then                                !exclude this cell - no events left
c               excluded(m) = 1
c               print '(/,''** No events remaining for IDcode '',I3,
c     +              '', so it''''s been excluded. **'')',m
         end if
      end do                    !CELLS1

*
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                       *
*       *       calculate phrenic overlay histogram and write   *
*       *         it to .qdt. file                              *
*       *                                                       *
*       *   *****   THIS IS RECORD #1 IN .qdt FILE      *****   *
*       *                                                       *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
      print '(/,''Calculating CTHs . . .'')'
      CALL calcCTH(SPIKETIMES,ITAL,BINW,IHIST,LAST_REF,LAST_TAR,
     +     E_pulse,phrenic,STIM_OFFSET,NUM_STIM,0,
     +     begin_E,end_E,TOTAL_SELECTED_CYCLES,IDs)
*
      NHW=INT(BINW*50.)
      ICN(1)=E_pulse
      ICN(2)=LAST_REF
      ICN(3)=phrenic
      ICN(4)=LAST_TAR
      ICN(5)=E_pulse
      ICN(6)=LAST_REF
      IF (NUM_STIM.NE.0) ICN(6)=NUM_STIM
*
*               ***** write phrenic CTH to *.qdt: ***** 
* 
      do i = 1,total_num_qdts
         write (i+200,c_format,rec=recnum(i)) 2,NHW,ICN,IHIST !0=cross,1=cont cross,2=auto,cth or pst
         recnum(i) = recnum(i) + 1
         histogram_number = histogram_number + 1 !update counter
      end do
*
*
*
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                       *
*       *       calculate NORMALIZED phrenic overlay histogram  *
*       *         and write it to .qdt. file                    *
*       *                                                       *
*       *   *****   THIS IS RECORD #2 IN .qdt FILE      *****   *
*       *                                                       *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
*
      call SCALCN (begin_E,end_E,middle_I,ebins,ibins,IDs,
     +     ITAL,NORM_BW,NORM_OFFSET,SPIKETIMES,TOTAL_SELECTED_CYCLES,
     +     phrenic,IHIST,spike_count)
*
      NHW=INT(NORM_BW*50.)
      ICN(1)=E_pulse
      ICN(2)=TOTAL_SELECTED_CYCLES           !# of E-pulses
      ICN(3)=phrenic               
      ICN(4)=spike_count           !# of normalized unit events
      ICN(5)=ebins
      ICN(6)=ibins
      IF (NUM_STIM.NE.0) then
         print *,'NUM_STIM is ',NUM_STIM
         stop
      end if
*
*               ***** write normalized phrenic CTH to *.qdt: *****      
* 
      do i = 1,total_num_qdts
         write (i+200,c_format,rec=recnum(i)) 3,NHW,ICN,IHIST !0=cross,1=cont cross,2=auto,cth or pst,3=phase norm CTH
         recnum(i) = recnum(i) + 1
         histogram_number = histogram_number + 1 !update counter
      end do
*
*
*
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                       *
*       *       calculate cardiac overlay histogram and write   *
*       *         it to .qdt file:                              *
*       *                                                       *
*       *   *****   THIS IS RECORD #3 IN .qdt FILE    *****     *
*       *                                                       *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
      if(C_pulse.eq.'y')then    !IF_THERE_IS_A_C_PULSE_2: 
         CALL calc_cardCCH(SPIKETIMES,ITAL,BINW2,IHIST,
     +        LAST_REF, LAST_TAR,cardiac_pls,
     +        cardiac_pls,IDs,
     +        selected_C,total_sel_C_pulses)
*
         NHW=INT(BINW2*50.)
         ICN(1)=cardiac_pls
         ICN(2)=LAST_REF
         ICN(3)=cardiac_pls
         ICN(4)=LAST_TAR
         ICN(5)=cardiac_pls
         ICN(6)=LAST_REF
      else
         NHW=0
         do i = 1, 6
            ICN(i) = 0
         end do
         do i = 1, NUM_BINS
            IHIST(i) = 0
         end do
      end if                    !IF_THERE_IS_A_C_PULSE_2
*
*               ***** write cardiac CCH to *.qdt: *****
*
      do i = 1,total_num_qdts
         write (i+200,c_format,rec=recnum(i)) 0,NHW,ICN,IHIST !0=cross,1=cont cross,2=auto,cth or pst
         recnum(i) = recnum(i) + 1
         histogram_number = histogram_number + 1 !update counter
      end do
*
*
*
*
*       **************************************************************************
*       **************************************************************************
*       **                                                                      **
*       **      BEGIN HUGE DO LOOP #1 that will CALCULATE and WRITE             **
*       **      the regular RESPIRATORY CYCLE-TRIGGERED HISTOGRAM, the          **
*       **      normalized RESPIRATORY CTH, and the CARDIAC                     **
*       **      CYCLE-TRIGGERED HISTOGRAM for each cell.                        **
*       **      Those cells which the user has designated as "excluded"         **
*       **      will not be used in these calculations, nor will they be used   **
*       **      later in computation of cross-correlograms.                     **
*       **                                                                      **
*       **************************************************************************
*       **************************************************************************
*
*
      do j=1,MAX_NUM_CODES      !CELL_TYPING:
         if(IDs(j).eq.0)cycle   !CELL_TYPING !this ID code not in input file
         if(excluded(j).eq.1)cycle !CELL_TYPING !j=excluded code, so get another
*
*       
*               * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*               *                                                               *
*               *       CALCULATE RESPIRATORY CTHs for                          *
*               *       each spike code (herein defined as j and TAR) and       *
*               *       write CTH to .qdt file.                                 *
*               *                                                               *
*               * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
         TAR=j                  !now have a valid, non-excluded cell ID code; 
!  define it as TAR
         
*
*
*               ***** calculate CTH for cell *****
*
*
         CALL calcCTH(SPIKETIMES,ITAL,BINW,IHIST,LAST_REF,LAST_TAR,
     +        E_pulse,TAR,STIM_OFFSET,NUM_STIM,0,
     +        begin_E,end_E,TOTAL_SELECTED_CYCLES,IDs)
*

         NHW=INT(BINW*50.)
         ICN(1)=E_pulse
         ICN(2)=LAST_REF
         ICN(3)=TAR
         ICN(4)=LAST_TAR
         ICN(5)=E_pulse
         ICN(6)=LAST_REF
         IF (NUM_STIM.NE.0) ICN(6)=NUM_STIM
*
*               ***** write cell's respiratory CTH to ".qdt" file: *****
*
         do i = 1,total_num_qdts
            write (i+200,c_format,rec=recnum(i)) 2,NHW,ICN,IHIST !0=cross,1=cont cross,2=auto,cth or pst
            recnum(i) = recnum(i) + 1
            histogram_number = histogram_number + 1 !update the counter
         end do
*
*
*
*
*       *****************************************************
         call respsig5(SPIKETIMES,TAR,E_pulse,
     +        mean_control_cycle,MAXCYCLE,
     +        control_max_value,user_max_value,
     +        mean_resp_cycle,begin_E,end_E,
     +        BNDRY_START_TIME,BNDRY_END_TIME,
     +        ETA2_1,zmodsig_1,zmodsig2_1,
     +        ETA2_2,zmodsig_2,zmodsig2_2,
     +        ETA2_3,zmodsig_3,zmodsig2_3,
     +        ETA2_4,zmodsig_4,zmodsig2_4,
     +        ETA2_5,zmodsig_5,zmodsig2_5,
     +        ETA2_6,zmodsig_6,zmodsig2_6,
     +        IDs,mode)
*
*
*       set max interspike interval allowed = mean resp cycle
*       length plus 2 sd.
         ZMAX = IMAXCYC         ! max interspike interval allowed ms
         call intvar(SPIKETIMES,ITAL(IDs(TAR)),TAR,ZMAX,coefvr,IDs)
         coefnum(IDs(TAR))=coefvr !fp number for acumsum subrtn
         write  (coef(IDs(TAR)),('(F5.2)')) coefvr
*       ****************************************************
*
*
*
*               ***** calculate normalized CTH for cell *****
*
*
         call SCALCN (begin_E,end_E,middle_I,ebins,ibins,IDs,
     +        ITAL,NORM_BW,NORM_OFFSET,SPIKETIMES,TOTAL_SELECTED_CYCLES,
     +        TAR,IHIST,spike_count)
*     
         NHW=INT(NORM_BW*50.)
         ICN(1)=E_pulse
         ICN(2)=TOTAL_SELECTED_CYCLES !# of E-pulses
         ICN(3)=phrenic               
         ICN(4)=spike_count     !# of normalized unit events
         ICN(5)=ebins
         ICN(6)=ibins
         IF (NUM_STIM.NE.0) then
            print *,'NUM_STIM is ',NUM_STIM
            stop
         end if
*
*               ***** write cell's norm respiratory CTH to ".qdt" file: *****
*
         do i = 1,total_num_qdts
            write (i+200,c_format,rec=recnum(i)) 3,NHW,ICN,IHIST !0=cross,1=cont cross,2=auto,cth or pst,3=phase norm CTH
            recnum(i) = recnum(i) + 1
            histogram_number = histogram_number + 1 !update the counter
         end do
*
*
*
*               * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*               *                                                       *
*               *   Calculate cardiac CCHs and write to *.qdt :         *
*               *                                                       *               
*               * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
         if(C_pulse.eq.'y')then !IF_THERE_IS_A_C_PULSE_3: 
            CALL calc_cardCCH(SPIKETIMES,ITAL,BINW2,IHIST,
     +           LAST_REF,LAST_TAR,cardiac_pls,
     +           TAR,IDs,
     +           selected_C,total_sel_C_pulses)
*
            NHW=INT(BINW2*50.)
            ICN(1)=cardiac_pls
            ICN(2)=LAST_REF
            ICN(3)=TAR
            ICN(4)=LAST_TAR
            ICN(5)=cardiac_pls
            ICN(6)=LAST_REF
         else
            NHW=0
            do i = 1, 6
               ICN(i) = 0
            end do
            do i = 1, NUM_BINS
               IHIST(i) = 0
            end do
         end if                 !IF_THERE_IS_A_C_PULSE_3
*
*               ***** write cell's cardiac CCH to ".qdt" file *****
*
         do i = 1,total_num_qdts
            write (i+200,c_format,rec=recnum(i)) 0,NHW,ICN,IHIST !0=cross,1=cont cross,2=auto,cth or pst
            recnum(i) = recnum(i) + 1
            histogram_number = histogram_number + 1 !update the counter
            call percent_completed(previous_percent,fildes, !recalculate % of histograms completed
     +           histogram_number,total_histograms,total_num_qdts)
         end do
*
*
*       *** GET REMAINDER OF CARDIAC INFORMATION FOR CURRENT CELL ***
*
*
*       ********************************************************
*
         if(C_pulse.eq.'y' .AND. ITAL(IDs(TAR)) .gt. 0)
     +        call cardsig_delta2(SPIKETIMES,TAR,cardiac_pls,ICMAXCYC,
     +        card_type,IDs,tedfactor,DELTA2)
         
         
*     ******* Calculate the discriminant scores to evaluate the *****
*     *******  "serotonergic-like-ness" of each cell: *****

c           print *,TAR,IDs(TAR),ITAL(IDs(TAR))
         call serotonergic_discr_score(SPIKETIMES,IDs,ITAL,TAR,
     +        meanISI(IDs(TAR)),
     +        sdISI(IDs(TAR)),fiveHT(IDs(TAR)),
     +        MAX_INT,num_rej_ISI(IDs(TAR)))


c           print '(I3,'': 5-HT = '',A10)',TAR,fiveHT(IDs(TAR))
c           print '(I3,'': meanISI = '',A10)',TAR,meanISI(IDs(TAR))
c           print '(I3,'': sdISI = '',A10)',TAR,sdISI(IDs(TAR))

         call calc_rISI(SPIKETIMES,IDs,ITAL,
     +        mean_rISI(IDs(TAR)),sd_rISI(IDs(TAR)),TAR,
     +        MAX_INT,num_rej_rISI(IDs(TAR)))
c           print '(I3,'': mean_rISI = '',A10)',TAR,mean_rISI(IDs(TAR))
c           print '(I3,'': sd_rISI = '',A10)',TAR,sd_rISI(IDs(TAR))
c           read (*,'(A)')
*
*
*
*

*
c        print *,'end do'
      end do                    !CELL_TYPING      
*
*
*       **************************************************************************
*       **************************************************************************
*       **                                                                      **
*       **  BEGIN DO LOOP #2 that will calculate ACHs at 4 binwidths for each   **
*       **      cell and write them to the .qdt files.                          **
*       **                                                                      **
*       **************************************************************************
*       **************************************************************************
*
*
*
*               * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*               *                                                               *
*               *       BUT FIRST, remove events from each spike train that     *
*               *               do not occur within an acceptable resp cycle    *
*               *           (unacceptable cycles are longer than MAXduration    *
*               *                and do not have an I pulse)                    *
*               *        -- use these modified spike trains to compute ACHs,    *
*               *               CCHs, and shift-control CCHs                    *
*               *                                                               *
*               * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
      print '(/,''Removing unacceptable cycles from spike trains'')'
      do i = 1, MAX_NUM_CODES
         if(IDs(i).eq.0) cycle  !be sure that this IDcode exists and is included
         if(excluded(i).eq.1)cycle
         ito = 0
         do ifrom = 1, ITAL(IDs(i)) !J_LOOP:                        !look at every event in this spike train
            do k = 1, num_acc_cycles !check to see if this event occurs within an acceptable cycle
               if((SPIKETIMES(ifrom,IDs(i)).ge.E_begin(k)).and.
     +              (SPIKETIMES(ifrom,IDs(i)).lt.E_end(k)))then
                  ito = ito+1
                  SPIKETIMES(ito,IDs(i)) = SPIKETIMES(ifrom,IDs(i)) 
                  sp_per_cycle(IDs(i),k)=sp_per_cycle(IDs(i),k)+1 !count the # of spikes in each acceptable cycle
c                    k1 = k                                      !don't want to go back to the beginning each time
                  exit          !cycle J_LOOP                                !this spike is OK, so look at the next spike
               end if
            end do
         end do                 !J_LOOP
         do jj = ito+1, MAX_NUM_EVENTS
            SPIKETIMES(jj,IDs(i)) = 0.0 !fill the rest of the array section with zeros
         end do
         ITAL(IDs(i)) = ito
      end do

      do i = 1, MAX_NUM_CODES
         if(IDs(i).eq.0)cycle
      end do

*       ***** OK  -- NOW we can compute the auto-correlograms: *****

      print '(/,''Calculating ACHs . . .'')'
      do k = 1,MAX_NUM_CODES
         if(IDs(k).eq.0) cycle  !invalid code
         if(excluded(k).eq.1) cycle !exclude this cell
*
         cell = k               !have a "good" cell
         
         do i = 1,total_num_qdts !for each new QDT file:
            do j= 1,4           !for each binwidth:
               BINWa = BWs(i,j)
               call calcACH(SPIKETIMES,ITAL,BINWa,IHIST,LAST_REF,
     +              LAST_TAR,cell,cell,0.0,0,0,IDs)
*
               NHW=INT(BINWa*50.)
               ICN(1)=cell
               ICN(2)=LAST_REF
               ICN(3)=cell
               ICN(4)=LAST_TAR
               ICN(5)=cell
               ICN(6)=LAST_REF

               write (i+200,c_format,rec=recnum(i)) 2,NHW,ICN,IHIST
               histogram_number=histogram_number+1 !update the counter
               recnum(i) = recnum(i) + 1
               call percent_completed(previous_percent,fildes, !recalculate % of histograms completed
     +              histogram_number,total_histograms,total_num_qdts)
            end do
         end do
      end do
*
*       **************************************************************************
*       **************************************************************************
*       **                                                                      **
*       **      BEGIN HUGE DO LOOP #3 that will procure a REFERENCE cell,       **
*       **      a TARGET cell, calculate CROSS-CORRELATION HISTOGRAMS at        **
*       **      4 binwidths for them, and write the CCHs to the *.qdt file.     **
*       **                                                                      **
*       **      Another pair of reference and target cells will then            **
*       **      be selected automatically, and CCHs will be calculated for      **
*       **      them.  This procedure continues until all possible pairings     **
*       **      have been made.                                                 **
*       **                                                                      **
*       **      Those cells which the user has designated as "excluded"         **
*       **      will not be used in these calculations.                         **
*       **                                                                      **
*       **      Remember -- spikes occurring in non-acceptable cycles           **
*       **      (i.e., those longer than MAXduration) have been removed         **
*       **      from the spike trains.                                          **
*       **************************************************************************
*       **************************************************************************
*
*
      call sc_write_edt (SPIKETIMES, size(spiketimes, 1), IDs, ITAL,
     +     excluded, total_num_qdts,BWs, DB_FILES)
      print '(/,''Calculating CCHs . . .'')'
      do j=1,MAX_NUM_CODES      !GET_REF:  
         if(IDs(j).eq.0)cycle   !GET_REF !this ID code not in input file
         if(excluded(j).eq.1)cycle !GET_REF !j = excluded code, so get another ref code
         REF=j                  !"good" code --> define as REFerence cell
*
         do k=j+1,MAX_NUM_CODES !GET_TAR: 
            if(IDs(k).eq.0)cycle !GET_TAR        !this ID code not in input file 
            if(excluded(k).eq.1)cycle !GET_TAR    !k = excluded code, so get another tar code
            TAR= k              !"good" code --> define as TARget cell
            if(j.eq.k) cycle    !GET_TAR            !ref = tar --> not a valid pair
*
*
*               ***** REF and TAR codes now defined for 1 pair *****
*
*               
            do i = 1, total_num_qdts
               do izz=1,4       !COMPUTE_CCHs:
                  BINWx = BWs(i,izz)
*
*
*               ***** calculate CCH and write it to .qdt file: *****
*
*
                  calc_shift = 'n'
                  CALL calcCCH_newshift(calc_shift,SPIKETIMES,
     +                 SPIKESHIFT,REF,TAR,ITAL,BINWx,IHIST,0.0,
     +                 LAST_REF,LAST_TAR,IDs,count,i)

                  NHW=INT(BINWx*50.)
                  ICN(1)=REF           
                  ICN(2)=LAST_REF
                  ICN(3)=TAR
                  ICN(4)=LAST_TAR
                  ICN(5)=E_pulse
                  ICN(6)=ITAL(IDs(E_pulse))

                  write (i+200,c_format,rec=recnum(i)) 0,NHW,ICN,IHIST !0=cross,1=cont cross,2=auto,cth or pst
                  recnum(i) = recnum(i) + 1
*
                  histogram_number = histogram_number + 1 !update counter
                  call percent_completed(previous_percent,fildes, !recalculate % of histograms completed
     +                 histogram_number,total_histograms,
     +                 total_num_qdts)
               end do           !COMPUTE_CCHs
            end do
*
*
         end do                 !GET_TAR          !get next TARget cell
         prevREFcode = REF
      end do                    !GET_REF            !get next REFerence cell



*       **************************************************************************
*       **************************************************************************
*       **                                                                      **
*       **      BEGIN  DO LOOP #4 that will procure a REFERENCE cell,           **
*       **      a TARGET cell, calculate CONTROL CROSS-CORRELATION HISTOGRAMS at**
*       **      4 binwidths for them, and write the control CCHs to the *.qdt   **
*       **      file.                                                           **
*       **                                                                      **
*       **      Default calculation of the control CCHs is done using Russ's    **
*       **      method of producing surrogate spike trains which preserve       **
*       **      the cell's firing rate and respiratory pattern.  The surrogate  **
*       **      spike trains are then used to calculate control CCHs that       **
*       **      will reflect correlations between the two cells that are due    **
*       **      to their firing rates but will not include short-time scale     **
*       **      correlations.                                                   **
*       **                                                                      **
*       **      The earlier method utilizing Lauren's normalized shift of       **
*       **      spikes into the next respiratory cycle (and whose philosophy is **
*       **      painstakingly and probably inappropriately explained below)     **
*       **      can be used by starting this program with xanalysis_old.        **
*       **                                                                      **
*       **                                                                      **
*       **      Explanation of earlier control CCH method:                      **
*       **      Each REFerence spike train                                     **
*       **      undergoes a shift process, whereby each event in the train is   **
*       **      proportionally shifted into appropriate phase of the next       **
*       **      respiratory cycle. This is repeated 19 more times for a total   **
*       **      of 20 shifts.  These 20                                         **      
*       **      interim shift-predictor control CCHs are then averaged together **
*       **      and saved as the official CONTROL CCH for that pair at that     **
*       **      binwidth.  SDs of each CONTROL CCH are stored in QDT.SAV.       **
*       **      Events occurring within the last respiratory cycle              **
*       **      are shifted to the now-empty first cycle.  [Events occurring    **
*       **      prior to the first or after the last cycle have already been    **
*       **      discarded from SPIKETIMES().]                                   **
*       **                                                                      **
*       **      Here's my thinking about this:                                  **
*       **                                                                      **
*       **          [relloc==relative location; absloc==absolute location]      **
*       **          [A represents the spike located in its original cycle;      **
*       **              B is the new cycle]                                     **
*       **                                                                      **
*       **      since rellocA = abslocA - E_begin(A)                            **
*       **      and   rellocA/cycle_length(A)  =  rellocB/cycle_length(B)       **
*       **      and   abslocB = rellocB + E_begin(B)                            **
*       **      then                                                            **
*       **                                                                      **
*       **      rellocB = (rellocA * cycle_length(B)) / cycle_length(A)         **
*       **                                                                      **
*       **      abslocB = ((rellocA * cycle_length(B))/cycle_length(A)) + E_begin(B)            **
*       **                                                                      **
*       **      abslocB = (((abslocA-E_begin(A))/cycle_length(A))*cycle_length(B)) + E_begin(B) **
*       **                                                                      **
*       **      abslocB = (((abslocA-E_begin(A))*cycle_length(B))/cycle_length(A)) + E_begin(B) **
*       **                                                                      **
*       **      Remember -- spikes occurring in non-acceptable cycles           **
*       **      (i.e., those longer than MAXduration) have been removed         **
*       **      from the spike trains.                                          **
*       **                                                                      **
*       **      Another pair of reference and target cells will then            **
*       **      be selected automatically, and control CCHs will be calculated  **
*       **      for them.  This procedure continues until all possible pairings **
*       **      have been made.                                                 **
*       **                                                                      **
*       **      Those cells which the user has designated as "excluded"         **
*       **      will not be used in these calculations.                         **
*       **
*       **
*       **      Code logic:                                                     **
*       **              1.  get the cell pair - keep track of the pair number   **
*       **          ^^  2.  shift the REF spike train by x number of cycles     **
*       **                      -- round the values to the nearest:             **
*       **                              -- .5 (for BDT files) or                **   
*       **                              -- .1 (for EDT files)                   **
*       **                      -- save it in SPIKESHIFT()                      **
*       **              3.  calc the SP for each BW in each QDT file            **
*       **                      -- define SP(MAX_NUM_QDTS,4)                    **
*       **              4.  add each SP to corresponding CONTROL                **
*       **                      -- CONTROL(MAX_NUM_QDTS,4)=CONTROL() + SP()     **
*       **              5.  repeat from ^^ until have calculated the SP         **
*       **                      a total of jTOTAL_SHIFTS times                   **
*       **              6.  calc the average SP at each BW for each QDT         **
*       **                      -- CONTROL(qdt#,bw)=IFIX(float(CONTROL(qdt#,bw))/float(jTOTAL_SHIFTS)))
*       **              7.  calc the sd for each CONTROL histogram              **
*       **                   -- define sd_CONTROL(MAX_NUM_QDTS,MAX_NUM_PAIRS,4) **
*       **                      (write to qdt.sav as sd_CONTROL(total_num_qdts,pair#,4)
*       **              8.  get the next cell pair and start over               **
*       **                                                                      **
*       **************************************************************************
*       **************************************************************************
*

*
      call strlength(BDT_FILE,LEN(BDT_FILE),l) !need to know what type of input file so can
      file_ext = BDT_FILE(l-2:l) ! handle shifted values correctly
      if(file_ext.eq.'bdt')factor=2.0D0
      if(file_ext.eq.'edt')factor=10.0D0
      calc_shift = 'y'

      call getenv ("XANALYSIS_CONTROL",cntrl)
      cycle_shift = .FALSE.
      cth_cch = .TRUE.
      surrogate = .FALSE.

      IF(cntrl.EQ.'cycle_shift')then
         cycle_shift = .TRUE.
         cth_cch = .FALSE.
      ELSE IF (cntrl.EQ.'surrogate')then
         surrogate = .TRUE.
         cth_cch = .FALSE.
      END IF

      IF (cth_cch) THEN
         IF(version.ne.'8')STOP ', BUG: cth_cch but version != 8'
         jTOTAL_SHIFTS = 1
      END IF

      IF(surrogate) THEN
         call tar_surrogate(SPIKETIMES,size(spiketimes,1),ITAL,IDs
     +        ,excluded,STP,SURR_ITAL, E_begin, E_end, num_acc_cycles)
      END IF

      do REF=1,MAX_NUM_CODES    !GET_REF_1:  
         if(IDs(REF).eq.0)cycle !GET_REF_1                             !this ID code not in input file
         if(excluded(REF).eq.1)cycle !GET_REF_1                        !j = excluded code, so get another ref code

*       *****  SHIFT THE REF TRAIN (but not if it's the last included cell): *****

         last_one = 'y'
         do i = REF+1, MAX_NUM_CODES !find out if REF is the last included cell
            if(included(i).ne.0)last_one = 'n'
         end do
         IF(cycle_shift)THEN
            if(last_one.eq.'n')then !SHIFT_LOOP:
               print "(/,'Shifting spike train (IDcode ',I3,')')",REF
               do ishift = 1,jTOTAL_SHIFTS+1
                  count = 0
                  k1 = num_acc_cycles-ishift+1
                  if(ITAL(IDs(REF)).gt.100)then
                     do j1 = ITAL(IDs(REF)),1,-ITAL(IDs(REF))/10 !don't waste time starting at the beginning of the data
                        if(SPIKETIMES(j1,IDs(REF)).lt.E_begin(k1))exit !jump backwards thru the data until find a point before the
                     end do     ! beginning of the first cycle in question
                  else
                     j1=1
                  end if
                  do j = j1, ITAL(IDs(REF)) !first, shift spike times from the back of SPIKETIMES
                     REFtime = SPIKETIMES(j,IDs(REF)) !  to the front of SPIKESHIFT
                     do k = k1, num_acc_cycles
                        k2 = MOD(k,num_acc_cycles-ishift) !# of cycle that will receive shifted times
                        if((REFtime.ge.E_begin(k))
     +                       .and.(REFtime.lt.E_end(k)))then
                           count = count + 1
                           if(REFtime.lt.I_assoc(k))then !spike occurs during E phase
                              SPIKESHIFT(count,ishift)=
     +                             (((REFtime-E_begin(k))*
     +                             Te(k2))/Te(k)) + E_begin(k2)
                           else !spike occurs during I phase
                              SPIKESHIFT(count,ishift) =
     +                             (((REFtime-I_assoc(k))*
     +                             Ti(k2))/Ti(k)) + I_assoc(k2)
                           endif
                           k1 = k
                           exit !get the next REF spike time
                        end if
                     end do
                  end do

                  k1 = 1
                  do j = 1, ITAL(IDs(REF)) !now can shift the rest of the spike times
                     REFtime = SPIKETIMES(j,IDs(REF))
                     do k = k1, num_acc_cycles - ishift
                        if((REFtime.ge.E_begin(k))
     +                       .and.(REFtime.lt.E_end(k))) then
                           count = count + 1
                           k2 = k + ishift

                           if(REFtime.lt.I_assoc(k))then !spike occurs during E phase
                              SPIKESHIFT(count,ishift)= (((REFtime
     +                             -E_begin(k))*Te(k2))/Te(k)) +
     +                             E_begin(k2)
                           else !spike occurs during I phase
                              SPIKESHIFT(count,ishift) = (((REFtime
     +                             -I_assoc(k))*Ti(k2))/Ti(k)) +
     +                             I_assoc(k2)
                           endif
                           k1 = k
                           exit !get the next REF spike time
                        end if
                     end do
                  end do

*       ***** Go through the shifted spike train and force all data points to:
*       *****  for a bdt file: the nearest half of a millisecond (each spike time will end in .0 or .5) 
*       *****  for an edt file: the nearest tenth of a millisecond (each spike time will end in .0, .1, .2, ... or .9) 

                  do j = 1, count !look at each shifted spike time, one at a time
                     shifted = SPIKESHIFT(j,ishift)
                     shift1 = SPIKESHIFT(j,ishift)*factor !multiply by the appropriate factor according to type of input data
*                                                               ! to express the spike time in units of:
*                                                               !  for bdt files: half milliseconds ( =2.0D0)
*                                                               !  for edt files: tenth milliseconds ( =10.0D0)
                     
*               ***** Implement Banker's Rounding -- round a value ending in .5 to the nearest even number.
*               *****  For example, 5.5 rounds to 6.0, but 12.5 rounds to 12.0.  Using this method ensures
*               *****  that we do not add a bias to a large set of numbers by rounding up more often than
*               *****  than we round down.

                     diff = shift1-INT(shift1) !look at the fractional part of the shifted spike time
                     if(diff.lt.0.5)shift2=INT(shift1) !< .5, so round down
                     if(diff.gt.0.5)shift2=INT(shift1)+1 !> .5, so round up
                     if(diff.eq.0.5)then != .5 --> now do the Banker's Rounding:
                        if(MOD(INT(shifted),2).ne.0)then !       whole number is odd, so round up
                           shift2=INT(shift1)+1
                        else    !       whole number is even, so round down
                           shift2=INT(shift1)
                        end if
                     end if
                     SPIKESHIFT(j,ishift)=shift2/factor !divide by appropriate factor to again express the spike time
*                                                               ! in terms of milliseconds
                  end do

*       ***** Here's what just happened:
*       *****
*       *****   1. Take, for example, a shifted value from a BDT file:  SPIKESHIFT = shifted = 25.6745
*       *****   2. Multiply it by the BDT factor:  25.6745 * 2 = 51.3490 = shift1
*       *****   3. Figure out what the fractional portion is:  shift1-INT(shift1) = 51.3490-51 = .3490 = diff
*       *****   4. Since diff < .5, shift2 = INT(shift1) = 51.0
*       *****   5. Divide by the BDT factor to obtain the corrected shifted value: 51.0 / 2 = 25.5
*       *****
*       ***** Now, in the case of diff = .5, we have to use the banker's rounding formulae:
*       *****
*       *****   1. Take a shifted value from a BDT file:  SPIKESHIFT = shifted = 25.25
*       *****   2. Multiply it by the BDT factor:  25.25 * 2 = 50.50 = shift1
*       *****   3. Figure out what the fractional portion is:  shift1-INT(shift1) = 50.50-50 = .50 = diff
*       *****   4. Since diff = .5 AND the whole portion of 25.25 is odd, shift2 = INT(shift1)+1 = INT(50.50)+1 = 51.0
*       *****   5. Divide by the BDT factor to obtain the corrected shifted value: 51.0 / 2 = 25.5

*       *****   1. Take a shifted value from a BDT file:  SPIKESHIFT = shifted = 24.25
*       *****   2. Multiply it by the BDT factor:  24.25 * 2 = 48.50 = shift1
*       *****   3. Figure out what the fractional portion is:  shift1-INT(shift1) = 48.50-48 = .50 = diff
*       *****   4. Since diff = .5 AND the whole portion of 24.25 is even, shift2 = INT(shift1)+1 = INT(48.50) = 48.0
*       *****   5. Divide by the BDT factor to obtain the corrected shifted value: 48.0 / 2 = 24.0

*       *****   1. Take a shifted value from a BDT file:  SPIKESHIFT = shifted = 25.75
*       *****   2. Multiply it by the BDT factor:  25.75 * 2 = 51.50 = shift1
*       *****   3. Figure out what the fractional portion is:  shift1-INT(shift1) = 51.50-51 = .50 = diff
*       *****   4. Since diff = .5 AND the whole portion of 25.75 is odd, shift2 = INT(shift1)+1 = INT(51.50)+1 = 52.0
*       *****   5. Divide by the BDT factor to obtain the corrected shifted value: 52.0 / 2 = 26.0

*       *****   1. Take a shifted value from a BDT file:  SPIKESHIFT = shifted = 24.75
*       *****   2. Multiply it by the BDT factor:  24.75 * 2 = 49.50 = shift1
*       *****   3. Figure out what the fractional portion is:  shift1-INT(shift1) = 49.50-49 = .50 = diff
*       *****   4. Since diff = .5 AND the whole portion of 24.75 is even, shift2 = INT(shift1)+1 = INT(49.50) = 49.0
*       *****   5. Divide by the BDT factor to obtain the corrected shifted value: 49.0 / 2 = 24.5


*       ***** OK - now have totally shifted, time-corrected REF trains. *****

               end do 
            end if              !SHIFT_LOOP
         ELSE IF (surrogate)then
            call ref_surrogate(SPIKETIMES(:,IDs(REF)),ITAL(IDs(REF)),
     +           IDs(REF),SPIKESHIFT,size(spikeshift,1),spikecount)
         END IF

         do TAR = REF+1,MAX_NUM_CODES !GET_TAR_1:         !now match the REF train up with its TAR trains
            if(IDs(TAR).eq.0)cycle !GET_TAR_1              !this ID code not in input file 
            if(excluded(TAR).eq.1)cycle !GET_TAR_1          !k = excluded code, so get another tar code
            if(REF.eq.TAR) cycle !GET_TAR_1               !ref = tar --> not a valid pair
            pair_num=pair_num+1 !keep a tally of how many good pairs encountered thus far

            Q_pos = 0           !calculate the position of this pair in the pairs queue
            if(included(REF).eq.1)goto 390
            do i2 = 1,included(REF)-1                    
               Q_pos = Q_pos + (total_num_cells-i2)              
            end do
 390        Q_pos = Q_pos + (included(TAR) - included(REF))
            rel_loc = (Q_pos*4)-3

            do k = 1, NUM_BINS
               do j = 1, 4
                  do i = 1, MAX_NUM_QDTS
                     CONTROL (i,j,k)= 0
                  end do
               end do
            end do

*       *****   Proceed to calculate the SPs:  *****
            if(cth_cch)then
               print '(T5,''calculating CTH-based control CCHs for '',
     +I3, '' vs '',I3,'' '')',REF,TAR
               call calc_cth_cch (SPIKETIMES,size(spiketimes,1), ITAL,
     +              IDs(REF),IDs(TAR),IDs(E_pulse),factor)
            else
               print '(T5,''calculating shift-control CCHs for '',I3,
     +'' vs '',I3,'' '')',REF,TAR
            end if

            do i = 1, jTOTAL_SHIFTS+1
               do m = 1, total_num_qdts
                  do izz=1,4
                     REC_num_control_1 = 3 + (5*total_num_cells)+ !record # of "shift by 1" control CCH @ bw = izz
     +                    (2*(total_num_cells**2)) + rel_loc + izz-1
                     REC_num_control_2 = 3 + (3*total_num_cells)+ !record # of "shift by 2" control CCH @ bw = izz
     +                    (4*(total_num_cells**2)) + rel_loc + izz-1
                     REC_num_control_avg = 3 + total_num_cells + !record # of "averaged shift" control CCH @ bw = izz
     +                    (6*(total_num_cells**2)) + rel_loc + izz-1
                     BINWs = BWs(m,izz)
                     do i2 = 1, NUM_BINS
                        SP (i2)= 0
                     end do
                     IF(cycle_shift)THEN
                        CALL calcCCH_newshift(calc_shift,SPIKETIMES,
     +                       SPIKESHIFT,REF,TAR,ITAL,BINWs,SP,
     +                       0.0,LAST_REF,LAST_TAR,IDs,count,i)
                     ELSE IF (surrogate)then
                        CALL calcCCH_newshift(calc_shift,SURR_TIMES,
     +                       SPIKESHIFT,REF,TAR,SURR_ITAL,BINWs,SP,
     +                       0.0,LAST_REF,LAST_TAR,IDs,spikecount(i),i)
                     ELSE IF (cth_cch)then
                        if(i.eq.1)call rebin_cth_cch  (BINWs, SP)
                        if(i.eq.2)call rebin_var      (BINWs, SP)
                     END IF                    
                     NHW=INT(BINWs*50.)
                     ICN(1)=REF
                     ICN(2)=LAST_REF
                     ICN(3)=TAR
                     ICN(4)=LAST_TAR
                     ICN(5)=E_pulse
                     ICN(6)=ITAL(IDs(E_pulse))
                     if(i.eq.1)then
                        write (m+200,c_format,REC=REC_num_control_1)
     +                       1, NHW,ICN,SP !0=cross,1=cont cross,2=auto,cth or pst
                        histogram_number = histogram_number + 1 !update counter
                        call percent_completed(previous_percent,fildes, !recalculate % of histograms completed
     +                       histogram_number,total_histograms,
     +                       total_num_qdts)
                     end if
                     if(i.eq.2)then
                        write (m+200,c_format,REC=REC_num_control_2) 
     +                       1, NHW,ICN,SP !0=cross,1=cont cross,2=auto,cth or pst
                        histogram_number = histogram_number + 1 !update counter
                        call percent_completed(previous_percent,fildes, !recalculate % of histograms completed
     +                       histogram_number,total_histograms,
     +                       total_num_qdts)
                     end if
                     if(i.gt.1)then
                        do n = 1, NUM_BINS
                           CONTROL(m,izz,n) = CONTROL(m,izz,n) + 
     +                          SP(n) !add this SP CCH to the others to create the average CONTROL SP
                        end do
                     end if
                     if(i.eq.jTOTAL_SHIFTS+1)then
                        do n = 1, NUM_BINS
                           CONTROL(m,izz,n)=
     +                          NINT(FLOAT(CONTROL(m,izz,n)) /
     +                          FLOAT(jTOTAL_SHIFTS)) !compute the average shift-control CCH
                        end do
                        write (m+200,c_format,REC=REC_num_control_avg) 
     +                       1, NHW,ICN,
     +                       (CONTROL(m,izz,n),n=1,NUM_BINS) !0=cross,1=cont cross,2=auto,cth or pst
                        histogram_number = histogram_number + 1 !update counter
                     end if
                     call percent_completed(previous_percent,fildes, !recalculate % of histograms completed
     +                    histogram_number,total_histograms,
     +                    total_num_qdts)
                  end do
               end do
            end do 
            if(cth_cch)call free_cth_cch ()
         end do                 !GET_TAR_1        !get next TARget cell
      end do                    !GET_REF_1          !get next REFerence cell

*
*         
*
*
*       ***** Write gamesave data to *.qdt.sav: *****
*
*
      SPC_array = MAX_NUM_CHAN*MAX_NUM_ACC_CYCLES

      print '(/,''Writing qdt save file(s) (version '',A,
     +'') . . .'')',version
      do i = 1, total_num_qdts
         QDTSAV=qdt_files(i)//'.sav'
         call remove_all_blanks(QDTSAV,LEN(QDTSAV))

         call write_QDTSAV(version,QDTSAV,ios,
     +        date_exp,recording,protocol,BDT_FILE,qdt_files(i),
     +        IDs,excluded,included,total_num_cells,
     +        I_pulse,E_pulse,phrenic,BNDRY,cardiac_pls,
     +        BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +        BWs(i,1),BWs(i,2),BWs(i,3),BWs(i,4),BOUNDARY,
     +        start_time,end_time,icycles,
     +        total_histograms,ETA2_1,ETA2_2,ETA2_3,ETA2_4,
     +        ETA2_5,ETA2_6,zmodsig_1,zmodsig_2,zmodsig_3,
     +        zmodsig_4,zmodsig_5,zmodsig_6,zmodsig2_1,
     +        zmodsig2_2,zmodsig2_3,zmodsig2_4,zmodsig2_5,
     +        zmodsig2_6,coef,coefnum,card_type,exp_name,
     +        DELTA2,tedfactor,meanISI,sdISI,fiveHT,mean_rISI,sd_rISI,
     +        num_rej_ISI,num_rej_rISI,c_MAX_INT,jTOTAL_SHIFTS,
     +        num_acc_cycles,sp_per_cycle,ITAL)

         QDTTXT = qdt_files(i)//'.txt'
         call remove_all_blanks(QDTTXT,LEN(QDTTXT))
         call write_QDTTXT(qdt_files(i),shortest,
     +        longest)
      end do

      call percent_completed(previous_percent,fildes, !100% of histograms completed
     +     histogram_number,total_histograms,total_num_qdts)
      call FDATE(TODAY)
      text='Finish time: '//TODAY
      call strlength(text,LEN(text),l)
      call text2d(fildes,600.,180.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      text='File(s) successfully generated:  '//
     +     qdt_files(1)(1:INDEX(qdt_files(1),' ')-1)//char(0)
      call text2d(fildes,150.,300.,text,
     +     ANNOTATION_TEXT,1)
      do i = 2, total_num_qdts
         text=', '//qdt_files(i)(1:INDEX(qdt_files(i),' ')-1)//char(0)
         call append_text(fildes,text,ANNOTATION_TEXT,1)
      end do
      text=' '//char(0)
      call append_text(fildes,text,ANNOTATION_TEXT,0)
      call make_picture_current(fildes)


*
*       *********************************************************************
*       *********************************************************************
*       *********************************************************************
      if(total_num_qdts.gt.1)then
         print '(20(/),''The following QDT files were ''
     +''successfully generated:'',/)'
         do i = 1, total_num_qdts
            print '(T10,A)',qdt_files(i)
         end do
      else
         print '(20(/),''The following QDT file was ''
     +''successfully generated:'',//,T10,A)',qdt_files(1)
      end if           
      print *,"Waiting for surrogate generation to complete"
      call sc_finish ()
      print '(//,''Press <cr> to continue  >> ''$)'
      read (*,'(A)')
      call clear_all(fildes)
      do i = 1,total_num_qdts
         CLOSE(UNIT=i+200)
      end do

*       ***** close the main window *****
 500  retval = gclose(fildes)
*

      RETURN
*
*       *** END END END END END END END END END END END END END ***
*
      END
*
*
*       ************************************************************
*       ************************************************************
*
*       flag = 1 ---> bad IDcode
*
      subroutine CHECK_CODE(i,flag,IDs,display_message)
*
*
      INCLUDE 'x2000parameter.defs'
      dimension IDs(MAX_NUM_CODES)
      integer flag
      character*1 display_message
*
*
      flag = 0                  !clear the flag
      if(IDs(i).eq.0) then
         flag=1
         if(display_message.ne.'n')
     +        print '(''***** '',I3,'' IS AN INVALID CODE ..... ''
     +''please re-enter *****'')',i
      end if
*
      return
      end
*
*
*       ************************************************************
*       ************************************************************
*
*
      subroutine percent_completed(previous_percent,fildes,
     +     histogram_number,total_histograms,total_num_qdts)
      use mod_clear_routines
      use mod_miscellaneous_subroutines
      use mod_new_draw_button
      

      include 'x2000parameter.defs'

      real x0,y0,x1,y1,x00
      INCLUDE 'gopen_type.defs'
      integer total_histograms, histogram_number,total_num_qdts,
     +     percent,previous_percent
      character*120 text
      character*30 perimeter
      character*6 c_total_hist
      character*6 c_hist_num
      character*3 c_percent
      real cputime
      real :: last_cputime = 0
      save last_cputime

      if(.false.)print *,fildes2 !suppress unused variable warning

      call cpu_time(cputime)
      if ((cputime - last_cputime < 1.).and.
     +     (histogram_number.lt.total_histograms*total_num_qdts))
     +     return
      last_cputime = cputime

      x0 = 50.0                 !establish dimensions of "percent completed" box
      y0 = 20.0
      x1 = 670.0
      y1 = 60.0
      perimeter='perimeter'
      call perimeter_color(fildes,0.,0.,0.)
      call draw_button(fildes,x0,y0,x1,y1,perimeter,0.,0.," ",0.,0.) !draw the "percent completed" box
      call character_height(fildes,.060)
      call character_width(fildes,.015)
      write (c_total_hist,'(I6)')total_histograms*total_num_qdts
      call strlength(c_total_hist,LEN(c_total_hist),l)
      text=c_total_hist(1:l)//' histograms will be calculated'
      call strlength(text,LEN(text),l)
      call text2d(fildes,30.,100.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)

      percent = IFIX((histogram_number)/ !recalculate % of histograms completed
     +     FLOAT(total_histograms*total_num_qdts)*100.)
      call clear_quarter(fildes)
      write (c_percent,'(I3)') percent
      
      text=c_percent//'%  complete'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x1+100.,y0+20.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      write (c_hist_num,'(I6)') histogram_number
      text = '('//c_hist_num
      call remove_all_blanks(text,LEN(text))
      call strlength(text,LEN(text),l)
      text=text(1:l)//' histograms have been calculated)'
      call strlength(text,LEN(text),l)
      call text2d(fildes,x1+100,y0,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call make_picture_current(fildes)
      do ipct=previous_percent+5,percent,5
         call interior_style(fildes,INT_SOLID,1)
         x00 = (x0+10.)+(((ipct/5.)-1.)*30.) ! calculate location of lower left corner
!   of block
         call rectangle(fildes,x00,(y0+5.),(x00+30.),(y1-5.)) !draw the block
         call interior_style(fildes,INT_HOLLOW,1)
      end do
      call make_picture_current(fildes)
      previous_percent=percent/5*5  !reset the comparison percent value

      return
      end

      SUBROUTINE phase_bins (begin_E, middle_I, end_E, NORM_BW,
     +     TOTAL_SELECTED_CYCLES, ebins, ibins)
      include 'x2000parameter.defs'
      double precision begin_E(:), end_E(:), middle_I(:)
      real NORM_BW
      INTEGER*4 TOTAL_SELECTED_CYCLES, ebins, ibins
      double precision esum, isum, emean, imean
      esum = 0
      isum = 0
      DO i = 1, TOTAL_SELECTED_CYCLES
         esum = esum + (middle_I(i) - begin_E(i))
         isum = isum + (end_E(i) - middle_I(i))
      END DO
      emean = esum / TOTAL_SELECTED_CYCLES
      imean = isum / TOTAL_SELECTED_CYCLES
      ebins = NINT (emean / NORM_BW)
      ibins = NINT (imean / NORM_BW)
      END SUBROUTINE
      end module mod_calc_and_write_newshift_cccs
