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

      module mod_stats6
      contains
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                                       *
*       *  filename: stats6.f                                                   *
*       *  date of last revision  14-Mar-2004 lss                               *
*       *                                                                       *
*       *  This subroutine of x2000 allows the user to indicate                 *
*       *       inclusive bin boundaries of primary features and also to        *
*       *       indicate the boundaries of the bins of background activity      *
*       *       that will be used to calculate k values and their               *
*       *       significances as well as some other indices of correlation.     *
*       *                                                                       *
*       *  New with version 2: "background" values used to calculate the k value*
*       *       are obtained from the                                           *
*       *       shift-control CCH ...OR... using a range of bins selected       *
*       *       by the user from a region of the original CCH.                  *
*       *                                                                       *
*       *                                                                       *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
      SUBROUTINE statistics(fildes,mouse,CCH,CONTROL,
     +     scaled_hist,BINVAL,
     +     ZK,P,det,vis,ZLAT,HALFWD,info,
     +     pixels_per_bin,region,NHW,REFcode,TARcode,
     +     REFtype,TARtype,
     +     QDT_FILENAME,BDT_FILE,
     +     WINDOW,mode,Q_pos,statcomm,
     +     czk,cprobk,cdet,cvis,chalfwd,czlat,
     +     rec_num,rec_num_control,
     +     show_control,show_conf_lim,show_single_shift,
     +     show_avg_shift,show_2_sd,show_3_sd,
     +     per_text_abbrev,
     +     per_results_REF,per_results_TAR,date,recording,
     +     ICN,title,CONTROL_cl,
     +     min_bin,max_bin,scaledup,rmagnify,IDs,ITAL,
     +     CCH_scaled,CONTROL_scaled,CONTROL_cl_scaled,
     +     stats_bw)
*
      use mod_clear_routines
      use mod_locate_region
      use mod_miscellaneous_subroutines
      use mod_new_draw_button
      use mod_print_and_write_routines
      INCLUDE 'x2000parameter.defs'

      real bkgrnd      !what was Lauren thinking?
      real pkgrnd,effect,contrib
      integer CCH(101),CONTROL(101),
     +     CCH_scaled(101),CONTROL_scaled(101),
     +     scaled_hist(101),CONTROL_cl(101),
     +     CONTROL_cl_scaled(101),
     +     REFcode,TARcode,NHW,Q_pos,
     +     rec_num,rec_num_control,ICN(6),IDs(MAX_NUM_CODES),
     +     ITAL(MAX_NUM_CHAN)
*
      real BINVAL,pixels_per_bin,x1,x2
*
      character*(*) BDT_FILE,QDT_FILENAME
      character*120 text,text2,title
      character*40 statcomm,bw_string
      character*15 REFtype,TARtype,screen
      character*12 per_text_abbrev(MAX_PERTURB)
      character*11 date
      character*5 WINDOW
      character*9 per_results_REF(MAX_PERTURB),
     +     per_results_TAR(MAX_PERTURB)
      character*8 czk,cprobk,cdet,cvis,czlat,chalfwd
      character*8 czk_tmp,cprobk_tmp,cdet_tmp,cvis_tmp,czlat_tmp,
     +     chalfwd_tmp
      character*8 stats_bw,stats_bw_tmp
      character*6 cancel,repeat
      character*3 CCHtype,OK,c_IPKF,c_IBKF
      character*2 mode,recording
      character*1 info,task,show_control,show_conf_lim,
     +     show_single_shift,
     +     show_avg_shift,show_2_sd,show_3_sd,scaledup
*
      integer*4 mouse,region,region1
      INCLUDE 'gopen_type.defs'

      real :: BKSUM = 0.0
      integer :: ISPBIN = 0
      integer STATS_CONTROL(101),status
*
*
*       ***** suppress unused variable warnings *****

      if(.false.)print *,mouse
      if(.false.)print *,q_pos
      if(.false.)print *,window
      if(.false.)print *,fildes2

      info = 'n'
      CCHtype = 'CCH'
      OK = 'OK'
      repeat = 'REPEAT'
      cancel = 'CANCEL'
      IPKS=0                    !initialize variables
      IPKF=0
*
*       ***** Have user select the range of bins which encompass the primary feature: *****
*        
      call clear_bottom(fildes)
      call character_height(fildes,.045)
      call character_width(fildes,.015)
      text='SELECT FIRST BIN OF PRIMARY FEATURE RANGE'
      call text2d(fildes,400.,55.,text//char(0),ANNOTATION_TEXT,0)
      call make_picture_current (fildes)
      call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
      call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc)
      IPKS=INT((xloc-500.)/pixels_per_bin)+1
      if(IPKS.lt.2)IPKS=2
      call clear_bottom(fildes) !clear bottom of window
*
      text='SELECT LAST BIN OF PRIMARY FEATURE RANGE'
      call text2d(fildes,400.,55.,text//char(0),ANNOTATION_TEXT,0)
      call make_picture_current(fildes)
      call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
      call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc)
      IPKF=INT((xloc-500.)/pixels_per_bin)+1
      call clear_bottom(fildes) !clear bottom of window
*
*       ***** Selection of "background" range of bins: *****
*

      call sc_hist(STATS_CONTROL,status)
      if (status .ne. 0) then
         STATS_CONTROL = CONTROL
         print *,"using cth control for stats"
      else
         print *,"using surrogate control for stats"
      end if

      if(.false.)then
*
      else if(region.eq.12)then !calculate stats using bin values in entire shift-control CCH 
         BKSUM = 0.0            ! use corresponding range of bins in the shift-control as "background"
         IBKS=IPKS
         IBKF=IPKF
         do i = IBKS,IBKF
            BKSUM = BKSUM + STATS_CONTROL(i)
         end do

      else if(region.eq.13)then !calculate stats using values in a user-selected range of 
!  bins within the original CCH
*
         text='SELECT FIRST BIN OF BACKGROUND RANGE'
         call character_height(fildes,.045)
         call character_width (fildes,.015)
         call text2d(fildes,400.,55.,text//char(0),ANNOTATION_TEXT,0)
         call make_picture_current (fildes)
         call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
         call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc)
         IBKS=INT((xloc-500.)/pixels_per_bin)+1
         if(IBKS.lt.2)IBKS=2
         call clear_bottom(fildes) !clear bottom of window
*
         text='SELECT LAST BIN OF BACKGROUND RANGE'
         call text2d(fildes,400.,55.,text//char(0),ANNOTATION_TEXT,0)
         call make_picture_current(fildes)
         call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
         call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc)
         IBKF=INT((xloc-500.)/pixels_per_bin)+1

         BKSUM=0.0
         DO I=IBKS,IBKF
            BKSUM=BKSUM+CCH(I)  !total of all background bins
         end do

      end if

*               ***** outline the bins used as feature range for statistics in red: *****
*
      call fill_color(fildes,1.0,0.0,0.0) !fill color = red
      call text_color (fildes,1.0,0.0,0.0)
      call interior_style(fildes,INT_SOLID,1)
      call line_type(fildes,SOLID)
      call rectangle(fildes,
     +     (500.+((IPKS-1)*pixels_per_bin)),
     +     350.,(500.+(IPKF*pixels_per_bin)),360.)
      write (text,'(I3)') IPKS
      write (c_IPKF,'(I3)') IPKF
      call character_height(fildes,.030)
      call character_width(fildes,.010)
      call text2d(fildes,(450.+((IPKS-1)*pixels_per_bin)),
     +     350.,text//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,(500.+(float(IPKF)*pixels_per_bin)),
     +     350.,c_IPKF//' (selected feature)'//
     +     char(0),ANNOTATION_TEXT,0)
*
      call line_color (fildes,1.0,0.0,0.0)
      x1 = 500. + ((IPKS-1)*pixels_per_bin)
      call move2d (fildes,x1,(scaled_hist(IPKS-1)+300.))
      call draw2d (fildes,x1,(scaled_hist(IPKS)+300.))
      do i = IPKS,IPKF          !outline selected feature bins in red
         x2 = x1 + pixels_per_bin
         call move2d (fildes,x1,(scaled_hist(i)+300.))
         call draw2d (fildes,x2,(scaled_hist(i)+300.))
         call draw2d (fildes,x2,(scaled_hist(i+1)+300.))
         x1 = x2
         print '(''number of events in feature bin #'',I3,
     +'' = '',I10)',i,CCH(i)
       end do
      call make_picture_current(fildes)
      call interior_style(fildes,INT_HOLLOW,1)
      call line_color (fildes,0.,0.,0.) !return line color to black
      call text_color (fildes,0.,0.,0.) !return text color to black


*               ***** outline the bins used as background range for statistics in blue: *****
*
      call fill_color(fildes,0.0,0.0,1.0) !fill color = blue
      call text_color(fildes,0.0,0.0,1.0)
      call line_color (fildes,0.0,0.0,1.0)
      call interior_style(fildes,INT_SOLID,1)
      call rectangle(fildes,
     +     (500.+((IBKS-1)*pixels_per_bin)),
     +     330.,(500.+(IBKF*pixels_per_bin)),340.)
      call make_picture_current(fildes)
      call interior_style(fildes,INT_HOLLOW,1)
      write (text,'(I3)') IBKS
      write (c_IBKF,'(I3)') IBKF
      call character_height(fildes,.030)
      call character_width(fildes,.010)
      call text2d(fildes,(450.+((IBKS-1)*pixels_per_bin)),
     +     330.,text//char(0),ANNOTATION_TEXT,0)
      call text2d(fildes,(500.+(float(IBKF)*pixels_per_bin)),
     +     330.,c_IBKF//' (background)'//
     +     char(0),ANNOTATION_TEXT,0)

      if(region.eq.13)then      !outline background bins in blue only if different from
         call line_color (fildes,0.0,0.0,1.0) ! region of interest (with the feature)
         x1 = 500. + ((IBKS-1)*pixels_per_bin)
         call move2d (fildes,x1,(scaled_hist(IBKS-1)+300.))
         call draw2d (fildes,x1,(scaled_hist(IBKS)+300.))
         do i = IBKS,IBKF       !outline selected feature bins in red
            x2 = x1 + pixels_per_bin
            call move2d (fildes,x1,(scaled_hist(i)+300.))
            call draw2d (fildes,x2,(scaled_hist(i)+300.))
            call draw2d (fildes,x2,(scaled_hist(i+1)+300.))
            x1 = x2
            print '(''number of events in background bin #'',I3,
     +'' = '',I10)',i,STATS_CONTROL(i)
         end do
         call make_picture_current(fildes)
         call interior_style(fildes,INT_HOLLOW,1)
      end if

      call line_color (fildes,0.,0.,0.) !return line color to black
      call text_color (fildes,0.,0.,0.) !return text color to black



*
*
*
*       *****  calculate the sum of the bins in "feature" range:  *****
*
*
      PKSUM=0.0
*
      DO i=IPKS,IPKF
         PKSUM = PKSUM + CCH(i) !total of all primary feature bins
      end do
*
      NMPKBN=IPKF-IPKS+1        !total number of bins in feature
      NMBKBN=IBKF-IBKS+1
*
      bkgrnd=BKSUM/NMBKBN      !mean value of background bins
      pkgrnd=PKSUM/NMPKBN      !mean value of feature bins
      BACKCT=bkgrnd*NMPKBN     !=(mean background bin)*(# of feature bins)
! = a calculated count of the background events that underlie the primary feature

      print '(''BKSUM = '',f10.2,''; NMBKBN = '',I5)', BKSUM,NMBKBN
      print '(''PKSUM = '',f10.2,''; NMPKBN = '',I5)', PKSUM,NMPKBN
      print '(''mean background (bkgrnd) = '',f10.2,
     +''; mean feature (pkgrnd) = '',f10.2)',bkgrnd,pkgrnd
      print '(''BACKCT = '',f10.2)',BACKCT

*
*
*
*       *****  find the greatest and least bin values in the primary feature: *****
*
*
      ILOW=999999               !load variables with dummy values
      IHI=0
*
      DO I=IPKS,IPKF
         IHI=MAX0(IHI,CCH(I))   !IHI = max bin value
         ILOW=MIN0(ILOW,CCH(I)) !ILOW = min bin value
      end do
*
*       define ISPBIN pointer
*
*
      if (pkgrnd.eq.bkgrnd) then !no difference between mean feature
! and background values - cannot analyze these ranges
!go back and try again
         call clear_bottom(fildes)
         call character_height(fildes,.054)
         call character_width(fildes,.018)
         text='MEAN VALUES OF FEATURE AND CONTROL RANGES ARE EQUAL'
         call text2d(fildes,100.,120.,text//char(0),ANNOTATION_TEXT,0)
         text =' TRY AGAIN'
         call character_height(fildes,.030)
         call character_width(fildes,.010)
         call draw_button (fildes,525.,5.,725.,80.,text,530.,35.,'',0.,
     +        0.) 
         text='  CANCEL'
         call draw_button (fildes,775.,5.,975.,80.,text,780.,35.,'',0.,
     +        0.)
 15      call make_picture_current (fildes)
         call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
         call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc)
         call locate_region(xloc,yloc,region1,'statsub')
*
         if(.false.)then
*                               !wait for a valid click
         else if(region1.eq.12)then !redraw histogram and try it again
            call clear_all(fildes)
            info='r'
            return
*
         else if(region1.eq.13)then !cancel
            call clear_all(fildes)
            return
*
         else
            goto 15
*
         end if
*
      end if
*
      if(pkgrnd.gt.bkgrnd) then !PEAK
         do i=IPKS,IPKF
            if (CCH(i).eq.IHI) ISPBIN=i !ISPBIN = # of the tallest bin
         end do
      end if
*
      if (pkgrnd.lt.bkgrnd) then !TROUGH
         do i= IPKS,IPKF
            if (CCH(i).eq.ILOW) ISPBIN=i !ISPBIN = # of the lowest bin
         end do
      end if


*
*       **** calculate detectability index (Aertsen & Gerstein                  ****
*       **** Brain Res. 340:341-354, 1985; Melssen & Epping, Biol. Cybernetics  ****
*       **** 57:403-414,1987):                                                  ****
*
*
      smalld=abs(CCH(ISPBIN)-bkgrnd) !difference between value of tallest or
! lowest bin and the mean background
      vis_tmp=smalld/bkgrnd    !visibility index = s
      det_tmp=smalld/sqrt(bkgrnd) !detectability index = d
*
*
*       K VALUE HERE = ZK
*
      ZK_tmp = PKSUM/BACKCT     !K-VALUE.... >1 FOR PEAK..<1 FOR TROUGH
*
      call clear_bottom(fildes) !clear bottom of window
*
*
*       *****  calculate the time lag from 0 to the primary feature:  *****
*
*
      if(((IPKS.le.51).and.(IPKF.gt.51)).or.(IPKF.eq.51))then !central feature
         ZLAT_tmp = 0.0 
      end if
*
      if (IPKF.lt.51) then      !offset to the left
         ZLAT_tmp=(IABS(51-IPKF))*BINVAL !=time from 0 to the final bin of the feature
      end if
*
      if (IPKS.gt.51) ZLAT_tmp=(IABS(IPKS-51))*BINVAL !offset to the right
!ZLAT=time from 0 to the beginning bin of the feature


*       **** calculate the measures of synaptic effectiveness and contribution ****
*       **** as evidenced by the primary feature as in                         ****
*       **** Levick et al., Invest Ophthalmol. 11:302-311, 1972;               ****
*       **** Lindsey and Gerstein, J Neurophysiol. 42:383-399, 1979.           ****
*       **** Aertsen and Gerstein, Brain Res. 340:341-354, 1985                ****
*       ****                                                                   ****
*       **** These values are reported/saved for offset peaks only.            ****

      if((ZLAT_tmp.gt.0.0).and.(pkgrnd.gt.bkgrnd))then !offset RIGHT peak
         if(NMBKBN.eq.NMPKBN)then !be sure using same number of background bins as in peak
            effect = (PKSUM-BKSUM)/ITAL(IDs(REFcode))
            contrib = (PKSUM-BKSUM)/ITAL(IDs(TARcode))
         else
            effect = (PKSUM-BACKCT)/ITAL(IDs(REFcode))
            contrib = (PKSUM-BACKCT)/ITAL(IDs(TARcode))
         end if
      elseif((ZLAT_tmp.lt.0.0).and.(pkgrnd.gt.bkgrnd))then !offset LEFT peak
         if(NMBKBN.eq.NMPKBN)then !be sure using same number of background bins as in peak
            effect = (PKSUM-BKSUM)/ITAL(IDs(TARcode))
            contrib = (PKSUM-BKSUM)/ITAL(IDs(REFcode))
         else
            effect = (PKSUM-BACKCT)/ITAL(IDs(TARcode))
            contrib = (PKSUM-BACKCT)/ITAL(IDs(REFcode))
         end if
         print '(''REF = '',I3,''; ITAL('',I3,'') = '',I10)',
     +     REFcode,IDs(REFcode),ITAL(IDs(REFcode))
         print '(''TAR = '',I3,''; ITAL('',I3,'') = '',I10)',
     +     TARcode,IDs(TARcode),ITAL(IDs(TARcode))
         print '(''effectiveness = '',f7.5,''; contribution = '',f7.5)',
     +     effect,contrib
      end if

*
*
*       ***** CALCULATE SIGNIFICANCE OF K VALUE AFTER SEARS,STAGG,KIRKWOOD: *****
*
*
      P_tmp=1.0                 !VALUE KEPT IF NOT SIG.
*
      ZMPK=PKSUM/NMPKBN         !MEAN FEATURE BIN COUNT
      ZMBK=BKSUM/NMBKBN         !MEAN BACKGROUND BIN COUNT
      SQ=SQRT(ZMBK)             !SQR ROOT OF BACKGROUND BIN COUNT
*
      P05=1.96*SQ               !probability values for 1 std dev
      P01=2.58*SQ               ! 2 std devs
      P001=3.28*SQ              ! 3 std devs
*
      IF((ZMPK.GT.(ZMBK+P05)).OR.(ZMPK.LT.(ZMBK-P05))) P_tmp=.05
      IF((ZMPK.GT.(ZMBK+P01)).OR.(ZMPK.LT.(ZMBK-P01))) P_tmp=.01
      IF((ZMPK.GT.(ZMBK+P001)).OR.(ZMPK.LT.(ZMBK-P001))) P_tmp=.001 !very large feature
*
*
*       ***** CALCULATE HALF WIDTH: *****
*
*
      ICRTCT=0                  !tally of # of bins that are at least half as tall/low as
! the maximum peak/trough bin minus background
*
      IF(ZK_tmp.GT.1.) then     !peak
         CRITZ=(IHI-ZMBK)*.5    != half of the value by which the peak is greater
!  than background
         DO I=IPKS,IPKF
            IF((CCH(I)-ZMBK).GE.CRITZ) ICRTCT=ICRTCT+1
            IF(((CCH(I)-ZMBK).LT.CRITZ).AND.(ICRTCT.GT.0)) exit !end of peak
         end do
         HALFWD_tmp=ICRTCT*BINVAL !HALF WIDTH FOR PEAK = duration of peak
      end if
*
*
      IF(ZK_tmp.LT.1.) then     !trough
         CRITZ=(ZMBK-ILOW)*.5   != half of the value by which the trough is less
!  than background
         DO I=IPKS,IPKF
            IF((ZMBK-CCH(I)).GE.CRITZ) ICRTCT=ICRTCT+1
            IF(((ZMBK-CCH(I)).LT.CRITZ).AND.(ICRTCT.GT.0)) exit !end of trough
         end do
         HALFWD_tmp=ICRTCT*BINVAL !HALF WIDTH FOR TROUGH = duration of trough
      end if
*
      IF(ZK_tmp.EQ.1.) HALFWD_tmp=0.0 !no significant feature in this histogram
*
 90   call clear1(fildes,2)     !clear bottom of window
      text='NOTE: +INF values not allowed!'
      call text2d(fildes,800.,125.,text//char(0),ANNOTATION_TEXT,0)
      text='k value'
      call text2d(fildes,5.,125.,text//char(0),ANNOTATION_TEXT,0)
      text='detectability index (d.i.)'
      call text2d(fildes,5.,105.,text//char(0),ANNOTATION_TEXT,0)
      text='visibility index'
      call text2d(fildes,5.,85.,text//char(0),ANNOTATION_TEXT,0)
      text='significance of k value'
      call text2d(fildes,5.,65.,text//char(0),ANNOTATION_TEXT,0)
      text='time lag to feature (ms)'
      call text2d(fildes,5.,45.,text//char(0),ANNOTATION_TEXT,0)
      text='half width of feature (ms)'
      call text2d(fildes,5.,25.,text//char(0),ANNOTATION_TEXT,0)
      text='binwidth (ms)'
      call text2d(fildes,5.,5.,text//char(0),ANNOTATION_TEXT,0)

*
*       ***** print interaction strength & detect indices: *****
*
*               ZK  = k value (F5.2)
*               det = detectability index = d (F4.1)
*               vis = visibility index = s (F5.2)
*
      write (czk_tmp,'(F6.2)')ZK_tmp    
      write (stats_bw_tmp,'(F6.2)')BINVAL
      write (cprobk_tmp,'(F6.3)')P_tmp
      write (cdet_tmp,'(F6.2)')det_tmp
      write (cvis_tmp,('(F6.2)'))vis_tmp
      write (czlat_tmp,('(F6.2)'))ZLAT_tmp
      write (chalfwd_tmp,('(F6.2)'))HALFWD_tmp
      call character_height(fildes,.030)
      call character_width(fildes,.010)
      if((mode.eq.'ed').or.(mode.eq.'cr'))then
         call text2d(fildes,215.,150.,'CURRENT vs. STORED VALUES:'
     +        //char(0),ANNOTATION_TEXT,0)
      end if
      
      call text2d(fildes,250.,125.,czk_tmp//char(0),
     +     ANNOTATION_TEXT,0)
      call text2d(fildes,250.,105.,cdet_tmp//char(0),
     +     ANNOTATION_TEXT,0)
      call text2d(fildes,250.,85.,cvis_tmp//char(0),
     +     ANNOTATION_TEXT,0)
      call text2d(fildes,250.,65.,cprobk_tmp//char(0),
     +     ANNOTATION_TEXT,0)
      call text2d(fildes,250.,45.,czlat_tmp//char(0),
     +     ANNOTATION_TEXT,0)
      call text2d(fildes,250.,25.,chalfwd_tmp//char(0),
     +     ANNOTATION_TEXT,0)
      call text2d(fildes,250.,5.,stats_bw_tmp//char(0),
     +     ANNOTATION_TEXT,0)

      if((mode.eq.'ed').or.(mode.eq.'cr'))then
         call text2d(fildes,350.,125.,czk//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,350.,105.,cdet//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,350.,85.,cvis//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,350.,65.,cprobk//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,350.,45.,czlat//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,350.,25.,chalfwd//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,350.,5.,stats_bw//char(0),ANNOTATION_TEXT,0)
      end if

*
      if((mode.eq.'ed').or.(mode.eq.'cr'))then
         text=   '  KEEP STORED'
         text2 = '    VALUES'
         call draw_button(fildes,525.,5.,725.,80.,
     +        text,530.,45.,text2,530.,25.)
         text=   '  KEEP CURRENT' 
         text2 = '    VALUES'
         call draw_button(fildes,525.,90.,725.,165.,
     +        text,530.,130.,text2,530.,110.)
      end if
      text='   PRINT'
      call draw_button(fildes,775.,5.,975.,80.,text,780.,35.,'',0.,0.)    
      text='   WRITE'
      call draw_button(fildes,1025.,5.,1225.,80.,text,1030.,35.,'',0.,
     +     0.) 
      text='  CONTINUE' 
      call draw_button(fildes,1295.,5.,1495.,80.,text,1300.,35.,'',0.,
     +     0.)

 100  call make_picture_current(fildes)
      call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
      call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc) !translate into fp coords.
      screen = 'statsub'
      call locate_region(xloc,yloc,region1,screen)
*
      if(region1.eq.12)then     !CANCEL - do not save the newly calculated statistical values (KEEP STORED VALUES)
         if((mode.eq.'ed').or.(mode.eq.'cr'))then
            call text_color(fildes,1.0,1.0,1.0) !white text - will cause a "flash" when this button pressed
            text=   '  KEEP STORED'
            text2 = '    VALUES'
            call draw_button(fildes,525.,5.,725.,80.,
     +           text,530.,45.,text2,530.,25.)
            call make_picture_current(fildes)
            do i = 1, 10000000
            end do
            call text_color(fildes,0.0,0.0,0.0) !back to black text
            call draw_button(fildes,525.,5.,725.,80.,
     +           text,530.,45.,text2,530.,25.)
            call make_picture_current(fildes)
            goto 100
         else
            goto 90
         end if
      else if(region1.eq.21.or.region1.eq.121.or.region1.eq.221)then !keep current values
         if((mode.eq.'ed').or.(mode.eq.'cr'))then
            call text_color(fildes,1.0,1.0,1.0) !white text - will cause a "flash" when this button pressed
            text=   '  KEEP CURRENT' 
            text2 = '    VALUES'
            call draw_button(fildes,525.,90.,725.,165.,
     +           text,530.,130.,text2,530.,110.)
            call make_picture_current(fildes)
            do i = 1, 10000000  !flash
            end do
            ZK = ZK_tmp
            det = det_tmp
            vis = vis_tmp
            P = P_tmp
            ZLAT = ZLAT_tmp
            HALFWD = HALFWD_tmp
            czk = czk_tmp
            cdet = cdet_tmp
            cvis = cvis_tmp
            czlat = czlat_tmp
            chalfwd = chalfwd_tmp
            cprobk = cprobk_tmp
            stats_bw = stats_bw_tmp
            call strlength(QDT_FILENAME,LEN(QDT_FILENAME),l)
            bw_string='bw='//stats_bw
            call remove_all_blanks(bw_string,LEN(bw_string))
            call strlength(bw_string,LEN(bw_string),m)
            statcomm=QDT_FILENAME(1:l)//
c     +                 ' @ bw = '//stats_bw
     +           ' @ '//bw_string(1:m)
            call text_color(fildes,0.0,0.0,0.0) !back to black text
            call draw_button(fildes,525.,90.,725.,165.,
     +           text,530.,130.,text2,530.,110.)
            call make_picture_current(fildes)
            goto 100
         else
            goto 90
         end if
      else if(region1.eq.13.or.region1.eq.14.or.region1.eq.213)then !PRINT or WRITE the graphics window
         if(region1.eq.13)task='p'
         if(region1.eq.14)task='w'
         if(region1.eq.213)task='v'
         if((mode.ne.'cr').and.(mode.ne.'ed'))then
            czk = czk_tmp
            cdet = cdet_tmp
            cvis = cvis_tmp
            czlat = czlat_tmp
            chalfwd = chalfwd_tmp
            cprobk = cprobk_tmp
            stats_bw = stats_bw_tmp
         end if
         if(scaledup.eq.'n')then
            call print_or_write_STATS(task,CCH,CONTROL,
     +           BDT_FILE,      !print the STATISTICS screen
     +           QDT_FILENAME,REFcode,TARcode,date,
     +           recording,NHW,
     +           rec_num,rec_num_control,title,ICN,
     +           show_control,show_conf_lim,
     +           REFtype,TARtype,
     +           show_single_shift,show_avg_shift,
     +           CONTROL_cl,show_2_sd,show_3_sd,
     +           czk,cdet,cvis,cprobk,czlat,
     +           chalfwd,
     +           IPKS,IPKF,IBKS,IBKF,per_text_abbrev,
     +           per_results_REF,per_results_TAR,
     +           min_bin,max_bin,scaledup,rmagnify,
     +           IDs,ITAL)
         else if(scaledup.eq.'y')then
            call print_or_write_STATS(task,CCH_scaled,
     +           CONTROL_scaled,BDT_FILE, !print the STATISTICS screen
     +           QDT_FILENAME,REFcode,TARcode,date,
     +           recording,NHW,
     +           rec_num,rec_num_control,title,ICN,
     +           show_control,show_conf_lim,
     +           REFtype,TARtype,
     +           show_single_shift,show_avg_shift,
     +           CONTROL_cl_scaled,show_2_sd,show_3_sd,
     +           czk,cdet,cvis,cprobk,czlat,
     +           chalfwd,
     +           IPKS,IPKF,IBKS,IBKF,per_text_abbrev,
     +           per_results_REF,per_results_TAR,
     +           min_bin,max_bin,scaledup,rmagnify,
     +           IDs,ITAL)
         end if

         call clear_bottom(fildes)
         goto 90
      else if(region1.eq.15)then !CONTINUE
         return
      else
         goto 100
      end if
*
*
      return
      end


      end module mod_stats6
