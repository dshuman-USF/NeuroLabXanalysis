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

      module mod_ps_plot_1
      contains
*
*       This is a new subroutine to draw histograms in a PS format for printing or saving.
*
*       filename = ps_plot_1.f
*
*       date of last revision = 16-Mar-2006     lss
*
*       link with x2000 code
*
*       The following parameters are passed into this subroutine:
*
*               the histograms to be plotted = IH, IH_overlay
*               (0,0) coordinates (x0,y0)
*               histogram type:  CCH, CTH, ACH, DIF, CCC (cardiac CCH)
*               overall dimensions of histogram (in pixels) -- height x width
*               histogram bin data (IH(101))
*               overlay bin data (IH_overlay(101)) -- (= IH(101) if no overlay)
*               whether or not the plot will contain labels and other info (for archiving and
*                       = info ('y' or 'n')                                 screen capture)
*               original *.bdt data file = BDT_FILE
*               *.qdt file that contains calculated histograms = QDT_FILENAME
*               date of experiment = date
*               which recording on that date = recording
*               number of this histogram within the *.qdt file = rec_num
*               record number of associated shift-control CCH = rec_num_control
*               title of histogram - = title
*               information re ref and tar codes and numbers of spikes = ICN()
*
*
      subroutine ps_plot(total_num_bins,
     +     IH,IH_overlay,IH1_scaled,
     +     hist_type,x0,y0,height,width,info,
     +     BDT_FILE,QDT_FILENAME,date,
     +     recording,NHW,rec_num,show_control,
     +     rec_num_control,title,ICN,max_bin,
     +     show_conf_lim,
     +     show_single_shift,show_avg_shift,
     +     IH_overlay_2,show_2_sd,show_3_sd,
     +     scaledup,min_bin,rmagnify,IDs,ITAL)
*
*
      use mod_miscellaneous_subroutines

c      use,intrinsic :: IEEE_ARITHMETIC
c     Quiet NaN, double precision
      REAL(8), PARAMETER :: D_QNAN =
     +     TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)
      
      include 'x2000parameter.defs'
      PARAMETER (iDOT=1,iSOLID=0)
      integer*4 total_num_bins
      integer IH(total_num_bins),IH_overlay(total_num_bins),
     +     IH_scaled(total_num_bins),
     +     ICN(:),IH1_scaled(total_num_bins),
     +     SP_DIFF(total_num_bins),
     +     IH_overlay_2(total_num_bins),IDs(MAX_NUM_CODES),
     +     ITAL(MAX_NUM_CHAN),sum_of_all_bins
*
      real x0,x1,y0,height,width,scale_factor,pixels_per_bin
      real binwd,tsd_plus(total_num_bins),tsd_minus(total_num_bins)
*
      integer rec_num,rec_num_control
*
      character*(*) title
      character*120 text
      character*(*) BDT_FILE,QDT_FILENAME
      character*(*) date
      character*10 c_rate,c_ITAL
      character*8 c_magnify
      character*7 c_ICN(6),c_binwd,c_max_bin_value,c_max,c_min,
     +     c_sum_of_all_bins
      character*5 c_rec_num,c_rec_num_control
      character*4 c_NHW
      character*(*) hist_type
      character*(*) recording
      character*(*) info,show_control,show_conf_lim,
     +     show_single_shift,show_avg_shift,show_2_sd,
     +     show_3_sd,scaledup

      real :: tsd = 0.0
      integer*4 int
      real*4 flt
      equivalence (int, flt)

*       ***** suppress unused variable warnings *****
      if(.false.)print *,show_avg_shift
      if(.false.)print *,show_single_shift

      call ps_color_set (1)
      IH_scaled = 0             !initialize arrays
      IH1_scaled = 0

      pixels_per_bin =width/FLOAT(total_num_bins) !calculate width of bin
      binwd = (float(NHW)/50.)
      if(hist_type.ne.'CCH')scaledup='n'
*
*       *****************************
*       *****   draw the axes:  *****
*       *****************************
*
      call ps_line_color(0.,0.,0.) !black
      call ps_line_type (iSOLID)
      call ps_move2d (x0,(y0+height)) !Y axis
      call ps_draw2d (x0,y0)
      if(hist_type.ne.'DIF')then
         call ps_draw2d((x0+width),y0) !X axis
      end if
*
      if((hist_type.eq.'CCH').or.(hist_type.eq.'DIF'))then !draw dotted line to mark middle bin
         call ps_line_type (iDOT)
         call ps_move2d ((x0+(50.*pixels_per_bin)
     +        +(pixels_per_bin/2.)),y0)
         call ps_draw2d ((x0+(50.*pixels_per_bin)
     +        +(pixels_per_bin/2.)),(y0+height))
         call ps_line_type (iSOLID)
      end if
*
c       if(scaled.eq.'y')goto 100                       !do not re-scale this CCH -- NOT USED?!?!?!?
*
*
*       **************************************************************
*       ***** handle the special case of a difference histogram: *****
*       **************************************************************
*
*
      if(hist_type.eq.'DIF')then
         call ps_move2d(x0,(y0+(height/2.0)))
         call ps_line_type(iDOT) !draw and label horizontal dotted line
         call ps_draw2d((x0+width),(y0+(height/2.0))) !  at zero for difference histogram
         call ps_line_type(iSOLID)
         text='0'
         call ps_text2d((x0-25.),(y0+150.),text//char(0),
     +        PS_ANNOTATION_TEXT)
*
         max_bin_value = 0
*
         do i = 1,total_num_bins
            max_bin_value = MAX0(ABS(IH(i)),max_bin_value)
         end do
*
         if(max_bin_value.eq.0)then
            scale_factor=1.0
            goto 50
         end if
*
         scale_factor = (height/2.0)/FLOAT(max_bin_value)
*
 50      IH_scaled = IFIX(FLOAT(IH)*scale_factor)
*
         x1 = x0
         do i = 1,total_num_bins-1
            x2 = x1 + pixels_per_bin
            call ps_rectangle(x1,y0+(height/2.0),
     +           x2,(IH_scaled(i)+y0+(height/2.0)))
            x1 = x2
         end do
*
         text='DIFFERENCE HISTOGRAM: CCH - SHIFT CONTROL'
         call ps_character_height(.033)
         call ps_character_width(.011)
         call ps_text2d (x0-50.,(y0+height+40.),
     +        text//char(0),PS_ANNOTATION_TEXT)
*
         call ps_character_height(.030)
         call ps_character_width(.010)
         text='-'//c_NHW
         call ps_text2d(x0,y0-20.,text//char(0),PS_ANNOTATION_TEXT)
*
         text=c_NHW//' msec'
         call ps_text2d(x0+width-50.,y0-20.,text//char(0),
     +        PS_ANNOTATION_TEXT)
*
*
         return
      end if
*
*
*       *************************************
*       ***** scale the main histogram: *****
*       *************************************
*
      max_bin_value = 0
      sum_of_all_bins = 0
      do i = 1,total_num_bins
         max_bin_value = MAX0(max_bin_value,(IH(i)))
         sum_of_all_bins = sum_of_all_bins + IH(i)
      end do

      if(max_bin_value.eq.0)then !avoid a "divide by zero" situation
         scale_factor = 1.0
         goto 100
      end if
*
      scale_factor = height/FLOAT(max_bin_value)
*
 100  do i = 1, total_num_bins
         IH_scaled(i) = IFIX(FLOAT(IH(i))*scale_factor) !scale the main histogram
         IH1_scaled(i)=IH_scaled(i) !store the scaled histogram
      end do

      rate = D_QNAN
      if (size(ICN).ge.2) then
         if((ICN(2).eq.sum_of_all_bins).and.(ICN(1).gt.0))then
c        In versions of xanalysis prior to 6.5.5, there was a bug that
c        made LAST_TAR equal to the sum of all bins, causing the rate to
c        be way off for mirrored plots.  This code detects that bug in
c        data generated by those old versions and uses ITAL(tar_id)
c        instead of LAST_TAR, and since those two numbers might be
c        slightly different, it prints a warning.  This code should
c        never run for data generated by newer versions.
            call ps_character_height(.030)
            call ps_character_width(.010)
            rate = (float(max_bin_value)*
     +           (1000.0/binwd))/ITAL(IDs(ICN(1))) !this plot is mirrored
            text = 'The rate was calculated '
            call strlength(text,LEN(text),l)
            call ps_text2d(5.,y0+height-50.,text(1:l)//char(0),
     +           PS_ANNOTATION_TEXT)
            text = 'using the total number of events'
            call strlength(text,LEN(text),l)
            call ps_text2d(5.,y0+height-70.,text(1:l)//char(0),
     +           PS_ANNOTATION_TEXT)
            text = 'in the reference spike train.'
            call strlength(text,LEN(text),l)
            call ps_text2d(5.,y0+height-90.,text(1:l)//char(0),
     +           PS_ANNOTATION_TEXT)
         else
            rate = (float(max_bin_value)*(1000.0/binwd))/ICN(2) !calculate the rate in spikes/sec
         end if
      end if



*
*
*
*       ************************************************************
*       *****   plot the data from the main histogram array:   *****
*       ************************************************************
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
*
*       ***** special case of respiratory cycle duration *****  !
*       *****   distribution histogram (hist_type=DIS):  *****  !
*                                                               !
      if(hist_type.eq.'DIS')then !
         call ps_fill_color(0.,0.,0.) !
         call ps_perimeter_color(0.,0.,0.) !
         call ps_interior_style(INT_SOLID,1) !
         x1 = x0                !               
         do i = 1,total_num_bins !
            x2 = x1 + pixels_per_bin !
            call ps_rectangle(x1,y0,x2,(IH_scaled(i)+y0)) !
            x1 = x2             !
         end do                 !
         goto 120               !
      end if                    !
*                                                               !
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
*
      x1 = x0 
      do i = 1,(total_num_bins-1)
         x2 = x1 + pixels_per_bin
         call ps_move2d (x1,(IH_scaled(i)+y0))
         call ps_draw2d (x2,(IH_scaled(i)+y0))
         call ps_draw2d (x2,(IH_scaled(i+1)+y0))
         x1 = x2
      end do

      if (hist_type.eq.'CCH') then
         call sc_plot_ps (fildes, x0, y0, pixels_per_bin, scale_factor)
      end if
*
*
*
*       *************************************************
*       *****   scale the overlay (if necessary):   *****
*       *************************************************
*
 120  if((hist_type.eq.'CTH').or.(show_control.eq.'y').or.
     +     (hist_type.eq.'CCC').or.
     +     (hist_type.eq.'DIS').or.
     +     (show_conf_lim.eq.'y'))then
         if((hist_type.eq.'CTH').or.(hist_type.eq.'CCC'))then !rescale the overlay for a CTH or cardiac pulse CCH
            max_bin_value = 0
            do i = 1,total_num_bins
               max_bin_value = MAX0(max_bin_value,IH_overlay(i))
            end do
         end if
*
         if(max_bin_value.eq.0)then
            scale_factor = 1.0
            goto 200
         end if
*
         scale_factor = height/FLOAT(max_bin_value)
*
 200     IH_scaled = IFIX(FLOAT(IH_overlay)*scale_factor) !scale the overlay
         do i = 1,total_num_bins
c             if((IH_scaled(i).gt.height).or.(IH_scaled(i).lt.0))
c     +            print '(''bin value error'')'
         end do
*
*
*       ********************************************************************
*       *****   plot the data from the overlay array (if necessary):   *****
*       ********************************************************************
*
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
*
*       ***** special case of respiratory cycle duration *****  !
*       *****   distribution histogram (hist_type=DIS):  *****  !
*                                                               !
         if(hist_type.eq.'DIS')then !
            call ps_perimeter_color(0,0.0,0.0) !outline selected bins in black
            call ps_interior_style(INT_SOLID,1) !
            call ps_fill_color(1.,1.,0.) !fill selected bins with yellow
            x1 = x0             !               
            do i = 1,total_num_bins !
               x2 = x1 + pixels_per_bin !
               if(IH_scaled(i).eq.0)goto 201 !bin value is 0 - don't draw it
               call ps_rectangle(x1,y0,x2,(IH_scaled(i)+y0)) !
 201           x1 = x2          !
            end do              !
            call ps_perimeter_color(0.,0.,0.) !default outline color = black
            call ps_line_type(iSOLID) !
            call ps_line_color(0.,0.,0.) !
            return              !
c         goto 206                                              !
         end if                 !
*                                                               !
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
*
         if(show_control.eq.'y')then
            call ps_line_type (iDOT)
            x1 = x0     
            do i = 1,(total_num_bins-1)
               x2 = x1 + pixels_per_bin
               call ps_move2d (x1,(IH_scaled(i)+y0))
               call ps_draw2d (x2,(IH_scaled(i)+y0))
               call ps_draw2d (x2,(IH_scaled(i+1)+y0))
               x1 = x2
            end do
            call ps_line_type (iSOLID)
            call ps_line_color(0.,0.,0.)
         end if
      end if
*
      if((show_control.eq.'y').and.(show_conf_lim.eq.'y'))then

*       ***** TO CALCULATE THE CONFIDENCE LIMITS:  *****
*       We were running into problems calculating the mean and standard deviation of data in which the cells were active
*         during different phases of the respiratory cycle, thereby producing a trend across the CCH ("tilts", "scoops", 
*         "hills", etc.).  Using a single shift-predictor (SP) control histogram to calculate the mean and sd wasn't working
*         because these trends obviously produced a mean which was more attuned to the trend than to shorter-term correlations
*         and the resulting confidence limits were unusable for statistical purposes.  So ...

*               1.  To remove possible long-term correlations or trends, we now use the histogram of the difference 
*                       between 2 SP CONTROL histograms (--> SP_DIFF = IH_overlay - IH_overlay_2)) to calculate the variance and 
*                       standard deviation of the control data.
*                       a. When viewing the single-shift control, use the "shift by 1" SP and the "shift by 2" SP to compute
*                               the difference histogram.
*                       b. When viewing the averaged-shift control, use the averaged-shift SP and the "shift by 1" SP.
*               2.  Because we're using a difference histogram to calculate the standard deviation of the control data,
*                       the mean of the difference histogram will be zero (or very close to it).  
*                       Therefore, technically speaking, the calculation of the mean could be left out of the 
*                       following calculations, but it is included because it is "expected" in the calculation of 
*                       the standard deviation.
*               3.  The standard deviation of the difference histogram is calculated across the 101 bins of the histogram.
*                       (one std dev value per difference histogram)
*               4.  The confidence limits (representing +/- 2 or 3 standard deviations) are then plotted around each bin of the
*                       appropriate SP CONTROL histogram (either the single-shift SP or the averaged-shift SP -- user's choice).

         if(show_2_sd.eq.'y')then
            imul = 2
            call ps_line_color(0.,0.,1.) !blue line
         else if(show_3_sd.eq.'y')then 
            imul = 3
            call ps_line_color(1.,0.,1.) !magenta line
         else
            STOP ', BUG: neither show_2_sd nor show_3_sd in ps_plot'
         end if

         if (qdt_file_version.eq.'8') then
            do i = 1, total_num_bins
               int = IH_overlay_2(i) !int and flt are equivalence'd
               tsd_plus (i)=ifix(IH_overlay(i)+imul*flt)*scale_factor
               tsd_minus(i)=ifix(IH_overlay(i)-imul*flt)*scale_factor
            end do
         else
            do i = 1, total_num_bins
               SP_DIFF(i) = IH_overlay(i) - IH_overlay_2(i)
            end do
            mean = 0
            sum = 0
            sum_sqr = 0.0
            sd_SP = 0.0
            do m = 1, total_num_bins
               sum = sum + SP_DIFF(m)
            end do
            mean = sum / total_num_bins !using a difference histogram, so mean is very close to zero
            do m = 1, total_num_bins
               sum_sqr = sum_sqr + float((mean-SP_DIFF(m))**2)
            end do
            sd_SP = SQRT(sum_sqr/(float(total_num_bins)-1.))
            tsd = imul*sd_SP

            do i = 1, total_num_bins
               tsd_plus(i) = ifix((IH_overlay(i) + tsd)*scale_factor)
               tsd_minus(i) = ifix((IH_overlay(i) - tsd)*scale_factor)
            end do
         end if

         call ps_line_type(iDOT)
         x1 = x0      
         do i = 1,(total_num_bins-1)
            x2 = x1 + pixels_per_bin
            call ps_move2d (x1,(tsd_plus(i)+y0))
            call ps_draw2d (x2,(tsd_plus(i)+y0))
            call ps_draw2d (x2,(tsd_plus(i+1)+y0))
            x1 = x2
         end do
         x1 = x0      
         do i = 1,(total_num_bins-1)
            x2 = x1 + pixels_per_bin
            call ps_move2d (x1,(tsd_minus(i)+y0))
            call ps_draw2d (x2,(tsd_minus(i)+y0))
            call ps_draw2d (x2,(tsd_minus(i+1)+y0))
            x1 = x2
         end do
         call ps_line_color(0.,0.,0.) !black

      end if
*
*       ********************************************************************************************************
*       ***** re-scale the main histogram and plot it again so that the overlay will be in the background: *****
*       ********************************************************************************************************
*
      max_bin_value = 0
      do i = 1,total_num_bins
         max_bin_value = MAX0(max_bin_value,(IH(i)))
      end do
      if(max_bin_value.eq.0)then !avoid a "divide by zero" situation
         scale_factor = 1.0
         goto 101
      end if
*
      scale_factor = height/FLOAT(max_bin_value)
*
 101  do i = 1, total_num_bins
         IH_scaled(i) = IFIX(FLOAT(IH(i))*scale_factor) !scale the main histogram
         IH1_scaled(i)=IH_scaled(i) !store the scaled histogram
      end do
      
      if (size(ICN).ge.2) then
         if((ICN(2).eq.sum_of_all_bins).and.(ICN(1).gt.0))then
            rate = (float(max_bin_value)*
     +           (1000.0/binwd))/ITAL(IDs(ICN(1))) !this plot is mirrored
         else
            rate = (float(max_bin_value)*(1000.0/binwd))/ICN(2) !calculate the rate in spikes/sec
         end if
      end if




c       rate = (float(max_bin_value)*(1000.0/binwd))/ICN(2)     !calculate the rate in spikes/sec
*
*
*
*       ************************************************************
*       *****   plot the data from the main histogram array:   *****
*       ************************************************************
*
      call ps_line_type(iSOLID)

*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
*
*       ***** special case of respiratory cycle duration *****  !
*       *****   distribution histogram (hist_type=DIS):  *****  !
*                                                               !
      if(hist_type.eq.'DIS')then !
         call ps_fill_color(0.,0.,0.) !
         call ps_perimeter_color(0.,0.,0.) !
         call ps_interior_style(INT_SOLID,1) !
         x1 = x0                !               
         do i = 1,total_num_bins !
            x2 = x1 + pixels_per_bin !
            call ps_rectangle(x1,y0,x2,(IH_scaled(i)+y0)) !
            x1 = x2             !
         end do                 !
         goto 120               !
      end if                    !
*                                                               !
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
*
      x1 = x0 
      do i = 1,(total_num_bins-1)
         x2 = x1 + pixels_per_bin
         call ps_move2d (x1,(IH_scaled(i)+y0))
         call ps_draw2d (x2,(IH_scaled(i)+y0))
         call ps_draw2d (x2,(IH_scaled(i+1)+y0))
         x1 = x2
      end do

      if (hist_type.eq.'CCH') then
         call sc_plot_ps (fildes, x0, y0, pixels_per_bin, scale_factor)
      end if
*
*       *************************************************************************************
*       *****   re-draw the axes, too (in case an overlay has obscured the x-axis):     *****
*       *************************************************************************************
*
      call ps_line_type (iSOLID)
      call ps_move2d (x0,(y0+height)) !Y axis
      call ps_draw2d (x0,y0)
      if(hist_type.ne.'DIF')then
         call ps_draw2d((x0+width),y0) !X axis
      end if
*
*
*       ***** if plotting an "enlarged" histogram, provide information on screen: *****
*
*
      if(info.eq.'y')then
*
         call ps_character_height(.030)
         call ps_character_width(.010)
*
         write (c_rec_num,'(I5)')rec_num
         if(show_control.eq.'y')
     +        write (c_rec_num_control,'(I5)')rec_num_control
         write (c_NHW,'(I4)')NHW
         write (c_binwd,'(F7.1)')binwd
         write (c_rate,'(F9.1)')rate
         write (c_max_bin_value,'(I7)')max_bin_value
         do i = 1,6             !convert to character format
            write (c_ICN(i),'(6I7)')ICN(i)
            call remove_all_blanks(c_ICN(i),LEN(c_ICN(i)))
         end do
         call strlength(c_ICN(1),LEN(c_ICN(1)),l_1)
         call strlength(c_ICN(2),LEN(c_ICN(2)),l_2)
         call strlength(c_ICN(3),LEN(c_ICN(3)),l_3)
         call strlength(c_ICN(4),LEN(c_ICN(4)),l_4)
         call strlength(c_ICN(5),LEN(c_ICN(5)),l_5)
         call strlength(c_ICN(6),LEN(c_ICN(1)),l_6)
*
         text='spike train filename:  '//BDT_FILE
         call strlength(text,LEN(text),l)
         call ps_text2d(5.,y0+height+160.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         call strlength(QDT_FILENAME,LEN(QDT_FILENAME),l)
         text='QDT filename:  '//QDT_FILENAME(1:l)
         call strlength(text,LEN(text),l)
         call ps_text2d(5.,y0+height+140.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         text='date of experiment:  '//date
         call strlength(text,LEN(text),l)
         call ps_text2d(5.,y0+height+120.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         text='recording #:  '//recording
         call strlength(text,LEN(text),l)
         call ps_text2d(5.,y0+height+100.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         if(show_control.eq.'y')then
            text='record # in .qdt file:  '//c_rec_num
            call strlength(text,LEN(text),l)
            call ps_text2d(5.,y0-40.,text(1:l)//char(0),
     +           PS_ANNOTATION_TEXT)
            text='record # of shift-control:  '//c_rec_num_control
            call strlength(text,LEN(text),l)
            call ps_text2d(x0+200.,y0-40.,text(1:l)//char(0),
     +           PS_ANNOTATION_TEXT)
         else
            text='record # in .qdt file:  '//c_rec_num
            call strlength(text,LEN(text),l)
            call ps_text2d(5.,y0-40.,text(1:l)//char(0),
     +           PS_ANNOTATION_TEXT)
         end if
*
         write (c_ITAL,'(I10)') ITAL(IDs(ICN(1)))
         call remove_all_blanks(c_ITAL,LEN(c_ITAL))
         call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)

         text='R: ID = '//c_ICN(1)(1:l_1)//
     +        ' (# events in spike train = '//
     +        c_ITAL(1:l_ITAL)//')'
         call strlength(text,LEN(text),l)
         call ps_text2d(5.,y0-60.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         call strlength(c_ICN(4),LEN(c_ICN(4)),l_4)
         write (c_sum_of_all_bins,'(I7)') sum_of_all_bins 
         call strlength(c_sum_of_all_bins,LEN(c_sum_of_all_bins),ls)
         write (c_ITAL,'(I10)') ITAL(IDs(ICN(3)))
         call remove_all_blanks(c_ITAL,LEN(c_ITAL))
         call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)

         text='T: ID = '//c_ICN(3)(1:l_3)//
     +        ' (# events in spike train = '//
     +        c_ITAL(1:l_ITAL)//')'
         call strlength(text,LEN(text),l)
         call ps_text2d(5.,y0-80.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         write (c_max,'(I7)') max_bin
         call remove_all_blanks(c_max,LEN(c_max))
         call strlength(c_max,LEN(c_max),l)
         write (c_min,'(I7)') min_bin
         call remove_all_blanks(c_min,LEN(c_min))
         call strlength(c_min,LEN(c_min),m)
         text='max bin = '//c_max(1:l)//'; min bin = '//
     +        c_min(1:m)//'; sum of all bins = '//
     +        c_sum_of_all_bins(1:ls)
         call strlength(text,LEN(text),l)
         call ps_text2d(5.,y0-100.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         text='binwidth:  '//c_binwd//' msec.'
         call strlength(text,LEN(text),l)
         call ps_text2d(5.,y0-120.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         text='-'//c_NHW
         call strlength(text,LEN(text),l)
         call ps_text2d(x0,y0-20.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         text=c_NHW//' msec'
         call strlength(text,LEN(text),l)
         call ps_text2d(x0+width-50.,y0-20.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         call ps_character_height(.060)
         call ps_character_width(.020)
         call strlength(title,LEN(title),l)
         call ps_text2d (x0-50.,(y0+height+60.),
     +        title(1:l)//char(0),PS_ANNOTATION_TEXT) !'REF>TAR'
*
         if(scaledup.eq.'n')then
            call ps_character_height(.036)
            call ps_character_width(.012)
            text=c_rate//' spikes/sec'
            call strlength(text,LEN(text),l)
            call ps_text2d(x0,(y0+height+20.),text(1:l)//char(0),
     +           PS_ANNOTATION_TEXT)
         else if(scaledup.eq.'y')then
            call ps_character_height(.030)
            call ps_character_width(.010)
            text = '# events/bin'
            call strlength(text,LEN(text),l)
            call ps_text2d(x0-150.,(y0+y0+height)/2.,
     +           text(1:l)//char(0),PS_ANNOTATION_TEXT)

            write (c_magnify,'(f5.2)') rmagnify
            call remove_all_blanks(c_magnify,LEN(c_magnify))
            call strlength(c_magnify,LEN(c_magnify),l)
            text = '(scaled up '//c_magnify(1:l)//'x)'
            call strlength(text,LEN(text),l)
            call ps_text2d(x0-45.,(y0+height+40.),
     +           text(1:l)//char(0),PS_ANNOTATION_TEXT)
            write (c_magnify,'(f5.2)') rmagnify
            call remove_all_blanks(c_magnify,LEN(c_magnify))
            call strlength(c_magnify,LEN(c_magnify),l)
            text = '(scaled up '//c_magnify(1:l)//'x)'
            call strlength(text,LEN(text),l)
            call ps_text2d(x0-45.,(y0+height+40.),
     +           text(1:l)//char(0),PS_ANNOTATION_TEXT)
         end if
      end if

      return
      end
      end module mod_ps_plot_1
